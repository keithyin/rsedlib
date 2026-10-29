//! Alignment path reconstruction.
//!
//! Port of `obtainAlignment`, `obtainAlignmentHirschberg` and
//! `obtainAlignmentTraceback` from edlib. Small problems are solved with a full
//! traceback table, large ones with Hirschberg's divide and conquer.

use crate::alphabet::{build_peq, EqualityDefinition};
use crate::myers::{
    ceil_div, read_block, read_block_reverse, Block, EDOP_DELETE, EDOP_INSERT, EDOP_MATCH,
    EDOP_MISMATCH, HIGH_BIT_MASK, WORD_SIZE,
};
use crate::search::{myers_calc_edit_distance_nw, AlignmentData};

/// Maximum estimated memory (in bytes) for the traceback table. Bigger problems
/// are split with Hirschberg's algorithm.
const MAX_ALIGNMENT_DATA_SIZE: i64 = 1024 * 1024;

/// Finds one optimal alignment between `query` and `target`.
///
/// `r_query` and `r_target` are the reversed copies of `query` and `target`.
/// The returned vector contains edlib edit operations (`0 = match`, `1 = insertion
/// to target`, `2 = insertion to query`, `3 = mismatch`).
pub(crate) fn obtain_alignment(
    query: &[u8],
    r_query: &[u8],
    target: &[u8],
    r_target: &[u8],
    equality: &EqualityDefinition,
    alphabet_length: usize,
    best_score: i32,
) -> Result<Vec<u8>, String> {
    // Handle the special case when one of the sequences has length 0.
    if query.is_empty() || target.is_empty() {
        let op = if query.is_empty() {
            EDOP_DELETE
        } else {
            EDOP_INSERT
        };
        return Ok(vec![op; query.len() + target.len()]);
    }

    let max_num_blocks = ceil_div(query.len(), WORD_SIZE);
    let w = max_num_blocks * WORD_SIZE - query.len();

    // If the estimated memory consumption of the traceback table is small
    // enough, use it, otherwise use Hirschberg's algorithm.
    let per_cell = (2 * std::mem::size_of::<u64>() + std::mem::size_of::<i32>()) as i64;
    let alignment_data_size = per_cell * max_num_blocks as i64 * target.len() as i64
        + 2 * std::mem::size_of::<i32>() as i64 * target.len() as i64;

    if alignment_data_size < MAX_ALIGNMENT_DATA_SIZE {
        let peq = build_peq(equality, query, alphabet_length);
        let (_score, _end, align_data) = myers_calc_edit_distance_nw(
            &peq,
            w,
            max_num_blocks,
            query.len(),
            target,
            best_score,
            true,
            -1,
        );
        let align_data = align_data
            .ok_or_else(|| "Failed to reconstruct alignment: no alignment data".to_string())?;
        obtain_alignment_traceback(query.len(), target.len(), best_score, &align_data)
    } else {
        obtain_alignment_hirschberg(
            query,
            r_query,
            target,
            r_target,
            equality,
            alphabet_length,
            best_score,
        )
    }
}

/// Hirschberg's algorithm: splits the DP matrix in half, finds the row where the
/// optimal alignment crosses the split and solves both halves recursively.
#[allow(clippy::too_many_arguments)]
fn obtain_alignment_hirschberg(
    query: &[u8],
    r_query: &[u8],
    target: &[u8],
    r_target: &[u8],
    equality: &EqualityDefinition,
    alphabet_length: usize,
    best_score: i32,
) -> Result<Vec<u8>, String> {
    let max_num_blocks = ceil_div(query.len(), WORD_SIZE);
    let w = max_num_blocks * WORD_SIZE - query.len();
    let target_length = target.len();
    let query_length = query.len();

    let peq = build_peq(equality, query, alphabet_length);
    let r_peq = build_peq(equality, r_query, alphabet_length);

    // Divide the dynamic matrix into two halves, left and right.
    let left_half_width = target_length / 2;
    let right_half_width = target_length - left_half_width;

    // Calculate the left half.
    let (_s, _e, align_data_left) = myers_calc_edit_distance_nw(
        &peq,
        w,
        max_num_blocks,
        query_length,
        target,
        best_score,
        false,
        left_half_width as i32 - 1,
    );
    // Calculate the right half.
    let (_s, _e, align_data_right) = myers_calc_edit_distance_nw(
        &r_peq,
        w,
        max_num_blocks,
        query_length,
        r_target,
        best_score,
        false,
        right_half_width as i32 - 1,
    );

    let align_data_left = align_data_left
        .ok_or_else(|| "Failed to reconstruct alignment: missing left half data".to_string())?;
    let align_data_right = align_data_right
        .ok_or_else(|| "Failed to reconstruct alignment: missing right half data".to_string())?;

    // Unwrap the left half.
    let first_block_left = align_data_left.first_blocks[0];
    let last_block_left = align_data_left.last_blocks[0];
    let mut scores_left_length = (last_block_left - first_block_left + 1) as usize * WORD_SIZE;
    let mut scores_left = vec![0i32; scores_left_length];
    for block_idx in first_block_left..=last_block_left {
        let block = Block::new(
            align_data_left.ps[block_idx as usize],
            align_data_left.ms[block_idx as usize],
            align_data_left.scores[block_idx as usize],
        );
        let offset = (block_idx - first_block_left) as usize * WORD_SIZE;
        read_block(&block, &mut scores_left[offset..offset + WORD_SIZE]);
    }
    let scores_left_start_idx = first_block_left * WORD_SIZE as i32;
    // If the last block contains padding, shorten the length of scores by it.
    if last_block_left == max_num_blocks as i32 - 1 {
        scores_left_length -= w;
    }

    // Unwrap the right half (also reversed while unwrapping).
    let first_block_right = align_data_right.first_blocks[0];
    let last_block_right = align_data_right.last_blocks[0];
    let mut scores_right_length = (last_block_right - first_block_right + 1) as usize * WORD_SIZE;
    let mut scores_right = vec![0i32; scores_right_length];
    for block_idx in first_block_right..=last_block_right {
        let block = Block::new(
            align_data_right.ps[block_idx as usize],
            align_data_right.ms[block_idx as usize],
            align_data_right.scores[block_idx as usize],
        );
        let offset = (last_block_right - block_idx) as usize * WORD_SIZE;
        read_block_reverse(&block, &mut scores_right[offset..offset + WORD_SIZE]);
    }
    let mut scores_right_start_idx =
        query_length as i32 - (last_block_right + 1) * WORD_SIZE as i32;
    // If there is padding at the beginning of scoresRight (because of the
    // reversing), skip it.
    let mut scores_right_offset = 0usize;
    if scores_right_start_idx < 0 {
        scores_right_offset = w;
        scores_right_start_idx += w as i32;
        scores_right_length -= w;
    }

    //--------------------- Find the best move ----------------//
    // Find the query/row index of the cell in the left column which, together
    // with its lower right neighbour from the right column, gives the best score.
    let query_idx_left_start = std::cmp::max(scores_left_start_idx, scores_right_start_idx - 1);
    let query_idx_left_end = std::cmp::min(
        scores_left_start_idx + scores_left_length as i32 - 1,
        scores_right_start_idx + scores_right_length as i32 - 2,
    );

    let mut left_score = -1;
    let mut right_score = -1;
    let mut query_idx_left_alignment: i32 = -1;
    let mut found = false;

    let mut query_idx = query_idx_left_start;
    while query_idx <= query_idx_left_end {
        left_score = scores_left[(query_idx - scores_left_start_idx) as usize];
        right_score =
            scores_right[scores_right_offset + (query_idx + 1 - scores_right_start_idx) as usize];
        if left_score + right_score == best_score {
            query_idx_left_alignment = query_idx;
            found = true;
            break;
        }
        query_idx += 1;
    }

    // Check boundary cells.
    if !found && scores_left_start_idx == 0 && scores_right_start_idx == 0 {
        left_score = left_half_width as i32;
        right_score = scores_right[scores_right_offset];
        if left_score + right_score == best_score {
            query_idx_left_alignment = -1;
            found = true;
        }
    }
    if !found
        && scores_left_start_idx + scores_left_length as i32 == query_length as i32
        && scores_right_start_idx + scores_right_length as i32 == query_length as i32
    {
        left_score = scores_left[scores_left_length - 1];
        right_score = right_half_width as i32;
        if left_score + right_score == best_score {
            query_idx_left_alignment = query_length as i32 - 1;
            found = true;
        }
    }

    if !found {
        // If there was no move that is part of an optimal alignment, then there
        // is no such alignment or the given bestScore is not correct.
        return Err("Failed to reconstruct alignment: optimal split not found".to_string());
    }
    //----------------------------------------------------------//

    let ul_height = (query_idx_left_alignment + 1) as usize;
    let lr_height = query_length - ul_height;
    let ul_width = left_half_width;
    let lr_width = right_half_width;

    // Upper left of the left half, and lower right of the right half.
    let mut alignment = obtain_alignment(
        &query[..ul_height],
        &r_query[lr_height..],
        &target[..ul_width],
        &r_target[lr_width..],
        equality,
        alphabet_length,
        left_score,
    )?;
    let lr_alignment = obtain_alignment(
        &query[ul_height..],
        &r_query[..lr_height],
        &target[ul_width..],
        &r_target[..lr_width],
        equality,
        alphabet_length,
        right_score,
    )?;
    alignment.extend_from_slice(&lr_alignment);

    Ok(alignment)
}

/// Reconstructs the alignment by walking back through the stored DP matrix.
fn obtain_alignment_traceback(
    query_length: usize,
    target_length: usize,
    best_score: i32,
    align_data: &AlignmentData,
) -> Result<Vec<u8>, String> {
    let max_num_blocks = ceil_div(query_length, WORD_SIZE);
    let w = max_num_blocks * WORD_SIZE - query_length;

    let mut alignment: Vec<u8> = Vec::with_capacity(query_length + target_length);

    let mut c: i32 = target_length as i32 - 1;
    let mut b: i32 = max_num_blocks as i32 - 1;
    let mut curr_score = best_score;
    let mut l_score: i32 = -1;
    let mut u_score: i32 = -1;
    let mut ul_score: i32 = -1;

    let mut curr_p = align_data.ps[align_data.idx(c as usize, b)];
    let mut curr_m = align_data.ms[align_data.idx(c as usize, b)];

    let mut there_is_left_block = c > 0
        && b >= align_data.first_blocks[(c - 1) as usize]
        && b <= align_data.last_blocks[(c - 1) as usize];
    let mut l_p = 0u64;
    let mut l_m = 0u64;
    if there_is_left_block {
        l_p = align_data.ps[align_data.idx((c - 1) as usize, b)];
        l_m = align_data.ms[align_data.idx((c - 1) as usize, b)];
    }

    curr_p <<= w;
    curr_m <<= w;
    let mut block_pos: i32 = WORD_SIZE as i32 - w as i32 - 1;

    loop {
        if c == 0 {
            there_is_left_block = true;
            l_score = b * WORD_SIZE as i32 + block_pos + 1;
            ul_score = l_score - 1;
        }

        //---------- Calculate scores ---------//
        if l_score == -1 && there_is_left_block {
            l_score = align_data.scores[align_data.idx((c - 1) as usize, b)];
            for _ in 0..WORD_SIZE as i32 - block_pos - 1 {
                if l_p & HIGH_BIT_MASK != 0 {
                    l_score -= 1;
                }
                if l_m & HIGH_BIT_MASK != 0 {
                    l_score += 1;
                }
                l_p <<= 1;
                l_m <<= 1;
            }
        }
        if ul_score == -1 {
            if l_score != -1 {
                ul_score = l_score;
                if l_p & HIGH_BIT_MASK != 0 {
                    ul_score -= 1;
                }
                if l_m & HIGH_BIT_MASK != 0 {
                    ul_score += 1;
                }
            } else if c > 0
                && b > align_data.first_blocks[(c - 1) as usize]
                && b - 1 <= align_data.last_blocks[(c - 1) as usize]
            {
                // The upper left cell is the last cell in the block and the block
                // to the left is not in the band, so lScore is -1.
                ul_score = align_data.scores[align_data.idx((c - 1) as usize, b - 1)];
            }
        }
        if u_score == -1 {
            u_score = curr_score;
            if curr_p & HIGH_BIT_MASK != 0 {
                u_score -= 1;
            }
            if curr_m & HIGH_BIT_MASK != 0 {
                u_score += 1;
            }
            curr_p <<= 1;
            curr_m <<= 1;
        }
        //-------------------------------------//

        //-------------- Move --------------//
        // Move up - insertion to target - deletion from query.
        if u_score != -1 && u_score + 1 == curr_score {
            curr_score = u_score;
            l_score = ul_score;
            u_score = -1;
            ul_score = -1;
            if block_pos == 0 {
                // If entering a new (upper) block.
                if b == 0 {
                    // If there are no cells above (only boundary cells).
                    alignment.push(EDOP_INSERT);
                    alignment.resize(alignment.len() + (c + 1) as usize, EDOP_DELETE);
                    break;
                } else {
                    block_pos = WORD_SIZE as i32 - 1;
                    b -= 1;
                    curr_p = align_data.ps[align_data.idx(c as usize, b)];
                    curr_m = align_data.ms[align_data.idx(c as usize, b)];
                    if c > 0
                        && b >= align_data.first_blocks[(c - 1) as usize]
                        && b <= align_data.last_blocks[(c - 1) as usize]
                    {
                        there_is_left_block = true;
                        l_p = align_data.ps[align_data.idx((c - 1) as usize, b)];
                        l_m = align_data.ms[align_data.idx((c - 1) as usize, b)];
                    } else {
                        there_is_left_block = false;
                    }
                }
            } else {
                block_pos -= 1;
                l_p <<= 1;
                l_m <<= 1;
            }
            alignment.push(EDOP_INSERT);
        }
        // Move left - deletion from target - insertion to query.
        else if l_score != -1 && l_score + 1 == curr_score {
            curr_score = l_score;
            u_score = ul_score;
            l_score = -1;
            ul_score = -1;
            c -= 1;
            if c == -1 {
                // If there are no cells to the left (only boundary cells).
                alignment.push(EDOP_DELETE);
                let num_up = b * WORD_SIZE as i32 + block_pos + 1;
                alignment.resize(alignment.len() + num_up as usize, EDOP_INSERT);
                break;
            }
            curr_p = l_p;
            curr_m = l_m;
            if c > 0
                && b >= align_data.first_blocks[(c - 1) as usize]
                && b <= align_data.last_blocks[(c - 1) as usize]
            {
                there_is_left_block = true;
                l_p = align_data.ps[align_data.idx((c - 1) as usize, b)];
                l_m = align_data.ms[align_data.idx((c - 1) as usize, b)];
            } else if c == 0 {
                // If there are no cells to the left (only boundary cells).
                there_is_left_block = true;
                l_score = b * WORD_SIZE as i32 + block_pos + 1;
                ul_score = l_score - 1;
            } else {
                there_is_left_block = false;
            }
            alignment.push(EDOP_DELETE);
        }
        // Move up left - (mis)match.
        else if ul_score != -1 {
            let move_code = if ul_score == curr_score {
                EDOP_MATCH
            } else {
                EDOP_MISMATCH
            };
            curr_score = ul_score;
            u_score = -1;
            l_score = -1;
            ul_score = -1;
            c -= 1;
            if c == -1 {
                // If there are no cells to the left (only boundary cells).
                alignment.push(move_code);
                let num_up = b * WORD_SIZE as i32 + block_pos;
                alignment.resize(alignment.len() + num_up as usize, EDOP_INSERT);
                break;
            }
            if block_pos == 0 {
                // If entering the upper left block.
                if b == 0 {
                    // If there are no more cells above (only boundary cells).
                    alignment.push(move_code);
                    alignment.resize(alignment.len() + (c + 1) as usize, EDOP_DELETE);
                    break;
                }
                block_pos = WORD_SIZE as i32 - 1;
                b -= 1;
                curr_p = align_data.ps[align_data.idx(c as usize, b)];
                curr_m = align_data.ms[align_data.idx(c as usize, b)];
            } else {
                // If entering the left block.
                block_pos -= 1;
                curr_p = l_p;
                curr_m = l_m;
                curr_p <<= 1;
                curr_m <<= 1;
            }
            // Set new left block.
            if c > 0
                && b >= align_data.first_blocks[(c - 1) as usize]
                && b <= align_data.last_blocks[(c - 1) as usize]
            {
                there_is_left_block = true;
                l_p = align_data.ps[align_data.idx((c - 1) as usize, b)];
                l_m = align_data.ms[align_data.idx((c - 1) as usize, b)];
            } else if c == 0 {
                // If there are no cells to the left (only boundary cells).
                there_is_left_block = true;
                l_score = b * WORD_SIZE as i32 + block_pos + 1;
                ul_score = l_score - 1;
            } else {
                there_is_left_block = false;
            }
            alignment.push(move_code);
        } else {
            // Reached the end - finished!
            break;
        }
        //----------------------------------//
    }

    alignment.reverse();
    Ok(alignment)
}
