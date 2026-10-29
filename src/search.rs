//! Banded bit-parallel Myers search.
//!
//! Faithful port of `myersCalcEditDistanceSemiGlobal` and
//! `myersCalcEditDistanceNW` from edlib, including the Ukkonen band reduction
//! heuristics.

use crate::myers::{
    all_block_cells_larger, block_cell_values, calculate_block, ceil_div, ceil_div_i32, Block,
    WORD_SIZE,
};

/// Alignment mode, mirrors `EdlibAlignMode`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Mode {
    /// Global alignment (Needleman-Wunsch).
    Nw,
    /// Prefix alignment, gaps after the query are free.
    Shw,
    /// Infix alignment, gaps before and after the query are free.
    Hw,
}

/// Each `STRONG_REDUCE_NUM` column is reduced in a more expensive but more
/// efficient way.
const STRONG_REDUCE_NUM: i32 = 2048;

/// Data needed to reconstruct an alignment (only filled when requested).
pub(crate) struct AlignmentData {
    pub(crate) ps: Vec<u64>,
    pub(crate) ms: Vec<u64>,
    pub(crate) scores: Vec<i32>,
    pub(crate) first_blocks: Vec<i32>,
    pub(crate) last_blocks: Vec<i32>,
    pub(crate) max_num_blocks: usize,
}

impl AlignmentData {
    pub(crate) fn new(max_num_blocks: usize, target_length: usize) -> Self {
        AlignmentData {
            ps: vec![0u64; max_num_blocks * target_length],
            ms: vec![0u64; max_num_blocks * target_length],
            scores: vec![0i32; max_num_blocks * target_length],
            first_blocks: vec![0i32; target_length],
            last_blocks: vec![0i32; target_length],
            max_num_blocks,
        }
    }

    #[inline]
    pub(crate) fn idx(&self, column: usize, block: i32) -> usize {
        column * self.max_num_blocks + block as usize
    }
}

#[inline]
fn k_plus_word(k: i32) -> i32 {
    k.saturating_add(WORD_SIZE as i32)
}

/// Myers' bit-vector algorithm for the semi-global methods (HW and SHW).
///
/// Returns the best score (`-1` if no score below `k` was found) and the
/// zero-based positions in `target` where that score was found.
pub(crate) fn myers_calc_edit_distance_semi_global(
    peq: &[u64],
    w: usize,
    max_num_blocks: usize,
    query_length: usize,
    target: &[u8],
    k_in: i32,
    mode: Mode,
) -> (i32, Vec<i32>) {
    let mut k = k_in;
    let max_blocks = max_num_blocks as i32;

    let mut first_block: i32 = 0;
    let k_for_band = if k < 0 { 0 } else { k };
    let mut last_block: i32 = std::cmp::min(
        ceil_div(k_for_band.saturating_add(1) as usize, WORD_SIZE) as i32,
        max_blocks,
    ) - 1;

    let mut blocks = vec![Block::default(); max_num_blocks];

    // For HW the solution is never larger than queryLength.
    if mode == Mode::Hw {
        k = std::cmp::min(query_length as i32, k);
    }

    let mut b = 0;
    while b <= last_block {
        let block = &mut blocks[b as usize];
        block.score = (b + 1) * WORD_SIZE as i32;
        block.p = u64::MAX;
        block.m = 0;
        b += 1;
    }

    let mut best_score: i32 = -1;
    let mut positions: Vec<i32> = Vec::new();
    let start_hout = if mode == Mode::Hw { 0 } else { 1 };
    let mut bl: i32 = 0;

    for (c, &target_char) in target.iter().enumerate() {
        let peq_c_base = target_char as usize * max_num_blocks;

        //----------------------- Calculate column -------------------------//
        let mut hout = start_hout;
        bl = first_block;
        let mut b = first_block;
        while b <= last_block {
            let block = &mut blocks[bl as usize];
            let (p, m, h) = calculate_block(block.p, block.m, peq[peq_c_base + b as usize], hout);
            block.p = p;
            block.m = m;
            block.score += h;
            hout = h;
            bl += 1;
            b += 1;
        }
        bl -= 1;
        //------------------------------------------------------------------//

        //---------- Adjust number of blocks according to Ukkonen ----------//
        if last_block < max_blocks - 1
            && blocks[bl as usize].score - hout <= k
            && ((peq[peq_c_base + (bl + 1) as usize] & 1) != 0 || hout < 0)
        {
            // If score of the left block is not too big, calculate one more block.
            last_block += 1;
            bl += 1;
            let previous_score = blocks[(bl - 1) as usize].score;
            let block = &mut blocks[bl as usize];
            block.p = u64::MAX;
            block.m = 0;
            let (p, m, h) = calculate_block(
                block.p,
                block.m,
                peq[peq_c_base + last_block as usize],
                hout,
            );
            block.p = p;
            block.m = m;
            block.score = previous_score - hout + WORD_SIZE as i32 + h;
        } else {
            while last_block >= first_block && blocks[bl as usize].score >= k_plus_word(k) {
                last_block -= 1;
                bl -= 1;
            }
        }

        // Every some columns, do some expensive but also more efficient block reducing.
        if c as i32 % STRONG_REDUCE_NUM == 0 {
            while last_block >= 0
                && last_block >= first_block
                && all_block_cells_larger(&blocks[bl as usize], k)
            {
                last_block -= 1;
                bl -= 1;
            }
        }
        // For HW there still may be a solution in the next column, so the
        // first block is always a candidate.
        if mode == Mode::Hw && last_block == -1 {
            last_block += 1;
            bl += 1;
        }

        // Reduce the band by increasing the first block if possible. Not applicable to HW.
        if mode != Mode::Hw {
            while first_block <= last_block && blocks[first_block as usize].score >= k_plus_word(k)
            {
                first_block += 1;
            }
            if c as i32 % STRONG_REDUCE_NUM == 0 {
                while first_block <= last_block
                    && all_block_cells_larger(&blocks[first_block as usize], k)
                {
                    first_block += 1;
                }
            }
        }

        // If the band stops to exist, finish.
        if last_block < first_block {
            return (best_score, positions);
        }
        //------------------------------------------------------------------//

        //------------------------- Update best score ----------------------//
        if last_block == max_blocks - 1 {
            let col_score = blocks[bl as usize].score;
            // Scores > k do not have correct values, but are certainly > k.
            if col_score <= k && (best_score == -1 || col_score <= best_score) {
                if col_score != best_score {
                    positions.clear();
                    best_score = col_score;
                    // Look only for equal or better scores than the best found so far.
                    k = best_score;
                }
                // NOTE: the score found in column c is actually the score of column c-W.
                positions.push(c as i32 - w as i32);
            }
        }
        //------------------------------------------------------------------//
    }

    // Obtain results for the last W columns from the last column.
    if last_block == max_blocks - 1 {
        let block_scores = block_cell_values(&blocks[bl as usize]);
        for i in 0..w {
            let col_score = block_scores[i + 1];
            if col_score <= k && (best_score == -1 || col_score <= best_score) {
                if col_score != best_score {
                    positions.clear();
                    k = col_score;
                    best_score = col_score;
                }
                positions.push(target.len() as i32 - w as i32 + i as i32);
            }
        }
    }

    (best_score, positions)
}

/// Myers' bit-vector algorithm for the global (NW) method.
///
/// Returns `(best_score, position, alignment_data)`; score and position are `-1`
/// when no score below `k` was found.
#[allow(clippy::too_many_arguments)]
pub(crate) fn myers_calc_edit_distance_nw(
    peq: &[u64],
    w: usize,
    max_num_blocks: usize,
    query_length: usize,
    target: &[u8],
    k_in: i32,
    find_alignment: bool,
    target_stop_position: i32,
) -> (i32, i32, Option<AlignmentData>) {
    if target_stop_position > -1 && find_alignment {
        // They can not be both set at the same time!
        return (-1, -1, None);
    }

    let target_length = target.len();
    let t_len = target_length as i32;
    let q_len = query_length as i32;
    let max_blocks = max_num_blocks as i32;

    if k_in < (t_len - q_len).abs() {
        return (-1, -1, None);
    }

    let mut k = std::cmp::min(k_in, std::cmp::max(q_len, t_len));

    let mut first_block: i32 = 0;
    let inner = std::cmp::min(k, (k + q_len - t_len) / 2);
    let mut last_block: i32 =
        std::cmp::min(max_blocks, ceil_div_i32(inner + 1, WORD_SIZE as i32)) - 1;

    let mut blocks = vec![Block::default(); max_num_blocks];

    let mut b = 0;
    while b <= last_block {
        let block = &mut blocks[b as usize];
        block.score = (b + 1) * WORD_SIZE as i32;
        block.p = u64::MAX;
        block.m = 0;
        b += 1;
    }

    let mut align_data = if find_alignment {
        Some(AlignmentData::new(max_num_blocks, target_length))
    } else if target_stop_position > -1 {
        Some(AlignmentData::new(max_num_blocks, 1))
    } else {
        None
    };

    let mut bl: i32;
    for (c, &target_char) in target.iter().enumerate() {
        let c_i32 = c as i32;
        let peq_c_base = target_char as usize * max_num_blocks;

        //----------------------- Calculate column -------------------------//
        let mut hout = 1;
        bl = first_block;
        let mut b = first_block;
        while b <= last_block {
            let block = &mut blocks[bl as usize];
            let (p, m, h) = calculate_block(block.p, block.m, peq[peq_c_base + b as usize], hout);
            block.p = p;
            block.m = m;
            block.score += h;
            hout = h;
            bl += 1;
            b += 1;
        }
        bl -= 1;
        //------------------------------------------------------------------//

        // Update k. Only done at the end of the column.
        // NOTE: W is added in the last block because the result is actually from
        // W cells to the left and W cells up.
        k = std::cmp::min(
            k,
            blocks[bl as usize].score
                + std::cmp::max(
                    t_len - c_i32 - 1,
                    q_len - ((1 + last_block) * WORD_SIZE as i32 - 1) - 1,
                )
                + if last_block == max_blocks - 1 {
                    w as i32
                } else {
                    0
                },
        );

        //---------- Adjust number of blocks according to Ukkonen ----------//
        //--- Adjust last block ---//
        if last_block + 1 < max_blocks
            && (last_block + 1) * WORD_SIZE as i32 - 1
                <= k - blocks[bl as usize].score + 2 * WORD_SIZE as i32 - 2 - t_len + c_i32 + q_len
        {
            last_block += 1;
            bl += 1;
            let previous_score = blocks[(bl - 1) as usize].score;
            let block = &mut blocks[bl as usize];
            block.p = u64::MAX;
            block.m = 0;
            let (p, m, new_hout) = calculate_block(
                block.p,
                block.m,
                peq[peq_c_base + last_block as usize],
                hout,
            );
            block.p = p;
            block.m = m;
            block.score = previous_score - hout + WORD_SIZE as i32 + new_hout;
        }

        // While block is out of band, move one block up.
        while last_block >= first_block
            && (blocks[bl as usize].score >= k_plus_word(k)
                || ((last_block + 1) * WORD_SIZE as i32 - 1
                    > k - blocks[bl as usize].score + 2 * WORD_SIZE as i32 - 2 - t_len
                        + c_i32
                        + q_len
                        + 1))
        {
            last_block -= 1;
            bl -= 1;
        }
        //-------------------------//

        //--- Adjust first block ---//
        while first_block <= last_block
            && (blocks[first_block as usize].score >= k_plus_word(k)
                || ((first_block + 1) * WORD_SIZE as i32 - 1
                    < blocks[first_block as usize].score - k - t_len + q_len + c_i32))
        {
            first_block += 1;
        }
        //--------------------------//

        // Every some columns do a more expensive but more efficient reduction.
        if c_i32 % STRONG_REDUCE_NUM == 0 {
            while last_block >= first_block {
                let scores = block_cell_values(&blocks[bl as usize]);
                let num_cells = if last_block == max_blocks - 1 {
                    WORD_SIZE as i32 - w as i32
                } else {
                    WORD_SIZE as i32
                };
                let mut r = last_block * WORD_SIZE as i32 + num_cells - 1;
                let mut reduce = true;
                for &score in scores.iter().skip((WORD_SIZE as i32 - num_cells) as usize) {
                    if score <= k && r <= k - score - t_len + c_i32 + q_len + 1 {
                        reduce = false;
                        break;
                    }
                    r -= 1;
                }
                if !reduce {
                    break;
                }
                last_block -= 1;
                bl -= 1;
            }

            while first_block <= last_block {
                let scores = block_cell_values(&blocks[first_block as usize]);
                let num_cells = if first_block == max_blocks - 1 {
                    WORD_SIZE as i32 - w as i32
                } else {
                    WORD_SIZE as i32
                };
                let mut r = first_block * WORD_SIZE as i32 + num_cells - 1;
                let mut reduce = true;
                for &score in scores.iter().skip((WORD_SIZE as i32 - num_cells) as usize) {
                    if score <= k && r >= score - k - t_len + c_i32 + q_len {
                        reduce = false;
                        break;
                    }
                    r -= 1;
                }
                if !reduce {
                    break;
                }
                first_block += 1;
            }
        }

        // If the band stops to exist, finish.
        if last_block < first_block {
            return (-1, -1, None);
        }
        //------------------------------------------------------------------//

        //---- Save column so it can be used for reconstruction ----//
        if find_alignment {
            if let Some(data) = align_data.as_mut() {
                let mut b = first_block;
                while b <= last_block {
                    let idx = data.idx(c, b);
                    let block = &blocks[b as usize];
                    data.ps[idx] = block.p;
                    data.ms[idx] = block.m;
                    data.scores[idx] = block.score;
                    b += 1;
                }
                data.first_blocks[c] = first_block;
                data.last_blocks[c] = last_block;
            }
        }
        //----------------------------------------------------------//

        //---- If this is the stop column, save it and finish ----//
        if c_i32 == target_stop_position {
            if let Some(data) = align_data.as_mut() {
                let mut b = first_block;
                while b <= last_block {
                    let block = &blocks[b as usize];
                    data.ps[b as usize] = block.p;
                    data.ms[b as usize] = block.m;
                    data.scores[b as usize] = block.score;
                    b += 1;
                }
                data.first_blocks[0] = first_block;
                data.last_blocks[0] = last_block;
            }
            return (-1, target_stop_position, align_data);
        }
        //----------------------------------------------------//
    }

    if last_block == max_blocks - 1 {
        // Obtain the best score from the block: the query is padded with W cells.
        let best_score = block_cell_values(&blocks[last_block as usize])[w];
        if best_score <= k {
            return (best_score, t_len - 1, align_data);
        }
    }

    (-1, -1, align_data)
}
