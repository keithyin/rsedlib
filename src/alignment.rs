//! Orchestration of a single alignment request.
//!
//! Port of the `edlibAlign` entry point: sequence transformation, k auto
//! adjustment, location lookup and (optional) alignment reconstruction.

use crate::alphabet::{build_peq, Alphabet};
use crate::myers::{ceil_div, WORD_SIZE};
use crate::param::{AlignMode, AlignTask, EdlibAlignParam};
use crate::search::{myers_calc_edit_distance_nw, myers_calc_edit_distance_semi_global, Mode};
use crate::traceback::obtain_alignment;

/// Result of an alignment, still in edlib's internal representation.
pub(crate) struct AlignOutcome {
    pub(crate) edit_distance: i32,
    pub(crate) alphabet_length: i32,
    /// Zero-based end positions in the target (`-1` is possible in HW mode).
    pub(crate) end_locations: Vec<i32>,
    /// Zero-based start positions in the target. Empty when not calculated
    /// (which the public API reports as `0`, mirroring the C wrapper).
    pub(crate) start_locations: Vec<i32>,
    /// Alignment operations, only present for the `Path` task.
    pub(crate) alignment: Option<Vec<u8>>,
}

/// Aligns `query` against `target` following edlib's semantics.
pub(crate) fn align(
    query: &[u8],
    target: &[u8],
    param: &EdlibAlignParam,
) -> Result<AlignOutcome, String> {
    if query.len() > i32::MAX as usize || target.len() > i32::MAX as usize {
        return Err("Sequence length exceeds i32::MAX".to_string());
    }

    let additional_equalities: Vec<(u8, u8)> = param
        .additional_eq_pairs()
        .iter()
        .map(|pair| (pair.first, pair.second))
        .collect();
    let alphabet = Alphabet::new(query, target, &additional_equalities);

    let mode = match param.mode() {
        AlignMode::Global => Mode::Nw,
        AlignMode::Prefix => Mode::Shw,
        AlignMode::Infix => Mode::Hw,
    };
    let task = param.task();

    let mut outcome = AlignOutcome {
        edit_distance: -1,
        alphabet_length: alphabet.alphabet_length as i32,
        end_locations: Vec::new(),
        start_locations: Vec::new(),
        alignment: None,
    };

    // Handle the special situation when at least one of the sequences is empty.
    if query.is_empty() || target.is_empty() {
        match mode {
            Mode::Nw => {
                outcome.edit_distance = std::cmp::max(query.len(), target.len()) as i32;
                outcome.end_locations = vec![target.len() as i32 - 1];
            }
            Mode::Shw | Mode::Hw => {
                outcome.edit_distance = query.len() as i32;
                outcome.end_locations = vec![-1];
            }
        }
        return Ok(outcome);
    }

    /*--------------------- INITIALIZATION ------------------*/
    let max_num_blocks = ceil_div(query.len(), WORD_SIZE);
    // Number of redundant cells in the last block.
    let w = max_num_blocks * WORD_SIZE - query.len();
    let peq = build_peq(
        &alphabet.equality,
        &alphabet.query,
        alphabet.alphabet_length,
    );
    /*-------------------------------------------------------*/

    /*------------------ MAIN CALCULATION -------------------*/
    let mut dynamic_k = false;
    let mut k = param.k();
    if k < 0 {
        // If a valid k is not given, auto-adjust k until a solution is found.
        dynamic_k = true;
        k = WORD_SIZE as i32;
    }

    loop {
        match mode {
            Mode::Nw => {
                let (best_score, _position, _align_data) = myers_calc_edit_distance_nw(
                    &peq,
                    w,
                    max_num_blocks,
                    query.len(),
                    &alphabet.target,
                    k,
                    false,
                    -1,
                );
                outcome.edit_distance = best_score;
            }
            Mode::Shw | Mode::Hw => {
                let (best_score, positions) = myers_calc_edit_distance_semi_global(
                    &peq,
                    w,
                    max_num_blocks,
                    query.len(),
                    &alphabet.target,
                    k,
                    mode,
                );
                outcome.edit_distance = best_score;
                outcome.end_locations = positions;
            }
        }
        k = k.saturating_mul(2);
        if !(dynamic_k && outcome.edit_distance == -1) {
            break;
        }
    }

    if outcome.edit_distance < 0 {
        // No solution within k: edlib reports edit distance -1 and no locations.
        outcome.end_locations.clear();
        return Ok(outcome);
    }

    // If NW mode, set the end location explicitly.
    if mode == Mode::Nw {
        outcome.end_locations = vec![target.len() as i32 - 1];
    }

    // Find starting locations.
    if task == AlignTask::Locations || task == AlignTask::Path {
        let num_locations = outcome.end_locations.len();
        outcome.start_locations = vec![0; num_locations];
        if mode == Mode::Hw {
            // For HW the start locations have to be calculated explicitly.
            let r_target: Vec<u8> = alphabet.target.iter().rev().copied().collect();
            let r_query: Vec<u8> = alphabet.query.iter().rev().copied().collect();
            // Peq for the reversed query.
            let r_peq = build_peq(&alphabet.equality, &r_query, alphabet.alphabet_length);
            for i in 0..num_locations {
                let end_location = outcome.end_locations[i];
                if end_location == -1 {
                    // Query starts before the target; edlib reports 0 here.
                    outcome.start_locations[i] = 0;
                } else {
                    let start = (target.len() as i32 - end_location - 1).max(0) as usize;
                    let end = (start + end_location as usize + 1).min(r_target.len());
                    let reversed_prefix = &r_target[start.min(end)..end];
                    let (_best_score, positions) = myers_calc_edit_distance_semi_global(
                        &r_peq,
                        w,
                        max_num_blocks,
                        query.len(),
                        reversed_prefix,
                        outcome.edit_distance,
                        Mode::Shw,
                    );
                    // Taking the last location as start ensures that the alignment
                    // does not start with insertions if it can start with
                    // mismatches instead.
                    let last_position = positions.last().copied().unwrap_or(0);
                    outcome.start_locations[i] = end_location - last_position;
                }
            }
        }
        // For NW and SHW the start location is always 0.
    }

    // Find the alignment: everything comes down to finding the alignment for NW.
    if task == AlignTask::Path {
        let start = outcome.start_locations[0];
        let end = outcome.end_locations[0];

        let aln_start = start.max(0) as usize;
        let aln_end_exclusive = (end + 1).max(0) as usize;
        let aln_end_exclusive = aln_end_exclusive.max(aln_start).min(alphabet.target.len());
        let aln_start = aln_start.min(aln_end_exclusive);

        let aln_target = &alphabet.target[aln_start..aln_end_exclusive];
        let r_aln_target: Vec<u8> = aln_target.iter().rev().copied().collect();
        let r_query: Vec<u8> = alphabet.query.iter().rev().copied().collect();

        let alignment = obtain_alignment(
            &alphabet.query,
            &r_query,
            aln_target,
            &r_aln_target,
            &alphabet.equality,
            alphabet.alphabet_length,
            outcome.edit_distance,
        )?;
        outcome.alignment = Some(alignment);
    }
    /*-------------------------------------------------------*/

    Ok(outcome)
}
