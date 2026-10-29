//! Conversion of an alignment (sequence of edit operations) into a CIGAR string.
//!
//! Port of `edlibAlignmentToCigar`.

use crate::myers::EDOP_MISMATCH;

/// Builds a CIGAR string from an alignment.
///
/// In the standard format matches and mismatches are both reported as `M`, in the
/// extended format they are reported as `=` and `X`. Returns `None` if the
/// alignment contains an invalid operation.
pub(crate) fn alignment_to_cigar(alignment: &[u8], extended: bool) -> Option<String> {
    if alignment.is_empty() {
        return Some(String::new());
    }

    // Maps a move code from the alignment to the char in the CIGAR.
    //                 0    1    2    3
    let mut code_to_char = *b"=IDX";
    if !extended {
        code_to_char[0] = b'M';
        code_to_char[3] = b'M';
    }

    let mut cigar = String::with_capacity(alignment.len() + 8);
    let mut last_move: u8 = 0; // 0 if there was no previous move.
    let mut num_of_same_moves: usize = 0;

    for i in 0..=alignment.len() {
        if i < alignment.len() && alignment[i] > EDOP_MISMATCH {
            // Invalid alignment value.
            return None;
        }
        // If a new sequence of same moves started.
        if i == alignment.len()
            || (code_to_char[alignment[i] as usize] != last_move && last_move != 0)
        {
            // Write the number of moves.
            cigar.push_str(&num_of_same_moves.to_string());
            // Write the code of the move.
            cigar.push(last_move as char);
            if i < alignment.len() {
                num_of_same_moves = 0;
            }
        }
        if i < alignment.len() {
            last_move = code_to_char[alignment[i] as usize];
            num_of_same_moves += 1;
        }
    }

    Some(cigar)
}
