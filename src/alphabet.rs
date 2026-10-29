//! Alphabet recognition, equality definition and query profile (`Peq`) building.
//!
//! Port of `transformSequences`, `EqualityDefinition` and `buildPeq` from edlib.

use crate::myers::{ceil_div, WORD_SIZE};

/// Largest byte value (`MAX_UCHAR` in edlib).
const ALPHABET_CAPACITY: usize = 256;

/// Defines the equality relation on alphabet symbols.
pub(crate) struct EqualityDefinition {
    alphabet_length: usize,
    /// Dense `alphabet_length x alphabet_length` matrix.
    matrix: Vec<bool>,
}

impl EqualityDefinition {
    fn new(alphabet_length: usize) -> Self {
        let mut matrix = vec![false; alphabet_length * alphabet_length];
        for i in 0..alphabet_length {
            matrix[i * alphabet_length + i] = true;
        }
        EqualityDefinition {
            alphabet_length,
            matrix,
        }
    }

    fn add_pair(&mut self, first: usize, second: usize) {
        if first < self.alphabet_length && second < self.alphabet_length {
            self.matrix[first * self.alphabet_length + second] = true;
            self.matrix[second * self.alphabet_length + first] = true;
        }
    }

    /// `a` and `b` are indexes into the transformed sequences.
    #[inline]
    pub(crate) fn are_equal(&self, a: u8, b: u8) -> bool {
        self.matrix[a as usize * self.alphabet_length + b as usize]
    }
}

/// Sequences transformed into alphabet indexes, plus the alphabet metadata.
pub(crate) struct Alphabet {
    pub(crate) alphabet_length: usize,
    pub(crate) equality: EqualityDefinition,
    pub(crate) query: Vec<u8>,
    pub(crate) target: Vec<u8>,
}

impl Alphabet {
    pub(crate) fn new(query: &[u8], target: &[u8], additional_equalities: &[(u8, u8)]) -> Alphabet {
        let mut in_alphabet = [false; ALPHABET_CAPACITY];
        let mut letter_idx = [0u8; ALPHABET_CAPACITY];
        let mut alphabet = Vec::with_capacity(ALPHABET_CAPACITY);

        let mut query_transformed = vec![0u8; query.len()];
        for (slot, &c) in query_transformed.iter_mut().zip(query.iter()) {
            let c_idx = c as usize;
            if !in_alphabet[c_idx] {
                in_alphabet[c_idx] = true;
                letter_idx[c_idx] = alphabet.len() as u8;
                alphabet.push(c);
            }
            *slot = letter_idx[c_idx];
        }

        let mut target_transformed = vec![0u8; target.len()];
        for (slot, &c) in target_transformed.iter_mut().zip(target.iter()) {
            let c_idx = c as usize;
            if !in_alphabet[c_idx] {
                in_alphabet[c_idx] = true;
                letter_idx[c_idx] = alphabet.len() as u8;
                alphabet.push(c);
            }
            *slot = letter_idx[c_idx];
        }

        let alphabet_length = alphabet.len();
        let mut equality = EqualityDefinition::new(alphabet_length);
        for &(first, second) in additional_equalities {
            if in_alphabet[first as usize] && in_alphabet[second as usize] {
                equality.add_pair(
                    letter_idx[first as usize] as usize,
                    letter_idx[second as usize] as usize,
                );
            }
        }

        Alphabet {
            alphabet_length,
            equality,
            query: query_transformed,
            target: target_transformed,
        }
    }
}

/// Builds the `Peq` query profile.
///
/// `Peq` has dimensions `(alphabetLength + 1) x maxNumBlocks`. Bit `i` of
/// `Peq[s * maxNumBlocks + b]` is 1 if the `i`-th symbol of block `b` of the
/// query equals symbol `s`. The last symbol is a wildcard that matches
/// everything and is used for the padding of the last block.
pub(crate) fn build_peq(
    equality: &EqualityDefinition,
    query: &[u8],
    alphabet_length: usize,
) -> Vec<u64> {
    let max_num_blocks = ceil_div(query.len(), WORD_SIZE);
    let mut peq = vec![0u64; (alphabet_length + 1) * max_num_blocks];

    for symbol in 0..=alphabet_length {
        for b in 0..max_num_blocks {
            let idx = symbol * max_num_blocks + b;
            if symbol < alphabet_length {
                let mut word = 0u64;
                for r in (b * WORD_SIZE..(b + 1) * WORD_SIZE).rev() {
                    word <<= 1;
                    // The query is padded at the end with wildcard symbols.
                    if r >= query.len() || equality.are_equal(query[r], symbol as u8) {
                        word += 1;
                    }
                }
                peq[idx] = word;
            } else {
                // Last symbol is the wildcard, so it is all 1s.
                peq[idx] = u64::MAX;
            }
        }
    }

    peq
}
