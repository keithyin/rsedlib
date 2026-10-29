//! Shared helpers for integration tests: a naive O(n*m) reference implementation
//! plus a small deterministic pseudo random generator (no external crates).

#![allow(dead_code)]

/// Deterministic linear congruential generator.
pub struct Lcg(u64);

impl Lcg {
    pub fn new(seed: u64) -> Self {
        Lcg(seed ^ 0x9E37_79B9_7F4A_7C15)
    }

    pub fn next_u32(&mut self) -> u32 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 33) as u32
    }

    pub fn next_below(&mut self, n: usize) -> usize {
        (self.next_u32() as usize) % n
    }

    pub fn sequence(&mut self, len: usize, alphabet: &[u8]) -> Vec<u8> {
        (0..len)
            .map(|_| alphabet[self.next_below(alphabet.len())])
            .collect()
    }
}

/// Mode of the reference implementation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RefMode {
    /// Global (gaps before and after the query are penalized).
    Nw,
    /// Prefix (gaps after the query are free).
    Shw,
    /// Infix (gaps before and after the query are free).
    Hw,
}

/// Naive dynamic programming reference implementation.
///
/// Returns the last DP row, i.e. `row[j]` is the edit distance between the whole
/// query and `target[..j]`.
pub fn dp_last_row(query: &[u8], target: &[u8], mode: RefMode) -> Vec<i32> {
    let t = target.len();
    let free_prefix = mode == RefMode::Hw;
    let mut prev: Vec<i32> = if free_prefix {
        vec![0; t + 1]
    } else {
        (0..=t as i32).collect()
    };

    for i in 1..=query.len() {
        let mut cur = vec![0i32; t + 1];
        cur[0] = i as i32;
        for j in 1..=t {
            let cost = if query[i - 1] == target[j - 1] { 0 } else { 1 };
            cur[j] = (prev[j - 1] + cost).min(prev[j] + 1).min(cur[j - 1] + 1);
        }
        prev = cur;
    }

    prev
}

/// Reference edit distance for the given mode.
pub fn ref_distance(query: &[u8], target: &[u8], mode: RefMode) -> i32 {
    let row = dp_last_row(query, target, mode);
    match mode {
        RefMode::Nw => row[target.len()],
        RefMode::Shw | RefMode::Hw => *row.iter().min().unwrap(),
    }
}

/// Reference distance between the full query and `target[start..=end]` (global).
pub fn ref_substring_distance(query: &[u8], target: &[u8], start: usize, end: usize) -> i32 {
    dp_last_row(query, &target[start..=end], RefMode::Nw)[end - start + 1]
}

/// Parses a CIGAR string into `(length, op)` pairs.
pub fn parse_cigar(cigar: &str) -> Vec<(usize, char)> {
    let mut result = Vec::new();
    let mut digits = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            digits.push(ch);
        } else {
            assert!(!digits.is_empty(), "CIGAR op without length: {cigar}");
            result.push((digits.parse::<usize>().unwrap(), ch));
            digits.clear();
        }
    }
    assert!(digits.is_empty(), "CIGAR ends with a number: {cigar}");
    result
}

/// Checks that the CIGAR consumes exactly the query and the reported target span
/// and that the edit operations add up to the reported edit distance.
pub fn check_cigar(
    cigar: &str,
    extended: bool,
    query_len: usize,
    target_span: usize,
    edit_distance: usize,
) {
    let ops = parse_cigar(cigar);
    let mut consumed_query = 0usize;
    let mut consumed_target = 0usize;
    let mut edits = 0usize;

    for (len, op) in ops {
        match op {
            '=' => {
                assert!(extended, "= is only valid in the extended CIGAR format");
                consumed_query += len;
                consumed_target += len;
            }
            'M' => {
                assert!(!extended, "M is only valid in the standard CIGAR format");
                consumed_query += len;
                consumed_target += len;
            }
            'X' => {
                assert!(extended, "X is only valid in the extended CIGAR format");
                consumed_query += len;
                consumed_target += len;
                edits += len;
            }
            'I' => {
                // Insertion to the target: consumes a query character.
                consumed_query += len;
                edits += len;
            }
            'D' => {
                // Deletion from the target: consumes a target character.
                consumed_target += len;
                edits += len;
            }
            other => panic!("unexpected CIGAR operation: {other}"),
        }
    }

    assert_eq!(
        consumed_query, query_len,
        "CIGAR {cigar} does not consume the query"
    );
    assert_eq!(
        consumed_target, target_span,
        "CIGAR {cigar} does not consume the target span"
    );
    if extended {
        assert_eq!(
            edits, edit_distance,
            "CIGAR {cigar} does not explain the edit distance"
        );
    }
}
