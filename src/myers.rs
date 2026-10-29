//! Bit-parallel Myers primitives.
//!
//! Faithful port of the low level helpers from edlib's `edlib.cpp`
//! (`calculateBlock`, `readBlock`, `readBlockReverse`, `getBlockCellValues`,
//! `allBlockCellsLarger`). All arithmetic that overflows in C is expressed with
//! wrapping operations so behaviour matches the reference implementation.

/// Number of bits in a `Word`.
pub(crate) const WORD_SIZE: usize = 64;

/// `100..00`
pub(crate) const HIGH_BIT_MASK: u64 = 1u64 << (WORD_SIZE - 1);

/// `00..01`
pub(crate) const WORD_1: u64 = 1;

/// Edit operation codes used by edlib's alignment representation.
pub(crate) const EDOP_MATCH: u8 = 0;
pub(crate) const EDOP_INSERT: u8 = 1;
pub(crate) const EDOP_DELETE: u8 = 2;
pub(crate) const EDOP_MISMATCH: u8 = 3;

/// Does ceiling division `x / y` for non-negative values.
#[inline]
pub(crate) fn ceil_div(x: usize, y: usize) -> usize {
    x.div_ceil(y)
}

/// Same as [`ceil_div`] but for `i32` with C style truncated division.
#[inline]
pub(crate) fn ceil_div_i32(x: i32, y: i32) -> i32 {
    if x % y != 0 {
        x / y + 1
    } else {
        x / y
    }
}

/// One 64-cell block (part of a DP column).
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Block {
    /// `Pv` of the last cell of the block.
    pub(crate) p: u64,
    /// `Mv` of the last cell of the block.
    pub(crate) m: u64,
    /// Score of the last cell of the block.
    pub(crate) score: i32,
}

impl Block {
    #[inline]
    pub(crate) fn new(p: u64, m: u64, score: i32) -> Self {
        Block { p, m, score }
    }
}

/// Advances one block of one column (Myers' `Advance_Block`).
///
/// The highest bit of the word is the most bottom cell of the block.
/// `hin` is one of `1`, `0`, `-1`.
#[inline]
pub(crate) fn calculate_block(pv: u64, mv: u64, eq: u64, hin: i32) -> (u64, u64, i32) {
    // 1 if hin == -1, 0 otherwise.
    let hin_is_neg = ((hin >> 2) as u64) & WORD_1;

    let xv = eq | mv;
    let eq = eq | hin_is_neg;
    let xh = ((eq & pv).wrapping_add(pv) ^ pv) | eq;

    let ph = mv | !(xh | pv);
    let mh = pv & xh;

    let mut hout: i32 = ((ph & HIGH_BIT_MASK) >> (WORD_SIZE - 1)) as i32;
    hout -= ((mh & HIGH_BIT_MASK) >> (WORD_SIZE - 1)) as i32;

    let mut ph = ph << 1;
    let mut mh = mh << 1;

    mh |= hin_is_neg;
    ph |= ((hin + 1) >> 1) as u64;

    let pv_out = mh | !(xv | ph);
    let mv_out = ph & xv;

    (pv_out, mv_out, hout)
}

/// Values of the cells in the block, starting with the bottom cell.
#[inline]
pub(crate) fn block_cell_values(block: &Block) -> [i32; WORD_SIZE] {
    let mut scores = [0i32; WORD_SIZE];
    let mut score = block.score;
    let mut mask = HIGH_BIT_MASK;
    for slot in scores.iter_mut().take(WORD_SIZE - 1) {
        *slot = score;
        if block.p & mask != 0 {
            score -= 1;
        }
        if block.m & mask != 0 {
            score += 1;
        }
        mask >>= 1;
    }
    scores[WORD_SIZE - 1] = score;
    scores
}

/// Writes values of the cells in the block into `dest`, starting with the top cell.
///
/// `dest` must have at least `WORD_SIZE` elements.
#[inline]
pub(crate) fn read_block(block: &Block, dest: &mut [i32]) {
    debug_assert!(dest.len() >= WORD_SIZE);
    let mut score = block.score;
    let mut mask = HIGH_BIT_MASK;
    for i in 0..WORD_SIZE - 1 {
        dest[WORD_SIZE - 1 - i] = score;
        if block.p & mask != 0 {
            score -= 1;
        }
        if block.m & mask != 0 {
            score += 1;
        }
        mask >>= 1;
    }
    dest[0] = score;
}

/// Writes values of the cells in the block into `dest`, starting with the bottom cell.
///
/// `dest` must have at least `WORD_SIZE` elements.
#[inline]
pub(crate) fn read_block_reverse(block: &Block, dest: &mut [i32]) {
    debug_assert!(dest.len() >= WORD_SIZE);
    let mut score = block.score;
    let mut mask = HIGH_BIT_MASK;
    for slot in dest.iter_mut().take(WORD_SIZE - 1) {
        *slot = score;
        if block.p & mask != 0 {
            score -= 1;
        }
        if block.m & mask != 0 {
            score += 1;
        }
        mask >>= 1;
    }
    dest[WORD_SIZE - 1] = score;
}

/// True if all cells in the block have value larger than `k`.
pub(crate) fn all_block_cells_larger(block: &Block, k: i32) -> bool {
    let scores = block_cell_values(block);
    scores.iter().all(|&score| score > k)
}
