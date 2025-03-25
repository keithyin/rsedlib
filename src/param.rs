use std::ffi::c_char;

use crate::edlib_sys::EdlibEqualityPair;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AlignMode {
    Global,
    Prefix,
    Infix,
}

impl Default for AlignMode {
    fn default() -> Self {
        AlignMode::Global
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AlignTask {
    Distance,
    Locations,
    Path,
}

impl Default for AlignTask {
    fn default() -> Self {
        AlignTask::Distance
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CigarFmt {
    NoCigar,
    Standard,
    Extended,
}

impl Default for CigarFmt {
    fn default() -> Self {
        CigarFmt::NoCigar
    }
}

pub struct EdlibAlignParam {
    /// Set k to non-negative value to tell edlib that edit distance is not larger than k
    ///  Smaller k can significantly improve speed of computation.
    ///  If edit distance is larger than k, edlib will set edit distance to -1.
    ///  Set k to negative value and edlib will internally auto-adjust k until score is found.
    k: i32,

    mode: AlignMode,

    ///  Alignment task - tells Edlib what to calculate. Less to calculate, faster it is.
    ///  AlignTask::Distance - find edit distance and end locations of optimal alignment paths in target.
    ///  AlignTask::Locations - find edit distance and start and end locations of optimal alignment paths in target.
    ///  AlignTask::Path - find edit distance, alignment path (and start and end locations of it in target).
    task: AlignTask,
    cigar_fmt: CigarFmt,

    ///  List of pairs of characters, where each pair defines two characters as equal.
    ///  This way you can extend edlib's definition of equality (which is that each character is equal only to itself).
    ///  This can be useful if you have some wildcard characters that should match multiple other characters,
    ///  or e.g. if you want edlib to be case insensitive.
    ///  Can be set to NULL if there are none.
    additional_eq_pairs: Vec<EdlibEqualityPair>,
}

impl EdlibAlignParam {
    pub fn new(k: i32, mode: AlignMode, task: AlignTask) -> Self {
        EdlibAlignParam {
            k,
            mode,
            task,
            cigar_fmt: CigarFmt::NoCigar,
            additional_eq_pairs: Vec::new(),
        }
    }
    pub fn set_eq_pairs(&mut self, eq_pairs: Vec<(u8, u8)>) {
        self.additional_eq_pairs = eq_pairs
            .into_iter()
            .map(|(a, b)| EdlibEqualityPair {
                first: a as c_char,
                second: b as c_char,
            })
            .collect();
    }

    pub fn add_eq_pair(&mut self, eq_pair: (u8, u8)) {
        self.additional_eq_pairs.push(EdlibEqualityPair {
            first: eq_pair.0 as c_char,
            second: eq_pair.1 as c_char,
        });
    }

    pub fn set_k(&mut self, k: i32) {
        self.k = k;
    }

    pub fn k(&self) -> i32 {
        self.k
    }

    pub fn set_mode(&mut self, mode: AlignMode) {
        self.mode = mode;
    }

    pub fn mode(&self) -> AlignMode {
        self.mode
    }

    pub fn set_task(&mut self, task: AlignTask) {
        self.task = task;
    }

    pub fn task(&self) -> AlignTask {
        self.task
    }

    pub fn set_cigar_fmt(&mut self, cigar_fmt: CigarFmt) {
        self.cigar_fmt = cigar_fmt;
    }

    pub fn cigar_fmt(&self) -> CigarFmt {
        self.cigar_fmt
    }

    pub fn additional_eq_pairs(&self) -> &[EdlibEqualityPair] {
        &self.additional_eq_pairs
    }
}

impl Default for EdlibAlignParam {
    fn default() -> Self {
        EdlibAlignParam {
            k: 0,
            mode: AlignMode::default(),
            task: AlignTask::default(),
            cigar_fmt: CigarFmt::default(),
            additional_eq_pairs: Vec::new(),
        }
    }
}
