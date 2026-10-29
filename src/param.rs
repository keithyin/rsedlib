/// Defines two given characters as equal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EdlibEqualityPair {
    pub first: u8,
    pub second: u8,
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum AlignMode {
    #[default]
    Global,
    Prefix,
    Infix,
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum AlignTask {
    #[default]
    Distance,
    Locations,
    Path,
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum CigarFmt {
    #[default]
    NoCigar,
    Standard,
    Extended,
}

pub struct EdlibAlignParam {
    ///  Set k to non-negative value to tell edlib that edit distance is not larger than k
    ///  Smaller k can significantly improve speed of computation.
    ///  If edit distance is larger than k, edlib will set edit distance to -1.
    ///  Set k to negative value and edlib will internally auto-adjust k until score is found.
    ///  Set to -1 (default) to have no limit on edit distance.
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
    ///  Can be empty if there are none.
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
                first: a,
                second: b,
            })
            .collect();
    }

    pub fn add_eq_pair(&mut self, eq_pair: (u8, u8)) {
        self.additional_eq_pairs.push(EdlibEqualityPair {
            first: eq_pair.0,
            second: eq_pair.1,
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
            k: -1,
            mode: AlignMode::default(),
            task: AlignTask::default(),
            cigar_fmt: CigarFmt::default(),
            additional_eq_pairs: Vec::new(),
        }
    }
}
