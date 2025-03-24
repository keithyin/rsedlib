use std::ffi::c_char;

use crate::edlib_sys::EdlibEqualityPair;

#[derive(Debug, Clone, Copy)]
pub enum AlignMode {
    Global,
    Prefix,
    Infix,
}

#[derive(Debug, Clone, Copy)]
pub enum AlignTask {
    Distance,
    Locations,
    Path,
}

pub struct EdlibAlignParam {
    k: i32,
    mode: AlignMode,
    task: AlignTask,
    additional_eq_pairs: Vec<EdlibEqualityPair>,
}

impl EdlibAlignParam {
    pub fn new(k: i32, mode: AlignMode, task: AlignTask) -> Self {
        EdlibAlignParam {
            k,
            mode,
            task,
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

    pub fn k(&self) -> i32 {
        self.k
    }

    pub fn mode(&self) -> AlignMode {
        self.mode
    }

    pub fn task(&self) -> AlignTask {
        self.task
    }

    pub fn additional_eq_pairs(&self) -> &[EdlibEqualityPair] {
        &self.additional_eq_pairs
    }
}

impl Default for EdlibAlignParam {
    fn default() -> Self {
        EdlibAlignParam {
            k: 0,
            mode: AlignMode::Global,
            task: AlignTask::Distance,
            additional_eq_pairs: Vec::new(),
        }
    }
}
