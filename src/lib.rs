use std::{ffi::CStr, ops::Deref};

use edlib_sys::{edlibAlignmentToCigar, EdlibAlignConfig, EDLIB_STATUS_OK};
use param::EdlibAlignParam;

pub mod edlib_sys;

pub mod param;
pub mod utils;

#[derive(Debug)]
struct AlignResultGuard(edlib_sys::EdlibAlignResult);
impl AlignResultGuard {
    fn target_start_ends(&self) -> Vec<(usize, usize)> {
        let mut start_ends = Vec::new();
        if self.0.numLocations > 0 {
            for i in 0..self.0.numLocations {
                unsafe {
                    let start = *self.0.startLocations.add(i as usize);
                    let end = *self.0.endLocations.add(i as usize);
                    start_ends.push((start as usize, end as usize));
                }
            }
        }

        start_ends
    }
}

impl Deref for AlignResultGuard {
    type Target = edlib_sys::EdlibAlignResult;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl From<edlib_sys::EdlibAlignResult> for AlignResultGuard {
    fn from(result: edlib_sys::EdlibAlignResult) -> Self {
        AlignResultGuard(result)
    }
}

struct AlignCigarGuard(*mut i8);

impl Deref for AlignCigarGuard {
    type Target = *mut i8;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl Drop for AlignCigarGuard {
    fn drop(&mut self) {
        unsafe {
            libc::free(self.0 as *mut std::ffi::c_void);
        }
    }
}

impl From<*mut i8> for AlignCigarGuard {
    fn from(cigar_str: *mut i8) -> Self {
        AlignCigarGuard(cigar_str)
    }
}

impl Drop for AlignResultGuard {
    fn drop(&mut self) {
        unsafe {
            edlib_sys::edlibFreeAlignResult(self.0);
        }
    }
}

#[derive(Debug)]
pub struct EdlibAlignResult {
    pub edit_distance: i32,
    pub alphabet_length: i32,
    pub locations: Vec<(usize, usize)>,
    pub cigar: Option<String>,
}

pub fn edlib_align(
    query: &[u8],
    target: &[u8],
    aln_param: &EdlibAlignParam,
) -> Result<EdlibAlignResult, String> {
    let config = EdlibAlignConfig {
        k: aln_param.k(),
        mode: match aln_param.mode() {
            param::AlignMode::Global => edlib_sys::EdlibAlignMode_EDLIB_MODE_NW,
            param::AlignMode::Prefix => edlib_sys::EdlibAlignMode_EDLIB_MODE_SHW,
            param::AlignMode::Infix => edlib_sys::EdlibAlignMode_EDLIB_MODE_HW,
        },
        task: match aln_param.task() {
            param::AlignTask::Distance => edlib_sys::EdlibAlignTask_EDLIB_TASK_DISTANCE,
            param::AlignTask::Locations => edlib_sys::EdlibAlignTask_EDLIB_TASK_LOC,
            param::AlignTask::Path => edlib_sys::EdlibAlignTask_EDLIB_TASK_PATH,
        },
        additionalEqualities: if aln_param.additional_eq_pairs().len() > 0 {
            aln_param.additional_eq_pairs().as_ptr()
        } else {
            std::ptr::null()
        },
        additionalEqualitiesLength: aln_param.additional_eq_pairs().len() as i32,
    };

    let edlib_raw_res: AlignResultGuard = unsafe {
        edlib_sys::edlibAlign(
            query.as_ptr() as *const i8,
            query.len() as i32,
            target.as_ptr() as *const i8,
            target.len() as i32,
            config,
        )
        .into()
    };

    // println!("{:?}", edlib_raw_res);

    let align_cigar_str = if aln_param.task() == param::AlignTask::Path
        && aln_param.cigar_fmt() != param::CigarFmt::NoCigar
        && edlib_raw_res.status == EDLIB_STATUS_OK as i32
    {
        let cigar_fmt = match aln_param.cigar_fmt() {
            param::CigarFmt::Standard => edlib_sys::EdlibCigarFormat_EDLIB_CIGAR_STANDARD,
            param::CigarFmt::Extended => edlib_sys::EdlibCigarFormat_EDLIB_CIGAR_EXTENDED,
            _ => panic!("Invalid cigar format"),
        };

        unsafe {
            let cigar_str_guard: AlignCigarGuard = edlibAlignmentToCigar(
                edlib_raw_res.alignment,
                edlib_raw_res.alignmentLength,
                cigar_fmt,
            )
            .into();
            Some(
                CStr::from_ptr(*cigar_str_guard)
                    .to_str()
                    .unwrap()
                    .to_string(),
            )
        }
    } else {
        None
    };

    if edlib_raw_res.status != EDLIB_STATUS_OK as i32 {
        Err(format!(
            "Edlib alignment failed with status: {}",
            edlib_raw_res.status
        ))
    } else {
        Ok(EdlibAlignResult {
            edit_distance: edlib_raw_res.editDistance,
            alphabet_length: edlib_raw_res.alphabetLength,
            locations: edlib_raw_res.target_start_ends(),
            cigar: align_cigar_str,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::reverse_complement;

    #[test]
    fn test_edlib_align() {
        let query = b"AAAGAGAGGAGATGCATCGAGAGTCAGAGCAGAGAGGGTCGATTATGGTCTCAAGAATGACGCTATGAAGTATGCATACACATGGCACAGTTGAACAAATATATTTTAATATAAGTAGCGAACGGGCGTAATATACAAGGAATTTTACCTTTAATTTGGATAACTGAGTAATTTCGGCTAATATACCAAATGACCTTGTGTTAAACGCAAGAAAGATTGTATCGTTTAAAGTAAGCTTGTACAGATTAAGTTGGAAGTAATAGCTTACCAGTCTATGATCCTTGAGCTAGTCTTGGGGATGGACGGATAGCCACACGAAGACTGAAGACAAGGTTCAGATTCTATGGAAGGAAGGCAGCAGTGTGATTTGGACAATGAGCGAAAGCTTGATCCACAAATACTTCGTGTGTGATCAGAAGAATAGGAGACTTTTGAAAGCACTAACAGTAAAGACGAAATTCGCTTATTTACAGAAGAAGCTCCGCAAATTTCCTGCAAGCAGCCGCGGTAATACGAAAGGAGCAGGTGTTATCAGATTAACTGGCGGTAAAGGGCACGTAGATGGTTTATAGTGTATACTATGAAATTATAAAATAAAGCATAACTTTTATAATGTGGTATAACACATATTAGACTTGAGTTTATAGGAGGGGTAGTAGAATTTTAATGTAAAGGTAAAATTTTTGATAAATAGAGAGATACCAAGGCGTAGGCAAACATCTAGACCTTAACTGACGATTGAGGTGCAAGGTATGGGAGCAATAGGATAGATACCCTAGTAGTCATACAGTAACGATGAATATTAAATTGTTGGGAATAAATCCAACAGTATACAGCTAACCGTAAAATATTCCCGCCTGGGAGTACAACCCAAGGTTGAAACTTAAAGGAATTGACGGGATCTAAACAAGCGGTGGAACATGTGGGTTTAATCCATACTACTGCGGTAAATCTACAGTTTTTGTATTATATAAGGTGTTGCAGGCTTCCGTCAGTTCGTGTTGTGAAAACGTGGTTACCCCTTTTTAACGGAACGCAACCCCCTAATTTTGTTACTAAGTCAAACCTTTAAAATAGTTTGTCAAGAACTCAAAAAATGCTTATGATCAATAAAGATGAAAGGTGGGGATGACGTCAAGTCGTCAGGCCCTTTATAAACTGGGCACACACGTGTTACAATAGTTACTACAATAAGAAGCAACGAACGTGAGTTAGAGGCAAATCTTGAAAAGTAATTTAAGTTCGAATTGTTCTTGTAACTCGAGAACATGAAGTTGAAATCGCTAGTAATCACAATTAGTAACTTTCGTTGTGGTGAATACGTTCTTAGATCTGTACACACCGCCCGTCACACCCGCGAAATTGGTTTTACTTTAAGTAGATTTAAAGAACCAAAGTTGGAGGCTTTTGAGAAAACAAGGACTTTTAGTACCAAAGTAGAATTAGTGACTGTGGTGAAGTCGTAACAAGGTAACTCTGTGTGCTCGCATGGATCATCTCTCT";

        let query_rc = &reverse_complement(query);

        let target = b"GAGACAGATGCATCCATAGCGACACACAGAGTTACCATTGTTACGACTTCACCACAGTCACTAATTCTACTTTGGTACTAAAAGTCCTAGTTTTGCTCAAAGCCTCCACTGGTTCTTTAAATCTTTTTTTAGCTACTTAAAGTAAACCAATTTCCCAGGGTGTGACGGGCGGTGTGTACAAGACTAATGAACGTATTCACCACAACGTAAACTAAATTGTGGATACTAGCGATTTCAACTTCATGTTCTCGAACTTACAAAACAATTCGAACTTAAATTACTTTTCAAGATTTGCTCAACTCACGGTTGTTGCTTCTTATTGTAGTAACTATGGTAACACGTGTGTATGCCCAGTTTATAAGGGCCATGACGACTTTGACGTCATCCCACTTCATCTTAATATCATAAGCAATTTTTTGAGTGCTTGACAAACTATTTTAAAGTTTGACTTAGTAACAAAATAGGGATGCGTTCTAAAAGGACTGAACCAAACGTTTCACAACACGAACTGACGACAGCCATGCAACACCTGTATATTATACAAAAACTGGTAAGATTTACGCGTAGTATCGGATTAAACCGCACAATGTTCCCACCGCTTGTATAGATCCCCGTCAATTCCTTTAAGTTTCAACCTTGCGGGTTGTACCTCCCCAGGCGGATATTGTAACGCGTTAGCTGTAATACTGAGATTTATCCCAACATTTAATATTCATCGTTTACTGTATGGACTACTAGGGTATCCTAATCCTATTTGTCCACCTTCGCACCTCCAATGTCATTAAGTCTAGTAGTTGCTACGCATTGGTAAATTCTCTCTTAATTTTATCAAATTTTACTTTACATTAGAAATTCTACACCCACCTACCTAACCCAAGCTAATATGAGTTTATACCAACATTATAATTGTTATGCTTTATAATTTCATAGTATACACTATAAACCATCTACGTGCCCTTTACCGCCCAGTAATCTGAATAACACCTGCTCCTTTCGTATTACCGCGGCTGCTTGCACGAAATTTGCCGGAGCTTCTTCTGTAAATAAAGTCAATTTCGTCTTACTGATAGTCTTTACAAATAGATCTCCCTATTCTTCTGATCACACACGGAAGTATTGCTGGATCAAGCTTTTCGCCTCATGTCCAAATTCCACACTTGCTGCTTCCAAAAGGTCTGAACCTTGTCTAGTTTCAGTGTGGCTATCCCGTCATCCCAACTAGCTAAGACATAGACTTGGTAAGCTATTACCTTACACAAACTATCTATCCTGTACAAAGCTTATCTTTTAACGATACAATTCCTTTCTTGCTTTAACACAAGTCATTTGGTATTAGCCGAAATACTCCAGTTTATCCCATATTAAAAGGTAAATTCTTGTATATTACGCACCCGTTCGCTACTTAATTAAATATATTTGTTCAACTTGCATGTGTTATGCATCTCCTAGCGTTCATTCTGAGCCATAATCGAACCCTCTCTGCTTGACTTCGATGCATCTCTCTC";

        let mut param = EdlibAlignParam::default();
        param.set_k((query.len() / 2) as i32);
        // param.set_mode(param::AlignMode::Global);
        param.set_task(param::AlignTask::Path);
        param.set_cigar_fmt(param::CigarFmt::Extended);
        let aln_res = edlib_align(query, target, &param);
        println!("{:?}", aln_res);

        let aln_res = edlib_align(query_rc, target, &param);
        println!("{:?}", aln_res);
    }

    #[test]
    fn test_edlib_align2() {
        let query = b"elephant";
        let target = b"telephone";

        let mut param = EdlibAlignParam::default();
        // param.set_mode(param::AlignMode::Global);
        param.set_task(param::AlignTask::Path);
        param.set_cigar_fmt(param::CigarFmt::Extended);
        let aln_res = edlib_align(query, target, &param);
        println!("{:?}", aln_res);
    }
}
