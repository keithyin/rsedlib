//! Pure Rust implementation of the [edlib](https://github.com/Martinsos/edlib)
//! edit distance algorithms.
//!
//! The public API mirrors the previous C/C++ binding, so existing callers keep
//! working unchanged, but no C++ toolchain, `cmake` or `libc` is needed anymore.
//!
//! ```no_run
//! use rsedlib::{edlib_align, param::{AlignTask, CigarFmt, EdlibAlignParam}};
//!
//! let mut param = EdlibAlignParam::default();
//! param.set_task(AlignTask::Path);
//! param.set_cigar_fmt(CigarFmt::Extended);
//! let res = edlib_align(b"elephant", b"telephone", &param).unwrap();
//! assert_eq!(res.edit_distance, 3);
//! assert_eq!(res.cigar.as_deref(), Some("1D5=1X1=1X"));
//! ```

#![deny(unsafe_code)]

use param::{AlignTask, CigarFmt, EdlibAlignParam};

mod alignment;
mod alphabet;
mod cigar;
mod myers;
mod search;
mod traceback;

pub mod param;
pub mod utils;

#[derive(Debug, Clone, PartialEq)]
pub struct EdlibAlignResult {
    pub edit_distance: i32,
    pub alphabet_length: i32,
    pub locations: Vec<(usize, usize)>,
    pub cigar: Option<String>,
}

/// Aligns `query` against `target` following edlib's semantics.
///
/// `Err` is returned for invalid input only (e.g. sequences longer than
/// `i32::MAX`). Note that, exactly like the C library, an edit distance of `-1`
/// is reported when it is larger than the configured `k` and therefore no
/// alignment could be found.
pub fn edlib_align(
    query: &[u8],
    target: &[u8],
    aln_param: &EdlibAlignParam,
) -> Result<EdlibAlignResult, String> {
    let outcome = alignment::align(query, target, aln_param)?;

    let mut locations = Vec::with_capacity(outcome.end_locations.len());
    for (i, &end) in outcome.end_locations.iter().enumerate() {
        // When start locations are not calculated the C wrapper reported 0.
        let start = outcome.start_locations.get(i).copied().unwrap_or(0) as usize;
        locations.push((start, end as usize));
    }

    let cigar = if aln_param.task() == AlignTask::Path && aln_param.cigar_fmt() != CigarFmt::NoCigar
    {
        let extended = aln_param.cigar_fmt() == CigarFmt::Extended;
        let alignment = outcome.alignment.as_deref().unwrap_or(&[]);
        cigar::alignment_to_cigar(alignment, extended)
            .map(Some)
            .ok_or_else(|| "Invalid alignment path".to_string())?
    } else {
        None
    };

    Ok(EdlibAlignResult {
        edit_distance: outcome.edit_distance,
        alphabet_length: outcome.alphabet_length,
        locations,
        cigar,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::reverse_complement;

    #[test]
    fn test_edlib_align() {
        let query = b"ATAATTTGTTCTCTGCGTTGACGTGCTTGAACCCTTAGCATCTTCTAAAATCTTTTTGACTTCATTCTTGTGTTTCCTTTAAGTTTTGTTATTTTCTTCTTCAAGTTCTGTGCAATTTAATAGCTTGTCTGCGTCATCGAATATCTCTGTTAATATCTCTTTCACGTTTATCCCATTACATCTTTCAGGCACCCCAATGCGAACTTTTTAAGTAACGCAAGAATAGCGATGAAATTAGGGACCC";

        let query_rc = &reverse_complement(query);

        let target = b"GGTCCTAACTTTCATCGAATTACTTGGCGTTACTTAAAAGTCGCATGGGTCCATTGAAAGATGTAATGGATAACGCTGAAAGAGATAATTAACAAGATATCGATGACGCAGAACAAGCTAAGTTAAATGGCACAGAAACCTTGAAGAAGAAATAAACCAAAAACTTAAGAAAGCCAAGAAAAGTTCAAAAAGATTTTTAAGAAGATGCTAAGGTTCAGCAGTCAACAGCAAGAAACAAATTTTAT";

        let param = EdlibAlignParam::default();
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
        let aln_res = edlib_align(query, target, &param).unwrap();
        println!("{:?}", aln_res);
        assert_eq!(aln_res.edit_distance, 3);
        assert_eq!(aln_res.alphabet_length, 8);
        assert_eq!(aln_res.locations, vec![(0, 8)]);
        assert_eq!(aln_res.cigar.as_deref(), Some("1D5=1X1=1X"));
    }
}
