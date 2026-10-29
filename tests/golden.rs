//! Golden values taken from the documented behaviour of the C implementation.

mod support;

use rsedlib::edlib_align;
use rsedlib::param::{AlignMode, AlignTask, CigarFmt, EdlibAlignParam};
use support::{check_cigar, ref_distance, RefMode};

/// Example from the README / the original crate docs.
#[test]
fn readme_example() {
    let query = b"elephant";
    let target = b"telephone";

    let mut param = EdlibAlignParam::default();
    param.set_task(AlignTask::Path);
    param.set_cigar_fmt(CigarFmt::Extended);
    let result = edlib_align(query, target, &param).unwrap();

    assert_eq!(
        result,
        rsedlib::EdlibAlignResult {
            edit_distance: 3,
            alphabet_length: 8,
            locations: vec![(0, 8)],
            cigar: Some("1D5=1X1=1X".to_string()),
        }
    );
}

#[test]
fn standard_cigar_of_readme_example() {
    let mut param = EdlibAlignParam::default();
    param.set_task(AlignTask::Path);
    param.set_cigar_fmt(CigarFmt::Standard);
    let result = edlib_align(b"elephant", b"telephone", &param).unwrap();

    assert_eq!(result.cigar.as_deref(), Some("1D8M"));
}

#[test]
fn no_cigar_without_path_task() {
    for task in [AlignTask::Distance, AlignTask::Locations] {
        let mut param = EdlibAlignParam::default();
        param.set_task(task);
        param.set_cigar_fmt(CigarFmt::Extended);
        let result = edlib_align(b"elephant", b"telephone", &param).unwrap();
        assert!(result.cigar.is_none());
        assert_eq!(result.edit_distance, 3);
    }
}

/// The DNA example used by the crate's own unit tests.
#[test]
fn dna_example() {
    let query = b"ATAATTTGTTCTCTGCGTTGACGTGCTTGAACCCTTAGCATCTTCTAAAATCTTTTTGACTTCATTCTTGTGTTTCCTTTAAGTTTTGTTATTTTCTTCTTCAAGTTCTGTGCAATTTAATAGCTTGTCTGCGTCATCGAATATCTCTGTTAATATCTCTTTCACGTTTATCCCATTACATCTTTCAGGCACCCCAATGCGAACTTTTTAAGTAACGCAAGAATAGCGATGAAATTAGGGACCC";
    let target = b"GGTCCTAACTTTCATCGAATTACTTGGCGTTACTTAAAAGTCGCATGGGTCCATTGAAAGATGTAATGGATAACGCTGAAAGAGATAATTAACAAGATATCGATGACGCAGAACAAGCTAAGTTAAATGGCACAGAAACCTTGAAGAAGAAATAAACCAAAAACTTAAGAAAGCCAAGAAAAGTTCAAAAAGATTTTTAAGAAGATGCTAAGGTTCAGCAGTCAACAGCAAGAAACAAATTTTAT";

    for (mode, ref_mode) in [
        (AlignMode::Global, RefMode::Nw),
        (AlignMode::Prefix, RefMode::Shw),
        (AlignMode::Infix, RefMode::Hw),
    ] {
        for (task, fmt) in [
            (AlignTask::Distance, CigarFmt::NoCigar),
            (AlignTask::Locations, CigarFmt::NoCigar),
            (AlignTask::Path, CigarFmt::Extended),
        ] {
            let mut param = EdlibAlignParam::default();
            param.set_mode(mode);
            param.set_task(task);
            param.set_cigar_fmt(fmt);
            let result = edlib_align(query, target, &param).unwrap();

            let expected = ref_distance(query, target, ref_mode);
            assert_eq!(result.edit_distance, expected, "mode {mode:?}");
            assert!(!result.locations.is_empty());

            if task == AlignTask::Path {
                let cigar = result.cigar.clone().unwrap();
                let (start, end) = result.locations[0];
                check_cigar(
                    &cigar,
                    true,
                    query.len(),
                    end - start + 1,
                    result.edit_distance as usize,
                );
            }
        }
    }
}

/// The query must fit into the target at the position reported by HW mode.
#[test]
fn infix_finds_substring() {
    let target = b"GGGGACGTACGTACGTGGGG";
    let query = b"ACGTACGTACGT";

    let mut param = EdlibAlignParam::default();
    param.set_mode(AlignMode::Infix);
    param.set_task(AlignTask::Locations);
    let result = edlib_align(query, target, &param).unwrap();

    assert_eq!(result.edit_distance, 0);
    assert_eq!(result.locations, vec![(4, 15)]);

    // A mutated copy still maps to the same span.
    let mutated = b"ACGTTCGTACGT";
    let result = edlib_align(mutated, target, &param).unwrap();
    assert_eq!(result.edit_distance, 1);
    assert_eq!(result.locations, vec![(4, 15)]);
}
