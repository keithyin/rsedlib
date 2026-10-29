//! Combination matrix of modes, tasks and k values, plus edge cases.

mod support;

use rsedlib::edlib_align;
use rsedlib::param::{AlignMode, AlignTask, CigarFmt, EdlibAlignParam};
use support::{check_cigar, ref_distance, RefMode};

const MODES: [(AlignMode, RefMode); 3] = [
    (AlignMode::Global, RefMode::Nw),
    (AlignMode::Prefix, RefMode::Shw),
    (AlignMode::Infix, RefMode::Hw),
];

const TASKS: [AlignTask; 3] = [AlignTask::Distance, AlignTask::Locations, AlignTask::Path];

#[test]
fn mode_task_k_matrix() {
    let query = b"elephant";
    let target = b"telephone";

    for (align_mode, ref_mode) in MODES {
        let expected = ref_distance(query, target, ref_mode);
        assert!(expected > 0);

        let mut k_values = vec![-1, 0, expected - 1, expected, expected + 5];
        k_values.sort_unstable();
        k_values.dedup();

        for task in TASKS {
            for &k in &k_values {
                let mut param = EdlibAlignParam::default();
                param.set_mode(align_mode);
                param.set_task(task);
                param.set_cigar_fmt(CigarFmt::Extended);
                param.set_k(k);

                let result = edlib_align(query, target, &param).unwrap();
                let label = format!("mode {align_mode:?}, task {task:?}, k {k}");

                if k >= 0 && k < expected {
                    assert_eq!(result.edit_distance, -1, "{label}");
                    assert!(result.locations.is_empty(), "{label}");
                    if task == AlignTask::Path {
                        // No alignment was computed; the C binding then produced
                        // an empty CIGAR string.
                        assert_eq!(result.cigar.as_deref(), Some(""), "{label}");
                    } else {
                        assert!(result.cigar.is_none(), "{label}");
                    }
                    continue;
                }

                assert_eq!(result.edit_distance, expected, "{label}");
                assert!(!result.locations.is_empty(), "{label}");

                if task == AlignTask::Path {
                    let cigar = result.cigar.clone().expect("cigar expected");
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
}

#[test]
fn distance_task_reports_start_zero() {
    let query = b"ACGTACGT";
    let target = b"TTACGTACGTGG";

    for (align_mode, _) in MODES {
        let mut param = EdlibAlignParam::default();
        param.set_mode(align_mode);
        param.set_task(AlignTask::Distance);
        let result = edlib_align(query, target, &param).unwrap();

        assert!(result.edit_distance >= 0);
        assert!(!result.locations.is_empty());
        // Start locations are not calculated for the distance task, the binding
        // reported 0 for them.
        assert!(result.locations.iter().all(|&(start, _)| start == 0));
        assert!(result.cigar.is_none());
    }
}

#[test]
fn identical_sequences() {
    let seq = b"ACGTACGTACGT";

    for task in TASKS {
        let mut param = EdlibAlignParam::default();
        param.set_task(task);
        param.set_cigar_fmt(CigarFmt::Extended);
        param.set_k(0);
        let result = edlib_align(seq, seq, &param).unwrap();

        assert_eq!(result.edit_distance, 0);
        assert_eq!(result.locations, vec![(0, seq.len() - 1)]);
        if task == AlignTask::Path {
            assert_eq!(result.cigar.as_deref(), Some("12="));
        }
    }
}

#[test]
fn empty_sequences() {
    let target = b"ACGT";

    // Global mode: all target characters have to be deleted.
    let mut param = EdlibAlignParam::default();
    param.set_mode(AlignMode::Global);
    param.set_task(AlignTask::Locations);
    let result = edlib_align(b"", target, &param).unwrap();
    assert_eq!(result.edit_distance, 4);
    assert_eq!(result.locations, vec![(0, 3)]);
    assert_eq!(result.alphabet_length, 4);

    // Prefix / infix: an empty query matches anywhere for free.
    for (align_mode, _) in MODES {
        if align_mode == AlignMode::Global {
            continue;
        }
        let mut param = EdlibAlignParam::default();
        param.set_mode(align_mode);
        param.set_task(AlignTask::Locations);
        let result = edlib_align(b"", target, &param).unwrap();
        assert_eq!(result.edit_distance, 0);
        // -1 end location is exposed as a wrapped usize by the binding.
        assert_eq!(result.locations, vec![(0, usize::MAX)]);
    }

    // Empty target with a non-empty query.
    for (align_mode, ref_mode) in MODES {
        let mut param = EdlibAlignParam::default();
        param.set_mode(align_mode);
        param.set_task(AlignTask::Path);
        param.set_cigar_fmt(CigarFmt::Extended);
        let result = edlib_align(target, b"", &param).unwrap();
        assert_eq!(result.edit_distance, ref_distance(target, b"", ref_mode));
        // edlib returns before calculating an alignment, so the binding reported
        // an empty CIGAR.
        assert_eq!(result.cigar.as_deref(), Some(""));
    }
}

#[test]
fn additional_equalities() {
    // Case insensitive comparison.
    let mut param = EdlibAlignParam::default();
    param.set_eq_pairs(vec![(b'A', b'a'), (b'C', b'c'), (b'G', b'g'), (b'T', b't')]);
    param.set_task(AlignTask::Path);
    param.set_cigar_fmt(CigarFmt::Extended);
    let result = edlib_align(b"ACGT", b"acgt", &param).unwrap();
    assert_eq!(result.edit_distance, 0);
    assert_eq!(result.alphabet_length, 8);
    assert_eq!(result.locations, vec![(0, 3)]);
    assert_eq!(result.cigar.as_deref(), Some("4="));

    // N matches G (wildcard).
    let mut param = EdlibAlignParam::default();
    param.add_eq_pair((b'N', b'G'));
    let result = edlib_align(b"ACNT", b"ACGT", &param).unwrap();
    assert_eq!(result.edit_distance, 0);

    // Equalities between characters that are not part of the alphabet are ignored.
    let mut param = EdlibAlignParam::default();
    param.add_eq_pair((b'X', b'A'));
    let result = edlib_align(b"ACGT", b"ACGT", &param).unwrap();
    assert_eq!(result.edit_distance, 0);
    assert_eq!(result.alphabet_length, 4);
}

#[test]
fn large_alignment_uses_hirschberg() {
    let mut rng = support::Lcg::new(1234);
    let alphabet = b"ACGT";
    let query = rng.sequence(3000, alphabet);
    // Target with roughly 5% substitutions + some indels.
    let mut target = Vec::with_capacity(query.len() + 60);
    for (i, &base) in query.iter().enumerate() {
        if i % 97 == 0 && target.len() > 4 {
            target.pop();
            continue;
        }
        if i % 131 == 0 {
            target.push(b'A');
        }
        target.push(if i % 173 == 0 {
            alphabet[rng.next_below(alphabet.len())]
        } else {
            base
        });
    }

    let mut param = EdlibAlignParam::default();
    param.set_task(AlignTask::Path);
    param.set_cigar_fmt(CigarFmt::Extended);
    let result = edlib_align(&query, &target, &param).unwrap();

    assert!(result.edit_distance > 0);
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
