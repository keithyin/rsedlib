//! Randomised comparison against a naive O(n*m) dynamic programming reference.

mod support;

use rsedlib::edlib_align;
use rsedlib::param::{AlignMode, AlignTask, CigarFmt, EdlibAlignParam};
use support::{check_cigar, ref_distance, ref_substring_distance, Lcg, RefMode};

fn mode_of(mode: RefMode) -> AlignMode {
    match mode {
        RefMode::Nw => AlignMode::Global,
        RefMode::Shw => AlignMode::Prefix,
        RefMode::Hw => AlignMode::Infix,
    }
}

#[test]
fn distance_matches_naive_dp() {
    let mut rng = Lcg::new(42);
    let alphabets: [&[u8]; 2] = [b"ACGT", b"AB"];

    for (case, alphabet) in alphabets.iter().enumerate() {
        for _ in 0..300 {
            let query_len = rng.next_below(40);
            let target_len = rng.next_below(40);
            let query = rng.sequence(query_len, alphabet);
            let target = rng.sequence(target_len, alphabet);

            for mode in [RefMode::Nw, RefMode::Shw, RefMode::Hw] {
                let mut param = EdlibAlignParam::default();
                param.set_mode(mode_of(mode));
                let result = edlib_align(&query, &target, &param).unwrap();
                let expected = ref_distance(&query, &target, mode);
                assert_eq!(
                    result.edit_distance, expected,
                    "case {case}, mode {mode:?}, query {query:?}, target {target:?}"
                );
            }
        }
    }
}

#[test]
fn locations_match_naive_dp() {
    let mut rng = Lcg::new(7);
    let alphabet = b"ACGT";

    for _ in 0..300 {
        let query_len = rng.next_below(30);
        let target_len = rng.next_below(30);
        let query = rng.sequence(query_len, alphabet);
        let target = rng.sequence(target_len, alphabet);

        for mode in [RefMode::Nw, RefMode::Shw, RefMode::Hw] {
            let mut param = EdlibAlignParam::default();
            param.set_mode(mode_of(mode));
            param.set_task(AlignTask::Locations);
            let result = edlib_align(&query, &target, &param).unwrap();

            if result.edit_distance < 0 {
                assert!(result.locations.is_empty());
                continue;
            }
            assert!(
                !result.locations.is_empty(),
                "no locations for a valid alignment (mode {mode:?})"
            );

            for &(start, end) in &result.locations {
                if target.is_empty() {
                    // end location -1 is exposed as usize::MAX by the binding.
                    assert_eq!((start, end), (0, usize::MAX));
                    continue;
                }
                if mode == RefMode::Nw {
                    assert_eq!((start, end), (0, target.len() - 1));
                } else if mode == RefMode::Shw {
                    assert_eq!(start, 0);
                } else {
                    assert!(start <= end);
                }
                if end >= target.len() {
                    // HW can report end == -1 (query starts before the target);
                    // the binding exposes that as a huge usize, so skip it here.
                    continue;
                }
                // The reported span must reproduce the reported edit distance.
                let span_distance = ref_substring_distance(&query, &target, start, end);
                assert_eq!(
                    span_distance, result.edit_distance,
                    "mode {mode:?}, query {query:?}, target {target:?}, span {start}..={end}"
                );
            }
        }
    }
}

#[test]
fn path_and_cigar_are_consistent() {
    let mut rng = Lcg::new(2024);
    let alphabet = b"ACGT";

    for _ in 0..200 {
        let query_len = rng.next_below(35);
        let target_len = rng.next_below(35);
        let query = rng.sequence(query_len, alphabet);
        let target = rng.sequence(target_len, alphabet);

        for mode in [RefMode::Nw, RefMode::Shw, RefMode::Hw] {
            let mut param = EdlibAlignParam::default();
            param.set_mode(mode_of(mode));
            param.set_task(AlignTask::Path);
            param.set_cigar_fmt(CigarFmt::Extended);
            let result = edlib_align(&query, &target, &param).unwrap();

            if result.edit_distance < 0 {
                assert!(result.cigar.is_none());
                continue;
            }
            if query.is_empty() || target.is_empty() {
                // edlib returns before computing an alignment; the C binding then
                // reported an empty CIGAR for the path task.
                assert_eq!(result.cigar.as_deref(), Some(""));
                continue;
            }
            let cigar = result.cigar.expect("cigar expected for the path task");
            let (start, end) = result.locations[0];
            if target.is_empty() || end < start || end >= target.len() {
                continue;
            }
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

#[test]
fn standard_and_extended_cigar_describe_the_same_alignment() {
    let mut rng = Lcg::new(99);
    let alphabet = b"ACGT";

    for _ in 0..100 {
        let query_len = 1 + rng.next_below(25);
        let target_len = 1 + rng.next_below(25);
        let query = rng.sequence(query_len, alphabet);
        let target = rng.sequence(target_len, alphabet);

        let mut param = EdlibAlignParam::default();
        param.set_task(AlignTask::Path);
        param.set_cigar_fmt(CigarFmt::Extended);
        let extended = edlib_align(&query, &target, &param).unwrap();

        param.set_cigar_fmt(CigarFmt::Standard);
        let standard = edlib_align(&query, &target, &param).unwrap();

        let (start, end) = extended.locations[0];
        let span = end - start + 1;

        let ext = extended.cigar.unwrap();
        let std = standard.cigar.unwrap();
        check_cigar(
            &ext,
            true,
            query.len(),
            span,
            extended.edit_distance as usize,
        );
        check_cigar(
            &std,
            false,
            query.len(),
            span,
            extended.edit_distance as usize,
        );

        // Both formats must describe the same amount of I/D operations, and the
        // standard format must not contain '=' / 'X' operations.
        let edit_ops = |cigar: &str| -> Vec<(char, usize)> {
            support::parse_cigar(cigar)
                .into_iter()
                .filter(|(_, op)| matches!(op, 'I' | 'D'))
                .map(|(len, op)| (op, len))
                .collect()
        };
        assert_eq!(edit_ops(&ext), edit_ops(&std));
    }
}

#[test]
fn medium_sequences_use_multi_block_traceback() {
    let mut rng = Lcg::new(555);
    let alphabet = b"ACGT";

    for _ in 0..12 {
        let query_len = 200 + rng.next_below(300);
        let target_len = 200 + rng.next_below(300);
        let query = rng.sequence(query_len, alphabet);
        let target = rng.sequence(target_len, alphabet);

        for mode in [RefMode::Nw, RefMode::Shw, RefMode::Hw] {
            let mut param = EdlibAlignParam::default();
            param.set_mode(mode_of(mode));
            param.set_task(AlignTask::Path);
            param.set_cigar_fmt(CigarFmt::Extended);
            let result = edlib_align(&query, &target, &param).unwrap();

            assert_eq!(result.edit_distance, ref_distance(&query, &target, mode));

            let cigar = result.cigar.clone().unwrap();
            let (start, end) = result.locations[0];
            if end < start || end >= target.len() {
                continue;
            }
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
