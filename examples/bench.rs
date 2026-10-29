//! Rough timing example, not meant as a CI gate.
//!
//! Run with `cargo run --release --example bench`.

use std::time::Instant;

use rsedlib::edlib_align;
use rsedlib::param::{AlignMode, AlignTask, CigarFmt, EdlibAlignParam};

struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Lcg(seed)
    }

    fn next_u32(&mut self) -> u32 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 33) as u32
    }

    fn next_below(&mut self, n: usize) -> usize {
        (self.next_u32() as usize) % n
    }
}

/// Random DNA sequence.
fn random_sequence(rng: &mut Lcg, len: usize) -> Vec<u8> {
    const ALPHABET: &[u8] = b"ACGT";
    (0..len).map(|_| ALPHABET[rng.next_below(4)]).collect()
}

/// Copy of `seq` with ~5% substitutions and ~1% indels.
fn mutate(rng: &mut Lcg, seq: &[u8]) -> Vec<u8> {
    const ALPHABET: &[u8] = b"ACGT";
    let mut result = Vec::with_capacity(seq.len() + seq.len() / 50);
    for &base in seq {
        match rng.next_below(100) {
            0 => {}                                        // deletion
            1 => result.push(ALPHABET[rng.next_below(4)]), // substitution
            2 => {
                result.push(base);
                result.push(ALPHABET[rng.next_below(4)]); // insertion
            }
            _ => result.push(base),
        }
    }
    result
}

fn bench(label: &str, query: &[u8], target: &[u8], param: &EdlibAlignParam) {
    let start = Instant::now();
    let result = edlib_align(query, target, param).unwrap();
    let elapsed = start.elapsed();

    println!(
        "{label:<28} q={:<6} t={:<6} -> distance={:<5} {:>10.3} ms",
        query.len(),
        target.len(),
        result.edit_distance,
        elapsed.as_secs_f64() * 1000.0
    );
}

fn main() {
    let mut rng = Lcg::new(0xDEAD_BEEF);

    for len in [1_000usize, 10_000] {
        let query = random_sequence(&mut rng, len);
        let target = mutate(&mut rng, &query);

        for (task, cigar_fmt) in [
            (AlignTask::Distance, CigarFmt::NoCigar),
            (AlignTask::Locations, CigarFmt::NoCigar),
            (AlignTask::Path, CigarFmt::Extended),
        ] {
            let mut param = EdlibAlignParam::default();
            param.set_mode(AlignMode::Global);
            param.set_task(task);
            param.set_cigar_fmt(cigar_fmt);
            bench(&format!("global/{task:?}"), &query, &target, &param);
        }

        let mut param = EdlibAlignParam::default();
        param.set_mode(AlignMode::Infix);
        param.set_task(AlignTask::Locations);
        param.set_k(50); // small k -> narrow band
        bench("infix/Locations k=50", &query, &target, &param);
    }
}
