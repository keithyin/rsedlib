

```Rust

let query = b"elephant";
let target = b"telephone";

let mut param = EdlibAlignParam::default();
// param.set_mode(param::AlignMode::Global);
param.set_task(param::AlignTask::Path);
param.set_cigar_fmt(param::CigarFmt::Extended);
let aln_res = edlib_align(query, target, &param);
println!("{:?}", aln_res);

/*
Ok(EdlibAlignResult { edit_distance: 3, alphabet_length: 8, locations: [(0, 8)], cigar: Some("1D5=1X1=1X") })
*/

```