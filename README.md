# rsedlib

[![CI](https://github.com/keithyin/rsedlib/actions/workflows/ci.yml/badge.svg)](https://github.com/keithyin/rsedlib/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/keithyin/rsedlib)](https://github.com/keithyin/rsedlib/releases)
[![crates.io](https://img.shields.io/crates/v/rsedlib)](https://crates.io/crates/rsedlib)
[![docs.rs](https://img.shields.io/docsrs/rsedlib)](https://docs.rs/rsedlib)
[![unsafe forbidden](https://img.shields.io/badge/unsafe-forbidden-success)](https://github.com/rust-secure-code/safety-dance/)
[![License](https://img.shields.io/badge/license-MIT-blue)](#license)
[![Rust](https://img.shields.io/badge/rust-2021-dea584?logo=rust)](https://www.rust-lang.org/)

[edlib](https://github.com/Martinsos/edlib) 的**纯 Rust 实现**（不是绑定）。

编辑距离（Levenshtein distance）计算的位并行 Myers 算法、Ukkonen 带状优化、k 限制与
Hirschberg 分治回溯路径还原全部用 Rust 重写，因此：

- 不需要 C++ 编译器、`cmake`、`libc` 或 git submodule；
- `cargo build` 开箱即用，可交叉编译到任意 Rust 目标平台；
- 全程 `#![deny(unsafe_code)]`，没有裸指针与手工内存管理。

## 用法

```rust
use rsedlib::edlib_align;
use rsedlib::param::{AlignTask, CigarFmt, EdlibAlignParam};

let query = b"elephant";
let target = b"telephone";

let mut param = EdlibAlignParam::default();
param.set_task(AlignTask::Path);
param.set_cigar_fmt(CigarFmt::Extended);
let aln_res = edlib_align(query, target, &param).unwrap();
println!("{:?}", aln_res);

/*
Ok(EdlibAlignResult { edit_distance: 3, alphabet_length: 8, locations: [(0, 8)], cigar: Some("1D5=1X1=1X") })
*/
```

## 支持的能力

| 模式 | 说明 |
| --- | --- |
| `AlignMode::Global` | 全局比对（Needleman-Wunsch），query 首尾 gap 都计分 |
| `AlignMode::Prefix` | 前缀比对（SHW），target 末尾多余的字符免费 |
| `AlignMode::Infix` | 中缀比对（HW），query 可匹配 target 中任意子串 |

| 任务 | 说明 |
| --- | --- |
| `AlignTask::Distance` | 只求编辑距离 + target 中的结束位置 |
| `AlignTask::Locations` | 额外求起始位置 |
| `AlignTask::Path` | 额外求对齐路径，可输出 CIGAR（Standard: `M/I/D`，Extended: `=/X/I/D`） |

- `k` 限制：`k < 0`（默认）自动倍增搜索；`k >= 0` 时若编辑距离大于 `k` 则返回
  `edit_distance == -1`，行为与 C 版一致。
- 自定义等价字符：`param.set_eq_pairs(vec![(b'A', b'a')])` 可实现大小写不敏感或 `N` 通配。

## 与旧版绑定的差异（迁移提示）

- 删除了 `build.rs`、`edlib/` 子模块、`libc` 依赖与 `edlib_sys` 模块（bindgen 生成的 C 声明）。
- `EdlibEqualityPair` 现在定义在 `param` 模块中，字段类型由 `c_char` 改为 `u8`；
  `EdlibAlignParam` 的公开方法签名保持不变。
- `edlib_align(query, target, &param)`、`EdlibAlignResult` 的字段与语义保持不变，
  调用方代码无需修改。

## 测试与性能

```bash
cargo test                       # 单元测试 + 与朴素 DP 的随机对拍 + 黄金用例
cargo run --release --example bench
```

`tests/` 中包含：与 O(n·m) 朴素动态规划的距离/位置对拍、模式 × 任务 × k 组合矩阵、
CIGAR 展开校验（消耗长度与编辑距离一致性）、空序列/等价字符等边界，以及触发
Hirschberg 分治路径的大序列用例。

参考计时（单线程，Release）：

```
global/Distance     q=10000 -> 1.5 ms
global/Path         q=10000 -> 7.6 ms
```

## License

MIT，与原 edlib（MIT）及旧绑定保持一致。
