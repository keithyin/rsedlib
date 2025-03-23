use std::env;
use std::process::Command;

fn main() {
    let build_output_dir = env::var("OUT_DIR").unwrap();
    Command::new("sh")
        .arg("-c")
        .arg(&format!("cp -r edlib {}/", build_output_dir))
        .output()
        .expect("copy edlib failed");

    let edlib_source_code_dir = format!("{build_output_dir}/edlib");

    Command::new("sh")
        .arg("-c")
        .arg("mkdir -p build && cd build && cmake -D CMAKE_BUILD_TYPE=Release .. && make")
        .current_dir(&edlib_source_code_dir)
        .output()
        .expect("cmake failed");
    //  try to build edlib in build.rs with edlib cloned in edlib (from which we removed )
    
    println!("cargo:rustc-link-search=native={}", format!("{}/build/lib", edlib_source_code_dir));

    // adapted fix from sam217pa. libstdc++ name varies...
    let lib_std: &str;
    if cfg!(target_os = "macos") {
        lib_std = "cargo:rustc-link-lib=c++";
    } else {
        lib_std = "cargo:rustc-link-lib=stdc++";
    }

    println!("cargo:rustc-link-lib=edlib");
    println!("{}", lib_std);

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    //
    println!("cargo:rerun-if-changed=edlib/");
}
