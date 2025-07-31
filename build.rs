fn main() {
    // println!("cargo:rustc-link-search=native=./cpp_src/");
    // println!("cargo:rustc-link-lib=dylib=segmentation"); // libx.so
    // println!("cargo:rustc-link-search=native=./cpp_src/");
    // println!("cargo:rustc-link-lib=static=segmentation"); // libx.so
    
    cc::Build::new()
        .cpp(true)
        .file("cpp_src/segmentation.cpp")
        .compile("segmentation"); // creates libx.a

    println!("cargo:rerun-if-changed=cpp_src/segmentation.cpp");
    
}
