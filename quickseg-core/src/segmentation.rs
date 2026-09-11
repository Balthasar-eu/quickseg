
#[inline]
pub fn set_bit(bits: &mut [u8], idx: usize) {
    bits[idx / 8] |= 1u8 << (idx % 8);
}

#[inline]
fn get_bit(bits: &[u8], idx: usize) -> bool {
    (bits[idx / 8] & (1u8 << (idx % 8))) != 0
}

pub fn segment_chromosome(
    values: &[f64],
    seg_values: &[f64],
    penalty: f64,
) -> (Vec<usize>, Vec<f64>) {

    assert!(
        !seg_values.is_empty(),
        "seg_values must not be empty"
    );

    assert!(
        !values.is_empty(),
        "Values must not be empty"
    );

    let mut seg_values = seg_values.to_vec();
    let padded_len = seg_values.len().next_multiple_of(8);
    // Pad with +inf so the SIMD implementations can process
    // full 8-wide vectors without affecting the argmin.
    seg_values.resize(padded_len, f64::INFINITY);

    let n = values.len();
    let s = seg_values.len();

    let num_bits = s * n;
    // this would need div_ceil(8), but we ensure above that s is a multiple of 8
    let mut backbool = vec![0u8; num_bits / 8];
    let mut backidx = vec![0usize; n];
    let mut breakidx = vec![0usize; n];

    segment_forward(
        values,
        &seg_values,
        penalty,
        &mut backbool,
        &mut backidx,
    );

    let mut b = 1;

    breakidx[0] = n;
    let mut state = backidx[n - 1];

    for i in (1..n).rev() {
        if get_bit(&backbool, i * s + state) {
            state = backidx[i - 1];
            breakidx[b] = i;
            b += 1;
        }
    }

    let mut out_index = vec![0usize; b];
    let mut out_values = vec![0.0; b];

    // breakidx stores breakpoints in reverse order with a sentinel at position 0.
    for i in 0..b {
        out_index[i] = breakidx[b - i];
        out_values[i] = seg_values[backidx[breakidx[b - i - 1] - 1]];
    }

    (out_index, out_values)
}

pub fn segment_forward(
    values: &[f64],
    seg_values: &[f64],
    penalty: f64,
    backbool: &mut [u8],
    backidx: &mut [usize],
) {
    #[cfg(target_arch = "x86_64")]
    {
        if let Ok(backend) = std::env::var("QUICKSEG_BACKEND") {
            match backend.as_str() {
                "scalar" => {
                    return segment_forward_base(
                        values,
                        seg_values,
                        penalty,
                        backbool,
                        backidx,
                    );
                }
                "avx2" => {
                    assert!(
                        std::arch::is_x86_feature_detected!("avx2"),
                        "QUICKSEG_BACKEND=avx2 but AVX2 is not supported by this CPU"
                    );

                    return unsafe {
                        segment_forward_asm_avx2(
                            values,
                            seg_values,
                            penalty,
                            backbool,
                            backidx,
                        )
                    };
                }
                "avx512" => {
                    assert!(
                        std::arch::is_x86_feature_detected!("avx512f"),
                        "QUICKSEG_BACKEND=avx512 but AVX512 is not supported by this CPU"
                    );

                    return unsafe {
                        segment_forward_asm_avx512(
                            values,
                            seg_values,
                            penalty,
                            backbool,
                            backidx,
                        )
                    };
                }
                _ => {}
            }
        }

        if std::arch::is_x86_feature_detected!("avx512f") {
            return unsafe {
                segment_forward_asm_avx512(
                    values,
                    seg_values,
                    penalty,
                    backbool,
                    backidx,
                )
            };
        }

        if std::arch::is_x86_feature_detected!("avx2") {
            return unsafe {
                segment_forward_asm_avx2(
                    values,
                    seg_values,
                    penalty,
                    backbool,
                    backidx,
                )
            };
        }
    }

    segment_forward_base(
        values,
        seg_values,
        penalty,
        backbool,
        backidx,
    )
}

fn segment_forward_base(
    values: &[f64],
    seg_values: &[f64],
    penalty: f64,
    backbool: &mut [u8],
    backidx: &mut [usize],
){

    let n = values.len();
    let s = seg_values.len();

    let mut score = vec![0.0; s];
    let mut minscore = f64::INFINITY;

    for i in 0..n {
        let v = values[i];
        let mut jmin = f64::INFINITY;

        for j in 0..s {
            if minscore < score[j] {
                set_bit(backbool, s * i + j);
                score[j] = minscore;
            }

            score[j] += (v - seg_values[j]).abs();

            if jmin > score[j] {
                jmin = score[j];
                backidx[i] = j;
            }
        }

        minscore = jmin + penalty;
    }

}

#[cfg(target_arch = "x86_64")]
unsafe fn segment_forward_asm_avx2(
    values: &[f64],
    seg_values: &[f64],
    penalty: f64,
    backbool: &mut [u8],
    backidx: &mut [usize],
) {
    let n = values.len();
    let s = seg_values.len();
    let mut score = vec![0.0; s];

    let values_ptr = values.as_ptr();
    let seg_ptr = seg_values.as_ptr();
    let score_ptr = score.as_mut_ptr();
    let bits_ptr = backbool.as_mut_ptr();
    let idx_ptr = backidx.as_mut_ptr();

    let inf = f64::INFINITY;

    unsafe {
    core::arch::asm!(
        // signbit mask
        "mov rax, 0x7FFFFFFFFFFFFFFF",
        "vmovq xmm2, rax",
        "vbroadcastsd ymm2, xmm2",

        // penalty broadcast
        "vbroadcastsd ymm0, xmm0",

        // minscore is infinite
        "vbroadcastsd ymm1, xmm15", // setup minscore

        "xor r13, r13", // initialize i = 0
        "xor r14, r14", // address for bitarray

        "2:",

        "xor r11, r11", // bitmask = 0
        "xor r12, r12", // j = 0
        "xor r15, r15", // current best idx
        "vbroadcastsd ymm8, xmm15", // setup jmin
        "vbroadcastsd ymm3, [rdi + r13*8]", // broadcast value

        "3:",

        // load score[j..j+3]
        "vmovupd ymm4, [rdx + r12*8]",

        // load score[j+4..j+7]
        "vmovupd ymm10, [rdx + r12*8 + 32]",

        // compare minscore < score
        "vcmppd ymm5, ymm1, ymm4, 1",
        "vcmppd ymm11, ymm1, ymm10, 1",

        // extract masks
        "vmovmskpd eax, ymm5",
        "vmovmskpd r11d, ymm11",

        // combine into one byte mask
        "shl r11d, 4",
        "or eax, r11d",

        // store bitmask
        "mov byte ptr [rcx + r14], al",

        // score = min(score, minscore)
        "vminpd ymm4, ymm4, ymm1",
        "vminpd ymm10, ymm10, ymm1",

        // load seg_values
        "vmovupd ymm6, [rsi + r12*8]",
        "vmovupd ymm12, [rsi + r12*8 + 32]",

        // abs(v - seg)
        "vsubpd ymm7, ymm3, ymm6",
        "vandpd ymm7, ymm7, ymm2",

        "vsubpd ymm13, ymm3, ymm12",
        "vandpd ymm13, ymm13, ymm2",

        // score += abs(...)
        "vaddpd ymm4, ymm4, ymm7",
        "vaddpd ymm10, ymm10, ymm13",

        // store scores
        "vmovupd [rdx + r12*8], ymm4",
        "vmovupd [rdx + r12*8 + 32], ymm10",

        // running minimum
        "vminpd ymm8, ymm8, ymm4",
        "vminpd ymm8, ymm8, ymm10",

        // reduce ymm8 -> xmm8[0]
        "vextractf128 xmm9, ymm8, 1",
        "vminpd xmm8, xmm8, xmm9",

        "vunpckhpd xmm9, xmm8, xmm8",
        "vminsd xmm8, xmm8, xmm9",

        // broadcast jmin
        "vbroadcastsd ymm8, xmm8",

        // check first block
        "vcmppd ymm12, ymm4, ymm8, 0",
        "vmovmskpd eax, ymm12",

        "test eax, eax",
        "jz 4f",

        "tzcnt eax, eax",
        "lea r15, [r12 + rax]",
        "jmp 5f",

        "4:",

        // check second block
        "vcmppd ymm13, ymm10, ymm8, 0",
        "vmovmskpd eax, ymm13",

        "test eax, eax",
        "jz 5f",

        "tzcnt eax, eax",
        "lea r15, [r12 + rax + 4]",

        "5:",

        "add r14, 1",
        "add r12, 8",
        "cmp r12, r10",
        "jl 3b",

        // backidx[i] = r15
        "mov [r8 + r13*8], r15",

        "vbroadcastsd ymm8, xmm8",
        "vaddpd ymm1, ymm8, ymm0",

        "add r13, 1",
        "cmp r13, r9",
        "jl 2b",

        "vextractf128 xmm1, ymm1, 0",
        "vextractf128 xmm8, ymm8, 0",


        in("rdi") values_ptr,
        in("rsi") seg_ptr,
        in("rdx") score_ptr,
        in("rcx") bits_ptr,
        in("r8") idx_ptr,
        in("r9") n,
        in("r10") s,
        in("xmm0") penalty,
        in("xmm15") inf,
        lateout("rax") _,
        lateout("r11") _,
        lateout("r12") _,
        lateout("r13") _,
        lateout("r14") _,
        lateout("r15") _,
        lateout("ymm1") _,
        lateout("ymm2") _,
        lateout("ymm3") _,
        lateout("ymm4") _,
        lateout("ymm5") _,
        lateout("ymm6") _,
        lateout("ymm7") _,
        lateout("ymm8") _,
        lateout("ymm9") _,
        lateout("ymm10") _,
        lateout("ymm11") _,
        lateout("ymm12") _,
        lateout("ymm13") _,
    );
    }

}


#[repr(C, align(64))]
struct AlignedIdxBase([u64; 8]);

static IDX_BASE: AlignedIdxBase =
    AlignedIdxBase([0, 1, 2, 3, 4, 5, 6, 7]);

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512dq,avx512vl")]
unsafe fn segment_forward_asm_avx512(
    values: &[f64],
    seg_values: &[f64],
    penalty: f64,
    backbool: &mut [u8],
    backidx: &mut [usize],
)  {
    let n = values.len();
    let s = seg_values.len();
    let mut score = vec![0.0; s];

    let values_ptr = values.as_ptr();
    let seg_ptr = seg_values.as_ptr();
    let score_ptr = score.as_mut_ptr();
    let bits_ptr = backbool.as_mut_ptr();
    let idx_ptr = backidx.as_mut_ptr();

    let inf = f64::INFINITY;


    unsafe {
    core::arch::asm!(

        "mov rax, 0x7FFFFFFFFFFFFFFF",
        "vmovq xmm2, rax",
        "vbroadcastsd zmm2, xmm2",

        // penalty broadcast
        "vbroadcastsd zmm0, xmm0",

        // minscore is infinite
        "vbroadcastsd zmm1, xmm15", // setup minscore

        "xor r13, r13", // initialize i = 0
        "xor r14, r14", // address for bitarray

        //index vector initialize
        "vmovupd zmm14, [rip + {idx}]",
        "mov rax, 8",
        "vmovq xmm12, rax",
        "vpbroadcastq zmm12, xmm12",


        "2:",

        "xor r11, r11", // bitmask = 0
        "xor r12, r12", // j = 0
        "vpxord zmm13, zmm13, zmm13", // current idx
        "vpaddq zmm13, zmm13, zmm14", // zmm10 = [j+0 .. j+7]
        "vpxord zmm9, zmm9, zmm9", // current best idx
        "vbroadcastsd zmm8, xmm15", // setup jmin
        "vbroadcastsd zmm3, [rdi + r13*8]", // broadcast value

        "3:",

        // load score[j..j+7]
        "vmovupd zmm4, [rdx + r12*8]",

        // compare minscore < score
        "vcmppd k1, zmm1, zmm4, 1",
        "kmovb eax, k1",

        // store bitmask
        "mov byte ptr [rcx + r14], al",

        // score = min(score, minscore)
        "vminpd zmm4, zmm4, zmm1",

        // load seg_values
        "vmovupd zmm6, [rsi + r12*8]",

        // abs(v - seg)
        "vsubpd zmm7, zmm3, zmm6",
        "vandpd zmm7, zmm7, zmm2",

        // score += abs(...)
        "vaddpd zmm4, zmm4, zmm7",

        // store scores
        "vmovupd [rdx + r12*8], zmm4",

        // k1 = score < best_score
        "vcmppd k1, zmm4, zmm8, 1",

        // update best scores
        "vminpd zmm8, zmm8, zmm4",

        // update best indices
        "vmovdqa64 zmm9{{k1}}, zmm13",
        "vpaddq zmm13, zmm13, zmm12",

        "add r14, 1",
        "add r12, 8",
        "cmp r12, r10",
        "jl 3b",

        // values
        "vextractf64x4 ymm10, zmm8, 1",
        // indices
        "vextracti64x4 ymm11, zmm9, 1",
        // mask = upper < lower
        "vcmppd k1, ymm10, ymm8, 1",
        // values
        "vminpd ymm8, ymm8, ymm10",
        // indices
        "vmovdqa64 ymm9{{k1}}, ymm11",

        "vextractf128 xmm10, ymm8, 1",
        "vextracti128 xmm11, ymm9, 1",

        "vcmppd k1, xmm10, xmm8, 1",

        "vminpd xmm8, xmm8, xmm10",
        "vmovdqa64 xmm9{{k1}}, xmm11",

        "vpermilpd xmm10, xmm8, 1",
        "pshufd xmm11, xmm9, 0x4E",

        "vcmppd k1, xmm10, xmm8, 1",

        "vminsd xmm8, xmm8, xmm10",
        "vmovdqa64 xmm9{{k1}}, xmm11",

        // backidx[i]
        "vmovq [r8 + r13*8], xmm9",

        "vbroadcastsd zmm8, xmm8",
        "vaddpd zmm1, zmm8, zmm0",

        "add r13, 1",
        "cmp r13, r9",
        "jl 2b",

        in("rdi") values_ptr,
        in("rsi") seg_ptr,
        in("rdx") score_ptr,
        in("rcx") bits_ptr,
        in("r8") idx_ptr,
        in("r9") n,
        in("r10") s,
        in("xmm0") penalty,
        in("xmm15") inf,
        lateout("rax") _,
        lateout("r11") _,
        lateout("r12") _,
        lateout("r13") _,
        lateout("r14") _,
        lateout("zmm1") _,
        lateout("zmm2") _,
        lateout("zmm3") _,
        lateout("zmm4") _,
        lateout("zmm5") _,
        lateout("zmm6") _,
        lateout("zmm7") _,
        lateout("zmm8") _,
        lateout("zmm9") _,
        lateout("zmm10") _,
        lateout("zmm11") _,
        lateout("zmm12") _,
        lateout("zmm13") _,
        lateout("k1") _,
        idx = sym IDX_BASE,
    );
    }

}
