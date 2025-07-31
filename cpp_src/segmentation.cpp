#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <vector>
#include <iostream>

#ifdef __AVX512DQ__

// Inline assembly function to perform the operation
void segmentation_asm(const double *values, const double *segments, uint64_t *decisionBitMask, uint32_t *decisionIdx, size_t asize, double penalty) {
    asm volatile(
        "vzeroall\n"

        "vmovupd (%1), %%zmm24 \n"     // Load seg_values into zmm8
        "vmovupd 64(%1), %%zmm25 \n"     // Load seg_values into zmm9
        "vmovupd 128(%1), %%zmm26 \n"     // Load seg_values into zmm10
        "vmovupd 192(%1), %%zmm27 \n"     // Load seg_values into zmm11
        "vmovupd 256(%1), %%zmm28 \n"     // Load seg_values into zmm12
        "vmovupd 320(%1), %%zmm29 \n"     // Load seg_values into zmm13
        "vmovupd 384(%1), %%zmm30 \n"     // Load seg_values into zmm14
        "vmovupd 448(%1), %%zmm31 \n"     // Load seg_values into zmm15
        "vpbroadcastq (%0), %%zmm7 \n" // Load val_values[0] into zmm7 and broadcast

        "movq $0x7FFFFFFFFFFFFFFF, %%rax \n"   // Load the sign bit pattern into a register
        "vpbroadcastq %%rax, %%zmm6 \n"        // Broadcast the sign bit pattern to all elements of %%zmm6
        "vpbroadcastq %5, %%zmm5 \n" // broadcast penalty

        "xor %%rax, %%rax \n" // Initialize the loop counter to one
        "jmp minelem\n"

        // start loop

        "mainloop:\n"
        "vbroadcastsd (%0, %%rax, 8), %%zmm7\n" // Broadcast the next element element of the array to segment to all elements of zmm7
        "vaddpd %%zmm5, %%zmm0, %%zmm0 \n" // add penalty (in zmm5)

        //
        // Get the min score between break and no break and save decision
        //

        "vcmppd $1, %%zmm8, %%zmm0, %%k1 \n" // check if zmm0 (minscore plus penalty) is smaller than current score zmm8
        "vblendmpd %%zmm0, %%zmm8, %%zmm8 %{%%k1%} \n" // then set the smaller elements to zmm0
        "kmovb %%k1, (%2, %%rax, 8) \n"
        "vcmppd $1, %%zmm9, %%zmm0, %%k1 \n" // check if zmm0 (minscore) is smaller than zmm9
        "vblendmpd %%zmm0, %%zmm9, %%zmm9 %{%%k1%} \n" // then set the smaller elements to zmm0
        "kmovb %%k1, 1(%2, %%rax, 8) \n"
        "vcmppd $1, %%zmm10, %%zmm0, %%k1 \n" // check if zmm0 (minscore) is smaller than zmm10
        "vblendmpd %%zmm0, %%zmm10, %%zmm10 %{%%k1%} \n" // then set the smaller elements to zmm0
        "kmovb %%k1, 2(%2, %%rax, 8) \n"
        "vcmppd $1, %%zmm11, %%zmm0, %%k1 \n" // check if zmm0 (minscore) is smaller than zmm11
        "vblendmpd %%zmm0, %%zmm11, %%zmm11 %{%%k1%} \n" // then set the smaller elements to zmm0
        "kmovb %%k1, 3(%2, %%rax, 8) \n"
        "vcmppd $1, %%zmm12, %%zmm0, %%k1 \n" // check if zmm0 (minscore) is smaller than zmm12
        "vblendmpd %%zmm0, %%zmm12, %%zmm12 %{%%k1%} \n" // then set the smaller elements to zmm0
        "kmovb %%k1, 4(%2, %%rax, 8) \n"
        "vcmppd $1, %%zmm13, %%zmm0, %%k1 \n" // check if zmm0 (minscore) is smaller than zmm13
        "vblendmpd %%zmm0, %%zmm13, %%zmm13 %{%%k1%} \n" // then set the smaller elements to zmm0
        "kmovb %%k1, 5(%2, %%rax, 8) \n"
        "vcmppd $1, %%zmm14, %%zmm0, %%k1 \n" // check if zmm0 (minscore) is smaller than zmm14
        "vblendmpd %%zmm0, %%zmm14, %%zmm14 %{%%k1%} \n" // then set the smaller elements to zmm0
        "kmovb %%k1, 6(%2, %%rax, 8) \n"
        "vcmppd $1, %%zmm15, %%zmm0, %%k1 \n" // check if zmm0 (minscore) is smaller than zmm15
        "vblendmpd %%zmm0, %%zmm15, %%zmm15 %{%%k1%} \n" // then set the smaller elements to zmm0
        "kmovb %%k1, 7(%2, %%rax, 8) \n"

        // This is the bitmask for each of the 64 elements. Are we switching or staying.
        "minelem:\n"

        "vsubpd %%zmm24, %%zmm7, %%zmm16 \n" // subtract %%zmm8 from %%zmm7 and store in zmm16
        "vsubpd %%zmm25, %%zmm7, %%zmm17 \n" // subtract %%zmm8 from %%zmm7 and store in zmm17
        "vsubpd %%zmm26, %%zmm7, %%zmm18 \n" // subtract %%zmm8 from %%zmm7 and store in zmm18
        "vsubpd %%zmm27, %%zmm7, %%zmm19 \n" // subtract %%zmm8 from %%zmm7 and store in zmm19
        "vsubpd %%zmm28, %%zmm7, %%zmm20 \n" // subtract %%zmm8 from %%zmm7 and store in zmm20
        "vsubpd %%zmm29, %%zmm7, %%zmm21 \n" // subtract %%zmm8 from %%zmm7 and store in zmm21
        "vsubpd %%zmm30, %%zmm7, %%zmm22 \n" // subtract %%zmm8 from %%zmm7 and store in zmm22
        "vsubpd %%zmm31, %%zmm7, %%zmm23 \n" // subtract %%zmm8 from %%zmm7 and store in zmm23

        "vandpd %%zmm6, %%zmm16, %%zmm16 \n" // abs
        "vandpd %%zmm6, %%zmm17, %%zmm17 \n" // abs
        "vandpd %%zmm6, %%zmm18, %%zmm18 \n" // abs
        "vandpd %%zmm6, %%zmm19, %%zmm19 \n" // abs
        "vandpd %%zmm6, %%zmm20, %%zmm20 \n" // abs
        "vandpd %%zmm6, %%zmm21, %%zmm21 \n" // abs
        "vandpd %%zmm6, %%zmm22, %%zmm22 \n" // abs
        "vandpd %%zmm6, %%zmm23, %%zmm23 \n" // abs

        "vaddpd %%zmm16, %%zmm8, %%zmm8 \n" // add the difference of the predefined median from the minscore
        "vaddpd %%zmm17, %%zmm9, %%zmm9 \n" //
        "vaddpd %%zmm18, %%zmm10, %%zmm10 \n" //
        "vaddpd %%zmm19, %%zmm11, %%zmm11 \n" //
        "vaddpd %%zmm20, %%zmm12, %%zmm12 \n" //
        "vaddpd %%zmm21, %%zmm13, %%zmm13 \n" //
        "vaddpd %%zmm22, %%zmm14, %%zmm14 \n" //
        "vaddpd %%zmm23, %%zmm15, %%zmm15 \n" //

        "vminpd %%zmm8, %%zmm9, %%zmm16 \n"
        "vminpd %%zmm10, %%zmm11, %%zmm17 \n"
        "vminpd %%zmm12, %%zmm13, %%zmm18 \n"
        "vminpd %%zmm14, %%zmm15, %%zmm19 \n"

        "vminpd %%zmm16, %%zmm17, %%zmm20 \n"
        "vminpd %%zmm18, %%zmm19, %%zmm21 \n"

        "vminpd %%zmm20, %%zmm21, %%zmm0 \n"

        "vextractf64x4 $1, %%zmm0, %%ymm1 \n"  // Extract upper 256 bits of zmm0 to ymm1
        "vminpd %%ymm0, %%ymm1, %%ymm0 \n" // Find minimum
        "vextractf64x2 $1, %%ymm0, %%xmm1 \n" // Extract upper 128 bits of ymm0 to xmm1
        "vminpd %%xmm0, %%xmm1, %%xmm0 \n" // Find minimum of lower and upper halves
        "movhlps %%xmm0, %%xmm1 \n" // Move upper half of xmm0 to xmm1
        "vminpd %%xmm0, %%xmm1, %%xmm0 \n" // Find minimum of lower and upper halves

        "vpbroadcastq %%xmm0, %%zmm0 \n" // broadcast to zmm0

        "vpcmpeqq %%zmm0, %%zmm21, %%k2 \n"
        "ktestb %%k2, %%k2 \n"
        "jnz upper \n"

        "vpcmpeqq %%zmm0, %%zmm17, %%k2 \n"
        "ktestb %%k2, %%k2 \n"
        "jnz lower_upper \n"

        "vpcmpeqq %%zmm0, %%zmm9, %%k2 \n"
        "ktestb %%k2, %%k2 \n"
        "jnz lower_lower_upper \n"

        // LOWER HALF

        "vpcmpeqq %%zmm0, %%zmm20, %%k2 \n"
        "kmovq %%k2, %%rcx \n" // move the bitmask to rcx
        "tzcnt %%rcx, %%rcx \n" // get the index of the min element
        "jmp end\n"

        "lower_lower_upper:"
        "kmovq %%k2, %%rcx \n" // move the bitmask to rcx
        "tzcnt %%rcx, %%rcx \n" // get the index of the min element
        "or $8, %%rcx \n"
        "jmp end\n"

        "lower_upper:"
        "kmovq %%k2, %%rcx \n" // move the bitmask to rcx
        "tzcnt %%rcx, %%rcx \n" // get the index of the min element

        "vpcmpeqq %%zmm0, %%zmm11, %%k1 \n"
        "ktestb %%k1, %%k1 \n"
        "jnz lower_upper_upper \n"

        "or $16, %%rcx \n"
        "jmp end\n"

        "lower_upper_upper:"

        "or $24, %%rcx \n"
        "jmp end\n"

        // UPPER HALF

        "upper:\n"

        "kmovq %%k2, %%rcx \n" // move the bitmask to rcx
        "tzcnt %%rcx, %%rcx \n" // get the index of the min element

        "vpcmpeqq %%zmm0, %%zmm19, %%k1 \n"
        "ktestb %%k1, %%k1 \n"
        "jnz upper_upper \n"

        "vpcmpeqq %%zmm0, %%zmm13, %%k1 \n"
        "ktestb %%k1, %%k1 \n"
        "jnz upper_lower_upper \n"

        "or $32, %%rcx \n"
        "jmp end \n"

        "upper_lower_upper:\n"

        "or $40, %%rcx \n"
        "jmp end \n"

        "upper_upper:\n"

        "vpcmpeqq %%zmm0, %%zmm15, %%k1 \n"
        "ktestb %%k1, %%k1 \n"
        "jnz upper_upper_upper \n"

        "or $48, %%rcx \n"
        "jmp end \n"

        "upper_upper_upper:\n"

        "or $56, %%rcx \n"

        "end:\n"

        "mov %%ecx, (%3, %%rax, 4) \n" // save the index in memory

        "incq %%rax\n" // Increment the loop counter
        "cmpq %4, %%rax\n" // Compare the loop counter with the number of elements
        "jl mainloop" // Jump back if the loop counter is less than the number of elements

        :
        : "r"(values), "r"(segments), "r"(decisionBitMask), "r"(decisionIdx), "r"(asize), "r"(penalty)
        : "rax", "rcx", "zmm0", "zmm1", "zmm2", "zmm3", "zmm5", "zmm6", "zmm7", "zmm8", "zmm9", "zmm10", "zmm11", "zmm12", "zmm13", "zmm14", "zmm15", "zmm16", "zmm17", "zmm18", "zmm19", "zmm20", "zmm21", "zmm22", "zmm23", "zmm24", "zmm25", "zmm26", "zmm27", "zmm28", "zmm29", "zmm30", "zmm31"
    );
}

# endif





extern "C" int segment(double* val_values,
                       size_t n,
                       double* seg_values,
                       size_t s,
                       double penalty,
                       int* out_index,
                       double* out_values) {

    std::vector<uint32_t> backidx (n);
    std::vector<size_t> breakidx (n);
    // b starts with 1 because no breaks means one segment
    int b = 1;

    if (s<=64) {

    std::vector<double> score (s);
    std::vector<uint64_t> backboolbit (n);


#ifdef __AVX512DQ__

    if (s==64) {

        segmentation_asm(val_values, seg_values, backboolbit.data(), backidx.data(), n, penalty);

    } else {

# endif


    // first element is min
    double minscore = score[0] = fabs(val_values[0] - seg_values[0]);
    for (size_t j=1;j<s;++j) {
        score[j] = fabs(val_values[0] - seg_values[j]);

        if (minscore > score[j]){
            minscore = score[j];
            backidx[0] = j;
        }
    }

    minscore += penalty;

    for (size_t i=1;i<n;++i) {
        uint64_t bitmask = (uint64_t)0;
        if (minscore < score[0]) {
            bitmask |= ((uint64_t)1 );
            score[0] = minscore;
        }
        double jmin = score[0] += fabs(val_values[i] - seg_values[0]);

        for (size_t j=1;j<s;++j) {

            // double diff = score[j] - minscore;
            if (minscore < score[j]) {
                bitmask |= (uint64_t)1 << j;
                score[j] = minscore;
            }
            score[j] = score[j] + fabs(val_values[i] - seg_values[j]);

            if (jmin > score[j]) {
                jmin = score[j];
                backidx[i] = j;
            }
        }

        backboolbit[i] |= bitmask;
        minscore = jmin + penalty;

    }

    #ifdef __AVX512DQ__
    }
    #endif

    breakidx[0] = n;
    uint64_t bitmask = (uint64_t)1 << backidx[n-1];
    for (int i=n-1;i>0;--i) {
        if ((backboolbit[i]) & bitmask) {
            bitmask = (uint64_t)1 << backidx[i-1];
            breakidx[b] = i;
            b++;
        }
    }

    } else {

    std::vector<double> score (s);
    std::vector<bool> backbool (s*n, false);

    // FORWARD PASS
    // first element is min
    double minscore = score[0] = fabs(val_values[0] - seg_values[0]);
    for (size_t j=1;j<s;++j) {
        score[j] = fabs(val_values[0] - seg_values[j]);

        if (minscore > score[j]){
            minscore = score[j];
            backidx[0] = j;
        }
    }

    minscore += penalty;

    for (size_t i=1;i<n;++i) {
        if (minscore < score[0]) {
            backbool[s*i] = true;
            score[0] = minscore;
        }
        double jmin = score[0] = score[0] + fabs(val_values[i] - seg_values[0]);

        for (size_t j=1;j<s;++j) {

            // double diff = score[j] - minscore;
            if (minscore < score[j]) {
                backbool[s*i+j] = true;
                score[j] = minscore;
            }
            score[j] = score[j] + fabs(val_values[i] - seg_values[j]);

            if (jmin > score[j]) {
                jmin = score[j];
                backidx[i] = j;
            }
        }

        minscore = jmin + penalty;
    }

    // BACKTRACK
    breakidx[0] = n;
    size_t maxixtmp = backidx[n-1];
    for (size_t i=n-1;i>0;--i) {
        if (backbool[i*s+maxixtmp]) {
            maxixtmp = backidx[i-1];
            breakidx[b] = i;
            b++;
        }
    }

    }

    out_index[0] = 0;
    out_values[0] = seg_values[backidx[breakidx[b-1]-1]];

    for (int i=b-2;i>=0;i--) {
        out_index[b-(i+1)] = breakidx[i+1];
        out_values[b-(i+1)] = seg_values[backidx[breakidx[i]-1]];
    }

    return b;


}

