#pragma once

typedef float wide_float;
typedef int wide_int;

// struct wide_float {
//     float V[8];
// };
//
// struct wide_int {
//     int V[8];
// };

wide_float fabs(wide_float d){
    return (d>=0) ?  d: -d;
};

wide_float max(wide_float l, wide_float r){
    return (l>r) ? l: r;
};

void ConditionalAssign(wide_int* dest, wide_int mask, wide_int source){
    mask = mask ? 0xFFFFFFFF : 0;
    *dest =  ( (~mask & *dest) | (mask & source) );
}

void ConditionalAssign(wide_float* dest, wide_int mask, wide_float source){
    ConditionalAssign( (wide_int*)dest, mask, *(wide_int*)&source);
}


float HorizontalAdd(wide_float A) {
    return A;
}

int HorizontalAdd(wide_int A) {
    return A;
}

bool IsEmpty(wide_int A) {
    return A == 0 ? true : false;
}

wide_int WideIndex(int i) {
    return i;
}


// __m256 fabs(const __m256& wide){
//     __m256i mask = _mm256_set1_epi32(0x7FFFFFFF);
//     mask = _mm256_and_si256(_mm256_castps_si256(wide), mask);
//     return _mm256_castsi256_ps(mask);
// };
