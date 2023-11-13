#include <immintrin.h>

#define LANE_WIDTH 8

// AVX512
#if (LANE_WIDTH==16)


//======
// AVX2
//======
#elif (LANE_WIDTH==8)
struct wide_float {
    __m256 V;
    wide_float(){};
    wide_float(float);
    wide_float(float*);
    wide_float& operator=(float);
};
struct wide_int {
    __m256i V;
    wide_int(){};
    wide_int(int);
    wide_int(int*);
    wide_int& operator=(int);
};


//-----------
// MAKE WIDE
inline wide_float
WideFloatFromFloat(float A){
    wide_float ret;
    ret.V = _mm256_set1_ps(A);
    return ret;
}
inline wide_int
WideIntFromInt(int A){
    wide_int ret;
    ret.V = _mm256_set1_epi32(A);
    return ret;
}
inline wide_float
WideFloatFromWideInt(wide_int A) {
    wide_float ret;
    ret.V = _mm256_cvtepi32_ps(A.V);
    return ret;
}
inline wide_float
WideFloatFromInt(int A) {
    wide_float ret = WideFloatFromFloat((float) A);
    return ret;
}

//--------------
// CONSTRUCTORS
wide_int::wide_int(int A) {
    this->V = WideIntFromInt(A).V;
}
wide_int::wide_int(int* A) {
    this->V = _mm256_loadu_si256((__m256i*)A);
}
wide_int& wide_int::operator=(int A) {
    this->V = WideIntFromInt(A).V;
    return *this;
}
wide_float::wide_float(float A) {
    this->V = WideFloatFromFloat(A).V;
}
wide_float::wide_float(float* A) {
    this->V = _mm256_loadu_ps(A);
}
wide_float& wide_float::operator=(float A) {
    this->V = WideFloatFromFloat(A).V;
    return *this;
}

//--------------
// LOAD - STORE
inline wide_int
LoadPackedWideInt(int* A) {
    wide_int ret;
    ret.V = _mm256_loadu_si256((__m256i*)A);
    return ret;
}
inline wide_float
LoadPackedWideFloat(float* A) {
    wide_float ret;
    ret.V = _mm256_loadu_ps(A);
    return ret;
}
inline void
StoreWideInt(int* dest, wide_int source) {
    _mm256_storeu_si256((__m256i*)dest, source.V);
}
inline void
StoreWideFloat(float* dest, wide_float source) {
    _mm256_storeu_ps(dest, source.V);
}

//----------
// WIDE_INT
inline wide_int
operator==(wide_int A, wide_int B) {
    wide_int ret;
    ret.V = _mm256_cmpeq_epi32(A.V, B.V);
    return ret;
}
inline wide_int
operator==(int A, wide_int B) {
    wide_int ret = (WideIntFromInt(A) == B);
    return ret;
}
inline wide_int
operator==(wide_int A, int B) {
    wide_int ret = (A == WideIntFromInt(B));
    return ret;
}

inline wide_int
operator!=(wide_int A, wide_int B) {
    wide_int ret;
    ret.V = _mm256_andnot_si256((A==B).V, _mm256_set1_epi32(0xFFFFFFFF));
    return ret;
}
inline wide_int
operator!=(int A, wide_int B) {
    wide_int ret = (WideIntFromInt(A)!=B);
    return ret;
}
inline wide_int
operator!=(wide_int A, int B) {
    wide_int ret = (A!=WideIntFromInt(B));
    return ret;
}

inline wide_int
operator|(wide_int A, wide_int B) {
    wide_int ret;
    ret.V = _mm256_or_si256(A.V, B.V);
    return ret;
}
inline wide_int
operator|(int A, wide_int B) {
    wide_int ret = WideIntFromInt(A) | B;
    return ret;
}
inline wide_int
operator|(wide_int A, int B) {
    wide_int ret = A | WideIntFromInt(B);
    return ret;
}
inline wide_int&
operator|=(wide_int& A, wide_int B) {
    A = A | B;
    return A;
}
inline wide_int&
operator|=(wide_int& A, int B) {
    A |= WideIntFromInt(B);
    return A;
}

inline wide_int
operator&(wide_int A, wide_int B) {
    wide_int ret;
    ret.V = _mm256_and_si256(A.V, B.V);
    return ret;
}
inline wide_int
operator&(int A, wide_int B) {
    wide_int ret = WideIntFromInt(A) & B;
    return ret;
}
inline wide_int
operator&(wide_int A, int B) {
    wide_int ret = A & WideIntFromInt(B);
    return ret;
}

inline wide_int
operator>(wide_int A, wide_int B) {
    wide_int ret;
    ret.V = _mm256_cmpgt_epi32(A.V, B.V);
    return ret;
}
inline wide_int
operator>(int A, wide_int B) {
    wide_int ret = WideIntFromInt(A)>B;
    return ret;
}
inline wide_int
operator>(wide_int A, int B) {
    wide_int ret = A>WideIntFromInt(B);
    return ret;
}

inline wide_int
operator>=(wide_int A, wide_int B) {
    wide_int ret = (A > B) | (A == B);
    return ret;
}
inline wide_int
operator>=(int A, wide_int B) {
    wide_int ret = (WideIntFromInt(A)>=B);
    return ret;
}
inline wide_int
operator>=(wide_int A, int B) {
    wide_int ret = (A>=WideIntFromInt(B));
    return ret;
}

inline wide_int
operator<(wide_int A, wide_int B) {
    return (B > A);
}
inline wide_int
operator<(int A, wide_int B) {
    return (B > WideIntFromInt(A));
}
inline wide_int
operator<(wide_int A, int B) {
    return (WideIntFromInt(B) > A);
}

inline wide_int
operator<=(wide_int A, wide_int B) {
    return (B >= A);
}
inline wide_int
operator<=(int A, wide_int B) {
    return (B >= WideIntFromInt(A));
}
inline wide_int
operator<=(wide_int A, int B) {
    return (WideIntFromInt(B) >= A);
}

inline wide_int
operator+(wide_int A, wide_int B) {
    wide_int ret;
    ret.V = _mm256_add_epi32(A.V, B.V);
    return ret;
}
inline wide_int
operator+(int A, wide_int B) {
    wide_int ret = WideIntFromInt(A) + B;
    return ret;
}
inline wide_int
operator+(wide_int A, int B) {
    wide_int ret = A + WideIntFromInt(B);
    return ret;
}

inline wide_int&
operator+=(wide_int& A, wide_int B) {
    A = A + B;
    return A;
}
inline wide_int&
operator+=(wide_int& A, int B) {
    A += WideIntFromInt(B);
    return A;
}

inline wide_int
operator-(wide_int A) {
    wide_int ret;
    ret.V = 0 - A.V;
    return ret;
}
inline wide_int
operator-(wide_int A, wide_int B) {
    wide_int ret;
    ret.V = _mm256_sub_epi32(A.V, B.V);
    return ret;
}
inline wide_int
operator-(int A, wide_int B) {
    wide_int ret = WideIntFromInt(A)-B;
    return ret;
}
inline wide_int
operator-(wide_int A, int B) {
    wide_int ret = A-WideIntFromInt(B);
    return ret;
}

inline wide_int&
operator-=(wide_int& A, wide_int B) {
    A = A - B;
    return A;
}
inline wide_int&
operator-=(wide_int& A, int B) {
    A -= WideIntFromInt(B);
    return A;
}

// NOTE(alex): This is darn slow, maybe we should never even do this on int
inline wide_int
operator*(wide_int A, wide_int B) {
    wide_int ret;
    ret.V = _mm256_mullo_epi32(A.V, B.V);
    return ret;
}
inline wide_int
operator*(int A, wide_int B) {
    wide_int ret = WideIntFromInt(A)*B;
    return ret;
}
inline wide_int
operator*(wide_int A, int B) {
    wide_int ret = A*WideIntFromInt(B);
    return ret;
}

inline wide_int&
operator*=(wide_int& A, wide_int B) {
    A = A * B;
    return A;
}
inline wide_int&
operator*=(wide_int& A, int B) {
    A *= WideIntFromInt(B);
    return A;
}


//------------
// WIDE_FLOAT
inline wide_int
operator>(wide_float A, wide_float B) {
    wide_int ret;
    ret.V = _mm256_castps_si256(_mm256_cmp_ps(A.V, B.V, _CMP_GT_OQ));
    return ret;
}
inline wide_int
operator>(wide_float A, float B) {
    return (A > WideFloatFromFloat(B));
}
inline wide_int
operator>(float A, wide_float B) {
    return (WideFloatFromFloat(A) > B);
}

inline wide_int
operator>=(wide_float A, wide_float B) {
    wide_int ret;
    ret.V = _mm256_castps_si256(_mm256_cmp_ps(A.V, B.V, _CMP_GE_OQ));
    return ret;
}
inline wide_int
operator>=(wide_float A, float B) {
    return (A >= WideFloatFromFloat(B));
}
inline wide_int
operator>=(float A, wide_float B) {
    return (WideFloatFromFloat(A) >= B);
}

inline wide_int
operator<(wide_float A, wide_float B) {
    wide_int ret;
    ret.V = _mm256_castps_si256(_mm256_cmp_ps(A.V, B.V, _CMP_LT_OQ));
    return ret;
}
inline wide_int
operator<(wide_float A, float B) {
    return (A < WideFloatFromFloat(B));
}
inline wide_int
operator<(float A, wide_float B) {
    return (WideFloatFromFloat(A) < B);
}

inline wide_int
operator<=(wide_float A, wide_float B) {
    wide_int ret;
    ret.V = _mm256_castps_si256(_mm256_cmp_ps(A.V, B.V, _CMP_LE_OQ));
    return ret;
}
inline wide_int
operator<=(wide_float A, float B) {
    return (A <= WideFloatFromFloat(B));
}
inline wide_int
operator<=(float A, wide_float B) {
    return (WideFloatFromFloat(A) <= B);
}

inline wide_int
operator==(wide_float A, wide_float B) {
    wide_int ret;
    ret.V = _mm256_castps_si256(_mm256_cmp_ps(A.V, B.V, _CMP_EQ_OQ));
    return ret;
}
inline wide_int
operator==(wide_float A, float B) {
    return (A == WideFloatFromFloat(B));
}
inline wide_int
operator==(float A, wide_float B) {
    return (WideFloatFromFloat(A) == B);
}

inline wide_int
operator!=(wide_float A, wide_float B) {
    wide_int ret;
    ret.V = _mm256_castps_si256(_mm256_cmp_ps(A.V, B.V, _CMP_NEQ_OQ));
    return ret;
}
inline wide_int
operator!=(wide_float A, float B) {
    return (A != WideFloatFromFloat(B));
}
inline wide_int
operator!=(float A, wide_float B) {
    return (WideFloatFromFloat(A) != B);
}

inline wide_float
operator+(wide_float A, wide_float B) {
    wide_float ret;
    ret.V = _mm256_add_ps(A.V, B.V);
    return ret;
}
inline wide_float
operator+(wide_float A, float B) {
    return A + WideFloatFromFloat(B);
}
inline wide_float
operator+(float A, wide_float B) {
    return WideFloatFromFloat(A) + B;
}

inline wide_float&
operator+=(wide_float& A, wide_float B) {
    A.V = _mm256_add_ps(A.V, B.V);
    return A;
}
inline wide_float&
operator+=(wide_float& A, float B) {
    A += WideFloatFromFloat(B);
    return A;
}

inline wide_float
operator-(wide_float A) {
    wide_float ret;
    ret.V = 0 - A.V;
    return ret;
}

inline wide_float
operator-(wide_float A, wide_float B) {
    wide_float ret;
    ret.V = _mm256_sub_ps(A.V, B.V);
    return ret;
}
inline wide_float
operator-(float A, wide_float B) {
    return WideFloatFromFloat(A) - B;
}
inline wide_float
operator-(wide_float A, float B) {
    return A - WideFloatFromFloat(B);
}

inline wide_float&
operator-=(wide_float& A, wide_float B) {
    A.V = _mm256_sub_ps(A.V, B.V);
    return A;
}
inline wide_float&
operator-=(wide_float& A, float B) {
    A -= WideFloatFromFloat(B);
    return A;
}

inline wide_float
operator*(wide_float A, wide_float B) {
    wide_float ret;
    ret.V = _mm256_mul_ps(A.V, B.V);
    return ret;
}
inline wide_float
operator*(float A, wide_float B) {
    wide_float ret = WideFloatFromFloat(A) * B;
    return ret;
}
inline wide_float
operator*(wide_float A, float B) {
    wide_float ret = A * WideFloatFromFloat(B);
    return ret;
}

inline wide_float&
operator*=(wide_float& A, wide_float B) {
    A.V = _mm256_mul_ps(A.V, B.V);
    return A;
}
inline wide_float&
operator*=(wide_float& A, float B) {
    A *= WideFloatFromFloat(B);
    return A;
}

inline wide_float
operator/(wide_float A, wide_float B) {
    wide_float ret;
    ret.V = _mm256_div_ps(A.V, B.V);
    return ret;
}
inline wide_float
operator/(float A, wide_float B) {
    wide_float ret = WideFloatFromFloat(A) / B;
    return ret;
}
inline wide_float
operator/(wide_float A, float B) {
    wide_float ret = A / WideFloatFromFloat(B);
    return ret;
}

inline wide_float&
operator/=(wide_float& A, wide_float B) {
    A.V = _mm256_div_ps(A.V, B.V);
    return A;
}
inline wide_float&
operator/=(wide_float& A, float B) {
    A /= WideFloatFromFloat(B);
    return A;
}

inline wide_float
operator|(wide_float A, wide_float B) {
    wide_float ret;
    ret.V = _mm256_or_ps(A.V, B.V);
    return ret;
}
inline wide_float
operator|(float A, wide_float B) {
    wide_float ret = WideFloatFromFloat(A) | B;
    return ret;
}
inline wide_float
operator|(wide_float A, float B) {
    wide_float ret = A | WideFloatFromFloat(B);
    return ret;
}
inline wide_float&
operator|=(wide_float& A, wide_float B) {
    A = A | B;
    return A;
}
inline wide_float&
operator|=(wide_float& A, float B) {
    A |= WideFloatFromFloat(B);
    return A;
}

inline wide_float
operator&(wide_float A, wide_float B) {
    wide_float ret;
    ret.V = _mm256_and_ps(A.V, B.V);
    return ret;
}
inline wide_float
operator&(float A, wide_float B) {
    wide_float ret = WideFloatFromFloat(A) & B;
    return ret;
}
inline wide_float
operator&(wide_float A, float B) {
    wide_float ret = A & WideFloatFromFloat(B);
    return ret;
}

// MIXED
inline wide_float
operator-(wide_float A, wide_int B) {
    wide_float ret;
    ret.V = A.V - _mm256_cvtepi32_ps(B.V);
    return ret;
}
inline wide_float
operator-(wide_int A, wide_float B) {
    wide_float ret;
    ret.V = _mm256_cvtepi32_ps(A.V) - B.V;
    return ret;
}

inline wide_float
operator+(wide_float A, wide_int B) {
    wide_float ret;
    ret.V = A.V + _mm256_cvtepi32_ps(B.V);
    return ret;
}
inline wide_float
operator+(wide_int A, wide_float B) {
    wide_float ret;
    ret.V = _mm256_cvtepi32_ps(A.V) + B.V;
    return ret;
}

inline wide_float
operator*(wide_float A, wide_int B) {
    wide_float ret;
    ret.V = A.V * _mm256_cvtepi32_ps(B.V);
    return ret;
}
inline wide_float
operator*(wide_int A, wide_float B) {
    wide_float ret;
    ret.V = _mm256_cvtepi32_ps(A.V) * B.V;
    return ret;
}

inline wide_float
operator/(wide_float A, wide_int B) {
    wide_float ret;
    ret.V = A.V / _mm256_cvtepi32_ps(B.V);
    return ret;
}
inline wide_float
operator/(wide_int A, wide_float B) {
    wide_float ret;
    ret.V = _mm256_cvtepi32_ps(A.V) / B.V;
    return ret;
}

// ---------
// Functions
// ---------
inline wide_int
Min(wide_int A, wide_int B) {
    wide_int ret;
    ret.V = _mm256_min_epi32(A.V, B.V);
    return ret;
}
inline wide_float
Min(wide_float A, wide_float B) {
    wide_float ret;
    A.V = _mm256_min_ps(A.V, B.V);
    return ret;
}

inline wide_int
Max(wide_int A, wide_int B) {
    wide_int ret;
    ret.V = _mm256_max_epi32(A.V, B.V);
    return ret;
}
inline wide_float
Max(wide_float A, wide_float B) {
    wide_float ret;
    ret.V = _mm256_max_ps(A.V, B.V);
    return ret;
}

inline wide_float
Abs(wide_float A) {
    wide_float ret;
    ret.V = _mm256_andnot_ps(WideFloatFromFloat(-0.0f).V, A.V);
    return ret;
}
inline wide_int
Abs(wide_int A) {
    wide_int ret;
    ret.V = _mm256_abs_epi32(A.V);
    return ret;
}

inline wide_int
Square(wide_int A){
    wide_int ret = A*A;
    return ret;
}
inline wide_float
Square(wide_float A){
    wide_float ret = A*A;
    return ret;
}

inline wide_int
WideIndex(int A) {
    wide_int ret;
    ret.V = _mm256_add_epi32(_mm256_set1_epi32(A), _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7));
    return ret;
};

inline void
ConditionalAssign(wide_float* A, wide_int mask, wide_float B){
    A->V = _mm256_blendv_ps(A->V, B.V, _mm256_castsi256_ps(mask.V));
}

inline void
ConditionalAssign(wide_int* A, wide_int mask, wide_int B){
    A->V = _mm256_blendv_epi8(A->V, B.V, mask.V);
}

inline float
HorizontalMax (wide_float A) {
    float ret = *(float*)&A.V;
    for (int i=1; i<8; i++) {
        float v = *((float*)&A.V+i);
        ret = ret>v ? ret : v;
    }
    return ret;
}
inline float
HorizontalMin (wide_float A) {
    float ret = *(float*)&A.V;
    for (int i=1; i<8; i++) {
        float v = *((float*)&A.V+i);
        ret = ret<v ? ret : v;
    }
    return ret;
}
inline int
HorizontalMax (wide_int A) {
    int ret = *(int*)&A.V;
    for (int i=1; i<8; i++) {
        int v = *((int*)&A.V+i);
        ret = ret>v ? ret : v;
    }
    return ret;
}
inline int
HorizontalMin (wide_int A) {
    int ret = *(int*)&A.V;
    for (int i=1; i<8; i++) {
        int v = *((int*)&A.V+i);
        ret = ret<v ? ret : v;
    }
    return ret;
}

// SSE
#elif (LANE_WIDTH==4)
// Single Scalar
#elif (LANE_WIDTH==1)
#else
#endif

#if (LANE_WIDTH!=1)
#endif
