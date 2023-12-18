#include <immintrin.h>

#define LANE_WIDTH 8

// AVX512
#if (LANE_WIDTH==16)


//======
// AVX2
//======
#elif (LANE_WIDTH==8)

// https://www.youtube.com/watch?v=qejTqnxQRcw
// https://github.com/CppCon/CppCon2020/blob/main/Presentations/adventures_in_simd_thinking_part_1/adventures_in_simd_thinking_part_1__bob_steagall__cppcon_2020.pdf

// typedef __m256 wfloat;
// typedef __m256i wint;

struct wide_float {
    __m256 V;
    wide_float(){};
    wide_float(float);
    wide_float(float*);
    wide_float(__m256);
    wide_float& operator=(float);
};
struct wide_int {
    __m256i V;
    wide_int(){};
    wide_int(int);
    wide_int(int*);
    wide_int(__m256i);
    wide_int& operator=(int);
};


//-----------
// MAKE WIDE
inline wide_float
WideFloatFromFloat(float A){
    return _mm256_set1_ps(A);
}
inline wide_int
WideIntFromInt(int A){
    return _mm256_set1_epi32(A);
}
inline wide_float
WideFloatFromWideInt(wide_int A) {
    return _mm256_cvtepi32_ps(A.V);
}
inline wide_float
WideFloatFromInt(int A) {
    return WideFloatFromFloat((float) A);
}

//--------------
// CONSTRUCTORS
wide_int::wide_int(int A) {
    this->V = WideIntFromInt(A).V;
}
wide_int::wide_int(int* A) {
    this->V = _mm256_loadu_si256((__m256i*)A);
}
wide_int::wide_int(__m256i A) {
    this->V = A;
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
wide_float::wide_float(__m256 A) {
    this->V = A;
}
wide_float& wide_float::operator=(float A) {
    this->V = WideFloatFromFloat(A).V;
    return *this;
}

//--------------
// LOAD - STORE
inline wide_int
LoadInts(int a, int b, int c, int d, int e, int f, int g, int h){
    return _mm256_setr_epi32(a,b,c,d,e,f,g,h);
}
inline wide_float
LoadFloats(float a, float b, float c, float d, float e, float f, float g, float h){
    return _mm256_setr_ps(a,b,c,d,e,f,g,h);
}
inline wide_int
LoadPackedWideInt(int* A) {
    return _mm256_loadu_si256((__m256i*)A);
}
inline wide_float
LoadPackedWideFloat(float* A) {
    return _mm256_loadu_ps(A);
}
inline wide_int
LoadMaskedPackedWideInt(int* A, wide_int mask) {
    return _mm256_maskload_epi32(A, mask.V);
}
inline wide_float
LoadMaskedPackedWideFloat(float* A, wide_int mask) {
    return _mm256_maskload_ps(A, mask.V);
}
inline void
StoreWideInt(int* dest, wide_int source) {
    _mm256_storeu_si256((__m256i*)dest, source.V);
}
inline void
StoreWideFloat(float* dest, wide_float source) {
    _mm256_storeu_ps(dest, source.V);
}
inline void
StoreMaskedWideFloat(float* dest, wide_int mask, wide_float source){
    // if mask -> then save
    _mm256_maskstore_ps(dest, mask.V, source.V);
}
inline void
StoreMaskedWideInt(int* dest, wide_int mask, wide_int source){
    // if mask -> then save
    _mm256_maskstore_epi32(dest, mask.V, source.V);
}

//----------
// WIDE_INT
inline wide_int
operator==(wide_int A, wide_int B) {
    return _mm256_cmpeq_epi32(A.V, B.V);
}
inline wide_int
operator==(int A, wide_int B) {
    return (WideIntFromInt(A) == B);
}
inline wide_int
operator==(wide_int A, int B) {
    return (A == WideIntFromInt(B));
}

inline wide_int
operator!=(wide_int A, wide_int B) {
    return _mm256_andnot_si256((A==B).V, _mm256_set1_epi32(0xFFFFFFFF));
}
inline wide_int
operator!=(int A, wide_int B) {
    return (WideIntFromInt(A)!=B);
}
inline wide_int
operator!=(wide_int A, int B) {
    return (A!=WideIntFromInt(B));
}

inline wide_int
operator|(wide_int A, wide_int B) {
    return _mm256_or_si256(A.V, B.V);
}
inline wide_int
operator|(int A, wide_int B) {
    return WideIntFromInt(A) | B;
}
inline wide_int
operator|(wide_int A, int B) {
    return A | WideIntFromInt(B);
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
    return _mm256_and_si256(A.V, B.V);
}
inline wide_int
operator&(int A, wide_int B) {
    return WideIntFromInt(A) & B;
}
inline wide_int
operator&(wide_int A, int B) {
    return A & WideIntFromInt(B);
}

inline wide_int
operator>(wide_int A, wide_int B) {
    return _mm256_cmpgt_epi32(A.V, B.V);
}
inline wide_int
operator>(int A, wide_int B) {
    return WideIntFromInt(A)>B;
}
inline wide_int
operator>(wide_int A, int B) {
    return A>WideIntFromInt(B);
}

inline wide_int
operator>=(wide_int A, wide_int B) {
    return (A > B) | (A == B);
}
inline wide_int
operator>=(int A, wide_int B) {
    return (WideIntFromInt(A)>=B);
}
inline wide_int
operator>=(wide_int A, int B) {
    return (A>=WideIntFromInt(B));
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
    return _mm256_add_epi32(A.V, B.V);
}
inline wide_int
operator+(int A, wide_int B) {
    return WideIntFromInt(A) + B;
}
inline wide_int
operator+(wide_int A, int B) {
    return A + WideIntFromInt(B);
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
    return 0 - A.V;
}
inline wide_int
operator-(wide_int A, wide_int B) {
    return _mm256_sub_epi32(A.V, B.V);
}
inline wide_int
operator-(int A, wide_int B) {
    return WideIntFromInt(A)-B;
}
inline wide_int
operator-(wide_int A, int B) {
    return A-WideIntFromInt(B);
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
    return _mm256_mullo_epi32(A.V, B.V);
}
inline wide_int
operator*(int A, wide_int B) {
    return WideIntFromInt(A)*B;
}
inline wide_int
operator*(wide_int A, int B) {
    return A*WideIntFromInt(B);
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
    return _mm256_castps_si256(_mm256_cmp_ps(A.V, B.V, _CMP_GT_OQ));
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
    return _mm256_castps_si256(_mm256_cmp_ps(A.V, B.V, _CMP_GE_OQ));
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
    return _mm256_castps_si256(_mm256_cmp_ps(A.V, B.V, _CMP_LT_OQ));
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
    return _mm256_castps_si256(_mm256_cmp_ps(A.V, B.V, _CMP_LE_OQ));
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
    return _mm256_castps_si256(_mm256_cmp_ps(A.V, B.V, _CMP_EQ_OQ));
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
    return _mm256_castps_si256(_mm256_cmp_ps(A.V, B.V, _CMP_NEQ_OQ));
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
    return _mm256_add_ps(A.V, B.V);
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
    return 0 - A.V;
}

inline wide_float
operator-(wide_float A, wide_float B) {
    return _mm256_sub_ps(A.V, B.V);
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
    return _mm256_mul_ps(A.V, B.V);
}
inline wide_float
operator*(float A, wide_float B) {
    return WideFloatFromFloat(A) * B;
}
inline wide_float
operator*(wide_float A, float B) {
    return A * WideFloatFromFloat(B);
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
    return _mm256_div_ps(A.V, B.V);
}
inline wide_float
operator/(float A, wide_float B) {
    return WideFloatFromFloat(A) / B;
}
inline wide_float
operator/(wide_float A, float B) {
    return A / WideFloatFromFloat(B);
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
    return _mm256_or_ps(A.V, B.V);
}
inline wide_float
operator|(float A, wide_float B) {
    return WideFloatFromFloat(A) | B;
}
inline wide_float
operator|(wide_float A, float B) {
    return A | WideFloatFromFloat(B);
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
    return _mm256_and_ps(A.V, B.V);
}
inline wide_float
operator&(float A, wide_float B) {
    return WideFloatFromFloat(A) & B;
}
inline wide_float
operator&(wide_float A, float B) {
    return A & WideFloatFromFloat(B);
}
inline wide_int
operator^(wide_int A, wide_int B){
    return _mm256_xor_si256(A.V, B.V);
}
inline wide_float
operator^(wide_float A, wide_float B){
    return _mm256_xor_ps(A.V, B.V);
}
inline wide_int
operator~(wide_int A){
    return A^wide_int(0xFFFFFFFF);
};

// MIXED
inline wide_float
operator-(wide_float A, wide_int B) {
    return A.V - _mm256_cvtepi32_ps(B.V);
}
inline wide_float
operator-(wide_int A, wide_float B) {
    return _mm256_cvtepi32_ps(A.V) - B.V;
}

inline wide_float
operator+(wide_float A, wide_int B) {
    return A.V + _mm256_cvtepi32_ps(B.V);
}
inline wide_float
operator+(wide_int A, wide_float B) {
    return _mm256_cvtepi32_ps(A.V) + B.V;
}

inline wide_float
operator*(wide_float A, wide_int B) {
    return A.V * _mm256_cvtepi32_ps(B.V);
}
inline wide_float
operator*(wide_int A, wide_float B) {
    return _mm256_cvtepi32_ps(A.V) * B.V;
}

inline wide_float
operator/(wide_float A, wide_int B) {
    return A.V / _mm256_cvtepi32_ps(B.V);
}
inline wide_float
operator/(wide_int A, wide_float B) {
    return _mm256_cvtepi32_ps(A.V) / B.V;
}

// ---------
// Functions
// ---------
inline wide_int
Min(wide_int A, wide_int B) {
    return _mm256_min_epi32(A.V, B.V);
}
inline wide_float
Min(wide_float A, wide_float B) {
    return _mm256_min_ps(A.V, B.V);
}

inline wide_int
Max(wide_int A, wide_int B) {
    return _mm256_max_epi32(A.V, B.V);
}
inline wide_float
Max(wide_float A, wide_float B) {
    return _mm256_max_ps(A.V, B.V);
}

inline wide_float
Abs(wide_float A) {
    return _mm256_andnot_ps(WideFloatFromFloat(-0.0f).V, A.V);
}
inline wide_int
Abs(wide_int A) {
    return _mm256_abs_epi32(A.V);
}

inline wide_float
Fmadd(wide_float a, wide_float b, wide_float c) {
    return _mm256_fmadd_ps(a.V, b.V, c.V);
}

inline wide_int
Square(wide_int A){
    return A*A;
}
inline wide_float
Square(wide_float A){
    return A*A;
}

inline wide_int
WideIndex(int A) {
    return _mm256_add_epi32(_mm256_set1_epi32(A), _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7));
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

template <int R>
inline wide_int
Rotate(wide_int IN) {
    if constexpr ((R % 8) == 0)
        return IN;

    constexpr int V0 = (R > 0) ? (8 - (R % 8)) : -R;
    constexpr int V1 = (V0 + 0) % 8;
    constexpr int V2 = (V0 + 1) % 8;
    constexpr int V3 = (V0 + 2) % 8;
    constexpr int V4 = (V0 + 3) % 8;
    constexpr int V5 = (V0 + 4) % 8;
    constexpr int V6 = (V0 + 5) % 8;
    constexpr int V7 = (V0 + 6) % 8;
    constexpr int V8 = (V0 + 7) % 8;

    return _mm256_permutevar8x32_epi32(IN.V,
            _mm256_setr_epi32(V1,V2,V3,V4,V5,V6,V7,V8));
}
template <int R>
inline wide_int
RotateRight(wide_int IN) {
    static_assert(R>=0);
    return Rotate<R>(IN);
}
template <int R>
inline wide_int
RotateLeft(wide_int IN) {
    static_assert(R>=0);
    return Rotate<-R>(IN);
}

template <int R>
inline wide_int
MakeShiftMask(){
    const __m256i idx = _mm256_setr_epi32(0,1,2,3,4,5,6,7);
    if (R>=0) {
        return _mm256_cmpgt_epi32(_mm256_set1_epi32(R), idx);          // checks for R > idx
    }
    return _mm256_cmpgt_epi32(idx, _mm256_set1_epi32(R+7));
}
template <int R>
inline wide_int
MakeRightShiftMask(){
    static_assert(R>=0);
    return MakeShiftMask<R>();
}
template <int R>
inline wide_int
MakeLeftShiftMask(){
    static_assert(R>=0);
    return MakeShiftMask<-R>();
}

template <int R>
inline wide_int
Shift(wide_int IN) {
    // rotates by R and set 0 instead of wrapping around
    wide_int rotated = Rotate<R>(IN);
    ConditionalAssign(&rotated, MakeShiftMask<R>(), 0);
    return rotated;
}
template <int R>
inline wide_int
LeftShift(wide_int IN) {
    static_assert(R>=0);
    return Shift<-R>(IN);
}
template <int R>
inline wide_int
RightShift(wide_int IN) {
    static_assert(R>=0);
    return Shift<R>(IN);
}

template <int R>
inline wide_int
ShiftWithCarry(wide_int A, wide_int B) {
    wide_int ret = Rotate<R>(A);
    ConditionalAssign(&ret, MakeShiftMask<R>(), Rotate<R>(B));
    return ret;
}
template <int R>
inline wide_int
ShiftLeftWithCarry(wide_int A, wide_int B) {            // Returns A with B shifted in it
    static_assert(R>=0);
    return ShiftWithCarry<-R>(A, B);
}
template <int R>
inline wide_int
ShiftRightWithCarry(wide_int A, wide_int B) {
    static_assert(R>=0);
    return ShiftWithCarry<R>(B, A);
}

// SSE
#elif (LANE_WIDTH==4)
// Single Scalar
#elif (LANE_WIDTH==1)
#else
#endif

#if (LANE_WIDTH!=1)
#endif
