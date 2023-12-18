// g++ -I../src/Math -march=native simd.cpp -o simd && ./simd
#include <assert.h>
#include <iostream>
#include <ostream>
#include "Utils/simd.hpp"

void print(wide_int a) {
    int v[8];
    StoreWideInt(v, a);
    for (int i=0; i<8; i++) {
        std::cout << v[i] << ", ";
    }
    std::cout << std::endl;
}

void print(wide_float a) {
    float v[8];
    StoreWideFloat(v, a);
    for (int i=0; i<8; i++) {
        std::cout << v[i] << ", ";
    }
    std::cout << std::endl;
}


void test_wide_int();
void test_wide_float();
void test_wide_functions();

int main() {
    test_wide_int();
    test_wide_float();
    test_wide_functions();
}

void test_wide_int() {
    // Construct
    if (1){
        wide_int A(1);
        wide_int B = WideIntFromInt(1);
        wide_int C(B);
        wide_int D = C;
    }
    // Basic Arithmetic
    // +
    if (1){
        wide_int A(1);
        wide_int B(2);
        wide_int C;
        C = A + 1;
        std::cout << ">> 2" << std::endl;
        print(C);
        C = 1 + A;
        std::cout << ">> 2" << std::endl;
        print(C);
        A += 1;
        std::cout << ">> 2" << std::endl;
        print(A);
        A += B;
        std::cout << ">> 4" << std::endl;
        print(A);
        C = A + B;
        std::cout << ">> 6" << std::endl;
        print(C);
    }
    // -
    if (1){
        wide_int A(1);
        wide_int B(2);
        wide_int C;
        C = A - 1;
        std::cout << ">> 0" << std::endl;
        print(C);
        C = 1 - A;
        std::cout << ">> 0" << std::endl;
        print(C);
        A -= 1;
        std::cout << ">> 0" << std::endl;
        print(A);
        A -= B;
        std::cout << ">> -2" << std::endl;
        print(A);
        C = A - B;
        std::cout << ">> -4" << std::endl;
        print(C);
    }
    // *
    if (1){
        std::cout << "\n*" << std::endl;
        wide_int A(1);
        wide_int B(2);
        wide_int C;
        C = A * 2;
        std::cout << ">> 2" << std::endl;
        print(C);
        C = 2 * A;
        std::cout << ">> 2" << std::endl;
        print(C);
        A *= 2;
        std::cout << ">> 2" << std::endl;
        print(A);
        A *= B;
        std::cout << ">> 4" << std::endl;
        print(A);
        C = A * B;
        std::cout << ">> 8" << std::endl;
        print(C);
    }
    // <
    if (1){
        std::cout << "\n<" << std::endl;
        int v[8] = {0, 1, 2, 3, 4, 5, 6, 7};
        int w[8] = {3, 3, 3, 3, 3, 3, 3, 3};
        wide_int A(v);
        wide_int B(w);
        wide_int C;

        C = A < B;
        std::cout << ">>\n-1, -1, -1, 0, 0, 0, 0, 0," << std::endl;
        print(C);

        C = A < 3;
        std::cout << ">>\n-1, -1, -1, 0, 0, 0, 0, 0," << std::endl;
        print(C);

        C = 3 < A;
        std::cout << ">>\n0, 0, 0, 0, -1, -1, -1, -1," << std::endl;
        print(C);

        C = A<=B;
        std::cout << ">>\n-1, -1, -1, -1, 0, 0, 0, 0," << std::endl;
        print(C);

        C = A<=3;
        std::cout << ">>\n-1, -1, -1, -1, 0, 0, 0, 0," << std::endl;
        print(C);

        C = 3<=A;
        std::cout << ">>\n0, 0, 0, -1, -1, -1, -1, -1," << std::endl;
        print(C);
    }
    // >
    if (1){
        std::cout << "\n<" << std::endl;
        int v[8] = {0, 1, 2, 3, 4, 5, 6, 7};
        int w[8] = {3, 3, 3, 3, 3, 3, 3, 3};
        wide_int A(v);
        wide_int B(w);
        wide_int C;

        C = A > B;
        std::cout << ">>\n0, 0, 0, 0, -1, -1, -1, -1," << std::endl;
        print(C);

        C = A > 3;
        std::cout << ">>\n0, 0, 0, 0, -1, -1, -1, -1," << std::endl;
        print(C);

        C = 3 > A;
        std::cout << ">>\n-1, -1, -1, 0, 0, 0, 0, 0," << std::endl;
        print(C);

        C = A>=B;
        std::cout << ">>\n0, 0, 0, -1, -1, -1, -1, -1," << std::endl;
        print(C);

        C = A>=3;
        std::cout << ">>\n0, 0, 0, -1, -1, -1, -1, -1," << std::endl;
        print(C);

        C = 3>=A;
        std::cout << ">>\n-1, -1, -1, -1, 0, 0, 0, 0," << std::endl;
        print(C);
    }
    // ==
    if (1){
        int v[8] = {1, 2, 1, 2, 1, 2, 1, 2};
        int w[8] = {2, 2, 2, 2, 2, 2, 2, 2};
        wide_int A(v);
        wide_int B(w);
        wide_int C;

        C = (A==B);
        std::cout << ">>\n0, -1, 0, -1, 0, -1, 0, -1," << std::endl;
        print(C);

        C = (A==A);
        std::cout << ">>\n-1, -1, -1, -1, -1, -1, -1, -1" << std::endl;
        print(C);

        C = (A==2);
        std::cout << ">>\n0, -1, 0, -1, 0, -1, 0, -1," << std::endl;
        print(C);

        C = (2==A);
        std::cout << ">>\n0, -1, 0, -1, 0, -1, 0, -1," << std::endl;
        print(C);
    }
    // !=
    if (1){
        int v[8] = {1, 2, 1, 2, 1, 2, 1, 2};
        int w[8] = {2, 2, 2, 2, 2, 2, 2, 2};
        wide_int A(v);
        wide_int B(w);
        wide_int C;

        C = (A!=B);
        std::cout << ">>\n-1, 0, -1, 0, -1, 0, -1, 0," << std::endl;
        print(C);

        C = (A!=A);
        std::cout << ">>\n0, 0, 0, 0, 0, 0, 0, 0," << std::endl;
        print(C);

        C = (A!=2);
        std::cout << ">>\n-1, 0, -1, 0, -1, 0, -1, 0," << std::endl;
        print(C);

        C = (2!=A);
        std::cout << ">>\n-1, 0, -1, 0, -1, 0, -1, 0," << std::endl;
        print(C);
    }
    // |
    if (1){
        int v[8] = {1, 2, 1, 2, 1, 2, 1, 2};
        int w[8] = {2, 2, 2, 2, 2, 2, 2, 2};
        wide_int A(v);
        wide_int B(w);
        wide_int C;

        C = A | B;
        std::cout << ">>\n3, 2, 3, 2, 3, 2, 3, 2," << std::endl;
        print(C);

        C = B | 1;
        std::cout << ">>\n3, 3, 3, 3, 3, 3, 3, 3," << std::endl;
        print(C);

        C = 1 | A;
        std::cout << ">>\n1, 3, 1, 3, 1, 3, 1, 3," << std::endl;
        print(C);

        C = B | 2;
        std::cout << ">>\n2, 2, 2, 2, 2, 2, 2, 2," << std::endl;
        print(C);

        C = 2 | B;
        std::cout << ">>\n2, 2, 2, 2, 2, 2, 2, 2," << std::endl;
        print(C);
    }
    // &
    if (1){
        int v[8] = {1, 2, 1, 2, 1, 2, 1, 2};
        int w[8] = {2, 2, 2, 2, 2, 2, 2, 2};
        wide_int A(v);
        wide_int B(w);
        wide_int C;

        C = (A & B);
        std::cout << ">>\n0, 2, 0, 2, 0, 2, 0, 2," << std::endl;
        print(C);

        C = (A & 2);
        std::cout << ">>\n0, 2, 0, 2, 0, 2, 0, 2," << std::endl;
        print(C);

        C = (2 & A);
        std::cout << ">>\n0, 2, 0, 2, 0, 2, 0, 2," << std::endl;
        print(C);

        C = (0xFFFFFFFF & A);
        std::cout << ">>\n";
        print(A);
        print(C);

        C = (0x0 & A);
        std::cout << ">>\n0, 0, 0, 0, 0, 0, 0, 0," << std::endl;
        print(C);
    }
}

void test_wide_float() {
    // Construct
    if (1){
        std::cout << "\nConstruct" << std::endl;
        wide_float A(1);
        wide_float B = WideFloatFromFloat(1);
        wide_float C(B);
        wide_float D = C;
    }
    // Basic Arithmetic
    // +
    if (1){
        std::cout << "\n+" << std::endl;
        wide_float A(1);
        wide_float B(2);
        wide_float C;
        C = A + 1;
        std::cout << ">> 2" << std::endl;
        print(C);
        C = 1 + A;
        std::cout << ">> 2" << std::endl;
        print(C);
        A += 1;
        std::cout << ">> 2" << std::endl;
        print(A);
        A += B;
        std::cout << ">> 4" << std::endl;
        print(A);
        C = A + B;
        std::cout << ">> 6" << std::endl;
        print(C);
    }
    // -
    if (1){
        std::cout << "\n-" << std::endl;
        wide_float A(1);
        wide_float B(2);
        wide_float C;
        C = A - 1;
        std::cout << ">> 0" << std::endl;
        print(C);
        C = 1 - A;
        std::cout << ">> 0" << std::endl;
        print(C);
        A -= 1;
        std::cout << ">> 0" << std::endl;
        print(A);
        A -= B;
        std::cout << ">> -2" << std::endl;
        print(A);
        C = A - B;
        std::cout << ">> -4" << std::endl;
        print(C);
    }
    // *
    if (1){
        std::cout << "\n*" << std::endl;
        wide_float A(1);
        wide_float B(2);
        wide_float C;
        C = A * 2;
        std::cout << ">> 2" << std::endl;
        print(C);
        C = 2 * A;
        std::cout << ">> 2" << std::endl;
        print(C);
        A *= 2;
        std::cout << ">> 2" << std::endl;
        print(A);
        A *= B;
        std::cout << ">> 4" << std::endl;
        print(A);
        C = A * B;
        std::cout << ">> 8" << std::endl;
        print(C);
    }
    // /
    if (1){
        std::cout << "\n/" << std::endl;
        wide_float A(1);
        wide_float B(2);
        wide_float C;
        C = A / 2;
        std::cout << ">> 0.5" << std::endl;
        print(C);
        C = 2 / A;
        std::cout << ">> 2" << std::endl;
        print(C);
        A /= 2;
        std::cout << ">> 0.5" << std::endl;
        print(A);
        A /= B;
        std::cout << ">> 0.25" << std::endl;
        print(A);
        C = A / B;
        std::cout << ">> 0.125" << std::endl;
        print(C);
    }
    // <
    if (1){
        std::cout << "\n<" << std::endl;
        float v[8] = {0, 1, 2, 3, 4, 5, 6, 7};
        float w[8] = {3, 3, 3, 3, 3, 3, 3, 3};
        wide_float A(v);
        wide_float B(w);
        wide_int C;

        C = A < B;
        std::cout << ">>\n-1, -1, -1, 0, 0, 0, 0, 0," << std::endl;
        print(C);

        C = A < 3;
        std::cout << ">>\n-1, -1, -1, 0, 0, 0, 0, 0," << std::endl;
        print(C);

        C = 3 < A;
        std::cout << ">>\n0, 0, 0, 0, -1, -1, -1, -1," << std::endl;
        print(C);

        C = A<=B;
        std::cout << ">>\n-1, -1, -1, -1, 0, 0, 0, 0," << std::endl;
        print(C);

        C = A<=3;
        std::cout << ">>\n-1, -1, -1, -1, 0, 0, 0, 0," << std::endl;
        print(C);

        C = 3<=A;
        std::cout << ">>\n0, 0, 0, -1, -1, -1, -1, -1," << std::endl;
        print(C);
    }
    // >
    if (1){
        std::cout << "\n>" << std::endl;
        float v[8] = {0, 1, 2, 3, 4, 5, 6, 7};
        float w[8] = {3, 3, 3, 3, 3, 3, 3, 3};
        wide_float A(v);
        wide_float B(w);
        wide_int C;

        C = A > B;
        std::cout << ">>\n0, 0, 0, 0, -1, -1, -1, -1," << std::endl;
        print(C);

        C = A > 3;
        std::cout << ">>\n0, 0, 0, 0, -1, -1, -1, -1," << std::endl;
        print(C);

        C = 3 > A;
        std::cout << ">>\n-1, -1, -1, 0, 0, 0, 0, 0," << std::endl;
        print(C);

        C = A>=B;
        std::cout << ">>\n0, 0, 0, -1, -1, -1, -1, -1," << std::endl;
        print(C);

        C = A>=3;
        std::cout << ">>\n0, 0, 0, -1, -1, -1, -1, -1," << std::endl;
        print(C);

        C = 3>=A;
        std::cout << ">>\n-1, -1, -1, -1, 0, 0, 0, 0," << std::endl;
        print(C);
    }
    // ==
    if (1){
        std::cout << "\n==" << std::endl;
        float v[8] = {1, 2, 1, 2, 1, 2, 1, 2};
        float w[8] = {2, 2, 2, 2, 2, 2, 2, 2};
        wide_float A(v);
        wide_float B(w);
        wide_int C;

        C = (A==B);
        std::cout << ">>\n0, -1, 0, -1, 0, -1, 0, -1," << std::endl;
        print(C);

        C = (A==A);
        std::cout << ">>\n-1, -1, -1, -1, -1, -1, -1, -1" << std::endl;
        print(C);

        C = (A==2);
        std::cout << ">>\n0, -1, 0, -1, 0, -1, 0, -1," << std::endl;
        print(C);

        C = (2==A);
        std::cout << ">>\n0, -1, 0, -1, 0, -1, 0, -1," << std::endl;
        print(C);
    }
    // !=
    if (1){
        std::cout << "\n!=" << std::endl;
        float v[8] = {1, 2, 1, 2, 1, 2, 1, 2};
        float w[8] = {2, 2, 2, 2, 2, 2, 2, 2};
        wide_float A(v);
        wide_float B(w);
        wide_int C;

        C = (A!=B);
        std::cout << ">>\n-1, 0, -1, 0, -1, 0, -1, 0," << std::endl;
        print(C);

        C = (A!=A);
        std::cout << ">>\n0, 0, 0, 0, 0, 0, 0, 0," << std::endl;
        print(C);

        C = (A!=2);
        std::cout << ">>\n-1, 0, -1, 0, -1, 0, -1, 0," << std::endl;
        print(C);

        C = (2!=A);
        std::cout << ">>\n-1, 0, -1, 0, -1, 0, -1, 0," << std::endl;
        print(C);
    }
    // |
    if (1){
        std::cout << "\n|" << std::endl;
        float v[8] = {1, 2, 1, 2, 1, 2, 1, 2};
        float w[8] = {2, 2, 2, 2, 2, 2, 2, 2};
        wide_float A(v);
        wide_float B(w);
        wide_float C;

        C = A | B;
        std::cout << ">>\ninf, 2, inf, 2, inf, 2, inf, 2," << std::endl;
        print(C);

        C = B | 1;
        std::cout << ">>\ninf, inf, inf, inf, inf, inf, inf, inf," << std::endl;
        print(C);

        C = 1 | A;
        std::cout << ">>\n1, inf, 1, inf, 1, inf, 1, inf," << std::endl;
        print(C);

        C = B | 2;
        std::cout << ">>\n2, 2, 2, 2, 2, 2, 2, 2," << std::endl;
        print(C);

        C = 2 | B;
        std::cout << ">>\n2, 2, 2, 2, 2, 2, 2, 2," << std::endl;
        print(C);
    }
    // &
    if (1){
        std::cout << "\n&" << std::endl;
        float v[8] = {1, 2, 1, 2, 1, 2, 1, 2};
        float w[8] = {2, 2, 2, 2, 2, 2, 2, 2};
        wide_float A(v);
        wide_float B(w);
        wide_float C;

        C = (A & B);
        std::cout << ">>\n0, 2, 0, 2, 0, 2, 0, 2," << std::endl;
        print(C);

        C = (A & 2);
        std::cout << ">>\n0, 2, 0, 2, 0, 2, 0, 2," << std::endl;
        print(C);

        C = (2 & A);
        std::cout << ">>\n0, 2, 0, 2, 0, 2, 0, 2," << std::endl;
        print(C);

        std::cout << "Following is wrong, I don't know how to set all the exponents bit 'up'" << std::endl;
        C = (0xFFFFFFFF & A);
        std::cout << ">>\n1, 2, 1, 2, 1, 2, 1, 2," << std::endl;
        print(C);

        C = ((float)0x0 & A);
        std::cout << ">>\n0, 0, 0, 0, 0, 0, 0, 0," << std::endl;
        print(C);

        C = ((float)1 & A);
        std::cout << ">>\n1, 0, 1, 0, 1, 0, 1, 0," << std::endl;
        print(C);
    }
}

void test_wide_functions() {
    if (1) {
        int v[8] = {1, 2, 3, 4, 5, 6, 7, 8};
        int w[8] = {8, 7, 6, 5, 4, 3, 2, 1};
        wide_int A(v);
        wide_int B(w);
        wide_int C;
        wide_int D(A);
        wide_int E(-A);

        C = Min(A, B);
        std::cout << "\n1, 2, 3, 4, 4, 3, 2, 1" << std::endl;
        print(C);

        ConditionalAssign(&D, (A<B), B);
        std::cout << std::endl;
        print(C);
        print(D);

        C = Max(A, B);
        std::cout << "\n8, 7, 6, 5, 5, 6, 7, 8" << std::endl;
        print(C);

        D = A;
        ConditionalAssign(&D, (A>B), B);
        std::cout << std::endl;
        print(C);
        print(D);

        C = Square(B);
        std::cout << "\n64, 49, 36, 25, 16, 9, 4, 1" << std::endl;
        print(C);

        C = WideIndex(3);
        std::cout << "\n3, 4, 5, 6, 7, 8, 9, 10" << std::endl;
        print(C);

        std::cout << "\n>>\tABS" << std::endl;
        print(E);
        print(Abs(E));

        std::cout << "\n>>\tHorizontalMax" << std::endl;
        print(E);
        std::cout << HorizontalMax(E) << std::endl;

        std::cout << "\n>>\tHorizontalMin" << std::endl;
        print(E);
        std::cout << HorizontalMin(E) << std::endl;
    }

    if (1) {
        float v[8] = {-1, 2, -3, 4, -5, 6, -7, 8};
        float w[8] = {-8, 7, -6, 5, -4, 3, -2, 1};
        float x[8] = {0, -1, -6, 5, -4, 3, -2, 1};
        wide_float A(v);
        wide_float B(w);
        wide_float C;
        wide_float D(x);

        C = Abs(A);
        print(A);
        print(C);

        std::cout << "\n>>\tABS" << std::endl;
        print(D);
        print(Abs(D));

        std::cout << "\n>>\tMax" << std::endl;
        print(Abs(D));
        print(Max(D, Abs(D)));

        std::cout << "\n>>\tHorizontalMax" << std::endl;
        print(D);
        std::cout << HorizontalMax(D) << std::endl;

        std::cout << "\n>>\tHorizontalMin" << std::endl;
        print(D);
        std::cout << HorizontalMin(D) << std::endl;
    }

    if (1) {
        int v[8] = {1,2,3,4,5,6,7,8};
        int w[8] = {9,10,11,12,13,14,15,16};
        wide_int A(v);
        wide_int B(w);

        print(Rotate<2>(A));
        print(RotateRight<2>(A));
        print(Rotate<-2>(A));
        print(RotateLeft<2>(A));

        print(MakeShiftMask<2>());
        print(MakeRightShiftMask<2>());
        print(MakeShiftMask<-2>());
        print(MakeLeftShiftMask<2>());
        print(Shift<2>(A));
        print(RightShift<2>(A));
        print(Shift<-2>(A));
        print(LeftShift<2>(A));

        print(ShiftWithCarry<2>(A, B));
        print(ShiftWithCarry<-2>(A, B));
        print(ShiftRightWithCarry<2>(A, B));
        print(ShiftLeftWithCarry<2>(A, B));

        print(LoadMaskedPackedWideInt(v, MakeShiftMask<7>()));
        print(LoadMaskedPackedWideInt(v, MakeShiftMask<-7>()));

        float f[8] = {1,2,3,4,5,6,7,8};
        print(LoadMaskedPackedWideFloat(f, MakeShiftMask<7>()));
        print(LoadMaskedPackedWideFloat(f, MakeShiftMask<-7>()));
    }
}
