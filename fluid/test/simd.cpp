// g++ -I../src/Math -march=native simd.cpp -o simd && ./simd
#include <assert.h>
#include <iostream>
#include <ostream>
#include "simd.h"

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

int main() {
    // Construct
    if (0){
        wide_int A(1);
        wide_int B = WideIntFromInt(1);
        wide_int C(B);
        wide_int D = C;
    }
    // Basic Arithmetic
    // +
    if (0){
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
    if (0){
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
    if (0){
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
    if (0){
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
    if (0){
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
    {
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
    if (0){
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
    if (0){
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


}
