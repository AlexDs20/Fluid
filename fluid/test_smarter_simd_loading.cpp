#include <iostream>
#include "src/Utils/simd.hpp"


void print(int* M, const int imax, const int jmax) {
    for (int j=0; j<jmax; j++) {
        for (int i=0; i<imax; i++) {
            const int index = i + j*imax;
            std::cout << M[index] << ", ";
        }
        std::cout << std::endl;
    }
}

void init(int* M, const int imax, const int jmax) {
    int n=0;
    for (int j=0; j<jmax; j++) {
        for (int i=0; i<imax; i++) {
            const int index = i + j*imax;
            M[index] = n++;
        }
    }
}

void print(wide_int A) {
    int v[8];
    StoreWideInt(v, A);
    for (int i=0; i<LANE_WIDTH; i++) {
        printf("%d, ", v[i]);
    }
    printf("\n");
}

void diff_with_previous(int* M, const int imax, const int jmax) {
    // assert imax > LANE_WIDTH

    wide_int previous;          // -> index-1
    wide_int current;           // -> index
    for (int j=0; j<jmax; j++) {
        for (int i=0; i<imax; i+=LANE_WIDTH) {
            int idx = i + j*imax;

            if (i == 0) {
                current = LoadPackedWideInt(&M[idx]);
                previous = RightShift<1>(current);
            } else {
                previous = current;         // Previous 8 packed
                current = LoadPackedWideInt(&M[idx]);
                previous = ShiftRightWithCarry<1>(previous, current);
            }

            StoreWideInt(&M[idx], current-previous);
        }
    }
}

void diff_with_next(int* M, const int imax, const int jmax) {
    wide_int current;           // -> index
    wide_int next;              // -> index+1
    wide_int next_8;            // -> index+LANE_WIDTH
    wide_int mask;
    wide_int next_mask;
    for (int j=0; j<jmax; j++) {
        for (int i=0; i<imax; i+=LANE_WIDTH) {
            int idx = i + j*imax;

            //          current  next_8
            // -------- -------- -------- -----
            //
            //                   current  next_8 // (which is not 8)
            // -------- -------- -------- -----

            if (i == 0) {
                if (i+LANE_WIDTH-1<imax) {
                    mask = -1;
                    current = LoadPackedWideInt(&M[idx]);
                } else {
                    mask = WideIndex(i) < imax;
                    current = LoadMaskedPackedWideInt(&M[idx], mask);
                }
            } else {
                mask = next_mask;
                current = next_8;
            }

            if (i+2*LANE_WIDTH-1<imax) {
                next_mask = -1;
                next_8 = LoadPackedWideInt(&M[idx+LANE_WIDTH]);
            } else {
                next_mask = WideIndex(i+LANE_WIDTH) < imax;
                next_8 = LoadMaskedPackedWideInt(&M[idx+LANE_WIDTH], next_mask);
            }
            next = ShiftLeftWithCarry<1>(current, next_8);

            StoreMaskedWideInt(&M[idx], mask, next-current);
        }
    }
}


int main(int , char**){

    const int imax = 31;
    const int jmax = 4;

    int M[imax*jmax];

    init(M, imax, jmax);
    // diff_with_previous(M, imax, jmax);
    diff_with_next(M, imax, jmax);
    print(M, imax, jmax);
}
