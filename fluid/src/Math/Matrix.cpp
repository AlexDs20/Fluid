#include "Math/Matrix.hpp"
#include <iostream>
#include <immintrin.h>


//====================
//  INT MATRIX
//====================
Matrixi::Matrixi(int width, int height, int value): _width(width), _height(height) {
    _size = width * height;

    _data = new int[_size];
    int i;
    for (i=0; i+8<_size; i+=8) {
        _mm256_storeu_si256((__m256i*)&_data[i], _mm256_set1_epi32(value));
    }
    for (; i<_size; i++)
        _data[i] = value;
};

int& Matrixi::operator()(int i, int j) {
    return _data[get_linear_index(i, j)];
};

const int& Matrixi::operator()(int i, int j) const {
    return _data[get_linear_index(i, j)];
};

int& Matrixi::operator[](int linear_index) {
    return _data[linear_index];
};

const int& Matrixi::operator[](int linear_index) const {
    return _data[linear_index];
};

const int* Matrixi::data() const {
    return _data;
};


int inline Matrixi::get_linear_index(int i, int j) const {
    return j * _width + i;
};


std::ostream& operator<<(std::ostream& os, const Matrixi& M) {
    for (int j=0; j<M._height; j++) {
        for (int i=0; i<M._width; i++) {
            os << M._data[M.get_linear_index(i, j)] << "   ";
        }
        os << "\n";
    }
    return os;
};

//====================
//  FLOAT MATRIX
//====================
Matrix::Matrix(int width, int height, float value): _width(width), _height(height) {
    _size = width * height;

    _data = new float[_size];
    int i;
    for (i=0; i+8<_size; i+=8) {
        _mm256_storeu_ps(&_data[i], _mm256_set1_ps(value));
    }
    for (; i<_size; i++)
        _data[i] = value;
};

float& Matrix::operator()(int i, int j) {
    return _data[get_linear_index(i, j)];
};

const float& Matrix::operator()(int i, int j) const {
    return _data[get_linear_index(i, j)];
};

float& Matrix::operator[](int linear_index) {
    return _data[linear_index];
};

const float& Matrix::operator[](int linear_index) const {
    return _data[linear_index];
};

const float* Matrix::data() const {
    return _data;
};


int inline Matrix::get_linear_index(int i, int j) const {
    return j * _width + i;
};


std::ostream& operator<<(std::ostream& os, const Matrix& M) {
    for (int j=0; j<M._height; j++) {
        for (int i=0; i<M._width; i++) {
            os << M._data[M.get_linear_index(i, j)] << "   ";
        }
        os << "\n";
    }
    return os;
};

