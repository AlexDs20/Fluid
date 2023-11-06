#pragma once
#include <ostream>

class Matrixi {
    friend std::ostream& operator<<(std::ostream& os, const Matrixi& T);

    public:
        Matrixi(int width, int height, int value);
        ~Matrixi() { delete[] _data; };

        int& operator()(int i, int j);
        const int& operator()(int i, int j) const;

        int& operator[](int linear_index);
        const int& operator[](int linear_index) const;

        const int* data() const;

    private:
        int inline get_linear_index(int i, int j) const;
        int _width;
        int _height;
        int _size;
        int* _data;
};


class Matrix {
    friend std::ostream& operator<<(std::ostream& os, const Matrix& M);

    public:
        Matrix(int width, int height, float value);
        ~Matrix() { delete[] _data; };

        float& operator()(int i, int j);
        const float& operator()(int i, int j) const;

        float& operator[](int linear_index);
        const float& operator[](int linear_index) const;

        const float* data() const;

    private:
        int get_linear_index(int i, int j) const;
        int _width;
        int _height;
        int _size;
        float* _data;
};
