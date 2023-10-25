#pragma once
#include <vector>
#include <ostream>

class Tensor {
    friend std::ostream& operator<<(std::ostream& os, const Tensor& T);

    public:
        Tensor(std::vector<int> shape, float value);
        ~Tensor();

        float& operator()(int i, int j);
        const float& operator()(int i, int j) const;

        float& operator()(int linear_index);
        const float& operator()(int linear_index) const;

        float max() const;
        float min() const;
        // Tensor normalize() const;
        int inline shape(int dim) const {
            if (dim>=0)
                return _shape[dim];
            else
                return _shape[_shape.size()+dim];
        };
        const float* data() const;

    private:
        int get_linear_index(int i, int j) const;
        std::vector<int> _shape;
        float* _data;
        // int block_size=64;
        // int b2;
        // int kimax;
};


// Uses bit field for the domain
struct Cell {       // could default initialize with false with C++20
    bool obstacle: 1;
    bool N: 1;
    bool S: 1;
    bool W: 1;
    bool E: 1;
};

class Domain {
    friend std::ostream& operator<<(std::ostream& os, const Domain& T);

    public:
        typedef std::vector<int>::size_type size_type;

        Domain(std::vector<int> shape, Cell value);
        Domain(std::vector<int> shape, std::vector<Cell> values);

        Cell& operator()(int i, int j);
        const Cell& operator()(int i, int j) const;

        // int& operator()(int linear_index);
        // const int& operator()(int linear_index) const;

        int inline shape(int dim) const {
            if (dim>=0)
                return _shape[dim];
            else
                return _shape[_shape.size()+dim];
        };
        const std::vector<Cell>& data() const;

    private:
        size_type get_linear_index(int i, int j) const;
        std::vector<int> _shape;
        std::vector<Cell> _data;
        // int block_size=64;
        // int b2;
        // int kimax;
};
