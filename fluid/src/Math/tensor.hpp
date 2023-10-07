#pragma once
#include <iostream>
#include <ostream>
#include <vector>
#include <bitset>

class Tensor {
    friend std::ostream& operator<<(std::ostream& os, const Tensor& T);

    public:
        typedef std::vector<float>::size_type size_type;

        Tensor(): _shape({0}) {};
        Tensor(std::vector<int> shape);
        Tensor(std::vector<int> shape, float value);
        Tensor(std::vector<int> shape, std::vector<float> values);

        float& operator()(const std::vector<int>& idx);
        const float& operator()(const std::vector<int>& idx) const;

        float max() const;
        float amax() const;
        float min() const;
        Tensor normalize() const;
        int inline shape(int dim) const {
            if (dim>=0)
                return _shape[dim];
            else
                return _shape[_shape.size()+dim];
        };
        std::vector<float> data();
        const std::vector<float>& data() const;
        int imax() const;
        int jmax() const;

    private:
        size_type inline get_linear_index(std::vector<int> idx) const;
        std::vector<int> _shape;
        std::vector<float> _data;
};


class TensorBit {
    friend std::ostream& operator<<(std::ostream& os, const TensorBit& T);

    public:
        typedef std::vector<std::bitset<5>>::size_type size_type;

        TensorBit(): _shape({0}) {};
        TensorBit(std::vector<int> shape);
        TensorBit(std::vector<int> shape, std::bitset<5> value);

        std::bitset<5>& operator()(const std::vector<int>& idx);
        const std::bitset<5>& operator()(const std::vector<int>& idx) const;

        int inline shape(int dim) const {
            if (dim>=0)
                return _shape[dim];
            else
                return _shape[_shape.size()+dim];
        };
        std::vector<std::bitset<5>> data();
        const std::vector<std::bitset<5>>& data() const;
        int imax() const;
        int jmax() const;

    private:
        size_type inline get_linear_index(std::vector<int> idx) const;
        std::vector<int> _shape;
        std::vector<std::bitset<5>> _data;
};
