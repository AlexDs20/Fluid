#pragma once
#include <iostream>
#include <ostream>
#include <vector>

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
