#pragma once
#include <ostream>
#include <vector>

class Tensor {
    friend std::ostream& operator<<(std::ostream& os, const Tensor& T);

    public:
        typedef std::vector<float>::size_type size_type;

        Tensor(): _shape({0}) {};
        Tensor(std::vector<int> shape);
        Tensor(std::vector<int> shape, float value);

        float& operator()(std::vector<int> idx);
        const float& operator()(std::vector<int> idx) const;

        float max() const;
        float amax() const;
        int shape(int dim) const;
        int imax() const;
        int jmax() const;

    private:
        size_type inline get_linear_index(std::vector<int> idx) const;
        std::vector<int> _shape;
        std::vector<float> _data;
};
