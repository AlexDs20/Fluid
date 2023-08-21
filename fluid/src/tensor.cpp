#include <math.h>
#include "tensor.hpp"

//------------------------------
//  CONSTRUCTORS
//------------------------------
Tensor::Tensor(std::vector<int> shape): _shape(shape) {
    std::vector<float>::size_type size = 1;

    for (std::vector<int>::const_iterator it=_shape.begin(); it!=_shape.end(); ++it)
        size *= *it;

    _data.reserve(size);
};

Tensor::Tensor(std::vector<int> shape, float value): _shape(shape) {
    std::vector<float>::size_type size = 1;

    for (std::vector<int>::const_iterator it=_shape.begin(); it!=_shape.end(); ++it)
        size *= *it;

    _data.assign(size, value);
};

//------------------------------
//  OPERATORS element access
//------------------------------
float& Tensor::operator()(std::vector<int> idx) {
    return _data[get_linear_index(idx)];
};
const float& Tensor::operator()(std::vector<int> idx) const {
    return _data[get_linear_index(idx)];
};

//------------------------------
//  MEMBER FUNCTIONS
//------------------------------
float Tensor::max() const {
    // Only take max in the domain, removes the extra edge cell
    float ret = 0.0f;
    for (int i=1; i!=shape(-2); ++i)
        for (int j=1; j!=shape(-1); ++j)
            ret = std::max(_data[get_linear_index({i, j})], ret);
    return ret;
};
float Tensor::amax() const {
    // Only take max in the domain, removes the extra edge cell
    float ret = 0.0f;
    for (int i=1; i!=shape(-2); ++i)
        for (int j=1; j!=shape(-1); ++j)
            ret = std::max(fabs(_data[get_linear_index({i, j})]), fabs(ret));
    return ret;
};

int Tensor::shape(int dim) const {
    if (dim>=0)
        return _shape[dim];
    else
        return _shape[_shape.size()+dim];
};
std::vector<float> Tensor::data() const {
    return _data;
};
int Tensor::imax() const {
    return shape(-2)-2;
};
int Tensor::jmax() const {
    return shape(-1)-2;
};

Tensor::size_type inline Tensor::get_linear_index(std::vector<int> idx) const {
    size_type linearIdx = 0;
    size_type stride = 1;

    // Packing technique 1: first dimension tightly packed
    for (size_type i=0; i!=_shape.size(); ++i) {
        linearIdx += idx[i] * stride;
        stride *= _shape[i];
    }

    // Packing technique 2: last dimension tightly packed
    // for (size_type i=_shape.size()-1; i!=-1; --i) {
    //     linearIdx += idx[i] * stride;
    //     stride *= _shape[i];
    // }
    return linearIdx;
}

//------------------------------
//  Friend function
//------------------------------
std::ostream& operator<<(std::ostream& os, const Tensor& T) {
    // TODO: Fix for all dimensions

    int size_x = T._shape[T._shape.size()-2];
    int size_y = T._shape[T._shape.size()-1];

    for (int j=size_y-1; j>=0; --j) {
        for (int i=0; i!=size_x; ++i)
            os << T({i, j}) << "   ";
        os << std::endl;
    }
    return os;
};
