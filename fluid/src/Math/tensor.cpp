#include <math.h>
#include <stdexcept>
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

Tensor::Tensor(std::vector<int> shape, std::vector<float> values): _shape(shape) {
    std::vector<float>::size_type size = 1;

    for (std::vector<int>::const_iterator it=_shape.begin(); it!=_shape.end(); ++it)
        size *= *it;

    if (values.size() != size)
        throw std::domain_error ("Wrong sizes");

    _data = values;
};

//------------------------------
//  OPERATORS element access
//------------------------------
float& Tensor::operator()(const std::vector<int>& idx) {
    return _data[get_linear_index(idx)];
};
const float& Tensor::operator()(const std::vector<int>& idx) const {
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
float Tensor::min() const {
    // Only take max in the domain, removes the extra edge cell
    float ret = 1.0 / 0.0f;
    for (int i=1; i!=shape(-2); ++i)
        for (int j=1; j!=shape(-1); ++j)
            ret = std::min(_data[get_linear_index({i, j})], ret);
    return ret;
};
Tensor Tensor::normalize() const {
    float min = this->min();
    float max = this->max();

    std::vector<float> vec = _data;
    for (std::vector<float>::iterator it=vec.begin(); it!=vec.end(); ++it)
        *it = (*it-min)/(max-min);
    return Tensor(_shape, vec);
};

std::vector<float> Tensor::data() {
    return _data;
};
const std::vector<float>& Tensor::data() const {
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



//--------------------------------------------------
//          TENSORBIT
//--------------------------------------------------
TensorBit::TensorBit(std::vector<int> shape): _shape(shape) {
    TensorBit::size_type size = 1;

    for (std::vector<int>::const_iterator it=_shape.begin(); it!=_shape.end(); ++it)
        size *= *it;

    _data.reserve(size);
};

TensorBit::TensorBit(std::vector<int> shape, std::bitset<5> value): _shape(shape) {
    TensorBit::size_type size = 1;

    for (std::vector<int>::const_iterator it=_shape.begin(); it!=_shape.end(); ++it)
        size *= *it;

    _data.assign(size, value);
};

std::bitset<5>& TensorBit::operator()(const std::vector<int>& idx) {
    return _data[get_linear_index(idx)];
};

const std::bitset<5>& TensorBit::operator()(const std::vector<int>& idx) const {
    return _data[get_linear_index(idx)];
};

std::vector<std::bitset<5>> TensorBit::data() {
    return _data;
};
const std::vector<std::bitset<5>>& TensorBit::data() const {
    return _data;
};

int TensorBit::imax() const {
    return shape(-2)-2;
};
int TensorBit::jmax() const {
    return shape(-1)-2;
};

TensorBit::size_type TensorBit::get_linear_index(std::vector<int> idx) const {
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
};

std::ostream& operator<<(std::ostream& os, const TensorBit& T) {
    int size_x = T._shape[T._shape.size()-2];
    int size_y = T._shape[T._shape.size()-1];

    for (int j=size_y-1; j>=0; --j) {
        for (int i=0; i!=size_x; ++i)
            os << T({i, j}) << "   ";
        os << std::endl;
    }
    return os;
};
