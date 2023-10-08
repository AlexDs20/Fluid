#include <math.h>
#include <stdexcept>
#include "tensor.hpp"

//------------------------------
//  CONSTRUCTORS
//------------------------------
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
// float& Tensor::operator()(int linear_index) {
//     return _data[linear_index];
// };
// const float& Tensor::operator()(int linear_index) const {
//     return _data[linear_index];
// };
float& Tensor::operator()(int i, int j) {
    return _data[get_linear_index(i, j)];
};
const float& Tensor::operator()(int i, int j) const {
    return _data[get_linear_index(i, j)];
};

//------------------------------
//  MEMBER FUNCTIONS
//------------------------------
float Tensor::max() const {
    // Only take max in the domain, removes the extra edge cell
    float ret = _data[0];
    for (int i=1; i!=shape(-2); ++i)
        for (int j=1; j!=shape(-1); ++j)
            ret = std::max(_data[get_linear_index(i, j)], ret);
    return ret;
};
float Tensor::min() const {
    // Only take min in the domain, removes the extra edge cell
    float ret = _data[0];
    for (int i=1; i!=shape(-2); ++i)
        for (int j=1; j!=shape(-1); ++j)
            ret = std::min(_data[get_linear_index(i, j)], ret);
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

const std::vector<float>& Tensor::data() const {
    return _data;
};

Tensor::size_type inline Tensor::get_linear_index(int i, int j) const {
    return i + j * _shape[0];
}

//------------------------------
//  Friend function
//------------------------------
std::ostream& operator<<(std::ostream& os, const Tensor& T) {
    int size_x = T._shape[T._shape.size()-2];
    int size_y = T._shape[T._shape.size()-1];

    for (int j=size_y-1; j>=0; --j) {
        for (int i=0; i!=size_x; ++i)
            os << T(i, j) << "   ";
        os << std::endl;
    }
    return os;
};


//------------------------------
//------------------------------
//  Domain CONSTRUCTORS
//------------------------------
Domain::Domain(std::vector<int> shape, Cell value): _shape(shape) {
    std::vector<float>::size_type size = 1;

    for (std::vector<int>::const_iterator it=_shape.begin(); it!=_shape.end(); ++it)
        size *= *it;

    _data.assign(size, value);
};

Domain::Domain(std::vector<int> shape, std::vector<Cell> values): _shape(shape) {
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
// float& Domain::operator()(int linear_index) {
//     return _data[linear_index];
// };
// const float& Domain::operator()(int linear_index) const {
//     return _data[linear_index];
// };
Cell& Domain::operator()(int i, int j) {
    return _data[get_linear_index(i, j)];
};
const Cell& Domain::operator()(int i, int j) const {
    return _data[get_linear_index(i, j)];
};

//------------------------------
//  MEMBER FUNCTIONS
//------------------------------
const std::vector<Cell>& Domain::data() const {
    return _data;
};

Domain::size_type inline Domain::get_linear_index(int i, int j) const {
    return i + j * _shape[0];
}

//------------------------------
//  Friend function
//------------------------------
std::ostream& operator<<(std::ostream& os, const Domain& T) {
    int size_x = T._shape[T._shape.size()-2];
    int size_y = T._shape[T._shape.size()-1];

    for (int j=size_y-1; j>=0; --j) {
        for (int i=0; i!=size_x; ++i)
            os << T(i, j).obstacle << "   ";
        os << std::endl;
    }
    return os;
};
