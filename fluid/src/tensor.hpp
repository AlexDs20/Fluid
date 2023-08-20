#include <ostream>
#include <cstddef>
#include <vector>
#include <math.h>

class Tensor {
    friend std::ostream& operator<<(std::ostream& os, const Tensor& T);
    public:
        typedef std::vector<float>::size_type size_type;

        Tensor(): _shape({0}) {};
        Tensor(std::vector<int> shape): _shape(shape) {
            std::vector<float>::size_type size = 1;

            for (std::vector<int>::const_iterator it=_shape.begin(); it!=_shape.end(); ++it)
                size *= *it;

            _data.reserve(size);
        };

        Tensor(std::vector<int> shape, float value): _shape(shape) {
            std::vector<float>::size_type size = 1;

            for (std::vector<int>::const_iterator it=_shape.begin(); it!=_shape.end(); ++it)
                size *= *it;

            _data.assign(size, value);
        };

        const float& operator()(std::vector<int> idx) const {
            // NO-CHECK
            return _data[get_linear_index(idx)];
        };
        float& operator()(std::vector<int> idx) {
            // NO-CHECK
            return _data[get_linear_index(idx)];
        };

        float max() const {
            float ret = 0.0f;
            for (const auto& val: _data )
                ret = std::max(val, ret);
            return ret;
        };
        float amax() const {
            float ret = 0.0f;
            for (const auto& val: _data )
                ret = std::max(fabs(val), fabs(ret));
            return ret;
        };

        int shape(int dim) const {
            return _shape[dim];
        };

    private:
        size_type inline get_linear_index(std::vector<int> idx) const {
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
        std::vector<int> _shape;
        std::vector<float> _data;
};

std::ostream& operator<<(std::ostream& os, const Tensor& T) {
    // TODO: Fix for all dimensions

    int size_x = T._shape[T._shape.size()-2];
    int size_y = T._shape[T._shape.size()-1];

    for (int i=0; i!=size_x; ++i) {
        for (int j=0; j!=size_y; ++j)
            os << T({i, j}) << "   ";
        os << std::endl;
    }
    return os;
};

// int main() {
//     std::vector<int> shape({5, 3, 4});
//     Tensor T(shape);
//
//     for (int i = 0; i !=shape[0]*shape[1]*shape[2] ; ++i) {
//         T._data[i] = i;
//     }
//
//     for (int i = 0; i != shape[0]; ++i)
//         for (int j = 0; j != shape[1]; ++j)
//             for (int k = 0; k != shape[2]; ++k)
//                 std::cout << "(" << i << ", " << j << ", " << k << "): " << T({i, j, k}) << std::endl;
//     return 0;
// };
