/**
\file my_vector.hpp
\author Tobias Olbrich <yourmomisamilf@uni-bonn.de>
\brief Contains a macro triggered range check vector, oh yeah!
*/
#ifndef MY_VECTOR_H_
#define MY_VECTOR_H_
#include <vector> 



#if defined(ENABLE_RANGE_CHECK)
#include <cstddef> //for size_t
namespace estd {


//! If ENABLE_RANGE_CHECK is defined, then \c vector_t is a rangechecked
//! variant of std::vector, else, it is plain std::vector
template<class T>
class vector_t : public std::vector<T> {
    public:
        using vector<T>::vector;

        T& operator[] (std::size_t i) {
            return vector<T>::at(i);
        }

        const T& operator[] (std::size_t i) const {
            return vector<T>::at(i);
        }
};
}//estd

#else

namespace estd {

//! vector_t is only a \c std::vector if ENABLE_RANGE_CHECK is not defined or true.
template<class T>
using vector_t = std::vector<T>;

#endif //Rangecheck?
}//estd
#endif //fileguard
