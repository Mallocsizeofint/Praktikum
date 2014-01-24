#include <cassert>
#include <cmath>
#include <numeric>
#include "vector_routines.hpp"

template <class T>
T estd::dot_product (const vector_t<T>& x, const vector_t<T>& y) {
	assert (x.size()==y.size());
	return std::inner_product(x.begin(),x.end(),y.begin(),0);
}

template <class T>
T estd::nrm2 (const vector_t<T>& x) {
	return std::sqrt(dot_product(x,x));
}

template <class T>
T estd::diffnrm2 (const vector_t<T>& x, const vector_t<T>& y){
	assert(x.size() == y.size());
	T ret {0};
	for (std::size_t i; i < x.size(); ++i)
		ret += (x[i]-y[i])*(x[i]-y[i]);
	return std::sqrt(ret);
}
