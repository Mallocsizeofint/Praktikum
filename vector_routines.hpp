/**
\file vector_routines.hpp
\author Tobias Olbrich <s6toolbr@uni-bonn.de>
\brief Some basic routines for vectors
*/
#include <algorithm>
#include <cmath>
#include <numeric>
#include "my_vector.hpp"

namespace estd {

	//! Calculates the dot-product of 2 vectors
	/*!
	  \param x first input-vector
	  \param y second input vector
	  \return the dot-product of x and y
	*/
	template <class T>
	inline T dot_product (const vector_t<T>& x, const vector_t<T>& y){
		assert (x.size()==y.size());
		return std::inner_product(x.begin(),x.end(),y.begin(),T{0});
	}

	//! Calculates the 2-Norm of a vector
	/*!
	  \param x vector to calculate norm from
	  \return the norm of x
	*/
	template <class T>
	inline T nrm2 (const vector_t<T>& x) {
		return std::sqrt(dot_product(x,x));
	}

	//! Calculates the 2-norm of the difference of to vectors
	/*!
	  \param x first input-vector
	  \param y second input vector
	  \return the 2-norm of x-y
	*/
	template <class T>
	inline T diffnrm2 (const vector_t<T>& x, const vector_t<T>& y){
		assert(x.size() == y.size());
		T ret {0};
		for (std::size_t i = 0; i < x.size(); ++i)
			ret += (x[i]-y[i])*(x[i]-y[i]);
		return std::sqrt(ret);
	}

	//! Multiplies a vector with a scalar
	/*!
	 Multiplication of type S and T must be defined,
	 result must implicitly convert to T
	 \param alpha the constant to multiply with
	 \param x the vector to be modified
	*/
	template <class S, class T>
	inline void scale (const S& alpha, vector_t<T>& x){
		std::transform(x.begin(), x.end(), x.begin(), [&](const T& v){return alpha*v;});
	}

	//! Adds a vector times a scalar to a vector
	/*!
	 Multiplication of type S and T must be defined,
	 result type must have addition with type U 
	 implicitly convertible to U
	 y = alpha * x + y after call
	*/
	template <class S, class T, class U>
	inline void axpy (const S& alpha, const vector_t<T>& x, vector_t<U>& y) {
		assert(x.size() == y.size());
		std::transform(y.begin(), y.end(),x.begin(), y.begin(), [&](const U& v1, const T& v2){return v1 + alpha * v2;});
	}
}