/**
\file vector_routines.hpp
\author Tobias Olbrich <s6toolbr@uni-bonn.de>
\brief Some basic routines for vectors
*/
#include "my_vector.hpp"

namespace estd {

	//! Calculates the dot-product of 2 vectors
	/*!
	  \param x first input-vector
	  \param y second input vector
	  \return the dot-product of x and y
	*/
	template <class T>
	T dot_product (const vector_t<T>& x, const vector_t<T>& y);

	//! Calculates the 2-Norm of a vector
	/*!
	  \param x vector to calculate norm from
	  \return the norm of x
	*/
	template <class T>
	T nrm2 (const vector_t<T>& x);

	//! Calculates the 2-norm of the difference of to vectors
	/*!
	  \param x first input-vector
	  \param y second input vector
	  \return the 2-norm of x-y
	*/
	template <class T>
	T diffnrm2 (const vector_t<T>& x, const vector_t<T>& y); 
}