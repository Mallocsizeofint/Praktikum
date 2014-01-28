#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cstddef>
#include "my_vector.hpp"

namespace FDM {

//! An abstract matrix class
/*!
 Handles a n x n matrix without explicitly storing it
*/
class Matrix {
private:
  //! matrix is n x n
  const std::size_t n_;
public:
  // constructor
  Matrix() : n_{0} { }
  Matrix(const std::size_t n) : n_{n} { }

  // destructor
  virtual ~Matrix();

  //! Returns size
  std::size_t size() const { return n_; }

  //! Multiplies a vector with the matrix
  /*!
   \param x input vector
   \param y output vector
  */
  virtual void multMatVec(const estd::vector_t<double>& x, estd::vector_t<double>& y) const = 0;       
  
  //! Abbplies a preconditioner to a vector
  virtual void applyPrecond(estd::vector_t<double>&) const { }

};

}//FDM
 

#endif // MATRIX_HPP



