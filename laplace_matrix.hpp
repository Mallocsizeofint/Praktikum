#ifndef LAPLACE_MATRIX_HPP
#define LAPLACE_MATRIX_HPP

#include "matrix.hpp"

namespace FDM {

//! Handles a discretized laplace operator without storing its entries
class LaplaceMatrix : public Matrix {
private:
  //!The number of nodes per row
  const std::size_t size_;
public:
  //! Constructor
  /*
   \param size number of nodes per row
  */
  LaplaceMatrix(const std::size_t size) : Matrix(size*size), size_{size} { }

  // destructor
  ~LaplaceMatrix() { }
  std::size_t operator_size() const { return size_; }
  
  // methods
  void multMatVec(const estd::vector_t<double>&, estd::vector_t<double>&) const;
};

}//FDM

#endif //fileguard
