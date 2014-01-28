#ifndef LAPLACE_MATRIX_JACOBI_HPP
#define LAPLACE_MATRIX_JACOBI_HPP

#include "laplace_matrix.hpp"

namespace FDM {
class LaplaceMatrixJacobi : public LaplaceMatrix {
public:
  // constructor
  LaplaceMatrixJacobi(const std::size_t size) 
    : LaplaceMatrix(size) { }

  // destructor
  ~LaplaceMatrixJacobi() { }

  // methods
  void applyPrecond(estd::vector_t<double>& x) const;
}; 

}//FDM

#endif //fileguard