#ifndef LAPLACE_MATRIX_SSOR_HPP
#define LAPLACE_MATRIX_SSOR_HPP

#include <cmath>
#include "laplace_matrix.hpp"

namespace FDM {

constexpr double PI {3.141592653589793238462643383279502884197169399375105820974944};


class LaplaceMatrixSSOR : public LaplaceMatrix {

	public:

	const double omega {2/(1+std::sqrt(1
		-std::cos(PI/operator_size())*std::cos(PI/operator_size())))};
  // constructor
  LaplaceMatrixSSOR(const std::size_t size) 
    : LaplaceMatrix(size) { }

  // destructor
  ~LaplaceMatrixSSOR() { }

  // methods
  void applyPrecond(estd::vector_t<double>& x) const;
};

} //FDM


#endif //fileguard