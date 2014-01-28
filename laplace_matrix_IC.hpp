#ifndef LAPLACE_MATRIX_IC_HPP
#define LAPLACE_MATRIX_IC_HPP

#include "laplace_matrix.hpp"
#include <iostream>

namespace FDM {
class LaplaceMatrixIC : public LaplaceMatrix {
private:
	mutable estd::vector_t<double> diag_;
	mutable estd::vector_t<double> subdiag_;
	mutable estd::vector_t<double> subsubdiag_;

	//Applies the upper and the lower cholesky matrix to a vector
	//private member because they must read the decomposition
	void apply_upper_cholesky (estd::vector_t<double>& x) const;
	void apply_lower_cholesky (estd::vector_t<double>& x) const;
public:
  // constructor
  LaplaceMatrixIC(const std::size_t size) 
    : LaplaceMatrix(size), diag_(size*size,4), subdiag_(size*size-1,-1),
    	subsubdiag_(size*size-size,-1)  { initPrecond();}

    void initPrecond() const;

  // destructor
  ~LaplaceMatrixIC() { }

  // methods
  void applyPrecond(estd::vector_t<double>& x) const;
}; 

}//FDM

#endif //fileguard