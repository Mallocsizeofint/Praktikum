#include "laplace_matrix_jacobi.hpp"
#include "vector_routines.hpp"

void FDM::LaplaceMatrixJacobi::applyPrecond(estd::vector_t<double>& x) const {
	estd::scale(.25,x);
}