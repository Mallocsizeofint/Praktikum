#include "laplace_matrix.hpp"
#include "laplace_matrix_routines.hpp"

void FDM::LaplaceMatrix::multMatVec(
	const estd::vector_t<double>& x,
	estd::vector_t<double>& y) const {
	apply_laplace(operator_size(),x,y);
}