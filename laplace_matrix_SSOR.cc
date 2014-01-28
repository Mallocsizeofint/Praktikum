#include <cassert>
#include <cstddef> //for size_t ffs
#include "laplace_matrix_SSOR.hpp"
#include "vector_routines.hpp"

namespace FDM {
	using std::size_t;
	//Some more helper functions

	//applies (D + \omega*L)^-1 to vector
	void apply_lower_matrix (estd::vector_t<double>& x, const LaplaceMatrixSSOR& L);
	//applies (D + \omega*L)^-1 to vector
	void apply_upper_matrix (estd::vector_t<double>& x, const LaplaceMatrixSSOR& L);
}

void FDM::LaplaceMatrixSSOR::applyPrecond (estd::vector_t<double>& x) const {
	assert(x.size() == size());
	apply_upper_matrix(x,*this);
	estd::scale(.25,x);
	apply_lower_matrix(x,*this);
}

void FDM::apply_lower_matrix (estd::vector_t<double>& x, const LaplaceMatrixSSOR& L) {
	//the top block:
	//first entry:
	x[0] /= 4;

	//rest of the top block:
	for (size_t i = 1; i < L.operator_size(); ++i)
		x[i] = (x[i] + L.omega * x[i-1]) / 4;

	//rest of them blocks:
	for (size_t i = 1; i < L.operator_size(); ++i){
		//First entry of block
		x[i*L.operator_size()] = (x[L.operator_size()*i] + 
			L.omega * (x[(i-1)*L.operator_size()])) / 4;
		for (size_t j = i * L.operator_size() + 1; j < (i+1)*L.operator_size(); ++j)
			x[j] = (x[j] + L.omega * (x[j-1] + x[j-L.operator_size()])) / 4;
	}
}
 
void FDM::apply_upper_matrix (estd::vector_t<double>& x, const LaplaceMatrixSSOR& L) {
	//the bottom block:
	//last entry:
	x[L.size()-1] /= 4;

	//rest of the bottom block:
	for (size_t i = L.size() - 2; i >= L.size() - L.operator_size(); --i)
		x[i] = (x[i] + L.omega * x[i+1]) / 4;

	//rest of them blocks:
	for (size_t i = L.operator_size() - 2; i > 0; --i) {
		//bottom entry of block
		x[(i+1)*L.operator_size()-1] = (x[(i+1)*L.operator_size()-1] + L.omega * (x[(i+2)*L.operator_size()-1])) / 4;
		//rest of block
		for (size_t j = (i+1)*L.operator_size() - 2; j >= i*L.operator_size(); --j)
			x[j] = (x[j] + L.omega * (x[j+1] + x[j + L.operator_size()])) / 4;	
	}
	
	//bottom entry of top block
	x[L.operator_size() - 1] = (x[L.operator_size() - 1] + L.omega * x[2*L.operator_size() - 1]) / 4;
	//rest of first block
	for (size_t i = L.operator_size() - 2; i > 0; --i)
		x[i] = (x[i] + L.omega * (x[i+1] + x[i + L.operator_size()])) / 4;	

	//first entry because reverse loops are evil
	x[0] = (x[0] + L.omega * (x[1] + x[L.operator_size()]));
	
}