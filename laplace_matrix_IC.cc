#include <cassert>
#include <cmath>
#include "laplace_matrix_IC.hpp"

namespace FDM {
	//And even more helper
	using std::size_t;

	//initializes the precondition on different types of blocks
	void init_top_block (estd::vector_t<double>& diag_,
		estd::vector_t<double>& subdiag_,
		estd::vector_t<double>& subsubdiag_,
		const LaplaceMatrixIC& L);
	void init_middle_blocks (estd::vector_t<double>& diag_,
		estd::vector_t<double>& subdiag_,
		estd::vector_t<double>& subsubdiag_,
		const LaplaceMatrixIC& L);
	void init_bottom_block (estd::vector_t<double>& diag_,
		estd::vector_t<double>& subdiag_,
		estd::vector_t<double>& subsubdiag_,
		const LaplaceMatrixIC& L);

	
}

void FDM::LaplaceMatrixIC::initPrecond () const {
	//Set the zero entries of subdiag_ to zero, they will stay that way.
	//The extra n-1 doubles to safe won't hurt memory or runtime
	for (size_t i = 1; i < operator_size(); ++i)
		subdiag_[i*operator_size()-1] = 0;
	//Now calculate the IC-decomposition
	init_top_block(diag_,subdiag_,subsubdiag_,*this);
	init_middle_blocks(diag_,subdiag_,subsubdiag_,*this);
	init_bottom_block(diag_,subdiag_,subsubdiag_,*this);
}

void FDM::LaplaceMatrixIC::applyPrecond (estd::vector_t<double>& x) const {
	assert( x.size() == size());
	apply_upper_cholesky(x);
	apply_lower_cholesky(x);
}

void FDM::init_top_block (estd::vector_t<double>& diag_,
	estd::vector_t<double>& subdiag_,
	estd::vector_t<double>& subsubdiag_,
	const LaplaceMatrixIC& L) {

	//First entry:
	diag_[0] = std::sqrt(diag_[0]);
	subdiag_[0] /= diag_[0];
	subsubdiag_[0] /= diag_[0];
	//Rest of first block
	for (size_t i = 1; i < L.operator_size() ; ++i) {
		//new diagonal
		diag_[i] -= subdiag_[i-1] * subdiag_[i-1];

		//i-1 only column with 2 non-zero entries
		//special diag update
		diag_[i] -= subdiag_[i-1]*subsubdiag_[i-1];
		diag_[i-1] -= subdiag_[i-1]*subsubdiag_[i-1];

		//Rest of the update
		diag_[i] = std::sqrt(diag_[i]);
		subdiag_[i] /= diag_[i];
		subsubdiag_[i] /= diag_[i];
	}	
}

void FDM::init_middle_blocks (estd::vector_t<double>& diag_,
		estd::vector_t<double>& subdiag_,
		estd::vector_t<double>& subsubdiag_,
		const LaplaceMatrixIC& L) {

	//For all middle blocks
	for (size_t i = 1; i < L.operator_size() - 1; ++i) {
		//special case first entry of block, no -1 left of diagonal
		diag_[i*L.operator_size()] -= subsubdiag_[(i-1)*L.operator_size()]
			* subsubdiag_[(i-1)*L.operator_size()];
		//No column with more than one entry below i left of that, update column
		diag_[i*L.operator_size()] = std::sqrt(diag_[i*L.operator_size()]);
		subdiag_[i*L.operator_size()] /= diag_[i*L.operator_size()];
		subsubdiag_[i*L.operator_size()] /= diag_[i*L.operator_size()];

		//Now for all the other columns of the block
		for (size_t j = i*L.operator_size() + 1; j < (i+1)*L.operator_size(); ++j) {
			//Update diagonal element
			diag_[j] -= (subsubdiag_[j-L.operator_size()] * subsubdiag_[j-L.operator_size()]
				+ subdiag_[j-1] * subdiag_[j-1]);

			//Still j-1-th column the only one with multiple entries
			//procede with special diag update
			diag_[j] -= subdiag_[j-1]*subsubdiag_[j-1];
			diag_[j-1] -= subdiag_[j-1]*subsubdiag_[j-1];

			//Rest of the update
			diag_[j] = std::sqrt(diag_[j]);
			subdiag_[j] /= diag_[j];
			subsubdiag_[j] /= diag_[j];
		}

	}
}


void FDM::init_bottom_block (estd::vector_t<double>& diag_,
		estd::vector_t<double>& subdiag_,
		estd::vector_t<double>& subsubdiag_,
		const LaplaceMatrixIC& L) {

	//First entry of last block
	diag_[L.size()-L.operator_size()] -= 
		subsubdiag_[L.size()-2*L.operator_size()]*
		subsubdiag_[L.size()-2*L.operator_size()];
	//Like above nothing special to do for first entry
	diag_[L.size()-L.operator_size()] = 
		std::sqrt(diag_[L.size()-L.operator_size()]);
	subdiag_[L.size()-L.operator_size()] /=
		diag_[L.size()-L.operator_size()];

	//The other entries of the block, except last one
	for (size_t i = L.size() - L.operator_size() + 1; i < L.size()-1; ++i) {
		//Update diagonal element
		diag_[i] -= (subsubdiag_[i-L.operator_size()] * subsubdiag_[i-L.operator_size()]
			+ subdiag_[i-1] * subdiag_[i-1]);
		
		//Here we don't even get the special diag-update, so we're almost done
		diag_[i] = std::sqrt(diag_[i]);
			subdiag_[i] /= diag_[i];
	}
	//And the last entry, no subdiag to update:
	diag_[L.size()-1] -= (subdiag_[L.size()-2] * subdiag_[L.size()-2]
		+ subsubdiag_[L.size()-L.operator_size()-1] * 
			subsubdiag_[L.size()-L.operator_size()-1]);
	diag_[L.size()-1] = std::sqrt(diag_[L.size()-1]);
}

void FDM::LaplaceMatrixIC::apply_upper_cholesky(estd::vector_t<double>& x) const {

	// bottom entry
	x[size()-1] /= diag_[size()-1];

	//bottom block, subdiag_ is zero where it should be
	for (size_t i = size() - 2; i >= size() - operator_size(); --i)
		x[i] = (x[i] - subdiag_[i]*x[i+1]) / diag_[i];

	//the rest
	for (size_t i = size() - operator_size() - 1; i > 0; --i)
		x[i] = (x[i] - subdiag_[i]*x[i+1] 
			- subsubdiag_[i]*x[i+operator_size()]) / diag_[i];
	//last entrie because reverse loop
	x[0] = (x[0] - subdiag_[0] * x[1] 
		- subsubdiag_[0]*x[operator_size()] ) / diag_[0]; 
}

void FDM::LaplaceMatrixIC::apply_lower_cholesky (
	estd::vector_t<double>& x) const {

	//first entry
	x[0] /= diag_[0];

	//top block
	for (size_t i = 1; i < operator_size(); ++i)
		x[i] = (x[i] - subdiag_[i-1]*x[i-1]) / diag_[i];

	//the rest
	for (size_t i = operator_size(); i < size(); ++i)
		x[i] = (x[i] - subdiag_[i-1]*x[i-1]
			- subsubdiag_[i-operator_size()] * x[i - operator_size()]) / diag_[i];
}