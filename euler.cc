#include <cassert>
#include <iostream>
#include "euler.hpp"
#include "laplace_matrix_routines.hpp"
#include "multigrid.hpp"
#include "vector_routines.hpp"

namespace FDM {
	using std::size_t;

	//Some helper functions

	//! Does one explicit euler step
	/*
	 \param n number of cols/rows of grid
	 \param u start vector
	 \param f RHS of system
	 \param ht stepsize
	 \return updated vector
	*/
	estd::vector_t<double> expl_euler_step (
		const size_t n,
		const estd::vector_t<double>& u,
		const estd::vector_t<double>& f,
		const double ht);
}

void FDM::EulerExpl (const std::size_t n,
		estd::vector_t<double>& u,
		const estd::vector_t<double>& f,
		const double ht,
		const std::size_t iter) {

	//Check dimensions
	assert (u.size() == n*n && f.size() == n*n);

	//Check for sane stepsize
	assert (ht > 0);

	u=f;

	//Do desired number of steps and write results to VTK-File
	for (size_t i = 0; i < iter; ++i) {
		//u = expl_euler_step(n,u,f,ht);
		estd::vector_t<double> temp (n*n);
		apply_laplace(n,u,temp);
		estd::scale(-1*ht,temp);
		estd::axpy(1,u,temp);
		estd::axpy(ht,f,temp);
		u=temp;
		createVTKFile(n,u.data());
	}

}

void FDM::EulerImpl (const std::size_t n,
		estd::vector_t<double>& u,
		const estd::vector_t<double>& f,
		const double ht,
		const std::size_t iter) {

	//Check dimensions
	assert (u.size() == n*n && f.size() == n*n);

	//Check for sane stepsize
	assert (ht > 0);

	u=f;

	for (size_t i = 0; i < iter; ++i) {
		estd::vector_t<double> rhs = u;
		estd::axpy(ht,f,u);

		MGM_euler(n,rhs,u,1e-4,3,1,ht);
		createVTKFile(n,u.data());
	}

}


estd::vector_t<double> FDM::expl_euler_step (
		const size_t n,
		const estd::vector_t<double>& u,
		const estd::vector_t<double>& f,
		const double ht) {

	//Updated vector
	estd::vector_t<double> ret(n*n);

	//First block
	//Special case first row
    ret[0] = u[0] - ht * (4*u[0] - u[1] - u[n]) + ht * f[0];
    for (size_t i = 1; i < n-1; ++i)
        ret[i] = u[i] - ht * (-u[i-1] + 4*u[i] - u[i+1] - u[i+n]) + ht * f[i];
    //n-th row
    ret[n-1] = u[n-1] - ht * (-u[n-2] + 4*u[n-1] - u[2*n-1]) + ht * f[n-1];

    //Middle blocks
    for(size_t i = 1; i < n - 1; ++i){ 
        ret[i*n] = u[i*n] - ht * (-u[i*n-n] + 4*u[i*n] - u[i*n+1] - u[i*n+n])
        			+ ht * f[i*n];
        for(size_t j = i*n + 1; j < (i+1)*n - 1; ++j)
            ret[j] = u[j] - ht * (-u[j-n] - u[j-1] + 4*u[j] -u[j+1] - u[j+n])
        				+ ht * f[j];
        ret[(i+1)*n-1] = u[(i+1)*n-1] - ht * (-u[i*n-1] - u[(i+1)*n-2] 
        	+ 4*u[(i+1)*n-1] - u[(i+2)*n-1]) + ht * f[(i+1)*n-1];
    }

    //bottom block
    ret[n*n-n] = u[n*n-n] - ht * (-u[n*n-2*n] + 4*u[n*n-n] - u[n*n-n+1]) 
    				+ ht * f[n*n-n];
    for(size_t i = n*n - n + 1; i < n*n - 1; ++i)
        ret[i] = u[i] - ht * (-u[i-n] - u[i-1] + 4*u[i] - u[i+1]) + ht * f[i];
    //Special case last row
    ret[n*n - 1] = u[n*n-1] - ht * (-u[n*n-n-1] - u[n*n-2] + 4*u[n*n-1])
    				+ ht * f[n*n-1];

    return ret;
}