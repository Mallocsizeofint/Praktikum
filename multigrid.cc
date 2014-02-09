#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream> //for debug
#include "laplace_matrix.hpp"
#include "laplace_matrix_SSOR.hpp"
#include "laplace_matrix_routines.hpp"
#include "multigrid.hpp"
#include "vector_routines.hpp"

namespace FDM {
	using std::size_t;

	//More helpers

	//For optimal omega
	//constexpr double PI = 3.141592653589793238462643383279502884197169399375105820974944;

	
	//returns the correspondenting index in the node-vector to index (i,j)
	size_t mesh_index (const size_t n, const size_t i, const size_t j);

	//! Does the prolongation everywhere but the last column and row
	/*
	 \param n number of rows in fine mesh
	 \param nc number of rows in coarse mesh
	 \param xH data of coarse mesh
	 \param xh to be set to data of fine mesh
	*/
	void prolong_no_special_case (const size_t n,
		const size_t nc,
		const estd::vector_t<double>& xH,
		estd::vector_t<double>& xh);

	//! Does the prolongation on the last column and row
	/*
	 \param n number of rows in fine mesh
	 \param nc number of rows in coarse mesh
	 \param xH data of coarse mesh
	 \param xh to be set to data of fine mesh
	*/
	void prolong_last_column_and_row (const size_t n,
		const size_t nc,
		const estd::vector_t<double>& xH,
		estd::vector_t<double>& xh);

	//! Does the restriction for inner nodes
	/*
	 \param n number of rows in fine mesh
	 \param nc number of rows in coarse mesh
	 \param xh data data of fine mesh
	 \param xH set to data of coarse mesh
	*/
	void restrict_inner (const size_t n,
		const size_t nc,
		const estd::vector_t<double>& xh,
		estd::vector_t<double>& xH); 

	//! Does the restrition on the border by injection
	/*
	 \param n number of rows in fine mesh
	 \param nc number of rows in coarse mesh
	 \param xh data data of fine mesh
	 \param xH set to data of coarse mesh
	*/
	void restrict_border (const size_t n,
		const size_t nc,
		const estd::vector_t<double>& xh,
		estd::vector_t<double>& xH);

	//! Does mu multigrid steps with depth depth
	void mgm_step (const size_t n,
		estd::vector_t<double>& u,
		const estd::vector_t<double>& f,
		const size_t depth,
		const size_t mu);

	//! Does mu multigrid steps with depth depth for explicit euler with stepsize ht
	void mgm_step_euler (const size_t n,
		estd::vector_t<double>& u,
		const estd::vector_t<double>& f,
		const size_t depth,
		const size_t mu,
		const double ht);
	
}

size_t FDM::mesh_index (const size_t n, const size_t i, const size_t j) {
	return n*i + j;
}

void FDM::prolong_no_special_case (const size_t n,
		const size_t nc,
		const estd::vector_t<double>& xH,
		estd::vector_t<double>& xh) {

	for (size_t i = 0; i < nc - 1; ++i) {
		for (size_t j = 0; j < nc - 1; ++j) {

			//The nodes of the coarse mesh
			xh[mesh_index(n,2*i,2*j)] = xH[mesh_index(nc,i,j)];

			//The  fine node below this
			xh[mesh_index(n,2*i+1,2*j)] = .5 * (xH[mesh_index(nc,i,j)]
											+ xH[mesh_index(nc,i+1,j)]);

			//The fine node right to this
			xh[mesh_index(n,2*i,2*j+1)] = .5 * (xH[mesh_index(nc,i,j)]
											+ xH[mesh_index(nc,i,j+1)]);

			//The fine node to the bottom right
			xh[mesh_index(n,2*i+1,2*j+1)] = .25 * (xH[mesh_index(nc,i,j)]
												+ xH[mesh_index(nc,i+1,j)]
												+ xH[mesh_index(nc,i,j+1)]
												+ xH[mesh_index(nc,i+1,j+1)]);
		}
	}
}

void FDM::prolong_last_column_and_row (const size_t n,
		const size_t nc,
		const estd::vector_t<double>& xH,
		estd::vector_t<double>& xh) {

	//everything but bottom right edge
	for (size_t i = 0; i < nc - 1; ++i) {

		//Right column

		//identity:
		xh[mesh_index(n,2*i,n-1)] = xH[mesh_index(nc,i,nc-1)];

		//Fine node below that
		xh[mesh_index(n,2*i+1,n-1)] = .5 * (xH[mesh_index(nc,i,nc-1)]
										+ xH[mesh_index(nc,i+1,nc-1)]);

		//Bottom row

		//identity:
		xh[mesh_index(n,n-1,2*i)] = xH[mesh_index(nc,nc-1,i)];

		//Fine node right to this
		xh[mesh_index(n,n-1,2*i+1)] = .5 * (xH[mesh_index(nc,nc-1,i)]
										+ xH[mesh_index(nc,nc-1,i+1)]);
	}

	//bottom right node
	xh[mesh_index(n,n-1,n-1)] = xH[mesh_index(nc,nc-1,nc-1)];
}

void FDM::prolong (const size_t n,
	const estd::vector_t<double>& xH,
	estd::vector_t<double>& xh) {

	//number of rows in the coarse mesh
	const size_t nc {(n-1)/2 + 1};

	assert( n*n == xh.size() && nc*nc == xH.size() );

	//prolong everything but the special cases 'right column' and 'bottom row'
	prolong_no_special_case(n,nc,xH,xh);

	//prolong the rest
	prolong_last_column_and_row(n,nc,xH,xh);	
}

void FDM::restrict_inner(const size_t n ,
	const size_t nc,
	const estd::vector_t<double>& xh,
	estd::vector_t<double>& xH) {

	for (size_t i = 1; i < nc - 1; ++i) {
		for (size_t j = 1; j < nc - 1; ++j) {
			xH[mesh_index(nc,i,j)] = (xh[mesh_index(n,2*i-1,2*j-1)]
									+ xh[mesh_index(n,2*i-1,2*j+1)]
									+ xh[mesh_index(n,2*i+1,2*j-1)]
									+ xh[mesh_index(n,2*i+1,2*j+1)]
									+ 2 * (xh[mesh_index(n,2*i,2*j-1)]
										+ xh[mesh_index(n,2*i,2*j+1)]
										+ xh[mesh_index(n,2*i-1,2*j)]
										+ xh[mesh_index(n,2*i+1,2*j)])
									+ 4 * xh[mesh_index(n,2*i,2*j)]) / 16;
		}
	}
}

void FDM::restrict_border (const size_t n,
		const size_t nc,
		const estd::vector_t<double>& xh,
		estd::vector_t<double>& xH) {

	for (size_t i = 0; i < nc; ++i) {

		//top row
		xH[mesh_index(nc,0,i)] = xh[mesh_index(n,0,2*i)];

		//bottom row
		xH[mesh_index(nc,nc-1,i)] = xh[mesh_index(n,n-1,2*i)];

		//left column
		xH[mesh_index(nc,i,0)] = xh[mesh_index(n,2*i,0)];

		//right column
		xH[mesh_index(nc,i,nc-1)] = xh[mesh_index(n,2*i,n-1)];
	}
}

void FDM::restrict (const size_t n,
		const estd::vector_t<double>& xh,
		estd::vector_t<double>& xH) {

	//number of rows in the coarse mesh
	const size_t nc {(n-1)/2 + 1};

	assert( n*n == xh.size() && nc*nc == xH.size() );

	//restrict to inner nodes
	restrict_inner(n,nc,xh,xH);

	//restrict on borders
	restrict_border(n,nc,xh,xH);
}


std::size_t FDM::TGM (const std::size_t n,
		estd::vector_t<double>& u,
		const estd::vector_t<double>& f,
		const double epsmin) {

	//number of rows in the coarse mesh
	const size_t nc {(n-1)/2 + 1};

	assert( n*n == u.size() && n*n == f.size() );

	//Temporary variables
	estd::vector_t<double> fine_residue (n*n);
	estd::vector_t<double> coarse_residue (nc*nc);
	LaplaceMatrixSSOR L {nc};
	estd::vector_t<double> coarse_update (nc*nc);
	estd::vector_t<double> fine_update (n*n);
	const double omega  {.9};

	size_t iterations{0};

	do {
		//smooth on fine mesh:
		SOR(n,f,u,1e-14,5,omega);

		
		//calculate residue
		apply_laplace(n,u,fine_residue);
		estd::scale(-1,fine_residue);
		estd::axpy(1,f,fine_residue);
		
		//get coarse residue
		restrict(n,fine_residue,coarse_residue);

		//Solve on coarse mesh;
		PCG(L,coarse_residue,coarse_update,1e-14,nc*nc);

		//Prolong update
		prolong(n,coarse_update,fine_update);

		//Update u
		estd::axpy(1,fine_update,u);

		//smooth on fine mesh:
		SOR(n,f,u,1e-14,5,omega);

		++iterations;
		
	} while (estd::nrm2(fine_residue) > epsmin);
	return iterations;
}


void FDM::mgm_step (const size_t n,
		estd::vector_t<double>& u,
		const estd::vector_t<double>& f,
		const size_t depth,
		const size_t mu) {

	assert( n*n == u.size() && n*n == f.size() );

	//If coarsest mesh
	if (depth == 0) { 		
		//Solve with SOR with optimal omega
		SOR (n,f,u,1e-20,5,2/(1+std::sqrt(1-std::cos(PI/n)*std::cos(PI/n))));
		return;
	}


	//Else solve with muligrid
	const size_t nc {(n-1)/2 + 1};

	//Temporary variables
	estd::vector_t<double> fine_residue (n*n);
	estd::vector_t<double> coarse_residue (nc*nc);
	estd::vector_t<double> coarse_update (nc*nc);
	estd::vector_t<double> fine_update (n*n);
	const double omega  {.9};

	//Smooth
	SOR(n,f,u,1e-20,5,omega);

	//calculate residue
	apply_laplace(n,u,fine_residue);
	estd::scale(-1,fine_residue);
	estd::axpy(1,f,fine_residue);
	
	//get coarse residue
	restrict(n,fine_residue,coarse_residue);

	//Solve on coarse mesh with mu multigrid steps
	for (size_t i = 0; i < mu; ++i)
		mgm_step(nc,coarse_update,coarse_residue,depth-1,mu);

	//Prolong update
	prolong(n,coarse_update,fine_update);

	//Update u
	estd::axpy(1,fine_update,u);

	//Smooth
	SOR(n,f,u,1e-20,5,omega);
}


void FDM::mgm_step_euler (const size_t n,
		estd::vector_t<double>& u,
		const estd::vector_t<double>& f,
		const size_t depth,
		const size_t mu,
		const double ht) {

	assert( n*n == u.size() && n*n == f.size() );

	//If coarsest mesh
	if (depth == 0) { 		
		//Solve with SOR with optimal omega
		SOR_euler (n,f,u,1e-20,5,1.5,ht);
		return;
	}


	//Else solve with multigrid
	const size_t nc {(n-1)/2 + 1};

	//Temporary variables
	estd::vector_t<double> fine_residue (n*n);
	estd::vector_t<double> coarse_residue (nc*nc);
	estd::vector_t<double> coarse_update (nc*nc);
	estd::vector_t<double> fine_update (n*n);
	const double omega  {.9};

	//Smooth
	SOR_euler(n,f,u,1e-20,5,omega,ht);

	//calculate residue
	apply_laplace(n,u,fine_residue);
	estd::scale(-1*ht,fine_residue);
	estd::axpy(-1,u,fine_residue);
	estd::axpy(1,f,fine_residue);
	
	//get coarse residue
	restrict(n,fine_residue,coarse_residue);

	//Solve on coarse mesh with mu multigrid steps
	for (size_t i = 0; i < mu; ++i)
		mgm_step_euler(nc,coarse_update,coarse_residue,depth-1,mu,ht);

	//Prolong update
	prolong(n,coarse_update,fine_update);

	//Update u
	estd::axpy(1,fine_update,u);

	//Smooth
	SOR_euler(n,f,u,1e-20,5,omega,ht);
}

std::size_t FDM::MGM (const std::size_t n,
	estd::vector_t<double>& u,
	const estd::vector_t<double>& f,
	const double epsmin,
	const std::size_t depth,
	const std::size_t mu) {


	//number of rows in the coarse mesh
	const size_t nc {(n-1)/2 + 1};

	assert( n*n == u.size() && n*n == f.size() );

	//depth == 0 is no multigrid, I consider this missuse
	assert(depth != 0);

	//Temporary variables
	estd::vector_t<double> fine_residue (n*n);
	estd::vector_t<double> coarse_residue (nc*nc);
	estd::vector_t<double> coarse_update (nc*nc);
	estd::vector_t<double> fine_update (n*n);
	const double omega  {.9};//{2/(1+std::sqrt(1-std::cos(PI/n)*std::cos(PI/n)))};

	size_t iterations{0};

	do {

		//First smoothing
		SOR(n,f,u,1e-20,5,omega);

				//calculate residue
		apply_laplace(n,u,fine_residue);
		estd::scale(-1,fine_residue);
		estd::axpy(1,f,fine_residue);
		
		//get coarse residue
		restrict(n,fine_residue,coarse_residue);

		//Solve on coarse mesh;
		std::fill(coarse_update.begin(),coarse_update.end(),0);

		for (size_t i = 0; i < mu; ++i)
			mgm_step(nc,coarse_update,coarse_residue,depth-1,mu);

		//Prolong update
		prolong(n,coarse_update,fine_update);

		//Update u
		estd::axpy(1,fine_update,u);

		//smooth on fine mesh:
		SOR(n,f,u,1e-14,5,omega);

		++iterations;

	} while (estd::nrm2(fine_residue) > epsmin);
	return iterations;
}


std::size_t FDM::MGM_euler (const std::size_t n,
	estd::vector_t<double>& u,
	const estd::vector_t<double>& f,
	const double epsmin,
	const std::size_t depth,
	const std::size_t mu,
	const double ht) {


	//number of rows in the coarse mesh
	const size_t nc {(n-1)/2 + 1};

	assert( n*n == u.size() && n*n == f.size() );

	//depth == 0 is no multigrid, I consider this missuse
	assert(depth != 0);

	//Temporary variables
	estd::vector_t<double> fine_residue (n*n);
	estd::vector_t<double> coarse_residue (nc*nc);
	estd::vector_t<double> coarse_update (nc*nc);
	estd::vector_t<double> fine_update (n*n);
	const double omega  {.9};//{2/(1+std::sqrt(1-std::cos(PI/n)*std::cos(PI/n)))};

	size_t iterations{0};

	do {

		//First smoothing
		SOR_euler(n,f,u,1e-20,5,omega,ht);

				//calculate residue
		apply_laplace(n,u,fine_residue);
		estd::scale(-1*ht,fine_residue);
		estd::axpy(-1,u,fine_residue);
		estd::axpy(1,f,fine_residue);
		
		//get coarse residue
		restrict(n,fine_residue,coarse_residue);

		//Solve on coarse mesh;
		std::fill(coarse_update.begin(),coarse_update.end(),0);

		for (size_t i = 0; i < mu; ++i)
			mgm_step_euler(nc,coarse_update,coarse_residue,depth-1,mu,ht);

		//Prolong update
		prolong(n,coarse_update,fine_update);

		//Update u
		estd::axpy(1,fine_update,u);

		//smooth on fine mesh:
		SOR_euler(n,f,u,1e-14,5,omega,ht);

		++iterations;

	} while (estd::nrm2(fine_residue) > epsmin);
	return iterations;
}

void FDM::FMGM (const std::size_t n,
		estd::vector_t<double>& u,
		const estd::vector_t<double>& f,
		const std::size_t depth) {

	assert( n*n == u.size() && n*n == f.size() );

	// if coarsest grid
	if (depth == 0) {
		//set u = 0
		std::fill(u.begin(),u.end(),0);

		//Solve with V-step
		mgm_step(n,u,f,depth,1);
		return;
	}

	//Else next step down
	//number of rows in the coarse mesh
	const size_t nc {(n-1)/2 + 1};

	estd::vector_t<double> coars_f (nc*nc);
	estd::vector_t<double> coars_u (nc*nc);


	//restrict f
	restrict(n,f,coars_f);

	//Do full multigrid step one level lower
	FMGM(nc,coars_u,coars_f,depth-1);

	//prolong u
	prolong(n,coars_u,u);

	//Do V-cycle
	mgm_step(n,u,f,depth,1);
}