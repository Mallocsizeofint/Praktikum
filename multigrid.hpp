#ifndef MULTIGRID_HPP_
#define MULTIGRID_HPP_

#include <cstddef>
#include "my_vector.hpp"

namespace FDM {
	//!the prolongation operator
	/*
	 \param n number of nodes per row of the finer mesh
	 \param xH data of coarse mesh
	 \param xh data of fine mesh
	*/
	void prolong (const std::size_t n,
		const estd::vector_t<double>& xH,
		estd::vector_t<double>& xh);

	//! the restricion operator
	/*
	 \param n number of rows in finer mesh
	 \param xh data of finer mesh
	 \param xH set to data of coarser mesh
	*/
	void restrict (const std::size_t n,
		const estd::vector_t<double>& xh,
		estd::vector_t<double>& xH);

	//! Solves the problem using the two-grid-method
	/*
	 \param n number of rows of finer mesh
	 \param u the solution vector
	 \param f the RHS of the system
	 \param epsmin minimal accuarcy to be achieved
	 \return number of iterations needed
	*/
	std::size_t TGM (const std::size_t n,
		estd::vector_t<double>& u,
		const estd::vector_t<double>& f,
		const double epsmin);

	//! Solves the problem using the multi-grid-method
	/*
	 \param n number of rows of finer mesh
	 \param u the solution vector
	 \param f the RHS of the system
	 \param epsmin minimal accuarcy to be achieved
	 \param depth the depth of the recursion
	 \param mu number of iterations per level
	 \return number of iterations needed
	*/
	std::size_t MGM (const std::size_t n,
		estd::vector_t<double>& u,
		const estd::vector_t<double>& f,
		const double epsmin,
		const std::size_t depth,
		const std::size_t mu);

	//! Solves the problem using the multi-grid-method, SOR changed for euler
	/*
	 \param n number of rows of finer mesh
	 \param u the solution vector
	 \param f the RHS of the system
	 \param epsmin minimal accuarcy to be achieved
	 \param depth the depth of the recursion
	 \param mu number of iterations per level
	 \param ht stepsize for euler
	 \return number of iterations needed
	*/
	std::size_t MGM_euler (const std::size_t n,
		estd::vector_t<double>& u,
		const estd::vector_t<double>& f,
		const double epsmin,
		const std::size_t depth,
		const std::size_t mu,
		const double ht);


	//! Calculates a good start vector for multi grid
	/*
	 \param n number of rows of mesh
	 \param u the solution vector
	 \param f the RHS of the system
	 \param depth the depth of the recursion
	*/
	void FMGM (const std::size_t n,
		estd::vector_t<double>& u,
		const estd::vector_t<double>& f,
		const std::size_t depth);

}

#endif //Fileguard