#ifndef EULER_HPP_
#define EULER_HPP_

#include <cstddef>
#include "my_vector.hpp"

namespace FDM {

	//! Solves the PDE using explicit euler method
	/*
	 \param n Number of rows/cols
	 \param u start/solution vector
	 \param f RHS of the system
	 \param ht size of time-step
	 \param iter number of time steps
	*/
	void EulerExpl (const std::size_t n,
		estd::vector_t<double>& u,
		const estd::vector_t<double>& f,
		const double ht,
		const std::size_t iter);

	//! Solves the PDE using implicit euler method
	/*
	 \param n Number of rows/cols
	 \param u start/solution vector
	 \param f RHS of the system
	 \param ht size of time-step
	 \param iter number of time steps
	*/
	void EulerImpl (const std::size_t n,
		estd::vector_t<double>& u,
		const estd::vector_t<double>& f,
		const double ht,
		const std::size_t iter);
}


#endif //fileguard