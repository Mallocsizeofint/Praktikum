#ifndef LAPLACE_MATRIX_ROUTINES_H
#define LAPLACE_MATRIX_ROUTINES_H


#include <cstddef>
#include <functional>
#include <utility>
#include "my_vector.hpp"


//! contains a finite-differences library
namespace FDM {


//! Applies the Laplace-Operator to a vector
/*!
 \param n number of nodes per row, Laplace-Operator is n x n
 \param x the vector to be applied to the Laplace-Operator
 \param y the vector the result is saved to
*/
 void apply_laplace(const std::size_t n,
        const estd::vector_t<double>& x, 
        estd::vector_t<double>& y);

//! Initialises the RHS for the discretation
/*!
 \param n number of nodes per row
 \param f vector to save result
 \param g function of type \c double(double,double). Returns the boundary values
*/
 void init_rhs(const std::size_t n,
                estd::vector_t<double>& f,
                const std::function<double(double,double)>& g);

//! Solves the linear equation system using the Gauss-Seidel-Method
/*!
 \param n the system must be n*n x n*n
 \param f the RHS of the system
 \param u start/solution vector
 \param epsmin the acuracy to achieve
 \param itmax maximum number of iterations to do
 \return a \c std::pair (acuracy reached, iterations needed)
*/
std::pair<double,size_t> gauss_seidel (
        const size_t n,             //The system must be n*nxn*n
        const estd::vector_t<double>& f,  //RHS of the system
        estd::vector_t<double>& u,        //The solution vector
        const double epsmin,        //The best acuracy to achieve
        const size_t itmax);       //The maximum number of iterations



}//FDM

//! Creates a VTK-data-file from array
/*!
 \param n length of vector
 \param u the data-vector
*/
extern void createVTKFile(const unsigned int n,
        const double * const u);

#endif //fileguard
