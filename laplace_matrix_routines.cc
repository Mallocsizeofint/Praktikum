
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include "laplace_matrix_routines.hpp"
#include "vector_routines.hpp"


//! Some helper-functions for FDM
namespace FDM {

    using std::size_t;
    

    //! helper 1 for \c FDM::apply_laplace
    /*! 
      Applies the top block of the operator
      \param n number of nodes per row
      \param x vector to be applied
      \param y vector to save result
    */
    void apply_top_block(const size_t n,
       const estd::vector_t<double>& x,
       estd::vector_t<double>& y);

    //! helper 2 for \c FDM::apply_laplace
    /*! 
      Applies the middle blocks of the operator
      \param n number of nodes per row
      \param x vector to be applied
      \param y vector to save result
    */
    void apply_middle_blocks(const size_t n,
       const estd::vector_t<double>& x,
       estd::vector_t<double>& y);

    //! helper 3 for \c FDM::apply_laplace
    /*! 
      Applies the bottom block of the operator
      \param n number of nodes per row
      \param x vector to be applied
      \param y vector to save result
    */
    void apply_bottom_block(const size_t n,
       const estd::vector_t<double>& x,
       estd::vector_t<double>& y);

    //calculates the j-th entry of the gauss-seidel vector
    double gauss_seidel_entry (
        const size_t n,             //The system is n*nxn*n
        const estd::vector_t<double>& f,   //RHS of the system
        const estd::vector_t<double>& u,   //entries 0 to j-1 contain entries
                                    //of the updated vector, j to n-1
                                    //entries of the old one
        const size_t j);           //The entry to calculate

    //calculates the j-th entry of the SOR vector
    double SOR_entry (
        const size_t n,             //The system is n*nxn*n
        const estd::vector_t<double>& f,   //RHS of the system
        const estd::vector_t<double>& u,   //entries 0 to j-1 contain entries
                                    //of the updated vector, j to n-1
                                    //entries of the old one
        const size_t j,              //The entry to calculate 
        const double omega);        //The relaxion parameter


    //returns infnorm of new - old, gets updated vec
    double gauss_seidel_top_block(
        size_t n,
        const estd::vector_t<double>& f,
        estd::vector_t<double>& u);
    double gauss_seidel_middle_blocks(
        size_t n,
        const estd::vector_t<double>& f,
        estd::vector_t<double>& u);
    double gauss_seidel_bottom_block(
        size_t n,
        const estd::vector_t<double>& f,
        estd::vector_t<double>& u);

}

double FDM::gauss_seidel_top_block(
        size_t n,
        const estd::vector_t<double>& f,
        estd::vector_t<double>& u) {
             
    //first line
    double new_entry = (f[0] + u[1] + u[n])/4;
    double infnorm = std::abs(new_entry - u[0]);
    u[0] = new_entry;
    for (size_t i = 1; i < n - 1; ++i) {           //middle line of first block
        new_entry = (f[i] + u[i-1] + u[i+1] + u[i+n]) / 4;
        infnorm = std::max(infnorm,
            std::abs(new_entry - u[i]));
        u[i] = new_entry;
    }
    //last line block 1
    new_entry = (f[n-1] + u[n-2]  + u[2*n-1]) / 4;
    infnorm = std::max(infnorm,
        std::abs(new_entry - u[n-1]));
    u[n-1] = new_entry;
    return infnorm;
}



double FDM::gauss_seidel_bottom_block(
        size_t n,
        const estd::vector_t<double>& f,
        estd::vector_t<double>& u) {
             
    //first line
    double new_entry = (f[n*n-n] + u[n*n-n+1] + u[n*n-2*n])/4;
    double infnorm = std::abs(new_entry - u[n*n-n]);
    u[n*n-n] = new_entry;
    for (size_t i = n*n-n+1; i < n*n - 2; ++i) {           //middle line of first block
        new_entry = (f[i] + u[i-1] + u[i+1] + u[i-n]) / 4;
        infnorm = std::max(infnorm,
            std::abs(new_entry - u[i]));
        u[i] = new_entry;
    }
    //last line block 1
    new_entry = (f[n*n-1] + u[n*n-2]  + u[n*n-n-1]) / 4;
    infnorm = std::max(infnorm,
        std::abs(new_entry - u[n*n-1]));
    u[n*n-1] = new_entry;
    return infnorm;
}
double FDM::gauss_seidel_middle_blocks(
        size_t n,
        const estd::vector_t<double>& f,
        estd::vector_t<double>& u) {
    double infnorm{0};
    //for all middle blocks
    for (size_t j = 1; j < n-1; ++j) {     
        //first line
        double new_entry = (f[j*n] + u[j*n-n] 
            + u[j*n+n] + u[j*n+1])/4;
        infnorm = std::max(infnorm,
            std::abs(new_entry - u[j*n]));
        u[j*n] = new_entry;
        for (size_t i = j*n+1; i < (j+1)*n - 2; ++i) {           //middle line of first block
            new_entry = (f[i] + u[i-1] + u[i+1] 
                + u[i-n] + u[i+n]) / 4;
            infnorm = std::max(infnorm,
                std::abs(new_entry - u[i]));
            u[i] = new_entry;
        }
        //last line block 1
        new_entry = (f[(j+1)*n-1] + u[(j+1)*n-2]  
            + u[(j)*n-1] + (j+2)*n-1) / 4;
        infnorm = std::max(infnorm,
            std::abs(new_entry - u[(j+1)*n-1]));
        u[(j+1)*n-1] = new_entry;
    }
    return infnorm;
}

//Left -1*unity-matrix block not present
void FDM::apply_top_block(const size_t n,
        const estd::vector_t<double>& x,
        estd::vector_t<double>& y) {

    //Special case first row
    y[0] = 4*x[0] - x[1] - x[n];
    for (size_t i = 1; i < n-1; ++i)
        y[i] = -x[i-1] + 4*x[i] - x[i+1] - x[i+n];
    //n-th row
    y[n-1] = -x[n-2] + 4*x[n-1] - x[2*n-1];
}

//All blocks present, no special cases
void FDM::apply_middle_blocks(const size_t n,
        const estd::vector_t<double>& x, 
        estd::vector_t<double>& y) {
    
    for(size_t i = 1; i < n - 1; ++i){ 
        y[i*n] = -x[i*n-n] + 4*x[i*n] - x[i*n+1] - x[i*n+n];
        for(size_t j = i*n + 1; j < (i+1)*n - 1; ++j)
            y[j] = -x[j-n] - x[j-1] + 4*x[j] -x[j+1] - x[j+n];
        y[(i+1)*n-1] = -x[i*n-1] - x[(i+1)*n-2] + 4*x[(i+1)*n-1] - x[(i+2)*n-1];
    }
}

//Right -1*unity-matrix block not present
void FDM::apply_bottom_block(const size_t n,
        const estd::vector_t<double>& x, 
        estd::vector_t<double>& y) {
    y[n*n-n] = -x[n*n-2*n] + 4*x[n*n-n] - x[n*n-n+1];
    for(size_t i = n*n - n; i < n*n - 1; ++i)
        y[i] = -x[i-n] - x[i-1] + 4*x[i] - x[i+1];
    //Special case last row
    y[n*n - 1] = -x[n*n-n-1] - x[n*n-2] + 4*x[n*n-1];
}


//calculates y=L_h*x
void FDM::apply_laplace(const size_t n, //number of nodes
        const estd::vector_t<double>& x, //input
        estd::vector_t<double>& y) {     //output

    //Dimension check
    assert( n*n==x.size() && x.size()==y.size());
    //Check for evil self assignment
    assert( &x != &y);

    //Handle special cases
    apply_top_block(n,x,y);
    apply_middle_blocks(n,x,y);
    apply_bottom_block(n,x,y);
}

//Initialise the RHS-vector for the discretisation
void FDM::init_rhs(const size_t n,   //Number of nodes per row
                estd::vector_t<double>& f,    //Target vector
                     //the boundary function w/o values for corners
                const std::function<double(double,double)>& g
                ) {
    //dimension check
    assert( n*n == f.size() );

    //Zero out f
    std::fill(f.begin(),f.end(),0); 

    //The first n nodes are the top-border-nodes:
    for(size_t i = 1; i < n-1; ++i)
        f[i] = g(static_cast<double>(i)/(n-1),1);

    //The left-border nodes: 
    for(size_t i = 1; i < n-1; ++i)
        f[i*n] = g(0,1-static_cast<double>(i)/(n-1));

    //The right-border nodes:
    for(size_t i = 2; i < n; ++i)
        f[i*n-1] = g(1,1-static_cast<double>(i-1)/(n-1));

    //The bottom-border nodes:
    for(size_t i = 1; i < n-1; ++i)
        f[n*n-n+i] = g(static_cast<double>(i)/(n-1),0);

    //Handle the edges
    f[0] = .5*(f[1] + f[n]);
    f[n-1] = .5*(f[n-2] + f[2*n-1]);
    f[n*n-n] = .5*(f[n*n-2*n] + f[n*n-n+1]);
    f[n*n-1] = .5*(f[n*n-n-1] + f[n*n-2]);
}

//calculates the j-th entry of the gauss-seidel vector
double FDM::gauss_seidel_entry (
        const size_t n,             //The system is n*nxn*n
        const estd::vector_t<double>& f,   //RHS of the system
        const estd::vector_t<double>& u,   //entries 0 to j-1 contain entries
                                    //of the updated vector, j to n-1
                                    //entries of the old one
        const size_t j) {           //The entry to calculate

        double sum{0};

        if (j < n) {                    //first block
            if (j == 0)                 //first line
                sum = -u[1] - u[n];
            else if (j < n-1)           //middle line of first block
                sum = -u[j-1] - u[j+1] - u[j+n];
            else                        //last line block 1
              sum = -u[j-1] - u[j+n];
        }
        else if (j < n*n-n) {           //middle block
            if (j%n == 0)               //first line of block
                sum = -u[j-n] - u[j+1] - u[j+n];
            else if ((j+1)%n == 0)        //last line of block
                sum = -u[j-n] - u[j-1] - u[j+n];
            else                        //middle line
                sum = -u[j-n] - u[j-1] - u[j+1] - u[j+n];
        }
        else  {                         //last block
            if (j == n*n -n)            //first line of block
                sum = -u[j-n] - u[j+1];
            else if(j == n*n - 1)       //last line
                sum = -u[j-n] - u[j-1];
            else                        //middle line of last block
                sum = -u[j-n] - u[j-1] - u[j+1];
        }

        return (f[j]-sum)/4;
}

//calculates the j-th entry of the gauss-seidel vector
double FDM::SOR_entry (
        const size_t n,             //The system is n*nxn*n
        const estd::vector_t<double>& f,   //RHS of the system
        const estd::vector_t<double>& u,   //entries 0 to j-1 contain entries
                                    //of the updated vector, j to n-1
                                    //entries of the old one
        const size_t j,              //The entry to calculate
        const double omega) {       //the relaxion parameter

        double sum{0};

        if (j < n) {                    //first block
            if (j == 0)                 //first line
                sum = -u[1] - u[n];
            else if (j < n-1)           //middle line of first block
                sum = -u[j-1] - u[j+1] - u[j+n];
            else                        //last line block 1
              sum = -u[j-1] - u[j+n];
        }
        else if (j < n*n-n) {           //middle block
            if (j%n == 0)               //first line of block
                sum = -u[j-n] - u[j+1] - u[j+n];
            else if ((j+1)%n == 0)        //last line of block
                sum = -u[j-n] - u[j-1] - u[j+n];
            else                        //middle line
                sum = -u[j-n] - u[j-1] - u[j+1] - u[j+n];
        }
        else  {                         //last block
            if (j == n*n -n)            //first line of block
                sum = -u[j-n] - u[j+1];
            else if(j == n*n - 1)       //last line
                sum = -u[j-n] - u[j-1];
            else                        //middle line of last block
                sum = -u[j-n] - u[j-1] - u[j+1];
        }

        return (f[j]-sum)/4*omega + (1-omega)*u[j];
}


//Solves L_h*u=x, returns (achieved acuracy, needed iterations)
std::pair<double,size_t> FDM::gauss_seidel (
        const size_t n,             //The system must be n*nxn*n
        const estd::vector_t<double>& f,  //RHS of the system
        estd::vector_t<double>& u,        //The solution vector
        const double epsmin,        //The best acuracy to achieve
        const size_t itmax) {       //The maximum number of iterations

    //dimension check:
    assert(f.size() == u.size() && u.size() == n*n);


    //u=f;
    //variable to store the current acuracy
    double infnorm;
    //The gauss-seidel-method
    for(size_t i = 1; i <= itmax; ++i) {
        infnorm = 0;        
        //calculate new iteration-vector
       
        infnorm = gauss_seidel_top_block(n,f,u);
        infnorm = std::max(infnorm,
            gauss_seidel_middle_blocks(n,f,u));
        infnorm = std::max(infnorm,
            gauss_seidel_top_block(n,f,u));
        //Check convergence
        //std::cout << i << " " << infnorm << std::endl;
        if (infnorm < epsmin)
            return {infnorm,i};
    }

    //itmax iterations done
    return {infnorm,itmax};
}


//Solves L_h*u=x, returns (achieved acuracy, needed iterations)
std::pair<double,size_t> FDM::SOR (
        const size_t n,
        const estd::vector_t<double>& f,
        estd::vector_t<double>& u,
        const double epsmin,
        const size_t itmax,
        const double omega) { 

    //dimension check:
    assert(f.size() == u.size() && u.size() == n*n);

    //convergence check
    assert((omega > 0) && (omega < 2));


    //u=f;
    //variable to store the current acuracy
    double infnorm;
    //The gauss-seidel-method
    for(size_t i = 1; i <= itmax; ++i) {
        infnorm = 0;        
        //calculate new iteration-vector
        for(size_t j = 0; j < n*n; ++j) {
            //if (!(j%n == 0 || (j+1)%n == 0)){
                double curr_entry = SOR_entry(n,f,u,j,omega);
                //
                //double curr_entry = (f[j] + u[j-n] + u[j+n] + u[j-1] + u[j+1])/4;
    
                //update the infinity norm
                infnorm = std::max(infnorm,
                    std::abs(curr_entry-u[j]));

                //update iteration vector entry
                u[j] = curr_entry;
            //}
        }

        //Check convergence
        //std::cout << i << " " << infnorm << std::endl;
        if (infnorm < epsmin)
            return {infnorm,i};
    }

    //itmax iterations done
    return {infnorm,itmax};
}


std::pair<double,std::size_t> FDM::CG (
        const Matrix& L,
        const estd::vector_t<double>& f,
        estd::vector_t<double>& u,
        const double epsmin,
        const std::size_t itmax) {

    //range checks
    assert(L.size() == f.size() && f.size() == u.size());

    //the initial residual
    estd::vector_t<double> residual = f;
    estd::vector_t<double> temp(L.size());
    L.multMatVec(u,temp);
    estd::axpy(-1,temp,residual);

    //the first search direction
    estd::vector_t<double> direction = residual;

    //The metod itself
    for (size_t i = 1; i <= itmax; ++i) {
        L.multMatVec(direction,temp);
        double alpha = estd::dot_product(residual,residual) 
            / estd::dot_product(direction,temp);

        //update solution vector
        estd::axpy(alpha,direction,u);

        //update residual, save dotproduct of old one for beta-update
        double beta_updater = estd::dot_product(residual,residual);
        estd::axpy(-1*alpha,temp,residual);

        //check convergence:
        double norm_res = estd::dot_product(residual,residual);
        //sstd::cout << norm_res << '\n';
        if (norm_res < epsmin*epsmin)
            return {estd::nrm2(residual),i};

        //update direction
        double beta = norm_res/beta_updater;
        estd::scale(beta,direction);
        estd::axpy(1,residual,direction);
    }

    return {estd::nrm2(residual),itmax};
}

std::pair<double,std::size_t> FDM::PCG (
        const Matrix& L,
        const estd::vector_t<double>& f,
        estd::vector_t<double>& u,
        const double epsmin,
        const std::size_t itmax) {

    //range checks
    assert(L.size() == f.size() && f.size() == u.size());

    //the initial residual
    estd::vector_t<double> residual = f;
    estd::vector_t<double> temp(L.size());
    L.multMatVec(u,temp);
    estd::axpy(-1,temp,residual);

    //the first search direction
    estd::vector_t<double> direction = residual;
    L.applyPrecond(direction);

    //PCG-specific helper vector
    estd:: vector_t<double> helper = direction;

    //The metod itself
    for (size_t i = 1; i <= itmax; ++i) {
        L.multMatVec(direction,temp);
        double alpha = estd::dot_product(residual,helper) 
            / estd::dot_product(direction,temp);

        //update solution vector
        estd::axpy(alpha,direction,u);

        //update residual, save dotproduct of old one for beta-update
        double beta_updater = estd::dot_product(residual,helper);
        estd::axpy(-1*alpha,temp,residual);

        //std::cout << norm_res << '\n';
        if (estd::nrm2(residual) < epsmin)
            return {estd::nrm2(residual),i};

        //update helper
        helper = residual;
        L.applyPrecond(helper);

        //update direction
        double beta = estd::dot_product(residual,helper)/beta_updater;
        estd::scale(beta,direction);
        estd::axpy(1,helper,direction);
    }

    return {estd::nrm2(residual),itmax};
}

