#include <cassert>
#include <iostream>
#include <utility>
#include "laplace_matrix_routines.hpp"
#include "laplace_matrix.hpp"
#include "laplace_matrix_jacobi.hpp"
#include "laplace_matrix_SSOR.hpp"
#include "laplace_matrix_IC.hpp"
#include "my_vector.hpp"
#include "vector_routines.hpp"

//Die Randwertfunktion von Blatt 1
double g_blatt1 (double x, double y) {

    //Randwert?
    assert( x==0 || y==0 || x==1 || y==1);

    //Ecke?
    assert( !(x==0 && y==0) &&
            !(x==1 && y==0) &&
            !(x==0 && y==1) &&
            !(x==1 && y==1) );

    if (y == 0)
        return 100;
    if (y == 1)
        return 20;
    return 40;
}



//Main Blatt 1
int main () {
    size_t n;
    double eps;
    std::cout << "Number of Nodes per row:\n";
    std::cin >> n;
    std::cout << "Required precision for Gauss-Seidel:\n";
    std::cin >> eps;
    estd::vector_t<double> rhs(n*n);
    estd::vector_t<double> solution(n*n);


    std::cout << "Setting RHS... " << std::flush;
    FDM::init_rhs(n,rhs,g_blatt1);
    //init_rhs(n,rhs,[](double x,double y){return x + y;});
    std::cout << "Done.\nSolving with Gauss-Seidel... "<< std::flush;
    auto x = FDM::gauss_seidel(n,rhs,solution,eps,n*n*n);
    std::cout << "done.\n";

    
    std::cout << "Precision reached: " << x.first
        << "\nIterations needed: " << x.second << "\n"; 
    createVTKFile(n,solution.data());

}


/*
//Main Blatt 2
int main () {
    size_t n;
    double eps, omega;
    std::cout << "Number of Nodes per row:\n";
    std::cin >> n;
    std::cout << "Required precision for SOR:\n";
    std::cin >> eps;
    std::cout << "Relaxion-parameter:\n";
    std::cin >> omega;
    estd::vector_t<double> rhs(n*n);
    estd::vector_t<double> solution(n*n);

    std::cout << "Setting RHS... " << std::flush;
    FDM::init_rhs(n,rhs,g_blatt1);
    //init_rhs(n,rhs,[](double x,double y){return x + y;});
    std::cout << "Done.\nSolving with SOR... "<< std::flush;
    auto x = FDM::SOR(n,rhs,solution,eps,n*n*n,omega);
    std::cout << "done.\n";

    
    std::cout << "Precision reached: " << x.first
        << "\nIterations needed: " << x.second << "\n"; 
    createVTKFile(n,solution.data());

}
*/

/*
//Main Blatt 3 CG
int main () {
    size_t n;
    double eps;
    std::cout << "Number of Nodes per row:\n";
    std::cin >> n;
    std::cout << "Required precision for CG:\n";
    std::cin >> eps;
    estd::vector_t<double> rhs(n*n);
    estd::vector_t<double> solution(n*n);

    std::cout << "Setting RHS... " << std::flush;
    FDM::init_rhs(n,rhs,g_blatt1);
    //init_rhs(n,rhs,[](double x,double y){return x + y;});
    std::cout << "Done.\nSolving with CG... "<< std::flush;
    FDM::LaplaceMatrix L{n};
    auto x = FDM::CG(L,rhs,solution,eps,n*n*n);
    std::cout << "done.\n";

    
    std::cout << "Precision reached: " << x.first
        << "\nIterations needed: " << x.second << "\n"; 
    createVTKFile(n,solution.data());

}
*/


/*
//Main Blatt 3 PCG
int main () {
    size_t n;
    double eps;
    std::cout << "Number of Nodes per row:\n";
    std::cin >> n;
    std::cout << "Required precision for PCG:\n";
    std::cin >> eps;
    estd::vector_t<double> rhs(n*n);
    estd::vector_t<double> solution(n*n);

    std::cout << "Setting RHS... " << std::flush;
    FDM::init_rhs(n,rhs,g_blatt1);
    //init_rhs(n,rhs,[](double x,double y){return x + y;});
    std::cout << "Done.\nSolving with PCG... "<< std::flush;
    FDM::LaplaceMatrixJacobi L{n};
    auto x = FDM::PCG(L,rhs,solution,eps,n*n*n);
    std::cout << "done.\n";

    
    std::cout << "Precision reached: " << x.first
        << "\nIterations needed: " << x.second << "\n"; 
    createVTKFile(n,solution.data());

}
*/

/*
//Main Blatt 3 SSOR
int main () {
    size_t n;
    double eps;
    std::cout << "Number of Nodes per row:\n";
    std::cin >> n;
    std::cout << "Required precision for SSOR:\n";
    std::cin >> eps;
    estd::vector_t<double> rhs(n*n);
    estd::vector_t<double> solution(n*n);

    std::cout << "Setting RHS... " << std::flush;
    FDM::init_rhs(n,rhs,g_blatt1);
    //init_rhs(n,rhs,[](double x,double y){return x + y;});
    std::cout << "Done.\nSolving with SSOR... "<< std::flush;
    FDM::LaplaceMatrixSSOR L{n};
    auto x = FDM::PCG(L,rhs,solution,eps,n*n*n);
    std::cout << "done.\n";

    
    std::cout << "Precision reached: " << x.first
        << "\nIterations needed: " << x.second << "\n"; 
    createVTKFile(n,solution.data());

}
*/

/*
//Main Blatt 3 IC
int main () {
    size_t n;
    double eps;
    std::cout << "Number of Nodes per row:\n";
    std::cin >> n;
    std::cout << "Required precision for IC-PCG:\n";
    std::cin >> eps;
    estd::vector_t<double> rhs(n*n);
    estd::vector_t<double> solution(n*n);

    std::cout << "Setting RHS... " << std::flush;
    FDM::init_rhs(n,rhs,g_blatt1);
    //init_rhs(n,rhs,[](double x,double y){return x + y;});
    std::cout << "Done.\nSolving with IC-PCG... "<< std::flush;
    FDM::LaplaceMatrixIC L{n};
    auto x = FDM::PCG(L,rhs,solution,eps,n*n*n);
    std::cout << "done.\n";

    
    std::cout << "Precision reached: " << x.first
        << "\nIterations needed: " << x.second << "\n"; 
    createVTKFile(n,solution.data());

}
*/

/*
//small test
int main () {
    estd::vector_t<double> test1(5,2), test2(5,4);
    estd::axpy(-1,test2,test1);

    for (auto v:test1)
        std::cout << v;
    std::cout << '\n' << estd::dot_product(test1,test1) << '\n';
}
*/