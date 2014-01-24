#include <cassert>
#include <iostream>
#include "my_vector.hpp"
#include <utility>
#include "laplace_matrix_routines.hpp"

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



//first testcase
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

    /*
    for (auto x:rhs)
        std::cout << x << "\n";

    for (size_t i = 0; i < n; ++i){
        for (size_t j = 0; j < n; ++j) 
            std::cout << solution[i*n+j] << " ";
        std::cout << "\n";
    }
    */
    std::cout << "Precision reached: " << x.first
        << "\nIterations needed: " << x.second << "\n"; 
    createVTKFile(n,solution.data());

}
