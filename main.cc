#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <utility>
#include "euler.hpp"
#include "laplace_matrix_routines.hpp"
#include "laplace_matrix.hpp"
#include "laplace_matrix_jacobi.hpp"
#include "laplace_matrix_SSOR.hpp"
#include "laplace_matrix_IC.hpp"
#include "multigrid.hpp"
#include "my_vector.hpp"
#include "vector_routines.hpp"

//for the timetests
extern "C" double cputime (const double);


//For best omega
constexpr double PI = 3.141592653589793238462643383279502884197169399375105820974944;



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


/*
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
*/


/*
//Main for time-tests:
int main () {

    std::ofstream out ("FMGM_Laufzeiten.txt");

    for (size_t i = 64; i < 512; i*=2) {

        size_t n =  i+1;

        estd::vector_t<double> rhs(n*n);

        FDM::init_rhs(n,rhs,g_blatt1);

        //for (size_t j = 0; j < 2; ++j, eps/=100) {
    
        
        estd::vector_t<double> solution(n*n);

        

        auto start = cputime(0.0);
        FDM::FMGM(n,solution,rhs,4);
        auto x = FDM::MGM(n,solution,rhs,1e-6,4,1);
        auto ende = cputime(start);
        std::fill(solution.begin(),solution.end(),0);

        out << n << '\t' << x << '\t' << ende 
            << "\n";

        start = cputime(0.0);
        FDM::FMGM(n,solution,rhs,4);
        x = FDM::MGM(n,solution,rhs,1e-6,4,2);
        ende = cputime(start);
        //std::fill(solution.begin(),solution.end(),0);

        out  << x << '\t' << ende 
            << "\n";
        

        
        FDM::LaplaceMatrix L{i};
        auto start = cputime(0.0);
        auto x = FDM::CG(L,rhs,solution,eps,i*i*i*i);
        auto ende = cputime(start);
        std::fill(solution.begin(),solution.end(),0);

        out << i <<  '\t' << eps <<'\t' << x.second << '\t' << ende 
            << "\t";

        FDM::LaplaceMatrixSSOR L2{i};
        start = cputime(0.0);
        x = FDM::PCG(L2,rhs,solution,eps,i*i*i*i);
        ende = cputime(start);
        std::fill(solution.begin(),solution.end(),0);

        out << i << '\t' << eps <<'\t' << x.second << '\t' << ende 
            << "\t";

        
        start = cputime(0.0);
        x = FDM::PCG(L3,rhs,solution,eps,i*i*i*i);
        ende = cputime(start);
        std::fill(solution.begin(),solution.end(),0);

        out << i <<  '\t' << eps <<'\t' << x.second << '\t' << ende 
            << "\n";

        }

        std::cerr  << i;
        
    }
}
*/

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
//Main Blatt 4
int main () {
    size_t n;
    double eps;
    std::cout << "Number of Nodes per row:\n";
    std::cin >> n;
    std::cout << "Required precision for TGM:\n";
    std::cin >> eps;
    estd::vector_t<double> rhs(n*n);
    estd::vector_t<double> solution(n*n);

    std::cout << "Setting RHS... " << std::flush;
    FDM::init_rhs(n,rhs,g_blatt1);
    //init_rhs(n,rhs,[](double x,double y){return x + y;});
    std::cout << "Done.\nSolving with TGM... "<< std::flush;
    auto x = FDM::TGM(n,solution,rhs,eps);
    std::cout << "done.\n";

    
    std::cout << "\nIterations needed: " << x << "\n"; 
    createVTKFile(n,solution.data());

}
*/

/*
//Main Blatt  5 MGM
int main () {
    size_t n,depth,mu;
    double eps;
    std::cout << "Number of Nodes per row:\n";
    std::cin >> n;
    std::cout << "Required precision for MGM:\n";
    std::cin >> eps;
    std::cout << "Required depth for multigrid:\n";
    std::cin >> depth;
    std::cout << "Required multigrid steps per level:\n";
    std::cin >> mu;
    estd::vector_t<double> rhs(n*n);
    estd::vector_t<double> solution(n*n);

    std::cout << "Setting RHS... " << std::flush;
    FDM::init_rhs(n,rhs,g_blatt1);
    //init_rhs(n,rhs,[](double x,double y){return x + y;});
    std::cout << "Done.\nSolving with MGM... "<< std::flush;
    auto x = FDM::MGM(n,solution,rhs,eps,depth,mu);
    std::cout << "done.\n";

    
    std::cout << "\nIterations needed: " << x << "\n"; 
    createVTKFile(n,solution.data());

}
*/

/*
//Main Blatt  5 FMGM
int main () {
    size_t n,depth,mu;
    double eps;
    std::cout << "Number of Nodes per row:\n";
    std::cin >> n;
    std::cout << "Required precision for MGM:\n";
    std::cin >> eps;
    std::cout << "Required depth for multigrid:\n";
    std::cin >> depth;
    std::cout << "Required multigrid steps per level:\n";
    std::cin >> mu;
    estd::vector_t<double> rhs(n*n);
    estd::vector_t<double> solution(n*n);

    std::cout << "Setting RHS... " << std::flush;
    FDM::init_rhs(n,rhs,g_blatt1);
    FDM::FMGM(n,solution,rhs,depth);
    //init_rhs(n,rhs,[](double x,double y){return x + y;});
    std::cout << "Done.\nSolving with MGM... "<< std::flush;
    auto x = FDM::MGM(n,solution,rhs,eps,depth,mu);
    std::cout << "done.\n";

    
    std::cout << "\nIterations needed: " << x << "\n"; 
    createVTKFile(n,solution.data());

}
*/


//expl euler
int main () {
    size_t n, num_steps;
    double stepsize;

    std::cout << "Number of Nodes per row:\n";
    std::cin >> n;
    std::cout << "Number of timesteps:\n";
    std::cin >> num_steps;
    std::cout << "Size of step:\n";
    std::cin >> stepsize;

    estd::vector_t<double> rhs(n*n);
    estd::vector_t<double> solution(n*n);
    estd::vector_t<double> solution2(n*n);

    std::cout << "Setting RHS... " << std::flush;
    FDM::init_rhs(n,rhs,g_blatt1);
    
    std::cout << "Done.\nSolving with explicit euler... "<< std::flush;
    FDM::EulerExpl(n,solution2,rhs,stepsize,num_steps);
    std::cout << "done.\n";
}



/*
//impl euler
int main () {
    size_t n, num_steps;
    double stepsize;

    std::cout << "Number of Nodes per row:\n";
    std::cin >> n;
    std::cout << "Number of timesteps:\n";
    std::cin >> num_steps;
    std::cout << "Size of step:\n";
    std::cin >> stepsize;

    estd::vector_t<double> rhs(n*n);
    estd::vector_t<double> solution(n*n);
    estd::vector_t<double> solution2(n*n);

    std::cout << "Setting RHS... " << std::flush;
    FDM::init_rhs(n,rhs,g_blatt1);
    std::cout << "Done.\nSolving with implicit euler... "<< std::flush;
    FDM::EulerImpl(n,solution2,rhs,stepsize,num_steps);
    std::cout << "done.\n";
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