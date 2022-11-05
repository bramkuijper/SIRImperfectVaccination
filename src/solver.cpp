#include <Rcpp.h>
#include "solver.h"

double Solver::beta(double const alpha)
{
    return(params.b1 + pow(alpha, params.b2));
}

double Solver::chi(double const alpha)
{
    return(params.c1 + pow(alpha, params.c1));
}

// h as described on p753 in Gandon et al 2001: h = beta * y + beta' * y'
double Solver::h(double const alpha)
{
    return(beta(alpha) * y + betaprime(alpha) * yprime);
}

Solver::Solver(Params const &parameters) :
    params{parameters}
{

} // end Solver


Rcpp::DataFrame Solver::run()
{
    double x = params.density_init * (1.0 - params.f);
    double xprime = params.density_init * params.f;
    double y = params.density_init * (1.0 - params.f);
    double yprime = params.density_init * params.f;

    double delta = params.delta;
    double lambda = params.lambda;
    
    double dxdt = 0.0;
    double dydt = 0.0;
    double dxprimedt = 0.0;
    double dyprimedt = 0.0;

    Rcpp::DataFrame df = Rcpp::DataFrame::create();

    for (int time_step = 0; time_step < params.maxt; ++time_step)
    {
        dxdt = (1.0 - params.f) * params.lambda - 
                    (params.delta + h(alpha)) * x + chi(alpha) * y;
    }


    return(df);
} // end solve_SIR_alpha
