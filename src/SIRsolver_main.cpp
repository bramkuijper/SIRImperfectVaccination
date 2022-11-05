#include <Rcpp.h>
#include "solver.h"

//' Run a single SIR simulation a la eqns. 6, 7 from Gandon
//' et al 2001 Nature to optimize virulence in the face of
//' different medical interventions
//' @param density_init initial numbers (equal for susceptible and infected) (floating-point)
//' @param f fraction vaccinated (floating-point)
//' @param lambda birth rate (floating-point)
//' @param delta parasite-independent mortality rate (floating-point)
//' @param b1 intercept of trade-off curve virulence-induced mortality vs transmission (floating-point)
//' @param b2 power of trade-off curve virulence-induced mortality vs transmission (floating-point)
//' @param c1 intercept of trade-off curve virulence-induced mortality vs recovery (floating-point)
//' @param c2 power of trade-off curve virulence-induced mortality vs recovery (floating-point)
//' @param r1 anti-infection immunity (floating-point value)
//' @param r2 anti-growth rate immunity (floating-point value)
//' @param r3 transmission-blocking immunity (floating-point value)
//' @param r4 anti-toxin immunity (floating-point value)
//' @return A \code{list} that contains the values of the equilibrium 
//'     densities x, x', y, y' and virulence alpha
//' @examples
//'     # test run of a simulation
//'     list.result <- SIRsolver(
//'                         density_init=100,
//'                         f=0.2,
//'                         lambda=25,
//'                         delta=1,
//'                         b1=1,
//'                         b2=0.5,
//'                         c1=1,
//'                         c2=0.2,
//'                         r1=0,
//'                         r2=0,
//'                         r3=0,
//'                         r4=0)
//'
//' @export 
// [[Rcpp::export]]
Rcpp::DataFrame SIRsolver(
         double density_init=100,
         double f=0.3,
         double lambda=25,
         double delta=1,
         double b1=1,
         double b2=0.5,
         double c1=1,
         double c2=0.2,
         double r1=0.0,
         double r2=0.0,
         double r3=0.0,
         double r4=0.0)
{
    // initialize a struct with parameters
    Params pars;

    pars.density_init = density_init;
    pars.f = f; 
    pars.lambda = lambda;
    pars.delta = delta;
    pars.b1 = b1;
    pars.b2 = b2;
    pars.c1 = c1;
    pars.c2 = c2;
    pars.r1 = r1;
    pars.r2 = r2;
    pars.r3 = r3;
    pars.r4 = r4;

    Solver simulateSIR(pars);

    Rcpp::DataFrame output = simulateSIR.run();

    return(output);
} // end SIRsolver
        
