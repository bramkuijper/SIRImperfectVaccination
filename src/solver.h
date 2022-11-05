#ifndef _SOLVER_H_
#define _SOLVER_H_

#include <Rcpp.h>

struct Params
{
    double density_init=100;
    double delta = 1;
    double lambda = 25;
    double f = 0.2;
    double b1 = 1.0;
    double b2 = 0.0;
    double c1 = 1.0;
    double c2 = 0.0;
    double r1 = 0.0;
    double r2 = 0.0;
    double r3 = 0.0;
    double r4 = 0.0;
    double eul = 0.01;
    int maxt = 1e08;
};

class Solver
{
    private:
        Params params;
        double y;
        double yprime;
        double x;
        double xprime;
        double alpha;

        double h(double const alpha);
        double beta(double const alpha);
        double betaprime(double const alpha);
        double chi(double const alpha);

    public:
        Solver(Params const &parameters);
        Rcpp::DataFrame run();
};

#endif
