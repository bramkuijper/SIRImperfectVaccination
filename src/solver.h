#ifndef _SOLVER_H_
#define _SOLVER_H_

#include <Rcpp.h>

struct Params
{
    double n_susceptible_init=100;
    double n_infected_init=100;
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
    double eul_eco = 0.0001;
    double eul_evo = 0.0001;
    double sigma = 0.0;
    int maxt_evo = 5;
    int maxt_eco = 5;
    double convergence_limit = 1e-07;
    double alpha_init = 1.0;
    bool all_data = false;
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
        double hprime(double const alpha);
        double alphaprime(double const alpha);
        double beta(double const alpha);
        double betaprime(double const alpha);

        double chi(double const alpha);
        double chiprime(double const alpha);

        double dalphaprime_dv(double const alpha);
        double dbeta_dv(double const alpha);
        double dbetaprime_dv(double const alpha);
        double dchi_dv(double const alpha);
        double dchiprime_dv(double const alpha);

        double dR0dalpha(double const alpha);

    public:
        Solver(Params const &parameters);
        void solve_ecol_dynamic();
        Rcpp::DataFrame run_sep_timescale();
};

#endif
