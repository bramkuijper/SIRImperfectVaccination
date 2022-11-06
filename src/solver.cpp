#include <Rcpp.h>
#include <math.h>
#include <iostream>
#include "solver.h"

double Solver::alphaprime(double const v)
{
    return((1.0 - params.r2) * (1.0 - params.r4) * v);
}

double Solver::dalphaprime_dv(double  const v)
{
    return((1.0 - params.r2) * (1.0 - params.r4));
}

double Solver::beta(double const v)
{
    return(params.b1 + pow(v, params.b2));
}

double Solver::dbeta_dv(double const v)
{
    double val = params.b2 * pow(v, (params.b2 - 1.0));

    return(val);
}

double Solver::betaprime(double const v)
{
    return((1.0 - params.r3) * beta((1.0 - params.r2)*v));
}

double Solver::dbetaprime_dv(double const v)
{
    double val = (1.0 - params.r3) * (1.0 - params.r2) * dbeta_dv((1.0 - params.r2)*v);

    assert(!std::isnan(val));

    return(val);
}

double Solver::chi(double const v)
{
    return(params.c1 + pow(v, params.c2));
}

double Solver::dchi_dv(double const v)
{
    return(params.c2 * pow(v, params.c2 - 1.0));
}

double Solver::chiprime(double const v)
{
    return(chi((1.0 - params.r2)*v));
}

double Solver::dchiprime_dv(double const v)
{
    return((1.0 - params.r2) * dchi_dv((1.0 - params.r2)*v));
}

// h as described on p753 in Gandon et al 2001: h = beta * y + beta' * y'
double Solver::h(double const v)
{
    return(beta(v) * y + betaprime(v) * yprime);
}

double Solver::hprime(double const v)
{
    return((1.0 - params.r1) * h(v));
}

double Solver::dR0dalpha(double const v)
{
    double dR0dv_part1 = (dbeta_dv(v) * (x + params.sigma * y) * (params.delta + v + chi(v) + params.sigma * h(v)) - beta(v) * (x + params.sigma * y) * (1 + dchi_dv(v))) /
        pow(params.delta + v + chi(v) + params.sigma * h(v),2.0);

    double dR0dv_part2 = (
            dbetaprime_dv(v) * (1.0 - params.r1) * (xprime + params.sigma * yprime) * (
                params.delta + alphaprime(v) + chiprime(v) + params.sigma * hprime(v)
                )
            - betaprime(v) * (1.0 - params.r1) * (xprime + params.sigma * yprime) * (
                dalphaprime_dv(v) + dchiprime_dv(v)
                ) ) /
        pow(params.delta + alphaprime(v) + chiprime(v) + params.sigma * hprime(v),2.0);

    assert(isnormal(dR0dv_part1) != 0);
    assert(isnormal(dR0dv_part2) != 0);

    return(dR0dv_part1 + dR0dv_part2);
} // end dR0dalpha

Solver::Solver(Params const &parameters) :
    params{parameters}
    ,x{params.n_susceptible_init * (1.0 - params.f)}
    ,xprime{params.n_susceptible_init * params.f}
    ,y{params.n_infected_init * (1.0 - params.f)}
    ,yprime{params.n_infected_init * params.f}
    ,alpha{params.alpha_init}
{} // end Solver

void Solver::solve_ecol_dynamic()
{
    double dxdt, dydt, dxprimedt, dyprimedt;
    double xtplus1, xprimetplus1, ytplus1, yprimetplus1;
    bool converged_ecol;

    for (int ecol_time_step = 0; ecol_time_step < params.maxt_eco;
            ++ecol_time_step)
    {
        // work out eq. 6 from Gandon et al 2001
        dxdt = (1.0 - params.f) * params.lambda -
                    (params.delta + h(alpha)) * x + chi(alpha) * y;


        dxprimedt = params.f * params.lambda -
                    (params.delta + hprime(alpha)) * xprime
                    + chiprime(alpha) * yprime;

        dydt = h(alpha) * x - (params.delta + alpha + chi(alpha)) * y;

        dyprimedt = hprime(alpha) * xprime -
            (params.delta + alphaprime(alpha) + chiprime(alpha) ) * yprime;

        xtplus1 = x + params.eul_eco * dxdt;

        if (xtplus1 < 0)
        {
            xtplus1 = 0;
        }

        ytplus1 = y + params.eul_eco * dydt;

        if (ytplus1 < 0)
        {
            ytplus1 = 0;
        }

        xprimetplus1 = xprime + params.eul_eco * dxprimedt;

        if (xprimetplus1 < 0)
        {
            xprimetplus1 = 0;
        }

        yprimetplus1 = yprime + params.eul_eco * dyprimedt;

        if (yprimetplus1 < 0)
        {
            yprimetplus1 = 0;
        }

//        std::cout << "ecol_time_step: " << ecol_time_step << " dxdt: " << dxdt << " dxprimedt: " << dxprimedt << " dydt: " << dydt << " yprime: " << dyprimedt << std::endl;

        converged_ecol = true;

        if (std::fabs(xtplus1 - x) > params.convergence_limit)
        {
            converged_ecol = false;
        }

        if (std::fabs(ytplus1 - y) > params.convergence_limit)
        {
            converged_ecol = false;
        }

        if (std::fabs(xprimetplus1 - xprime) > params.convergence_limit)
        {
            converged_ecol = false;
        }

        if (std::fabs(yprimetplus1 - yprime) > params.convergence_limit)
        {
            converged_ecol = false;
        }

        x = xtplus1;
        xprime = xprimetplus1;
        y = ytplus1;
        yprime = yprimetplus1;

        if (converged_ecol)
        {
            break;
        }

        if (ecol_time_step % 1000 == 0)
        {
            Rcpp::checkUserInterrupt();
        }

//        std::cout << "ecol_time_step: " << ecol_time_step << " x: " << x << " xprime: " << xprime << " y: " << y << " yprime: " << yprime << std::endl;
    }// end for int ecol time
}

Rcpp::DataFrame Solver::run_sep_timescale()
{
    bool converged_evol;

    double alphatplus1;

    int evol_time_step;

    int vector_size = params.all_data ? params.maxt_evo : 1;

    // allocate vectors to fill the dataframe
    Rcpp::NumericVector time_steps_vec(vector_size);
    Rcpp::NumericVector x_vec(vector_size);
    Rcpp::NumericVector y_vec(vector_size);
    Rcpp::NumericVector xprime_vec(vector_size);
    Rcpp::NumericVector yprime_vec(vector_size);
    Rcpp::NumericVector alpha_vec(vector_size);

    for (evol_time_step = 0; evol_time_step < params.maxt_evo; ++evol_time_step)
    {
        solve_ecol_dynamic();

        alphatplus1 = alpha + params.eul_evo * dR0dalpha(alpha);

//        std::cout << "a(t+1) " << alphatplus1 << " " << alpha << " " << params.eul_evo * dR0dalpha(alpha) << " x: " << x << " y: " << y << " xprime: " << xprime << " yprime: " << yprime <<  std::endl;

        if (alphatplus1 < 0)
        {
            alphatplus1 = 0.0;
        }

        converged_evol = true;

        if (std::fabs(alphatplus1 - alpha) > params.convergence_limit)
        {
            converged_evol = false;
        }

        alpha = alphatplus1;

        if (converged_evol)
        {
            break;
        }

        if (evol_time_step % 1000 == 0)
        {
            Rcpp::checkUserInterrupt();
        }
        if (params.all_data)
        {
            time_steps_vec[evol_time_step] = evol_time_step;
            x_vec[evol_time_step] = x;
            xprime_vec[evol_time_step] = xprime;
            y_vec[evol_time_step] = y;
            yprime_vec[evol_time_step] = yprime;
            alpha_vec[evol_time_step] = alpha;
        }
    } // end for int evol time

    if (!params.all_data)
    {
        time_steps_vec[0] = evol_time_step;
        x_vec[0] = x;
        xprime_vec[0] = xprime;
        y_vec[0] = y;
        yprime_vec[0] = yprime;
        alpha_vec[0] = alpha;
    }

    Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named("t_converge") = time_steps_vec
            ,Rcpp::Named("x") = x_vec
            ,Rcpp::Named("xprime") = xprime_vec
            ,Rcpp::Named("y") = y_vec
            ,Rcpp::Named("yprime") = yprime_vec
            ,Rcpp::Named("alpha") = alpha_vec
            );

    return(df);

} // end solve_SIR_alpha
