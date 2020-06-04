//
// Created by diegorenner on 2/14/20.
//
#include <complex>

#ifndef ADVNUMCSE_GENERATE_SOLUTION_H
#define ADVNUMCSE_GENERATE_SOLUTION_H
namespace sol {
    std::complex<double> hn(int n,
                            double x);

    std::complex<double> fund_sol(double k, double x1, double x2, double ipt1, double ipt2);
    std::complex<double> fund_sol_N(double k, double x1, double x2, double ipt1, double ipt2);

    std::complex<double> hn_der(int n,
                                double x);

    std::complex<double> jn_der(int n,
                  double x);

    std::complex<double> r_coeff(int n,
                                 double r,
                                 double k,
                                 double n_i);

    std::complex<double> t_coeff(int n,
                                 double eps,
                                 double k,
                                 double n_i);

    std::complex<double> u_i(double x1,
                             double x2,
                             int l,
                             double eps,
                             double *a_n,
                             double k,
                             double n_i);

    std::complex<double> u_s(double x1,
                             double x2,
                             int l,
                             double eps,
                             double *a_n,
                             double k,
                             double n_i);

    std::complex<double> u_t(double x1,
                             double x2,
                             int l,
                             double eps,
                             double *a_n,
                             double k,
                             double n_i);

    std::complex<double> u_i_N(double x1,
                               double x2,
                               int l,
                               double eps,
                               double *a_n,
                               double k,
                               double n_i);

    std::complex<double> u_s_N(double x1,
                               double x2,
                               int l,
                               double eps,
                               double *a_n,
                               double k,
                               double n_i);

    std::complex<double> u_t_N(double x1,
                               double x2,
                               int l,
                               double eps,
                               double *a_n,
                               double k,
                               double n_i);
}
#endif //ADVNUMCSE_GENERATE_SOLUTION_H
