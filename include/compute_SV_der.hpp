/**
 * \file dirichlet.hpp
 * \brief This file defines lowest order indirect/direct BVP solvers to solve a
 * Dirichlet Boundary Value problem of the form given in \f$\eqref{eq:dirbvp}\f$
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#ifndef COMPUTESVDERHPP
#define COMPUTESVDERHPP

#include "parametrized_mesh.hpp"

namespace parametricbem2d {
    std::complex<double> compute_SV_der(const ParametrizedMesh &mesh,
                                        unsigned order,
                                        const double c,
                                        const double k_o,
                                        const double k_i);
    std::complex<double> compute_SV_der_debug(const ParametrizedMesh &mesh,
                                              unsigned order,
                                              const double c,
                                              const double k_o,
                                              const double k_i);
    std::complex<double> eval_SV_der(const ParametrizedMesh &mesh,
                                     unsigned order,
                                     const double c,
                                     const double k_o,
                                     const double k_i);
    double eval_sv(const ParametrizedMesh &mesh,
                   unsigned order,
                   const double k_o,
                   const double k_i);
    double eval_sv_debug(const ParametrizedMesh &mesh,
                         unsigned order,
                         const double k_o,
                         const double k_i);
    double diffex ( std::function<double(double)> f , const double x , const double h0 ,
                    const double rtol , const double atol );

} // namespace parametricbem2d
#endif // DIRICHLETHPP
