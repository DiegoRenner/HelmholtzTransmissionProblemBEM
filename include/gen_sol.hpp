/**
 * \file gen_sol.hpp
 * \brief This file contains functions that generate solutions
 * for the Helmholtz Equation and the Helmholtz transmission problem.
 *
 * This File is a part of the HelmholtzTransmissionProblemBEM library.
 */


#ifndef GEN_SOLH
#define GEN_SOLH

#include <complex>

/**
 * \namespace sol
 * \brief This namesapce nontains all functions for generating solutions.
 */
namespace sol {

    /**
     * This function computes the derivative of the bessel function of order \p n at \p x.
     * @param n order of bessel function for which to compute derivative
     * @param x point at which to evaluate
     * @return value of derivative of bessel function of order \p n at point \p x
     */
    std::complex<double> jn_der(int n,
                                double x);

    /**
     * This function takes the order \p n and a point \p x on the real axis
     * and returns the value of the first kind Hankel function of
     * <tt>n</tt>-th order at the point \p x.
     * @param n order of Hankel function to be computed
     * @param x point at which to evaluate
     * @return value of first kind Hankel function of order n at point x
     */
    std::complex<double> hn(int n,
                            double x);

    /**
     * This function takes the order \p n and a point \p x on the real axis
     * and returns the value of the derivative of the first kind Hankel function of
     * <tt>n</tt>-th order at the point \p x.
     * @param n order of Hankel function to be computed
     * @param x point at which to evaluate
     * @return value of the derivative of the first kind Hankel function of order n at point x
     */
    std::complex<double> hn_der(int n,
                                double x);

    /**
     * This function takes a wavenumber \p k, two points in \f$ \mathbb{R}^2 \f$ defined by \p x1, \p x2
     * and \p ipt1, \p ipt2 and returns the value of the Dirichlet data of the fundamental
     * solution to the Helmholtz equation with the interior evaluation point defined
     * by \p ipt1, \p ipt2 at the point defined by \p x1, \p x2.
     * @param k wavenumber
     * @param x1 x-value of point at which to evaluate
     * @param x2 y-value of point at which to evaluate
     * @param ipt1 x-value of interior evaluation point
     * @param ipt2 y-value of interior evaluation point
     * @return Dirichlet data of fundamental solution to the Helmholtz equation
     * with interior evaluation point defined by \p ipt1, \p ipt2 at point definde by \p x1, \p x2
     */
    std::complex<double> fund_sol_dir(double k,
            double x1,
            double x2,
            double ipt1,
            double ipt2);

    /**
     * This function takes a wavenumber \p k, two points in \f$ \mathbb{R}^2 \f$ defined by \p x1, \p x2
     * and \p ipt1, \p ipt2 and returns the value of the Neumann data of the fundamental
     * solution to the Helmholtz equation with the interior evaluation point defined
     * by \p ipt1, \p ipt2 at the point defined by \p x1,\p x2.
     * @param k wavenumber
     * @param x1 x-value of point at which to evaluate
     * @param x2 y-value of point at which to evaluate
     * @param ipt1 x-value of interior evaluation point
     * @param ipt2 y-value of interior evaluation point
     * @return Neumann data of fundamental solution to the Helmholtz equation
     * with interior evaluation point defined by \p ipt1, \p ipt2 at point definde by \p x1, \p x2
     */
    std::complex<double> fund_sol_neu(double k,
            double x1,
            double x2,
            double ipt1,
            double ipt2);

    /**
     * This function computes the reflection coefficient of the Helmholtz transmission problem
     * with a refraction index of the inner domain of \p n_i on a 2-dimensional ball of radius \p r.
     * The wavenumber is given by \p k and \p n is the index of the coefficient in the series expansion of the solution.
     * @param n index of coefficient in series expansion of solution
     * @param r radius of ball on which solution is computed
     * @param k wavenuber
     * @param n_i refraction index of inner domain
     * @return reflection coefficient of <tt>n</tt>-th product in series expansion of the solution
     * given radius \p r, wavenumber \p k and refraction index \p n_i
     */
    std::complex<double> r_coeff(int n,
                                 double r,
                                 double k,
                                 double n_i);

    /**
     * This function computes the transmission coefficient of the Helmholtz transmission problem
     * with a refraction index of the inner domain of \p n_i on a 2-dimensional ball of radius \p r.
     * The wavenumber is given by \p k and \p n is the index of the coefficient in the series expansion of the solution.
     * @param n index of coefficient in series expansion of solution
     * @param eps radius of ball on which solution is computed
     * @param k wavenuber
     * @param n_i refraction index of inner domain
     * @return transmission coefficient of <tt>n</tt>-th product in series expansion of the solution
     * given radius \p r, wavenumber \p k and refraction index \p n_i
     */
    std::complex<double> t_coeff(int n,
                                 double eps,
                                 double k,
                                 double n_i);

    /**
     * Evaluates an incoming wave for the wavenumber \p k at point \p x1, \p x2 with the
     * \f$ 2\verb|l| + 1 \f$ series expansion coefficients stored in a_n.
     * @param x1 x-value of evaluation point
     * @param x2 y-value of evaluation point
     * @param l number of coefficients for series expansion
     * @param a_n coefficients fo rseries expansion
     * @param k wavenumber
     * @return value of incoming wave given the wavenumber \p k and the series expansion
     * defined by \p a_n and \p l at the point \p x1 \p x2
     */
    std::complex<double> u_i(double x1,
                             double x2,
                             int l,
                             double *a_n,
                             double k);

    /**
     * Evaluates a wave scattered by a 2-dimensional ball of radius \p eps
     * for the wavenumber \p k at point \p x1, \p x2 with the
     * \f$ 2\verb|l| + 1 \f$ series expansion coefficients stored in a_n.
     * The refraction index of the ball is given by n_i.
     * @param x1 x-value of evaluation point
     * @param x2 y-value of evaluation point
     * @param l number of coefficients for series expansion
     * @param eps radius of ball
     * @param a_n coefficients fo rseries expansion
     * @param k wavenumber
     * @param n_i refraction index of ball
     * @return value of scattered wave given the wavenumber \p k, the series expansion
     * defined by \p a_n and \p l, and the scattering ball of radius \p eps at the point \p x1 \p x2
     */
    std::complex<double> u_s(double x1,
                             double x2,
                             int l,
                             double eps,
                             double *a_n,
                             double k,
                             double n_i);

    /**
     * Evaluates a wave transmitted through a 2-dimensional ball of radius \p eps
     * for the wavenumber \p k at point \p x1, \p x2 with the
     * \f$ 2\verb|l|+1 \f$ series expansion coefficients stored in a_n.
     * The refraction index of the ball is given by n_i.
     * @param x1 x-value of evaluation point
     * @param x2 y-value of evaluation point
     * @param l number of coefficients for series expansion
     * @param eps radius of ball
     * @param a_n coefficients fo rseries expansion
     * @param k wavenumber
     * @param n_i refraction index of ball
     * @return value of transmitted wave given the wavenumber \p k, the series expansion
     * defined by \p a_n and \p l, and the scattering ball of radius \p eps at the point \p x1 \p x2
     */
    std::complex<double> u_t(double x1,
                             double x2,
                             int l,
                             double eps,
                             double *a_n,
                             double k,
                             double n_i);

    /**
     * Evaluates Neumann data of an incoming wave for the wavenumber \p k at point \p x1, \p x2 with the
     * \f$ 2\verb|l| + 1 \f$ series expansion coefficients stored in a_n.
     * @param x1 x-value of evaluation point
     * @param x2 y-value of evaluation point
     * @param l number of coefficients for series expansion
     * @param a_n coefficients fo rseries expansion
     * @param k wavenumber
     * @return Neumann data of incoming wave given the wavenumber \p k and the series expansion
     * defined by \p a_n and \p l at the point \p x1 \p x2
     */
    std::complex<double> u_i_neu(double x1,
                                 double x2,
                                 int l,
                                 double *a_n,
                                 double k);

    /**
     * Evaluates Neumann data of a wave scattered by a 2-dimensional ball of radius \p eps
     * for the wavenumber \p k at point \p x1, \p x2 with the
     * \f$ 2\verb|l| + 1 \f$ series expansion coefficients stored in a_n.
     * The refraction index of the ball is given by n_i.
     * @param x1 x-value of evaluation point
     * @param x2 y-value of evaluation point
     * @param l number of coefficients for series expansion
     * @param eps radius of ball
     * @param a_n coefficients fo rseries expansion
     * @param k wavenumber
     * @param n_i refraction index of ball
     * @return Neumann data of scattered wave given the wavenumber \p k, the series expansion
     * defined by \p a_n and \p l, and the scattering ball of radius \p eps at the point \p x1 \p x2
     */
    std::complex<double> u_s_neu(double x1,
                                 double x2,
                                 int l,
                                 double eps,
                                 double *a_n,
                                 double k,
                                 double n_i);

    /**
     * Evaluates Neumann data of a wave transmitted through a 2-dimensional ball of radius \p eps
     * for the wavenumber \p k at point \p x1, \p x2 with the
     * \f$ 2\verb|l|+1 \f$ series expansion coefficients stored in a_n.
     * The refraction index of the ball is given by n_i.
     * @param x1 x-value of evaluation point
     * @param x2 y-value of evaluation point
     * @param l number of coefficients for series expansion
     * @param eps radius of ball
     * @param a_n coefficients fo rseries expansion
     * @param k wavenumber
     * @param n_i refraction index of ball
     * @return Neumann data of transmitted wave given the wavenumber \p k, the series expansion
     * defined by \p a_n and \p l, and the scattering ball of radius \p eps at the point \p x1 \p x2
     */
    std::complex<double> u_t_neu(double x1,
                                 double x2,
                                 int l,
                                 double eps,
                                 double *a_n,
                                 double k,
                                 double n_i);
}
#endif //ADVNUMCSE_GENERATE_SOLUTION_H
