/**
 * \file cbessel.hpp
 *
 * \brief This is a C++ library for Bessel functions of complex argument and
 *        real order, intended as a replacement for complex_bessel library
 *        which uses the Fortran code by Donald E. Amos.
 *
 * The same theoretical basis is used, as explained in the two papers by
 * Amos (1983). This code includes optimizations for order 0 and 1, making it
 * slightly faster than the Fortran library. The routines below will throw
 * exceptions if CBESSEL_EXCEPT is defined (by default it is not defined,
 * the respective line in the CPP file is commented out).
 *
 * The following exceptions are thrown:
 *  - std::overflow_error in the case of overflow
 *  - std::underflow_error in the case of underflow (the result is 0)
 *  - std::invalid_argument in the case of undefined result (NAN)
 *  - std::length_error in the case of convergence failure
 *  - std::domain_error in the case of invalid parameter value
 * Throwing exceptions is more informative but slightly slower.
 * Convergence problems should not happen in practice.
 *
 * (c) 2023 Luka Marohnić
 */

#ifndef CBESSEL_HPP
#define CBESSEL_HPP

#include <complex>

using namespace std;

namespace complex_bessel {

    /**
     * typedef for double precision real type
     */
    typedef double Real;

    /**
     * typedef for double precision complex type
     */
    typedef complex<Real> Cplx;

    /**
     * workspace for Olver expansions for large order
     */
    typedef struct olver_data {
        bool is_valid;
        Cplx S1,S2,xi,phi; // xi = 2/3*zeta^(3/2)
        olver_data();
    } OlverData;

    /**
     * This function computes the value of the modified
     * Bessel function of the first kind.
     *
     * @param v real order
     * @param z complex argument
     * @param scaled if true, the value is scaled by exp(-real(z))
     */
    Cplx I  (Real v,const Cplx &z,bool scaled=false);

    /**
     * This function computes the nth derivative of the modified
     * Bessel function of the first kind.
     *
     * @param v real order
     * @param z complex argument
     * @param n positive integer (the order of derivative)
     */
    Cplx Ip (Real v,const Cplx &z,int n=1);

    /**
     * This function computes the value of the
     * Bessel function of the first kind.
     *
     * @param v real order
     * @param z complex argument
     * @param scaled if true, the value is scaled by exp(-imag(z))
     */
    Cplx J  (Real v,const Cplx &z,bool scaled=false);

    /**
     * This function computes the nth derivative of the
     * Bessel function of the first kind.
     *
     * @param v real order
     * @param z complex argument
     * @param n positive integer (the order of derivative)
     */
    Cplx Jp (Real v,const Cplx &z,int n=1);

    /**
     * This function computes the value of the modified
     * Bessel function of the second kind.
     *
     * @param v real order
     * @param z complex argument
     * @param scaled if true, the value is scaled by exp(z)
     */
    Cplx K  (Real v,const Cplx &z,bool scaled=false);

    /**
     * This function computes the nth derivative of the modified
     * Bessel function of the second kind.
     *
     * @param v real order
     * @param z complex argument
     * @param n positive integer (the order of derivative)
     */
    Cplx Kp (Real v,const Cplx &z,int n=1);

    /**
     * This function computes the value of the
     * Bessel function of the second kind.
     *
     * @param v real order
     * @param z complex argument
     * @param scaled if true, the value is scaled by exp(-imag(z))
     */
    Cplx Y  (Real v,const Cplx &z,bool scaled=false);

    /**
     * This function computes the nth derivative of the
     * Bessel function of the second kind.
     *
     * @param v real order
     * @param z complex argument
     * @param n positive integer (the order of derivative)
     */
    Cplx Yp (Real v,const Cplx &z,int n=1);

    /**
     * This function computes the value of the
     * Hankel function of the first kind.
     *
     * @param v real order
     * @param z complex argument
     * @param scaled if true, the value is scaled by exp(-i*z)
     */
    Cplx H1 (Real v,const Cplx &z,bool scaled=false);

    /**
     * This function computes the nth derivative of the
     * Hankel function of the first kind.
     *
     * @param v real order
     * @param z complex argument
     * @param n positive integer (the order of derivative)
     */
    Cplx H1p(Real v,const Cplx &z,int n=1);

    /**
     * This function computes the value of the
     * Hankel function of the second kind.
     *
     * @param v real order
     * @param z complex argument
     * @param scaled if true, the value is scaled by exp(i*z)
     */
    Cplx H2 (Real v,const Cplx &z,bool scaled=false);

    /**
     * This function computes the nth derivative of the
     * Hankel function of the second kind.
     *
     * @param v real order
     * @param z complex argument
     * @param n positive integer (the order of derivative)
     */
    Cplx H2p(Real v,const Cplx &z,int n=1);

    /**
     * This function computes the value of the
     * Airy function of the first kind.
     *
     * @param z complex argument
     * @param scaled if true, the result is scaled by exp(2/3*z*sqrt(z))
     */
    Cplx Ai (       const Cplx &z,bool scaled=false);

    /**
     * This function computes the first derivative of the
     * Airy function of the first kind.
     *
     * @param z complex argument
     * @param scaled if true, the result is scaled by exp(2/3*z*sqrt(z))
     */
    Cplx Aip(       const Cplx &z,bool scaled=false);

    /**
     * This function computes the value of the
     * Airy function of the second kind.
     *
     * @param z complex argument
     * @param scaled if true, the result is scaled by exp(-abs(real(2/3*z*sqrt(z))))
     */
    Cplx Bi (       const Cplx &z,bool scaled=false);

    /**
     * This function computes the first derivative of the
     * Airy function of the second kind.
     *
     * @param z complex argument
     * @param scaled if true, the result is scaled by exp(-abs(real(2/3*z*sqrt(z))))
     */
    Cplx Bip(       const Cplx &z,bool scaled=false);

}

#endif // CBESSEL_HPP
