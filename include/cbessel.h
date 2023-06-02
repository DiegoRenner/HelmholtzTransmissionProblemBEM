// Compute Bessel functions of real order and complex argument.
// SPDX-FileCopyrightText: 2023 Luka Marohnić <lmarohni@tvz.hr>
// SPDX-License-Identifier: LGPL-3-or-later

#ifndef BESSEL_H
#define BESSEL_H

#include <complex>

using namespace std;

/**
 * This is a C++ library for Bessel functions of complex argument and
 * real order, intended as a replacement for complex_bessel library
 * using the Fortran code by Donald E. Amos.
 *
 * The same theoretical basis is used, as explained in the two papers
 * by Amos (1983). This code includes optimizations for order 0 and 1,
 * making it slightly faster than the Fortran library.
 *
 * The routines below will throw exceptions if CBESSEL_EXCEPT is defined
 * (by default it is not defined, the respective line in the CPP file
 * is commented out). The following exceptions are thrown:
 *
 *  - std::overflow_error in the case of overflow
 *  - std::underflow_error in the case of underflow (the result is 0)
 *  - std::invalid_argument in the case of undefined result (NAN)
 *  - std::length_error in the case of convergence failure
 *  - std::domain_error in the case of invalid parameter value
 */
namespace ComplexBessel
{

    typedef double Real;
    typedef complex<Real> Cplx;
    typedef struct olver_data {
        bool is_valid;
        Cplx S1,S2,xi,phi; // xi = 2/3*zeta^(3/2)
        olver_data();
    } OlverData;

    /* Bessel functions and their derivatives */
    Cplx I  (Real v,const Cplx &z,bool scaled=false);
    Cplx Ip (Real v,const Cplx &z,int n=1);
    Cplx J  (Real v,const Cplx &z,bool scaled=false);
    Cplx Jp (Real v,const Cplx &z,int n=1);
    Cplx K  (Real v,const Cplx &z,bool scaled=false);
    Cplx Kp (Real v,const Cplx &z,int n=1);
    Cplx Y  (Real v,const Cplx &z,bool scaled=false);
    Cplx Yp (Real v,const Cplx &z,int n=1);

    /* Hankel functions and their derivatives */
    Cplx H1 (Real v,const Cplx &z,bool scaled=false);
    Cplx H1p(Real v,const Cplx &z,int n=1);
    Cplx H2 (Real v,const Cplx &z,bool scaled=false);
    Cplx H2p(Real v,const Cplx &z,int n=1);

    /* Airy functions and their first derivatives */
    Cplx Ai (       const Cplx &z,bool scaled=false);
    Cplx Aip(       const Cplx &z,bool scaled=false);
    Cplx Bi (       const Cplx &z,bool scaled=false);
    Cplx Bip(       const Cplx &z,bool scaled=false);

}

#endif // BESSEL_H
