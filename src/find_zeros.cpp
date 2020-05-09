//
// Created by diegorenner on 4/29/20.
//

#include <limits>
#include <functional>
#include <iostream>
#include <cmath>
#include <fstream>
#include "find_zeros.hpp"

#define ITMAX 100
#define EPS std::numeric_limits<double>::epsilon()


namespace parametricbem2d {
    using namespace std;
    // Using Brent’s method, find the root of a function func known to lie between x1 and x2.
    // The root, returned as zbrent , will be refined until its accuracy is tol .
    double zbrent( const function<double(double)> f,
                    double x1,
                    double x2,
                    double tol){
        int iter;
        double a=x1, b=x2, c=x2, d, e, min1, min2;
        double fa=f(a), fb=f(b), fc, p, q, r, s, tol1, xm;
        if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
            cout << "Root must be bracketed in zbrent" << endl;
            return 0.0;
        }
        fc=fb;
        for (iter=1; iter <= ITMAX; iter++) {
            if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
                c = a; //Rename a, b, c, and adjust bounding interval d.
                fc=fa;
                e=d=b-a;
            }
            if (fabs(fc) < fabs(fb)) {
                a=b;
                b=c;
                c=a;
                fa=fb;
                fb=fc;
                fc=fa;
            }
            tol1=2.0*EPS*fabs(b)+0.5*tol; //Convergence Check
            xm=0.5*(c-b);
            if (fabs(xm) <= tol1 || fb == 0.0){
                std::ofstream filename;
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/SV_analysis_roots.dat", std::ios::app);
                filename << 0.0 << " " << b << " " << f(b) << std::endl;
                filename.close();
                return b;
            }
            if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
                s=fb/fa; //Attempt inverse quadratic interpolation.
                if (a == c) {
                    p=2.0*xm*s;
                    q=1.0-s;
                } else {
                    q=fa/fc;
                    r=fb/fc;
                    p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                    q=(q-1.0)*(r-1.0)*(s-1.0);
                }
                if (p > 0.0) q = -q; //Check whether in bounds
                p=fabs(p);
                min1=3.0*xm*q-fabs(tol1*q);
                min2=fabs(e*q);
                if(2.0*p < (min1 < min2 ? min1 : min2)) {
                    e=d; //Accept interpolation
                    d=p/q;
                } else {
                    d=xm; //Interpolation failed, use bisection.
                    e=d;
                }
            } else {
                d = xm; //Bounds decreasing too slowly, use bisection.
                e = d;
            }
            a=b; //Move last best guess to a.
            fa=fb;
            if (fabs(d) > tol1) //Evaluate new trial root.
                b += d;
            else
                b += (xm >= 0) ? fabs(tol1) : -fabs(tol1);
            fb = f(b);
        }
        cout << "Maximum number of iterations exceeded in zbrent" << endl;
        return 0.0; //Never get here.
    }

    double secant_method(const function<double(double)> f,
                         double x1,
                         double x2,
                         const double tol,
                         const unsigned maxIter) {
        unsigned n = 0;
        double xm = 0.0;
        double x0 = 0.0;
        double c = 0.0;
        if (f(x1) * f(x2) < 0) {
            do {
                // calculate the intermediate value
                x0 = (x1 * f(x2) - x2 * f(x1)) / (f(x2) - f(x1));

                // check if x0 is root of equation or not
                c = f(x1) * f(x0);

                // update the value of interval
                x1 = x2;
                x2 = x0;

                // update number of iteration
                n++;

                // if x0 is the root of equation then break the loop
                if (abs(c) < tol)
                    break;
                xm = (x1 * f(x2) - x2 * f(x1)) / (f(x2) - f(x1));
            } while (abs(xm - x0) >= tol); // repeat the loop
            // until the convergence

            if (abs(c) < tol) {
                cout << "Root of the given equation = " << x0 << endl;
                cout << "No. of iterations = " << n << endl;
                std::ofstream filename;
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/SV_analysis_roots.dat", std::ios::app);
                filename << 0.0 << " " << x0 << std::endl;
                filename.close();
            } else {
                cout << "No root found after " << n << " iterations." << endl;
            }
        } else
            cout << "There might not be a root in this interval." << endl;
        return x0;
    }
}