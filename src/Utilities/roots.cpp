//
// Created by diegorenner on 4/29/20.
//

#include <limits>
#include <functional>
#include <iostream>
#include <cmath>
#include <fstream>
#include "roots.hpp"

#define ITMAX 100
#define MAXIT 100
#define EPS std::numeric_limits<double>::epsilon()


    using namespace std;
    // Using Brent’s method, find the root of a function func known to lie between x1 and x2.
    // The root, returned as zbrent , will be refined until its accuracy is tol .
    double zbrent( const function<double(double)> f,
                    double x1,
                    double x2,
                    double tol,
                    bool &root_found,
                    unsigned &num_iter){
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
                root_found = true;
                num_iter = iter;
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

    double rtsafe( std::function<Eigen::MatrixXd(double)> fct,
                   double x1,
                   double x2,
                   double tol,
                   bool &root_found,
                   unsigned &num_iter){
        int j;
        double df,dx,dxold,f,fh,fl;
        double temp,xh,xl,rts;
        Eigen::MatrixXd tempM = fct(x1);
        fl = tempM(0,0);
        df = tempM(0,1);
        tempM = fct(x2);
        fh = tempM(0,0);
        df = tempM(0,1);
        if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
            std::cout << "Root must be bracketed in rtsafe" << std::endl;
            if (fl > fh){
                return x1;
            }else{
                return x2;
            };
        }
        if (fl == 0.0) {
            root_found = true;
            return x1;
        }
        if (fh == 0.0) {
            root_found = true;
            return x2;
        }
        if (fl < 0.0) {
            //Orient the search so that f (xl) < 0.
            xl=x1;
            xh=x2;
        } else {
            xh=x1;
            xl=x2;
        }
        rts=0.5*(x1+x2);
        // Initialize the guess for root,the “stepsize before last,”and the last step.
        dxold=fabs(x2-x1);
        dx=dxold;
        tempM = fct(rts);
        f = tempM(0,0);
        df = tempM(0,1);
        for (j=1;j<=MAXIT;j++) {
            //Loop over allowed iterations. || Bisect if Newton out of range,or not decreasing fast enough.
            if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) || (fabs(2.0*f) > fabs(dxold*df))) {
                dxold=dx;
                dx=0.5*(xh-xl);
                rts=xl+dx;
                //Change in root is negligible.
                if (xl == rts) {
                    root_found = true;
                    num_iter = j;
                    return rts;
                }
            } else {
                //Newton step acceptable. Take it.
                dxold=dx;
                dx=f/df;
                temp=rts;
                rts -= dx;
                if (temp == rts) {
                    root_found = true;
                    num_iter = j;
                    return rts;
                }
            }
            if (fabs(dx) < tol) {
                root_found = true;
                num_iter = j;
                return rts;
            }
            //Convergence criterion.
            tempM = fct(rts);
            f = tempM(0,0);
            df = tempM(0,1);
            //The one new function evaluation per iteration.
            //Maintain the bracket on the root.
            if (f < 0.0)
                xl=rts;
            else
                xh=rts;
        }
        std::cout << "Maximum number of iterations exceeded in rtsafe" << std::endl;
        return 0.0;
        //Never get here.
    };
    Eigen::VectorXd parabolic_approximation(const std::function<Eigen::VectorXd(double)> f,
                                            const std::function<Eigen::VectorXd(double)> f_der,
                                            const std::function<Eigen::VectorXd(double)> f_der2,
                                            const double x0,
                                            double step){
        Eigen::VectorXd vals = f(x0);
        Eigen::VectorXd ders = f_der(x0);
        Eigen::VectorXd ders2 = f_der2(x0);
        unsigned numsvs = vals.size();
        Eigen::VectorXd exts(numsvs);
        bool adjust_step;
        step = step/(1+ders2.cwiseAbs().maxCoeff());
        Eigen::VectorXd res(4);
        do {
            adjust_step = false;
            for (unsigned i = 0; i < numsvs; i++) {
                exts[i] = (-ders[i]/ders2[i] + x0);
                if ( exts[i] - x0 > step) {
                    exts[i] = x0 + step;
                } else if( exts[i] - x0 < 0 ){
                    if (ders2[i]<0){
                        exts[i] = x0+step;
                    } else if (exts[i] -x0 < -step){
                        exts[i] = x0 - 0.5 * step;
                    }
                }
            }
            res[0] = exts[0];
            res[1] = f(res[0])[0];
            double temp;
            unsigned index = 0;
            for (unsigned i = 1; i < numsvs; i++) {
                temp = f(exts[i])[i];
                if (temp < res[1]) {
                    res[1] = temp;
                    index = i;
                }
            }
            res[0] = exts[index];
            res[1] = vals[index];
            res[2] = ders[index];
            res[3] = ders2[index];
            if (res[1] < 0) {
                adjust_step = true;
                step *= 0.5;
            }
        } while(adjust_step);
       return res;

    }

    double secant_method(const function<double(double)> f,
                         double x1,
                         double x2,
                         const double tol,
                         const unsigned maxIter,
                         bool &root_found,
                         unsigned &num_iter) {
        unsigned n = 0;
        double xm = 0.0;
        double x0 = 0.0;
        double c = 0.0;
        double f2 = f(x2);
        double f1 = f(x1);
        if (f1 * f2 < 0) {
            do {
                // calculate the intermediate value
                x0 = (x1 * f2 - x2 * f1) / (f2 - f1);
                double f0 = f(x0);

                // check if x0 is root of equation or not
                c = f1 * f0;

                // update the value of interval
                x1 = x2;
                x2 = x0;
                f1 = f(x1);
                f2 = f(x2);
                // update number of iteration
                n++;
                if (abs(c) < tol)
                    break;

                // if x0 is the root of equation then break the loop
                xm = (x1 * f2 - x2 * f1) / (f2 - f1);
            } while (abs(xm - x0) >= tol); // repeat the loop
            // until the convergence

            if (abs(c) < tol) {
                cout << "Root of the given equation = " << x0 << endl;
                cout << "No. of iterations = " << n << endl;
                root_found = true;
                num_iter = n;
                return x0;
            } else {
                cout << "No root found after " << n << " iterations." << endl;
            }
        } else
            cout << "There might not be a root in this interval." << endl;
        return 0.0;
    }
