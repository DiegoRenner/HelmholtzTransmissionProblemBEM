#include <limits>
#include <functional>
#include <iostream>
#include <ostream>
#include <cmath>
#include "find_roots.hpp"

#define MAXIT 10
#define EPS std::numeric_limits<double>::epsilon()


using namespace std;
double zbrent( const function<double(double)> f,
               double x1,
               double x2,
               double tol,
               bool &root_found,
               unsigned &num_iter){
    // counter for iterations
    int iter;
    // initialize function values and boundaries
    double a=x1, b=x2, c=x2, d, e, min1, min2;
    double fa=f(a), fb=f(b), fc, p, q, r, s, tol1, xm;
    // sanity checks
    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
		#ifdef CMDL
        cout << "Root must be bracketed in zbrent" << endl;
		#endif
        return 0.0;
    }
    fc=fb;
    for (iter=1; iter <= MAXIT; iter++) {
        // reorient boundary for next interpolation
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c = a;
            fc=fa;
            e=d=b-a;
        }
        if (fabs(fc) < fabs(fb)){
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        }
        // check if converged
        tol1=2.0*EPS*fabs(b)+0.5*tol;
        xm=0.5*(c-b);
        if (fabs(xm) <= tol1 || fb == 0.0){
            root_found = true;
            num_iter = iter;
            return b;
        }
        // try quadratic interpolation if
        // bounds are decreasing
        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
            // compute parameters for quadratic interpolation
            s=fb/fa;
            if (a == c) {
                p=2.0*xm*s;
                q=1.0-s;
            } else {
                q=fa/fc;
                r=fb/fc;
                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                q=(q-1.0)*(r-1.0)*(s-1.0);
            }
            // check if quadratic interpolation
            // would fall within bounds
            if (p > 0.0) q = -q;
            p=fabs(p);
            min1=3.0*xm*q-fabs(tol1*q);
            min2=fabs(e*q);
            if(2.0*p < (min1 < min2 ? min1 : min2)) {
                // interpolation accepted
                e=d;
                d=p/q;
            } else {
                // interpolation rejected, bisect instead
                d=xm;
                e=d;
            }
        } else {
            // convergence too slow, bounds not collapsing
            // fast enough, bisect
            d = xm;
            e = d;
        }
        // store previous best guess before
        // computing new best guess
        a=b;
        fa=fb;
        // compute new best guess
        // if step taken is too small, take
        // minimum accepted step towards 0
        if (fabs(d) > tol1)
            b += d;
        else
            b += (xm >= 0) ? fabs(tol1) : -fabs(tol1);
        // one new function evaluation per iteration
        fb = f(b);
    }
	#ifdef CMDL
    cout << "Maximum number of iterations exceeded in zbrent" << endl;
	#endif
    // should never be reached
    return 0.0;
}

double rtsafe( std::function<double(double)> fct,
               std::function<Eigen::MatrixXd(double)> fct_both,
               double x1,
               double x2,
               double tol,
               bool &root_found,
               unsigned &num_iter){
    // initialize counter, function values and boundaries
    int j;
    double df,dx,dxold,f,fh,fl;
    double temp,xh,xl,rts;
    fl = fct(x1);
    fh = fct(x2);
    // sanity checks
    if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
		#ifdef CMDL
        std::cout << "Root must be bracketed in rtsafe" << std::endl;
		#endif
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
    // assign boundaries s.t f(xl) < 0
    if (fl < 0.0) {
        xl=x1;
        xh=x2;
    } else {
        xh=x1;
        xl=x2;
    }
    // set first guess for root, lest step, and "step before last"
    rts=0.5*(x1+x2);
    dxold=fabs(x2-x1);
    dx=dxold;
    // initialize functions and derivative value
    Eigen::MatrixXd tempM = fct_both(rts);
    f = tempM(0,0);
    df = tempM(0,1);
    // loop over max number of allowed iterations
    for (j=1;j<=MAXIT;j++) {
        // check if newton step is out of range or too slow
        // if it is, bisect
        if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) || (fabs(2.0*f) > fabs(dxold*df))) {
            dxold=dx;
            dx=0.5*(xh-xl);
            rts=xl+dx;
            // if change in root is negligible assume convergence
            if (xl == rts) {
                root_found = true;
                num_iter = j;
                return rts;
            }
        } else {
            // newton step accepted
            dxold=dx;
            dx=f/df;
            temp=rts;
            rts -= dx;
            // if change in root is negligible assume convergence
            if (temp == rts) {
                root_found = true;
                num_iter = j;
                return rts;
            }
        }
        // if change in root is small enough assume convergence
        if (fabs(dx) < tol) {
            root_found = true;
            num_iter = j;
            return rts;
        }
        // one new function/derivative evaluation per iteration
        tempM = fct_both(rts);
        f = tempM(0,0);
        df = tempM(0,1);
        // assign new boundary making sure root stay bracketed
        if (f < 0.0)
            xl=rts;
        else
            xh=rts;
    }
	#ifdef CMDL
    std::cout << "Maximum number of iterations exceeded in rtsafe" << std::endl;
	#endif
    // should never be reached
    return 0.0;
};
Eigen::VectorXd parabolic_approximation(const std::function<Eigen::VectorXd(double)> f,
                                        const std::function<Eigen::VectorXd(double)> f_der,
                                        const std::function<Eigen::VectorXd(double)> f_der2,
                                        const double x0,
                                        double step){
    // initialize function values and derivatives
    Eigen::VectorXd vals = f(x0);
    Eigen::VectorXd ders = f_der(x0);
    Eigen::VectorXd ders2 = f_der2(x0);
    // initialize vecotor for minima approximations
    unsigned numsvs = vals.size();
    Eigen::VectorXd exts(numsvs);
    //bool adjust_step;
    //step = step/(1+ders2.cwiseAbs().maxCoeff());
    Eigen::VectorXd res(4);
    // do {
    //adjust_step = false;
    // compute approximations for next minima
    for (unsigned i = 0; i < numsvs; i++) {
        // find minima of parabola
        exts[i] = (-ders[i]/ders2[i] + x0);
        // try and stop convergence to maxima
        if ( exts[i] - x0 > step) {
            exts[i] = x0 + step;
        }else if (exts[i] -x0 < -step){
            exts[i] = x0 - step;
        }/* else if( exts[i] - x0 < 0 ){
            if (ders2[i]<0){
                exts[i] = x0+step;
            } else if (exts[i] -x0 < -step){
                exts[i] = x0 - 0.5 * step;
            }
        }*/
    }
    // find best approximation for minima
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
    // return best approx for minima and
    // values used for parabolic approximation
    res[0] = exts[index];
    res[1] = vals[index];
    res[2] = ders[index];
    res[3] = ders2[index];

    /*if (res[1] < 0) {
        adjust_step = true;
        step *= 0.5;
    }
} while(adjust_step);*/
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
