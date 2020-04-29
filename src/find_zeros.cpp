//
// Created by diegorenner on 4/29/20.
//

#include <limits>
#include <functional>
#include "find_zeros.hpp"


namespace parametricbem2d {
    using namespace std;

    double epsilon = numeric_limits<double>::epsilon();
    double secant_method(const function<double(double)> f,
                         double x1,
                         double x2,
                         const double tol,
                         const unsigned maxIter) {
        unsigned n;
        double xm, x0, c;
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
                if (abs(c) < epsilon)
                    break;
                xm = (x1 * f(x2) - x2 * f(x1)) / (f(x2) - f(x1));
            } while (fabs(xm - x0) >= E); // repeat the loop
            // until the convergence

            if (abs(c) < epsilon) {
                cout << "Root of the given equation=" << x0 << endl;
                cout << "No. of iterations = " << n << endl;
            } else {
                cout << "No root found after" << n << " iterations." << endl;
            }
        } else
            cout << "There might not be a root in this interval." << endl;
    }
}