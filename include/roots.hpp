//
// Created by diegorenner on 4/29/20.
//

#ifndef HELMHOLTZTRANSMISSIONPROBLEM_FIND_ZEROS_HPP
#define HELMHOLTZTRANSMISSIONPROBLEM_FIND_ZEROS_HPP

#include <Eigen/Dense>

namespace parametricbem2d {
    double zbrent( const std::function<double(double)> f,
                   double x1,
                   double x2,
                   double tol);

    double rtsafe( const std::function<Eigen::MatrixXd(double)> f,
                   double x1,
                   double x2,
                   double tol);

    double secant_method( const std::function<double(double)> f,
            double x0,
            double x1,
            const double tol,
            const unsigned maxIter);

}
#endif //HELMHOLTZTRANSMISSIONPROBLEM_FIND_ZEROS_HPP
