//
// Created by diegorenner on 4/29/20.
//

#ifndef HELMHOLTZTRANSMISSIONPROBLEM_ROOTS_HPP
#define HELMHOLTZTRANSMISSIONPROBLEM_ROOTS_HPP

#include <Eigen/Dense>

    double zbrent( const std::function<double(double)> f,
                   double x1,
                   double x2,
                   double tol);

    double rtsafe( const std::function<Eigen::MatrixXd(double)> f,
                   double x1,
                   double x2,
                   double tol,
                   bool &root_found,
                   unsigned &num_iter);

    Eigen::VectorXd parabolic_approximation(const std::function<Eigen::VectorXd(double)> f,
                                            const std::function<Eigen::VectorXd(double)> f_der,
                                            const std::function<Eigen::VectorXd(double)> f_der2,
                                            const double x0,
                                            double step);

    double secant_method( std::function<double(double)> f,
            double x0,
            double x1,
            const double tol,
            const unsigned maxIter,
            bool &root_found,
            unsigned &num_iter);

#endif //HELMHOLTZTRANSMISSIONPROBLEM_ROOTS_HPP
