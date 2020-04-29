//
// Created by diegorenner on 4/29/20.
//

#ifndef HELMHOLTZTRANSMISSIONPROBLEM_FIND_ZEROS_HPP
#define HELMHOLTZTRANSMISSIONPROBLEM_FIND_ZEROS_HPP

namespace parametricbem2d {
    double secant_method( const std::function<double(double)> f,
            double x0,
            double x1,
            const double tol,
            const unsigned maxIter);
}
#endif //HELMHOLTZTRANSMISSIONPROBLEM_FIND_ZEROS_HPP
