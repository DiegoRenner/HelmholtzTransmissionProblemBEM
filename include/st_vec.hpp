//
// Created by diego on 11/26/20.
//

#ifndef HELMHOLTZTRANSMISSIONPROBLEM_ST_VEC_HPP
#define HELMHOLTZTRANSMISSIONPROBLEM_ST_VEC_HPP


#include <complex>
#include <Eigen/Dense>

class StVec {
public:
    void set_k(std::complex<double> k);
    void set_vec(std::complex<double>* vec);
    std::complex<double> get_k();
    std::complex<double>* get_vec();
    StVec(int dim,
          std::complex<double> k = std::complex<double>(0.0,0.0),
          std::complex<double>* vec = NULL);
    bool operator < (const StVec& str) const
    {
        return ( k.real() < str.k.real());
    }
private:
    std::complex<double> k;
    std::complex<double>* vec = NULL;
    int dim;

};


#endif //HELMHOLTZTRANSMISSIONPROBLEM_ST_VEC_HPP
