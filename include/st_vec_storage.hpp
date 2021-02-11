//
// Created by diego on 11/26/20.
//

#ifndef HELMHOLTZTRANSMISSIONPROBLEM_ST_VEC_STORAGE_HPP
#define HELMHOLTZTRANSMISSIONPROBLEM_ST_VEC_STORAGE_HPP


#include <complex>
#include <vector>
#include "st_vec.hpp"


class StVecStorage {
public:
    StVecStorage(int dim,
                 std::complex<double> k = std::complex<double>(0.0,0.0),
                 std::complex<double>* vec = NULL);
    void add_vec(std::complex<double> k,
                 std::complex<double>* vec);
    std::complex<double>* get_closest(std::complex<double> k);
    std::vector<StVec> get_st_vecs();
    void status();
    ~StVecStorage();

private:
    std::vector<StVec> st_vecs;
    int dim = 0;
};

#endif //HELMHOLTZTRANSMISSIONPROBLEM_ST_VEC_STORAGE_HPP
