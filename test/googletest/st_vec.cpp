/**
 * \file hypersingular_test.cpp
 * \brief This test file compares the Hypersingular BIO for the
 * Helmholtz kernel to a precomputed known solution from
 * a file.
 */

#include <complex>
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include "st_vec_storage.hpp"


typedef std::complex<double> complex_t;


TEST(ST_VEC_TEST, check_consistency) {

    unsigned size = 1000;
    StVecStorage storage(size, 0, NULL);
    complex_t* vector_in = new complex_t[size];
    for (unsigned i = 0; i<size; i++){
        vector_in[i] = complex_t(std::rand(),std::rand());
    }

    storage.add_vec(0.0, vector_in);
    complex_t* vector_out = new complex_t[size];
    vector_out = storage.get_closest(complex_t(0.0,0.0));
    for (unsigned i = 0; i<size; i++){
        ASSERT_EQ(vector_in[i],vector_out[i]);
    }

    for (unsigned i = 0; i<size; i++){
        vector_out[i] = complex_t(std::rand(),std::rand());
    }
    for (unsigned i = 0; i<size; i++){
        ASSERT_NE(vector_in[i], vector_out[i]);
    }

};

