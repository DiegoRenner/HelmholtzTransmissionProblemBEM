// // Created by diego on 11/26/20. //

#include "st_vec.hpp"
#include <Eigen/Dense>
#include <vector>
typedef std::complex<double> complex_t;

void StVec::set_k(complex_t k){
    this->k = k;
}
void StVec::set_vec(complex_t* vec){
    this->vec = vec;
}
complex_t StVec::get_k(){
    return this->k;
}
complex_t* StVec::get_vec(){
    return this->vec;
}
StVec::StVec(int dim, complex_t k, complex_t* vec){
    this->dim = dim;
    this->k = k;
    this->vec = new complex_t [dim];
    if(vec != NULL) {
        for (int i = 0; i < dim; i++) {
            this->vec[i] = vec[i];
        }
    }
}


