//
// Created by diego on 11/26/20.
//

#include <iostream>
#include "st_vec_storage.hpp"
typedef std::complex<double> complex_t ;

StVecStorage::StVecStorage(int dim, complex_t k , complex_t* vec){
    this->dim = dim;
    st_vecs = std::vector<StVec>();
    StVec first_vec(dim, k);

    double norm = 0;
    for (int i = 0; i < this->dim; i++) {
        first_vec.get_vec()[i] = complex_t(2*(std::rand()%10000)/10000.0-1,2*(std::rand()%10000)/10000.0-1);
        norm += (first_vec.get_vec()[i]*std::conj(first_vec.get_vec()[i])).real();
    }
    norm = std::sqrt(norm);
    for (int i = 0; i < this->dim; i++){
        first_vec.get_vec()[i] /= norm;
    }
    st_vecs.push_back(first_vec);
}

void StVecStorage::add_vec(complex_t k, complex_t* vec){
    StVec new_vec(this->dim,k,vec);
    st_vecs.push_back(new_vec);
    std::sort(st_vecs.begin(), st_vecs.end());
}

std::complex<double>* StVecStorage::get_closest(std::complex<double> k) {
    int iter = 0;
    for ( int i = 0; i < st_vecs.size(); i++){
        if ( this->st_vecs[i].get_k().real() > k.real()){
            iter = i+1;
            break;
        }
        if (i == st_vecs.size()-1){
            iter = i;
            break;
        }
    }
    return this->st_vecs[iter].get_vec();
}

void StVecStorage::status(){
    for ( int i = 0; i < st_vecs.size(); i++){
        std::cout << st_vecs[i].get_k() << std::endl;
        if (st_vecs[i].get_vec() != NULL){
            std::cout << "************************************" << std::endl;
            for( unsigned j = 0; j < this->dim; j++){
                std::cout << st_vecs[i].get_vec()[j];
            }
            std::cout << std::endl;
        }
    }
}
std::vector<StVec> StVecStorage::get_st_vecs(){
    return st_vecs;
}

StVecStorage::~StVecStorage() {
    for (unsigned i = 0; i < this->st_vecs.size(); i++){
        delete this->st_vecs[i].get_vec();
    }
}
