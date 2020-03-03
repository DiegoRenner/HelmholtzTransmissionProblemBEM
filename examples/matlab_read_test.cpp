//
// Created by diegorenner on 2/20/20.
//
#include <complex>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>

typedef std::complex<double> complex_t;
int main(){
    Eigen::VectorXcd K_expected;
    std::ifstream fp_data;
    double real, imag;
    char sign;
    std::cout << "test" << std::endl;
    fp_data.open("/home/diegorenner/Uni/Thesis/HelmholtzBEM/double_layer_50.dat");
    while(fp_data >> real >> sign >> imag >> sign) {
        std::cout << real << " " << imag << std::endl;
        //K_expected <<  complex_t(real, imag);
    }




}
