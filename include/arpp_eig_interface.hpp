/**
 * \file arpp_eig_interface.hpp
 * \brief This file contains ARPACK to/from EIGEN converters.
 *
 * This File is a part of the HelmholtzTransmissionProblemBEM library
 */
#ifndef ARPP_EIG_INTERFACEHPP
#define ARPP_EIG_INTERFACEHPP

#include "../arpackpp_lib/include/arrscomp.h"

void arpp_to_eig(std::complex<double>* in, Eigen::VectorXcd& out);
void eig_to_arpp(Eigen::VectorXcd& in, std::complex<double>* out);
void arpp_to_eig(ARrcCompStdEig<double>& in, Eigen::VectorXd& out_vals, Eigen::MatrixXcd& out_vectors);

#endif
