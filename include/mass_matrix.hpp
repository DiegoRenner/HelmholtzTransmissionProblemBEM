/**
 * \file mass_matrix.hpp
 * \brief This file declares the functions to evaluate the entries of
 *        Galerkin matrices based on the bilinear form induced by the
 *        Identity BIO,
 *        using common and composite Gauss-Legendre quadrature rules.
 *
 * This File is a part of the HelmholtzTransmissionProblemBEM library.
 */

#ifndef MASSMATRIXRHPP
#define MASSMATRIXRHPP

#include "abstract_bem_space.hpp"
#include "abstract_parametrized_curve.hpp"
#include "gauleg.hpp"

/**
 * \namespace mass_matrix
 * \brief This namespace contains all the functions for evaluating the Identity
 * Galerkin Matrix using common and composite Gauss-Legendre quadrature rules.
 */
namespace mass_matrix{

    Eigen::MatrixXd GalerkinMatrix(const ParametrizedMesh &mesh,
                                   const AbstractBEMSpace &trial_space,
                                   const AbstractBEMSpace &test_space,
                                   const QuadRule &GaussQR);
} //namespace mass_matrix

#endif // SINGLELAYERHPP
