/**
 * \file adj_double_layer.cpp
 * \brief This file defines the functions to evaluate the entries of
 *        Galerkin matrices based on the bilinear form induced by the
 *        Adjoint Double Layer BIO, internally using the Double Layer
 *        implementation
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#include "adj_double_layer.hpp"

#include <math.h>
#include <vector>
#include <limits>

#include <Eigen/Dense>
#include "abstract_bem_space.hpp"
#include "abstract_parametrized_curve.hpp"
#include "discontinuous_space.hpp"
#include "double_layer.hpp"
#include "gauleg.hpp"
#include "logweight_quadrature.hpp"
#include "parametrized_mesh.hpp"

namespace parametricbem2d {
    namespace adj_double_layer_helmholtz {
        Eigen::MatrixXcd GalerkinMatrix(const ParametrizedMesh mesh,
                                       const AbstractBEMSpace &trial_space,
                                       const AbstractBEMSpace &test_space,
                                       const unsigned int &N,
                                       const double k)
                                       {
            // Getting the adjoint double layer matrix by calculating the double layer
            // Galerkin matrix and transposing it
            return parametricbem2d::double_layer_helmholtz::GalerkinMatrix(mesh, test_space,
                                                                 trial_space, N, k);

        }
    } // namespace parametricbem2d
}
