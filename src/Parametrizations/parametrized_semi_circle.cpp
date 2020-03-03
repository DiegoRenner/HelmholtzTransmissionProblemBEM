/**
 * \file parametrized_semi_circle.cpp
 * \brief This file defines the class representing parametrization
 *        of a semi circle.
 * @see parametrized_semi_circle.hpp
 *
 *  This File is a part of the 2D-Parametric BEM package
 */

#include "parametrized_semi_circle.hpp"

#include <math.h>
#include <assert.h>
#include <iostream>
#include <utility>
#include <memory>

#include <Eigen/Dense>
#include "parametrized_circular_arc.hpp"

#define _USE_MATH_DEFINES //for Pi

namespace parametricbem2d {
  ParametrizedSemiCircle::ParametrizedSemiCircle(double r) : radius_(r) {}

  Eigen::Vector2d ParametrizedSemiCircle::operator() (double t) const {
    assert(IsWithinParameterRange(t));
    // Parametrization using polar coordinates based on parameter t
    Eigen::Vector2d point(radius_*cos(M_PI*t/2.),radius_*sin(M_PI*t/2.));
    return point;
  }

  Eigen::Vector2d ParametrizedSemiCircle::Derivative(double t) const {
    assert(IsWithinParameterRange(t));
    // Derivative of the polar coordinaties used in the function operator()
    Eigen::Vector2d derivative(-radius_*M_PI*sin(M_PI*t/2.)/2.,
                                M_PI*radius_*cos(M_PI*t/2.)/2.);
    return derivative;
  }

  Eigen::Vector2d ParametrizedSemiCircle::DoubleDerivative(double t) const {
    assert(IsWithinParameterRange(t));
    // Double Derivative of the polar coordinaties in the function operator()
    return -M_PI*M_PI/4.*operator() (t);
  }

  PanelVector ParametrizedSemiCircle::split(unsigned int N) const {
    // PanelVector for storing the part parametrizations
    PanelVector parametrization_parts;
    for (int i = 0 ; i < N ; ++i) {
      double phi_min = 0.;
      double phi_max = M_PI/2.;
      Eigen::Vector2d center(2); center << 0,0;
      double phi_start = phi_min + i*(phi_max-phi_min)/N;
      double phi_end = phi_min + (i+1)*(phi_max-phi_min)/N;
      if (i==N-1)
        phi_end = phi_max;
      parametrization_parts.push_back(std::make_shared<ParametrizedCircularArc>
                                      (center,radius_,phi_start,phi_end)) ;
    }
    return parametrization_parts;
  }
} // namespace parametricbem2d
