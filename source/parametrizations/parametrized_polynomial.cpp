/**
 * \file parametrized_polynomial.cpp
 * \brief This file defines the class for representing a polynomial
 *        parametrization
 * @see parametrized_polynomial.cpp
 *
 *  This File is a part of the 2D-Parametric BEM package
 */

#include "parametrized_polynomial.hpp"

#include <assert.h>
#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

#include <Eigen/Dense>

using CoefficientsList = typename ParametrizedPolynomial::CoefficientsList;

ParametrizedPolynomial::ParametrizedPolynomial(CoefficientsList coeffs,
                                               double tmin, double tmax)
    : coeffs_(coeffs), tmin_(tmin), tmax_(tmax) {
  assert(coeffs.cols() >= 3); // Asserting degree is at least 2
}

Eigen::Vector2d ParametrizedPolynomial::operator()(double t) const {
  assert(IsWithinParameterRange(t));
  t = t * (tmax_ - tmin_) / 2 +
      (tmax_ + tmin_) / 2; // converting to the range [tmin,tmax]
  // Number of terms in the polynomial
  int N = coeffs_.cols();
  // The output vector
  Eigen::VectorXd point = Eigen::VectorXd::Zero(2);
  // Filling the cosine and sine vectors
  for (int i = 0; i < N; ++i) {
    point += std::pow(t, i) * coeffs_.col(i);
  }
  return point;
}

Eigen::Vector2d ParametrizedPolynomial::Derivative(double t) const {
  assert(IsWithinParameterRange(t));
  t = t * (tmax_ - tmin_) / 2 +
      (tmax_ + tmin_) / 2; // converting to the range [tmin,tmax]
  // Number of terms in the polynomial
  int N = coeffs_.cols();
  // The output vector
  Eigen::VectorXd derivative = Eigen::VectorXd::Zero(2);
  // Filling the cosine and sine vectors
  for (int i = 1; i < N; ++i) {
    derivative += i * std::pow(t, i - 1) * coeffs_.col(i);
  }
  return derivative * (tmax_ - tmin_) / 2;
}

Eigen::Vector2d ParametrizedPolynomial::DoubleDerivative(double t) const {
  assert(IsWithinParameterRange(t));
  t = t * (tmax_ - tmin_) / 2 +
      (tmax_ + tmin_) / 2; // converting to the range [tmin,tmax]
  // Number of terms in the polynomial
  int N = coeffs_.cols();
  // The output vector
  Eigen::VectorXd double_derivative = Eigen::VectorXd::Zero(2);
  // Filling the cosine and sine vectors
  for (int i = 2; i < N; ++i) {
    double_derivative += i * (i - 1) * std::pow(t, i - 2) * coeffs_.col(i);
  }
  return double_derivative * (tmax_ - tmin_) / 2 * (tmax_ - tmin_) / 2;
}

PanelVector ParametrizedPolynomial::split(unsigned int N) const {
  // PanelVector for storing the part parametrizations
  PanelVector parametrization_parts;
  // Generating the parts
  for (int i = 0; i < N; ++i) {
    // Partitioning by splitting the parameter range [tmin,tmax]
    double tmin = tmin_ + i * (tmax_ - tmin_) / N;
    double tmax = tmin_ + (i + 1) * (tmax_ - tmin_) / N;
    if (i==N-1)
      tmax = tmax_;
    // Adding the part parametrization to the vector with a shared pointer
    parametrization_parts.push_back(std::make_shared<ParametrizedPolynomial>(coeffs_, tmin, tmax));
  }
  return parametrization_parts;
}
