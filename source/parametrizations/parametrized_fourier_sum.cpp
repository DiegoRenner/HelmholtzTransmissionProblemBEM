/**
 * \file parametrized_fourier_sum.cpp
 * \brief This file defines the class representing Fourier
 *        Sum based parametrization
 * @see parametrized_fourier_sum.hpp
 *
 *  This File is a part of the 2D-Parametric BEM package
 */

#include "parametrized_fourier_sum.hpp"

#include <assert.h>
#include <iostream>
#include <math.h>
#include <utility>
#include <vector>

#include <Eigen/Dense>

using CoefficientsList = typename ParametrizedFourierSum::CoefficientsList;

ParametrizedFourierSum::ParametrizedFourierSum(Eigen::Vector2d center,
                                               CoefficientsList cos_list,
                                               CoefficientsList sin_list,
                                               double tmin, double tmax)
: cosine_(cos_list), sine_(sin_list), tmin_(tmin), tmax_(tmax), center_(center) {
  // Checking consistency
  assert(cosine_.cols() == sine_.cols());
}

Eigen::Vector2d ParametrizedFourierSum::operator()(double t) const {
  assert(IsWithinParameterRange(t));
  t = t * (tmax_ - tmin_) / 2 + (tmax_ + tmin_) / 2; // converting to the range [tmin,tmax]
  // Number of cosine/sine terms in the sum
  int N = cosine_.cols();
  // Vectors for storing the cosine and sine values
  Eigen::VectorXd cos_theta(N);
  Eigen::VectorXd sin_theta(N);
  // Filling the cosine and sine vectors
  for (int i = 0; i < N; ++i) {
    // Computing cosine and sine values at the parameter value t
    cos_theta(i) = cos((i + 1) * t);
    sin_theta(i) = sin((i + 1) * t);
  }
  // Matrix multiplication to create the Fourier Sum (in vector form)
  // [2 x N] X [N x 1] = [2 x 1]
  Eigen::Vector2d point = cosine_ * cos_theta + sine_ * sin_theta;
  return point + center_;
}

Eigen::Vector2d ParametrizedFourierSum::Derivative(double t) const {
  assert(IsWithinParameterRange(t));
  t = t * (tmax_ - tmin_) / 2 + (tmax_ + tmin_) / 2; // converting to the range [tmin,tmax]
  // Number of cosine/sine terms in the sum
  int N = cosine_.cols();
  // Vectors for storing the cosine and sine derivatives
  Eigen::VectorXd cos_theta_dot(N);
  Eigen::VectorXd sin_theta_dot(N);
  // Filling the vectors with derivatives
  for (int i = 0; i < N; ++i) {
    // Computing the derivatives for cosine and sine terms
    cos_theta_dot(i) = -(i + 1) * (tmax_ - tmin_) / 2 * sin((i + 1) * t);
    sin_theta_dot(i) = (i + 1) * (tmax_ - tmin_) / 2 * cos((i + 1) * t);
  }
  // Matrix multiplication to create the derivative of Fourier Sum (in vector
  // form) [2 x N] X [N x 1] = [2 x 1]
  Eigen::Vector2d derivative = cosine_ * cos_theta_dot + sine_ * sin_theta_dot;
  return derivative;
}

Eigen::Vector2d ParametrizedFourierSum::DoubleDerivative(double t) const {
  assert(IsWithinParameterRange(t));
  t = t * (tmax_ - tmin_) / 2 + (tmax_ + tmin_) / 2; // converting to the range [tmin,tmax]
  // Number of cosine/sine terms in the sum
  int N = cosine_.cols();
  // Vectors for storing the cosine and sine double derivatives
  Eigen::VectorXd cos_theta_ddot(N);
  Eigen::VectorXd sin_theta_ddot(N);
  // Filling the vectors with double derivatives
  for (int i = 0; i < N; ++i) {
    // Computing the double derivatives for cosine and sine terms
    cos_theta_ddot(i) = -(i + 1) * (i + 1) * (tmax_ - tmin_) / 2 *
    (tmax_ - tmin_) / 2 * cos((i + 1) * t);
    sin_theta_ddot(i) = -(i + 1) * (i + 1) * (tmax_ - tmin_) / 2 *
    (tmax_ - tmin_) / 2 * sin((i + 1) * t);
  }
  // Matrix multiplication to create the double derivative of Fourier Sum
  // [2 x N] X [N x 1] = [2 x 1]
  Eigen::Vector2d double_derivative =
  cosine_ * cos_theta_ddot + sine_ * sin_theta_ddot;
  return double_derivative;
}

PanelVector ParametrizedFourierSum::split(unsigned int N) const {
  // PanelVector for storing the part parametrizations
  PanelVector parametrization_parts;
  // Generating the parts
  for (int i = 0; i < N; ++i) {
    // Partitioning by splitting the parameter range [tmin,tmax]
    double tmin = tmin_ + i * (tmax_ - tmin_) / N;
    double tmax = tmin_ + (i + 1) * (tmax_ - tmin_) / N;
    if (i == N - 1)
      tmax = tmax_;
    // Adding the part parametrization to the vector with a shared pointer
    parametrization_parts.push_back(std::make_shared<ParametrizedFourierSum>(center_, cosine_, sine_, tmin, tmax));
  }
  return parametrization_parts;
}
