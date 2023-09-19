/**
 * \file parametrized_circular_arc.cpp
 * \brief This file defines the class for representing parametrization
 *        of a circular arc.
 *
 *  This File is a part of the 2D-Parametric BEM package
 */

#include "parametrized_circular_arc.hpp"
#include <math.h>
#include <utility>
#include <Eigen/Dense>

#define _USE_MATH_DEFINES // for Pi

using namespace std::complex_literals;

ParametrizedCircularArc::ParametrizedCircularArc()
: center_(Eigen::Vector2d(0.0, 0.0)), radius_(0.25), phi_start_(0.0), phi_end_(2 * M_PI) { }

ParametrizedCircularArc::ParametrizedCircularArc(Eigen::Vector2d center,
                                                 double r,
                                                 double phi_start,
                                                 double phi_end)
: center_(center), radius_(r), phi_start_(phi_start), phi_end_(phi_end) {
    // Asserting that the circular arc does not overlap itself
    assert(fabs(phi_end - phi_start) <= 2 * M_PI);
}

Eigen::Vector2d ParametrizedCircularArc::operator()(double t) const {
    assert(IsWithinParameterRange(t));
    // Parametrization using polar coordinates based on parameter t
    double mean = (phi_start_ + phi_end_) / 2.;
    double difference = (phi_end_ - phi_start_) / 2.;
    Eigen::Vector2d point(center_(0) + radius_ * cos(t * difference + mean),
                          center_(1) + radius_ * sin(t * difference + mean));
    return point;
}
Eigen::ArrayXXcd ParametrizedCircularArc::operator()(const Eigen::ArrayXXd &t) const {
    assert(IsWithinParameterRange(t));
    // Parametrization using polar coordinates based on parameter t
    double mean = (phi_start_ + phi_end_) / 2.;
    double difference = (phi_end_ - phi_start_) / 2.;
    return center_(0) + 1i * center_(1) + radius_ * (1i * (t * difference + mean)).exp();
}
Eigen::Vector2d ParametrizedCircularArc::operator[](double t) const {
    assert(IsWithinParameterRange(t));
    // Parametrization using polar coordinates based on parameter t
    double difference = (phi_end_ - phi_start_);
    Eigen::Vector2d point(center_(0) + radius_ * cos(t * difference + phi_start_),
                          center_(1) + radius_ * sin(t * difference + phi_start_));
    return point;
}
Eigen::ArrayXXcd ParametrizedCircularArc::operator[](const Eigen::ArrayXXd &t) const {
    assert(IsWithinParameterRange(t));
    // Parametrization using polar coordinates based on parameter t
    double difference = (phi_end_ - phi_start_);
    return center_(0) + 1i * center_(1) + radius_ * (1i * (t * difference + phi_start_)).exp();
}
Eigen::Vector2d ParametrizedCircularArc::swapped_op(double t) const {
    assert(IsWithinParameterRange(t));
    // Parametrization using polar coordinates based on parameter t
    double difference = (phi_start_ - phi_end_);
    Eigen::Vector2d point(center_(0) + radius_ * cos(t * difference + phi_end_),
                          center_(1) + radius_ * sin(t * difference + phi_end_));
    return point;
}
Eigen::ArrayXXcd ParametrizedCircularArc::swapped_op(const Eigen::ArrayXXd &t) const {
    assert(IsWithinParameterRange(t));
    // Parametrization using polar coordinates based on parameter t
    double difference = (phi_start_ - phi_end_);
    return center_(0) + 1i * center_(1) + radius_ * (1i * (t * difference + phi_end_)).exp();
}

Eigen::Vector2d ParametrizedCircularArc::Derivative(double t) const {
    assert(IsWithinParameterRange(t));
    // Derivative of the polar coordinaties used in the function operator()
    double mean = (phi_start_ + phi_end_) / 2.;
    double difference = (phi_end_ - phi_start_) / 2.;
    Eigen::Vector2d derivative(-radius_ * difference * sin(t * difference + mean),
                                radius_ * difference * cos(t * difference + mean));
    return derivative;
}
void ParametrizedCircularArc::Derivative(const Eigen::ArrayXXd &t, Eigen::ArrayXXcd &res, Eigen::ArrayXXd &norm) const {
    assert(IsWithinParameterRange(t));
    // Derivative of the polar coordinaties used in the function operator()
    double mean = (phi_start_ + phi_end_) / 2.;
    double difference = (phi_end_ - phi_start_) / 2.;
    res = 1i * radius_ * difference * (1i * (t * difference + mean)).exp();
    norm.setConstant(radius_ * std::abs(difference));
}
Eigen::Vector2d ParametrizedCircularArc::Derivative_01(double t) const {
    assert(IsWithinParameterRange(t));
    // Derivative of the polar coordinaties used in the function operator()
    double difference = (phi_end_ - phi_start_);
    Eigen::Vector2d derivative(-radius_ * difference* sin(t * difference + phi_start_),
                               radius_ * difference* cos(t * difference + phi_start_));
    return derivative;
}
void ParametrizedCircularArc::Derivative_01(const Eigen::ArrayXXd &t, Eigen::ArrayXXcd &res, Eigen::ArrayXXd &norm) const {
    assert(IsWithinParameterRange(t));
    // Derivative of the polar coordinaties used in the function operator()
    double difference = (phi_end_ - phi_start_);
    res = 1i * radius_ * difference * (1i * (t * difference + phi_start_)).exp();
    norm.setConstant(radius_ * std::abs(difference));
}
Eigen::Vector2d ParametrizedCircularArc::Derivative_01_swapped(double t) const {
    assert(IsWithinParameterRange(t));
    // Derivative of the polar coordinaties used in the function operator()
    double difference = (phi_start_ - phi_end_);
    Eigen::Vector2d derivative(-radius_ * difference* sin(t * difference + phi_end_),
                               radius_ * difference* cos(t * difference + phi_end_));
    return derivative;
}
void ParametrizedCircularArc::Derivative_01_swapped(const Eigen::ArrayXXd &t, Eigen::ArrayXXcd &res, Eigen::ArrayXXd &norm, bool neg) const {
    assert(IsWithinParameterRange(t));
    // Derivative of the polar coordinaties used in the function operator()
    double difference = (phi_start_ - phi_end_);
    res = (neg ? -1. : 1.) * 1i * radius_ * difference * (1i * (t * difference + phi_end_)).exp();
    norm.setConstant(radius_ * std::abs(difference));
}

Eigen::Vector2d ParametrizedCircularArc::DoubleDerivative(double t) const {
    assert(IsWithinParameterRange(t));
    // Double derivative of the polar coordinaties used in the function operator()
    double mean = (phi_start_ + phi_end_) / 2.;
    double diff = (phi_end_ - phi_start_) / 2.;
    Eigen::Vector2d double_derivative(
        -radius_ * diff * diff * cos(t * diff + mean),
                                      -radius_ * diff * diff * sin(t * diff + mean));
    return double_derivative;
}

PanelVector ParametrizedCircularArc::split(unsigned int N) const {
    // PanelVector for storing the part parametrizations
    PanelVector parametrization_parts;
    // Generating the parts
    for (int i = 0; i < N; ++i) {
        // Partitioning by splitting the polar angle
        double phi_start = phi_start_ + i * (phi_end_ - phi_start_) / N;
        double phi_end = phi_start_ + (i + 1) * (phi_end_ - phi_start_) / N;
        if (i==N-1)
            phi_end = phi_end_;
        // Adding the part parametrization to the vector with a shared pointer
        parametrization_parts.push_back(std::make_shared<ParametrizedCircularArc>(center_, radius_, phi_start, phi_end));
    }
    return parametrization_parts;
}

double ParametrizedCircularArc::length() const {
    return radius_ * (phi_end_ - phi_start_);
}

