/**
 * \file parametrized_semi_circle.hpp
 * \brief This file declares the class for representing parametrization
 *        of a semi circle.
 *
 *  This File is a part of the 2D-Parametric BEM package
 */

#ifndef PARAMETRIZEDSEMICIRCLEHPP
#define PARAMETRIZEDSEMICIRCLEHPP

#include "abstract_parametrized_curve.hpp"

#include <utility>

#include <Eigen/Dense>

namespace parametricbem2d {
/**
 * \class ParametrizedSemiCircle
 * \brief This class represents the parametrization
 *        of a semi circle. It is centered at (0,0) and lies in
 *        positive x half plane. This class inherits from the
 *        Abstract base class representing parametrized curves
 * @see abstract_parametrized_curve.hpp
 */
class ParametrizedSemiCircle : public AbstractParametrizedCurve {
public:
  /**
   * Constructor with specified radius; default value = 1.
   *
   * @param r radius of the semi circle
   */
  ParametrizedSemiCircle(double r = 1.);

  /**
   * See documentation in AbstractParametrizedCurve
   */
  Eigen::Vector2d operator()(double) const;
    Eigen::Vector2d operator[](double) const;
    Eigen::Vector2d swapped_op(double) const;

  /**
   * See documentation in AbstractParametrizedCurve
   */
  Eigen::Vector2d Derivative(double) const;
    Eigen::Vector2d Derivative_01(double) const;
    Eigen::Vector2d Derivative_01_swapped(double) const;

  /**
   * See documentation in AbstractParametrizedCurve
   */
  Eigen::Vector2d DoubleDerivative(double) const;

  /**
   * See documentation in AbstractParametrizedCurve
   */
  PanelVector split(unsigned int) const;

private:
  /**
   * Private const field for storing the radius
   */
  const double radius_;
}; // class ParametrizedSemiCircle
} // namespace parametricbem2d

#endif // PARAMETRIZEDSEMICIRCLEHPP
