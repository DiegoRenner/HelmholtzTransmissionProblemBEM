/**
 * \file parametrized_circular_arc.hpp
 * \brief This file declares the class for representing parametrization
 *        of a circular arc.
 *
 *  This File is a part of the 2D-Parametric BEM package
 */

#ifndef PARAMETRIZEDCIRCULARARCHPP
#define PARAMETRIZEDCIRCULARARCHPP

#include "abstract_parametrized_curve.hpp"

#include <utility>

#include <Eigen/Dense>

namespace parametricbem2d {
/**
 * \class ParametrizedCircularArc
 * \brief This class represents the parametrization
 *        of a circular arc. This class inherits from the
 *        Abstract base class representing parametrized curves
 * @see abstract_parametrized_curve.hpp
 */
class ParametrizedCircularArc : public AbstractParametrizedCurve {
public:
  /**
   * Constructor with specified center,radius and the starting and ending
   * polar angles
   *
   * @param center A 2-D vector containing the center for the arc
   * @param r radius of the circular arc
   * @param phi_start starting polar angle for the arc
   * @param phi_end ending polar angle for the arc
   */
  ParametrizedCircularArc(Eigen::Vector2d center, double r, double phi_start,
                          double phi_end);

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

  /**
   * See documentation in AbstractParametrizedCurve
   */
  double length() const;

private:
  /**
   * Private const field for storing the center
   */
  const Eigen::Vector2d center_;
  /**
   * Private const field for storing the radius
   */
  const double radius_;
  /**
   * Private const field for storing the phi_start
   */
  const double phi_start_;
  /**
   * Private const field for storing the phi_end
   */
  const double phi_end_;
}; // class ParametrizedCircularArc
} // namespace parametricbem2d

#endif // PARAMETRIZEDCIRCULARARCHPP
