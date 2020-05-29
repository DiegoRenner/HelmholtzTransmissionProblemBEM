/**
 * \file parametrized_line.hpp
 * \brief This file declares the class for representing parametrization
 *        of a line segment in 2D.
 *
 *  This File is a part of the 2D-Parametric BEM package
 */

#ifndef PARAMETRIZEDLINEHPP
#define PARAMETRIZEDLINEHPP

#include "abstract_parametrized_curve.hpp"
#include <utility>
#include <Eigen/Dense>

/**
 * \class ParametrizedLine
 * \brief This class represents a Parametrized Line
 *        segment in 2-D and inherits from the
 *        Abstract base class representing parametrized curves
 * @see abstract_parametrized_curve.hpp
 */
class ParametrizedLine : public AbstractParametrizedCurve {
public:
  /**
   * Defining the type for the points on a line segment
   */
  using Point = Eigen::Vector2d;

  /**
   * Constructor with start and end points of the 2D line
   *
   * @param first Starting point for the line segment
   * @param second Ending point for the line segment
   */
  ParametrizedLine(Point first, Point second);

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
   * private const fields for storing the starting and
   * ending point of the line.
   */
  const Point start_;
  const Point end_;
}; // class ParametrizedLine

#endif // PARAMETRIZEDLINEHPP
