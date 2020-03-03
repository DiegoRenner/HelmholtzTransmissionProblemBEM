/**
 * \file abstract_parametrized_curve.hpp
 * \brief This file declares an abstract class for representing
 *        Parametrized curves.
 *
 *  This File is a part of the 2D-Parametric BEM package
 */

#ifndef PARAMETRIZEDCURVEHPP
#define PARAMETRIZEDCURVEHPP

#include <iostream>
#include <iterator>
#include <memory>
#include <utility>
#include <vector>

#include <Eigen/Dense>
/**
 * \namespace parametricbem2d
 * \brief This namespace contains all the classes and functions for
 *        the 2D-Parametric BEM package.
 */
namespace parametricbem2d {
// Forward declaration of the class for using it in the typedef for PanelVector
class AbstractParametrizedCurve;

/**
 * This typedef is used for defining a vector of AbstractParametrizedCurve
 * pointers which can store any sequence of parametrized curves derived from the
 * abstract base class AbstractParametrizedCurve. The PanelVector is also used
 * to store the components of a mesh by using additional constraints. Smart
 * pointers are used to free up memory whenever the pointer is destroyed.
 */
using PanelVector = std::vector<std::shared_ptr<AbstractParametrizedCurve>>;

/**
 * \class AbstractParametrizedCurve
 * \brief This abstract class provides the interface for a parametric
 *        curve. The interface follows the one outlined in
 *        \f$\ref{cpp:crv}\f$ in the lecture material for
 *        Advanced Numerical Methods for CSE.
 */
class AbstractParametrizedCurve {
public:
  /**
   * This function is used for querying the parameter interval.
   * The standard parameter interval [-1,1] is used and it can't
   * be overriden in the inherited classes as the function is non
   * virtual. The function is declared static as it is independent
   * of the concrete object.
   *
   * @return A std::pair<double,double> object containing
   *          the valid parameter range for parametrization
   *          that is [-1,1]
   */
  static std::pair<double, double> ParameterRange(void) {
    // Parameter range : [-1,1]
    return std::make_pair(-1., 1.);
  }

  /**
   * This function is used for accessing a point on the parametrized
   * curve with parameter "t" for the parametrization \f$\gamma\f$(t).
   * This is a pure virtual function which has to be implemented
   * in the inherited classes.
   *
   * @param t Double type argument for the parameter value.
   * @return A 2-D vector of type Eigen::Vector2d
   *         containing the coordinates of the
   *         parametrized point.
   */
  virtual Eigen::Vector2d operator()(double t) const = 0;
  virtual Eigen::Vector2d operator[](double t) const = 0;
  virtual Eigen::Vector2d swapped_op(double t) const = 0;

  /**
   * This function is used for retrieving the derivative \f$\dot{\gamma}\f$(t)
   * at a point on the parametrized curve \f$\gamma\f$(t).
   * This is a pure virtual function which has to be implemented
   * in the inherited classes.
   *
   * @param t Parameter value of Double type
   *          at which the derivative for
   *          parametrization is evaluated
   * @return A 2-D vector of type Eigen::Vector2d
   *         containing the parametrized
   *         derivative at point 't'
   */
  virtual Eigen::Vector2d Derivative(double t) const = 0;
  virtual Eigen::Vector2d Derivative_01(double t) const = 0;
  virtual Eigen::Vector2d Derivative_01_swapped(double t) const = 0;

  /**
   * This function is used for retrieving the double derivative
   * \f$\ddot{\gamma}\f$(t)
   * at a point on the parametrized curve \f$\gamma\f$(t).
   * This is a pure virtual function which has to be implemented
   * in the inherited classes.
   *
   * @param t Parameter value of Double type
   *          at which the double derivative for
   *          the parametrization is evaluated
   * @return A 2-D vector of type Eigen::Vector2d
   *         containing the parametrized double
   *         derivative at point 't'
   */
  virtual Eigen::Vector2d DoubleDerivative(double t) const = 0;

  /**
   * This function is used for checking whether a value t is within the
   * valid parameter range. This function is non virtual to prevent it
   * from being overriden as the parameter interval is fixed. It is
   * declared static because the check is independent of the concrete
   * implementation.
   *
   * @param t The value to be checked
   * @return boolean indicating result of the performed check
   */
  static bool IsWithinParameterRange(double t) {
    double a, b;
    // Getting the parameter range
    std::tie(a, b) = ParameterRange();
    // Checking if the value is within parameter range
    return (t >= a && t <= b);
  }

  /**
   * This function is used for splitting a parametrized curve into several
   * self-similar part curves. It is useful to make a parametrized mesh as it
   * returns all the part-parametrizations stored in a sequence in a PanelVector
   * such that the end point of a part is the starting point of the next part.
   *
   * @param N Integer indicating the number of parts for split
   * @return PanelVector containing all the part parametrizations
   */
  virtual PanelVector split(unsigned int N) const = 0;

  /**
   * This function is used for calculating distance to another parametrized
   * curve which is passed as an input to the function. This distance
   * computation is done approximately by evaluating distances on a sample and
   * finding the minimum out of them. A sample of size 5 X 5 is hard coded into
   * the function. Please note that this is an approximate computation which is
   * intended for evaluating the admissibility in the case of not intersecting
   * panels
   *
   * @param curve Parametrized curve, distance to which is evaluated
   * @return double Distance between the parametric curves
   */
  double distanceTo(const AbstractParametrizedCurve &curve) const {
    // Number of points sampled in one dimension for naive distance calculation
    unsigned N = 5;
    double tmin, tmax, min_dist;
    // Getting the parameter range
    std::tie(tmin, tmax) = ParameterRange();
    // Sample points over the domain
    Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(N, tmin, tmax);
    // Evaluating distances over all possible pairs
    for (unsigned i = 0; i < N; ++i) {
      for (unsigned j = 0; j < N; ++j) {
        double temp = (this->operator()(t(i)) - curve.operator()(t(j))).norm();
        min_dist =
            i == 0 && j == 0 ? temp : (temp < min_dist ? temp : min_dist);
      }
    }
    return min_dist;
  }

  /**
   * This function is used for calculating the length of a panel. The function
   * in this abstract base class naively calculates distance by dividing the
   * curve into line segments. The function is virtual and can be overriden in
   * derived classes for a more accurate implementation.
   *
   * @return double Approximate length of the parametrized curve
   */
  virtual double length() const {
    // Number of points sampled for naive length calculation
    unsigned N = 10;
    double tmin, tmax, length = 0.;
    // Getting the parameter range
    std::tie(tmin, tmax) = ParameterRange();
    // Sample points over the domain
    Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(N, tmin, tmax);
    for (unsigned i = 0; i < N - 1; ++i) {
      length += (this->operator()(t(i)) - this->operator()(t(i + 1))).norm();
    }
    return length;
  }
}; // class AbstractParametrizedCurve

/**
 * This function is used for calculating the admissibility of two disjoint
 * parametrized curves. It is defined as \f$\rho(\Pi,\Pi') = \frac
 * {max(|\Pi|,|\Pi'|)} {dist(\Pi,\Pi')}\f$
 *
 * @param pi First parametrized curve
 * @param pi_p Second parametrized curve
 * @return Admissibility value for the given pair of panels
 */
inline double rho(const AbstractParametrizedCurve &pi,
                  const AbstractParametrizedCurve &pi_p) {
  // Numerator: longer panel length
  double num = std::max(pi.length(), pi_p.length());
  // Denominator: distance between panels
  double denom = pi.distanceTo(pi_p);
  // Admissibility formula
  return num / denom;
}

} // namespace parametricbem2d

#endif // PARAMETRIZEDCURVEHPP
