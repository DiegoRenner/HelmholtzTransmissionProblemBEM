/**
 * \file parametrized_fourier_sum.hpp
 * \brief This file declares the class for representing Fourier Sum
 *        based parametrization.
 *
 *  This File is a part of the 2D-Parametric BEM package
 */

#ifndef PARAMETRIZEDFOURIERSUMHPP
#define PARAMETRIZEDFOURIERSUMHPP

#include "abstract_parametrized_curve.hpp"

/**
 * \class ParametrizedFourierSum
 * \brief This class represents Fourier Sum based parametrization
 *        of the form \f$\gamma\f$(t) = \f$c+\sum_{n=1}^{N} a_n cos(nt)+b_n
 *        sin(nt)\f$ and inherits from the Abstract base class representing
 *        parametrized curves
 * @see abstract_parametrized_curve.hpp
 */
class ParametrizedFourierSum : public AbstractParametrizedCurve {
public:
  /**
   * Defining the type for coefficient list of sine and cosine terms for 2D
   * Fourier Sum based parametrizations.
   */
  using CoefficientsList = typename Eigen::Matrix<double, 2, Eigen::Dynamic>;

  /**
   * Constructor with two coefficient list type parameters
   * that is Eigen::MatrixXd with size: 2 X N.
   * Here N is the number of sine and cosine terms in the sum
   *
   * @param center Constant term 'c' in the parametrization
   * @param cos_list Coefficient list for cosine terms
   * @param sin_list Coefficient list for sine terms
   * @param tmin Lower end of the actual parameter interval used which is
   *             linearly mapped to the standard interval
   * @param tmax Upper end of the actual parameter interval used which is
   *             linearly mapped to the standard interval
   */
  ParametrizedFourierSum(Eigen::Vector2d center, CoefficientsList cos_list,
                         CoefficientsList sin_list, double tmin = -1.,
                         double tmax = 1.);

  /**
   * See documentation in AbstractParametrizedCurve
   */
  Eigen::Vector2d operator()(double) const;
  Eigen::Vector2d operator[](double) const;
  Eigen::Vector2d swapped_op(double) const;
  Eigen::ArrayXXcd operator()(const Eigen::ArrayXXd &) const;
  Eigen::ArrayXXcd operator[](const Eigen::ArrayXXd &) const;
  Eigen::ArrayXXcd swapped_op(const Eigen::ArrayXXd &) const;

  /**
   * See documentation in AbstractParametrizedCurve
   */
  Eigen::Vector2d Derivative(double) const;
  Eigen::Vector2d Derivative_01(double) const;
  Eigen::Vector2d Derivative_01_swapped(double) const;
  void Derivative(const Eigen::ArrayXXd &, Eigen::ArrayXXcd &, Eigen::ArrayXXd &) const;
  void Derivative_01(const Eigen::ArrayXXd &, Eigen::ArrayXXcd &, Eigen::ArrayXXd &) const;
  void Derivative_01_swapped(const Eigen::ArrayXXd &, Eigen::ArrayXXcd &, Eigen::ArrayXXd &, bool) const;

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
   * List of coefficients for the cosine terms in Fourier Sum based
   * parametrization
   */
  const CoefficientsList cosine_;

  /**
   * List of coefficients for the sine terms in Fourier Sum based
   * parametrization
   */
  const CoefficientsList sine_;

  /**
   * Storing the actual range of parameter within the parameter range
   * By default, it is exactly equal to the parameter range. It is used
   * for the split functionality to make part Fourier Sum parameterizations.
   */
  const double tmin_, tmax_;

  /**
   * The constant term in the parametrization.
   */
  Eigen::Vector2d center_;
}; // class ParametrizedFourierSum

#endif // PARAMETRIZEDFOURIERSUMHPP
