/**
 * \file find_roots.hpp
 * \brief This file defines different root finding methods applied
 * to find the minimas in the smallest singular value
 * of the Helmholtz transmission problem solutions operator.
 *
 * This File is a part of the HelmholtzTransmissionProblemBEM library.
 * It incoporates some functions taken from the book
 * "Numerical Recipes in C".
 * These have been marked in their documentation.
 */
#ifndef FIND_ROOTSHPP
#define FIND_ROOTSHPP

#include <Eigen/Dense>
#include <vector>
#include <iomanip>
#include <ostream>
#include "operators.hpp"

using namespace std;
/**
 * Search for roots using Van Wijngaarden-Dekker-Brent method.
 * If there is a root of the function \p f within the interval defined by [\p x1, \p x2]
 * then this function will set the \p root_found flag to true
 * and return the value of the root with a precision of \p tol.
 * This function was taken from the book "Numerical Recipes in C"
 * @param f function to search for root
 * @param x1 initial interval, left bound
 * @param x2 initial interval, right bound
 * @param tol tolerance at which to stop search
 * @param root_found indicator if root was actually foudn
 * @param num_iter numbers of iterations taken to find root
 * @return value of root if root was found
 */
double zbrent( const std::function<double(double)> f,
               double x1,
               double x2,
               double tol,
               bool &root_found,
               unsigned &num_iter);

/**
 * Search for roots using Newton Raphson method.
 * Here the function \p f returns the function of which to find a root and
 * \p f_der returns the derivative of \p f.
 * If there is a root of the function \p f[0] within the interval defined by [\p x1, \p x2]
 * then this function will set the \p root_found flag to true
 * and return the value of the root with a precision of \p tol.
 * This function was taken from the book "Numerical Recipes in C"
 * @param f function to search for root
 * @param f_both function and derivative of function to search for root
 * @param x1 initial interval, left bound
 * @param x2 initial interval, right bound
 * @param tol tolerance at which to stop search
 * @param root_found indicator if root was actually foudn
 * @param num_iter numbers of iterations taken to find root
 * @return value of root if root was found
 */
double rtsafe( const std::function<double(double)> f,
               const std::function<Eigen::MatrixXd(double)> f_both,
               double x1,
               double x2,
               double tol,
               bool &root_found,
               unsigned &num_iter);

/**
 * Search for roots using Newton Raphson method.
 * Here the function \p f returns the function of which to find a root and
 * \p f_der returns the derivative of \p f.
 * If there is a root of the function \p f[0] within the interval defined by [\p a, \p b]
 * then this function will set the \p root_found flag to true
 * and return the value of the root with a precision of \p tol.
 * This function was taken from the book "Numerical Recipes in C"
 * @param f function to search for root
 * @param f_both function and derivative of function to search for root
 * @param a initial interval, left bound
 * @param b initial interval, right bound
 * @param tol tolerance at which to stop search
 * @param root_found indicator if root was actually foudn
 * @param num_iter numbers of iterations taken to find root
 * @return value of root if root was found
 */
std::vector<double> findZeros( const std::function<Eigen::MatrixXd(double)> f_both,
               double a,
               double b,
               double init_len);

/**
 * Search for roots using Newton Raphson method.
 * Here the function \p f returns the function of which to find a root and
 * \p f_der returns the derivative of \p f.
 * If there is a root of the function \p f[0] within the interval defined by [\p a, \p b]
 * then this function will set the \p root_found flag to true
 * and return the value of the root with a precision of \p tol.
 * This function was taken from the book "Numerical Recipes in C"
 * @param f function to search for root
 * @param f_both function and derivative of function to search for root
 * @param a initial interval, left bound
 * @param b initial interval, right bound
 * @param tol tolerance at which to stop search
 * @param root_found indicator if root was actually foudn
 * @param num_iter numbers of iterations taken to find root
 * @return value of root if root was found
 */
template <typename RECORDER = std::function<void(std::vector<grid_data>)>>
std::vector<double> findZeros_seq(const std::function<Eigen::MatrixXd(double)> f_both,
                               double a,
                               double b,
                               unsigned int m,
                               RECORDER rec = [](std::vector<grid_data>)->void{});

/**
 * This function computes the value of the function \p f and it's two
 * first derivatives \p f_der, \p f_der2 at the point \p x0 which are then
 * used to find the next closest minima in the function \p f by parabolic
 * approximation.
 * The function \p f can be vector valued.
 * This will be interpreted as multiple functions being passed and the
 * function that returns the smallest minima is chosen as the reult.
 * @param f function for which to find minima
 * @param f_der first derivative of \p f
 * @param f_der2 second derivative of \p f
 * @param x0 initial guess for minima
 * @param step range in which to seach for minimas
 * @return best approximation of minima and value of function and
 * it's derivatives used to find minima
 */
Eigen::VectorXd parabolic_approximation(const std::function<Eigen::VectorXd(double)> f,
                                        const std::function<Eigen::VectorXd(double)> f_der,
                                        const std::function<Eigen::VectorXd(double)> f_der2,
                                        const double x0,
                                        const double step);

/**
 * Search for root in function \p f given the initial interval [\p x0, \p x1]
 * using the secant .
 * The root does not stay bracketed.
 *
 * @param f function of which to find root
 * @param x0 initial interval, left boundary
 * @param x1 initial interval, right boundary
 * @param tol tolerance with which to determine root
 * @param maxIter maximum number of allowed iterations
 * @param root_found flag to signal if root was found
 * @param num_iter number of iterations taken to find root
 * @return value of found root if a root was found
 */
double secant_method( std::function<double(double)> f,
                      double x0,
                      double x1,
                      const double tol,
                      const unsigned maxIter,
                      bool &root_found,
                      unsigned &num_iter);

std::vector<double> general_cubic_formula(double a, double b, double c, double d, double x0, double x1);

std::vector<double> zerosquadpolstab( double b, double c, double x0, double x1);

void min_shrinkage(double mu, std::vector<double> &pot_zeros, double init_len);
#endif //FIND_ROOTSHPP
