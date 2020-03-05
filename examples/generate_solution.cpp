//
// Created by diegorenner on 2/14/20.
//

#include <Eigen/Dense>
#include <iostream>
#include "generate_solution.hpp"
#include "math.h"

namespace sol {
    typedef std::complex<double> complex_t;
    complex_t ii = complex_t(0., 1.);

    complex_t hn(int n, double x) {
    return jn(n,x) + ii*yn(n,x);
    }

    complex_t fund_sol(double k, double x1, double x2, double ipt1, double ipt2) {
        double r = sqrt((x1-ipt1)*(x1-ipt1) + (x2-ipt2)*(x2-ipt2));
        return ii/4.*hn(0,k*r);
    }

    complex_t fund_sol_N(double k, double x1, double x2, double ipt1, double ipt2) {
        Eigen::Vector2d x(x1,x2);
        Eigen::Vector2d normal = x.normalized();
        Eigen::Vector2d ipt(ipt1,ipt2);
        double r = (x-ipt).norm();

        return k*ii/4.*hn_der(0,k*r)*(x-ipt).normalized().dot(normal);
    }

    complex_t hn_der(int n, double x) {
        return (hn(n-1,x) - hn(n+1,x))/2.;
    }

    complex_t jn_der(int n, double x) {
        return (jn(n-1,x) - jn(n+1,x))/2.;
    }

    complex_t t_coeff(int n,
                      double eps,
                      double k,
                      double n_i) {
        double k_eps = k * eps;
        double lambda = sqrt(n_i);
        complex_t result = (2. * ii /( M_PI * k_eps)) / (hn_der(n, k_eps) * jn(n, lambda * k_eps) - lambda * jn_der(n, lambda * k_eps) * hn(n, k_eps));
        return result;
    };

    complex_t r_coeff(int n,
                      double eps,
                      double k,
                      double n_i) {
        double k_eps = k * eps;
        double lambda = sqrt(n_i);
        return -(hn_der(n, k_eps) * jn(n, lambda * k_eps) - lambda * jn_der(n, lambda * k_eps) * hn(n, k_eps)).real() /
               (hn_der(n, k_eps) * jn(n, lambda * k_eps) - lambda * jn_der(n, lambda * k_eps) * hn(n, k_eps));
    }

    complex_t u_i(double x1,
                  double x2,
                  int l,
                  double eps,
                  double *a_n,
                  double k,
                  double n_i) {
        complex_t result = complex_t(0.0, 0.0);
        double r = sqrt(x1 * x1 + x2 * x2);
        double eta = atan2(x2 / r, x1 / r);
        if (eta < 0 ) eta += M_PI;
        for (int i = 0; i < 2 * l + 1; i++) {
            result += a_n[i] * jn(i - l, k * r) * complex_t(cos((i - l) * eta), sin((i - l) * eta));
        }
        return result;
    }

    complex_t u_s(double x1,
                  double x2,
                  int l,
                  double eps,
                  double *a_n,
                  double k,
                  double n_i) {
        complex_t result = complex_t(0.0, 0.0);
        double r = sqrt(x1 * x1 + x2 * x2);
        double eta = atan2(x2 / r, x1 / r);
        if (eta < 0 ) eta += M_PI;
        for (int i = 0; i < 2 * l + 1; i++) {
            result += a_n[i] * r_coeff(i - l, eps, k, n_i) * hn(i - l, k * r) *
                      complex_t(cos((i - l) * eta), sin((i - l) * eta));
        }
        return result;
    }

    complex_t u_t(double x1,
                  double x2,
                  int l,
                  double eps,
                  double *a_n,
                  double k,
                  double n_i) {
        complex_t result = complex_t(0.0, 0.0);
        double lambda = sqrt(n_i);
        double r = sqrt(x1 * x1 + x2 * x2);
        double eta = atan2(x2 / r, x1 / r);
        if (eta < 0 ) eta += M_PI;
        for (int i = 0; i < 2 * l + 1; i++) {
            result += a_n[i] * t_coeff(i - l, eps, k, n_i) * jn(i - l, lambda * k * r) *
                      complex_t(cos((i - l) * eta), sin((i - l) * eta));
        }
        return result;
    }

    complex_t u_i_N(double x1,
                    double x2,
                    int l,
                    double eps,
                    double *a_n,
                    double k,
                    double n_i) {
        Eigen::Vector2cd result;
        double r = sqrt(x1 * x1 + x2 * x2);
        double eta = atan2(x2 / r, x1 / r);
        if (eta < 0 ) eta += M_PI;
        Eigen::Vector2d e_r;
        e_r << cos(eta), sin(eta);
        Eigen::Vector2d e_eta;
        e_eta << -sin(eta), cos(eta);
        for (int i = 0; i < 2 * l + 1; i++) {
            result += a_n[i] * complex_t(cos((i - l) * eta), sin((i - l) * eta)) *
                      (k * jn_der(i - l, k * r) * e_r + ii * double(i - l) * jn(i - l, k * r) / r * e_eta);
        }
        return result.dot(e_r);
    }

    complex_t u_s_N(double x1,
                    double x2,
                    int l,
                    double eps,
                    double *a_n,
                    double k,
                    double n_i) {
        Eigen::Vector2cd result;
        double r = sqrt(x1 * x1 + x2 * x2);
        double eta = atan2(x2 / r, x1 / r);
        if (eta < 0 ) eta += M_PI;
        Eigen::Vector2d e_r;
        e_r << cos(eta), sin(eta);
        Eigen::Vector2d e_eta;
        e_eta << -sin(eta), cos(eta);
        for (int i = 0; i < 2 * l + 1; i++) {
            result += a_n[i] * r_coeff(i - l, eps, k, n_i) * complex_t(cos((i - l) * eta), sin((i - l) * eta)) *
                      (k * hn_der(i - l, k * r) * e_r + ii * double(i - l) * hn(i - l, k * r) * 1. / r * e_eta);
        }
        return result.dot(e_r);
    }

    complex_t u_t_N(double x1,
                    double x2,
                    int l,
                    double eps,
                    double *a_n,
                    double k,
                    double n_i) {
        Eigen::Vector2cd result;
        double lambda = sqrt(n_i);
        double r = sqrt(x1 * x1 + x2 * x2);
        double eta = atan2(x2 / r, x1 / r);
        if (eta < 0 ) eta += M_PI;
        Eigen::Vector2d e_r;
        e_r << cos(eta), sin(eta);
        Eigen::Vector2d e_eta;
        e_eta << -sin(eta), cos(eta);
        for (int i = 0; i < 2 * l + 1; i++) {
            result += a_n[i] * (t_coeff(i - l, eps, k, n_i).real() - ii*t_coeff(i - l, eps, k, n_i).imag()) * complex_t(cos((i - l) * eta), sin((i - l) * eta)) *
                      (lambda * k * jn_der(i - l, lambda * k * r) * e_r) ;//+ ii * double(i - l) * jn(i - l, lambda * k * r) * 1. / r * e_eta);
            std::cout << eta << std::endl;
        }
        return result.dot(e_r);
    }
}
