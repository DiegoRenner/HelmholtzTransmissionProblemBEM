#include <limits>
#include <functional>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "find_roots.hpp"

#define MAXIT 1000
#define EPS std::numeric_limits<double>::epsilon()

//using namespace std;
typedef std::complex<double> complex_t;
double zbrent(const function<double(double)> f,
              double x1,
              double x2,
              double tol,
              bool &root_found,
              unsigned &num_iter){
    // counter for iterations
    int iter;
    // initialize function values and boundaries
    double a=x1, b=x2, c=x2, d, e, min1, min2;
    double fa=f(a), fb=f(b), fc, p, q, r, s, tol1, xm;
    // sanity checks
    if (fa * fb > 0.) {
#ifdef CMDL
        cout << "Root must be bracketed in zbrent" << endl;
#endif
        return 0.0;
    }
    fc=fb;
    for (iter=1; iter <= MAXIT; iter++) {
        // reorient boundary for next interpolation
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c = a;
            fc=fa;
            e=d=b-a;
        }
        if (fabs(fc) < fabs(fb)){
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        }
        // check if converged
        tol1=2.0*EPS*fabs(b)+0.5*tol;
        xm=0.5*(c-b);
        if (fabs(xm) <= tol1 || fb == 0.0){
            root_found = true;
            num_iter = iter;
            return b;
        }
        // try quadratic interpolation if
        // bounds are decreasing
        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
            // compute parameters for quadratic interpolation
            s=fb/fa;
            if (a == c) {
                p=2.0*xm*s;
                q=1.0-s;
            } else {
                q=fa/fc;
                r=fb/fc;
                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                q=(q-1.0)*(r-1.0)*(s-1.0);
            }
            // check if quadratic interpolation
            // would fall within bounds
            if (p > 0.0) q = -q;
            p=fabs(p);
            min1=3.0*xm*q-fabs(tol1*q);
            min2=fabs(e*q);
            if(2.0*p < (min1 < min2 ? min1 : min2)) {
                // interpolation accepted
                e=d;
                d=p/q;
            } else {
                // interpolation rejected, bisect instead
                d=xm;
                e=d;
            }
        } else {
            // convergence too slow, bounds not collapsing
            // fast enough, bisect
            d = xm;
            e = d;
        }
        // store previous best guess before
        // computing new best guess
        a=b;
        fa=fb;
        // compute new best guess
        // if step taken is too small, take
        // minimum accepted step towards 0
        if (fabs(d) > tol1)
            b += d;
        else
            b += (xm >= 0) ? fabs(tol1) : -fabs(tol1);
        // one new function evaluation per iteration
        fb = f(b);
    }
#ifdef CMDL
    cout << "Maximum number of iterations exceeded in zbrent" << endl;
#endif
    // should never be reached
    return 0.0;
}

double rtsafe(std::function<double(double)> fct,
              std::function<Eigen::MatrixXd(double)> fct_both,
              double x1,
              double x2,
              double tol,
              bool &root_found,
              unsigned &num_iter){
    // initialize counter, function values and boundaries
    int j;
    double df,dx,dxold,f,fh,fl;
    double temp,xh,xl,rts;
    fl = fct(x1);
    fh = fct(x2);
    // sanity checks
    if (fl * fh > 0.) {
#ifdef CMDL
        std::cout << "Root must be bracketed in rtsafe" << std::endl;
#endif
        if (fl > fh)
            return x1;
        return x2;
    }
    if (fl == 0.0) {
        root_found = true;
        return x1;
    }
    if (fh == 0.0) {
        root_found = true;
        return x2;
    }
    // assign boundaries s.t f(xl) < 0
    if (fl < 0.0) {
        xl=x1;
        xh=x2;
    } else {
        xh=x1;
        xl=x2;
    }
    // set first guess for root, lest step, and "step before last"
    rts=0.5*(x1+x2);
    dxold=fabs(x2-x1);
    dx=dxold;
    // initialize functions and derivative value
    Eigen::MatrixXd tempM = fct_both(rts);
    f = tempM(0,0);
    df = tempM(0,1);
    // loop over max number of allowed iterations
    for (j=1;j<=MAXIT;j++) {
        // check if newton step is out of range or too slow
        // if it is, bisect
        if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) || (fabs(2.0*f) > fabs(dxold*df))) {
            dxold=dx;
            dx=0.5*(xh-xl);
            rts=xl+dx;
            // if change in root is negligible assume convergence
            if (xl == rts) {
                root_found = true;
                num_iter = j;
                return rts;
            }
        } else {
            // newton step accepted
            dxold=dx;
            dx=f/df;
            temp=rts;
            rts -= dx;
            // if change in root is negligible assume convergence
            if (temp == rts) {
                root_found = true;
                num_iter = j;
                return rts;
            }
        }
        // if change in root is small enough assume convergence
        if (fabs(dx) < tol) {
            root_found = true;
            num_iter = j;
            return rts;
        }
        // one new function/derivative evaluation per iteration
        tempM = fct_both(rts);
        f = tempM(0,0);
        df = tempM(0,1);
        // assign new boundary making sure root stay bracketed
        if (f < 0.0)
            xl=rts;
        else
            xh=rts;
    }
#ifdef CMDL
    std::cout << "Maximum number of iterations exceeded in rtsafe" << std::endl;
#endif
    // should never be reached
    return 0.0;
};



std::vector<double> findZeros( std::function<Eigen::MatrixXd(double)> fct_both,
                               double a,
                               double b,
                               double init_len ){

    std::cout << std::endl;
    std::cout << "findZeros called" << std::endl;
    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
    // initialize storage for all real zeros found
    std::vector<double> zeros;

    // evaluate function and it's derivative at boundaries
    Eigen::MatrixXd tempM = fct_both(a);
    double f_a = tempM(0,0);
    double df_a = tempM(0,1);
    tempM = fct_both(b);
    double f_b = tempM(0,0);
    double df_b = tempM(0,1);
    std::cout << "f(a) = " << f_a << std::endl;
    std::cout << "f(b) = " << f_b << std::endl;
    std::cout << "df(a) = " << df_a << std::endl;
    std::cout << "df(b) = " << df_b << std::endl;

    // set tolerances
    // jump detection
    double sigma = 10.0;
    // interpolation tolerance
    double eta = 0.1;
    // interval tolerances
    double tau_abs = 0.00001;
    double tau_rel = 0.0001;
    // early stopping, unused
    double tau = 1e-6;
    // minimal shrinkage parameter
    double mu = 0.1;
    double mu_splitting = 0.1 / (1 + std::abs(a));

    // early stopping, might lose zeros when checking like this here
    // std::cout << f_a << " " << f_b << std::endl;
    /*if (std::abs(f_a) < tau){
        zeros.push_back(a);
        return zeros;
    }
    if (std::abs(f_b) < tau){
        zeros.push_back(b);
        return zeros;
    }*/

    // check for sign change
    if (f_a*f_b <= 0){

        // jump detection
        double val = std::abs((f_b-f_a)/(b - a));
        double tol = sigma*std::max(std::abs(df_a),std::abs(df_b));
        if (val > tol){
            std::cout << "jump detected" << std::endl;
            return zeros;
        }

        // minima detecting
        double val1 = std::abs(b - a);
        double tol1 = tau_rel*std::min(std::abs(a), std::abs(b));
        double val2  = std::abs(b - a);
        double tol2  = tau_abs;
        if (val1 <= tol1 || val2 <= tol2){
            std::cout << "sign change + small interval + no jump detected" << std::endl;
            zeros.push_back(0.5*(a + b));
            return zeros;
        }
    }

    // hermite polynom using h functions
    auto h_00 = [&] (double t) {
        return (1+2*t)*std::pow(1-t,2);
    };
    auto h_10 = [&] (double t) {
        return t*std::pow(1-t,2);
    };
    auto h_01 = [&] (double t) {
        return std::pow(t,2)*(3-2*t);
    };
    auto h_11 = [&] (double t) {
        return std::pow(t,2)*(t-1);
    };

    double f_l = f_a - df_a;
    double f_r = f_b + df_b;

    auto p = [&] (double x) {
        double t = (x - a) / (b - a);
        return h_00(t)*f_a + h_10(t) * (b - a) * df_a
               + h_01(t)*f_b
               + h_11(t) * (b - a) * df_b;
    };

    // hermite polynom using explicit coefficients
    double d_p = f_a;
    double c_p = (b-a)*df_a;

    Eigen::Matrix2d A;
    A << 1, 1, 3, 2;
    Eigen::Vector2d rhs;
    rhs << f_b-f_a-(b-a)*df_a, (b-a)*df_b - (b-a)*df_a;

    Eigen::Vector2d res = A.householderQr().solve(rhs);
    double a_p = res[0];
    double b_p = res[1];
    std::cout << "Coefficients of Hermite Polynomial: " << a_p << " " << b_p << " " << c_p << " " << d_p << std::endl;

    // compute zeros of hermite interpolating function
    std::vector<double> pot_zeros;
    if ( std::abs(a_p) > 1000*EPS){
        pot_zeros = general_cubic_formula(a_p, b_p, c_p, d_p, a, b);
    } else {
        if (std::abs(b_p) > 1000*EPS){
            double a_p2 = 1.0;
            double b_p2 = c_p/(b_p);
            double c_p2 = d_p/(b_p);
            std::cout << "fall back to quadratic case" << std::endl;
            std::cout << "New coefficients: " << std::endl;
            std::cout << a_p2 << " " << b_p2 << " " << c_p2 << std::endl;
            pot_zeros = zerosquadpolstab(b_p2,c_p2,a,b);

        } else {
            if (std::abs(c_p) > 1000*EPS){
                std::cout << "fall back to linear case" << std::endl;
                if (-d_p/c_p > 0 && -d_p/c_p < 1){
                    std::cout << "Zero found: " << std::endl;
                    std::cout << -d_p/c_p << std::endl;
                    pot_zeros.push_back((1+d_p/c_p)*a - (d_p/c_p)*b);
                }
            }
        }
    }
    std::cout << "Zeros of Hermite Polynomial: " << pot_zeros << std::endl;

    // If there are potential zeros, use them to split the interval and rerun the algorithm.
    if (pot_zeros.end() != pot_zeros.begin()){
        pot_zeros.insert(pot_zeros.begin(), a);
        pot_zeros.push_back(b);
        std::sort(pot_zeros.begin(), pot_zeros.end());

        // another point for early stopping, maybe more reasonable but dangerous as well
        /*int size_prev = pot_zeros.size();
        pot_zeros.erase( unique( pot_zeros.begin(), pot_zeros.end() ), pot_zeros.end() );
        int size = pot_zeros.size();
        if ( size_prev > size ) {
            for (auto it = pot_zeros.begin(); it != (pot_zeros.end());) {
                if (std::abs(fct(*it)) < 10*EPS){
                    zeros.push_back(*it);
                    pot_zeros.erase(it);
                    //std::remove(pot_zeros.begin(), pot_zeros.end(), it);
                    //std::cout << "size" << " " << pot_zeros.size() << std::endl;
                    if( pot_zeros.size() < 2){
                        return zeros;
                    }
                } else {
                    ++it;
                }

            }
        }*/
        min_shrinkage(mu, pot_zeros, init_len);
        for (auto it = pot_zeros.begin(); it != (pot_zeros.end()-1); ++it){
            std::vector<double> temp = findZeros(fct_both, *it, *(it + 1), init_len);
            zeros.insert(zeros.end(),temp.begin(),temp.end());
        }

    } else {
        if ( std::abs(a_p) > 1000*EPS){
            double a_p2 = 1.0;
            double b_p2 = 2*b_p/(3*a_p);
            double c_p2 = c_p/(3*a_p);
            pot_zeros = zerosquadpolstab(b_p2,c_p2,a,b);
        } else {
            if (std::abs(b_p) > 1000 * EPS) {
                double b_p2 = 1.0;
                double c_p2 = c_p / (2 * b_p);
                double pot_zero = -c_p2;
                std::cout << "fall back to linear case" << std::endl;
                if (pot_zero > 0 && pot_zero < 1) {
                    std::cout << "Zero found: " << std::endl;
                    std::cout << pot_zero << std::endl;
                    pot_zeros.push_back((1 - pot_zero) * a + (pot_zero) * b);
                }
            }
        }

        std::cout << "Extrema of Hermite Polynomial: " << pot_zeros << std::endl;

        for (auto it = pot_zeros.begin(); it != pot_zeros.end();){
            double val_01 = (*it-a)/(b-a);
            double val_at_min = a_p*pow(val_01,3) + b_p*pow(val_01,2) + c_p*(val_01) + d_p;
            double der2_at_min = 6*a_p*(val_01) + 2*b_p;
            if(!(der2_at_min*val_at_min >= 0)){
                std::cout << "Removed Value " << *it << " due to being positive Max/ negative Min." << std::endl;
                pot_zeros.erase(it);
            } else if (!(std::abs(val_at_min)<eta*std::min(std::abs(f_a),std::abs(f_b)))) {
                std::cout << "Removed Value " << *it << " due to not being close enough to zero." << std::endl;
                pot_zeros.erase(it);
            } else {
                ++it;
            }
        }

        if (pot_zeros.end() != pot_zeros.begin()){
            pot_zeros.insert(pot_zeros.begin(), a);
            pot_zeros.push_back(b);
            std::sort(pot_zeros.begin(), pot_zeros.end());

            // only needed for early stopping
            /*int size_prev = pot_zeros.size();
            pot_zeros.erase( unique( pot_zeros.begin(), pot_zeros.end() ), pot_zeros.end() );*/

            // check if min or max, WRONG!
            //for (auto it = pot_zeros.begin(); it != pot_zeros.end();){
            //    if (fct_both(*it)(0,1) < 0.0){
            //        std::cout << *it << std::endl;
            //        std::cout << fct_both(*it)(0,1) << std::endl;
            //        pot_zeros.erase(it);
            //    } else {
            //        ++it;
            //    }
            //}
            //another potentially dangerous point for early stopping
            /*int size = pot_zeros.size();
            if ( size_prev > size ) {
                for (auto it = pot_zeros.begin(); it != (pot_zeros.end());) {
                    //std::cout << *it << std::endl;
                    if (std::abs(fct(*it)) < 10*EPS){
                        zeros.push_back(*it);

                    } else  {
                        ++it;
                    }

                }
            }*/
        } else {
            pot_zeros.push_back(a);
            pot_zeros.push_back(b);
        }

        if (pot_zeros.size() > 2){
            min_shrinkage(mu, pot_zeros, init_len);
            for (auto it = pot_zeros.begin(); it != (pot_zeros.end()-1); ++it){
                std::vector<double> temp = findZeros(fct_both, *it, *(it + 1), init_len);
                zeros.insert(zeros.end(),temp.begin(),temp.end());
            }
        } else if (std::abs(*std::prev(pot_zeros.end())-*pot_zeros.begin()) > mu_splitting * init_len){
            double midpoint = *pot_zeros.begin() + std::abs(*std::prev(pot_zeros.end())-*pot_zeros.begin())/2.0;
            std::vector<double> temp1 = findZeros(fct_both, pot_zeros.front(), midpoint, init_len);
            std::vector<double> temp2 = findZeros(fct_both, midpoint, pot_zeros.back(), init_len);
            zeros.insert(zeros.end(),temp1.begin(),temp1.end());
            zeros.insert(zeros.end(),temp2.begin(),temp2.end());
        }
    }

    return zeros;
};

template <typename RECORDER>
std::vector<double> findZeros_seq(std::function<Eigen::MatrixXd(double)> fct_both,
                                  double a,
                                  double b,
                                  unsigned int m,
                                  RECORDER rec) {

    std::cout << std::endl;
    std::cout << "findZeros_seq called" << std::endl;

    // set tolerances
    // jump detection
    double sigma = 40;
    // interpolation tolerance
    double eta = 0.1;
    // interval tolerances
    double tau_abs = 0.0001;
    double tau_rel = 0.001;
    // early stopping, unused
    double tau = 1e-6;
    // minimal shrinkage parameter
    double mu = 0.1;
    double mu_splitting = 0.1;
    // interval parameter for early detection (relation to tau_abs important inorder to guarantee possible detection)
    double gamma_rel_0 = 0.0001;
    double gamma_rel_1 = tau_abs / 3;

    // initialize search grid, function values, derivatives and flags
    std::vector<grid_data> S;
    double interval_len = std::abs(b - a);
    double sub_interval_len = interval_len / m;
    for (unsigned int i = 0; i <= m; i++) {
        Eigen::MatrixXd tempM = fct_both(a + i * sub_interval_len);
        grid_data data_field = {a + i * sub_interval_len, tempM(0, 0), tempM(0, 1), active};
        S.push_back(data_field);
    }
    S.back().flag_val = nozero;
    int N_active = S.size() - 1;
    std::cout << "Data initialized as:" << std::endl;
    std::cout << S << std::endl;
    rec(S);

    do {
        // initialize storage for all real zeros found
        std::vector<double> zeros;
        std::vector<double> abs_values;
        for (auto it = S.begin(); it != S.end(); ++it) {
            abs_values.push_back(std::abs((*it).value));
        }
        std::cout << "absolute values" << std::endl;
        std::cout << abs_values << std::endl;
        double M = *std::max_element(abs_values.begin(), abs_values.end());
        std::vector<int> Z;

        for (auto it = abs_values.begin(); it != abs_values.end(); ++it) {
            if (*it < M * tau) {
                int i = it - abs_values.begin();
                if (i == 0 && S.at(i).flag_val == active) {
                    Z.push_back(i);
                } else if (i == abs_values.size() - 1 && S.at(i - 1).flag_val == active) {
                    Z.push_back(i);
                } else if (i > 0 && i < abs_values.size()) {
                    if (S.at(i - 1).flag_val == active && S.at(i).flag_val == active) {
                        Z.push_back(i);
                    }
                }
            }
        }
        std::cout << "Zeros found through early checking: " << std::endl;
        std::cout << Z << std::endl;


        for (auto it = Z.end() - 1; it != Z.begin() - 1; --it) {
            double temp_zero = S[*it].grid_point;
            double kappa_l;
            double kappa_r;
            if (*it == 0) {
                kappa_l = temp_zero;
                kappa_r = temp_zero +
                          std::min(0.8*tau_abs,
                                   0.5*std::abs(S[*it + 1].grid_point - temp_zero));
                Eigen::MatrixXd tempM = fct_both(kappa_r);
                double f_kappa_l = S[*it].value;
                double f_kappa_r = tempM(0, 0);
                double df_kappa_l = S[*it].derivative;
                double df_kappa_r = tempM(0, 1);
                grid_data data_field = {kappa_r, f_kappa_r, df_kappa_r, active};
                S.insert(S.begin() + *it + 1, data_field);
                S[0].flag_val = active;
                N_active += 1;

            } else if (*it == S.size() - 1) {
                kappa_l = temp_zero -
                          std::min(0.8*tau_abs,
                                   0.5*std::abs(temp_zero - S[*it - 1].grid_point));
                kappa_r = temp_zero;
                Eigen::MatrixXd tempM = fct_both(kappa_l);
                double f_kappa_l = tempM(0, 0);
                double f_kappa_r = S[*it].value;
                double df_kappa_l = tempM(0, 1);
                double df_kappa_r = S[*it].derivative;
                grid_data data_field = {kappa_l, f_kappa_l, df_kappa_l, active};
                S.insert(S.begin() + *it, data_field);
                N_active += 1;

            } else {
                kappa_l = temp_zero -
                          std::min(0.4*tau_abs,
                                   0.5*std::abs(temp_zero - S[*it - 1].grid_point));
                kappa_r = temp_zero +
                          std::min(0.4*tau_abs,
                                   0.5*std::abs(S[*it + 1].grid_point - temp_zero));
                Eigen::MatrixXd tempM = fct_both(kappa_l);
                double f_kappa_l = tempM(0, 0);
                double df_kappa_l = tempM(0, 1);
                tempM = fct_both(kappa_r);
                double f_kappa_r = tempM(0, 0);
                double df_kappa_r = tempM(0, 1);
                grid_data data_field = {kappa_l, f_kappa_l, df_kappa_l, active};
                S.erase(S.begin() + *it);
                S.insert(S.begin() + *it, data_field);
                data_field = {kappa_r, f_kappa_r, df_kappa_r, active};
                S.insert(S.begin() + *it + 1, data_field);
                N_active += 1;

            }
        }
        std::cout << "Data after modifications: " << std::endl;
        std::cout << S << std::endl;


        for ( unsigned int i = 1; i < S.size(); i++ ){
            std::vector<double> pot_zeros;
            if ( S[i-1].flag_val == active ) {

                // gather function and it's derivative at boundaries
                double alpha = S[i-1].grid_point;
                double beta = S[i].grid_point;
                double f_alpha = S[i-1].value;
                double f_beta = S[i].value;
                double df_alpha = S[i-1].derivative;
                double df_beta = S[i].derivative;
                std::cout << "a = " << alpha << std::endl;
                std::cout << "b = " << beta << std::endl;
                std::cout << "f(a) = " << f_alpha << std::endl;
                std::cout << "f(b) = " << f_beta << std::endl;
                std::cout << "df(a) = " << df_alpha << std::endl;
                std::cout << "df(b) = " << df_beta << std::endl;

                // check for sign change
                double val = std::abs((f_beta - f_alpha) / (beta - alpha));
                double tol = sigma * std::max(std::abs(df_alpha), std::abs(df_beta));
                double val1 = std::abs(beta - alpha);
                double tol1 = tau_rel * std::min(std::abs(alpha), std::abs(beta));
                double val2 = std::abs(beta - alpha);
                double tol2 = tau_abs;
                // jump detection
                if (f_alpha * f_beta <= 0 && val > tol) {
                    std::cout << "jump detected" << std::endl;
                    std::cout << std::endl;
                    S[i - 1].flag_val = nozero;
                    N_active -= 1;
                // minima detecting
                }else if (f_alpha * f_beta <= 0 && (val1 <= tol1 || val2 <= tol2)) {
                    std::cout << "sign change + small interval + no jump detected" << std::endl;
                    std::cout << std::endl;
                    S[i - 1].flag_val = zerofound;
                    N_active -= 1;
                } else {
                    // hermite polynom using explicit coefficients
                    double d_p = f_alpha;
                    double c_p = (beta - alpha) * df_alpha;

                    Eigen::Matrix2d A;
                    A << 1, 1, 3, 2;
                    Eigen::Vector2d rhs;
                    rhs << f_beta - f_alpha - (beta - alpha) * df_alpha, (beta - alpha) * df_beta -
                                                                         (beta - alpha) * df_alpha;

                    Eigen::Vector2d res = A.householderQr().solve(rhs);
                    double a_p = res[0];
                    double b_p = res[1];
                    std::cout << "Coefficients of Hermite Polynomial: " << std::endl;
                    std::cout << a_p << " " << b_p << " " << c_p << " " << d_p << std::endl;

                    // compute zeros of hermite interpolating function
                    Eigen::VectorXd order_temp = Eigen::VectorXd::Zero(4);
                    order_temp << a_p, b_p, c_p, d_p;
                    order_temp = order_temp.cwiseAbs();
                    std::cout << order_temp.maxCoeff() << std::endl;
                    double order = order_temp.maxCoeff();
                    std::cout << order << std::endl;
                    if (std::abs(a_p) > order * EPS) {
                        pot_zeros = general_cubic_formula(a_p, b_p, c_p, d_p, alpha, beta);
                    } else {
                        if (std::abs(b_p) > order * EPS) {
                            double a_p2 = 1.0;
                            double b_p2 = c_p / (b_p);
                            double c_p2 = d_p / (b_p);
                            std::cout << "fall back to quadratic case" << std::endl;
                            std::cout << "New coefficients: " << std::endl;
                            std::cout << a_p2 << " " << b_p2 << " " << c_p2 << std::endl;
                            pot_zeros = zerosquadpolstab(b_p2, c_p2, alpha, beta);

                        } else {
                            if (std::abs(c_p) > order * EPS) {
                                std::cout << "fall back to linear case" << std::endl;
                                if (-d_p / c_p > 0 && -d_p / c_p < 1) {
                                    std::cout << "Zero found: " << std::endl;
                                    std::cout << -d_p / c_p << std::endl;
                                    pot_zeros.push_back((1 + d_p / c_p) * alpha - (d_p / c_p) * beta);
                                }
                            }
                        }
                    }
                    std::cout << "Zeros of Hermite Polynomial: " << pot_zeros << std::endl;

                    if (pot_zeros.end() == pot_zeros.begin()) {
                        if (std::abs(a_p) > order * EPS) {
                            double a_p2 = 1.0;
                            double b_p2 = 2 * b_p / (3 * a_p);
                            double c_p2 = c_p / (3 * a_p);
                            pot_zeros = zerosquadpolstab(b_p2, c_p2, alpha, beta);
                        } else {
                            if (std::abs(b_p) > order * EPS) {
                                double b_p2 = 1.0;
                                double c_p2 = c_p / (2 * b_p);
                                double pot_zero = -c_p2;
                                std::cout << "fall back to linear case" << std::endl;
                                if (pot_zero > 0 && pot_zero < 1) {
                                    std::cout << "Zero found: " << std::endl;
                                    std::cout << pot_zero << std::endl;
                                    pot_zeros.push_back((1 - pot_zero) * alpha + (pot_zero) * beta);
                                }
                            }
                        }
                        std::cout << "Extrema of Hermite Polynomial: " << pot_zeros << std::endl;

                        for (auto it = pot_zeros.begin(); it != pot_zeros.end();) {
                            double val_01 = (*it - alpha) / (beta - alpha);
                            double val_at_min = a_p * pow(val_01, 3) + b_p * pow(val_01, 2) + c_p * (val_01) + d_p;
                            double der2_at_min = 6 * a_p * (val_01) + 2 * b_p;
                            if (!(der2_at_min * val_at_min >= 0)) {
                                std::cout << "Removed Value " << *it << " due to being positive Max/ negative Min."
                                          << std::endl;
                                pot_zeros.erase(it);
                            } else if (!(std::abs(val_at_min) <
                                         eta * std::min(std::abs(f_alpha), std::abs(f_beta)))) {
                                std::cout << "Removed Value " << *it << " due to not being close enough to zero."
                                          << std::endl;
                                pot_zeros.erase(it);
                            } else {
                                ++it;
                            }
                        }
                        if (pot_zeros.end() == pot_zeros.begin()) {
                            if (std::abs(beta - alpha) > mu_splitting/(1.0+std::abs(alpha)) * sub_interval_len) {
                                std::cout << "Adding midpoint" << std::endl;
                                double midpoint = alpha + std::abs(alpha - beta) / 2.0;
                                pot_zeros.push_back(midpoint);
                            } else if (pot_zeros.end() == pot_zeros.begin()) {
                                S[i - 1].flag_val = nozero;
                                N_active -= 1;
                            }
                        }
                    }
                    std::cout << "Potential zeros: " << std::endl;
                    std::cout << pot_zeros << std::endl;
                    std::cout << std::endl;
                }
            }
            zeros.insert(zeros.end(), pot_zeros.begin(), pot_zeros.end());
        }
        for ( auto it = zeros.begin(); it != zeros.end(); ++it){
            Eigen::MatrixXd tempM = fct_both(*it);
            grid_data data_field = {*it, tempM(0,0), tempM(0,1), active};
            S.push_back(data_field);
            N_active += 1;
        }
        std::sort(S.begin(), S.end(), ordering());
        rec(S);
        std::cout << "active intervals" << std::endl;
        std::cout << N_active << std::endl;
    } while (N_active > 0);

    std::vector<double> final_zeros;

    for ( auto it = S.begin(); it != S.end()-1; ++it ){
        if ( (*it).flag_val == zerofound ){
            final_zeros.push_back((0.5*((*it).grid_point+(*(it+1)).grid_point)));
        }
    }
    rec(S);

    return final_zeros;
};
template
std::vector<double> findZeros_seq<>(std::function<Eigen::MatrixXd(double)> fct_both,
                                    double a,
                                    double b,
                                    unsigned int m,
                                    std::function<void(std::vector<grid_data>)> rec);

void min_shrinkage(double mu, std::vector<double> &pot_zeros, double init_len){
    std::cout << "Applying minimal shrinkage condition to vector: " << pot_zeros << std::endl;
    double min_len = std::abs(pot_zeros.back() - pot_zeros.front());
    for (auto it = pot_zeros.begin(); it != (pot_zeros.end()-1); ++it){
        if (std::abs(*(it+1)-*it) < mu*min_len) {
            if ((it + 1) != std::prev(pot_zeros.end())){
                *(it + 1) = *it + mu * min_len;
            } else {
                *(it) = *(it + 1) - mu * min_len;
            }
        }
    }
    std::cout << "Result: " << pot_zeros << std::endl;
}

std::vector<double> zerosquadpolstab( double b, double c, double x0, double x1) {
    std::vector<double> pot_zeros01(0);
    double D = std::pow(b,2) - 4*c ; // discriminant
    std::vector<double> pot_zeros(0);
    if ( D < 0){
        return pot_zeros;
    } else {
        // two cases for avoiding double zeros
        double wD = std::sqrt(D) ;
        if (std::abs(D) <10*EPS){
            // Use discriminant formula only for zero far away from 0
            // in order to avoid cancellation. For the other zero
            // use Vieta’s formula.
            if ( b >= 0 ) {
                double t = 0.5 * (-b - wD) ;
                pot_zeros01.push_back(t);
            } else {
                double t = 0.5 * (-b + wD);
                pot_zeros01.push_back(c / t);
            }

        } else {
            // Use discriminant formula only for zero far away from 0
            // in order to avoid cancellation. For the other zero
            // use Vieta’s formula.
            if ( b >= 0 ) {
                double t = 0.5 * (-b - wD) ;
                pot_zeros01.push_back(t);
                pot_zeros01.push_back(c / t);
            } else {
                double t = 0.5 * (-b + wD);
                pot_zeros01.push_back(c / t);
                pot_zeros01.push_back(t);
            }
        }
    }
    for (auto it = pot_zeros01.begin(); it != pot_zeros01.end(); ++it){
        double pot_zero01 = *it;
        double pot_zero = (1-pot_zero01) * x0 + (pot_zero01) * x1;
        if (pot_zero <= x1 && pot_zero >= x0 && (std::find(pot_zeros01.begin(), pot_zeros01.end(), pot_zero) == pot_zeros01.end())){
            //std::cout << pot_zero << std::endl;
            //std::cout << std::abs(a*std::pow(pot_zero,3) + b*std::pow(pot_zero,2) + c*pot_zero + d) << std::endl;
            // Catch close from below as well!
            //double val_at_zero = std::pow(pot_zero01,2) + b*pot_zero01 + c;
            //if (std::abs(val_at_zero)< 1e-4 &&
            //       (((2*pot_zero01 + b >= 0) && (val_at_zero > 0))|| ((2*pot_zero01 + b <= 0) && (val_at_zero <= 0)))){
            pot_zeros.push_back(pot_zero);
            //}
        }
    }
    return pot_zeros;
}

std::vector<double> general_cubic_formula(double a, double b, double c, double d, double x0, double x1){
    complex_t epsilon = (complex_t(-1) + std::sqrt(complex_t(-3.0,0.0)))/2.0;
    complex_t delta_0 = std::pow(b, 2) - 3 * a * c;
    complex_t delta_1 = 2*std::pow(b, 3) - 9 * a * b * c
                        + 27 * std::pow(a, 2) * d;
    auto C = [&] (){
        complex_t temp = complex_t(std::pow((delta_1 + std::pow(std::pow(delta_1,2) - complex_t(4.0)*std::pow(delta_0,3),0.5))/complex_t(2.0),1.0/3.0));
        if (std::abs(temp) > 1e-16){
            return temp;
        }
        else {
            return complex_t(std::pow((delta_1 - std::pow(std::pow(delta_1,2) - complex_t(4.0)*std::pow(delta_0,3),0.5))/complex_t(2.0),1.0/3.0));
        }
    };


    // can't handle a = 0
    auto zeros_fct = [&] (int k) {
        return complex_t(-1/(3 * a)) * (complex_t(b) + std::pow(epsilon, k) * C() + complex_t(delta_0) / (C() * complex_t(std::pow(epsilon, k))));
    };


    // find potential candidates for actual zeros
    std::vector<double> pot_zeros;
    for (unsigned  k=0; k<3; k++){
        complex_t pot_zero01_complex = zeros_fct(k);
        // unstable detection!!
        if (pot_zero01_complex.imag() < 1e-10){
            double pot_zero01 = zeros_fct(k).real();
            double pot_zero = (1-pot_zero01) * x0 + (pot_zero01) * x1;
            if (pot_zero <= x1 && pot_zero >= x0 && (std::find(pot_zeros.begin(), pot_zeros.end(), pot_zero) == pot_zeros.end())){
                // unstable detection!!
                //if (std::abs(a*std::pow(pot_zero01,3) + b*std::pow(pot_zero01,2) + c*pot_zero01 + d) < 1e-4){
                pot_zeros.push_back(pot_zero);
                //}
            }
        }
    }
    return pot_zeros;
}
Eigen::VectorXd parabolic_approximation(const std::function<Eigen::VectorXd(double)> f,
                                        const std::function<Eigen::VectorXd(double)> f_der,
                                        const std::function<Eigen::VectorXd(double)> f_der2,
                                        const double x0,
                                        double step){
    // initialize function values and derivatives
    Eigen::VectorXd vals = f(x0);
    Eigen::VectorXd ders = f_der(x0);
    Eigen::VectorXd ders2 = f_der2(x0);
    // initialize vecotor for minima approximations
    unsigned numsvs = vals.size();
    Eigen::VectorXd exts(numsvs);
    //bool adjust_step;
    //step = step/(1+ders2.cwiseAbs().maxCoeff());
    Eigen::VectorXd res(4);
    // do {
    //adjust_step = false;
    // compute approximations for next minima
    for (unsigned i = 0; i < numsvs; i++) {
        // find minima of parabola
        exts[i] = (-ders[i]/ders2[i] + x0);
        // try and stop convergence to maxima
        if ( exts[i] - x0 > step) {
            exts[i] = x0 + step;
        }else if (exts[i] -x0 < -step){
            exts[i] = x0 - step;
        }/* else if( exts[i] - x0 < 0 ){
            if (ders2[i]<0){
                exts[i] = x0+step;
            } else if (exts[i] -x0 < -step){
                exts[i] = x0 - 0.5 * step;
            }
        }*/
    }
    // find best approximation for minima
    res[0] = exts[0];
    res[1] = f(res[0])[0];
    double temp;
    unsigned index = 0;
    for (unsigned i = 1; i < numsvs; i++) {
        temp = f(exts[i])[i];
        if (temp < res[1]) {
            res[1] = temp;
            index = i;
        }
    }
    // return best approx for minima and
    // values used for parabolic approximation
    res[0] = exts[index];
    res[1] = vals[index];
    res[2] = ders[index];
    res[3] = ders2[index];

    /*if (res[1] < 0) {
        adjust_step = true;
        step *= 0.5;
    }
} while(adjust_step);*/
    return res;

}

double secant_method(const function<double(double)> f,
                     double x1,
                     double x2,
                     const double tol,
                     const unsigned maxIter,
                     bool &root_found,
                     unsigned &num_iter) {
    unsigned n = 0;
    double xm = 0.0;
    double x0 = 0.0;
    double c = 0.0;
    double f2 = f(x2);
    double f1 = f(x1);
    if (f1 * f2 < 0) {
        do {
            // calculate the intermediate value
            x0 = (x1 * f2 - x2 * f1) / (f2 - f1);
            double f0 = f(x0);
            // check if x0 is root of equation or not
            c = f1 * f0;
            // update the value of interval
            x1 = x2;
            x2 = x0;
            f1 = f(x1);
            f2 = f(x2);
            // update number of iteration
            n++;
            if (abs(c) < tol)
                break;
            // if x0 is the root of equation then break the loop
            xm = (x1 * f2 - x2 * f1) / (f2 - f1);
        } while (abs(xm - x0) >= tol); // repeat the loop
        // until the convergence
        if (abs(c) < tol) {
            cout << "Root of the given equation = " << x0 << endl;
            cout << "No. of iterations = " << n << endl;
            root_found = true;
            num_iter = n;
            return x0;
        } else {
            cout << "No root found after " << n << " iterations." << endl;
        }
    } else
        cout << "There might not be a root in this interval." << endl;
    return 0.0;
}
