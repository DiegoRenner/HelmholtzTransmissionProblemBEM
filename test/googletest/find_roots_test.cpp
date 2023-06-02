/**
 * \file double_layer_test.cpp
 * \brief This test file compares Hypersingular BIO for the
 * Helmholtz kernel to a precomputed known solution from
 * a file.
 */
#include <gtest/gtest.h>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include "find_roots.hpp"

TEST(find_roots_test, findZeros_seq_initalization) {

    auto f = [&] (double x) {
        return x;
    };
    auto df = [&] (double x){
        Eigen::MatrixXd df(1,2);
        df << x,1.0;
        return df;
    };

    double a = 0.0;
    double b = 1.0;
    int init_resolution = 10;
    Eigen::VectorXd zeros = stdToEig(findZeros_seq(df, a, b, init_resolution));
    std::cout << std::endl;
    std::cout << zeros;

}
TEST(find_roots_test, findZeros_seq_f1) {

    auto f1 = [&] (double x) {
        return std::min(std::cos(x),sin(x));
    };
    auto df1 = [&] (double x){
        double df = std::cos(x) < std::sin(x) ? -std::sin(x) : std::cos(x);
        return df;
    };
    auto ddf1 = [&] (double x){
        Eigen::MatrixXd both(1,2);
        double df = std::cos(x) < std::sin(x) ? -std::sin(x) : std::cos(x);
        double ddf = std::cos(x) < std::sin(x) ? -std::cos(x) : -std::sin(x);
        both << df, ddf;
        return both;
    };

    double a = -2*M_PI;
    double b = 2*M_PI;
    int init_resolution = 10;

    Eigen::VectorXd zeros = stdToEig(findZeros_seq(ddf1, a, b, init_resolution));

    int num_expected_zeros = 4;
    Eigen::VectorXd expected_zeros(num_expected_zeros);
    expected_zeros << -3.14159,  -1.5708,  3.14159,  4.71239;

    if (zeros.size() == num_expected_zeros){
        ASSERT_NEAR((zeros-expected_zeros).norm(),0.0, 1e-5);
    } else {
        ASSERT_EQ(zeros.size(),num_expected_zeros);
    }

    std::cout << std::endl;
    std::cout << zeros.transpose();

}

TEST(find_roots_test, findZeros_seq_f2) {

    auto f2 = [&] (double x) {
        return std::min(std::cos(2*x+2),2*sin(x));
    };
    auto df2 = [&] (double x){
        double df = std::cos(2*x+2) < 2*sin(x) ? -2*sin(2*x+2) : 2*cos(x);
        return df;
    };
    auto ddf2 = [&] (double x){
        Eigen::MatrixXd both(1,2);
        double df = std::cos(2*x+2) < 2*sin(x) ? -2*sin(2*x+2) : 2*cos(x);
        double ddf = std::cos(2*x+2) < 2*sin(x) ? -4*cos(2*x+2) : -2*sin(x);
        both << df, ddf;
        return both;
    };

    double a = -2*M_PI;
    double b = 2*M_PI;
    int init_resolution = 10;

    Eigen::VectorXd zeros = stdToEig(findZeros_seq(ddf2, a, b, init_resolution));

    int num_expected_zeros = 6;
    Eigen::VectorXd expected_zeros(num_expected_zeros);
    expected_zeros << -5.71025, -4.14161, -1.5708, 0.570816, 2.14157, 4.71239;

    if (zeros.size() == num_expected_zeros){
        ASSERT_NEAR((zeros-expected_zeros).norm(),0.0, 1e-5);
    } else {
        ASSERT_EQ(zeros.size(),num_expected_zeros);
    }

    std::cout << std::endl;
    std::cout << zeros.transpose();

}
TEST(find_roots_test, findZeros_seq_f3) {

    auto f3 = [&] (double x) {
        return std::min(std::cos(2*x+2),2*sin(x+0.1));
    };
    auto df3 = [&] (double x){
        double df = std::cos(2*x+2) < 2*sin(x+0.1) ? -2*sin(2*x+2) : 2*cos(x+0.1);
        return df;
    };
    auto ddf3 = [&] (double x){
        Eigen::MatrixXd both(1,2);
        double df = std::cos(2*x+2) < 2*sin(x+0.1) ? -2*sin(2*x+2) : 2*cos(x+0.1);
        double ddf = std::cos(2*x+2) < 2*sin(x+0.1) ? -4*cos(2*x+2) : -2*sin(x+0.1);
        both << df, ddf;
        return both;
    };

    double a = -2*M_PI;
    double b = 2*M_PI;
    double init_len = std::abs(b-a);
    int init_resolution = 10;

    Eigen::VectorXd zeros = stdToEig(findZeros_seq(ddf3, a, b, init_resolution));

    int num_expected_zeros = 6;
    Eigen::VectorXd expected_zeros(num_expected_zeros);
    expected_zeros << -5.71025, -4.14161, -1.6708, 0.570816, 2.14157, 4.61239;

    if (zeros.size() == num_expected_zeros){
        ASSERT_NEAR((zeros-expected_zeros).norm(),0.0, 1e-5);
    } else {
        ASSERT_EQ(zeros.size(),num_expected_zeros);
    }

    std::cout << std::endl;
    std::cout << zeros.transpose();
}

TEST(find_roots_test, findZeros_seq_recording) {

    unsigned N_fct_calls = 0;
    auto f3 = [&] (double x) {
        return std::min(std::cos(2*x+2),2*sin(x+0.1));
    };
    auto df3 = [&N_fct_calls] (double x){
        double df = std::cos(2*x+2) < 2*sin(x+0.1) ? -2*sin(2*x+2) : 2*cos(x+0.1);
        return df;
    };
    auto ddf3 = [&] (double x){
        N_fct_calls += 1;
        Eigen::MatrixXd both(1,2);
        double df = std::cos(2*x+2) < 2*sin(x+0.1) ? -2*sin(2*x+2) : 2*cos(x+0.1);
        double ddf = std::cos(2*x+2) < 2*sin(x+0.1) ? -4*cos(2*x+2) : -2*sin(x+0.1);
        both << df, ddf;
        return both;
    };

    double a = -2*M_PI;
    double b = 2*M_PI;
    int init_resolution = 10;

    std::vector<std::vector<grid_data>> records;
    std::function<void(std::vector<grid_data>)> recorder = [&records](std::vector<grid_data> entry)->void{records.push_back(entry);};

    std::vector<double> zeros = findZeros_seq(ddf3, a, b, init_resolution, recorder);

    std::ofstream file_out;
    file_out.open("test.dat", std::ios_base::app);
    file_out << "Initial grid: " << std::endl;
    file_out << *records.begin() << std::endl;
    file_out << std::endl;
    for (auto it = records.begin()+1; it != records.end(); ++it){
        int i = it-records.begin();
        file_out << "Grid after " << i << " iteration(s): " << std::endl;
        file_out << *it;
        file_out << std::endl;
    }
    file_out << "#function calls: "  << N_fct_calls << std::endl;
    file_out << "Roots found: " << zeros << std::endl;

}

TEST(find_roots_test, findZeros_f1) {

    auto f1 = [&] (double x) {
        return std::min(std::cos(x),sin(x));
    };
    auto df1 = [&] (double x){
        double df = std::cos(x) < std::sin(x) ? -std::sin(x) : std::cos(x);
        return df;
    };
    auto ddf1 = [&] (double x){
        Eigen::MatrixXd both(1,2);
        double df = std::cos(x) < std::sin(x) ? -std::sin(x) : std::cos(x);
        double ddf = std::cos(x) < std::sin(x) ? -std::cos(x) : -std::sin(x);
        both << df, ddf;
        return both;
    };

    double a = -2*M_PI;
    double b = 2*M_PI;
    double init_len = std::abs(b-a);

    std::vector<double> zeros = findZeros(ddf1, a, b, init_len);

    ASSERT_TRUE(zeros.size() == 4);
    for ( auto it = zeros.begin(); it != zeros.end(); ++it){
        ASSERT_NEAR(df1(*it), 0.0, 1e-4);
    }

    std::cout << std::endl;
    std::cout << "Roots found: " << zeros << std::endl;

}
TEST(find_roots_test, findZeros_f2) {

    auto f2 = [&] (double x) {
        return std::min(std::cos(2*x+2),2*sin(x));
    };
    auto df2 = [&] (double x){
        double df = std::cos(2*x+2) < 2*sin(x) ? -2*sin(2*x+2) : 2*cos(x);
        return df;
    };
    auto ddf2 = [&] (double x){
        Eigen::MatrixXd both(1,2);
        double df = std::cos(2*x+2) < 2*sin(x) ? -2*sin(2*x+2) : 2*cos(x);
        double ddf = std::cos(2*x+2) < 2*sin(x) ? -4*cos(2*x+2) : -2*sin(x);
        both << df, ddf;
        return both;
    };

    double a = -2*M_PI;
    double b = 2*M_PI;
    double init_len = std::abs(b-a);

    std::vector<double> zeros = findZeros(ddf2, a, b, init_len);

    ASSERT_TRUE(zeros.size() == 6);
    for ( auto it = zeros.begin(); it != zeros.end(); ++it){
        ASSERT_NEAR(df2(*it), 0.0, 1e-3);
    }

    std::cout << std::endl;
    std::cout << "Roots found: " << zeros << std::endl;
}
TEST(find_roots_test, findZeros_f3) {

    auto f3 = [&] (double x) {
        return std::min(std::cos(2*x+2),2*sin(x+0.1));
    };
    auto df3 = [&] (double x){
        double df = std::cos(2*x+2) < 2*sin(x+0.1) ? -2*sin(2*x+2) : 2*cos(x+0.1);
        return df;
    };
    auto ddf3 = [&] (double x){
        Eigen::MatrixXd both(1,2);
        double df = std::cos(2*x+2) < 2*sin(x+0.1) ? -2*sin(2*x+2) : 2*cos(x+0.1);
        double ddf = std::cos(2*x+2) < 2*sin(x+0.1) ? -4*cos(2*x+2) : -2*sin(x+0.1);
        both << df, ddf;
        return both;
    };

    double a = -2*M_PI;
    double b = 2*M_PI;
    double init_len = std::abs(b-a);

    std::vector<double> zeros = findZeros(ddf3, a, b, init_len);

    ASSERT_TRUE(zeros.size() == 6);
    for ( auto it = zeros.begin(); it != zeros.end(); ++it){
        ASSERT_NEAR(df3(*it), 0.0, 1e-3);
    }

    std::cout << std::endl;
    std::cout << "Roots found: " << zeros << std::endl;
}
