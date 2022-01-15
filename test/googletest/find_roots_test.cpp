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
    //findZeros_seq(f, df, 0, 1, 1, 1);
    double a = 0.0;
    double b = 1.0;
    double init_len = std::abs(b-a);
    int init_resolution = 10;
    std::cout << std::endl;
    std::cout << findZeros_seq(df, a, b, init_resolution);


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
    double init_len = std::abs(b-a);
    int init_resolution = 10;
    std::cout << findZeros_seq(ddf1, a, b, init_resolution);


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
    double init_len = std::abs(b-a);
    int init_resolution = 10;
    std::vector<double> zeros = findZeros_seq(ddf2, a, b, init_resolution);
    std::cout << std::endl;
    std::cout << "Roots found: " << zeros << std::endl;

    ASSERT_TRUE(zeros.size() == 6);
    for ( auto it = zeros.begin(); it != zeros.end(); ++it){
        ASSERT_NEAR(df2(*it), 0.0, 1e-3);
    }
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
    std::vector<double> zeros = findZeros_seq(ddf3, a, b, init_resolution);
    std::cout << std::endl;
    std::cout << "Roots found: " << zeros << std::endl;

    ASSERT_TRUE(zeros.size() == 6);
    for ( auto it = zeros.begin(); it != zeros.end(); ++it){
        ASSERT_NEAR(df3(*it), 0.0, 1e-3);
    }
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
    double init_len = std::abs(b-a);
    int init_resolution = 10;
    std::vector<std::vector<data>> records;
    //auto recorder = [](std::vector<data> entry)->void {
    //    //records.push_back(entry);
    //};
    std::function<void(std::vector<data>)> recorder = [&records](std::vector<data> entry)->void{records.push_back(entry);};
    std::vector<double> zeros = findZeros_seq(ddf3, a, b, init_resolution, recorder);
    std::cout << std::endl;
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


    ASSERT_TRUE(zeros.size() == 6);
    for ( auto it = zeros.begin(); it != zeros.end(); ++it){
        ASSERT_NEAR(df3(*it), 0.0, 1e-3);
    }

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
    std::cout << std::endl;
    std::cout << "Roots found: " << zeros << std::endl;

    ASSERT_TRUE(zeros.size() == 4);
    for ( auto it = zeros.begin(); it != zeros.end(); ++it){
        ASSERT_NEAR(df1(*it), 0.0, 1e-4);
    }
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
    std::cout << std::endl;
    std::cout << "Roots found: " << zeros << std::endl;

    ASSERT_TRUE(zeros.size() == 6);
    for ( auto it = zeros.begin(); it != zeros.end(); ++it){
        ASSERT_NEAR(df2(*it), 0.0, 1e-3);
    }
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
    //findZeros_seq(f, df, 0, 1, 1, 1);
    double a = -2*M_PI;
    double b = 2*M_PI;
    double init_len = std::abs(b-a);
    std::vector<double> zeros = findZeros(ddf3, a, b, init_len);
    std::cout << std::endl;
    std::cout << "Roots found: " << zeros << std::endl;

    ASSERT_TRUE(zeros.size() == 6);
    for ( auto it = zeros.begin(); it != zeros.end(); ++it){
        ASSERT_NEAR(df3(*it), 0.0, 1e-3);
    }
}
