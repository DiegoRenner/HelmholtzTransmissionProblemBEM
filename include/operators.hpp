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
#ifndef OPERATORS_HPP
#define OPERATORS_HPP

#include <vector>
#include <iomanip>
#include <ostream>
#include "data.hpp"

using namespace std;

inline ostream& operator<<(ostream& os, grid_data data_val){
    os << data_val.grid_point << " ";
    os << data_val.value << " ";
    os << data_val.derivative << " ";
    os << data_val.flag_val;
    return os;
};
template<typename type> ostream& operator<<(ostream& os, std::vector<type> in){
    for (auto it = in.begin(); it != in.end(); ++it){
        os << *it << " ";
    }
    return os;
};

class formatted_output
{
private:
    int width;
    ostream& stream_obj;

public:
    formatted_output(ostream& obj, int w): width(w), stream_obj(obj) {};

    template<typename T>
    formatted_output& operator<<(const T& output)
    {
        stream_obj << setw(width) << output;

        return *this;
    }

    formatted_output& operator<<(ostream& (*func)(ostream&))
    {
        func(stream_obj);
        return *this;
    }
};

inline ostream& operator<<(ostream& os, std::vector<grid_data> in){

    formatted_output output(os, 15);
    output << "Grid Points" << "Value" << "Derivative" << "Flag Value" << std::endl;
    for (auto it = in.begin(); it != in.end(); ++it){
        output << (*it).grid_point;
        output << (*it).value;
        output << (*it).derivative;
        output << (*it).flag_val;
        output << std::endl;
    }
    return os;
};

struct ordering {
    bool operator()(grid_data one, grid_data two) {
        if (std::abs(one.grid_point-two.grid_point) < std::max(std::abs(one.grid_point),std::abs(two.grid_point))*1e-16){
            return one.flag_val > two.flag_val;
        } else {
            return one.grid_point < two.grid_point;
        }
    };
};

inline Eigen::VectorXd stdToEig(std::vector<double> in ){
    Eigen::VectorXd out = Eigen::VectorXd::Zero(in.size());
    int i = 0;
    for ( auto it = in.begin(); it != in.end(); ++it){
       out[i] = *it;
       ++i;
    }
    return out;
}


#endif //OPERATORS
