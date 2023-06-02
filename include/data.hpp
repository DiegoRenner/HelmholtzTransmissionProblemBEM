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
#ifndef DATA_HPP
#define DATA_HPP

#include <vector>
#include <iomanip>
#include <ostream>

using namespace std;

enum flag {
    active,
    nozero,
    zerofound
};

struct grid_data {
public:
    double grid_point;
    double value;
    double derivative;
    flag flag_val;

};
#endif //OPERATORS
