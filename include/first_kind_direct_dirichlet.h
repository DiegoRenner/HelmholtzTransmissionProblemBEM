//
// Created by diego on 7/8/21.
//

#ifndef HELMHOLTZTRANSMISSIONPROBLEM_FIRST_KIND_DIRECT_DIRICHLET_H
#define HELMHOLTZTRANSMISSIONPROBLEM_FIRST_KIND_DIRECT_DIRICHLET_H

#include <complex>
#include <iostream>
#include <json/value.h>
#include <json/reader.h>
#include "parametrized_circular_arc.hpp"

class FirstKindDirectDirichlet {
public:
    FirstKindDirectDirichlet(Json::Value job_config);
    void run();

private:
    double eps = 0.25;
    ParametrizedCircularArc boundary;
    double k  = 1.0;
    unsigned order = 11;
    unsigned n_runs = 7;
    unsigned init_panels = 50;
    std::vector<double> numpanels{50, 100, 200, 400, 800, 1600, 3200};
    Eigen::Vector2d ipt;

};


#endif //HELMHOLTZTRANSMISSIONPROBLEM_FIRST_KIND_DIRECT_DIRICHLET_H
