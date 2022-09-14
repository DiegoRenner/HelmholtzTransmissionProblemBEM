/**
 * \file run_config.cpp
 * \brief This target builds a script that runs a configuration
 * as defined in config.json. It is meant to make running jobs
 * in this library more straightforward.In the future one might be able
 * to just setup a config file instead of writing your own script.
 * This is still in development.
 *
 * This File is a part of the HelmholtzTransmissionProblemBEM library.
 */

#include <json/value.h>
#include <json/reader.h>
#include <json/writer.h>
#include <fstream>
#include <iostream>
#include "first_kind_direct_dirichlet.hpp"

int main(int argc, char** argv){
    std::ifstream config_file("../examples/config.json");
    //char* a;
    //config_file.read(a,100);
    //std::cout << a << std::endl;
    Json::Value config;
    config_file >> config;
    unsigned n_jobs = config["jobs"].size();
    for (int i = 0; i < n_jobs; i ++){
        Json::Value job_config = config["jobs"][i];
        FirstKindDirectDirichlet firstKindDirectDirichlet(job_config);
        firstKindDirectDirichlet.run();
    }
    std::cout << config["jobs"][0].size(); //This will print the entire json object.


    //The following lines will let you access the indexed objects.
    //std::cout << config["Anna"]; //Prints the value for "Anna"
    //std::cout<<config["ben"]; //Prints the value for "Ben"
    //std::cout<<config["Anna"]["profession"]; //Prints the value corresponding to "profession" in the json for "Anna"

    //std::cout<<config["profession"]; //NULL! There is no element with key "profession". Hence a new empty element will be created.
    return 0;
}