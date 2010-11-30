/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __Geochemistry_Test_hpp__
#define __Geochemistry_Test_hpp__

#include <string>

#include "Beaker.hpp"
#include "Verbosity.hpp"

int CommandLineOptions(int argc, char **argv,
                       Verbosity* verbosity, int* test, std::string* model);
void ModelSpecificParameters(const std::string model,
                             Beaker::BeakerParameters* parameters);
void PrintDoubleVector(const std::vector<double> &total);

void calcite_kinetics(const Verbosity& verbosity,
                      std::string* thermo_database_file,
                      std::string* activity_model_name,
                      Beaker::BeakerComponents* components,
                      double* delta_time,
                      int* num_time_steps,
                      int* output_interval);

void calcite_kinetics_large_time_steps(const Verbosity& verbosity,
                                       std::string* thermo_database_file,
                                       std::string* activity_model_name,
                                       Beaker::BeakerComponents* components,
                                       double* delta_time,
                                       int* num_time_steps,
                                       int* output_interval);

void fbasin_initial_all(const Verbosity& verbosity,
                        std::string* thermo_database_file,
                        std::string* activity_model_name,
                        Beaker::BeakerComponents* components);
void fbasin_initial_speciation(const Verbosity& verbosity,
                               std::string* thermo_database_file,
                               std::string* activity_model_name,
                               Beaker::BeakerComponents* components);
void fbasin_minerals(const Verbosity& verbosity,
                     std::string* thermo_database_file,
                     std::string* activity_model_name,
                     Beaker::BeakerComponents* components);


#endif 
