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

void surface_complexation(const Verbosity& verbosity,
                          std::string* thermo_database_file,
                          std::string* activity_model_name,
                          Beaker::BeakerComponents* components,
                          double* delta_time,
                          int* num_time_steps,
                          int* output_interval);

void surface_complexation_full(const Verbosity& verbosity,
                          std::string* thermo_database_file,
                          std::string* activity_model_name,
                          Beaker::BeakerComponents* components,
                          double* delta_time,
                          int* num_time_steps,
                          int* output_interval);

void fbasin_initial(const Verbosity& verbosity,
                    std::string* thermo_database_file,
                    std::string* activity_model_name,
                    Beaker::BeakerComponents* components,
                    double* delta_time,
                    int* num_time_steps,
                    int* output_interval);
void fbasin_infiltration(const Verbosity& verbosity,
                         std::string* thermo_database_file,
                         std::string* activity_model_name,
                         Beaker::BeakerComponents* components,
                         double* delta_time,
                         int* num_time_steps,
                         int* output_interval);
void fbasin_source(const Verbosity& verbosity,
                   std::string* thermo_database_file,
                   std::string* activity_model_name,
                   Beaker::BeakerComponents* components,
                   double* delta_time,
                   int* num_time_steps,
                   int* output_interval);
void fbasin_aqueous_initial(Beaker::BeakerComponents* components);
void fbasin_aqueous_infiltration(Beaker::BeakerComponents* components);
void fbasin_aqueous_source(Beaker::BeakerComponents* components);
void fbasin_free_ions(Beaker::BeakerComponents* components);
void fbasin_minerals(Beaker::BeakerComponents* components);
void fbasin_sorbed(Beaker::BeakerComponents* components);

void uo2_5_component_initial(const Verbosity& verbosity,
                    std::string* thermo_database_file,
                    std::string* activity_model_name,
                    Beaker::BeakerComponents* components,
                    double* delta_time,
                    int* num_time_steps,
                    int* output_interval);
void uo2_5_component_outlet(const Verbosity& verbosity,
                         std::string* thermo_database_file,
                         std::string* activity_model_name,
                         Beaker::BeakerComponents* components,
                         double* delta_time,
                         int* num_time_steps,
                         int* output_interval);
void uo2_5_component_source(const Verbosity& verbosity,
                   std::string* thermo_database_file,
                   std::string* activity_model_name,
                   Beaker::BeakerComponents* components,
                   double* delta_time,
                   int* num_time_steps,
                   int* output_interval);
void uo2_5_component_aqueous_initial(Beaker::BeakerComponents* components);
void uo2_5_component_aqueous_outlet(Beaker::BeakerComponents* components);
void uo2_5_component_aqueous_source(Beaker::BeakerComponents* components);
void uo2_5_component_free_ions(Beaker::BeakerComponents* components);
void uo2_5_component_minerals(Beaker::BeakerComponents* components);
void uo2_5_component_sorbed(Beaker::BeakerComponents* components);


#endif 
