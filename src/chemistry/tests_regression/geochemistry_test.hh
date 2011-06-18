/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_GEOCHEMISTRY_TEST_HH_
#define AMANZI_CHEMISTRY_GEOCHEMISTRY_TEST_HH_

#include <string>
#include <vector>

#include "beaker.hh"
#include "verbosity.hh"

int CommandLineOptions(int argc, char** argv,
                       amanzi::chemistry::Verbosity* verbosity,
                       int* test, std::string* model);
void ModelSpecificParameters(const std::string model,
                             amanzi::chemistry::Beaker::BeakerParameters* parameters);
void PrintDoubleVector(const std::vector<double> &total);

void calcite_kinetics(const amanzi::chemistry::Verbosity& verbosity,
                      std::string* thermo_database_file,
                      std::string* activity_model_name,
                      amanzi::chemistry::Beaker::BeakerComponents* components,
                      double* delta_time,
                      int* num_time_steps,
                      int* output_interval);

void calcite_kinetics_large_time_steps(const amanzi::chemistry::Verbosity& verbosity,
                                       std::string* thermo_database_file,
                                       std::string* activity_model_name,
                                       amanzi::chemistry::Beaker::BeakerComponents* components,
                                       double* delta_time,
                                       int* num_time_steps,
                                       int* output_interval);

void surface_complexation(const amanzi::chemistry::Verbosity& verbosity,
                          std::string* thermo_database_file,
                          std::string* activity_model_name,
                          amanzi::chemistry::Beaker::BeakerComponents* components,
                          double* delta_time,
                          int* num_time_steps,
                          int* output_interval);

void fbasin_initial(const amanzi::chemistry::Verbosity& verbosity,
                    std::string* thermo_database_file,
                    std::string* activity_model_name,
                    amanzi::chemistry::Beaker::BeakerComponents* components,
                    double* delta_time,
                    int* num_time_steps,
                    int* output_interval);
void fbasin_infiltration(const amanzi::chemistry::Verbosity& verbosity,
                         std::string* thermo_database_file,
                         std::string* activity_model_name,
                         amanzi::chemistry::Beaker::BeakerComponents* components,
                         double* delta_time,
                         int* num_time_steps,
                         int* output_interval);
void fbasin_source(const amanzi::chemistry::Verbosity& verbosity,
                   std::string* thermo_database_file,
                   std::string* activity_model_name,
                   amanzi::chemistry::Beaker::BeakerComponents* components,
                   double* delta_time,
                   int* num_time_steps,
                   int* output_interval);
void fbasin_aqueous_initial(amanzi::chemistry::Beaker::BeakerComponents* components);
void fbasin_aqueous_infiltration(amanzi::chemistry::Beaker::BeakerComponents* components);
void fbasin_aqueous_source(amanzi::chemistry::Beaker::BeakerComponents* components);
void fbasin_free_ions(amanzi::chemistry::Beaker::BeakerComponents* components);
void fbasin_minerals(amanzi::chemistry::Beaker::BeakerComponents* components);
void fbasin_sorbed(amanzi::chemistry::Beaker::BeakerComponents* components);

void uo2_5_component_initial(const amanzi::chemistry::Verbosity& verbosity,
                             std::string* thermo_database_file,
                             std::string* activity_model_name,
                             amanzi::chemistry::Beaker::BeakerComponents* components,
                             double* delta_time,
                             int* num_time_steps,
                             int* output_interval);
void uo2_5_component_outlet(const amanzi::chemistry::Verbosity& verbosity,
                            std::string* thermo_database_file,
                            std::string* activity_model_name,
                            amanzi::chemistry::Beaker::BeakerComponents* components,
                            double* delta_time,
                            int* num_time_steps,
                            int* output_interval);
void uo2_5_component_source(const amanzi::chemistry::Verbosity& verbosity,
                            std::string* thermo_database_file,
                            std::string* activity_model_name,
                            amanzi::chemistry::Beaker::BeakerComponents* components,
                            double* delta_time,
                            int* num_time_steps,
                            int* output_interval);
void uo2_5_component_aqueous_initial(amanzi::chemistry::Beaker::BeakerComponents* components);
void uo2_5_component_aqueous_outlet(amanzi::chemistry::Beaker::BeakerComponents* components);
void uo2_5_component_aqueous_source(amanzi::chemistry::Beaker::BeakerComponents* components);
void uo2_5_component_free_ions(amanzi::chemistry::Beaker::BeakerComponents* components);
void uo2_5_component_minerals(amanzi::chemistry::Beaker::BeakerComponents* components);
void uo2_5_component_sorbed(amanzi::chemistry::Beaker::BeakerComponents* components);


#endif  /* AMANZI_CHEMISTRY_GEOCHEMISTRY_TEST_HH_ */
