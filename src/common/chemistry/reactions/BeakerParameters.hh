/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Supporting structs for beaker class.
*/

#ifndef AMANZI_CHEMISTRY_BEAKER_PARAMETERS_HH_
#define AMANZI_CHEMISTRY_BEAKER_PARAMETERS_HH_

#include <string>


namespace Amanzi {
namespace AmanziChemistry {

struct BeakerParameters {
  BeakerParameters()
    : tolerance(1.0e-12),
      max_iterations(250),
      update_activity_newton(false),
      activity_model_name("unit") {};

  // solver parameters
  double tolerance;
  int max_iterations;
  bool update_activity_newton;

  // models
  std::string activity_model_name;

  // Name of the Pitzer virial coefficients database
  std::string pitzer_database;
  // Name of the approach for J's functions for the Pitzer model
  std::string pitzer_jfunction;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi

#endif
