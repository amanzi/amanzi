/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_MULTIPHASE_TI_SPECS_HH_
#define AMANZI_MULTIPHASE_TI_SPECS_HH_

#include "MultiphaseDefs.hh"

namespace Amanzi {
namespace Multiphase {

class TI_Specs {
 public:
  TI_Specs() { 
    ti_method = MULTIPHASE_TIME_INTEGRATION_BDF1;
    ti_method_name = "none";
    preconditioner_name = "undefined";
    num_itrs = max_itrs = 0;
    error_control_options = 0;
    dT_method = 0;
    T0 = dT0 = dTmax = 0.0;
    T1 = 100.0;  // used in unit tests
    dTfactor = 1.0;
    atol = 1.0; 
    rtol = 1e-5;
    residual_tol = 0.0;
    initialize_with_darcy = false;
    clip_saturation = -1.0;
    clip_pressure = -1e+10;
    pressure_lambda_constraints = true;
    inflow_krel_correction = false;
    ti_list_ptr_ = NULL;
  }
  ~TI_Specs() {};

 public:
  int ti_method;
  std::string ti_method_name;

  std::string preconditioner_name;
  std::string solver_name;
  int num_itrs, max_itrs;

  int dT_method, error_control_options;
  double T0, T1, dT0, dTmax, dTfactor;
  double atol, rtol, residual_tol; 

  bool initialize_with_darcy;  // initialization options
  double clip_saturation, clip_pressure;
  std::string solver_name_ini;
  std::string preconditioner_name_ini;

  bool pressure_lambda_constraints, inflow_krel_correction;
  std::string solver_name_constraint;
  std::string preconditioner_name_constraint;

  std::vector<std::pair<double, double> > dT_history;  // statistics (relocate to debug?)
  Teuchos::ParameterList* ti_list_ptr_;
};

}  // namespace Flow
}  // namespace Amanzi

#endif

