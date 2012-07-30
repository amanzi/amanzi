/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#ifndef __TI_Specs_HPP__
#define __TI_Specs_HPP__

#include "Flow_constants.hpp"

namespace Amanzi {
namespace AmanziFlow {

class TI_Specs {
 public:
  TI_Specs() { 
    ti_method = 0;
    preconditioner_name = "Trilinos ML";
    preconditioner_method = FLOW_PRECONDITIONER_TRILINOS_ML;
    num_itrs = max_itrs = 0;
    T0 = T1 = dT0 = dTmax = 0.0;
    dTfactor = 1.0;
    residual_tol = 0.0;
    initialize_with_darcy = false;
    clip_saturation = 0.6;
  }

 public:
  int ti_method;
  std::string preconditioner_name;
  int preconditioner_method;
  int num_itrs, max_itrs;

  double T0, T1, dT0, dTmax, dTfactor;
  double residual_tol; 

  bool initialize_with_darcy;  // initialization options
  double clip_saturation;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif

