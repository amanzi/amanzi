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

namespace Amanzi {
namespace AmanziFlow {

class TI_Specs {
 public:
  TI_Specs() { 
    ti_method = 0;
    preconditioner_name = "";
    num_itrs = max_itrs = 0;
    T0 = T1 = dT0 = dTmax = 0.0;
    atol = 1.0; 
    rtol = residual_tol = 0.0;
  }

 public:
  int ti_method;
  std::string preconditioner_name;
  int num_itrs, max_itrs;

  double T0, T1, dT0, dTmax;
  double atol, rtol;  // obsolete options
  double residual_tol; 
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif

