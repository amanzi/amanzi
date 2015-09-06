/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_MULTISCALE_POROSITY_HH_
#define AMANZI_MULTISCALE_POROSITY_HH_

namespace Amanzi {
namespace Flow {

class MultiscalePorosity {
 public:
  virtual ~MultiscalePorosity() {};

  // Filed here is the water content.
  // There is no need to use evaluators for this task.
  virtual double ComputeField(double phi, double n_l, double pcm) = 0;

  // local (cell-based) solver returns water content and capilalry
  // pressure in the matrix. 
  virtual double WaterContentMatrix(
      double dt, double phi, double n_l, double wcm0, double pcf0,
      double& pcm, int& max_itrs) = 0;
};

}  // namespace Flow
}  // namespace Amanzi
  
#endif
  
