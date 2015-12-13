/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Base class for continuum multiscale porosity models. 
*/

#ifndef MULTISCALE_FLOW_POROSITY_HH_
#define MULTISCALE_FLOW_POROSITY_HH_

namespace Amanzi {
namespace Flow {

class MultiscaleFlowPorosity {
 public:
  virtual ~MultiscaleFlowPorosity() {};

  // Field here is the water content.
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
  
