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

#include "DenseVector.hh"

namespace Amanzi {
namespace Flow {

class MultiscaleFlowPorosity {
 public:
  virtual ~MultiscaleFlowPorosity() {};

  // Field here is the water storage.
  // There is no need to use evaluators for this task.
  virtual double ComputeField(double phi, double n_l, double pcm) = 0;

  // local (cell-based) solver returns water storage and capilalry
  // pressure in the matrix. 
  virtual double WaterContentMatrix(
      double pcf0, WhetStone::DenseVector& pcm,
      double wcm0, double dt, double phi, double n_l, int& max_itrs) = 0;

  // Number of matrix nodes
  virtual int NumberMatrixNodes() = 0;
};

}  // namespace Flow
}  // namespace Amanzi
  
#endif
  
