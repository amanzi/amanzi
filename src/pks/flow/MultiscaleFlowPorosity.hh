/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

The list *multiscale models* is the place for various multiscale models.
The list is extension of the list *water retention models*. 
Its ordered by soil regions and includes parameters for the multiscale,
capillary pressure, and relative permebility models.
This list is optional. 

* `"multiscale model`" [string] is the model name. Available options are `"dual porosity`"
  and `"generalized dual porosity`".

* `"xxx parameters`" [sublist] provides parameters for the model specified by variable `"multiscale model`".

* `"water retention model`" [string] specifies a model for the soil.
  The available models are `"van Genuchten`" and `"Brooks Corey`". 
  Parameters for each model are described above.

* `"relative permeability model`" [string] The available options are `"Mualem`" (default) 
  and `"Burdine`".

*/

#ifndef MULTISCALE_FLOW_POROSITY_HH_
#define MULTISCALE_FLOW_POROSITY_HH_

#include "DenseVector.hh"

namespace Amanzi {
namespace Flow {

class MultiscaleFlowPorosity {
 public:
  virtual ~MultiscaleFlowPorosity(){};

  // Field here is the water storage.
  // There is no need to use evaluators for this task.
  virtual double ComputeField(double phi, double n_l, double prm) = 0;

  // local (cell-based) solver returns water storage and pressure in matrix
  virtual 
  WhetStone::DenseVector WaterContentMatrix(double prf0,
                                            WhetStone::DenseVector& prm,
                                            WhetStone::DenseVector& wcm0,
                                            double dt,
                                            double phi,
                                            double n_l,
                                            double mu_l,
                                            double atm_pressure_,
                                            int& max_itrs) = 0;

  // Number of matrix nodes
  virtual int NumberMatrixNodes() = 0;
};

} // namespace Flow
} // namespace Amanzi

#endif
