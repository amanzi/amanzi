/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

A two-scale porosity model (fracture + matrix) aka dual porosity
model. Current naming convention is that the fields used in the
single-porosity model correspond now to the fracture continuum.
Example: pressure = pressure in the fracture continuum;
pressure_msp = pressure in the matrix continuum.

* `"mass transfer coefficient`" [double] is the mass transfer coefficient.

* `"tolerance`" [double] defines tolerance for iterative methods used to solve
  secondary equations. Default is 1e-8.

*/

#ifndef MULTISCALE_FLOW_POROSITY_DPM_HH_
#define MULTISCALE_FLOW_POROSITY_DPM_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"

#include "MultiscaleFlowPorosity.hh"
#include "WRM.hh"

namespace Amanzi {
namespace Flow {

class MultiscaleFlowPorosity_DPM : public MultiscaleFlowPorosity {
 public:
  MultiscaleFlowPorosity_DPM(Teuchos::ParameterList& plist);
  ~MultiscaleFlowPorosity_DPM(){};

  // Calculate field water storage assuming pressure equilibrium
  virtual double ComputeField(double phi, double n_l, double pcm) override;

  // local (cell-based) solver returns water storage and capilalry
  // pressure in the matrix. max_itrs is input/output parameter
  virtual double WaterContentMatrix(double pcf0,
                                    WhetStone::DenseVector& pcm,
                                    double wcm0,
                                    double dt,
                                    double phi,
                                    double n_l,
                                    int& max_itrs) override;

  // Number of matrix nodes
  virtual int NumberMatrixNodes() override { return 1; }

 private:
  Teuchos::RCP<WRM> wrm_;
  double alpha_;
  double tol_;

  static Utils::RegisteredFactory<MultiscaleFlowPorosity, MultiscaleFlowPorosity_DPM> factory_;
};

} // namespace Flow
} // namespace Amanzi

#endif
