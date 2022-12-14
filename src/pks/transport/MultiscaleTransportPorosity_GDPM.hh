/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Transport PK

  A generalized dual porosity model, fracture + multiple matrix nodes.
  Current naming convention is that the fields used in the single-porosity
  model correspond now to the fracture continuum.
  Example: tcc = total component concentration in the fracture continuum;
           tcc_matrix = total component concentration in the matrix continuum.
*/

#ifndef MULTISCALE_TRANSPORT_POROSITY_GDPM_HH_
#define MULTISCALE_TRANSPORT_POROSITY_GDPM_HH_

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "Mini_Diffusion1D.hh"

#include "MultiscaleTransportPorosity.hh"

namespace Amanzi {
namespace Transport {

class MultiscaleTransportPorosity_GDPM : public MultiscaleTransportPorosity {
 public:
  MultiscaleTransportPorosity_GDPM(Teuchos::ParameterList& plist);
  ~MultiscaleTransportPorosity_GDPM(){};

  // Compute solute flux: icomp - component id, phi - matrix porosity,
  // tcc_m_aux - vector of concentration values in secondary nodes,
  // wfm[0|1] - fracture water content at initial and final time moments,
  // wcm[0|1] - water content at initial and final time moments
  virtual double ComputeSoluteFlux(double flux_liquid,
                                   double& tcc_f,
                                   WhetStone::DenseVector& tcc_m,
                                   int icomp,
                                   double dt,
                                   double wcf0,
                                   double wcf1,
                                   double wcm0,
                                   double wcm1,
                                   double phi) override;

  // Number of matrix nodes
  virtual int NumberMatrixNodes() override { return matrix_nodes_; }

 private:
  static Utils::RegisteredFactory<MultiscaleTransportPorosity, MultiscaleTransportPorosity_GDPM>
    factory_;
  std::vector<Operators::Mini_Diffusion1D> op_diff_;

  int matrix_nodes_;
  double depth_, tau_;
  std::vector<double> mol_diff_;
};

} // namespace Transport
} // namespace Amanzi

#endif
