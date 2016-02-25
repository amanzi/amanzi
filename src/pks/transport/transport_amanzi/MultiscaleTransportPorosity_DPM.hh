/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  A two-scale porosity model (fracture + matrix) aka dual porosity
  model. Current naming convention is that the fields used in the 
  single-porosity model correspond now to the fracture continuum.
  Example: tcc = total component concentration in the fracture continuum;
           tcc_matrix = total component concentration in the matrix continuum.
*/

#ifndef MULTISCALE_TRANSPORT_POROSITY_DPM_HH_
#define MULTISCALE_TRANSPORT_POROSITY_DPM_HH_

#include "Teuchos_ParameterList.hpp"

#include "factory.hh"

#include "MultiscaleTransportPorosity.hh"

namespace Amanzi {
namespace Transport {

class MultiscaleTransportPorosity_DPM : public MultiscaleTransportPorosity {
 public:
  MultiscaleTransportPorosity_DPM(Teuchos::ParameterList& plist);
  ~MultiscaleTransportPorosity_DPM() {};

  // Advances concentrations in the matrix continuum to the next time level
  double ComputeSoluteFlux(double flux_liquid, double tcc_f, double tcc_m);

  // Modify outflux used in the stability estimate.
  void UpdateStabilityOutflux(double flux_liquid, double* outflux);

 private:
  double omega_;

  static Utils::RegisteredFactory<MultiscaleTransportPorosity, MultiscaleTransportPorosity_DPM> factory_;
};

}  // namespace Transport
}  // namespace Amanzi
  
#endif
  
