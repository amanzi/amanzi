/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Base class for continuum multiscale porosity models. 
*/

#ifndef MULTISCALE_TRANSPORT_POROSITY_HH_
#define MULTISCALE_TRANSPORT_POROSITY_HH_

#include <string>
#include <vector>

namespace Amanzi {
namespace Transport {

class MultiscaleTransportPorosity {
 public:
  virtual ~MultiscaleTransportPorosity() {};

  // Compute solute flux: icomp - component id, phi - matrix porosity phi,
  // tcc_m_aux - vector of concentration values in secondary nodes
  virtual double ComputeSoluteFlux(double flux_liquid, double tcc_f, double tcc_m, 
                 int icomp, double phi, std::vector<double>* tcc_m_aux) = 0;

  // Modify outflux used in the stability estimate.
  virtual void UpdateStabilityOutflux(double flux_liquid, double* outflux) = 0;

  // Number of matrix nodes
  virtual int NumberMatrixNodes() = 0;
};

}  // namespace Transport
}  // namespace Amanzi
  
#endif
  
