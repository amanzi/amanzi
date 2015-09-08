/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <string>
#include "MultiscaleTransportPorosity_DPM.hh"

namespace Amanzi {
namespace Transport {

/* ******************************************************************
* This model is minor extension of the WRM.
****************************************************************** */
MultiscaleTransportPorosity_DPM::MultiscaleTransportPorosity_DPM(Teuchos::ParameterList& plist)
{
  omega_ = plist.get<double>("solute transfer coefficient", 0.0);
}


/* ******************************************************************
* It should be called only once; otherwise, create an evaluator.
****************************************************************** */
double MultiscaleTransportPorosity_DPM::ComputeSoluteFlux(
    double flux_liquid, double tcc_f, double tcc_m)
{
  double tmp = (flux_liquid > 0.0) ? tcc_f : tcc_m; 
  return flux_liquid * tmp + omega_ * (tcc_f - tcc_m);
}


/* ******************************************************************
* It should be called only once; otherwise, create an evaluator.
****************************************************************** */
void MultiscaleTransportPorosity_DPM::UpdateStabilityOutflux(
    double flux_liquid, double* outflux)
{ 
  *outflux += std::max(0.0, flux_liquid) + omega_;
}

}  // namespace Transport
}  // namespace Amanzi
  
  
