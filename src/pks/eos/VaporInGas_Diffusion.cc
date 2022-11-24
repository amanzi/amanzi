/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Binary diffusion coefficient.
*/

#include "errors.hh"
#include "VaporInGas_Diffusion.hh"

namespace Amanzi {
namespace AmanziEOS {

VaporInGas_Diffusion::VaporInGas_Diffusion(Teuchos::ParameterList& plist) : EOS_Diffusion(plist)
{
  dref_ = plist.get<double>("reference binary diffusion", 0.282e-4);
  Tref_ = plist.get<double>("reference temperature", 298.0);
}


double
VaporInGas_Diffusion::Diffusion(double T, double p)
{
  return dref_ * std::pow(T / Tref_, 1.8);
}


double
VaporInGas_Diffusion::DDiffusionDT(double T, double p)
{
  return dref_ * (1.8 / Tref_) * std::pow(T / Tref_, 0.8);
}


double
VaporInGas_Diffusion::DDiffusionDp(double T, double p)
{
  return 0.0;
}

} // namespace AmanziEOS
} // namespace Amanzi
