/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <string>

#include "boost/math/tools/roots.hpp"

#include "MultiscalePorosity_DPM.hh"
#include "WRMFactory.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* This model is minor extension of the WRM.
****************************************************************** */
MultiscalePorosity_DPM::MultiscalePorosity_DPM(Teuchos::ParameterList& plist)
{
  WRMFactory factory;
  wrm_ = factory.Create(plist);

  alpha_ = plist.get<double>("mass transfer coefficient", 0.0);
}


/* ******************************************************************
* Main capability: cell-based Newton solver.
****************************************************************** */
double MultiscalePorosity_DPM::WaterContentMatrix(
    double dt, double phi, double n_l, double wcm0, double pcf0, double& pcm)
{
  double guess(pcm);
  double left = pcm - 1e+5; 
  double right = pcm + 1e+5; 

  // setup internal parameters for operator()
  pcf0_ = pcf0;
  sat0_ = wcm0 / (phi * n_l);
  alpha_mod_ = alpha_ * dt / (phi * n_l);

  pcm = boost::math::tools::newton_raphson_iterate<MultiscalePorosity_DPM, double>(*this, guess, left, right, 20);
  return wrm_->saturation(pcm) * phi * n_l;
}


/* ******************************************************************
* Funtion evaluation for Newton solver.
****************************************************************** */
boost::math::tuple<double, double> MultiscalePorosity_DPM::operator() (double pc)
{
  double ds = wrm_->saturation(pc) - sat0_;
  double dp = pc - pcf0_;
  double dsdp = wrm_->dSdPc(pc);
  return boost::math::make_tuple(ds - alpha_mod_ * dp, dsdp - alpha_mod_);
}

}  // namespace Flow
}  // namespace Amanzi
  
  
