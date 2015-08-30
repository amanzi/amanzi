/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_MULTISCALE_POROSITY_DPM_HH_
#define AMANZI_MULTISCALE_POROSITY_DPM_HH_

#include "boost/math/tools/tuple.hpp"
#include "Teuchos_ParameterList.hpp"

#include "factory.hh"

#include "MultiscalePorosity.hh"
#include "WRM.hh"

namespace Amanzi {
namespace Flow {

class MultiscalePorosity_DPM : public MultiscalePorosity {
 public:
  MultiscalePorosity_DPM(Teuchos::ParameterList& plist);
  ~MultiscalePorosity_DPM() {};

  // local (cell-based) solver
  double WaterContentMatrix(
      double dt, double phi, double n_l, double wcm0, double pcf0, double& pcm);

  // interface to boost
  boost::math::tuple<double, double> operator() (double pc);

 private:
  Teuchos::RCP<WRM> wrm_;
  double alpha_;

  // support of the Newton-Raphson solver
  double sat0_, pcf0_, alpha_mod_;

  static Utils::RegisteredFactory<MultiscalePorosity, MultiscalePorosity_DPM> factory_;
};

}  // namespace Flow
}  // namespace Amanzi
  
#endif
  
