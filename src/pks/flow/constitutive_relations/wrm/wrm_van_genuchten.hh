/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! WRMVanGenuchten : water retention model using van Genuchten's parameterization

/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!
  <ul>Native Spec Example</>
    <ParameterList name="moss" type="ParameterList">
      <Parameter name="region" type="string" value="moss" />
      <Parameter name="WRM Type" type="string" value="van Genuchten" />
      <Parameter name="van Genuchten alpha" type="double" value="0.002" />
      <Parameter name="van Genuchten m" type="double" value="0.2" />
      <Parameter name="residual saturation" type="double" value="0.0" />
      <Parameter name="smoothing interval width" type="double" value="100.0" />
    </ParameterList>

*/

#ifndef ATS_FLOWRELATIONS_WRM_VAN_GENUCHTEN_
#define ATS_FLOWRELATIONS_WRM_VAN_GENUCHTEN_

#include "Teuchos_ParameterList.hpp"
#include "wrm.hh"
#include "factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class WRMVanGenuchten : public WRM {

public:
  explicit WRMVanGenuchten(Teuchos::ParameterList& plist);

  // required methods from the base class
  double k_relative(double pc);
  double d_k_relative(double pc);
  double saturation(double pc);
  double d_saturation(double pc);
  double capillaryPressure(double saturation);
  double d_capillaryPressure(double saturation);
  double residualSaturation() { return sr_; }

 private:
  void InitializeFromPlist_();

  Teuchos::ParameterList& plist_;

  double m_;  // van Genuchten parameters: m, n, alpha
  double n_;
  double l_;
  double alpha_;
  double sr_;  // van Genuchten residual saturation

  int function_;
  double pc0_;  // regularization threshold (usually 0 to 500 Pa)
  double a_, b_, factor_dSdPc_;  // frequently used constant

  static Utils::RegisteredFactory<WRM,WRMVanGenuchten> factory_;
};

} //namespace
} //namespace
} //namespace

#endif
