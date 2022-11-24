/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_BROOKS_COREY_MODEL_HH_
#define AMANZI_BROOKS_COREY_MODEL_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"

#include "WRM.hh"

namespace Amanzi {
namespace Flow {

class WRM_BrooksCorey : public WRM {
 public:
  explicit WRM_BrooksCorey(Teuchos::ParameterList& plist);
  ~WRM_BrooksCorey(){};

  // required methods from the base class
  double k_relative(double pc) const;
  double saturation(double pc) const;
  double dSdPc(double pc) const;
  double capillaryPressure(double saturation) const;
  double residualSaturation() const { return sr_; }
  double dKdPc(double pc) const;

 private:
  void
  Init_(double lambda, double l, double alpha, double sr, std::string& krel_function, double pc0);

 private:
  double lambda_, l_, alpha_; // Brooks and Corey parameters: lambda, alpha
  double sr_;                 // residual saturation

  double pc0_;                        // regularization threshold (usually 0 to 500 Pa)
  double a_, b_, factor_, pc_bubble_; // frequently used constant

  static Utils::RegisteredFactory<WRM, WRM_BrooksCorey> factory_;
};

} // namespace Flow
} // namespace Amanzi

#endif
