/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_POROSITY_MODEL_COMPRESSIBLE_HH_
#define AMANZI_POROSITY_MODEL_COMPRESSIBLE_HH_

#include "Teuchos_ParameterList.hpp"

#include "PorosityModel.hh"

namespace Amanzi {
namespace Flow {

class PorosityModel_Compressible : public PorosityModel {
 public:
  explicit PorosityModel_Compressible(Teuchos::ParameterList& plist)
  {
    porosity_ = plist.get<double>("undeformed soil porosity");
    p_ref_ = plist.get<double>("reference pressure");
    c_ = plist.get<double>("pore compressibility");

    factor_ = porosity_ * c_;
  }
  ~PorosityModel_Compressible(){};

  // required methods from the base class
  inline double Porosity(double p)
  {
    double dp = p - p_ref_;
    return (dp <= 0.0) ? porosity_ : porosity_ * std::exp(c_ * dp);
  }
  inline double dPorositydPressure(double p)
  {
    double dp = p - p_ref_;
    return (dp <= 0.0) ? 0.0 : factor_ * std::exp(c_ * dp);
  }

 private:
  double porosity_, p_ref_, c_;
  double factor_;
};

} // namespace Flow
} // namespace Amanzi

#endif
