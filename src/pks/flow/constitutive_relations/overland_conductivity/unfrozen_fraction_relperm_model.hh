
/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluates the Kr associated with the unfrozen fraction of water.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_UNFROZEN_FRACTION_KR_MODEL_
#define AMANZI_FLOWRELATIONS_UNFROZEN_FRACTION_KR_MODEL_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class UnfrozenFractionRelPermModel {
 public:
  UnfrozenFractionRelPermModel(Teuchos::ParameterList& list);

  double UnfrozenFractionRelPerm(double uf, double h);
  double DUnfrozenFractionRelPermDUnfrozenFraction(double uf, double h);


 protected:
  Teuchos::ParameterList plist_;
  const double pi_;

};

} // namespace
} // namespace
} // namespace

#endif
