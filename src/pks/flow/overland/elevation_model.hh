/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  The elevation model gets the surface elevation, slope, and updates pres + elev.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_ELEVATION_MODEL_
#define AMANZI_FLOWRELATIONS_ELEVATION_MODEL_

#include "secondary_variables_field_model.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class ElevationModel : public SecondaryVariablesFieldModel {

 public:
  ElevationModel();

  // Required methods from SecondaryVariableFieldModel
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results);

  virtual void EvaluateElevationAndSlope_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results) = 0;

 private:
  bool updated_once_;
  Key pres_key_;
};

} //namespace
} //namespace
} //namespace

#endif
