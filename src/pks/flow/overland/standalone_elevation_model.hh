/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  An elevation model getting values from the volumetric mesh.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_STANDALONE_ELEVATION_MODEL_
#define AMANZI_FLOWRELATIONS_STANDALONE_ELEVATION_MODEL_

#include "composite_vector_function.hh"
#include "elevation_model.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class StandaloneElevationModel : public ElevationModel {

 public:
  StandaloneElevationModel(Teuchos::ParameterList& elev_plist);
  StandaloneElevationModel(const StandaloneElevationModel& other);

  Teuchos::RCP<FieldModel> Clone() const;

  virtual void EvaluateElevationAndSlope_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results);

 protected:
  Teuchos::ParameterList elev_plist_;
  Teuchos::RCP<Functions::CompositeVectorFunction> elevation_function_;
  Teuchos::RCP<Functions::CompositeVectorFunction> slope_function_;

};

} //namespace
} //namespace
} //namespace

#endif
