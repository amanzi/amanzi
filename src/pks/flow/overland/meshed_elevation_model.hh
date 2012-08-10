/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  An elevation model getting values from the volumetric mesh.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_MESHED_ELEVATION_MODEL_
#define AMANZI_FLOWRELATIONS_MESHED_ELEVATION_MODEL_

#include "elevation_model.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class MeshedElevationModel : public ElevationModel {

 public:
  MeshedElevationModel();
  MeshedElevationModel(const MeshedElevationModel& other);

  Teuchos::RCP<FieldModel> Clone() const;

  virtual void EvaluateElevationAndSlope_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results);

};

} //namespace
} //namespace
} //namespace

#endif
