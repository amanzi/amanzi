/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  An evaluator for pulling the darcy flux, at the surface, from the
  subsurface field and putting it into a surface field.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_OVERLAND_SOURCE_FROM_SUBSURFACE_FLUX_EVALUATOR_HH_
#define AMANZI_FLOWRELATIONS_OVERLAND_SOURCE_FROM_SUBSURFACE_FLUX_EVALUATOR_HH_

#include "field_evaluator_factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class OverlandSourceFromSubsurfaceFluxEvaluator :
    public SecondaryVariableFieldEvaluator {

 public:
  explicit
  OverlandSourceFromSubsurfaceFluxEvaluator(Teuchos::ParameterList& plist);

  OverlandSourceFromSubsurfaceFluxEvaluator(const OverlandSourceFromSubsurfaceFluxEvaluator& other);

  Teuchos::RCP<FieldEvaluator> Clone() const;

  // custom ensure compatibility as all data is not just on the same components
  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);

protected:
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  void IdentifyFaceAndDirection_(const Teuchos::Ptr<State>& S);

  typedef std::pair<int, double> FaceDir;
  Teuchos::RCP<std::vector<FaceDir> > face_and_dirs_;

  Key flux_key_;
  Key dens_key_;
  bool volume_basis_;

  Key surface_mesh_key_;
  Key subsurface_mesh_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,OverlandSourceFromSubsurfaceFluxEvaluator> fac_;

};

} //namespace
} //namespace
} //namespace

#endif
