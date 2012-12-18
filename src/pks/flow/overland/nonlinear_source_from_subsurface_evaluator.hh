/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  An evaluator for pulling the darcy flux, at the surface, from the
  subsurface field and putting it into a surface field.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_NONLINEAR_SOURCE_FROM_SUBSURFACE_EVALUATOR_HH_
#define AMANZI_FLOWRELATIONS_NONLINEAR_SOURCE_FROM_SUBSURFACE_EVALUATOR_HH_

#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {

// forward declaration
namespace Operators { class MatrixMFD; }

namespace Flow {
namespace FlowRelations {

class NonlinearSourceFromSubsurfaceEvaluator :
    public SecondaryVariableFieldEvaluator {

 public:
  explicit
  NonlinearSourceFromSubsurfaceEvaluator(Teuchos::ParameterList& plist);

  NonlinearSourceFromSubsurfaceEvaluator(const NonlinearSourceFromSubsurfaceEvaluator& other);

  Teuchos::RCP<FieldEvaluator> Clone() const;

  Teuchos::RCP<Operators::MatrixMFD> get_operator() { return op_; }
  void set_operator(const Teuchos::RCP<Operators::MatrixMFD> op) { op_ = op; }

protected:
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  void IdentifyFaceAndDirection_(const Teuchos::Ptr<State>& S);

  typedef std::pair<int, int> FaceDir;
  Teuchos::RCP<std::vector<FaceDir> > face_and_dirs_;
  Teuchos::RCP<Operators::MatrixMFD> op_;

  Key pressure_key_;
  Key density_key_;
  Key height_key_;

  Key surface_mesh_key_;
  Key subsurface_mesh_key_;

};

} //namespace
} //namespace
} //namespace

#endif
