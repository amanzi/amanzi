/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluator for determining height( rho, head )

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "height_evaluator.hh"


namespace Amanzi {
namespace Flow {
namespace FlowRelations {


HeightEvaluator::HeightEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  // my keys are for saturation and rel perm.
  my_key_ = plist_.get<string>("height key", "ponded_depth");
  setLinePrefix(my_key_+std::string(" evaluator"));

  // my dependencies
  dens_key_ = plist_.get<string>("mass density key", "surface_mass_density_liquid");
  dependencies_.insert(dens_key_);

  pres_key_ = plist_.get<string>("head key", "surface_pressure");
  dependencies_.insert(pres_key_);

  gravity_key_ = plist_.get<string>("gravity key", "gravity");
  patm_key_ = plist_.get<string>("atmospheric pressure key", "atmospheric_pressure");

}


HeightEvaluator::HeightEvaluator(const HeightEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    dens_key_(other.dens_key_),
    pres_key_(other.pres_key_),
    gravity_key_(other.gravity_key_),
    patm_key_(other.patm_key_) {}


Teuchos::RCP<FieldEvaluator>
HeightEvaluator::Clone() const {
  return Teuchos::rcp(new HeightEvaluator(*this));
}


void HeightEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  // Ensure my field exists.  Requirements should be already set.
  ASSERT(my_key_ != std::string(""));
  Teuchos::RCP<CompositeVectorFactory> my_fac = S->RequireField(my_key_, my_key_);

  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, true);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);

  // If my requirements have not yet been set, we'll have to hope they
  // get set by someone later.  For now just defer.
  if (my_fac->Mesh() != Teuchos::null) {
    // Create an unowned factory to check my dependencies.  This is done
    // manually here because we do NOT want faces, despite having faces in
    // my_key.  The faces will get updated directly from the mixed field.
    Teuchos::RCP<CompositeVectorFactory> dep_fac =
        Teuchos::rcp(new CompositeVectorFactory());
    dep_fac->set_owned(false);
    dep_fac->SetGhosted(my_fac->Ghosted());
    dep_fac->SetMesh(my_fac->Mesh());
    dep_fac->AddComponent("cell", AmanziMesh::CELL, 1);

    // Loop over my dependencies, ensuring they meet the requirements.
    for (KeySet::const_iterator key=dependencies_.begin();
         key!=dependencies_.end(); ++key) {
      Teuchos::RCP<CompositeVectorFactory> fac = S->RequireField(*key);
      fac->Update(*dep_fac);
    }

    // Recurse into the tree to propagate info to leaves.
    for (KeySet::const_iterator key=dependencies_.begin();
         key!=dependencies_.end(); ++key) {
      S->RequireFieldEvaluator(*key)->EnsureCompatibility(S);
    }
  }


}


void HeightEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pres_key_);

  // this is rather hacky.  surface_pressure is a mixed field vector -- it has
  // pressure on cells and ponded depth on faces.
  // -- copy the faces over directly
  *result->ViewComponent("face",false) = *pres->ViewComponent("face",false);

  // -- cells need the function eval
  const Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);
  const Epetra_MultiVector& pres_c = *pres->ViewComponent("cell",false);
  const Epetra_MultiVector& rho = *S->GetFieldData(dens_key_)
      ->ViewComponent("cell",false);

  const double& p_atm = *S->GetScalarData(patm_key_);
  const Epetra_Vector& gravity = *S->GetConstantVectorData(gravity_key_);
  double gz = -gravity[2];  // check this

  int ncells = res_c.MyLength();
  for (int c=0; c!=ncells; ++c) {
    res_c[0][c] = pres_c[0][c] < p_atm ? 0. : (pres_c[0][c] - p_atm) / (rho[0][c]*gz);
  }
}


void HeightEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  // this is rather hacky.  surface_pressure is a mixed field vector -- it has
  // pressure on cells and ponded depth on faces.
  // -- faces are 1
  result->ViewComponent("face",false)->PutScalar(1.0);

  // -- cells need the function eval
  const Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);
  const Epetra_MultiVector& pres_c = *S->GetFieldData(pres_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& rho = *S->GetFieldData(dens_key_)
      ->ViewComponent("cell",false);

  const double& p_atm = *S->GetScalarData(patm_key_);
  const Epetra_Vector& gravity = *S->GetConstantVectorData(gravity_key_);
  double gz = -gravity[2];  // check this

  if (wrt_key == pres_key_) {
    int ncells = res_c.MyLength();
    for (int c=0; c!=ncells; ++c) {
      //      res_c[0][c] = pres_c[0][c] < p_atm ? 0. : 1. / (rho[0][c]*gz);
      res_c[0][c] = 1. / (rho[0][c]*gz);
    }
  } else if (wrt_key == dens_key_) {
    int ncells = res_c.MyLength();
    for (int c=0; c!=ncells; ++c) {
      //      res_c[0][c] = pres_c[0][c] < p_atm ? 0. : -(pres_c[0][c] - p_atm)
      //          / (rho[0][c] * rho[0][c] * gz);
      res_c[0][c] = -(pres_c[0][c] - p_atm) / (rho[0][c] * rho[0][c] * gz);
    }
  } else {
    ASSERT(0);
  }
}



} //namespace
} //namespace
} //namespace
