/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow subgrid model.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "independent_variable_field_evaluator.hh"
#include "overland_conductivity_subgrid_evaluator.hh"
#include "manning_conductivity_model.hh"

namespace Amanzi {
namespace Flow {

OverlandConductivitySubgridEvaluator::OverlandConductivitySubgridEvaluator(Teuchos::ParameterList& plist)
    : SecondaryVariableFieldEvaluator(plist)
{
  Key domain = Keys::getDomain(my_key_);

  if (plist_.isParameter("height key") ||
      plist_.isParameter("ponded depth key") ||
      plist_.isParameter("height key suffix") ||
      plist_.isParameter("ponded depth key suffix")) {
    Errors::Message message("OverlandConductivitySubgrid: only use \"depth key\" or \"depth key suffix\", not \"height key\" or \"ponded depth key\".");
    Exceptions::amanzi_throw(message);
  }

  mobile_depth_key_ = Keys::readKey(plist_, domain, "mobile depth", "mobile_depth");
  dependencies_.insert(mobile_depth_key_);

  slope_key_ = Keys::readKey(plist_, domain, "slope", "slope_magnitude");
  dependencies_.insert(slope_key_);

  coef_key_ = Keys::readKey(plist_, domain, "coefficient", "manning_coefficient");
  dependencies_.insert(coef_key_);

  dens_key_ = Keys::readKey(plist_, domain, "molar density liquid", "molar_density_liquid");
  dependencies_.insert(dens_key_);

  frac_cond_key_ = Keys::readKey(plist_, domain, "fractional conductance", "fractional_conductance");
  dependencies_.insert(frac_cond_key_);

  drag_exp_key_ = Keys::readKey(plist_, domain, "drag exponent", "drag_exponent");
  dependencies_.insert(drag_exp_key_);

  // create the model
  Teuchos::ParameterList& sublist = plist_.sublist("overland conductivity model");
  model_ = Teuchos::rcp(new ManningConductivityModel(sublist));
}


Teuchos::RCP<FieldEvaluator>
OverlandConductivitySubgridEvaluator::Clone() const {
  return Teuchos::rcp(new OverlandConductivitySubgridEvaluator(*this));
}


// Required methods from SecondaryVariableFieldEvaluator
void OverlandConductivitySubgridEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  Teuchos::RCP<const CompositeVector> mobile_depth = S->GetFieldData(mobile_depth_key_);
  Teuchos::RCP<const CompositeVector> slope = S->GetFieldData(slope_key_);
  Teuchos::RCP<const CompositeVector> coef = S->GetFieldData(coef_key_);

  Teuchos::RCP<const CompositeVector> frac_cond = S->GetFieldData(frac_cond_key_);
  Teuchos::RCP<const CompositeVector> drag = S->GetFieldData(drag_exp_key_);

  Teuchos::RCP<const CompositeVector> dens = S->GetFieldData(dens_key_);

  for (const auto& comp : *result) {
    const Epetra_MultiVector& mobile_depth_v = *mobile_depth->ViewComponent(comp,false);
    const Epetra_MultiVector& slope_v = *slope->ViewComponent(comp,false);
    const Epetra_MultiVector& coef_v = *coef->ViewComponent(comp,false);

    const Epetra_MultiVector& frac_cond_v = *frac_cond->ViewComponent(comp,false);
    const Epetra_MultiVector& drag_v = *drag->ViewComponent(comp,false);
    const Epetra_MultiVector& dens_v = *S->GetFieldData(dens_key_)->ViewComponent(comp,false);
    Epetra_MultiVector& result_v = *result->ViewComponent(comp,false);

    int ncomp = result->size(comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = model_->Conductivity(mobile_depth_v[0][i], slope_v[0][i], coef_v[0][i]);
      result_v[0][i] *= dens_v[0][i] * std::pow(frac_cond_v[0][i], drag_v[0][i] + 1);
    }
  }
}


void OverlandConductivitySubgridEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  Teuchos::RCP<const CompositeVector> mobile_depth = S->GetFieldData(mobile_depth_key_);
  Teuchos::RCP<const CompositeVector> slope = S->GetFieldData(slope_key_);
  Teuchos::RCP<const CompositeVector> coef = S->GetFieldData(coef_key_);

  Teuchos::RCP<const CompositeVector> frac_cond = S->GetFieldData(frac_cond_key_);
  Teuchos::RCP<const CompositeVector> drag = S->GetFieldData(drag_exp_key_);

  Teuchos::RCP<const CompositeVector> dens = S->GetFieldData(dens_key_);

  if (wrt_key == mobile_depth_key_) {
    for (const auto& comp : *result) {
      const Epetra_MultiVector& mobile_depth_v = *mobile_depth->ViewComponent(comp,false);
      const Epetra_MultiVector& slope_v = *slope->ViewComponent(comp,false);
      const Epetra_MultiVector& coef_v = *coef->ViewComponent(comp,false);

      const Epetra_MultiVector& frac_cond_v = *frac_cond->ViewComponent(comp,false);
      const Epetra_MultiVector& drag_v = *drag->ViewComponent(comp,false);
      const Epetra_MultiVector& dens_v = *S->GetFieldData(dens_key_)->ViewComponent(comp,false);
      Epetra_MultiVector& result_v = *result->ViewComponent(comp,false);

      int ncomp = result->size(comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DConductivityDDepth(mobile_depth_v[0][i], slope_v[0][i], coef_v[0][i]);
        result_v[0][i] *= dens_v[0][i] * std::pow(frac_cond_v[0][i], drag_v[0][i] + 1);
      }
    }

  } else if (wrt_key == dens_key_) {
    for (const auto& comp : *result) {
      const Epetra_MultiVector& mobile_depth_v = *mobile_depth->ViewComponent(comp,false);
      const Epetra_MultiVector& slope_v = *slope->ViewComponent(comp,false);
      const Epetra_MultiVector& coef_v = *coef->ViewComponent(comp,false);

      const Epetra_MultiVector& frac_cond_v = *frac_cond->ViewComponent(comp,false);
      const Epetra_MultiVector& drag_v = *drag->ViewComponent(comp,false);
      const Epetra_MultiVector& dens_v = *S->GetFieldData(dens_key_)->ViewComponent(comp,false);
      Epetra_MultiVector& result_v = *result->ViewComponent(comp,false);

      int ncomp = result->size(comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->Conductivity(mobile_depth_v[0][i], slope_v[0][i], coef_v[0][i]);
        result_v[0][i] *= std::pow(frac_cond_v[0][i], drag_v[0][i] + 1);
      }
    }

  } else if (wrt_key == frac_cond_key_) {
    for (const auto& comp : *result) {
      const Epetra_MultiVector& mobile_depth_v = *mobile_depth->ViewComponent(comp,false);
      const Epetra_MultiVector& slope_v = *slope->ViewComponent(comp,false);
      const Epetra_MultiVector& coef_v = *coef->ViewComponent(comp,false);

      const Epetra_MultiVector& frac_cond_v = *frac_cond->ViewComponent(comp,false);
      const Epetra_MultiVector& drag_v = *drag->ViewComponent(comp,false);
      const Epetra_MultiVector& dens_v = *S->GetFieldData(dens_key_)->ViewComponent(comp,false);
      Epetra_MultiVector& result_v = *result->ViewComponent(comp,false);

      int ncomp = result->size(comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->Conductivity(mobile_depth_v[0][i], slope_v[0][i], coef_v[0][i]);
        result_v[0][i] *= dens_v[0][i] * (drag_v[0][i] + 1) *
                          std::pow(frac_cond_v[0][i], drag_v[0][i]);
      }
    }
  } else {
    result->PutScalar(0.);
  }
}


// void OverlandConductivitySubgridEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
//   // Ensure my field exists.  Requirements should be already set.
//   AMANZI_ASSERT(my_key_ != std::string(""));
//   Teuchos::RCP<CompositeVectorSpace> my_fac = S->RequireField(my_key_, my_key_);

//   // check plist for vis or checkpointing control
//   bool io_my_key = plist_.get<bool>("visualize", true);
//   S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
//   bool checkpoint_my_key = plist_.get<bool>("checkpoint", false);
//   S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);

//   // If my requirements have not yet been set, we'll have to hope they
//   // get set by someone later.  For now just defer.
//   if (my_fac->Mesh() != Teuchos::null) {
//     // Create an unowned factory to check my dependencies.
//     //
//     // Note only doing cells here, some other things would need to be fixed...
//     Teuchos::RCP<CompositeVectorSpace> dep_fac =
//         Teuchos::rcp(new CompositeVectorSpace());
//     dep_fac->SetMesh(my_fac->Mesh());
//     dep_fac->AddComponent("cell", AmanziMesh::CELL, 1);
//     dep_fac->SetGhosted(true);

//     // Loop over my dependencies, ensuring they meet the requirements.
//     for (const auto& key : dependencies_) {
//       if (key == my_key_) {
//         Errors::Message msg;
//         msg << "Evaluator for key \"" << my_key_ << "\" depends upon itself.";
//         Exceptions::amanzi_throw(msg);
//       }
//       Teuchos::RCP<CompositeVectorSpace> fac = S->RequireField(key);
//       fac->Update(*dep_fac);
//     }

//     // Recurse into the tree to propagate info to leaves.
//     for (const auto& key : dependencies_) {
//       S->RequireFieldEvaluator(key)->EnsureCompatibility(S);
//     }
//   }
// }



} // namespace Flow
} // namespace Amanzi

