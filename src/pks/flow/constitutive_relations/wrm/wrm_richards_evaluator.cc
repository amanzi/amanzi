/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  This WRM evaluator calls saturation and rel perm using a capillary pressure p_atm - pc.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "wrm_richards_evaluator.hh"


namespace Amanzi {
namespace Flow {
namespace FlowRelations {

WRMRichardsEvaluator::WRMRichardsEvaluator(Teuchos::ParameterList& plist) :
    WRMEvaluator(plist),
    calc_other_sat_(false) {
  InitializeFromPlist_();
}

WRMRichardsEvaluator::WRMRichardsEvaluator(Teuchos::ParameterList& plist,
        const Teuchos::RCP<WRMRegionPairList>& wrms) :
    WRMEvaluator(plist, wrms),
    calc_other_sat_(false) {
  InitializeFromPlist_();
}

WRMRichardsEvaluator::WRMRichardsEvaluator(const WRMRichardsEvaluator& other) :
    WRMEvaluator(other),
    calc_other_sat_(other.calc_other_sat_),
    pres_key_(other.pres_key_) {}

Teuchos::RCP<FieldEvaluator>
WRMRichardsEvaluator::Clone() const {
  return Teuchos::rcp(new WRMRichardsEvaluator(*this));
}


void WRMRichardsEvaluator::InitializeFromPlist_() {
  // my keys are for saturation and rel perm.
  my_keys_.push_back(plist_.get<string>("saturation key", "saturation_liquid"));

  calc_other_sat_ = plist_.get<bool>("calculate minor saturation", true);
  if (calc_other_sat_) {
    my_keys_.push_back(plist_.get<string>("other saturation key", "saturation_gas"));
  }

  // my dependencies are just pressure.
  pres_key_ = plist_.get<string>("pressure key", "pressure");
  dependencies_.insert(pres_key_);

  // set up the verbose object
  setLinePrefix(my_keys_[0]+std::string(" evaluator"));
}


void WRMRichardsEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const std::vector<Teuchos::Ptr<CompositeVector> >& results) {
  Epetra_MultiVector& sat = *results[0]->ViewComponent("cell",false);
  const Epetra_MultiVector& pres = *S->GetFieldData(pres_key_)->ViewComponent("cell",false);
  const double& p_atm = *(S->GetScalarData("atmospheric_pressure"));

  // Loop over names in the target and then owned entities in that name,
  // evaluating the evaluator to calculate sat.
  if (calc_other_sat_) {
    Epetra_MultiVector& sat_g = *results[1]->ViewComponent("cell",false);

    for (WRMRegionPairList::iterator region=wrms_->begin();
         region!=wrms_->end(); ++region) {
      std::string name = region->first;
      int ncells = results[0]->mesh()->get_set_size(name,
              AmanziMesh::CELL, AmanziMesh::OWNED);
      AmanziMesh::Entity_ID_List cells(ncells);
      results[0]->mesh()->get_set_entities(name,
              AmanziMesh::CELL, AmanziMesh::OWNED, &cells);

      // use the wrm to evaluate saturation on each cell in the region
      for (AmanziMesh::Entity_ID_List::iterator c=cells.begin(); c!=cells.end(); ++c) {
        double s = region->second->saturation(p_atm - pres[0][*c]);
        sat[0][*c] = s;
        sat_g[0][*c] = 1.0 - s;
      }
    }
  } else {
    for (WRMRegionPairList::iterator region=wrms_->begin();
         region!=wrms_->end(); ++region) {
      std::string name = region->first;
      int ncells = results[0]->mesh()->get_set_size(name,
              AmanziMesh::CELL, AmanziMesh::OWNED);
      AmanziMesh::Entity_ID_List cells(ncells);
      results[0]->mesh()->get_set_entities(name,
              AmanziMesh::CELL, AmanziMesh::OWNED, &cells);

      // use the wrm to evaluate saturation on each cell in the region
      for (AmanziMesh::Entity_ID_List::iterator c=cells.begin(); c!=cells.end(); ++c) {
        sat[0][*c] = region->second->saturation(p_atm - pres[0][*c]);
      }
    }
  }
}


void WRMRichardsEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results) {
  ASSERT(wrt_key == pres_key_);
  Epetra_MultiVector& dsat = *results[0]->ViewComponent("cell",false);
  const Epetra_MultiVector& pres = *S->GetFieldData(pres_key_)->ViewComponent("cell",false);
  const double& p_atm = *(S->GetScalarData("atmospheric_pressure"));

  // Loop over names in the target and then owned entities in that name,
  // evaluating the evaluator to calculate sat and rel perm.
  if (calc_other_sat_) {
    Epetra_MultiVector& dsat_g = *results[1]->ViewComponent("cell",false);

    for (WRMRegionPairList::iterator region=wrms_->begin();
         region!=wrms_->end(); ++region) {
      std::string name = region->first;
      int ncells = results[0]->mesh()->get_set_size(name,
              AmanziMesh::CELL, AmanziMesh::OWNED);
      AmanziMesh::Entity_ID_List cells(ncells);
      results[0]->mesh()->get_set_entities(name,
              AmanziMesh::CELL, AmanziMesh::OWNED, &cells);

      // use the wrm to evaluate saturation on each cell in the region
      for (AmanziMesh::Entity_ID_List::iterator c=cells.begin(); c!=cells.end(); ++c) {
        double ds = region->second->d_saturation(p_atm - pres[0][*c]);
        dsat[0][*c] = -ds;
        dsat_g[0][*c] = ds;
      }
    }
  } else {
    for (WRMRegionPairList::iterator region=wrms_->begin();
         region!=wrms_->end(); ++region) {
      std::string name = region->first;
      int ncells = results[0]->mesh()->get_set_size(name,
              AmanziMesh::CELL, AmanziMesh::OWNED);
      AmanziMesh::Entity_ID_List cells(ncells);
      results[0]->mesh()->get_set_entities(name,
              AmanziMesh::CELL, AmanziMesh::OWNED, &cells);

      // use the wrm to evaluate saturation on each cell in the region
      for (AmanziMesh::Entity_ID_List::iterator c=cells.begin(); c!=cells.end(); ++c) {
        double ds = region->second->d_saturation(p_atm - pres[0][*c]);
        dsat[0][*c] = -ds;
      }
    }
  }
}


} // namespace
} // namespace
} // namespace



