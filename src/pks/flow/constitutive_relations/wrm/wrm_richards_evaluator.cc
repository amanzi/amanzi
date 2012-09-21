/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  This WRM evaluator calls saturation and rel perm using a capillary pressure p_atm - pc.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "wrm_richards_evaluator.hh"


namespace Amanzi {
namespace Flow {
namespace FlowRelations {

WRMRichardsEvaluator::WRMRichardsEvaluator(Teuchos::ParameterList& wrm_plist) :
    WRMEvaluator(wrm_plist) {
  InitializeFromPlist_();
}

WRMRichardsEvaluator::WRMRichardsEvaluator(Teuchos::ParameterList& wrm_plist,
        const Teuchos::RCP<WRMRegionPairList>& wrms) :
    WRMEvaluator(wrm_plist, wrms) {
  InitializeFromPlist_();
}

WRMRichardsEvaluator::WRMRichardsEvaluator(const WRMRichardsEvaluator& other) :
    WRMEvaluator(other),
    pres_key_(other.pres_key_) {}

Teuchos::RCP<FieldEvaluator>
WRMRichardsEvaluator::Clone() const {
  return Teuchos::rcp(new WRMRichardsEvaluator(*this));
}


void WRMRichardsEvaluator::InitializeFromPlist_() {
  // my keys are for saturation and rel perm.
  my_keys_.push_back(wrm_plist_.get<string>("saturation key", "saturation_liquid"));

  calc_other_sat_ = wrm_plist_.get<bool>("calculate minor saturation", true);
  if (calc_other_sat_) {
    my_keys_.push_back(wrm_plist_.get<string>("other saturation key", "saturation_gas"));
  }

  // my dependencies are just pressure.
  pres_key_ = wrm_plist_.get<string>("pressure key", "pressure");
  dependencies_.insert(pres_key_);
}


void WRMRichardsEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const std::vector<Teuchos::Ptr<CompositeVector> >& results) {
  Teuchos::Ptr<CompositeVector> sat = results[0];

  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pres_key_);
  Teuchos::RCP<const double> p_atm = S->GetScalarData("atmospheric_pressure");

  // Loop over names in the target and then owned entities in that name,
  // evaluating the evaluator to calculate sat.
  if (calc_other_sat_) {
    Teuchos::Ptr<CompositeVector> sat_g = results[1];
    for (WRMRegionPairList::iterator region=wrms_->begin();
         region!=wrms_->end(); ++region) {
      std::string name = region->first;
      int ncells = sat->mesh()->get_set_size(name,
              AmanziMesh::CELL, AmanziMesh::OWNED);
      std::vector<int> cells(ncells);
      sat->mesh()->get_set_entities(name,
              AmanziMesh::CELL, AmanziMesh::OWNED, &cells);

      // use the wrm to evaluate saturation on each cell in the region
      for (std::vector<int>::iterator c=cells.begin(); c!=cells.end(); ++c) {
        double pc = *p_atm - (*pres)("cell", *c);
        (*sat)("cell",*c) = region->second->saturation(pc);
        (*sat_g)("cell", *c) = 1.0 - (*sat)("cell", *c);
      }
    }
  } else {
    for (WRMRegionPairList::iterator region=wrms_->begin();
         region!=wrms_->end(); ++region) {
      std::string name = region->first;
      int ncells = sat->mesh()->get_set_size(name,
              AmanziMesh::CELL, AmanziMesh::OWNED);
      std::vector<int> cells(ncells);
      sat->mesh()->get_set_entities(name,
              AmanziMesh::CELL, AmanziMesh::OWNED, &cells);

      // use the wrm to evaluate saturation on each cell in the region
      for (std::vector<int>::iterator c=cells.begin(); c!=cells.end(); ++c) {
        double pc = *p_atm - (*pres)("cell", *c);
        (*sat)("cell",*c) = region->second->saturation(pc);
      }
    }
  }
}


void WRMRichardsEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results) {
  ASSERT(wrt_key == pres_key_);
  Teuchos::Ptr<CompositeVector> dsat = results[0];

  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pres_key_);
  Teuchos::RCP<const double> p_atm = S->GetScalarData("atmospheric_pressure");

  // Loop over names in the target and then owned entities in that name,
  // evaluating the evaluator to calculate sat and rel perm.
  if (calc_other_sat_) {
    Teuchos::Ptr<CompositeVector> dsat_g = results[1];
    for (WRMRegionPairList::iterator region=wrms_->begin();
         region!=wrms_->end(); ++region) {
      std::string name = region->first;
      int ncells = dsat->mesh()->get_set_size(name,
              AmanziMesh::CELL, AmanziMesh::OWNED);
      std::vector<int> cells(ncells);
      dsat->mesh()->get_set_entities(name,
              AmanziMesh::CELL, AmanziMesh::OWNED, &cells);

      // use the wrm to evaluate saturation on each cell in the region
      for (std::vector<int>::iterator c=cells.begin(); c!=cells.end(); ++c) {
        double pc = *p_atm - (*pres)("cell", *c);
        (*dsat)("cell",*c) = -region->second->d_saturation(pc);
        (*dsat_g)("cell", *c) = -(*dsat)("cell", *c);
      }
    }
  } else {
    for (WRMRegionPairList::iterator region=wrms_->begin();
         region!=wrms_->end(); ++region) {
      std::string name = region->first;
      int ncells = dsat->mesh()->get_set_size(name,
              AmanziMesh::CELL, AmanziMesh::OWNED);
      std::vector<int> cells(ncells);
      dsat->mesh()->get_set_entities(name,
              AmanziMesh::CELL, AmanziMesh::OWNED, &cells);

      // use the wrm to evaluate saturation on each cell in the region
      for (std::vector<int>::iterator c=cells.begin(); c!=cells.end(); ++c) {
        double pc = *p_atm - (*pres)("cell", *c);
        (*dsat)("cell",*c) = -region->second->d_saturation(pc);
      }
    }
  }
}


} // namespace
} // namespace
} // namespace



