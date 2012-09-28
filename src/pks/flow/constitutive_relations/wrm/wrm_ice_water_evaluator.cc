/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  This WRM evaluator calls saturation and rel perm using a capillary pressure p_atm - pc.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "wrm_evaluator.hh"
#include "wrm_ice_water_evaluator.hh"


namespace Amanzi {
namespace Flow {
namespace FlowRelations {

WRMIceWaterEvaluator::WRMIceWaterEvaluator(Teuchos::ParameterList& plist) :
    WRMEvaluator(plist) {
  InitializeFromPlist_();
}

WRMIceWaterEvaluator::WRMIceWaterEvaluator(Teuchos::ParameterList& plist,
        const Teuchos::RCP<WRMRegionPairList>& wrms) :
    WRMEvaluator(plist, wrms) {
  InitializeFromPlist_();
}

WRMIceWaterEvaluator::WRMIceWaterEvaluator(const WRMIceWaterEvaluator& other) :
    WRMEvaluator(other),
    temp_key_(other.temp_key_),
    dens_key_(other.dens_key_),
    calc_other_sat_(other.calc_other_sat_),
    pc_(other.pc_) {}

Teuchos::RCP<FieldEvaluator>
WRMIceWaterEvaluator::Clone() const {
  return Teuchos::rcp(new WRMIceWaterEvaluator(*this));
}


void WRMIceWaterEvaluator::InitializeFromPlist_() {
  // my keys are for saturation and rel perm.
  my_keys_.push_back(plist_.get<string>("saturation key", "saturation_liquid"));
  setLinePrefix(my_keys_[0]+std::string(" evaluator"));

  calc_other_sat_ = plist_.get<bool>("calculate minor saturation", true);
  if (calc_other_sat_) {
    ASSERT(0);  // this is not implemented, as it physically doesn't make sense.
    //my_keys_.push_back(plist_.get<string>("other saturation key", "one_minus_sat_liq"));
  }

  // my dependencies are temperature and density
  temp_key_ = plist_.get<string>("temperature key", "temperature");
  dependencies_.insert(temp_key_);

  // get the physical evaluator for P_{c-il}(T)
  ASSERT(plist_.isSublist("capillary pressure of ice-water"));
  Teuchos::ParameterList pc_plist = plist_.sublist("capillary pressure of ice-water");
  pc_ = Teuchos::rcp(new PCIceWater(pc_plist));

  if (pc_->IsMolarBasis()) {
    dens_key_ = plist_.get<std::string>("molar density key", "molar_density_liquid");
  } else {
    dens_key_ = plist_.get<std::string>("mass density key", "mass_density_liquid");
  }
  dependencies_.insert(dens_key_);
}


void WRMIceWaterEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const std::vector<Teuchos::Ptr<CompositeVector> >& results) {
  Teuchos::Ptr<CompositeVector> sat = results[0];

  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  Teuchos::RCP<const CompositeVector> dens = S->GetFieldData(dens_key_);

  // Evaluate the model to calculate sat.
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
      double pc = pc_->CapillaryPressure((*temp)("cell", *c), (*dens)("cell", *c));
      (*sat)("cell",*c) = region->second->saturation(pc);
    }
  }
}


void WRMIceWaterEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results) {
  Teuchos::Ptr<CompositeVector> dsat = results[0];

  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  Teuchos::RCP<const CompositeVector> dens = S->GetFieldData(dens_key_);

  // Evaluate the model to calculate sat and rel perm.
  if (wrt_key == temp_key_) {
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
        double pc = pc_->CapillaryPressure((*temp)("cell", *c), (*dens)("cell", *c));
        (*dsat)("cell",*c) = region->second->d_saturation(pc)
            * pc_->DCapillaryPressureDT((*temp)("cell", *c), (*dens)("cell", *c));
      }
    }
  } else if (wrt_key == dens_key_) {
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
        double pc = pc_->CapillaryPressure((*temp)("cell", *c), (*dens)("cell", *c));
        (*dsat)("cell",*c) = region->second->d_saturation(pc)
            * pc_->DCapillaryPressureDRho((*temp)("cell", *c), (*dens)("cell", *c));
      }
    }
  } else {
    ASSERT(0);
  }
}

} // namespace
} // namespace
} // namespace



