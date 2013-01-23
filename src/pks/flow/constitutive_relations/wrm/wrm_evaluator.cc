/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  The WRM Evaluator simply calls the WRM with the correct arguments.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "wrm_evaluator.hh"
#include "wrm_factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,WRMEvaluator> WRMEvaluator::factory_("WRM");

WRMEvaluator::WRMEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariablesFieldEvaluator(plist),
    calc_other_sat_(true) {
  ASSERT(plist_.isSublist("WRM parameters"));

  WRMFactory fac;

  Teuchos::ParameterList region_list = plist_.sublist("WRM parameters");
  wrms_ = Teuchos::rcp(new WRMRegionPairList());

  for (Teuchos::ParameterList::ConstIterator lcv=region_list.begin();
       lcv!=region_list.end(); ++lcv) {
    std::string name = lcv->first;
    if (region_list.isSublist(name)) {
      Teuchos::ParameterList sublist = region_list.sublist(name);
      std::string region = sublist.get<std::string>("region");
      wrms_->push_back(std::make_pair(region, fac.createWRM(sublist)));
    } else {
      ASSERT(0);
    }
  }

  InitializeFromPlist_();
}

WRMEvaluator::WRMEvaluator(Teuchos::ParameterList& plist,
                           const Teuchos::RCP<std::vector<WRMRegionPair> >& wrms) :
    SecondaryVariablesFieldEvaluator(plist),
    wrms_(wrms) {
  InitializeFromPlist_();
}

WRMEvaluator::WRMEvaluator(const WRMEvaluator& other) :
    SecondaryVariablesFieldEvaluator(other),
    calc_other_sat_(other.calc_other_sat_),
    cap_pres_key_(other.cap_pres_key_),
    wrms_(other.wrms_) {}


Teuchos::RCP<FieldEvaluator> WRMEvaluator::Clone() const {
  return Teuchos::rcp(new WRMEvaluator(*this));
}

void WRMEvaluator::InitializeFromPlist_() {
  // my keys are for saturation
  my_keys_.push_back(plist_.get<string>("saturation key", "saturation_liquid"));

  calc_other_sat_ = plist_.get<bool>("calculate minor saturation", true);
  if (calc_other_sat_) {
    my_keys_.push_back(plist_.get<string>("other saturation key", "saturation_gas"));
  }

  // my dependencies are capillary pressure.
  cap_pres_key_ = plist_.get<string>("capillary pressure key", "capillary_pressure_gas_liq");
  dependencies_.insert(cap_pres_key_);

  // set up the verbose object
  setLinePrefix(std::string("saturation evaluator"));
}


void WRMEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const std::vector<Teuchos::Ptr<CompositeVector> >& results) {
  Epetra_MultiVector& sat = *results[0]->ViewComponent("cell",false);
  const Epetra_MultiVector& pres = *S->GetFieldData(cap_pres_key_)->ViewComponent("cell",false);

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
        double s = region->second->saturation(pres[0][*c]);
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
        sat[0][*c] = region->second->saturation(pres[0][*c]);
      }
    }
  }
}


void WRMEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results) {
  ASSERT(wrt_key == cap_pres_key_);
  Epetra_MultiVector& dsat = *results[0]->ViewComponent("cell",false);
  const Epetra_MultiVector& pres = *S->GetFieldData(cap_pres_key_)->ViewComponent("cell",false);
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
        double ds = region->second->d_saturation(pres[0][*c]);
        dsat[0][*c] = ds;
        dsat_g[0][*c] = -ds;
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
        dsat[0][*c] = region->second->d_saturation(pres[0][*c]);
      }
    }
  }
}


} //namespace
} //namespace
} //namespace
