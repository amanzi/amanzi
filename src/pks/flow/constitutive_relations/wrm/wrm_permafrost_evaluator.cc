/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  This WRM evaluator evaluates saturation of gas, liquid, and ice from
  capillary pressures for the ice-liquid and liquid-gas pairs.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "wrm_permafrost_evaluator.hh"
#include "wrm_factory.hh"
#include "wrm_permafrost_factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,WRMPermafrostEvaluator> WRMPermafrostEvaluator::factory_("permafrost WRM");


/* --------------------------------------------------------------------------------
  Constructor from just a ParameterList, reads WRMs and permafrost models from list.
 -------------------------------------------------------------------------------- */
WRMPermafrostEvaluator::WRMPermafrostEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariablesFieldEvaluator(plist) {

  // get the WRMs
  ASSERT(plist_.isSublist("WRM parameters"));
  WRMFactory wrm_fac;
  Teuchos::ParameterList region_list = plist_.sublist("WRM parameters");
  Teuchos::RCP<WRMRegionPairList> wrms = Teuchos::rcp(new WRMRegionPairList());

  for (Teuchos::ParameterList::ConstIterator lcv=region_list.begin();
       lcv!=region_list.end(); ++lcv) {
    std::string name = lcv->first;
    if (region_list.isSublist(name)) {
      Teuchos::ParameterList sublist = region_list.sublist(name);
      std::string region = sublist.get<std::string>("region");
      wrms->push_back(std::make_pair(region, wrm_fac.createWRM(sublist)));
    } else {
      ASSERT(0);
    }
  }

  // for each WRM create a permfrost_model
  WRMPermafrostFactory fac;
  permafrost_models_ = Teuchos::rcp(new WRMPermafrostModelRegionPairList());
  Teuchos::ParameterList pmodel_list = plist_.sublist("permafrost model parameters");

  for (WRMRegionPairList::const_iterator regionwrm=wrms->begin();
       regionwrm!=wrms->end(); ++regionwrm) {
    permafrost_models_->push_back(std::make_pair(regionwrm->first,
            fac.createWRMPermafrostModel(pmodel_list, regionwrm->second)));
  }

  InitializeFromPlist_();
}


/* --------------------------------------------------------------------------------
  Constructor with WRMs.
 -------------------------------------------------------------------------------- */
WRMPermafrostEvaluator::WRMPermafrostEvaluator(Teuchos::ParameterList& plist,
        const Teuchos::RCP<WRMRegionPairList>& wrms) :
    SecondaryVariablesFieldEvaluator(plist) {

  // for each WRM create a permfrost_model
  WRMPermafrostFactory fac;
  permafrost_models_ = Teuchos::rcp(new WRMPermafrostModelRegionPairList());
  Teuchos::ParameterList pmodel_list = plist_.sublist("permafrost model parameters");

  for (WRMRegionPairList::const_iterator regionwrm=wrms->begin();
       regionwrm!=wrms->end(); ++regionwrm) {
    permafrost_models_->push_back(std::make_pair(regionwrm->first,
            fac.createWRMPermafrostModel(pmodel_list, regionwrm->second)));
  }

  InitializeFromPlist_();
}


/* --------------------------------------------------------------------------------
  Constructor with Permafrost models.
 -------------------------------------------------------------------------------- */
WRMPermafrostEvaluator::WRMPermafrostEvaluator(Teuchos::ParameterList& plist,
        const Teuchos::RCP<WRMPermafrostModelRegionPairList>& models) :
    SecondaryVariablesFieldEvaluator(plist),
    permafrost_models_(models) {

  InitializeFromPlist_();
}


/* --------------------------------------------------------------------------------
  Copy constructor
 -------------------------------------------------------------------------------- */
WRMPermafrostEvaluator::WRMPermafrostEvaluator(const WRMPermafrostEvaluator& other) :
    SecondaryVariablesFieldEvaluator(other),
    pc_liq_key_(other.pc_liq_key_),
    pc_ice_key_(other.pc_ice_key_),
    s_l_key_(other.s_l_key_),
    permafrost_models_(other.permafrost_models_) {}


/* --------------------------------------------------------------------------------
  Virtual opy constructor as a FieldEvaluator.
 -------------------------------------------------------------------------------- */
Teuchos::RCP<FieldEvaluator>
WRMPermafrostEvaluator::Clone() const {
  return Teuchos::rcp(new WRMPermafrostEvaluator(*this));
}


/* --------------------------------------------------------------------------------
  Initialization of keys.
 -------------------------------------------------------------------------------- */
void WRMPermafrostEvaluator::InitializeFromPlist_() {
  // my keys are for saturation -- order matters... gas -> liq -> ice
  my_keys_.push_back(plist_.get<string>("gas saturation key", "saturation_gas"));
  s_l_key_ = plist_.get<string>("liquid saturation key", "saturation_liquid");
  my_keys_.push_back(s_l_key_);
  my_keys_.push_back(plist_.get<string>("ice saturation key", "saturation_ice"));
  setLinePrefix(my_keys_[0]+std::string(" evaluator"));

  // liquid-gas capillary pressure
  pc_liq_key_ = plist_.get<string>("liquid-gas capillary pressure key",
          "capillary_pressure_liq_gas");
  dependencies_.insert(pc_liq_key_);

  // liquid-gas capillary pressure
  pc_ice_key_ = plist_.get<string>("ice-liquid capillary pressure key",
          "capillary_pressure_ice_liq");
  dependencies_.insert(pc_ice_key_);

  // set up the verbose object
  setLinePrefix(std::string("saturation evaluator"));
}



void WRMPermafrostEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const std::vector<Teuchos::Ptr<CompositeVector> >& results) {
  Epetra_MultiVector& satg = *results[0]->ViewComponent("cell",false);
  Epetra_MultiVector& satl = *results[1]->ViewComponent("cell",false);
  Epetra_MultiVector& sati = *results[2]->ViewComponent("cell",false);

  const Epetra_MultiVector& pc_liq = *S->GetFieldData(pc_liq_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& pc_ice = *S->GetFieldData(pc_ice_key_)
      ->ViewComponent("cell",false);

  // Loop over names in the target and then owned entities in that name,
  // evaluating the evaluator to calculate sat.
  double sats[3];
  for (WRMPermafrostModelRegionPairList::iterator regionmodel=permafrost_models_->begin();
       regionmodel!=permafrost_models_->end(); ++regionmodel) {
    std::string region = regionmodel->first;
    int ncells = results[0]->mesh()->get_set_size(region,
            AmanziMesh::CELL, AmanziMesh::OWNED);
    AmanziMesh::Entity_ID_List cells(ncells);
    results[0]->mesh()->get_set_entities(region,
            AmanziMesh::CELL, AmanziMesh::OWNED, &cells);

    // use the model to evaluate saturations on each cell in the region
    for (AmanziMesh::Entity_ID_List::iterator c=cells.begin(); c!=cells.end(); ++c) {
      regionmodel->second->saturations(pc_liq[0][*c], pc_ice[0][*c], sats);
      satg[0][*c] = sats[0];
      satl[0][*c] = sats[1];
      sati[0][*c] = sats[2];
    }
  }
}


void WRMPermafrostEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results) {

  Epetra_MultiVector& dsatg = *results[0]->ViewComponent("cell",false);
  Epetra_MultiVector& dsatl = *results[1]->ViewComponent("cell",false);
  Epetra_MultiVector& dsati = *results[2]->ViewComponent("cell",false);

  const Epetra_MultiVector& pc_liq = *S->GetFieldData(pc_liq_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& pc_ice = *S->GetFieldData(pc_ice_key_)
      ->ViewComponent("cell",false);

  // Loop over names in the target and then owned entities in that name,
  // evaluating the evaluator to calculate sat and rel perm.
  double dsats[3];
  if (wrt_key == pc_liq_key_) {
    for (WRMPermafrostModelRegionPairList::iterator regionmodel=permafrost_models_->begin();
         regionmodel!=permafrost_models_->end(); ++regionmodel) {
      std::string region = regionmodel->first;
      int ncells = results[0]->mesh()->get_set_size(region,
              AmanziMesh::CELL, AmanziMesh::OWNED);
      AmanziMesh::Entity_ID_List cells(ncells);
      results[0]->mesh()->get_set_entities(region,
              AmanziMesh::CELL, AmanziMesh::OWNED, &cells);

      // use the model to evaluate saturation derivs on each cell in the region
      for (AmanziMesh::Entity_ID_List::iterator c=cells.begin(); c!=cells.end(); ++c) {
        regionmodel->second->dsaturations_dpc_liq(pc_liq[0][*c], pc_ice[0][*c], dsats);
        dsatg[0][*c] = dsats[0];
        dsatl[0][*c] = dsats[1];
        dsati[0][*c] = dsats[2];
      }
    }
  } else if (wrt_key == pc_ice_key_) {
    for (WRMPermafrostModelRegionPairList::iterator regionmodel=permafrost_models_->begin();
         regionmodel!=permafrost_models_->end(); ++regionmodel) {
      std::string region = regionmodel->first;
      int ncells = results[0]->mesh()->get_set_size(region,
              AmanziMesh::CELL, AmanziMesh::OWNED);
      AmanziMesh::Entity_ID_List cells(ncells);
      results[0]->mesh()->get_set_entities(region,
              AmanziMesh::CELL, AmanziMesh::OWNED, &cells);

      // use the model to evaluate saturation derivs on each cell in the region
      for (AmanziMesh::Entity_ID_List::iterator c=cells.begin(); c!=cells.end(); ++c) {
        regionmodel->second->dsaturations_dpc_ice(pc_liq[0][*c], pc_ice[0][*c], dsats);
        dsatg[0][*c] = dsats[0];
        dsatl[0][*c] = dsats[1];
        dsati[0][*c] = dsats[2];
      }
    }
  }
}



} // namespace
} // namespace
} // namespace



