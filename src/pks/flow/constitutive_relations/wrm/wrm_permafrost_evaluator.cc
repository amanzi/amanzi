/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  This WRM evaluator evaluates saturation of gas, liquid, and ice from
  capillary pressures for the ice-liquid and liquid-gas pairs.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "wrm_permafrost_evaluator.hh"
#include "wrm_partition.hh"

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
  Teuchos::ParameterList wrm_plist = plist_.sublist("WRM parameters");
  wrms_ = createWRMPartition(wrm_plist);

  // and the permafrost models
  ASSERT(plist_.isSublist("permafrost model parameters"));
  Teuchos::ParameterList perm_plist = plist_.sublist("permafrost model parameters");
  permafrost_models_ = createWRMPermafrostModelPartition(perm_plist, wrms_);

  InitializeFromPlist_();
}


/* --------------------------------------------------------------------------------
  Constructor with WRMs.
 -------------------------------------------------------------------------------- */
WRMPermafrostEvaluator::WRMPermafrostEvaluator(Teuchos::ParameterList& plist,
        const Teuchos::RCP<WRMPartition>& wrms) :
    SecondaryVariablesFieldEvaluator(plist),
    wrms_(wrms) {

  // and the permafrost models
  ASSERT(plist_.isSublist("permafrost model parameters"));
  Teuchos::ParameterList perm_plist = plist_.sublist("permafrost model parameters");
  permafrost_models_ = createWRMPermafrostModelPartition(perm_plist, wrms_);

  InitializeFromPlist_();
}


/* --------------------------------------------------------------------------------
  Constructor with Permafrost models.
 -------------------------------------------------------------------------------- */
WRMPermafrostEvaluator::WRMPermafrostEvaluator(Teuchos::ParameterList& plist,
        const Teuchos::RCP<WRMPermafrostModelPartition>& models) :
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
  pc_liq_key_ = plist_.get<string>("gas-liquid capillary pressure key",
          "capillary_pressure_gas_liq");
  dependencies_.insert(pc_liq_key_);

  // liquid-gas capillary pressure
  pc_ice_key_ = plist_.get<string>("liquid-ice capillary pressure key",
          "capillary_pressure_liq_ice");
  dependencies_.insert(pc_ice_key_);

  // set up the verbose object
  setLinePrefix(std::string("saturation evaluator"));
}



void WRMPermafrostEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const std::vector<Teuchos::Ptr<CompositeVector> >& results) {
  // Initialize the MeshPartition
  if (!permafrost_models_->first->initialized())
    permafrost_models_->first->Initialize(results[0]->mesh());

  // Cell values
  Epetra_MultiVector& satg_c = *results[0]->ViewComponent("cell",false);
  Epetra_MultiVector& satl_c = *results[1]->ViewComponent("cell",false);
  Epetra_MultiVector& sati_c = *results[2]->ViewComponent("cell",false);

  const Epetra_MultiVector& pc_liq_c = *S->GetFieldData(pc_liq_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& pc_ice_c = *S->GetFieldData(pc_ice_key_)
      ->ViewComponent("cell",false);

  double sats[3];
  int ncells = satg_c.MyLength();
  for (AmanziMesh::Entity_ID c=0; c!=ncells; ++c) {
    int i = (*permafrost_models_->first)[c];
    permafrost_models_->second[i]->saturations(pc_liq_c[0][c], pc_ice_c[0][c], sats);
    satg_c[0][c] = sats[0];
    satl_c[0][c] = sats[1];
    sati_c[0][c] = sats[2];
  }

  // Potentially do face values as well, though only for saturation_liquid?
  if (results[1]->has_component("boundary_face")) {
    ASSERT(!results[0]->has_component("boundary_face"));
    ASSERT(!results[2]->has_component("boundary_face"));

    Epetra_MultiVector& sat_bf = *results[1]->ViewComponent("boundary_face",false);
    const Epetra_MultiVector& pc_liq_bf = *S->GetFieldData(pc_liq_key_)
        ->ViewComponent("boundary_face",false);
    const Epetra_MultiVector& pc_ice_bf = *S->GetFieldData(pc_ice_key_)
        ->ViewComponent("boundary_face",false);

    // Need to get boundary face's inner cell to specify the WRM.
    Teuchos::RCP<const AmanziMesh::Mesh> mesh = results[0]->mesh();
    const Epetra_Map& vandelay_map = mesh->exterior_face_epetra_map();
    AmanziMesh::Entity_ID_List cells;

    // calculate boundary face values
    int nbfaces = sat_bf.MyLength();
    for (int bf=0; bf!=nbfaces; ++bf) {
      // given a boundary face, we need the internal cell to choose the right WRM
      AmanziMesh::Entity_ID f = vandelay_map.GID(bf);
      mesh->face_get_cells(f, AmanziMesh::USED, &cells);
      ASSERT(cells.size() == 1);

      int i = (*permafrost_models_->first)[cells[0]];
      permafrost_models_->second[i]
          ->saturations(pc_liq_bf[0][bf], pc_ice_bf[0][bf], sats);
      sat_bf[0][bf] = sats[1];
    }
  }
}


void
WRMPermafrostEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results) {

  // Cell values
  Epetra_MultiVector& satg_c = *results[0]->ViewComponent("cell",false);
  Epetra_MultiVector& satl_c = *results[1]->ViewComponent("cell",false);
  Epetra_MultiVector& sati_c = *results[2]->ViewComponent("cell",false);

  const Epetra_MultiVector& pc_liq_c = *S->GetFieldData(pc_liq_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& pc_ice_c = *S->GetFieldData(pc_ice_key_)
      ->ViewComponent("cell",false);

  double dsats[3];
  if (wrt_key == pc_liq_key_) {
    int ncells = satg_c.MyLength();
    for (AmanziMesh::Entity_ID c=0; c!=ncells; ++c) {
      int i = (*permafrost_models_->first)[c];
      permafrost_models_->second[i]->dsaturations_dpc_liq(
          pc_liq_c[0][c], pc_ice_c[0][c], dsats);

      satg_c[0][c] = dsats[0];
      satl_c[0][c] = dsats[1];
      sati_c[0][c] = dsats[2];
    }

  } else if (wrt_key == pc_ice_key_) {
    int ncells = satg_c.MyLength();
    for (AmanziMesh::Entity_ID c=0; c!=ncells; ++c) {
      int i = (*permafrost_models_->first)[c];
      permafrost_models_->second[i]->dsaturations_dpc_ice(
          pc_liq_c[0][c], pc_ice_c[0][c], dsats);

      satg_c[0][c] = dsats[0];
      satl_c[0][c] = dsats[1];
      sati_c[0][c] = dsats[2];
    }
  } else {
    ASSERT(0);
  }

  // Potentially do face values as well, though only for saturation_liquid?
  if (results[1]->has_component("boundary_face")) {
    ASSERT(!results[0]->has_component("boundary_face"));
    ASSERT(!results[2]->has_component("boundary_face"));

    Epetra_MultiVector& sat_bf = *results[1]->ViewComponent("boundary_face",false);
    const Epetra_MultiVector& pc_liq_bf = *S->GetFieldData(pc_liq_key_)
        ->ViewComponent("boundary_face",false);
    const Epetra_MultiVector& pc_ice_bf = *S->GetFieldData(pc_ice_key_)
        ->ViewComponent("boundary_face",false);

    // Need to get boundary face's inner cell to specify the WRM.
    Teuchos::RCP<const AmanziMesh::Mesh> mesh = results[0]->mesh();
    const Epetra_Map& vandelay_map = mesh->exterior_face_epetra_map();
    AmanziMesh::Entity_ID_List cells;

    if (wrt_key == pc_liq_key_) {
      // calculate boundary face values
      int nbfaces = sat_bf.MyLength();
      for (int bf=0; bf!=nbfaces; ++bf) {
        // given a boundary face, we need the internal cell to choose the right WRM
        AmanziMesh::Entity_ID f = vandelay_map.GID(bf);
        mesh->face_get_cells(f, AmanziMesh::USED, &cells);
        ASSERT(cells.size() == 1);

        int i = (*permafrost_models_->first)[cells[0]];
        permafrost_models_->second[i]->dsaturations_dpc_liq(
            pc_liq_bf[0][bf], pc_ice_bf[0][bf], dsats);
        sat_bf[0][bf] = dsats[1];
      }

    } else if (wrt_key == pc_ice_key_) {
      // calculate boundary face values
      int nbfaces = sat_bf.MyLength();
      for (int bf=0; bf!=nbfaces; ++bf) {
        // given a boundary face, we need the internal cell to choose the right WRM
        AmanziMesh::Entity_ID f = vandelay_map.GID(bf);
        mesh->face_get_cells(f, AmanziMesh::USED, &cells);
        ASSERT(cells.size() == 1);

        int i = (*permafrost_models_->first)[cells[0]];
        permafrost_models_->second[i]->dsaturations_dpc_ice(
            pc_liq_bf[0][bf], pc_ice_bf[0][bf], dsats);
        sat_bf[0][bf] = dsats[1];
      }
    } else {
      ASSERT(0);
    }
  }
}



} // namespace
} // namespace
} // namespace



