/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Rel perm( pc ( sat ) ).

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "rel_perm_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

RelPermEvaluator::RelPermEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist),
    min_val_(0.) {

  ASSERT(plist_.isSublist("WRM parameters"));
  Teuchos::ParameterList sublist = plist_.sublist("WRM parameters");
  wrms_ = createWRMPartition(sublist);
  InitializeFromPlist_();
}

RelPermEvaluator::RelPermEvaluator(Teuchos::ParameterList& plist,
        const Teuchos::RCP<WRMPartition>& wrms) :
    SecondaryVariableFieldEvaluator(plist),
    wrms_(wrms),
    min_val_(0.) {
  InitializeFromPlist_();
}

RelPermEvaluator::RelPermEvaluator(const RelPermEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    wrms_(other.wrms_),
    sat_key_(other.sat_key_),
    dens_key_(other.dens_key_),
    visc_key_(other.visc_key_),
    perm_scale_(other.perm_scale_),
    min_val_(other.min_val_) {}


Teuchos::RCP<FieldEvaluator>
RelPermEvaluator::Clone() const {
  return Teuchos::rcp(new RelPermEvaluator(*this));
}


void RelPermEvaluator::InitializeFromPlist_() {
  // my keys are for saturation and rel perm.
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<string>("rel perm key", "relative_permeability");
  }

  // dependencies
  sat_key_ = plist_.get<string>("saturation key", "saturation_liquid");
  dependencies_.insert(sat_key_);

  dens_key_ = plist_.get<string>("density key", "molar_density_liquid");
  dependencies_.insert(dens_key_);

  visc_key_ = plist_.get<string>("viscosity key", "viscosity_liquid");
  dependencies_.insert(visc_key_);

  // cutoff above 0?
  min_val_ = plist_.get<double>("minimum rel perm cutoff", 0.);
  perm_scale_ = plist_.get<double>("permeability rescaling", 1.);
}


void RelPermEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  // initialize the MeshPartition if needed to make sure all cells are covered
  if (!wrms_->first->initialized()) wrms_->first->Initialize(result->Mesh());

  const Epetra_MultiVector& sat_c = *S->GetFieldData(sat_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& dens_c = *S->GetFieldData(dens_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& visc_c = *S->GetFieldData(visc_key_)
      ->ViewComponent("cell",false);

  Epetra_MultiVector& res_v = *result->ViewComponent("cell",false);

  // Evaluate the evaluator to calculate sat.
  int ncells = res_v.MyLength();
  for (int c=0; c!=ncells; ++c) {
    int index = (*wrms_->first)[c];
    double pc = wrms_->second[index]->capillaryPressure(sat_c[0][c]);
    res_v[0][c] = dens_c[0][c] * std::max(wrms_->second[index]->k_relative(pc), min_val_)
      / visc_c[0][c];
  }

  // Potentially do face values as well.
  if (result->HasComponent("boundary_face")) {
    const Epetra_MultiVector& sat_bf = *S->GetFieldData(sat_key_)
        ->ViewComponent("boundary_face",false);
    const Epetra_MultiVector& dens_bf = *S->GetFieldData(dens_key_)
        ->ViewComponent("boundary_face",false);
    const Epetra_MultiVector& visc_bf = *S->GetFieldData(visc_key_)
        ->ViewComponent("boundary_face",false);
    Epetra_MultiVector& res_v = *result->ViewComponent("boundary_face",false);

    Teuchos::RCP<const AmanziMesh::Mesh> mesh = result->Mesh();
    const Epetra_Map& vandelay_map = mesh->exterior_face_epetra_map();
    const Epetra_Map& face_map = mesh->face_epetra_map(false);

    // Evaluate the evaluator to calculate sat.
    AmanziMesh::Entity_ID_List cells;
    int nbfaces = res_v.MyLength();
    for (int bf=0; bf!=nbfaces; ++bf) {
      // given a boundary face, we need the internal cell to choose the right WRM
      AmanziMesh::Entity_ID f = face_map.LID(vandelay_map.GID(bf));
      mesh->face_get_cells(f, AmanziMesh::USED, &cells);
      ASSERT(cells.size() == 1);

      int index = (*wrms_->first)[cells[0]];
      double pc = wrms_->second[index]->capillaryPressure(sat_bf[0][bf]);
      res_v[0][bf] = dens_bf[0][bf] * std::max(wrms_->second[index]->k_relative(pc), min_val_)
        / visc_bf[0][bf];
    }
  }

  result->Scale(1./perm_scale_);
}


void RelPermEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  // initialize the MeshPartition if needed to make sure all cells are covered
  if (!wrms_->first->initialized()) wrms_->first->Initialize(result->Mesh());
  const Epetra_MultiVector& sat_c = *S->GetFieldData(sat_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& dens_c = *S->GetFieldData(dens_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& visc_c = *S->GetFieldData(visc_key_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& res_v = *result->ViewComponent("cell",false);


  if (wrt_key == sat_key_) {
    // Evaluate the evaluator to calculate sat.
    int ncells = res_v.MyLength();
    for (int c=0; c!=ncells; ++c) {
      int index = (*wrms_->first)[c];
      double pc = wrms_->second[index]->capillaryPressure(sat_c[0][c]);
      if (pc <= 0.) {
        res_v[0][c] = 0.;
      } else {
        double dpc_ds = wrms_->second[index]->d_capillaryPressure(sat_c[0][c]);
        res_v[0][c] = dens_c[0][c] * wrms_->second[index]->d_k_relative(pc) * dpc_ds
          / visc_c[0][c];
      }
    }

  } else if (wrt_key == dens_key_) {
    // Evaluate the evaluator to calculate sat.
    int ncells = res_v.MyLength();
    for (int c=0; c!=ncells; ++c) {
      int index = (*wrms_->first)[c];
      double pc = wrms_->second[index]->capillaryPressure(sat_c[0][c]);
      res_v[0][c] = std::max(wrms_->second[index]->k_relative(pc), min_val_) / visc_c[0][c];
    }

  } else if (wrt_key == visc_key_) {
    // Evaluate the evaluator to calculate sat.
    int ncells = res_v.MyLength();
    for (int c=0; c!=ncells; ++c) {
      int index = (*wrms_->first)[c];
      double pc = wrms_->second[index]->capillaryPressure(sat_c[0][c]);
      res_v[0][c] = - dens_c[0][c] * std::max(wrms_->second[index]->k_relative(pc), min_val_) 
        / (std::pow(visc_c[0][c],2));
    }
  } else {
    ASSERT(0);
  }

  // Potentially do face values as well.
  if (result->HasComponent("boundary_face")) {
    const Epetra_MultiVector& sat_bf = *S->GetFieldData(sat_key_)
      ->ViewComponent("boundary_face",false);
    const Epetra_MultiVector& dens_bf = *S->GetFieldData(dens_key_)
      ->ViewComponent("boundary_face",false);
    const Epetra_MultiVector& visc_bf = *S->GetFieldData(visc_key_)
      ->ViewComponent("boundary_face",false);
    Epetra_MultiVector& res_v = *result->ViewComponent("boundary_face",false);

    Teuchos::RCP<const AmanziMesh::Mesh> mesh = result->Mesh();
    const Epetra_Map& vandelay_map = mesh->exterior_face_epetra_map();
    const Epetra_Map& face_map = mesh->face_epetra_map(false);

    if (wrt_key == sat_key_) {
      // Evaluate the evaluator to calculate sat.
      AmanziMesh::Entity_ID_List cells;
      int nbfaces = res_v.MyLength();
      for (int bf=0; bf!=nbfaces; ++bf) {
        // given a boundary face, we need the internal cell to choose the right WRM
        AmanziMesh::Entity_ID f = face_map.LID(vandelay_map.GID(bf));
        mesh->face_get_cells(f, AmanziMesh::USED, &cells);
        ASSERT(cells.size() == 1);

        int index = (*wrms_->first)[cells[0]];
        double pc = wrms_->second[index]->capillaryPressure(sat_bf[0][bf]);
        if (pc <= 0.) {
          res_v[0][bf] = 0.;
        } else {
          double dpc_ds = wrms_->second[index]->d_capillaryPressure(sat_bf[0][bf]);
          res_v[0][bf] = dens_bf[0][bf] * wrms_->second[index]->d_k_relative(pc) * dpc_ds
            / visc_bf[0][bf];
        }
      }
    } else if (wrt_key == dens_key_) {
      // Evaluate the evaluator to calculate sat.
      // Evaluate the evaluator to calculate sat.
      AmanziMesh::Entity_ID_List cells;
      int nbfaces = res_v.MyLength();
      for (int bf=0; bf!=nbfaces; ++bf) {
        // given a boundary face, we need the internal cell to choose the right WRM
        AmanziMesh::Entity_ID f = face_map.LID(vandelay_map.GID(bf));
        mesh->face_get_cells(f, AmanziMesh::USED, &cells);
        ASSERT(cells.size() == 1);

        int index = (*wrms_->first)[cells[0]];
        double pc = wrms_->second[index]->capillaryPressure(sat_bf[0][bf]);
        res_v[0][bf] = std::max(wrms_->second[index]->k_relative(pc), min_val_) / visc_bf[0][bf];
      }

    } else if (wrt_key == visc_key_) {
      // Evaluate the evaluator to calculate sat.
      // Evaluate the evaluator to calculate sat.
      AmanziMesh::Entity_ID_List cells;
      int nbfaces = res_v.MyLength();
      for (int bf=0; bf!=nbfaces; ++bf) {
        // given a boundary face, we need the internal cell to choose the right WRM
        AmanziMesh::Entity_ID f = face_map.LID(vandelay_map.GID(bf));
        mesh->face_get_cells(f, AmanziMesh::USED, &cells);
        ASSERT(cells.size() == 1);

        int index = (*wrms_->first)[cells[0]];
        double pc = wrms_->second[index]->capillaryPressure(sat_bf[0][bf]);
        res_v[0][bf] = - dens_bf[0][bf] * std::max(wrms_->second[index]->k_relative(pc), min_val_) 
          / (std::pow(visc_bf[0][bf],2));

      }
    } else {
      ASSERT(0);
    }
  }

  result->Scale(1./perm_scale_);
}



} //namespace
} //namespace
} //namespace
