/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Jeffrey Johnson (jnjohnson@lbl.gov)
*/

#include "TransportBoundaryFunction_Alquimia.hh"

#ifdef ALQUIMIA_ENABLED

namespace Amanzi {
namespace Transport {

/* ******************************************************************
* Constructor of BCs for Alquimia.
****************************************************************** */
TransportBoundaryFunction_Alquimia::TransportBoundaryFunction_Alquimia(
    const Teuchos::ParameterList& plist,
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
    Teuchos::RCP<AmanziChemistry::Alquimia_PK> chem_pk,
    Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine)
    : mesh_(mesh),
      chem_pk_(chem_pk),
      chem_engine_(chem_engine)
{
  // Check arguments.
  if (chem_engine_ != Teuchos::null) {
    chem_engine_->InitState(alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_);
    chem_engine_->GetPrimarySpeciesNames(tcc_names_);
  } else {
    Errors::Message msg;
    msg << "Geochemistry is off, but a geochemical condition was requested.";
    Exceptions::amanzi_throw(msg); 
  }

  // Get the regions assigned to this geochemical condition. If these regions have 
  // already been covered, we don't add a new condition.
  std::vector<std::string> regions = plist.get<Teuchos::Array<std::string> >("regions").toVector();
  cond_names_ = plist.get<Teuchos::Array<std::string> >("geochemical conditions").toVector();
  times_ = plist.get<Teuchos::Array<double> >("times").toVector();

  // Associate it with the given regions.
  Init_(regions);
}


/* ******************************************************************
* Delegating destructor.
****************************************************************** */
TransportBoundaryFunction_Alquimia::~TransportBoundaryFunction_Alquimia()
{
  chem_engine_->FreeState(alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_);
}


/* ******************************************************************
* Internal subroutine that defines a boundary function.
****************************************************************** */
void TransportBoundaryFunction_Alquimia::Init_(const std::vector<std::string> &regions)
{
  for (int i = 0; i < regions.size(); ++i) {
    // Get the faces that belong to this region (since boundary conditions
    // are applied on faces).
    assert(mesh_->valid_set_name(regions[i], AmanziMesh::FACE));
    unsigned int num_faces = mesh_->get_set_size(regions[i], AmanziMesh::FACE, AmanziMesh::OWNED);

    AmanziMesh::Entity_ID_List face_indices;
    mesh_->get_set_entities(regions[i], AmanziMesh::FACE, AmanziMesh::OWNED, &face_indices);

    // Now get the cells that are attached to these faces.
    for (int n = 0; n < num_faces; ++n) {
      int f = face_indices[n];
      value_[f] = WhetStone::DenseVector(chem_engine_->NumPrimarySpecies());

      AmanziMesh::Entity_ID_List cells_for_face;
      mesh_->face_get_cells(face_indices[f], AmanziMesh::OWNED, &cells_for_face);

      cell_for_face_[f] = cells_for_face[0];
    }
  }
}


/* ******************************************************************
* Evaluate values at time.
****************************************************************** */
void TransportBoundaryFunction_Alquimia::Compute(double t_old, double t_new) 
{
  // Find the condition that corresponds to the given time.
  int time_index = 0;
  while (time_index < (times_.size()-1)) {
    if (times_.at(time_index+1) > t_new)
      break;
    ++time_index;
  }
  std::string cond_name = cond_names_.at(time_index);

  // Loop over sides and evaluate values.
  for (TransportBoundaryFunction::Iterator it = begin(); it != end(); ++it) {
    int f = it->first; 
    // Find the index of the cell we're in.
    int cell = cell_for_face_[f];

    // Dump the contents of the chemistry state into our Alquimia containers.
    chem_pk_->CopyToAlquimia(cell, alq_mat_props_, alq_state_, alq_aux_data_);

    // Enforce the condition.
    chem_engine_->EnforceCondition(cond_name, t_new, alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_);

    // Move the concentrations into place.
    WhetStone::DenseVector& values = it->second;
    for (int i = 0; i < values.NumRows(); i++) {
      values(i) = alq_state_.total_mobile.data[i];
    }
  }
}


}  // namespace Transport
}  // namespace Amanzi

#endif

