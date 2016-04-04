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
    const std::vector<double>& times,
    const std::vector<std::string>& cond_names,
    const Teuchos::RCP<const AmanziMesh::Mesh> &mesh,
    Teuchos::RCP<AmanziChemistry::Alquimia_PK> chem_pk,
    Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine) :
    TransportBoundaryFunction(mesh),
    mesh_(mesh),
    times_(times),
    cond_names_(cond_names),
    chem_pk_(chem_pk),
    chem_engine_(chem_engine)
{
  // Check arguments.
  // NOTE: The times are always sorted in ascending order.
  if (times_.size() != cond_names_.size()) {
    Errors::Message msg;
    msg << "times and conditions arrays must be of equal size.";
    Exceptions::amanzi_throw(msg); 
  }
  if (times_.size() < 2) {
    Errors::Message msg;
    msg << "times and conditions arrays must contain at least two elements.";
    Exceptions::amanzi_throw(msg); 
  }

  if (chem_engine_ != Teuchos::null) {
    chem_engine_->InitState(alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_);
    chem_engine_->GetPrimarySpeciesNames(tcc_names_);
  } else {
    Errors::Message msg;
    msg << "Geochemistry is off, but a geochemical condition was requested.";
    Exceptions::amanzi_throw(msg); 
  }
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
void TransportBoundaryFunction_Alquimia::Define(const std::vector<std::string> &regions)
{
  for (size_t i = 0; i < regions.size(); ++i) {
    // Get the faces that belong to this region (since boundary conditions
    // are applied on faces).
    assert(mesh_->valid_set_name(regions[i], AmanziMesh::FACE));
    unsigned int num_faces = mesh_->get_set_size(regions[i], AmanziMesh::FACE, AmanziMesh::OWNED);

    AmanziMesh::Entity_ID_List face_indices;
    std::vector<double> vofs;
    mesh_->get_set_entities(regions[i], AmanziMesh::FACE, AmanziMesh::OWNED, &face_indices, &vofs);

    // Now get the cells that are attached to these faces.
    faces_.resize(face_indices.size());
    values_.resize(face_indices.size());
    for (unsigned int f = 0; f < num_faces; ++f) {
      faces_[f] = face_indices[f];
      values_[f].resize(chem_engine_->NumPrimarySpecies());
      AmanziMesh::Entity_ID_List cells_for_face;
      mesh_->face_get_cells(face_indices[f], AmanziMesh::OWNED, &cells_for_face);
      assert(cells_for_face.size() == 1); // Only one cell per boundary face, right?
      cell_for_face_[face_indices[f]] = cells_for_face[0];
    }
  }
}


/* ******************************************************************
* Internal subroutine that defines a boundary function.
****************************************************************** */
void TransportBoundaryFunction_Alquimia::Define(std::string region) 
{
  std::vector<std::string> regions(1, region);
  Define(regions);
}


/* ******************************************************************
* Evaluate values at time.
****************************************************************** */
void TransportBoundaryFunction_Alquimia::Compute(double time) 
{
  // Find the condition that corresponds to the given time.
  int time_index = 0;
  while (time_index < (times_.size()-1)) {
    if (times_.at(time_index+1) > time)
      break;
    ++time_index;
  }
  std::string cond_name = cond_names_.at(time_index);

  // Loop over sides and evaluate values.
  for (int n = 0; n < faces_.size(); n++) {
    // Find the index of the cell we're in.
    int cell = cell_for_face_[faces_[n]];

    // Dump the contents of the chemistry state into our Alquimia containers.
    chem_pk_->CopyToAlquimia(cell, alq_mat_props_, alq_state_, alq_aux_data_);

    // Enforce the condition.
    chem_engine_->EnforceCondition(cond_name, time, alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_);

    // Move the concentrations into place.
    for (int i = 0; i < values_[n].size(); i++) {
      values_[n][i] = alq_state_.total_mobile.data[i];
    }
  }
}

}  // namespace Transport
}  // namespace Amanzi

#endif

