/*
  This is the transport component of the Amanzi code. 

  License: see $AMANZI_DIR/COPYRIGHT
  Author (v1): Neil Carlson
         (v2): Ethan Coon
*/

#include "TransportBoundaryFunction_Alquimia.hh"

#ifdef ALQUIMIA_ENABLED

namespace Amanzi {
namespace AmanziTransport {

TransportBoundaryFunction_Alquimia::TransportBoundaryFunction_Alquimia(const std::string& cond_name,
                                                                       const Teuchos::RCP<const AmanziMesh::Mesh> &mesh,
                                                                       Teuchos::RCP<AmanziChemistry::Chemistry_State> chem_state,
                                                                       Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine):
  TransportBoundaryFunction(mesh), mesh_(mesh), cond_name_(cond_name), chem_state_(chem_state), chem_engine_(chem_engine)
{
  if (chem_engine_ != Teuchos::null)
  {
    chem_engine_->InitState(alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_);
    chem_engine_->GetPrimarySpeciesNames(tcc_names_);
  }
  else
  {
    Errors::Message msg;
    msg << "Geochemistry is off, but a geochemical condition was requested.";
    Exceptions::amanzi_throw(msg); 
  }
}

TransportBoundaryFunction_Alquimia::~TransportBoundaryFunction_Alquimia()
{
  chem_engine_->FreeState(alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_);
}

/* ******************************************************************
* Internal subroutine that defines a boundary function.
****************************************************************** */
void TransportBoundaryFunction_Alquimia::Define(const std::vector<std::string> &regions)
{
  for (size_t i = 0; i < regions.size(); ++i)
  {
    // Get the faces that belong to this region (since boundary conditions
    // are applied on faces).
    assert(mesh_->valid_set_name(regions[i], AmanziMesh::FACE));
    unsigned int num_faces = mesh_->get_set_size(regions[i], AmanziMesh::FACE, AmanziMesh::OWNED);
    AmanziMesh::Entity_ID_List face_indices;
    mesh_->get_set_entities(regions[i], AmanziMesh::FACE, AmanziMesh::OWNED, &face_indices);

    // Now get the cells that are attached to these faces.
    faces_.resize(face_indices.size());
    values_.resize(face_indices.size());
    for (unsigned int f = 0; f < num_faces; ++f)
    {
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
  // Loop over sides and evaluate values.
  for (int n = 0; n < faces_.size(); n++) {
    // Find the index of the cell we're in.
    int cell = cell_for_face_[faces_[n]];

    // Dump the contents of the chemistry state into our Alquimia containers.
    chem_state_->CopyToAlquimia(cell, alq_mat_props_, alq_state_, alq_aux_data_);

    // Enforce the condition.
    chem_engine_->EnforceCondition(cond_name_, time, alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_);

    // Move the concentrations into place.
    for (int i = 0; i < values_[n].size(); i++) {
      values_[n][i] = alq_state_.total_mobile.data[i];
std::cout << "[" << tcc_names_[i] << "] -> " << values_[n][i] << std::endl;
    }
  }
}

}  // namespace AmanziTransport
}  // namespace Amanzi

#endif

