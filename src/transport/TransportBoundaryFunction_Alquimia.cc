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
  TransportBoundaryFunction(mesh), cond_name_(cond_name), chem_state_(chem_state), chem_engine_(chem_engine)
{
  chem_engine_->InitState(alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_);
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
  // FIXME
}


/* ******************************************************************
* Internal subroutine that defines a boundary function.
****************************************************************** */
void TransportBoundaryFunction_Alquimia::Define(std::string region) 
{
  RegionList regions(1,region);
  // FIXME
}


/* ******************************************************************
* Evaluate values at time.
****************************************************************** */
void TransportBoundaryFunction_Alquimia::Compute(double time) 
{
  // Loop over sides and evaluate values.
  for (int n = 0; n < faces_.size(); n++) {
    // FIXME: Find out the index of the cell we're in.
    int cell;

    // Dump the contents of the chemistry state into our Alquimia containers.
    chem_state_->CopyToAlquimia(cell, alq_mat_props_, alq_state_, alq_aux_data_);

    // Enforce the condition.
    chem_engine_->EnforceCondition(cond_name_, time, alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_);

    // Move the concentrations into place.
    for (int i = 0; i < tcc_index_.size(); i++) {
      values_[n][i] = alq_state_.total_mobile.data[i];
    }
  }
}

}  // namespace AmanziTransport
}  // namespace Amanzi

#endif

