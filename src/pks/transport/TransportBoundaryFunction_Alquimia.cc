/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Jeffrey Johnson (jnjohnson@lbl.gov)
*/

/*
  Transport PK

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
  Teuchos::RCP<AmanziChemistry::Alquimia_PK> alquimia_pk,
  Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine)
  : name_("alquimia bc"), mesh_(mesh), alquimia_pk_(alquimia_pk), chem_engine_(chem_engine)
{
  // Check arguments.
  if (chem_engine_ != Teuchos::null) {
    chem_engine_->InitState(
      beaker_.properties, beaker_.state, beaker_.aux_data, beaker_.aux_output);
    chem_engine_->GetPrimarySpeciesNames(tcc_names_);
  } else {
    Errors::Message msg;
    msg << "Geochemistry is off, but a geochemical condition was requested.";
    Exceptions::amanzi_throw(msg);
  }

  // Get the regions assigned to this geochemical condition. We do not
  // check for region overlaps here, since a better way is to derive from
  // the generic mesh function.
  std::vector<std::string> regions = plist.get<Teuchos::Array<std::string>>("regions").toVector();
  std::vector<double> times = plist.get<Teuchos::Array<double>>("times").toVector();
  std::vector<std::string> conditions =
    plist.get<Teuchos::Array<std::string>>("geochemical conditions").toVector();

  // Function of geochemical conditions and the associates regions.
  f_ = Teuchos::rcp(new FunctionTabularString(times, conditions));
  Init_(regions);
}


/* ******************************************************************
* Delegating destructor.
****************************************************************** */
TransportBoundaryFunction_Alquimia::~TransportBoundaryFunction_Alquimia()
{
  chem_engine_->FreeState(beaker_.properties, beaker_.state, beaker_.aux_data, beaker_.aux_output);
}


/* ******************************************************************
* Internal subroutine that defines a boundary function.
****************************************************************** */
void
TransportBoundaryFunction_Alquimia::Init_(const std::vector<std::string>& regions)
{
  for (int i = 0; i < regions.size(); ++i) {
    // Get the faces that belong to this region (since boundary conditions
    // are applied on faces).
    assert(mesh_->isValidSetName(regions[i], AmanziMesh::Entity_kind::FACE));

    auto block = mesh_->getSetEntities(
      regions[i], AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
    int nblock = block.size();

    // Now get the cells that are attached to these faces.
    for (int n = 0; n < nblock; ++n) {
      int f = block[n];
      value_[f].resize(chem_engine_->NumPrimarySpecies());
    }
  }
}


/* ******************************************************************
* Evaluate values at time.
****************************************************************** */
void
TransportBoundaryFunction_Alquimia::Compute(double t_old, double t_new)
{
  std::string cond_name = (*f_)(t_new);
  // Loop over sides and evaluate values.
  for (auto it = begin(); it != end(); ++it) {
    // Find the index of the cell we're in.
    int f = it->first;
    int cell = AmanziMesh::getFaceOnBoundaryInternalCell(*mesh_, f);

    // Dump the contents of the chemistry state into our Alquimia containers.
    alquimia_pk_->copyToAlquimia(cell, beaker_);

    // Enforce the condition.
    chem_engine_->EnforceCondition(
      cond_name, t_new, beaker_.properties, beaker_.state, beaker_.aux_data, beaker_.aux_output);

    // Move the concentrations into place.
    std::vector<double>& values = it->second;
    for (int i = 0; i < values.size(); i++) {
      values[i] = beaker_.state.total_mobile.data[i];
    }
  }
}

} // namespace Transport
} // namespace Amanzi

#endif
