/*
  This is the data structures component of the Amanzi code.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

#include "SuperMapLumpedSurface.hh"

namespace Amanzi {
namespace Operators {

// Constructor
SuperMapLumpedSurface::SuperMapLumpedSurface(const SuperMapLumped& map,
        const Teuchos::RCP<const AmanziMesh::Mesh>& surf_mesh):
    surf_mesh_(surf_mesh),
    SuperMapLumped(map) {}

const std::vector<int>&
SuperMapLumpedSurface::CreateIndices_(const std::string& surf_compname, int dofnum, bool ghosted) const
{
  if (surf_compname == std::string("cell")) {
    std::string compname("face"); // surface cells correspond to subsurface faces

    if (ghosted) {
      int nentities_owned = counts_.at(compname);
      int nentities = nentities_owned + ghosted_counts_.at(compname);

      int offset = offsets_.at(compname);
      int num_dof = num_dofs_.at(compname);

      int surf_nentities_owned = surf_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
      int surf_nentities = surf_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

      std::vector<int> indices(surf_nentities, -1);
      for (int i=0; i!=surf_nentities_owned; ++i) {
        indices[i] = offset + dofnum + surf_mesh_->entity_get_parent(AmanziMesh::CELL, i) * num_dof;
      }
      
      int ghosted_offset = ghosted_offsets_.at(compname);
      for (int i=surf_nentities_owned; i!=surf_nentities; ++i) {
        indices[i] = ghosted_offset + dofnum
            + (surf_mesh_->entity_get_parent(AmanziMesh::CELL, i)-nentities_owned)*num_dof;
      }

      // assign
      ghosted_indices_[compname][dofnum] = indices;
      return ghosted_indices_[compname][dofnum];

    } else {
      int nentities_owned = counts_.at(compname);
      int nentities = nentities_owned + ghosted_counts_.at(compname);

      int offset = offsets_.at(compname);
      int num_dof = num_dofs_.at(compname);

      int surf_nentities = surf_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

      std::vector<int> indices(surf_nentities, -1);
      for (int i=0; i!=surf_nentities; ++i) {
        indices[i] = offset + dofnum + surf_mesh_->entity_get_parent(AmanziMesh::CELL, i) * num_dof;
      }

      // assign
      indices_[compname][dofnum] = indices;
      return indices_[compname][dofnum];
    }
    
  } else {
    Errors::Message msg("SuperMapLumpedSurface only provides CELL entities.");
    Exceptions::amanzi_throw(msg);
  }
}

} // namespace Operators
} // namespace Amanzi
