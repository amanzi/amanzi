/*
  This is the operators component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

#include "OperatorAccumulation.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Adds time derivative ss * (u - u0) / dT.
****************************************************************** */
void OperatorAccumulation::UpdateMatrices(
    const CompositeVector& u0, const CompositeVector& ss, double dT)
{
  AmanziMesh::Entity_ID_List nodes;

  std::string name;
  CompositeVector entity_volume(ss);

  if (ss.HasComponent("cell")) {
    name = "cell";
    Epetra_MultiVector& volume = *entity_volume.ViewComponent(name); 

    for (int c = 0; c < ncells_owned; c++) {
      volume[0][c] = mesh_->cell_volume(c); 
    }
  } else if (ss.HasComponent("face")) {
    name = "face";
    // Missing code.
  } else if (ss.HasComponent("node")) {
    name = "node";
    Epetra_MultiVector& volume = *entity_volume.ViewComponent(name, true); 
    volume.PutScalar(0.0);

    for (int c = 0; c < ncells_owned; c++) {
      mesh_->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();

      for (int i = 0; i < nnodes; i++) {
        volume[0][nodes[i]] += mesh_->cell_volume(c) / nnodes; 
      }
    }

    entity_volume.GatherGhostedToMaster(name);
  }

  const Epetra_MultiVector& u0c = *u0.ViewComponent(name);
  const Epetra_MultiVector& ssc = *ss.ViewComponent(name);

  Epetra_MultiVector& volume = *entity_volume.ViewComponent(name); 
  Epetra_MultiVector& diag = *diagonal_->ViewComponent(name);
  Epetra_MultiVector& rhs = *rhs_->ViewComponent(name);

  int n = u0c.MyLength();
  for (int i = 0; i < n; i++) {
    double factor = volume[0][i] * ssc[0][i] / dT;
    diag[0][i] += factor;
    rhs[0][i] += factor * u0c[0][i];
  }
}

}  // namespace Operators
}  // namespace Amanzi
