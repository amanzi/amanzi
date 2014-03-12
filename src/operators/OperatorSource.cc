/*
  This is the operators component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

#include "OperatorSource.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Adds time derivative ss * (u - u0) / dT.
****************************************************************** */
void OperatorSource::UpdateMatrices(const CompositeVector& source)
{
  AmanziMesh::Entity_ID_List nodes;

  std::string name;
  std::vector<double> volume;

  if (source.HasComponent("cell")) {
    name = "cell";
    volume.resize(ncells_owned);

    for (int c = 0; c < ncells_owned; c++) {
      volume[c] = mesh_->cell_volume(c); 
    }
  } else if (source.HasComponent("face")) {
    name = "face";
  } else if (source.HasComponent("node")) {
    name = "node";
    volume.resize(nnodes_owned, 0.0);

    for (int c = 0; c < ncells_owned; c++) {
      mesh_->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();

      for (int i = 0; i < nnodes; i++) {
        volume[nodes[i]] += mesh_->cell_volume(c) / nnodes; 
      }
    }
  }

  const Epetra_MultiVector& src = *source.ViewComponent(name);

  Epetra_MultiVector& diag = *diagonal_->ViewComponent(name);
  Epetra_MultiVector& rhs = *rhs_->ViewComponent(name);

  int n = src.MyLength();
  for (int i = 0; i < n; i++) {
    rhs[0][i] += volume[i] * src[0][i];
  }
}

}  // namespace Operators
}  // namespace Amanzi
