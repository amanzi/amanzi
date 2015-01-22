/*
  This is the operators component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

#include "Epetra_Vector.h"

#include "CompositeVector.hh"
#include "SuperMap.hh"
#include "OperatorUtils.hh"

namespace Amanzi {
namespace Operators {

// composite vector to/from super vector
int
CopyCompositeVectorToSuperVector(const SuperMap& smap, const CompositeVector& cv,
        Epetra_Vector& sv, int dofnum) {
  if (cv.HasComponent("face")) {
    const std::vector<int>& face_inds = smap.Indices("face", dofnum);
    const Epetra_MultiVector& data = *cv.ViewComponent("face");
    for (int f=0; f!=data.MyLength(); ++f) sv[face_inds[f]] = data[0][f];
  } 
  if (cv.HasComponent("cell")) {
    const std::vector<int>& cell_inds = smap.Indices("cell", dofnum);
    const Epetra_MultiVector& data = *cv.ViewComponent("cell");
    for (int c=0; c!=data.MyLength(); ++c) sv[cell_inds[c]] = data[0][c];
  } 
  if (cv.HasComponent("node")) {
    const std::vector<int>& node_inds = smap.Indices("node", dofnum);
    const Epetra_MultiVector& data = *cv.ViewComponent("node");
    for (int v=0; v!=data.MyLength(); ++v) sv[node_inds[v]] = data[0][v];
  }
  return 0;
}

int
CopyCompositeVectorToSuperVector(const SuperMap& smap, const CompositeVector& cv,
        Teuchos::RCP<Epetra_Vector>& sv, int dofnum) {
  if (sv == Teuchos::null) {
    sv = Teuchos::rcp(new Epetra_Vector(*smap.Map()));
  }
  return CopyCompositeVectorToSuperVector(smap, cv, *sv, dofnum);
}

int
CopySuperVectorToCompositeVector(const SuperMap& smap, const Epetra_Vector& sv,
        CompositeVector& cv, int dofnum) {

  if (cv.HasComponent("face")) {
    const std::vector<int>& face_inds = smap.Indices("face", dofnum);
    Epetra_MultiVector& data = *cv.ViewComponent("face");
    for (int f=0; f!=data.MyLength(); ++f) data[0][f] = sv[face_inds[f]];
  } 
  if (cv.HasComponent("cell")) {
    const std::vector<int>& cell_inds = smap.Indices("cell", dofnum);
    Epetra_MultiVector& data = *cv.ViewComponent("cell");
    for (int c=0; c!=data.MyLength(); ++c) data[0][c] = sv[cell_inds[c]];
  } 
  if (cv.HasComponent("node")) {
    const std::vector<int>& node_inds = smap.Indices("node", dofnum);
    Epetra_MultiVector& data = *cv.ViewComponent("node");
    for (int v=0; v!=data.MyLength(); ++v) data[0][v] = sv[node_inds[v]];
  } 
  return 0;
}


} // namespace Operators
} // namespace Amanzi
