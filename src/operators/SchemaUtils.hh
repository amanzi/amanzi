/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  local <-> global communications of DOFs using schema
*/

#ifndef AMANZI_OPERATOR_SCHEMA_UTILS_HH_
#define AMANZI_OPERATOR_SCHEMA_UTILS_HH_

#include "DenseVector.hh"
#include "CompositeVector.hh"
#include "Mesh.hh"

#include "Schema.hh"

namespace Amanzi {
namespace Operators {

void
ExtractVectorFaceOp(int f,
                    const AmanziMesh::Mesh& mesh,
                    const Schema& schema,
                    WhetStone::DenseVector& v,
                    const CompositeVector& X);
void
AssembleVectorFaceOp(int f,
                     const AmanziMesh::Mesh& mesh,
                     const Schema& schema,
                     const WhetStone::DenseVector& v,
                     CompositeVector& X);

void
ExtractVectorNodeOp(int n,
                    const AmanziMesh::Mesh& mesh,
                    const Schema& schema,
                    WhetStone::DenseVector& v,
                    const CompositeVector& X);
void
AssembleVectorNodeOp(int n,
                     const AmanziMesh::Mesh& mesh,
                     const Schema& schema,
                     const WhetStone::DenseVector& v,
                     CompositeVector& X);

void
ExtractVectorCellOp(int c,
                    const AmanziMesh::Mesh& mesh,
                    const Schema& schema,
                    WhetStone::DenseVector& v,
                    const CompositeVector& X);
void
AssembleVectorCellOp(int c,
                     const AmanziMesh::Mesh& mesh,
                     const Schema& schema,
                     const WhetStone::DenseVector& v,
                     CompositeVector& X);

} // namespace Operators
} // namespace Amanzi

#endif
