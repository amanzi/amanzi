/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Operator those domain and range are defined by two schemas and
  respected CVSs.
*/

#include "SchemaUtils.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Assemble local vector to the global CV when the base is cell.
****************************************************************** */
void
AssembleVectorCellOp(int c,
                     const AmanziMesh::Mesh& mesh,
                     const Schema& schema,
                     const WhetStone::DenseVector& v,
                     CompositeVector& X)
{
  int m(0);
  for (auto it = schema.begin(); it != schema.end(); ++it) {
    int num;
    AmanziMesh::Entity_kind kind;
    std::tie(kind, std::ignore, num) = *it;

    if (kind == AmanziMesh::Entity_kind::NODE) {
      Epetra_MultiVector& Xn = *X.ViewComponent("node", true);

      auto nodes = mesh.getCellNodes(c);
      int nnodes = nodes.size();

      for (int n = 0; n != nnodes; ++n) {
        for (int k = 0; k < num; ++k) { Xn[k][nodes[n]] += v(m++); }
      }
    }

    else if (kind == AmanziMesh::Entity_kind::FACE) {
      Epetra_MultiVector& Xf = *X.ViewComponent("face", true);

      const auto& faces = mesh.getCellFaces(c);
      int nfaces = faces.size();

      for (int n = 0; n != nfaces; ++n) {
        for (int k = 0; k < num; ++k) { Xf[k][faces[n]] += v(m++); }
      }
    }

    else if (kind == AmanziMesh::Entity_kind::EDGE) {
      Epetra_MultiVector& Xe = *X.ViewComponent("edge", true);

      auto edges = mesh.getCellEdges(c);
      int nedges = edges.size();

      for (int n = 0; n != nedges; ++n) {
        for (int k = 0; k < num; ++k) { Xe[k][edges[n]] += v(m++); }
      }
    }

    else if (kind == AmanziMesh::Entity_kind::CELL) {
      Epetra_MultiVector& Xc = *X.ViewComponent("cell", true);

      for (int k = 0; k < num; ++k) { Xc[k][c] += v(m++); }
    }

    else {
      AMANZI_ASSERT(false);
    }
  }
}


/* ******************************************************************
* Assemble local vector to the global CV when the base is face.
****************************************************************** */
void
AssembleVectorFaceOp(int f,
                     const AmanziMesh::Mesh& mesh,
                     const Schema& schema,
                     const WhetStone::DenseVector& v,
                     CompositeVector& X)
{
  int m(0);
  for (auto it = schema.begin(); it != schema.end(); ++it) {
    int num;
    AmanziMesh::Entity_kind kind;
    std::tie(kind, std::ignore, num) = *it;

    if (kind == AmanziMesh::Entity_kind::CELL) {
      Epetra_MultiVector& Xf = *X.ViewComponent("cell", true);

      auto cells = mesh.getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
      int ncells = cells.size();

      for (int n = 0; n != ncells; ++n) {
        for (int k = 0; k < num; ++k) { Xf[k][cells[n]] += v(m++); }
      }
    }
  }
}


/* ******************************************************************
* Assemble local vector to the global CV when the base is node.
****************************************************************** */
void
AssembleVectorNodeOp(int n,
                     const AmanziMesh::Mesh& mesh,
                     const Schema& schema,
                     const WhetStone::DenseVector& v,
                     CompositeVector& X)
{
  int m(0);
  for (auto it = schema.begin(); it != schema.end(); ++it) {
    int num;
    AmanziMesh::Entity_kind kind;
    std::tie(kind, std::ignore, num) = *it;

    if (kind == AmanziMesh::Entity_kind::CELL) {
      Epetra_MultiVector& Xc = *X.ViewComponent("cell", true);

      auto cells = mesh.getNodeCells(n, AmanziMesh::Parallel_kind::ALL);
      int ncells = cells.size();

      for (int i = 0; i != ncells; ++i) {
        for (int k = 0; k < num; ++k) { Xc[k][cells[i]] += v(m++); }
      }
    }
  }
}


/* ******************************************************************
* Extract local vector from the global CV when the base is cell.
****************************************************************** */
void
ExtractVectorCellOp(int c,
                    const AmanziMesh::Mesh& mesh,
                    const Schema& schema,
                    WhetStone::DenseVector& v,
                    const CompositeVector& X)
{
  int m(0);
  for (auto it = schema.begin(); it != schema.end(); ++it) {
    int num;
    AmanziMesh::Entity_kind kind;
    std::tie(kind, std::ignore, num) = *it;

    if (kind == AmanziMesh::Entity_kind::NODE) {
      const Epetra_MultiVector& Xn = *X.ViewComponent("node", true);

      auto nodes = mesh.getCellNodes(c);
      int nnodes = nodes.size();

      for (int n = 0; n != nnodes; ++n) {
        for (int k = 0; k < num; ++k) { v(m++) = Xn[k][nodes[n]]; }
      }
    }

    else if (kind == AmanziMesh::Entity_kind::FACE) {
      const Epetra_MultiVector& Xf = *X.ViewComponent("face", true);

      const auto& faces = mesh.getCellFaces(c);
      int nfaces = faces.size();

      for (int n = 0; n != nfaces; ++n) {
        for (int k = 0; k < num; ++k) { v(m++) = Xf[k][faces[n]]; }
      }
    }

    else if (kind == AmanziMesh::Entity_kind::EDGE) {
      const Epetra_MultiVector& Xe = *X.ViewComponent("edge", true);

      auto edges = mesh.getCellEdges(c);
      int nedges = edges.size();

      for (int n = 0; n != nedges; ++n) {
        for (int k = 0; k < num; ++k) { v(m++) = Xe[k][edges[n]]; }
      }
    }

    else if (kind == AmanziMesh::Entity_kind::CELL) {
      const Epetra_MultiVector& Xc = *X.ViewComponent("cell", true);

      for (int k = 0; k < num; ++k) { v(m++) = Xc[k][c]; }
    }

    else {
      AMANZI_ASSERT(false);
    }
  }
}


/* ******************************************************************
* Extract local vector from the global CV when the base is face.
****************************************************************** */
void
ExtractVectorFaceOp(int f,
                    const AmanziMesh::Mesh& mesh,
                    const Schema& schema,
                    WhetStone::DenseVector& v,
                    const CompositeVector& X)
{
  int m(0);
  for (auto it = schema.begin(); it != schema.end(); ++it) {
    int num;
    AmanziMesh::Entity_kind kind;
    std::tie(kind, std::ignore, num) = *it;

    if (kind == AmanziMesh::Entity_kind::CELL) {
      const Epetra_MultiVector& Xf = *X.ViewComponent("cell", true);

      auto cells = mesh.getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
      int ncells = cells.size();

      for (int n = 0; n != ncells; ++n) {
        for (int k = 0; k < num; ++k) { v(m++) = Xf[k][cells[n]]; }
      }
    }
  }
}


/* ******************************************************************
* Extract local vector from the global CV when the base is node.
****************************************************************** */
void
ExtractVectorNodeOp(int n,
                    const AmanziMesh::Mesh& mesh,
                    const Schema& schema,
                    WhetStone::DenseVector& v,
                    const CompositeVector& X)
{
  int m(0);
  for (auto it = schema.begin(); it != schema.end(); ++it) {
    int num;
    AmanziMesh::Entity_kind kind;
    std::tie(kind, std::ignore, num) = *it;

    if (kind == AmanziMesh::Entity_kind::CELL) {
      const Epetra_MultiVector& Xc = *X.ViewComponent("cell", true);

      auto cells = mesh.getNodeCells(n, AmanziMesh::Parallel_kind::ALL);
      int ncells = cells.size();

      for (int i = 0; i != ncells; ++i) {
        for (int k = 0; k < num; ++k) { v(m++) = Xc[k][cells[i]]; }
      }
    }
  }
}

} // namespace Operators
} // namespace Amanzi
