/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Remap of 2-cochains and 3-cochains with supporting tools 
  for physics packages.
*/

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "dbc.hh"
#include "LimiterCell.hh"
#include "MFD3D_Diffusion.hh"
#include "MFD3D_Electromagnetics.hh"
#include "Point.hh"
#include "ReconstructionCellLinear.hh"

#include "RemapCochains.hh"

namespace Amanzi {
namespace Operators {

/* *****************************************************************
* Extract degree of cochain.
***************************************************************** */
RemapCochains::RemapCochains(const Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
                             const Teuchos::RCP<const AmanziMesh::Mesh> mesh1)
  : mesh0_(mesh0), mesh1_(mesh1)
{
  d_ = mesh0_->getSpaceDimension();
}


/* *****************************************************************
* Low-order (donor) remap of 1-cochains
***************************************************************** */
template<>
void RemapCochains::LowOrder<1>(const CompositeVector& in, CompositeVector& out) 
{
  int nedges_owned = mesh0_->getNumEntities(AmanziMesh::EDGE, AmanziMesh::Parallel_kind::OWNED);
  int nnodes_wghost = mesh0_->getNumEntities(AmanziMesh::NODE, AmanziMesh::Parallel_kind::ALL);

  const auto& in_e = *in.ViewComponent("edge", true);
  auto& out_e = *out.ViewComponent("edge");

  // compute cell-centered vector field from co-chains
  const auto& vc = getCellFieldFromEdges(mesh0_, in);
  const auto& vc_c = *vc->ViewComponent("cell", true);

  // compute donor cells 
  double eps = 1e-8;
  AmanziGeometry::Point xprobe(d_), val(d_);
  std::vector<int> donor_cell(nnodes_wghost);  // FIXME in parallel

  for (int v = 0; v < nnodes_wghost; ++v) {
    const AmanziGeometry::Point& p0 = mesh0_->getNodeCoordinate(v);
    const AmanziGeometry::Point& p1 = mesh1_->getNodeCoordinate(v);

    xprobe = p0 + eps * (p1 - p0);

    const auto& cells = mesh0_->getNodeCells(v);
    int ncells = cells.size();

    donor_cell[v] = cells[0];
    for (int n = 0; n < ncells; ++n) {
      if (mesh0_->isPointInCell(xprobe, cells[n])) {
        donor_cell[v] = cells[n];
        break; 
      }
    }
  }

  // compute 1-cochains on mesh1
  for (int e = 0; e < nedges_owned; ++e) {
    int dir = 1;
    double sum = in_e[0][e];

    const auto& nodes = mesh0_->getEdgeNodes(e);
    for (int n = 0; n < 2; ++n) {
      int v = nodes[n];
      int c = donor_cell[v];

      const AmanziGeometry::Point& p0 = mesh0_->getNodeCoordinate(v);
      const AmanziGeometry::Point& p1 = mesh1_->getNodeCoordinate(v);
      for (int i = 0; i < d_; ++i) val[i] = vc_c[i][c];
      sum -= val * (p1 - p0) * dir;
      dir = -dir;
    }
    out_e[0][e] = sum;
  }
}


/* *****************************************************************
* Low-order (donor) remap of 2-cochains
***************************************************************** */
template<>
void RemapCochains::LowOrder<2>(const CompositeVector& in, CompositeVector& out) 
{
  int nfaces_owned = mesh0_->getNumEntities(AmanziMesh::FACE, AmanziMesh::Parallel_kind::OWNED);
  int ncells_wghost = mesh0_->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::ALL);
  int nedges_wghost = mesh0_->getNumEntities(AmanziMesh::EDGE, AmanziMesh::Parallel_kind::ALL);

  const auto& in_f = *in.ViewComponent("face", true);
  auto& out_f = *out.ViewComponent("face");

  // compute cell-centered vector field
  const auto& vc = getCellFieldFromFaces(mesh0_, in);
  const auto& vc_c = *vc->ViewComponent("cell", true);

  // compute swept areas and donor cells 
  double eps = 1e-8;
  AmanziGeometry::Point xa(d_), xprobe(d_), val(d_);

  std::vector<AmanziGeometry::Point> swept_area(nedges_wghost);
  std::vector<int> donor_cell(nedges_wghost);  // FIXME in parallel

  for (int e = 0; e < nedges_wghost; ++e) {
    const auto& nodes0 = mesh0_->getEdgeNodes(e);
    const auto& nodes1 = mesh1_->getEdgeNodes(e);

    const AmanziGeometry::Point& a0 = mesh0_->getNodeCoordinate(nodes0[0]);
    const AmanziGeometry::Point& a1 = mesh0_->getNodeCoordinate(nodes0[1]);

    const AmanziGeometry::Point& b0 = mesh1_->getNodeCoordinate(nodes1[0]);
    const AmanziGeometry::Point& b1 = mesh1_->getNodeCoordinate(nodes1[1]);

    swept_area[e] = 0.5 * ((b1 - a0)^(b0 - a1));

    const AmanziGeometry::Point& tau = mesh0_->getEdgeVector(e);
    const AmanziGeometry::Point& xe = mesh0_->getEdgeCentroid(e);

    xa = 0.25 * (a0 + a1 + b0 + b1);
    xa = xa - (xa * tau) / L22(tau) * tau;
    xprobe = xe + eps * xa;

    const auto& cells = mesh0_->getEdgeCells(e);
    int ncells = cells.size();

    donor_cell[e] = cells[0];
    for (int n = 0; n < ncells; ++n) {
      if (mesh0_->isPointInCell(xprobe, cells[n])) {
        donor_cell[e] = cells[n];
        break; 
      }
    }
  }

  // compute 2-cochains on mesh1
  for (int f = 0; f < nfaces_owned; ++f) {
    double sum = in_f[0][f];

    const auto& [edges, dirs] = mesh0_->getFaceEdgesAndDirections(f);
    int nedges = edges.size();
    for (int n = 0; n < nedges; ++n) {
      int e = edges[n];
      int c = donor_cell[e];
      for (int i = 0; i < d_; ++i) val[i] = vc_c[i][c];
      sum -= val * swept_area[e] * dirs[n];
    }
    out_f[0][f] = sum;
  }
}


/* *****************************************************************
* High-order remap of 2-cochains
***************************************************************** */
template<>
void RemapCochains::HighOrder<2>(const CompositeVector& in, CompositeVector& out) 
{
  int nfaces_owned = mesh0_->getNumEntities(AmanziMesh::FACE, AmanziMesh::Parallel_kind::OWNED);
  int ncells_wghost = mesh0_->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::ALL);
  int nedges_wghost = mesh0_->getNumEntities(AmanziMesh::EDGE, AmanziMesh::Parallel_kind::ALL);

  // projection to constant vectors in cells
  const auto& in_f = *in.ViewComponent("face", true);
  auto& out_f = *out.ViewComponent("face");

  const auto& vc = getCellFieldFromFaces(mesh0_, in);
  auto vc_c = vc->ViewComponent("cell");

  // component-wise linear reconstruction in cells
  Teuchos::ParameterList plist;
  std::vector<Teuchos::RCP<Operators::ReconstructionCellLinear>> grad(d_);

  auto limiter = Teuchos::rcp(new Operators::LimiterCell(mesh0_));
  limiter->Init(plist);

  for (int i = 0; i < d_; ++i) {
    grad[i] = Teuchos::rcp(new Operators::ReconstructionCellLinear(mesh0_));
    grad[i]->Init(plist);
    grad[i]->Compute(vc_c, i);
    // limiter->ApplyLimiter(vc, 0, grad, bc_model_vector, bc_value_qx);
    grad[i]->data()->ScatterMasterToGhosted("cell");
  }

  // compute swept areas and donor cells 
  double eps = 1e-8;
  AmanziGeometry::Point xa(d_), xprobe(d_), val(d_);

  std::vector<AmanziGeometry::Point> swept_area(nedges_wghost);
  std::vector<int> donor_cell(nedges_wghost);  // FIXME in parallel

  for (int e = 0; e < nedges_wghost; ++e) {
    const auto& nodes0 = mesh0_->getEdgeNodes(e);
    const auto& nodes1 = mesh1_->getEdgeNodes(e);

    const AmanziGeometry::Point& a0 = mesh0_->getNodeCoordinate(nodes0[0]);
    const AmanziGeometry::Point& a1 = mesh0_->getNodeCoordinate(nodes0[1]);

    const AmanziGeometry::Point& b0 = mesh1_->getNodeCoordinate(nodes1[0]);
    const AmanziGeometry::Point& b1 = mesh1_->getNodeCoordinate(nodes1[1]);

    swept_area[e] = 0.5 * ((b1 - a0)^(b0 - a1));

    const AmanziGeometry::Point& tau = mesh0_->getEdgeVector(e);
    const AmanziGeometry::Point& xe = mesh0_->getEdgeCentroid(e);

    xa = 0.25 * (a0 + a1 + b0 + b1);
    xa = xa - (xa * tau) / L22(tau) * tau;
    xprobe = xe + eps * xa;

    const auto& cells = mesh0_->getEdgeCells(e);
    int ncells = cells.size();

    donor_cell[e] = cells[0];
    for (int n = 0; n < ncells; ++n) {
      if (mesh0_->isPointInCell(xprobe, cells[n])) {
        donor_cell[e] = cells[n];
        break; 
      }
    }
  }

  // compute 2-cochains on mesh1
  for (int f = 0; f < nfaces_owned; ++f) {
    double sum = in_f[0][f];

    const auto& [edges, dirs] = mesh0_->getFaceEdgesAndDirections(f);
    int nedges = edges.size();
    for (int n = 0; n < nedges; ++n) {
      int e = edges[n];
      int c = donor_cell[e];

      const auto& p = getSweptEdgeCentroid(e);

      for (int i = 0; i < d_; ++i) {
        val[i] = grad[i]->getValue(c, p);
      }
      sum -= val * swept_area[e] * dirs[n];
    }
    out_f[0][f] = sum;
  }
}


/* *****************************************************************
* Cell reconstruction from face-based data
***************************************************************** */
Teuchos::RCP<CompositeVector>
RemapCochains::getCellFieldFromFaces(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                     const CompositeVector& in)
{
  in.ScatterMasterToGhosted();
  const auto& in_f = *in.ViewComponent("face", true);

  WhetStone::MFD3D_Diffusion mfd(mesh);
  WhetStone::Polynomial linear(d_, 1);

  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh0_)->SetGhosted(true)->AddComponent("cell", AmanziMesh::Entity_kind::CELL, d_);
  auto vc = Teuchos::rcp(new CompositeVector(cvs));
  auto& vc_c = *vc->ViewComponent("cell", true);

  // to same on communications, we reconstruct in ghost cells
  int ncells_wghost = mesh->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::ALL);
  for (int c = 0; c < ncells_wghost; ++c) {
    const auto& faces = mesh->getCellFaces(c);
    int nfaces = faces.size();

    std::vector<WhetStone::Polynomial> dofs(nfaces);
    for (int n = 0; n < nfaces; ++n) {
      dofs[n].Reshape(d_, 0);
      dofs[n](0) = in_f[0][faces[n]];
    }
    mfd.L2Cell(c, dofs, dofs, nullptr, linear);
    for (int i = 0; i < d_; ++i) vc_c[i][c] = -linear(i + 1);
  }

  return vc;
}


/* *****************************************************************
* Cell reconstruction from edge-based data
***************************************************************** */
Teuchos::RCP<CompositeVector>
RemapCochains::getCellFieldFromEdges(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                     const CompositeVector& in)
{
  in.ScatterMasterToGhosted();
  const auto& in_e = *in.ViewComponent("edge", true);

  WhetStone::MFD3D_Electromagnetics mfd(mesh);
  WhetStone::Polynomial linear(d_, 1);

  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh0_)->SetGhosted(true)->AddComponent("cell", AmanziMesh::Entity_kind::CELL, d_);
  auto vc = Teuchos::rcp(new CompositeVector(cvs));
  auto& vc_c = *vc->ViewComponent("cell", true);

  // to same on communications, we reconstruct in ghost cells
  int ncells_wghost = mesh->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::ALL);
  for (int c = 0; c < ncells_wghost; ++c) {
    const auto& edges = mesh->getCellEdges(c);
    int nedges = edges.size();

    std::vector<WhetStone::Polynomial> dofs(nedges);
    for (int n = 0; n < nedges; ++n) {
      int e = edges[n];
      dofs[n].Reshape(d_, 0);
      dofs[n](0) = in_e[0][e] / mesh->getEdgeLength(e);
    }
    mfd.L2Cell(c, dofs, dofs, nullptr, linear);
    for (int i = 0; i < d_; ++i) vc_c[i][c] = linear(i + 1);
  }

  return vc;
}


/* *****************************************************************
* Supporting tools for swept regions
***************************************************************** */
AmanziGeometry::Point
RemapCochains::getSweptEdgeCentroid(int e) const
{
  const auto& nodes0 = mesh0_->getEdgeNodes(e);
  const auto& nodes1 = mesh1_->getEdgeNodes(e);

  const AmanziGeometry::Point& v0 = mesh0_->getNodeCoordinate(nodes0[0]);
  const AmanziGeometry::Point& v1 = mesh0_->getNodeCoordinate(nodes0[1]);

  const AmanziGeometry::Point& v2 = mesh1_->getNodeCoordinate(nodes1[0]);
  const AmanziGeometry::Point& v3 = mesh1_->getNodeCoordinate(nodes1[1]);

  return (v0 + v1 + v2 + v3) / 4;

  /*
  auto b = v1 - v0;
  auto c = v3 - v0;
  auto d = v2 - v0;
 
  double s1 = norm(b ^ c);
  double s2 = norm(d ^ c);

  if (s1 + s2 <= 1.0e-12) return 0.5 * (v0 + v1);
  else return v0 + (s1 * (c + b) + s2 * (c + d)) / (3 * (s1 + s2));
  */
}

} // namespace Operators
} // namespace Amanzi

