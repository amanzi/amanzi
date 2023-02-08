/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Collection of mesh deformation tools.
*/

#ifndef AMANZI_OPERATOR_DEFORM_MESH_HH_
#define AMANZI_OPERATOR_DEFORM_MESH_HH_

#include "Mesh.hh"
#include "MeshCurved.hh"
#include "CompositeVector.hh"

#include "OperatorDefs.hh"

namespace Amanzi {

// collection of routines for mesh deformation
void
DeformMesh(const Teuchos::RCP<AmanziMesh::Mesh>& mesh1,
           int deform,
           double t,
           const Teuchos::RCP<const AmanziMesh::Mesh>& mesh0 = Teuchos::null);

void
DeformMeshCurved(const Teuchos::RCP<AmanziMesh::Mesh>& mesh1,
                 int deform,
                 double t,
                 const Teuchos::RCP<const AmanziMesh::Mesh>& mesh0,
                 int order);

AmanziGeometry::Point
MovePoint(double t, const AmanziGeometry::Point& xv, int deform);

AmanziGeometry::Point
TaylorGreenVortex(double t, const AmanziGeometry::Point& xv);
AmanziGeometry::Point
Rotation2D(double t, const AmanziGeometry::Point& xv);
AmanziGeometry::Point
CompressionExpansion(double t, const AmanziGeometry::Point& xv);
AmanziGeometry::Point
BubbleFace3D(double t, const AmanziGeometry::Point& xv);
AmanziGeometry::Point
Unused(double t, const AmanziGeometry::Point& xv);
AmanziGeometry::Point
SineProduct(double t, const AmanziGeometry::Point& xv);


/* *****************************************************************
* Deform mesh1 using coordinates given by mesh0
***************************************************************** */
inline void
DeformMesh(const Teuchos::RCP<AmanziMesh::Mesh>& mesh1,
           int deform,
           double t,
           const Teuchos::RCP<const AmanziMesh::Mesh>& mesh0)
{
  if (mesh1->getComm()->MyPID() == 0) std::cout << "Deforming mesh...\n";

  // create distributed random vector
  int d = mesh1->getSpaceDimension();

  // consistent parallel data are needed for the random mesh deformation
  AmanziMesh::Entity_ID_View bnd_ids;
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh1)->SetGhosted(true)->AddComponent("node", AmanziMesh::Entity_kind::NODE, d);
  CompositeVector random(cvs);
  Epetra_MultiVector& random_n = *random.ViewComponent("node", true);

  int gid = mesh1->getMap(AmanziMesh::Entity_kind::NODE, false).MaxAllGID();
  double scale = 0.2 * std::pow(gid, -2.0 / d);

  if (deform == 7) {
    random_n.Random();
    random_n.Scale(scale);
    random.ScatterMasterToGhosted();

    std::vector<double> vofs;
    mesh1->getSetEntitiesAndVolumeFractions(
      "Boundary", AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::ALL, &bnd_ids, &vofs);
  }

  // relocate mesh nodes
  AmanziGeometry::Point xv(d), yv(d), uv(d);
  AmanziMesh::Entity_ID_View nodeids;
  AmanziGeometry::Point_View new_positions, final_positions;

  int nnodes = mesh1->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::ALL);

  for (int v = 0; v < nnodes; ++v) {
    if (mesh0.get())
      xv = mesh0->getNodeCoordinate(v);
    else
      xv = mesh1->getNodeCoordinate(v);

    nodeids.push_back(v);

    if (deform == 7) {
      yv = xv;
      if (std::find(bnd_ids.begin(), bnd_ids.end(), v) == bnd_ids.end()) {
        for (int i = 0; i < d; ++i) yv[i] += random_n[i][v];
      }
    } else {
      yv = MovePoint(t, xv, deform);
    }

    new_positions.push_back(yv);
  }
  mesh1->deform(nodeids, new_positions, false, &final_positions);
}


/* *****************************************************************
* Deform high-order (curved) mesh1
***************************************************************** */
inline void
DeformMeshCurved(const Teuchos::RCP<AmanziMesh::Mesh>& mesh1,
                 int deform,
                 double t,
                 const Teuchos::RCP<const AmanziMesh::Mesh>& mesh0,
                 int order)
{
  DeformMesh(mesh1, deform, t, mesh0);

  int dim = mesh1->getSpaceDimension();
  if (order > 1) {
    int nfaces =
      mesh0->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);
    auto ho_nodes0f = std::make_shared<std::vector<AmanziGeometry::Point_View>>(nfaces);
    auto ho_nodes1f = std::make_shared<std::vector<AmanziGeometry::Point_View>>(nfaces);

    for (int f = 0; f < nfaces; ++f) {
      const AmanziGeometry::Point& xf = mesh0->getFaceCentroid(f);
      (*ho_nodes0f)[f].push_back(xf);

      auto yv = MovePoint(t, xf, deform);
      (*ho_nodes1f)[f].push_back(yv);
    }
    auto tmp = Teuchos::rcp_const_cast<AmanziMesh::Mesh>(mesh0);
    Teuchos::rcp_static_cast<AmanziMesh::MeshCurved>(tmp)->set_face_ho_nodes(ho_nodes0f);
    Teuchos::rcp_static_cast<AmanziMesh::MeshCurved>(mesh1)->set_face_ho_nodes(ho_nodes1f);

    if (dim == 3) {
      int nedges =
        mesh0->getNumEntities(AmanziMesh::Entity_kind::EDGE, AmanziMesh::Parallel_kind::ALL);
      auto ho_nodes0e = std::make_shared<std::vector<AmanziGeometry::Point_View>>(nedges);
      auto ho_nodes1e = std::make_shared<std::vector<AmanziGeometry::Point_View>>(nedges);

      for (int e = 0; e < nedges; ++e) {
        const AmanziGeometry::Point& xe = mesh0->getEdgeCentroid(e);
        (*ho_nodes0e)[e].push_back(xe);

        auto yv = MovePoint(t, xe, deform);
        (*ho_nodes1e)[e].push_back(yv);
      }
      Teuchos::rcp_static_cast<AmanziMesh::MeshCurved>(tmp)->set_edge_ho_nodes(ho_nodes0e);
      Teuchos::rcp_static_cast<AmanziMesh::MeshCurved>(mesh1)->set_edge_ho_nodes(ho_nodes1e);
    }
  }
}


/* *****************************************************************
* Factory of point relocation methods
***************************************************************** */
inline AmanziGeometry::Point
MovePoint(double t, const AmanziGeometry::Point& xv, int deform)
{
  AmanziGeometry::Point yv(xv);

  switch (deform) {
  case 1:
    yv = TaylorGreenVortex(t, xv);
    break;
  case 2:
    yv = Unused(t, xv);
    break;
  case 3:
    yv = SineProduct(t, xv);
    break;
  case 4:
    AMANZI_ASSERT(false);
    break;
  case 5:
    yv = CompressionExpansion(t, xv);
    break;
  case 6:
    yv = Rotation2D(t, xv);
    break;
  case 8:
    yv = BubbleFace3D(t, xv);
    break;
  default:
    break;
  }

  return yv;
}


/* *****************************************************************
* Taylor-Green vortex
***************************************************************** */
inline AmanziGeometry::Point
TaylorGreenVortex(double t, const AmanziGeometry::Point& xv)
{
  int d = xv.dim();
  AmanziGeometry::Point uv(d), yv(xv);

  double ds(0.0001);
  int n = t / ds;
  for (int i = 0; i < n; ++i) {
    if (d == 2) {
      uv[0] = 0.2 * std::sin(M_PI * yv[0]) * std::cos(M_PI * yv[1]);
      uv[1] = -0.2 * std::cos(M_PI * yv[0]) * std::sin(M_PI * yv[1]);
    } else {
      uv[0] = 0.2 * std::sin(M_PI * yv[0]) * std::cos(M_PI * yv[1]) * std::sin(M_PI * yv[2] / 2);
      uv[1] = -0.2 * std::cos(M_PI * yv[0]) * std::sin(M_PI * yv[1]) * std::sin(M_PI * yv[2] / 2);
      uv[2] = 0.0;
    }
    yv += uv * ds;
  }

  return yv;
}


/* *****************************************************************
* Rotation
***************************************************************** */
inline AmanziGeometry::Point
Rotation2D(double t, const AmanziGeometry::Point& xv)
{
  double phi = t * 2 * M_PI;
  double cs(std::cos(phi)), sn(std::sin(phi));

  AmanziGeometry::Point yv(xv);
  yv[0] = cs * xv[0] - sn * xv[1];
  yv[1] = sn * xv[0] + cs * xv[1];
  return yv;
}


/* *****************************************************************
* Compression/Expansion
***************************************************************** */
inline AmanziGeometry::Point
CompressionExpansion(double t, const AmanziGeometry::Point& xv)
{
  AmanziGeometry::Point yv(xv);
  int d = xv.dim();

  double factor(t);
  for (int i = 0; i < d; ++i) factor *= yv[i];

  for (int i = 0; i < d; ++i) yv[i] += factor * (1.0 - yv[i]) / 2;
  return yv;
}


/* *****************************************************************
* Bubble face in Z-direction
***************************************************************** */
inline AmanziGeometry::Point
BubbleFace3D(double t, const AmanziGeometry::Point& xv)
{
  AmanziGeometry::Point yv(xv);
  yv[2] += 8 * xv[0] * xv[1] * xv[2] * (1.0 - xv[0]) * (1.0 - xv[1]) * (1.0 - xv[2]);
  return yv;
}


/* *****************************************************************
* Unused deformations
***************************************************************** */
inline AmanziGeometry::Point
Unused(double t, const AmanziGeometry::Point& xv)
{
  AmanziGeometry::Point yv(2);
  yv[0] = xv[0] * xv[1] + (1.0 - xv[1]) * std::pow(xv[0], 0.8);
  yv[1] = xv[1] * xv[0] + (1.0 - xv[0]) * std::pow(xv[1], 0.8);
  return yv;
}


/* *****************************************************************
* Sine-type
***************************************************************** */
inline AmanziGeometry::Point
SineProduct(double t, const AmanziGeometry::Point& xv)
{
  int d = xv.dim();
  double phi = 2 * M_PI;

  AmanziGeometry::Point yv(xv);
  double tmp = t * 0.1;
  for (int i = 0; i < d; ++i) tmp *= std::sin(xv[i] * phi);

  for (int i = 0; i < d; ++i) yv[i] = xv[i] + tmp;
  return yv;
}

} // namespace Amanzi

#endif
