/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Collection of mesh deformation tools.
*/

#ifndef AMANZI_OPERATOR_DEFORM_MESH_HH_
#define AMANZI_OPERATOR_DEFORM_MESH_HH_

#include "Mesh.hh"
#include "CompositeVector.hh"

#include "OperatorDefs.hh"

namespace Amanzi{

// collection of routines for mesh deformation
void DeformMesh(const Teuchos::RCP<AmanziMesh::Mesh>& mesh1, int deform, double t,
                const Teuchos::RCP<const AmanziMesh::Mesh>& mesh0 = Teuchos::null);

AmanziGeometry::Point TaylorGreenVortex(double t, const AmanziGeometry::Point& xv);
AmanziGeometry::Point Rotation2D(double t, const AmanziGeometry::Point& xv);
AmanziGeometry::Point CompressionExpansion2D(double t, const AmanziGeometry::Point& xv);
AmanziGeometry::Point CompressionExpansion3D(double t, const AmanziGeometry::Point& xv);
AmanziGeometry::Point Unused(double t, const AmanziGeometry::Point& xv);


/* *****************************************************************
* Deform mesh1 using coordinates given by mesh0
***************************************************************** */
inline
void DeformMesh(const Teuchos::RCP<AmanziMesh::Mesh>& mesh1, int deform, double t,
                const Teuchos::RCP<const AmanziMesh::Mesh>& mesh0)
{
  // create distributed random vector
  int d = mesh1->space_dimension();
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh1)->SetGhosted(true)->AddComponent("node", AmanziMesh::NODE, d);
  CompositeVector random(cvs);

  int gid = mesh1->node_map(false).MaxAllGID();
  double scale = 0.2 * std::pow(gid, -1.0 / d);
  Epetra_MultiVector& random_n = *random.ViewComponent("node", true);

  random_n.Random();
  random_n.Scale(scale);
  random.ScatterMasterToGhosted();

  // relocate mesh nodes
  AmanziGeometry::Point xv(d), yv(d), uv(d);
  AmanziMesh::Entity_ID_List nodeids;
  AmanziGeometry::Point_List new_positions, final_positions;

  int nnodes = mesh1->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);

  for (int v = 0; v < nnodes; ++v) {
    if (mesh0.get()) 
      mesh0->node_get_coordinates(v, &xv);
    else
      mesh1->node_get_coordinates(v, &xv);
      
    nodeids.push_back(v);

    if (deform == 1) {
      yv = TaylorGreenVortex(t, xv);
    } else if (deform == 2) {
      yv = Unused(t, xv);
    } else if (deform == 4) {
      yv = CompressionExpansion2D(t, xv);
    } else if (deform == 5) {
      yv = CompressionExpansion3D(t, xv);
    } else if (deform == 6) {
      yv = Rotation2D(t, xv);
    } else if (deform == 7) {
      for (int i = 0; i < d; ++i) yv[i] = random_n[i][v];
    }

    new_positions.push_back(yv);
  }
  mesh1->deform(nodeids, new_positions, false, &final_positions);
}


/* *****************************************************************
* Taylor-Green vortex
***************************************************************** */
inline
AmanziGeometry::Point TaylorGreenVortex(double t, const AmanziGeometry::Point& xv)
{
  int d = xv.dim();
  AmanziGeometry::Point uv(d), yv(xv);

  double ds(0.0001);
  int n = t / ds;
  for (int i = 0; i < n; ++i) {
    if (d == 2) {
      uv[0] = 0.2 * std::sin(M_PI * yv[0]) * std::cos(M_PI * yv[1]);
      uv[1] =-0.2 * std::cos(M_PI * yv[0]) * std::sin(M_PI * yv[1]);
    } else {
      uv[0] = 0.2 * std::sin(M_PI * yv[0]) * std::cos(M_PI * yv[1]) * std::cos(M_PI * yv[2]);
      uv[1] =-0.1 * std::cos(M_PI * yv[0]) * std::sin(M_PI * yv[1]) * std::cos(M_PI * yv[2]);
      uv[2] =-0.1 * std::cos(M_PI * yv[0]) * std::cos(M_PI * yv[1]) * std::sin(M_PI * yv[2]);
    }
    yv += uv * ds;
  }

  return yv;
}


/* *****************************************************************
* Rotation
***************************************************************** */
inline
AmanziGeometry::Point Rotation2D(double t, const AmanziGeometry::Point& xv)
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
inline
AmanziGeometry::Point CompressionExpansion2D(double t, const AmanziGeometry::Point& xv)
{
  AmanziGeometry::Point yv(xv);
  yv[0] += t * xv[0] * xv[1] * (1.0 - xv[0]) / 2;
  yv[1] += t * xv[0] * xv[1] * (1.0 - xv[1]) / 2;
  return yv;
}


inline
AmanziGeometry::Point CompressionExpansion3D(double t, const AmanziGeometry::Point& xv)
{
  AmanziGeometry::Point yv(xv);
  yv[0] += t * xv[0] * xv[1] * xv[2] * (1.0 - xv[0]) / 2;
  yv[1] += t * xv[0] * xv[1] * xv[2] * (1.0 - xv[1]) / 2;
  yv[2] += t * xv[0] * xv[1] * xv[2] * (1.0 - xv[2]) / 2;
  return yv;
}


/* *****************************************************************
* Unused deformations
***************************************************************** */
inline
AmanziGeometry::Point Unused(double t, const AmanziGeometry::Point& xv)
{
  AmanziGeometry::Point yv(2);
  yv[0] = xv[0] * xv[1] + (1.0 - xv[1]) * std::pow(xv[0], 0.8);
  yv[1] = xv[1] * xv[0] + (1.0 - xv[0]) * std::pow(xv[1], 0.8);
  return yv;
}

}  // namespace Amanzi

#endif
