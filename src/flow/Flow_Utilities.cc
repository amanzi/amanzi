/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include <string>
#include <vector>

#include "tensor.hh"

#include "Mesh.hh"
#include "Flow_PK.hh"


namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Calculating an extended vector of Darcy velocities. The velocity
* is evaluated at cell-center and at boundary points.
****************************************************************** */
void Flow_PK::CalculateDarcyVelocity(std::vector<AmanziGeometry::Point>& xyz, 
                                     std::vector<AmanziGeometry::Point>& velocity)
{
  xyz.clear();
  velocity.clear();

  int nnodes_owned = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);
  int nnodes_wghost = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::USED);

  // set markers for boundary nodes and faces
  std::vector<int> node_marker(nnodes_wghost);  
  std::vector<int> face_marker(nfaces_wghost);  

  AmanziMesh::Entity_ID_List nodes, faces, cells;
  std::vector<int> dirs;

  for (int f = 0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    if (ncells == 1) {
      face_marker[f] = 1;

      mesh_->face_get_nodes(f, &nodes);
      int nnodes = nodes.size();

      for (int n = 0; n < nnodes; n++) node_marker[nodes[n]] = 1;
    }
  }

  // STEP 1: recover velocity in cell centers
  Epetra_Vector& flux = FS->ref_darcy_flux();
  
#ifdef HAVE_MPI
  Epetra_Vector flux_wghost(mesh_->face_map(true));
  FS->CopyMasterFace2GhostFace(flux, flux_wghost);
#else
  Epetra_Vector& flux_wghost = flux;
#endif

  int d = mesh_->space_dimension();
  AmanziGeometry::Point local_velocity(d);

  Teuchos::LAPACK<int, double> lapack;
  Teuchos::SerialDenseMatrix<int, double> matrix(d, d);
  double rhs[d];

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int i = 0; i < d; i++) rhs[i] = 0.0;
    matrix.putScalar(0.0);

    for (int n = 0; n < nfaces; n++) {  // populate least-square matrix
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      double area = mesh_->face_area(f);

      for (int i = 0; i < d; i++) {
        rhs[i] += normal[i] * flux_wghost[f];
        matrix(i, i) += normal[i] * normal[i];
        for (int j = i+1; j < d; j++) {
          matrix(j, i) = matrix(i, j) += normal[i] * normal[j];
        }
      }
    }

    int info;
    lapack.POSV('U', d, 1, matrix.values(), d, rhs, d, &info);

    for (int i = 0; i < d; i++) local_velocity[i] = rhs[i];
    velocity.push_back(local_velocity);

    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
    xyz.push_back(xc);
  }

  // STEP 2: recover velocity at boundary nodes
  WhetStone::Tensor N(d, 2);
  AmanziGeometry::Point tmp(d);

  for (int v = 0; v < nnodes_owned; v++) {
    if (node_marker[v] > 0) {
      mesh_->node_get_cells(v, AmanziMesh::USED, &cells);
      int ncells = cells.size();

      local_velocity.set(0.0);
      for (int n = 0; n < ncells; n++) {
        int c = cells[n];
        mesh_->node_get_cell_faces(v, c, AmanziMesh::USED, &faces);
        int nfaces = faces.size();

        if (nfaces > d) {  // Move boundary face to the top of the list. 
          for (int i = d; i < nfaces; i++) {
            int f = faces[i];
            if (face_marker[f] > 0) {
              faces[0] = f;
              break;
            }
          }
        }

        for (int i = 0; i < d; i++) {
          int f = faces[i];
          N.add_row(i, mesh_->face_normal(f));
          tmp[i] = flux_wghost[f];
        }
        N.inverse();
        local_velocity += N * tmp;
      }
      local_velocity /= (double)d;
      velocity.push_back(local_velocity);

      mesh_->node_get_coordinates(v, &tmp);
      xyz.push_back(tmp);
    }
  }
}


/* ******************************************************************
* Calculating an extended vector of Darcy velocities. The velocity
* is evaluated at cell-center and at boundary points.
****************************************************************** */
void Flow_PK::CalculatePoreVelocity(std::vector<AmanziGeometry::Point>& xyz, 
                                    std::vector<AmanziGeometry::Point>& velocity)
{
  CalculateDarcyVelocity(xyz, velocity);

  Epetra_Vector& porosity = FS->ref_porosity();

  int n = xyz.size();
  for (int i = 0; i< n; i++) velocity[i] /= porosity[i];  
}

}  // namespace AmanziFlow
}  // namespace Amanzi

