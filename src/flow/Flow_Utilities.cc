/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <string>
#include <vector>

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "DenseMatrix.hh"
#include "Mesh.hh"
#include "tensor.hh"

#include "Flow_PK.hh"

namespace Amanzi {
namespace Flow {

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
  const Epetra_MultiVector& flux = *S_->GetFieldData("darcy_flux")->ViewComponent("face", true);
  
  int d(dim);
  AmanziGeometry::Point local_velocity(d);

  Teuchos::LAPACK<int, double> lapack;
  Teuchos::SerialDenseMatrix<int, double> matrix(d, d);
  double rhs[d];

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    for (int i = 0; i < d; i++) rhs[i] = 0.0;
    matrix.putScalar(0.0);

    for (int n = 0; n < nfaces; n++) {  // populate least-square matrix
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      double area = mesh_->face_area(f);

      for (int i = 0; i < d; i++) {
        rhs[i] += normal[i] * flux[0][f];
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
          N.AddRow(i, mesh_->face_normal(f));
          tmp[i] = flux[0][f];
        }
        N.Inverse();
        local_velocity += N * tmp;
      }
      local_velocity /= static_cast<double>(ncells);
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
void Flow_PK::CalculatePoreVelocity(
    std::vector<AmanziGeometry::Point>& xyz, 
    std::vector<AmanziGeometry::Point>& velocity,
    std::vector<double>& porosity, std::vector<double>& saturation)
{
  S_->GetFieldData("porosity")->ScatterMasterToGhosted();
  S_->GetFieldData("water_saturation")->ScatterMasterToGhosted();

  const Epetra_MultiVector& flux = *(S_->GetFieldData("darcy_flux")->ViewComponent("face", true));
  const Epetra_MultiVector& phi = *(S_->GetFieldData("porosity")->ViewComponent("cell", true));
  const Epetra_MultiVector& ws = *(S_->GetFieldData("water_saturation")->ViewComponent("cell", true));

  CalculateDarcyVelocity(xyz, velocity);

  // set markers for boundary nodes and faces
  int nnodes_owned = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);
  int nnodes_wghost = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::USED);

  std::vector<int> node_marker(nnodes_wghost);  
  AmanziMesh::Entity_ID_List nodes, cells;

  for (int f = 0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    if (ncells == 1) {
      mesh_->face_get_nodes(f, &nodes);
      int nnodes = nodes.size();

      for (int n = 0; n < nnodes; n++) node_marker[nodes[n]] = 1;
    }
  }

  // STEP 1: populate porosity and saturations in cell centers
  porosity.clear();
  saturation.clear();

  for (int c = 0; c < ncells_owned; c++) {
    porosity.push_back(phi[0][c]);
    saturation.push_back(ws[0][c]);
  }
    
  // STEP 2: recover porosity and saturation at boundary nodes
  double local_phi, local_ws;

  for (int v = 0; v < nnodes_owned; v++) {
    if (node_marker[v] > 0) {
      mesh_->node_get_cells(v, AmanziMesh::USED, &cells);
      int ncells = cells.size();

      local_phi = 0.0;
      local_ws = 0.0;
      for (int n = 0; n < ncells; n++) {
        int c = cells[n];
        local_phi += phi[0][c];
        local_ws += ws[0][c];
      }
      local_phi /= static_cast<double>(ncells);
      porosity.push_back(local_phi);

      local_ws /= static_cast<double>(ncells);
      saturation.push_back(local_ws);
    }
  }

  int n = xyz.size();
  for (int i = 0; i < n; i++) velocity[i] /= porosity[i] * saturation[i];  
}


/* ******************************************************************
* Write walkabout data
****************************************************************** */
void Flow_PK::WriteWalkabout(const Teuchos::Ptr<Checkpoint>& wlk)
{
  if ( !wlk->is_disabled() ) {

    wlk->CreateFile(S_->cycle());

    std::vector<AmanziGeometry::Point> xyz;
    std::vector<AmanziGeometry::Point> velocity;
    std::vector<double> porosity;
    std::vector<double> saturation;
    CalculatePoreVelocity(xyz, velocity, porosity, saturation);
    
    int n_loc = xyz.size();
    int n_glob;

    // create an epetra block map that we can use to create an appropriate
    // epetra multi vector
    const AmanziMesh::Mesh& mesh = *S_->GetMesh();
    mesh.get_comm()->SumAll(&n_loc, &n_glob, 1);
    
    Epetra_BlockMap map(n_glob, n_loc, 1, 0, *mesh.get_comm());
    
    int dim = mesh.space_dimension();
    // create an auxiliary vector that will hold the centrod and velocity, this is a cell based vector
    Teuchos::RCP<Epetra_MultiVector> aux =
        Teuchos::rcp(new Epetra_MultiVector(map, dim));
    
    std::vector<AmanziGeometry::Point>::const_iterator it;
    int i;
    for (it = xyz.begin(), i = 0; it != xyz.end(); ++it, ++i) {      
      (*(*aux)(0))[i]  = (*it)[0];
      if (dim > 1) (*(*aux)(1))[i]  = (*it)[1];
      if (dim > 2) (*(*aux)(2))[i]  = (*it)[2];
    }
    std::vector<std::string>  name;
    name.resize(0);
    name.push_back("x");
    if (dim > 1) name.push_back("y");
    if (dim > 2) name.push_back("z");
    wlk->WriteVector(*aux, name);
       

    for (it = velocity.begin(), i = 0; it != velocity.end(); ++it, ++i) {
      (*(*aux)(0))[i]  = (*it)[0];
      if (dim > 1) (*(*aux)(1))[i]  = (*it)[1];
      if (dim > 2) (*(*aux)(2))[i]  = (*it)[2];
    }
    name.resize(0);
    name.push_back("pore velocity x");
    if (dim > 1) name.push_back("pore velocity y");
    if (dim > 2) name.push_back("pore velocity z");
    wlk->WriteVector(*aux, name);
    
    // reallocate a new aux vector only if we need to 
    if (dim != 2) aux = Teuchos::rcp(new Epetra_MultiVector(map, 2));
    
    std::vector<double>::const_iterator it0, it1;
    
    for (it0 = saturation.begin(), it1 = porosity.begin(), i = 0; it0 != saturation.end(); ++it0, ++it1, ++i) {    
      (*(*aux)(0))[i]  = *it0;
      (*(*aux)(1))[i]  = *it1;
    }
    name.resize(0);
    name.push_back("saturation");    
    name.push_back("porosity");
    wlk->WriteVector(*aux, name);

    wlk->WriteAttributes(S_->time(), S_->cycle());    
    wlk->Finalize();
  }
}



}  // namespace Flow
}  // namespace Amanzi

