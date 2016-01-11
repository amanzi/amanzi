/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Amanzi

License:


Walkabout

------------------------------------------------------------------------- */

#include "Walkabout_observations.hh"
#include "Epetra_MpiComm.h"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include <iostream>
#include <iomanip>

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "DenseMatrix.hh"
#include "Mesh.hh"
#include "Tensor.hh"

namespace Amanzi {

/* ******************************************************************
* Calculating an extended vector of Darcy velocities. The velocity
* is evaluated at cell-center and at boundary points.
****************************************************************** */
void Walkabout_observations::CalculateDarcyVelocity(Teuchos::RCP<State>& S,
				     std::vector<AmanziGeometry::Point>& xyz, 
                                     std::vector<AmanziGeometry::Point>& velocity)
{
  xyz.clear();
  velocity.clear();
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = S->GetMesh();

  int nnodes_owned  = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);
  int nnodes_wghost = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::USED);
  int nfaces_owned  = mesh->num_entities(AmanziMesh::FACE,  AmanziMesh::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE,  AmanziMesh::USED);
  int ncells_owned  = mesh->num_entities(AmanziMesh::CELL,  AmanziMesh::OWNED);

  // set markers for boundary nodes and faces
  std::vector<int> node_marker(nnodes_wghost);  
  std::vector<int> face_marker(nfaces_wghost);  

  AmanziMesh::Entity_ID_List nodes, faces, cells;

  for (int f = 0; f < nfaces_owned; f++) {
    mesh->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    if (ncells == 1) {
      face_marker[f] = 1;

      mesh->face_get_nodes(f, &nodes);
      int nnodes = nodes.size();

      for (int n = 0; n < nnodes; n++) node_marker[nodes[n]] = 1;
    }
  }

  // STEP 1: recover velocity in cell centers
  const Epetra_MultiVector& flux = *S->GetFieldData("darcy_flux")->ViewComponent("face", true);
  
  int dim = mesh->space_dimension();
  int d(dim);
  AmanziGeometry::Point local_velocity(d);

  Teuchos::LAPACK<int, double> lapack;
  Teuchos::SerialDenseMatrix<int, double> matrix(d, d);
  double rhs[d];

  for (int c = 0; c < ncells_owned; c++) {
    mesh->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    for (int i = 0; i < d; i++) rhs[i] = 0.0;
    matrix.putScalar(0.0);

    for (int n = 0; n < nfaces; n++) {  // populate least-square matrix
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh->face_normal(f);
      double area = mesh->face_area(f);

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

    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    xyz.push_back(xc);
  }

  // STEP 2: recover velocity at boundary nodes
  WhetStone::Tensor N(d, 2);
  AmanziGeometry::Point tmp(d);

  for (int v = 0; v < nnodes_owned; v++) {
    if (node_marker[v] > 0) {
      mesh->node_get_cells(v, AmanziMesh::USED, &cells);
      int ncells = cells.size();

      local_velocity.set(0.0);
      for (int n = 0; n < ncells; n++) {
        int c = cells[n];
        mesh->node_get_cell_faces(v, c, AmanziMesh::USED, &faces);
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
          N.SetRow(i, mesh->face_normal(f));
          tmp[i] = flux[0][f];
        }
        N.Inverse();
        local_velocity += N * tmp;
      }
      local_velocity /= static_cast<double>(ncells);
      velocity.push_back(local_velocity);

      mesh->node_get_coordinates(v, &tmp);
      xyz.push_back(tmp);
    }
  }
}


/* ******************************************************************
* Calculating an extended vector of Darcy velocities. The velocity
* is evaluated at cell-center and at boundary points.
****************************************************************** */
void Walkabout_observations::CalculatePoreVelocity( Teuchos::RCP<State>& S,
						    std::vector<AmanziGeometry::Point>& xyz, 
						    std::vector<AmanziGeometry::Point>& velocity,
						    std::vector<double>& porosity, std::vector<double>& saturation,
						    std::vector<double>& pressure, std::vector<double>& isotherm_kd)
{
  const AmanziMesh::Mesh& mesh = *S->GetMesh();
  S->GetFieldData("porosity")->ScatterMasterToGhosted();
  S->GetFieldData("saturation_liquid")->ScatterMasterToGhosted();

  const Epetra_MultiVector& flux = *S->GetFieldData("darcy_flux")->ViewComponent("face", true);
  const Epetra_MultiVector& phi = *S->GetFieldData("porosity")->ViewComponent("cell", true);
  const Epetra_MultiVector& ws = *S->GetFieldData("saturation_liquid")->ViewComponent("cell", true);
  const Epetra_MultiVector& p = *S->GetFieldData("pressure")->ViewComponent("cell", true);

  // process non-flow state variables
  bool flag(false);
  Teuchos::RCP<const Epetra_MultiVector> kd;
  if (S->HasField("isotherm_kd")) {
    flag = true;
    S->GetFieldData("isotherm_kd")->ScatterMasterToGhosted();
    kd = S->GetFieldData("isotherm_kd")->ViewComponent("cell", true);
  }

  CalculateDarcyVelocity(S, xyz, velocity);

  // set markers for boundary nodes and faces
  int nnodes_owned = mesh.num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);
  int nnodes_wghost = mesh.num_entities(AmanziMesh::NODE, AmanziMesh::USED);
  int nfaces_owned  = mesh.num_entities(AmanziMesh::FACE,  AmanziMesh::OWNED);
  int ncells_owned  = mesh.num_entities(AmanziMesh::CELL,  AmanziMesh::OWNED);

  std::vector<int> node_marker(nnodes_wghost);  
  AmanziMesh::Entity_ID_List nodes, cells;

  for (int f = 0; f < nfaces_owned; f++) {
    mesh.face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    if (ncells == 1) {
      mesh.face_get_nodes(f, &nodes);
      int nnodes = nodes.size();

      for (int n = 0; n < nnodes; n++) node_marker[nodes[n]] = 1;
    }
  }

  // STEP 1: populate porosity and saturations in cell centers
  porosity.clear();
  saturation.clear();
  pressure.clear();
  isotherm_kd.clear();

  for (int c = 0; c < ncells_owned; c++) {
    porosity.push_back(phi[0][c]);
    saturation.push_back(ws[0][c]);
    pressure.push_back(p[0][c]);
    if (flag) isotherm_kd.push_back((*kd)[0][c]);
  }
    
  // STEP 2: recover porosity and saturation at boundary nodes
  double local_phi, local_ws, local_p, local_kd;

  for (int v = 0; v < nnodes_owned; v++) {
    if (node_marker[v] > 0) {
      mesh.node_get_cells(v, AmanziMesh::USED, &cells);
      int ncells = cells.size();

      local_phi = 0.0;
      local_ws = 0.0;
      local_p = 0.0;
      local_kd = 0.0;
      for (int n = 0; n < ncells; n++) {
        int c = cells[n];
        local_phi += phi[0][c];
        local_ws += ws[0][c];
        local_p += p[0][c];
        if (flag) local_kd += (*kd)[0][c];
      }
      local_phi /= ncells;
      porosity.push_back(local_phi);

      local_ws /= ncells;
      saturation.push_back(local_ws);

      local_p /= ncells;
      pressure.push_back(local_p);

      if (flag) {
        local_kd /= ncells;
        isotherm_kd.push_back(local_kd);
      }
    }
  }

  int n = xyz.size();
  for (int i = 0; i < n; i++) velocity[i] /= porosity[i] * saturation[i];  
}


/* ******************************************************************
* Write walkabout data
****************************************************************** */
void Walkabout_observations::WriteWalkabout(Teuchos::RCP<State>& S)
{
  if (!is_disabled()) {
    CreateFile(S->cycle());

    std::vector<AmanziGeometry::Point> xyz;
    std::vector<AmanziGeometry::Point> velocity;
    std::vector<double> porosity, saturation;
    std::vector<double> pressure, isotherm_kd;
    CalculatePoreVelocity(S, xyz, velocity, porosity, saturation, pressure, isotherm_kd);
    
    int n_loc = xyz.size();
    int n_glob;

    // create a block map that we can use to create an Epetra_MultiVector
    const AmanziMesh::Mesh& mesh = *S->GetMesh();
    mesh.get_comm()->SumAll(&n_loc, &n_glob, 1);
    
    Epetra_BlockMap map(n_glob, n_loc, 1, 0, *mesh.get_comm());
    
    // create an auxiliary vector that will hold the centroid and velocity
    int dim = mesh.space_dimension();
    Teuchos::RCP<Epetra_MultiVector> aux = Teuchos::rcp(new Epetra_MultiVector(map, dim));
    
    int ndata = xyz.size();
    for (int n = 0; n < ndata; n++) {      
      for (int i = 0; i < dim; i++) {
        (*(*aux)(i))[n] = xyz[n][i];
      }
    }
    std::vector<std::string> name;
    name.resize(0);
    name.push_back("x");
    name.push_back("y");
    if (dim > 2) name.push_back("z");
    WriteVector(*aux, name);
       
    for (int n = 0; n < ndata; n++) {
      for (int i = 0; i < dim; i++) {
        (*(*aux)(i))[n] = velocity[n][i];
      }
    }
    name.resize(0);
    name.push_back("pore velocity x");
    name.push_back("pore velocity y");
    if (dim > 2) name.push_back("pore velocity z");
    WriteVector(*aux, name);
    
    // reallocate a new aux vector only if we need to 
    if (dim != 3) aux = Teuchos::rcp(new Epetra_MultiVector(map, 3));
    
    for (int n = 0; n < ndata; n++) {
      (*(*aux)(0))[n] = saturation[n];
      (*(*aux)(1))[n] = porosity[n];
      (*(*aux)(2))[n] = pressure[n];
    }
    name.resize(0);
    name.push_back("saturation");    
    name.push_back("porosity");
    name.push_back("pressure");
    WriteVector(*aux, name);

    // dump other parameters: "bulk density" and "isotherm kd"
    if (isotherm_kd.size() > 0) {
      aux = Teuchos::rcp(new Epetra_MultiVector(map, 1));
      // aux = Teuchos::rcp(new Epetra_MultiVector(map, 2));

      for (int n = 0; n < ndata; n++) {    
        (*(*aux)(0))[n] = isotherm_kd[n];
        // (*(*aux)(1))[n] =  particle_density[n] * (1.0 - porosity[n]);
      }

      name.resize(0);
      name.push_back("isotherm kd");
      // name.push_back("bulk density");    
      WriteVector(*aux, name);
    } 

    WriteAttributes(S->time(), S->cycle());    
    Finalize();
  }
}


} // namespace
