/*
  MPC

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Checkpointing Walkabout data.
*/

#include <iomanip>
#include <iostream>

// TPLs
#include "Epetra_MpiComm.h"
#include "Epetra_Export.h"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

// Amanzi
#include "DenseMatrix.hh"
#include "errors.hh"
#include "Darcy_PK.hh"
#include "Mesh.hh"
#include "OperatorDefs.hh"
#include "ParallelCommunication.hh"
#include "ReconstructionCell.hh"
#include "Tensor.hh"
#include "WhetStoneMeshUtils.hh"

#include "WalkaboutCheckpoint.hh"

namespace Amanzi {

/* ******************************************************************
* Calculate full vectors of Darcy velocities. The velocity is 
* evaluated at mesh nodes.
****************************************************************** */
void WalkaboutCheckpoint::CalculateDarcyVelocity(
    Teuchos::RCP<State>& S,
    std::vector<AmanziGeometry::Point>& xyz, 
    std::vector<AmanziGeometry::Point>& velocity) const
{
  xyz.clear();
  velocity.clear();
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = S->GetMesh();

  int nnodes_owned  = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);
  int nnodes_wghost = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  double rho = *S->GetScalarData("const_fluid_density");
  S->GetFieldData("darcy_flux")->ScatterMasterToGhosted();
  const Epetra_MultiVector& flux = *S->GetFieldData("darcy_flux")->ViewComponent("face", true);
  
  int d = mesh->space_dimension();
  AmanziGeometry::Point node_velocity(d);

  // optional flux constraint at nodes
  // -- flag
  bool projection(false);
  if (pk_ != Teuchos::null && pk_->name() == "flow") projection = true;

  // least-square recovery at mesh nodes 
  AmanziMesh::Entity_ID_List faces;
  AmanziGeometry::Point xv(d);
  WhetStone::DenseVector rhs(d), sol(d);
  WhetStone::DenseMatrix matrix(d, d);

  for (int v = 0; v < nnodes_owned; ++v) {
    mesh->node_get_faces(v, AmanziMesh::Parallel_type::ALL, &faces);
    int nfaces = faces.size();

    rhs.PutScalar(0.0);
    matrix.PutScalar(0.0);

    for (int n = 0; n < nfaces; n++) {  // populate least-square matrix
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh->face_normal(f);
      double area = mesh->face_area(f);

      for (int i = 0; i < d; i++) {
        rhs(i) += normal[i] * flux[0][f];
        matrix(i, i) += normal[i] * normal[i];
        for (int j = i+1; j < d; j++) {
          matrix(j, i) = matrix(i, j) += normal[i] * normal[j];
        }
      }
    }

    matrix.Inverse();
    matrix.Multiply(rhs, sol, false);
    
    // enforce constraint: formulas follow from solution of saddle-point problem
    if (projection) {
      const auto pk_flow = dynamic_cast<Flow::Flow_PK*>(pk_.get());
      const auto& bc = pk_flow->op_bc();
      const auto& bc_model = bc->bc_model();
      const auto& bc_value = bc->bc_value();

      // -- sort normals by angle
      int dir;
      std::vector<AmanziGeometry::Point> node_basis;
      std::vector<double> node_area, node_flux;

      for (int n = 0; n < nfaces; ++n) {
        int f = faces[n];
        if (bc_model[f] == Operators::OPERATOR_BC_NEUMANN) {
          double area = mesh->face_area(f);
          AmanziGeometry::Point normal = WhetStone::face_normal_exterior(*mesh, f, &dir);

          bool flag(false);
          for (int m = 0; m < node_basis.size(); ++m) {
            double angle = (node_basis[m] * normal) / norm(node_basis[m]) / area;
            if (angle > 0.9) {
              node_basis[m] += normal;
              node_area[m] += area;
              node_flux[m] += bc_value[f] / rho * area;
              flag = true;
              break;
            }
          }
          if (! flag) {
            node_basis.push_back(normal);
            node_area.push_back(area);
            node_flux.push_back(bc_value[f] / rho * area);
          }
        } 
      }

      int nbasis = std::min((int)node_basis.size(), d);
      if (nbasis > 0) { 
        WhetStone::DenseMatrix N(d, nbasis), MN(d, nbasis), sigma(nbasis, nbasis);
        WhetStone::DenseVector f(nbasis), v1(nbasis), v2(nbasis), v3(d);

        for (int m = 0; m < nbasis; ++m) {
          node_basis[m] /= node_area[m];
          node_flux[m] /= node_area[m];

          for (int k = 0; k < d; ++k) {
            N(k, m) = node_basis[m][k];
          }
          f(m) = node_flux[m];
        }

        MN.Multiply(matrix, N, false);
        sigma.Multiply(N, MN, true);

        N.Multiply(sol, v1, true);
        v1 -= f;
          
        sigma.Inverse();
        sigma.Multiply(v1, v2, false);

        MN.Multiply(v2, v3, false);
        sol -= v3;
      }
    }

    for (int i = 0; i < d; i++) node_velocity[i] = sol(i);
    velocity.push_back(node_velocity);

    mesh->node_get_coordinates(v, &xv);
    xyz.push_back(xv);
  }
}


/* ******************************************************************
* Calculating an extended vector of Darcy velocities. The velocity
* is evaluated at cell-center and at boundary points.
****************************************************************** */
void WalkaboutCheckpoint::CalculateData(
    Teuchos::RCP<State>& S,
    std::vector<AmanziGeometry::Point>& xyz, 
    std::vector<AmanziGeometry::Point>& velocity,
    std::vector<double>& porosity, std::vector<double>& saturation,
    std::vector<double>& pressure, std::vector<double>& isotherm_kd,
    std::vector<int>& material_ids)
{
  const auto& mesh = S->GetMesh();
  S->GetFieldData("porosity")->ScatterMasterToGhosted();
  S->GetFieldData("saturation_liquid")->ScatterMasterToGhosted();
  S->GetFieldData("pressure")->ScatterMasterToGhosted();

  const Epetra_MultiVector& flux = *S->GetFieldData("darcy_flux")->ViewComponent("face", true);
  const Epetra_MultiVector& phi = *S->GetFieldData("porosity")->ViewComponent("cell", true);
  const Epetra_MultiVector& ws = *S->GetFieldData("saturation_liquid")->ViewComponent("cell", true);
  auto p = S->GetFieldData("pressure")->ViewComponent("cell", true);

  int nnodes_owned = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);
  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  // process non-flow state variables
  bool flag(false);
  Teuchos::RCP<const Epetra_MultiVector> kd;
  if (S->HasField("isotherm_kd")) {
    flag = true;
    S->GetFieldData("isotherm_kd")->ScatterMasterToGhosted();
    kd = S->GetFieldData("isotherm_kd")->ViewComponent("cell", true);
  }

  CalculateDarcyVelocity(S, xyz, velocity);

  // collect material information
  AmanziMesh::Entity_ID_List cells;
  Epetra_IntVector cell_ids(mesh->cell_map(true));

  cell_ids.PutValue(-1);

  if (plist_.isSublist("write regions")) {
    const Teuchos::ParameterList& tmp = plist_.sublist("write regions");
    std::vector<std::string> regs = tmp.get<Teuchos::Array<std::string> >("region names").toVector();
    std::vector<int> ids = tmp.get<Teuchos::Array<int> >("material ids").toVector();

    for (int n = 0; n < regs.size(); ++n) {
      mesh->get_set_entities(regs[n], AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &cells);

      for (auto it = cells.begin(); it != cells.end(); ++it) {
        cell_ids[*it] = ids[n];
      }
    }
  }

  // verify material ids
  if (cell_ids.MinValue() == -1) {
    Errors::Message msg("Negative material id: check materials attribute \"id\"");
    Exceptions::amanzi_throw(msg);
  } 

  ParallelCommunication pp(mesh);
  pp.CopyMasterCell2GhostCell(cell_ids);

  // prepare reconstrction data
  Teuchos::ParameterList plist;
  plist.set<std::string>("limiter", "tensorial");

  Operators::ReconstructionCell lifting(mesh);
  lifting.Init(p, plist);
  lifting.ComputeGradient();

  // Populate state data at mesh nodes
  porosity.clear();
  saturation.clear();
  pressure.clear();
  isotherm_kd.clear();
  material_ids.clear();

  int dim = mesh->space_dimension();
  AmanziGeometry::Point xv(dim);

  double local_phi, local_ws, local_p, local_kd;
  int local_id;
  for (int v = 0; v < nnodes_owned; v++) {
    mesh->node_get_coordinates(v, &xv);
    mesh->node_get_cells(v, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();

    local_id = -1;
    local_phi = 0.0;
    local_ws = 0.0;
    local_p = 0.0;
    local_kd = 0.0;
    for (int n = 0; n < ncells; n++) {
      int c = cells[n];
      local_phi += phi[0][c];
      local_ws += ws[0][c];
      local_p += lifting.getValue(c, xv);
      if (flag) local_kd += (*kd)[0][c];
      local_id = std::max(local_id, cell_ids[c]);
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

    material_ids.push_back(local_id);
  }

  int n = xyz.size();
  for (int i = 0; i < n; i++) velocity[i] /= porosity[i] * saturation[i];  
}


/* ******************************************************************
* Write walkabout data
****************************************************************** */
void WalkaboutCheckpoint::WriteDataFile(
    Teuchos::RCP<State>& S, Teuchos::RCP<PK> pk)
{
  if (!is_disabled()) {
    pk_ = pk;
    CreateFile(S->cycle());

    std::vector<AmanziGeometry::Point> xyz;
    std::vector<AmanziGeometry::Point> velocity;
    std::vector<double> porosity, saturation;
    std::vector<double> pressure, isotherm_kd;
    std::vector<int> material_ids;

    CalculateData(S, xyz, velocity, porosity,
                  saturation, pressure, isotherm_kd, material_ids);
    
    // create a block map that we can use to create a new Epetra_MultiVector
    const AmanziMesh::Mesh& mesh = *S->GetMesh();
    const auto& map = mesh.node_map(false); 

    const auto& source_map = mesh.node_map(true); 
    Epetra_Map target_map(-1, map.NumMyElements(), 0, map.Comm());
    Epetra_Export exporter(source_map, target_map);
    
    // create an auxiliary vector that will hold the centroid and velocity
    int dim = mesh.space_dimension();
    Teuchos::RCP<Epetra_MultiVector> vs = Teuchos::rcp(new Epetra_MultiVector(source_map, dim));
    Teuchos::RCP<Epetra_MultiVector> vt = Teuchos::rcp(new Epetra_MultiVector(target_map, dim));
    
    // redistribute data
    // -- coordinates of mesh nodes
    int ndata = xyz.size();
    for (int n = 0; n < ndata; n++) {      
      for (int i = 0; i < dim; i++) {
        (*(*vs)(i))[n] = xyz[n][i];
      }
    }
    std::vector<std::string> name;
    name.resize(0);
    name.push_back("x");
    name.push_back("y");
    if (dim > 2) name.push_back("z");

    vt->Export(*vs, exporter, Add);
    WriteVector(*vt, name);
       
    // -- pore velocity
    vs->PutScalar(0.0);
    vt->PutScalar(0.0);
    for (int n = 0; n < ndata; n++) {
      for (int i = 0; i < dim; i++) {
        (*(*vs)(i))[n] = velocity[n][i];
      }
    }
    name.resize(0);
    name.push_back("pore velocity x");
    name.push_back("pore velocity y");
    if (dim > 2) name.push_back("pore velocity z");

    vt->Export(*vs, exporter, Add);
    WriteVector(*vt, name);
    
    // -- saturation, porosity, pressure
    //    reallocate temporary vectors only if we need to
    if (dim != 3) { 
      vs = Teuchos::rcp(new Epetra_MultiVector(source_map, 3));
      vt = Teuchos::rcp(new Epetra_MultiVector(target_map, 3));
    }
    vs->PutScalar(0.0);
    vt->PutScalar(0.0);
    
    for (int n = 0; n < ndata; n++) {
      (*(*vs)(0))[n] = saturation[n];
      (*(*vs)(1))[n] = porosity[n];
      (*(*vs)(2))[n] = pressure[n];
    }
    name.resize(0);
    name.push_back("saturation");    
    name.push_back("porosity");
    name.push_back("pressure");

    vt->Export(*vs, exporter, Add);
    WriteVector(*vt, name);

    // -- bulk density and isotherm kd
    if (isotherm_kd.size() > 0) {
      vs = Teuchos::rcp(new Epetra_MultiVector(source_map, 1));
      vt = Teuchos::rcp(new Epetra_MultiVector(target_map, 1));

      for (int n = 0; n < ndata; n++) {    
        (*(*vs)(0))[n] = isotherm_kd[n];
      }

      name.resize(0);
      name.push_back("isotherm kd");
      vt->Export(*vs, exporter, Add);
      WriteVector(*vt, name);
    } 

    // -- material ids
    vs = Teuchos::rcp(new Epetra_MultiVector(source_map, 1));
    vt = Teuchos::rcp(new Epetra_MultiVector(target_map, 1));

    for (int n = 0; n < ndata; n++) {    
      (*(*vs)(0))[n] = material_ids[n];
    }

    name.resize(0);
    name.push_back("material ids");

    vt->Export(*vs, exporter, Add);
    WriteVector(*vt, name);

    // dump pairs: material ids, material name
    if (plist_.isSublist("write regions")) {
      const Teuchos::ParameterList& tmp = plist_.sublist("write regions");
      std::vector<std::string> names = tmp.get<Teuchos::Array<std::string> >("material names").toVector();
      std::vector<int> ids = tmp.get<Teuchos::Array<int> >("material ids").toVector();
 
      int nnames = names.size();

      if (nnames > 0) {
        int *tmp_ids; 
        char **tmp_names;

        tmp_ids = (int*)malloc(nnames * sizeof(int));
        tmp_names = (char**)malloc(nnames * sizeof(char*));

        for (int i = 0; i < nnames; ++i) {
          tmp_ids[i] = ids[i];
          tmp_names[i] = (char*)malloc((names[i].size() + 1) * sizeof(char));
          strcpy(tmp_names[i], names[i].c_str());
        }

        checkpoint_output_->writeAttrInt(tmp_ids, nnames, "material_labels");
        checkpoint_output_->writeDataString(tmp_names, nnames, "material_names");

        for (int i = 0; i < nnames; ++i) free(tmp_names[i]);
        free(tmp_names);
        free(tmp_ids);
      }
    }

    // timestamp and cycle number 
    WriteAttributes(S->time(), S->cycle());    
    Finalize();
  }
}

} // namespace Amanzi
