/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (version 1) 
           Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"

#include "GMVMesh.hh"
#include "mfd3d.hh"

#include "Mesh.hh"
#include "State.hh"
#include "Flow_PK.hh"


namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Initiazition of fundamental flow sturctures.                                              
****************************************************************** */
void Flow_PK::Init()
{
  T_physics = dT = 0.0;

  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);

  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  nseepage_prev = 0;
  ti_phase_counter = 0;

  // Fundamental physical quantities
  double* gravity_data;
  S_->GetConstantVectorData("gravity")->ExtractView(&gravity_data);
  gravity_.init(dim);
  for (int k = 0; k < dim; k++) gravity_[k] = gravity_data[k];
  g_ = fabs(gravity_[dim - 1]);

  // Other constant (temporarily) physical quantaties
  rho_ = *(S_->GetScalarData("fluid_density"));
  mu_ = *(S_->GetScalarData("fluid_viscosity"));

  // parallel execution data
  MyPID = 0;
#ifdef HAVE_MPI
  MyPID = mesh_->cell_map(false).Comm().MyPID();
#endif
}


/* ******************************************************************
* Populate data needed by submodels.
****************************************************************** */
void Flow_PK::ProcessBCs()
{
  for (int f = 0; f < nfaces_owned; f++) {
    if (bc_submodel[f] & FLOW_BC_SUBMODEL_RAINFALL) {
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      rainfall_factor[f] = fabs(normal[dim - 1]) / norm(normal);
    }
  }
}


/* ******************************************************************
* Add a boundary marker to used faces.
* WARNING: we can skip update of ghost boundary faces, b/c they 
* should be always owned. 
****************************************************************** */
void Flow_PK::ComputeBCs(const CompositeVector& pressure)
{
  const Epetra_MultiVector& pressure_cells = *(pressure.ViewComponent("cell"));
  const Epetra_MultiVector& pressure_faces = *(pressure.ViewComponent("face"));
  
  int flag_essential_bc = 0;
  bc_tuple zero = {0.0, 0.0};
  for (int n = 0; n < bc_model.size(); n++) {
    bc_model[n] = FLOW_BC_FACE_NULL;
    bc_values[n] = zero;
  }

  Functions::FlowBoundaryFunction::Iterator bc;
  for (bc = bc_pressure->begin(); bc != bc_pressure->end(); ++bc) {
    int f = bc->first;
    bc_model[f] = FLOW_BC_FACE_PRESSURE;
    bc_values[f][0] = bc->second;
    flag_essential_bc = 1;
  }

  for (bc = bc_head->begin(); bc != bc_head->end(); ++bc) {
    int f = bc->first;
    bc_model[f] = FLOW_BC_FACE_PRESSURE;
    bc_values[f][0] = bc->second;
    flag_essential_bc = 1;
  }

  for (bc = bc_flux->begin(); bc != bc_flux->end(); ++bc) {
    int f = bc->first;
    bc_model[f] = FLOW_BC_FACE_FLUX;
    bc_values[f][0] = bc->second * rainfall_factor[f];
  }

  int nseepage = 0;
  double area_seepage = 0.0;
  for (bc = bc_seepage->begin(); bc != bc_seepage->end(); ++bc) {
    int f = bc->first;

    if (bc_submodel[f] & FLOW_BC_SUBMODEL_SEEPAGE_PFLOTRAN) {  // Model I.
      if (pressure_faces[0][f] < atm_pressure_) {
        bc_model[f] = FLOW_BC_FACE_FLUX;
        bc_values[f][0] = bc->second * rainfall_factor[f];
      } else {
        int c = BoundaryFaceGetCell(f);
        if (pressure_cells[0][c] < pressure_faces[0][f]) {
          bc_model[f] = FLOW_BC_FACE_FLUX;
          bc_values[f][0] = bc->second * rainfall_factor[f];
        } else {
          bc_model[f] = FLOW_BC_FACE_PRESSURE;
          bc_values[f][0] = atm_pressure_;
          flag_essential_bc = 1;
          nseepage++;
          area_seepage += mesh_->face_area(f);
        }
      }

    } else if (bc_submodel[f] & FLOW_BC_SUBMODEL_SEEPAGE_FACT) {  // Model II.
      double I = FLOW_BC_SEEPAGE_FACE_IMPEDANCE;
      double influx = bc->second * rainfall_factor[f];
      double pcreg = influx / I;
      double pcmin = 3 * pcreg / 2;
      double pcmax = pcreg / 2;

      double pc = pressure_faces[0][f] - atm_pressure_;
      if (pc < pcmin) {
        bc_model[f] = FLOW_BC_FACE_FLUX;
        bc_values[f][0] = influx;
      } else if (pc >= pcmax) {
        bc_model[f] = FLOW_BC_FACE_MIXED;
        bc_values[f][0] = I * atm_pressure_;
        bc_values[f][1] = -I;  // Impedance I should be positive.
        flag_essential_bc = 1;
        nseepage++;
        area_seepage += mesh_->face_area(f);
      } else {
        bc_model[f] = FLOW_BC_FACE_FLUX;
        double f = 2 * (pcreg - pc) / pcreg;
        double q = (7 - 2 * f - f * f) / 8;
        bc_values[f][0] = q * influx; 
      }

    } else if (bc_submodel[f] & FLOW_BC_SUBMODEL_SEEPAGE_AMANZI) {  // Model III.
      double influx = bc->second * rainfall_factor[f];
      double pcreg = -FLOW_BC_SEEPAGE_FACE_REGULARIZATION;

      double pc = pressure_faces[0][f] - atm_pressure_;
      if (pc < pcreg) {
        bc_model[f] = FLOW_BC_FACE_FLUX;
        bc_values[f][0] = influx;
      } else if (pc >= 0.0) {
        int c = BoundaryFaceGetCell(f);
        if (pressure_cells[c] < pressure_faces[f]) {
          bc_model[f] = FLOW_BC_FACE_FLUX;
          bc_values[f][0] = influx;
        } else {
          bc_model[f] = FLOW_BC_FACE_PRESSURE;
          bc_values[f][0] = atm_pressure_;
          flag_essential_bc = 1;
          nseepage++;
          area_seepage += mesh_->face_area(f);
        }
      } else {
        bc_model[f] = FLOW_BC_FACE_FLUX;
        double f = pc / pcreg;
        double q = f * f * (3 - 2 * f);
        bc_values[f][0] = q * influx;
      }
    }
  }

  // mark missing boundary conditions as zero flux conditions
  AmanziMesh::Entity_ID_List cells;
  missed_bc_faces_ = 0;
  for (int f = 0; f < nfaces_owned; f++) {
    if (bc_model[f] == FLOW_BC_FACE_NULL) {
      cells.clear();
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      int ncells = cells.size();

      if (ncells == 1) {
        bc_model[f] = FLOW_BC_FACE_FLUX;
        bc_values[f][0] = 0.0;
        missed_bc_faces_++;
      }
    }
  }

  // verify that the algebraic problem is consistent
#ifdef HAVE_MPI
  int flag = flag_essential_bc;
  mesh_->get_comm()->MaxAll(&flag, &flag_essential_bc, 1);  // find the global maximum
#endif
  if (! flag_essential_bc && vo_->getVerbLevel() >= Teuchos::VERB_LOW) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *(vo_->os()) << "WARNING: no essential boundary conditions, solver may fail" << endl;
  }

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
#ifdef HAVE_MPI
    int nseepage_tmp = nseepage;
    double area_tmp = area_seepage;
    mesh_->get_comm()->SumAll(&area_tmp, &area_seepage, 1);
    mesh_->get_comm()->SumAll(&nseepage_tmp, &nseepage, 1);
#endif
    if (MyPID == 0 && nseepage > 0 && nseepage != nseepage_prev) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *(vo_->os()) << "seepage face has changed: " << area_seepage << " [m^2]" << endl;
    }
  }
  nseepage_prev = nseepage;
}


/* ******************************************************************
*  Temporary convertion from double to tensor.                                               
****************************************************************** */
void Flow_PK::SetAbsolutePermeabilityTensor()
{
  const Epetra_MultiVector& perm = *(S_->GetFieldData("permeability")->ViewComponent("cell"));
  if (dim == 2) {
    for (int c = 0; c < K.size(); c++) {
      if (perm[0][c] == perm[1][c]) {
	K[c].init(dim, 1);
	K[c](0, 0) = perm[0][c];
      } else {
	K[c].init(dim, 2);
	K[c](0, 0) = perm[0][c];
	K[c](1, 1) = perm[1][c];
      }
    }    
  } else if (dim == 3) {
    for (int c = 0; c < K.size(); c++) {
      if (perm[0][c] == perm[1][c] && perm[0][c] == perm[2][c]) {
	K[c].init(dim, 1);
	K[c](0, 0) = perm[0][c];
      } else {
	K[c].init(dim, 2);
	K[c](0, 0) = perm[0][c];
	K[c](1, 1) = perm[1][c];
	K[c](2, 2) = perm[2][c];
      }
    }        
  }
}


/* ******************************************************************
*  Calculate inner product e^T K e in each cell.                                               
****************************************************************** */
void Flow_PK::CalculatePermeabilityFactorInWell()
{
  for (int c = 0; c < ncells_owned; c++) {
    (*Kxy)[c] = 0.0;
    int idim = std::max(1, K[c].size() - 1);
    for (int i = 0; i < idim; i++) (*Kxy)[c] += K[c](i, i);
    (*Kxy)[c] /= idim;
  }
}


/* ******************************************************************
* Add source and sink terms.                                   
****************************************************************** */
void Flow_PK::AddSourceTerms(CompositeVector& rhs)
{
  Epetra_MultiVector& rhs_cells = *rhs.ViewComponent("cell");
  Functions::FlowDomainFunction::Iterator src;

  for (src = src_sink->begin(); src != src_sink->end(); ++src) {
    int c = src->first;
    rhs_cells[0][c] += mesh_->cell_volume(c) * src->second;
  }
}


/* ******************************************************************
*                                             
****************************************************************** */
/*
void Flow_PK::AddNewtonFluxes_MFD(
    RelativePermeability& rel_perm,
    const Epetra_Vector& pressure_cells, const Epetra_Vector& flux,
    Epetra_Vector& rhs, Matrix_MFD_PLambda* matrix_operator)
{
  double rho = FS->ref_fluid_density();
  std::vector<Teuchos::SerialDenseMatrix<int, double> >& Acc_faces = matrix_operator->Acc_faces();
  Acc_faces.clear();

  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  Epetra_Vector rhs_wghost(cmap_wghost);
  for (int c = 0; c < ncells_owned; c++) rhs_wghost[c] = rhs[c];

  AmanziMesh::Entity_ID_List cells, faces;
  std::vector<int> dirs;

  const Epetra_Vector& Krel_faces = rel_perm.Krel_faces();
  const Epetra_Vector& dKdP_faces = rel_perm.dKdP_faces();

  for (int f = 0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    Teuchos::SerialDenseMatrix<int, double> Bcc(ncells, ncells);
    double factor = rho * flux[f] * dKdP_faces[f] / Krel_faces[f];

    int c1 = cells[0];
    mesh_->cell_get_faces_and_dirs(c1, &faces, &dirs);

    int k = FindPosition(f, faces);
    factor *= dirs[k];

    if (ncells == 1) {
      if (factor > 0.0) { 
        Bcc(0, 0) = factor;
        rhs_wghost[c1] += factor * pressure_cells[c1];
      }
    } else {
      int i = (factor > 0.0) ? 0 : 1;
      factor = fabs(factor);
      Bcc(i, i) = factor;
      Bcc(1 - i, i) = -factor;

      c1 = cells[i];
      int c2 = cells[1 - i];
      rhs_wghost[c1] += factor * pressure_cells[c1];
      rhs_wghost[c2] -= factor * pressure_cells[c1];
    }

    Acc_faces.push_back(Bcc);
  }

  FS->CombineGhostCell2MasterCell(rhs_wghost, Add);
  for (int c = 0; c < ncells_owned; c++) rhs[c] = rhs_wghost[c];
}
*/


/* ******************************************************************
* Updates global Darcy vector calculated by a discretization method.
* simplified implementation for a single phase flow.                                            
****************************************************************** */
void Flow_PK::AddGravityFluxes_DarcyFlux(Epetra_MultiVector& darcy_mass_flux)
{
  double rho = *(S_->GetScalarData("fluid_density"));

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  std::vector<int> flag(nfaces_wghost, 0);

  for (int c = 0; c < ncells_owned; c++) {
    AmanziGeometry::Point Kg = K[c] * gravity_;
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);

      if (f < nfaces_owned && !flag[f]) {
        darcy_mass_flux[0][f] += rho * (Kg * normal);
        flag[f] = 1;
      }
    }
  }
}


/* ******************************************************************
* Updates global Darcy vector calculated by a discretization method.                                             
****************************************************************** */
void Flow_PK::AddGravityFluxes_DarcyFlux(Epetra_MultiVector& darcy_mass_flux,
                                         RelativePermeability& rel_perm) 
{
  double rho = *(S_->GetScalarData("fluid_density"));

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  std::vector<int> flag(nfaces_wghost, 0);

  Epetra_Vector& Krel_cells = rel_perm.Krel_cells();
  Epetra_Vector& Krel_faces = rel_perm.Krel_faces();
  int method = rel_perm.method();

  for (int c = 0; c < ncells_owned; c++) {
    AmanziGeometry::Point Kg = K[c] * gravity_;
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      std::vector<double>& krel = rel_perm.Krel_amanzi()[c];

      if (f < nfaces_owned && !flag[f]) {
        if (method == FLOW_RELATIVE_PERM_NONE) {
          darcy_mass_flux[0][f] += rho * (Kg * normal);
        } else if (method == FLOW_RELATIVE_PERM_CENTERED) {
          darcy_mass_flux[0][f] += rho * (Kg * normal) * Krel_cells[c];
        } else if (method == FLOW_RELATIVE_PERM_AMANZI) {
          darcy_mass_flux[0][f] += rho * (Kg * normal) * krel[n];
        } else {
          darcy_mass_flux[0][f] += rho * (Kg * normal) * Krel_faces[f];
        }
        flag[f] = 1;
      }
    }
  }
}


/* ******************************************************************
* BDF methods need a good initial guess.
* This method gives a less smoother solution than in Flow 1.0.
* WARNING: Each owned face must have at least one owned cell. 
* Probability that this assumption is violated is close to zero. 
* Even when it happens, the code will not crash.
****************************************************************** */
void Flow_PK::DeriveFaceValuesFromCellValues(
    const Epetra_MultiVector& ucells, Epetra_MultiVector& ufaces)
{
  AmanziMesh::Entity_ID_List cells;

  for (int f = 0; f < nfaces_owned; f++) {
    cells.clear();
    mesh_->face_get_cells(f, AmanziMesh::OWNED, &cells);
    int ncells = cells.size();

    double face_value = 0.0;
    for (int n = 0; n < ncells; n++) face_value += ucells[0][cells[n]];
    ufaces[0][f] = face_value / ncells;
  }
}


/* ******************************************************************
* Calculate change of water volume per second due to boundary flux.                                          
****************************************************************** */
double Flow_PK::WaterVolumeChangePerSecond(std::vector<int>& bc_model,
                                           Epetra_MultiVector& darcy_flux)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> fdirs;

  double volume = 0.0;
  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);

    for (int i = 0; i < faces.size(); i++) {
      int f = faces[i];
      if (bc_model[f] != FLOW_BC_FACE_NULL && f < nfaces_owned) {
        if (fdirs[i] >= 0) {
          volume -= darcy_flux[0][f];
        } else {
          volume += darcy_flux[0][f];
        }
      }
    }
  }
  return volume;
}


/* ******************************************************************
* Returns the first cell attached to a boundary face.   
****************************************************************** */
int Flow_PK::BoundaryFaceGetCell(int f)
{
  AmanziMesh::Entity_ID_List cells;
  mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
  return cells[0];
}


/* ******************************************************************
* Returns position of face f in the list faces.  
****************************************************************** */
int Flow_PK::FindPosition(int f, AmanziMesh::Entity_ID_List faces)
{
  for (int i = 0; i < faces.size(); i++) {
    if (faces[i] == f) return i;
  }
  return -1;
}


/* ****************************************************************
* DEBUG: creating GMV file 
**************************************************************** */
void Flow_PK::WriteGMVfile(Teuchos::RCP<State> FS) const
{
  GMV::open_data_file(*mesh_, (std::string)"flow.gmv");
  GMV::start_data();
  GMV::write_cell_data(*(S_->GetFieldData("pressure")->ViewComponent("cell")), 0, "pressure");
  GMV::write_cell_data(*(S_->GetFieldData("water_saturation")->ViewComponent("cell")), 0, "saturation");
  GMV::close_data_file();
}

}  // namespace AmanziFlow
}  // namespace Amanzi

