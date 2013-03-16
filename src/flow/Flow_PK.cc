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
#include "Flow_PK.hh"
#include "Flow_State.hh"


namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Initiazition of fundamental flow sturctures.                                              
****************************************************************** */
void Flow_PK::Init(Teuchos::RCP<Flow_State> FS_MPC)
{
  flow_status_ = FLOW_STATUS_NULL;

  FS = Teuchos::rcp(new Flow_State(*FS_MPC));
  FS_aux = Teuchos::rcp(new Flow_State(*FS_MPC), AmanziFlow::FLOW_STATE_COPY);

  mesh_ = FS->mesh();
  dim = mesh_->space_dimension();
  MyPID = 0;

  T_physics = dT = 0.0;

  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);

  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  nseepage_prev = 0;
  ti_phase_counter = 0;
}


/* ******************************************************************
* Super-map combining cells and faces.                                                  
****************************************************************** */
Epetra_Map* Flow_PK::CreateSuperMap()
{
  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& fmap = mesh_->face_map(false);

  int ndof = ncells_owned + nfaces_owned;
  int ndof_global_cell = cmap.NumGlobalElements();
  int ndof_global = ndof_global_cell + fmap.NumGlobalElements();

  int* gids = new int[ndof];
  cmap.MyGlobalElements(&(gids[0]));
  fmap.MyGlobalElements(&(gids[ncells_owned]));

  for (int f = 0; f < nfaces_owned; f++) gids[ncells_owned + f] += ndof_global_cell;
  Epetra_Map* map = new Epetra_Map(ndof_global, ndof, gids, 0, cmap.Comm());

  delete [] gids;
  return map;
}


/* ******************************************************************
* Populate data needed by submodels.
****************************************************************** */
void Flow_PK::ProcessStaticBCsubmodels(const std::vector<int>& bc_submodel,
                                       std::vector<double>& rainfall_factor)
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
void Flow_PK::ProcessBoundaryConditions(
    Functions::BoundaryFunction* bc_pressure, Functions::BoundaryFunction* bc_head,
    Functions::BoundaryFunction* bc_flux, Functions::BoundaryFunction* bc_seepage,
    const Epetra_Vector& pressure_cells, 
    const Epetra_Vector& pressure_faces, const double atm_pressure,
    const std::vector<double>& rainfall_factor,
    const std::vector<int>& bc_submodel,
    std::vector<int>& bc_model, std::vector<bc_tuple>& bc_values)
{
  int flag_essential_bc = 0;
  bc_tuple zero = {0.0, 0.0};
  for (int n = 0; n < bc_model.size(); n++) {
    bc_model[n] = FLOW_BC_FACE_NULL;
    bc_values[n] = zero;
  }

  Functions::BoundaryFunction::Iterator bc;
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
      if (pressure_faces[f] < atm_pressure) {
        bc_model[f] = FLOW_BC_FACE_FLUX;
        bc_values[f][0] = bc->second * rainfall_factor[f];
      } else {
        int c = BoundaryFaceGetCell(f);
        if (pressure_cells[c] < pressure_faces[f]) {
          bc_model[f] = FLOW_BC_FACE_FLUX;
          bc_values[f][0] = bc->second * rainfall_factor[f];
        } else {
          bc_model[f] = FLOW_BC_FACE_PRESSURE;
          bc_values[f][0] = atm_pressure;
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

      double pc = pressure_faces[f] - atm_pressure;
      if (pc < pcmin) {
        bc_model[f] = FLOW_BC_FACE_FLUX;
        bc_values[f][0] = influx;
      } else if (pc >= pcmax) {
        bc_model[f] = FLOW_BC_FACE_MIXED;
        bc_values[f][0] = I * atm_pressure;
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

      double pc = pressure_faces[f] - atm_pressure;
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
          bc_values[f][0] = atm_pressure;
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
  if (! flag_essential_bc && MyPID == 0 && verbosity >= FLOW_VERBOSITY_LOW) {
    std::printf("Flow PK: WARNING: no essential boundary conditions, solver may fail\n");
  }

  // verbose output
  if (verbosity >= FLOW_VERBOSITY_HIGH) {
#ifdef HAVE_MPI
    int nseepage_tmp = nseepage;
    double area_tmp = area_seepage;
    mesh_->get_comm()->SumAll(&area_tmp, &area_seepage, 1);
    mesh_->get_comm()->SumAll(&nseepage_tmp, &nseepage, 1);
#endif
    if (MyPID == 0 && nseepage > 0 && nseepage != nseepage_prev) {
      std::printf("Flow PK: seepage face has changed: %9.4e [m^2]\n", area_seepage);
    }
  }
  nseepage_prev = nseepage;
}


/* ******************************************************************
*  Calculate inner product e^T K e in each cell.                                               
****************************************************************** */
void Flow_PK::CalculatePermeabilityFactorInWell(const std::vector<WhetStone::Tensor>& K, Epetra_Vector& Kxy)
{
  for (int c = 0; c < ncells_owned; c++) {
    Kxy[c] = 0.0;
    int idim = std::max(1, K[c].size() - 1);
    for (int i = 0; i < idim; i++) Kxy[c] += K[c](i, i);
    Kxy[c] /= idim;
  }
}


/* ******************************************************************
* Add source and sink terms. We use a simplified algorithms than for
* boundary conditions.                                          
****************************************************************** */
void Flow_PK::AddSourceTerms(Functions::UniqueMeshFunction* src_sink, Epetra_Vector& rhs)
{
  // Functions::UniqueMeshFunction::Iterator src;
  // for (src = src_sink->begin(); src != src_sink->end(); ++src) {
  //   int c = src->first;
  //   rhs[c] += mesh_->cell_volume(c) * src->second;
  // }
}


/* ******************************************************************
* Routine updates elemental discretization matrices and must be 
* called before applying boundary conditions and global assembling.                                             
****************************************************************** */
void Flow_PK::AddGravityFluxes_MFD(
    std::vector<WhetStone::Tensor>& K,
    const Epetra_Vector& Krel_cells, const Epetra_Vector& Krel_faces, int method,
    Matrix_MFD* matrix_operator)
{
  double rho = FS->ref_fluid_density();
  AmanziGeometry::Point gravity(dim);
  for (int k = 0; k < dim; k++) gravity[k] = (*(FS->gravity()))[k] * rho;

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    Epetra_SerialDenseVector& Ff = matrix_operator->Ff_cells()[c];
    double& Fc = matrix_operator->Fc_cells()[c];

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);

      double outward_flux = ((K[c] * gravity) * normal) * dirs[n]; 
      if (method == FLOW_RELATIVE_PERM_CENTERED) {
        outward_flux *= Krel_cells[c];
      } else {
        outward_flux *= Krel_faces[f];
      }
      Ff[n] += outward_flux;
      Fc -= outward_flux;  // Nonzero-sum contribution when flag_upwind = false.
    }
  }
}


/* ******************************************************************
*                                             
****************************************************************** */
void Flow_PK::AddNewtonFluxes_MFD(
    const Epetra_Vector& dKdP_faces, const Epetra_Vector& Krel_faces,
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


/* ******************************************************************
* Updates global Darcy vector calculated by a discretization method.                                             
****************************************************************** */
void Flow_PK::AddGravityFluxes_DarcyFlux(
    std::vector<WhetStone::Tensor>& K,
    const Epetra_Vector& Krel_cells, const Epetra_Vector& Krel_faces, int method,
    Epetra_Vector& darcy_mass_flux)
{
  double rho = FS->ref_fluid_density();
  AmanziGeometry::Point gravity(dim);
  for (int k = 0; k < dim; k++) gravity[k] = (*(FS->gravity()))[k] * rho;

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  std::vector<int> flag(nfaces_wghost, 0);

  for (int c = 0; c < ncells_owned; c++) {
    // AmanziGeometry::Point Kg = K[c] * gravity;
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);

      if (f < nfaces_owned && !flag[f]) {
        if (method == FLOW_RELATIVE_PERM_NONE) {
           darcy_mass_flux[f] += ((K[c] * gravity) * normal);
        } else if (method == FLOW_RELATIVE_PERM_CENTERED) {
          darcy_mass_flux[f] += ((K[c] * gravity) * normal) * Krel_cells[c];
        } else {
          darcy_mass_flux[f] += ((K[c] * gravity) * normal) * Krel_faces[f];
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
void Flow_PK::DeriveFaceValuesFromCellValues(const Epetra_Vector& ucells, Epetra_Vector& ufaces)
{
  AmanziMesh::Entity_ID_List cells;

  for (int f = 0; f < nfaces_owned; f++) {
    cells.clear();
    mesh_->face_get_cells(f, AmanziMesh::OWNED, &cells);
    int ncells = cells.size();

    double face_value = 0.0;
    for (int n = 0; n < ncells; n++) face_value += ucells[cells[n]];
    ufaces[f] = face_value / ncells;
  }
}


/* *******************************************************************
* Identify flux direction based on orientation of the face normal 
* and sign of the Darcy velocity. 
* WARNING: It is *not* used now.                              
******************************************************************* */
void Flow_PK::IdentifyUpwindCells(Epetra_IntVector& upwind_cell, Epetra_IntVector& downwind_cell)
{
  for (int f = 0; f < nfaces_owned; f++) {
    upwind_cell[f] = -1;  // negative value is indicator of a boundary
    downwind_cell[f] = -1;
  }

  Epetra_Vector& darcy_flux = FS->ref_darcy_flux();
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> fdirs;

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);

    for (int i = 0; i < faces.size(); i++) {
      int f = faces[i];
      if (darcy_flux[f] * fdirs[i] >= 0)
        upwind_cell[f] = c;
      else
        downwind_cell[f] = c;
    }
  }
}


/* ******************************************************************
* Calculate change of water volume per second due to boundary flux.                                          
****************************************************************** */
double Flow_PK::WaterVolumeChangePerSecond(std::vector<int>& bc_model, Epetra_Vector& darcy_flux)
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
          volume -= darcy_flux[f];
        } else {
          volume += darcy_flux[f];
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
void Flow_PK::WriteGMVfile(Teuchos::RCP<Flow_State> FS) const
{
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = FS->mesh();

  GMV::open_data_file(*mesh, (std::string)"flow.gmv");
  GMV::start_data();
  GMV::write_cell_data(FS->ref_pressure(), "pressure");
  GMV::write_cell_data(FS->ref_water_saturation(), "saturation");
  GMV::write_cell_data(FS->ref_darcy_velocity(), 0, "velocity_h");
  GMV::write_cell_data(FS->ref_darcy_velocity(), dim-1, "velocity_v");
  GMV::close_data_file();
}

}  // namespace AmanziFlow
}  // namespace Amanzi

