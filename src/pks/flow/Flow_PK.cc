/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
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
#include "Mesh.hh"
#include "mfd3d.hh"
#include "OperatorDefs.hh"
#include "State.hh"

#include "Flow_PK.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* default constructor that initializes all pointers to NULL
****************************************************************** */
Flow_PK::Flow_PK() :
    bc_pressure(NULL),
    bc_flux(NULL),
    bc_head(NULL),
    bc_seepage(NULL),
    src_sink(NULL),
    ti_specs(NULL),
    vo_(NULL),
    passwd_("state")
{
}

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
  gravity_.set(dim,&(gravity_data[0])); // do it in complicated way because we
                                        // are not sure if gravity_data is an
                                        // array or vector
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
void Flow_PK::ComputeBCs(const CompositeVector& u)
{
  const Epetra_MultiVector& u_cell = *u.ViewComponent("cell");
  
  int flag_essential_bc = 0;
  dirichlet_bc_faces_ = 0;
  for (int n = 0; n < bc_model.size(); n++) {
    bc_model[n] = Operators::OPERATOR_BC_NONE;
    bc_value[n] = 0.0;
    bc_coef[n] = 0.0;
  }

  Functions::FlowBoundaryFunction::Iterator bc;
  for (bc = bc_pressure->begin(); bc != bc_pressure->end(); ++bc) {
    int f = bc->first;
    bc_model[f] = Operators::OPERATOR_BC_FACE_DIRICHLET;
    bc_value[f] = bc->second;
    flag_essential_bc = 1;
    dirichlet_bc_faces_++;
  }

  for (bc = bc_head->begin(); bc != bc_head->end(); ++bc) {
    int f = bc->first;
    if (bc_submodel[f] & FLOW_BC_SUBMODEL_NOFLOW_ABOVE_WATER_TABLE) {
      if (bc->second < FLOW_PRESSURE_ATMOSPHERIC) {
        bc_model[f] = Operators::OPERATOR_BC_FACE_NEUMANN;
        bc_value[f] = 0.0;
        continue;
      }
    }
    bc_model[f] = Operators::OPERATOR_BC_FACE_DIRICHLET;
    bc_value[f] = bc->second;
    flag_essential_bc = 1;
    dirichlet_bc_faces_++;
  }

  for (bc = bc_flux->begin(); bc != bc_flux->end(); ++bc) {
    int f = bc->first;
    bc_model[f] = Operators::OPERATOR_BC_FACE_NEUMANN;
    bc_value[f] = bc->second * rainfall_factor[f];
  }

  // Seepage face BC is implemented for p-lambda discretization only.
  int nseepage = 0;
  double area_seepage = 0.0;

  if (u.HasComponent("face")) {
    const Epetra_MultiVector& u_face = *u.ViewComponent("face");
    double ref_pressure = bc_seepage->reference_pressure();
    double tol = ref_pressure * 1e-14;

    for (bc = bc_seepage->begin(); bc != bc_seepage->end(); ++bc) {
      int f = bc->first;

      // Model I. Due to round-off errors, we have to use tolerances.
      if (bc_submodel[f] & FLOW_BC_SUBMODEL_SEEPAGE_PFLOTRAN) {
        if (u_face[0][f] < ref_pressure - tol) {
          bc_model[f] = Operators::OPERATOR_BC_FACE_NEUMANN;
          bc_value[f] = bc->second * rainfall_factor[f];
        } else {
          int c = BoundaryFaceGetCell(f);
          if (u_cell[0][c] < u_face[0][f]) {
            bc_model[f] = Operators::OPERATOR_BC_FACE_NEUMANN;
            bc_value[f] = bc->second * rainfall_factor[f];
          } else {
            bc_model[f] = Operators::OPERATOR_BC_FACE_DIRICHLET;
            bc_value[f] = ref_pressure;
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

        double pc = u_face[0][f] - ref_pressure;
        if (pc < pcmin) {
          bc_model[f] = Operators::OPERATOR_BC_FACE_NEUMANN;
          bc_value[f] = influx;
        } else if (pc >= pcmax) {
          bc_model[f] = Operators::OPERATOR_BC_FACE_MIXED;
          bc_value[f] = I * ref_pressure;
          bc_coef[f] = -I;  // Impedance I should be positive.
          flag_essential_bc = 1;
          nseepage++;
          area_seepage += mesh_->face_area(f);
        } else {
          bc_model[f] = Operators::OPERATOR_BC_FACE_NEUMANN;
          double f = 2 * (pcreg - pc) / pcreg;
          double q = (7 - 2 * f - f * f) / 8;
          bc_value[f] = q * influx; 
        }

      } else if (bc_submodel[f] & FLOW_BC_SUBMODEL_SEEPAGE_AMANZI) {  // Model III.
        double influx = bc->second * rainfall_factor[f];
        double pcreg = -FLOW_BC_SEEPAGE_FACE_REGULARIZATION;

        double pc = u_face[0][f] - ref_pressure;
        if (pc < pcreg) {
          bc_model[f] = Operators::OPERATOR_BC_FACE_NEUMANN;
          bc_value[f] = influx;
        } else if (pc >= 0.0) {
          int c = BoundaryFaceGetCell(f);
          if (u_cell[0][c] < u_face[0][f]) {
            bc_model[f] = Operators::OPERATOR_BC_FACE_NEUMANN;
            bc_value[f] = influx;
          } else {
            bc_model[f] = Operators::OPERATOR_BC_FACE_DIRICHLET;
            bc_value[f] = ref_pressure;
            flag_essential_bc = 1;
            nseepage++;
            area_seepage += mesh_->face_area(f);
          }
        } else {
          bc_model[f] = Operators::OPERATOR_BC_FACE_NEUMANN;
          double f = pc / pcreg;
          double q = f * f * (3 - 2 * f);
          bc_value[f] = q * influx;
        }
      }
    }
  }
  else {
    double face_value;
    double ref_pressure = bc_seepage->reference_pressure();
    double tol = ref_pressure * 1e-14;

    for (bc = bc_seepage->begin(); bc != bc_seepage->end(); ++bc) {
      int f = bc->first;
      int c = BoundaryFaceGetCell(f);
      face_value = BoundaryFaceValue(f, u);

      if (bc_submodel[f] & FLOW_BC_SUBMODEL_SEEPAGE_PFLOTRAN) {   // Model I
    
	if (face_value < ref_pressure - tol) {
	  bc_model[f] = Operators::OPERATOR_BC_FACE_NEUMANN;
	  bc_value[f] = bc->second * rainfall_factor[f];
	} else {
	  int c = BoundaryFaceGetCell(f);
	  if (u_cell[0][c] < face_value) {
	    bc_model[f] = Operators::OPERATOR_BC_FACE_NEUMANN;
	    bc_value[f] = bc->second * rainfall_factor[f];
	  } else {
	    bc_model[f] = Operators::OPERATOR_BC_FACE_DIRICHLET;
	    bc_value[f] = ref_pressure;
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

        double pc = face_value - ref_pressure;
        if (pc < pcmin) {
          bc_model[f] = Operators::OPERATOR_BC_FACE_NEUMANN;
          bc_value[f] = influx;
        } else if (pc >= pcmax) {
          bc_model[f] = Operators::OPERATOR_BC_FACE_MIXED;
          bc_value[f] = I * ref_pressure;
          bc_coef[f] = -I;  // Impedance I should be positive.
          flag_essential_bc = 1;
          nseepage++;
          area_seepage += mesh_->face_area(f);
        } else {
          bc_model[f] = Operators::OPERATOR_BC_FACE_NEUMANN;
          double f = 2 * (pcreg - pc) / pcreg;
          double q = (7 - 2 * f - f * f) / 8;
          bc_value[f] = q * influx; 
        }
      } else if (bc_submodel[f] & FLOW_BC_SUBMODEL_SEEPAGE_AMANZI) {  // Model III.
        double influx = bc->second * rainfall_factor[f];
        double pcreg = -FLOW_BC_SEEPAGE_FACE_REGULARIZATION;

        double pc = face_value - ref_pressure;
        if (pc < pcreg) {
          bc_model[f] = Operators::OPERATOR_BC_FACE_NEUMANN;
          bc_value[f] = influx;
        } else if (pc >= 0.0) {
          int c = BoundaryFaceGetCell(f);
          if (u_cell[0][c] < face_value) {
            bc_model[f] = Operators::OPERATOR_BC_FACE_NEUMANN;
            bc_value[f] = influx;
          } else {
            bc_model[f] = Operators::OPERATOR_BC_FACE_DIRICHLET;
            bc_value[f] = ref_pressure;
            flag_essential_bc = 1;
            nseepage++;
            area_seepage += mesh_->face_area(f);
          }
        } else {
          bc_model[f] = Operators::OPERATOR_BC_FACE_NEUMANN;
          double f = pc / pcreg;
          double q = f * f * (3 - 2 * f);
          bc_value[f] = q * influx;
        }
      }

    }
  }

  // mark missing boundary conditions as zero flux conditions
  AmanziMesh::Entity_ID_List cells;
  missed_bc_faces_ = 0;
  for (int f = 0; f < nfaces_owned; f++) {
    if (bc_model[f] == Operators::OPERATOR_BC_NONE) {
      cells.clear();
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      int ncells = cells.size();

      if (ncells == 1) {
        bc_model[f] = Operators::OPERATOR_BC_FACE_NEUMANN;
        bc_value[f] = 0.0;
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
    *vo_->os() << "WARNING: no essential boundary conditions, solver may fail" << std::endl;
  }

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
#ifdef HAVE_MPI
    int nseepage_tmp = nseepage;
    double area_tmp = area_seepage;
    mesh_->get_comm()->SumAll(&area_tmp, &area_seepage, 1);
    mesh_->get_comm()->SumAll(&nseepage_tmp, &nseepage, 1);
#endif
    if (MyPID == 0 && nseepage > 0 && nseepage != nseepage_prev) {
    //if (MyPID == 0 && nseepage > 0) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "seepage face has changed: " << area_seepage << " [m^2] "<<nseepage<<" "<<nseepage_prev<< std::endl;
    }
  }
  nseepage_prev = nseepage;
}


/* ******************************************************************
*  Temporary convertion from double to tensor.                                               
****************************************************************** */
void Flow_PK::SetAbsolutePermeabilityTensor()
{
  const CompositeVector& cv = *S_->GetFieldData("permeability");
  cv.ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& perm = *cv.ViewComponent("cell", true);

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
  Epetra_MultiVector& rhs_cell = *rhs.ViewComponent("cell");
  Functions::FlowDomainFunction::Iterator src;

  for (src = src_sink->begin(); src != src_sink->end(); ++src) {
    int c = src->first;
    rhs_cell[0][c] += mesh_->cell_volume(c) * src->second;
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
  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

  for (int f = 0; f < nfaces; f++) {
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
      if (bc_model[f] != Operators::OPERATOR_BC_NONE && f < nfaces_owned) {
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
* Returns approximation of a solution on a boundary face   
****************************************************************** */

double Flow_PK::BoundaryFaceValue(int f, const CompositeVector& pressure){

  const Epetra_MultiVector& u_cell = *pressure.ViewComponent("cell");
  int c = BoundaryFaceGetCell(f);
  return u_cell[0][c];

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

}  // namespace Flow
}  // namespace Amanzi

