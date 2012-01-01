/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include "Teuchos_ParameterList.hpp"

#include "mfd3d.hpp"

#include "Flow_PK.hpp"
#include "Flow_State.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Initiazition of fundamental flow sturctures.                                              
****************************************************************** */
void Flow_PK::Init(Teuchos::RCP<Flow_State> FS_MPC)
{
  status = FLOW_NULL;

  FS = Teuchos::rcp(new Flow_State(*FS_MPC));
  mesh_ = FS->get_mesh();
  dim = mesh_->space_dimension();

  T_internal = T_physical = dT = 0.0;
  standalone_mode = false;

  FS_nextBIG = Teuchos::rcp(new Flow_State(*FS, CopyMemory) );  
  FS_nextMPC = Teuchos::rcp(new Flow_State(*FS_nextBIG, ViewMemory));

  const Epetra_Map& cmap = mesh_->cell_map(true);
  cmin = cmap.MinLID();
  cmax = cmap.MaxLID();

  number_owned_cells = mesh_->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  cmax_owned = cmin + number_owned_cells - 1;

  const Epetra_Map& fmap = mesh_->face_map(true);
  fmin = fmap.MinLID();
  fmax = fmap.MaxLID(); 
 
  number_owned_faces = mesh_->count_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  fmax_owned = fmin + number_owned_faces - 1;
}


/* ******************************************************************
* Super-map combining cells and faces.                                                  
****************************************************************** */
Epetra_Map* Flow_PK::create_super_map()
{
  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& fmap = mesh_->face_map(false);
  
  int ndof = number_owned_cells + number_owned_faces;
  int ndof_global_cell = cmap.NumGlobalElements();
  int ndof_global = ndof_global_cell + fmap.NumGlobalElements();

  int *gids = new int[ndof];
  cmap.MyGlobalElements(&(gids[0]));
  fmap.MyGlobalElements(&(gids[number_owned_cells]));

  for (int f=0; f<number_owned_faces; f++) gids[number_owned_cells + f] += ndof_global_cell;
  Epetra_Map *map = new Epetra_Map(ndof_global, ndof, gids, 0, cmap.Comm());

  delete [] gids;
  return map;
}


/* ******************************************************************
* Add a boundary marker to used faces.                                          
****************************************************************** */
void Flow_PK::updateBoundaryConditions(BoundaryFunction* bc_pressure,
                                       BoundaryFunction* bc_head, 
                                       BoundaryFunction* bc_flux,
                                       std::vector<int>& bc_markers,
                                       std::vector<double>& bc_values)
{
  for (int n=0; n<bc_markers.size(); n++) {
    bc_markers[n] = FLOW_BC_FACE_NULL;
    bc_values[n] = 0.0;
  }

  BoundaryFunction::Iterator bc;
  for (bc=bc_pressure->begin(); bc!=bc_pressure->end(); ++bc) {
    int f = bc->first;
    bc_markers[f] = FLOW_BC_FACE_PRESSURE;
    bc_values[f] = bc->second;
  }

  for (bc=bc_head->begin(); bc!=bc_head->end(); ++bc) {
    int f = bc->first;    
    bc_markers[f] = FLOW_BC_FACE_HEAD;
    bc_values[f] = bc->second;
  }

  for (bc=bc_flux->begin(); bc!=bc_flux->end(); ++bc) {
    int f = bc->first;    
    bc_markers[f] = FLOW_BC_FACE_FLUX;
    bc_values[f] = bc->second;
  }
}


/* ******************************************************************
* Gravity fluxes are calculated w.r.t. actual normals.                                                  
****************************************************************** */
void Flow_PK::calculateGravityFluxes(
  int c, WhetStone::Tensor& K, std::vector<double>& gravity_flux) 
{
  AmanziMesh::Entity_ID_List faces;
  mesh_->cell_get_faces(c, &faces);
  int nfaces = faces.size();

  double rho = FS->ref_fluid_density();
  AmanziGeometry::Point gravity(dim); 
  for (int k=0; k<dim; k++) gravity[k] = (*(FS->get_gravity()))[k] * rho * rho;

  gravity_flux.clear();
  for (int n=0; n<nfaces; n++) {
    int f = faces[n];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    gravity_flux.push_back((K * gravity) * normal);    
  }
}


/* ******************************************************************
* Updates global Darcy vector calculated by a discretization method.                                             
****************************************************************** */
void Flow_PK::addGravityFluxes(Epetra_Vector& darcy_flux)
{
  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

  double rho = FS->ref_fluid_density();
  AmanziGeometry::Point gravity(dim); 
  for (int k=0; k<dim; k++) gravity[k] = (*(FS->get_gravity()))[k] * rho * rho;

  for (int f=0; f<nfaces; f++) {
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    darcy_flux[f] -= gravity * normal;
  }
}


/* *******************************************************************
 * Identify flux direction based on orientation of the face normal 
 * and sign of the  Darcy velocity.                               
 ****************************************************************** */
void Flow_PK::identify_upwind_cells(Epetra_Vector& upwind_cell, Epetra_Vector& downwind_cell)
{
  for (int f=fmin; f<=fmax; f++) {
    upwind_cell[f] = -1;  // negative value is indicator of a boundary
    downwind_cell[f] = -1;
  }

  AmanziMesh::Entity_ID_List faces; 
  std::vector<int> fdirs;
  Epetra_Vector& darcy_flux = FS_nextBIG->ref_darcy_flux();

  for (int c=cmin; c<=cmax; c++) {
    mesh_->cell_get_faces(c, &faces);
    mesh_->cell_get_face_dirs(c, &fdirs);

    for (int i=0; i<faces.size(); i++) {
      int f = faces[i];
      if (darcy_flux[f] * fdirs[i] >= 0) { 
        upwind_cell[f] = c; 
      } else { 
        downwind_cell[f] = c; 
      }
    }
  }
}


}  // namespace AmanziFlow
}  // namespace Amanzi

