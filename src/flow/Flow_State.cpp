/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"

#include "Point.hh"
#include "errors.hh"
#include "Mesh.hh"

#include "State.hpp"
#include "Flow_State.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* *******************************************************************
* Flow state is build from scratch and filled with zeros.         
******************************************************************* */
Flow_State::Flow_State(Teuchos::RCP<AmanziMesh::Mesh> mesh)
{
  const Epetra_BlockMap& cmap = mesh->cell_map(false);
  const Epetra_BlockMap& fmap = mesh->face_map(false);

  porosity = Teuchos::rcp(new Epetra_Vector(cmap));
  fluid_density = Teuchos::rcp(new double);
  fluid_viscosity = Teuchos::rcp(new double);
  gravity = Teuchos::rcp(new double*);
  *gravity = new double[3];
  vertical_permeability = Teuchos::rcp(new Epetra_Vector(cmap));
  horizontal_permeability = Teuchos::rcp(new Epetra_Vector(cmap));
  pressure = Teuchos::rcp(new Epetra_Vector(cmap));
  darcy_mass_flux = Teuchos::rcp(new Epetra_Vector(fmap));
  mesh_ = mesh;

  S_ = NULL;  
}


/* *******************************************************************
* Flow state is build from state S.        
******************************************************************* */
Flow_State::Flow_State(Teuchos::RCP<State> S)
{
  porosity = S->get_porosity();
  fluid_density = S->get_density();
  fluid_viscosity = S->get_viscosity();
  gravity = S->get_gravity();
  vertical_permeability = S->get_vertical_permeability();
  horizontal_permeability = S->get_horizontal_permeability();
  pressure = S->get_pressure();
  darcy_mass_flux = S->get_darcy_flux();
  mesh_ = S->get_mesh_maps();

  S_ = &*S;
}


Flow_State::Flow_State(State& S)
{
  porosity = S.get_porosity();
  fluid_density = S.get_density();
  fluid_viscosity = S.get_viscosity();
  gravity = S.get_gravity();
  vertical_permeability = S.get_vertical_permeability();
  horizontal_permeability = S.get_horizontal_permeability();
  pressure = S.get_pressure();
  darcy_mass_flux = S.get_darcy_flux();
  mesh_ = S.get_mesh_maps();

  S_ = &S;
}


/* *******************************************************************
* mode = CopyPointers (default) a trivial copy of the given state           
* mode = ViewMemory   creates the flow state from internal one   
*                     as the MPC expected                     
******************************************************************* */
Flow_State::Flow_State(Flow_State& S, FlowCreateMode mode)
{
  if (mode == CopyPointers) {
    porosity = S.get_porosity();
    fluid_density = S.get_fluid_density();
    fluid_viscosity = S.get_fluid_viscosity();
    gravity = S.get_gravity();
    vertical_permeability = S.get_vertical_permeability();
    horizontal_permeability = S.get_horizontal_permeability();
    pressure = S.get_pressure();
    darcy_mass_flux = S.get_darcy_mass_flux();
    mesh_ = S.get_mesh();
  } 
  else if (mode == CopyMemory ) { 
    porosity = S.get_porosity();
    fluid_density = S.get_fluid_density();
    fluid_viscosity = S.get_fluid_viscosity();
    gravity = S.get_gravity();
    vertical_permeability = S.get_vertical_permeability();
    horizontal_permeability = S.get_horizontal_permeability();
    mesh_ = S.get_mesh();

    // allocate memory for the next state
    pressure = Teuchos::rcp(new Epetra_Vector(S.ref_pressure()));
    darcy_mass_flux = Teuchos::rcp(new Epetra_Vector(S.ref_darcy_mass_flux()));
  }

  S_ = S.S_;
}


/* *******************************************************************
* Copy cell-based data from master to ghost positions.              
* WARNING: vector v must contain ghost cells.              
******************************************************************* */
void Flow_State::copyMasterCell2GhostCell(Epetra_Vector& v)
{
#ifdef HAVE_MPI
  const Epetra_BlockMap& source_cmap = mesh_->cell_map(false);
  const Epetra_BlockMap& target_cmap = mesh_->cell_map(true);
  Epetra_Import importer(target_cmap, source_cmap);

  double* vdata;
  v.ExtractView(&vdata);
  Epetra_Vector vv(View, source_cmap, vdata);

  v.Import(vv, importer, Insert);
#endif
}


/* *******************************************************************
* Transfer face-based data from master to ghost positions and perform
* operation mode there. 
* WARNING: Vector v must contain ghost faces.              
******************************************************************* */
void Flow_State::combineGhostFace2MasterFace(Epetra_Vector& v, Epetra_CombineMode mode)
{
#ifdef HAVE_MPI
  const Epetra_BlockMap& target_fmap = mesh_->face_map(true);
  const Epetra_BlockMap& source_fmap = mesh_->face_map(false);
  Epetra_Import importer(target_fmap, source_fmap);

  double* vdata;
  v.ExtractView(&vdata);
  Epetra_Vector vv(View, source_fmap, vdata);

  vv.Export(v, importer, mode);
#endif
}


/* *******************************************************************
* Copy cell-based data from master to ghost positions.              
* WARNING: MultiVector v must contain ghost cells.              
******************************************************************* */
void Flow_State::copyMasterMultiCell2GhostMultiCell(Epetra_MultiVector& v)
{
#ifdef HAVE_MPI
  const Epetra_BlockMap& source_cmap = mesh_->cell_map(false);
  const Epetra_BlockMap& target_cmap = mesh_->cell_map(true);
  Epetra_Import importer(target_cmap, source_cmap);

  double** vdata;
  v.ExtractView(&vdata);
  Epetra_MultiVector vv(View, source_cmap, vdata, v.NumVectors());

  v.Import(vv, importer, Insert);
#endif
}


/* *******************************************************************
* L2 norms (USED only ONCE).              
******************************************************************* */
double Flow_State::normL2cell(Epetra_Vector& v1, Epetra_Vector& v2)
{
  int ncells = (mesh_->cell_map(false)).NumMyElements();

  double L2error = 0.0;
  for (int c=0; c<ncells; c++) {
    double volume = mesh_->cell_volume(c);
    L2error += volume * pow(v1[c] - v2[c], 2.0);
  }
  return sqrt(L2error);
}


double Flow_State::normL2cell(Epetra_Vector& v1)
{
  int ncells = (mesh_->cell_map(false)).NumMyElements();

  double L2norm = 0.0;
  for (int c=0; c<ncells; c++) {
    double volume = mesh_->cell_volume(c);
    L2norm += volume * pow(v1[c], 2.0);
  }
  return sqrt(L2norm);
}


/* *******************************************************************
* Extract cells from a supervector             
******************************************************************* */
Epetra_Vector* Flow_State::createCellView(const Epetra_Vector& u) const
{
  double* data;
  u.ExtractView(&data);
  return new Epetra_Vector(View, mesh_->cell_map(false), data);
}


/* *******************************************************************
* Extract faces from a supervector             
******************************************************************* */
Epetra_Vector* Flow_State::createFaceView(const Epetra_Vector& u) const
{
  double* data;
  u.ExtractView(&data);
  int ncells = (mesh_->cell_map(false)).NumMyElements();
  return new Epetra_Vector(View, mesh_->face_map(false), data+ncells);
}


/* *******************************************************************
* DEBUG: create constant fluid density    
******************************************************************* */
void Flow_State::set_fluid_density(double rho)
{
  *fluid_density = rho;  // verify that it is positive (lipnikov@lanl.gov)
}


/* *******************************************************************
 * DEBUG: create constant fluid viscosity
 ****************************************************************** */
void Flow_State::set_fluid_viscosity(double mu)
{
  *fluid_viscosity = mu;  // verify that it is positive (lipnikov@lanl.gov)
}


/* *******************************************************************
 * DEBUG: create constant porosity
 ****************************************************************** */
void Flow_State::set_porosity(double phi)
{
  porosity->PutScalar(phi);
}


/* *******************************************************************
 * DEBUG: create hydrostatic pressure with p0 at height z0.
 ****************************************************************** */
void Flow_State::set_pressure_hydrostatic(double z0, double p0)
{
  int dim = mesh_->space_dimension();
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);

  double rho = *fluid_density;
  double g = (*gravity)[dim - 1];

  for (int c=0; c<ncells; c++) {
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
    (*pressure)[c] = p0 + rho * g * (xc[dim - 1] - z0);
  }
}


/* *******************************************************************
 * DEBUG: create diagonal permeability
 ****************************************************************** */
void Flow_State::set_permeability(double Kh, double Kv)
{
  horizontal_permeability->PutScalar(Kh);
  vertical_permeability->PutScalar(Kv);
}

void Flow_State::set_permeability(double Kh, double Kv, const string region)
{
  std::vector<unsigned int> block;
  mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);
  int ncells = block.size();
      
  for (int i=0; i<ncells; i++) {
    int c = block[i];
    (*horizontal_permeability)[c] = Kh;
    (*vertical_permeability)[c] = Kv;
  }
}


/* *******************************************************************
 * DEBUG: create constant gravity
 ****************************************************************** */
void Flow_State::set_gravity(double g)
{
  int dim = mesh_->space_dimension();
  (*gravity)[dim-1] = g;
}



}  // namespace AmanziTransport
}  // namespace Amanzi

