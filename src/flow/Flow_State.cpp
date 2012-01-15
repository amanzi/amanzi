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

#include "State.hpp"
#include "Flow_State.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* **************************************************************** */
Flow_State::Flow_State(Teuchos::RCP<State> S)
{
  porosity = S->get_porosity();
  fluid_density = S->get_density();
  fluid_viscosity = S->get_viscosity();
  gravity = S->get_gravity();
  absolute_permeability = S->get_permeability();
  pressure = S->get_pressure();
  darcy_flux = S->get_darcy_flux();
  mesh_ = S->get_mesh_maps();

  S_ = &*S;
}


Flow_State::Flow_State(State& S)
{
  porosity = S.get_porosity();
  fluid_density = S.get_density();
  fluid_viscosity = S.get_viscosity();
  gravity = S.get_gravity();
  absolute_permeability = S.get_permeability();
  pressure = S.get_pressure();
  darcy_flux = S.get_darcy_flux();
  mesh_ = S.get_mesh_maps();

  S_ = &S;
}


/* *******************************************************************
 * mode = CopyPointers (default) a trivial copy of the given state           
 * mode = ViewMemory   creates the flow state from internal one   
 *                     as the MPC expected                     
 * mode = CopyMemory   creates internal flow state based on the
 *                     ovelapped mesh maps                       
 * **************************************************************** */
Flow_State::Flow_State(Flow_State& S, FlowCreateMode mode)
{
  if (mode == CopyPointers) {
    porosity = S.get_porosity();
    fluid_density = S.get_fluid_density();
    fluid_viscosity = S.get_fluid_viscosity();
    gravity = S.get_gravity();
    absolute_permeability = S.get_absolute_permeability();
    pressure = S.get_pressure();
    darcy_flux = S.get_darcy_flux();
    mesh_ = S.get_mesh();
  }
  else if (mode == CopyMemory ) { 
    porosity = S.get_porosity();
    fluid_density = S.get_fluid_density();
    fluid_viscosity = S.get_fluid_viscosity();
    gravity = S.get_gravity();
    absolute_permeability = S.get_absolute_permeability();
    mesh_ = S.get_mesh();

    // allocate memory for internal state
    const Epetra_Map& cmap = mesh_->cell_map(true);
    const Epetra_Map& fmap = mesh_->face_map(true);

    pressure = Teuchos::rcp(new Epetra_Vector(cmap));
    copyMemoryVector(S.ref_pressure(), *pressure);

    darcy_flux = Teuchos::rcp(new Epetra_Vector(fmap));
    copyMemoryVector(S.ref_darcy_flux(), *darcy_flux);
  }

  else if (mode == ViewMemory) {
    porosity = S.get_porosity(); 
    fluid_density = S.get_fluid_density();
    fluid_viscosity = S.get_fluid_viscosity();
    gravity = S.get_gravity();
    absolute_permeability = S.get_absolute_permeability();
    mesh_ = S.get_mesh();

    double *data_dp, *data_ddf;
    const Epetra_Map& cmap = mesh_->cell_map(false);
    const Epetra_Map& fmap = mesh_->face_map(false);

    Epetra_Vector& dp = S.ref_pressure();
    dp.ExtractView(&data_dp);     
    pressure = Teuchos::rcp(new Epetra_Vector(View, cmap, data_dp));

    Epetra_Vector& ddf = S.ref_darcy_flux();
    ddf.ExtractView(&data_ddf);     
    darcy_flux = Teuchos::rcp(new Epetra_Vector(View, fmap, data_ddf));
  }

  S_ = S.S_;
}


/* *******************************************************************
 * Routine imports a short multivector to a parallel overlaping vector.
 ****************************************************************** */
void Flow_State::copyMemoryMultiVector(Epetra_MultiVector& source, 
                                       Epetra_MultiVector& target)
{
  const Epetra_BlockMap& source_cmap = source.Map();
  const Epetra_BlockMap& target_cmap = target.Map();

  int cmin, cmax, cmax_s, cmax_t;
  cmin = source_cmap.MinLID();
  cmax_s = source_cmap.MaxLID();
  cmax_t = target_cmap.MaxLID();
  cmax = std::min(cmax_s, cmax_t);

  int number_vectors = source.NumVectors();
  for (int c=cmin; c<=cmax; c++) {
    for (int i=0; i<number_vectors; i++) target[i][c] = source[i][c];
  }

#ifdef HAVE_MPI
  if (cmax_s > cmax_t) {
    Errors::Message msg;
    msg << "Source map (in copy_multivector) is larger than target map.\n";
    Exceptions::amanzi_throw(msg);
  }

  Epetra_Import importer(target_cmap, source_cmap);
  target.Import(source, importer, Insert);
#endif
}


/* *******************************************************************
 * Routine imports a short vector to a parallel overlaping vector.                
 ****************************************************************** */
void Flow_State::copyMemoryVector(Epetra_Vector& source, Epetra_Vector& target)
{
  const Epetra_BlockMap& source_fmap = source.Map();
  const Epetra_BlockMap& target_fmap = target.Map();

  int fmin, fmax, fmax_s, fmax_t;
  fmin   = source_fmap.MinLID();
  fmax_s = source_fmap.MaxLID();
  fmax_t = target_fmap.MaxLID();
  fmax   = std::min(fmax_s, fmax_t);

  for (int f=fmin; f<=fmax; f++) target[f] = source[f];

#ifdef HAVE_MPI
  if (fmax_s > fmax_t)  {
    Errors::Message msg;
    msg << "Source map (in copy_vector) is larger than target map.\n";
    Exceptions::amanzi_throw(msg);
  }

  Epetra_Import importer(target_fmap, source_fmap);
  target.Import(source, importer, Insert);
#endif
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
 * DEBUG: create hydrostatic pressure with p0 at height z0.
 ****************************************************************** */
void Flow_State::set_pressure_head(double z0, double p0, Epetra_Vector& pressure)
{
  int dim = mesh_->space_dimension();
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);

  double rho = *fluid_density;
  double g = (*gravity)[dim - 1];

  for (int c=0; c<ncells; c++) {
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
    pressure[c] = p0 + rho * g * (xc[dim - 1] - z0);
  }

  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  if (pressure.MyLength() == ncells + nfaces) {
    for (int f=0; f<nfaces; f++) {
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
      pressure[ncells + f] = p0 + rho * g * (xf[dim - 1] - z0);
    }
  }
}



}  // namespace AmanziTransport
}  // namespace Amanzi

