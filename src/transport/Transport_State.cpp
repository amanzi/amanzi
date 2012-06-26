/*
This is the transport component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <algorithm>

#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"

#include "Point.hh"
#include "errors.hh"

#include "State.hpp"
#include "Transport_State.hpp"

namespace Amanzi {
namespace AmanziTransport {

/* *******************************************************************
* Create Flow state from a state.                     
******************************************************************* */
Transport_State::Transport_State(State& S)
{
  total_component_concentration_ = S.get_total_component_concentration();
  porosity_ = S.get_porosity();
  darcy_flux_ = S.get_darcy_flux();
  water_saturation_ = S.get_water_saturation();
  prev_water_saturation_ = S.get_prev_water_saturation();
  water_density_ = S.get_water_density();
  mesh_ = S.get_mesh_maps();

  S_ = &S;
}


/* *******************************************************************
* mode = CopyPointers (default) a trivial copy of the given state           
* mode = ViewMemory   creates the transport from internal one   
*                     as the MPC expected                     
* mode = CopyMemory   creates internal transport state based on 
*                     ovelapped mesh maps                       
******************************************************************* */
Transport_State::Transport_State(Transport_State& S, TransportCreateMode mode)
{
  if (mode == CopyPointers) {
    total_component_concentration_ = S.total_component_concentration();
    porosity_ = S.porosity();
    darcy_flux_ = S.darcy_flux();
    water_saturation_ = S.water_saturation();
    prev_water_saturation_ = S.prev_water_saturation();
    water_density_ = S.water_density();
    mesh_ = S.mesh();

  } else if (mode == CopyMemory) {
    porosity_ = S.porosity();
    water_saturation_ = S.water_saturation();
    prev_water_saturation_ = S.prev_water_saturation();
    water_density_ = S.water_density();
    mesh_ = S.mesh();

    // allocate memory for internal state
    const Epetra_Map& cmap = mesh_->cell_map(true);
    const Epetra_Map& fmap = mesh_->face_map(true);

    int number_vectors = S.total_component_concentration()->NumVectors();

    total_component_concentration_ = Teuchos::rcp(new Epetra_MultiVector(cmap, number_vectors));
    darcy_flux_ = Teuchos::rcp(new Epetra_Vector(fmap));

    copymemory_multivector(S.ref_total_component_concentration(), *total_component_concentration_);
    copymemory_vector(S.ref_darcy_flux(), *darcy_flux_);

  } else if (mode == ViewMemory) {
    porosity_ = S.porosity();
    water_saturation_ = S.water_saturation();
    prev_water_saturation_ = S.prev_water_saturation();
    water_density_ = S.water_density();
    mesh_ = S.mesh();

    double* data_df;
    double** data_tcc;
    const Epetra_Map& cmap = mesh_->cell_map(false);
    const Epetra_Map& fmap = mesh_->face_map(false);

    Epetra_Vector& df = S.ref_darcy_flux();
    df.ExtractView(&data_df);
    darcy_flux_ = Teuchos::rcp(new Epetra_Vector(View, fmap, data_df));

    Epetra_MultiVector & tcc = S.ref_total_component_concentration();
    tcc.ExtractView(&data_tcc);
    total_component_concentration_ = Teuchos::rcp(new Epetra_MultiVector(View, cmap, data_tcc, tcc.NumVectors()));
  }
  S_ = S.S_;
}


/* *******************************************************************
* Routine imports a short multivector to a parallel overlaping vector.
* No parallel communications are performed if target_is_parallel = 0.
******************************************************************* */
void Transport_State::copymemory_multivector(Epetra_MultiVector& source,
                                             Epetra_MultiVector& target,
                                             int target_is_parallel)
{
  const Epetra_BlockMap& source_cmap = source.Map();
  const Epetra_BlockMap& target_cmap = target.Map();

  int cmin, cmax, cmax_s, cmax_t;
  cmin = source_cmap.MinLID();
  cmax_s = source_cmap.MaxLID();
  cmax_t = target_cmap.MaxLID();
  cmax = std::min(cmax_s, cmax_t);

  int num_vectors = source.NumVectors();
  for (int c = cmin; c <= cmax; c++) {
    for (int i = 0; i < num_vectors; i++) target[i][c] = source[i][c];
  }

#ifdef HAVE_MPI
  if (target_is_parallel == 1) {
    if (cmax_s > cmax_t) {
      Errors::Message msg;
      msg << "The source map is bigger than the target map.\n";
      Exceptions::amanzi_throw(msg);
    }

    Epetra_Import importer(target_cmap, source_cmap);
    target.Import(source, importer, Insert);
  }
#endif
}


/* *******************************************************************
* Routine imports a short vector to a parallel overlaping vector.                
******************************************************************* */
void Transport_State::copymemory_vector(Epetra_Vector& source, Epetra_Vector& target)
{
  const Epetra_BlockMap& source_fmap = source.Map();
  const Epetra_BlockMap& target_fmap = target.Map();

  int fmin, fmax, fmax_s, fmax_t;
  fmin   = source_fmap.MinLID();
  fmax_s = source_fmap.MaxLID();
  fmax_t = target_fmap.MaxLID();
  fmax   = std::min(fmax_s, fmax_t);

  for (int f = fmin; f <= fmax; f++) target[f] = source[f];

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
* WARNING: MultiVector v must contain ghost cells.                
******************************************************************* */
void Transport_State::CopyMasterCell2GhostCell(Epetra_Vector& v)
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
* Copy data in ghost positions for a multivector.              
******************************************************************* */
void Transport_State::CopyMasterMultiCell2GhostMultiCell(Epetra_MultiVector& v)
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
* Interpolate linearly in time between two states v0 and v1. The time 
* is measuared relative to state v0; so that v1 is at time dT. The
* interpolated data are at time dT_int.            
******************************************************************* */
void Transport_State::InterpolateCellVector(
    const Epetra_Vector& v0, const Epetra_Vector& v1, double dT_int, double dT, Epetra_Vector& v_int)
{
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  double a = dT_int / dT;
  double b = 1.0 - a;
  for (int c = 0; c < ncells; c++) v_int[c] = v0[c] * b + v1[c] * a;
}


/* *******************************************************************
 * DEBUG: create constant analytical Darcy velocity fieldx u     
 ****************************************************************** */
void Transport_State::AnalyticDarcyFlux(const AmanziGeometry::Point& u)
{
  const Epetra_BlockMap& fmap = (*darcy_flux_).Map();

  for (int f = fmap.MinLID(); f <= fmap.MaxLID(); f++) {
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    (*darcy_flux_)[f] = u * normal;
  }
}
void Transport_State::AnalyticDarcyFlux(
    AmanziGeometry::Point f_vel(const AmanziGeometry::Point&, double), double t)
{
  const Epetra_BlockMap& fmap = (*darcy_flux_).Map();

  for (int f = fmap.MinLID(); f <= fmap.MaxLID(); f++) {
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    const AmanziGeometry::Point& fc = mesh_->face_centroid(f);
    (*darcy_flux_)[f] = f_vel(fc, t) * normal;
  }
}


/* *******************************************************************
 * DEBUG: create analytical concentration C = f(x, t)       
 ****************************************************************** */
void Transport_State::AnalyticTotalComponentConcentration(double f(const AmanziGeometry::Point&, double), double t)
{
  const Epetra_BlockMap& cmap = (*total_component_concentration_).Map();

  for (int c = cmap.MinLID(); c <= cmap.MaxLID(); c++) {
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
    (*total_component_concentration_)[0][c] = f(xc, t);
  }
}
void Transport_State::AnalyticTotalComponentConcentration(double tcc)
{
  const Epetra_BlockMap& cmap = (*total_component_concentration_).Map();

  for (int c = cmap.MinLID(); c <= cmap.MaxLID(); c++) {
    (*total_component_concentration_)[0][c] = tcc;
  }
}


/* **************************************************************** */
void Transport_State::error_total_component_concentration(
    double f(const AmanziGeometry::Point&, double), double t, double* L1, double* L2)
{
  int i, j, c;
  double d;
  const Epetra_BlockMap& cmap = (*total_component_concentration_).Map();

  *L1 = *L2 = 0.0;
  for (c = cmap.MinLID(); c <= cmap.MaxLID(); c++) {
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
    d = (*total_component_concentration_)[0][c] - f(xc, t);

    double volume = mesh_->cell_volume(c);
    *L1 += fabs(d) * volume;
    *L2 += d * d * volume;
  }

  *L2 = sqrt(*L2);
}


/* *******************************************************************
 * DEBUG: create constant analytical porosity                    
 ****************************************************************** */
void Transport_State::AnalyticPorosity(double phi)
{
  const Epetra_BlockMap& cmap = (*porosity_).Map();

  for (int c = cmap.MinLID(); c <= cmap.MaxLID(); c++) {
    (*porosity_)[c] = phi;  // default is 0.2
  }
}


/* ******************************************************************
 * DEBUG: create constant analytical water saturation            
 ***************************************************************** */
void Transport_State::AnalyticWaterSaturation(double ws)
{
  const Epetra_BlockMap& cmap = (*water_saturation_).Map();

  for (int c = cmap.MinLID(); c <= cmap.MaxLID(); c++) {
    (*water_saturation_)[c] = ws;  // default is 1.0
    (*prev_water_saturation_)[c] = ws;
  }
}


/* *****************************************************************
 * DEBUG: create constant analytical water density               
 **************************************************************** */
void Transport_State::AnalyticWaterDensity(double wd)
{
  const Epetra_BlockMap& cmap = (*water_density_).Map();

  for (int c = cmap.MinLID(); c <= cmap.MaxLID(); c++) {
    (*water_density_)[c] = wd;  // default is 1000.0
  }
}

}  // namespace AmanziTransport
}  // namespace Amanzi
