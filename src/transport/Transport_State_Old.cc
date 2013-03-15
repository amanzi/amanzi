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

#include "State_Old.hh"
#include "Transport_State.hh"

namespace Amanzi {
namespace AmanziTransport {

/* *******************************************************************
* Create Flow state from a state.                     
******************************************************************* */
Transport_State::Transport_State(State_Old& S)
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

    CopyMasterMultiCell2GhostMultiCell(
        S.ref_total_component_concentration(), *total_component_concentration_);
    CopyMasterFace2GhostFace(S.ref_darcy_flux(), *darcy_flux_);

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
* Copy cell-based data from master to ghost positions.              
* WARNING: Vector v must contain ghost cells.                
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
* Routine imports a short vector to a parallel overlapping vector.                
******************************************************************* */
void Transport_State::CopyMasterCell2GhostCell(const Epetra_Vector& v, Epetra_Vector& vv)
{
#ifdef HAVE_MPI
  const Epetra_BlockMap& source_cmap = mesh_->cell_map(false);
  const Epetra_BlockMap& target_cmap = mesh_->cell_map(true);
  Epetra_Import importer(target_cmap, source_cmap);

  double* vdata;
  v.ExtractView(&vdata);
  Epetra_Vector vcells(View, source_cmap, vdata);

  vv.Import(vcells, importer, Insert);
#else
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells_owned; c++) vv[c] = v[c];
#endif
}


/* *******************************************************************
* Copy cell-based data from master to ghost positions.              
* WARNING: MultiVector v must contain ghost cells.                
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
* Routine imports a short multivector to a parallel overlapping multivector.                
******************************************************************* */
void Transport_State::CopyMasterMultiCell2GhostMultiCell(const Epetra_MultiVector& v, 
                                                         Epetra_MultiVector& vv,
                                                         int parallel_comm)
{
#ifdef HAVE_MPI
  if (parallel_comm == 1) {
    const Epetra_BlockMap& source_cmap = mesh_->cell_map(false);
    const Epetra_BlockMap& target_cmap = mesh_->cell_map(true);
    Epetra_Import importer(target_cmap, source_cmap);

    vv.Import(v, importer, Insert);
  } else {
    int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    int num_vectors = v.NumVectors();
    for (int c = 0; c < ncells_owned; c++) {
      for (int i = 0; i < num_vectors; i++) vv[i][c] = v[i][c];
    }
  }
#else
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int num_vectors = v.NumVectors();
  for (int c = 0; c < ncells_owned; c++) {
    for (int i = 0; i < num_vectors; i++) vv[i][c] = v[i][c];
  }
#endif
}


/* *******************************************************************
* Routine imports a short vector to a parallel overlapping vector.                
******************************************************************* */
void Transport_State::CopyMasterFace2GhostFace(const Epetra_Vector& v, Epetra_Vector& vv)
{
#ifdef HAVE_MPI
  const Epetra_BlockMap& source_cmap = mesh_->face_map(false);
  const Epetra_BlockMap& target_cmap = mesh_->face_map(true);
  Epetra_Import importer(target_cmap, source_cmap);

  double* vdata;
  v.ExtractView(&vdata);
  Epetra_Vector vcells(View, source_cmap, vdata);

  vv.Import(vcells, importer, Insert);
#else
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f = 0; f < nfaces_owned; f++) vv[f] = v[f];
#endif
}


/* *******************************************************************
* Calculate minimum values in a multivector using master cells.            
******************************************************************* */
void Transport_State::MinValueMasterCells(Epetra_MultiVector& v, double* vmin)
{
#ifdef HAVE_MPI
  const Epetra_BlockMap& cmap = mesh_->cell_map(false);

  double** vdata;
  v.ExtractView(&vdata);
  Epetra_MultiVector vv(View, cmap, vdata, v.NumVectors());
  
  vv.MinValue(vmin);
#else
  v.MinValues(vmin);
#endif
}


/* *******************************************************************
* Calculate miximum values in a multivector using master cells.             
******************************************************************* */
void Transport_State::MaxValueMasterCells(Epetra_MultiVector& v, double* vmax)
{
#ifdef HAVE_MPI
  const Epetra_BlockMap& cmap = mesh_->cell_map(false);

  double** vdata;
  v.ExtractView(&vdata);
  Epetra_MultiVector vv(View, cmap, vdata, v.NumVectors());
  
  vv.MaxValue(vmax);
#else
  v.MaxValues(vmax);
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
