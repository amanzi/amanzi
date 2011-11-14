#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"

#include "Point.hh"

#include "State.hpp"
#include "Transport_State.hpp"

using namespace Teuchos;

namespace Amanzi {
namespace AmanziTransport {

/* **************************************************************** */
Transport_State::Transport_State(State& S)
{
  total_component_concentration = S.get_total_component_concentration();
  porosity = S.get_porosity();
  darcy_flux = S.get_darcy_flux();
  water_saturation = S.get_water_saturation();
  water_density = S.get_water_density();
  mesh_maps = S.get_mesh_maps();
}


/* *******************************************************************
 * mode = CopyPointers (default) a trivial copy of the given state           
 * mode = ViewMemory   creates the transport from internal one   
 *                     as the MPC expected                     
 * mode = CopyMemory   creates internal transport state based on 
 *                     ovelapped mesh maps                       
/* **************************************************************** */
Transport_State::Transport_State(Transport_State& S, TransportCreateMode mode)
{
  if (mode == CopyPointers) {
    total_component_concentration = S.get_total_component_concentration();
    porosity = S.get_porosity();
    darcy_flux = S.get_darcy_flux();
    water_saturation = S.get_water_saturation();
    water_density = S.get_water_density();
    mesh_maps = S.get_mesh_maps();
  }
  else if (mode == CopyMemory ) { 
    porosity = S.get_porosity(); 
    water_saturation = S.get_water_saturation(); 
    water_density = S.get_water_density();
    mesh_maps = S.get_mesh_maps();

    // allocate memory for internal state
    const Epetra_Map& cmap = mesh_maps->cell_map(true);
    const Epetra_Map& fmap = mesh_maps->face_map(true);

    int number_vectors = S.get_total_component_concentration()->NumVectors();

    total_component_concentration = rcp(new Epetra_MultiVector(cmap, number_vectors));
    darcy_flux = rcp(new Epetra_Vector(fmap));

    copymemory_multivector(S.ref_total_component_concentration(), *total_component_concentration);
    copymemory_vector(S.ref_darcy_flux(), *darcy_flux);
  }

  else if (mode == ViewMemory) {
    porosity = S.get_porosity(); 
    water_saturation = S.get_water_saturation(); 
    water_density = S.get_water_density();
    mesh_maps = S.get_mesh_maps();

    double* data_df;
    double** data_tcc;
    const Epetra_Map& cmap = mesh_maps->cell_map(false);
    const Epetra_Map& fmap = mesh_maps->face_map(false);

    Epetra_Vector& df = S.ref_darcy_flux();
    df.ExtractView(&data_df);     
    darcy_flux = rcp(new Epetra_Vector(View, fmap, data_df));

    Epetra_MultiVector & tcc = S.ref_total_component_concentration();
    tcc.ExtractView(&data_tcc);     
    total_component_concentration = rcp(new Epetra_MultiVector(View, cmap, data_tcc, tcc.NumVectors()));
  }
}


/* *******************************************************************
 * import concentrations to internal Transport state             
 ****************************************************************** */
void Transport_State::copymemory_multivector(Epetra_MultiVector& source, 
                                             Epetra_MultiVector& target)
{
  const Epetra_BlockMap& source_cmap = source.Map();
  const Epetra_BlockMap& target_cmap = target.Map();

  int cmin, cmax, cmax_s, cmax_t;
  cmin   = source_cmap.MinLID();
  cmax_s = source_cmap.MaxLID();
  cmax_t = target_cmap.MaxLID();
  cmax   = std::min(cmax_s, cmax_t);

  int number_vectors = source.NumVectors();
  for (int c=cmin; c<=cmax; c++) {
     for (int i=0; i<number_vectors; i++) target[i][c] = source[i][c];
  }

#ifdef HAVE_MPI
  if (cmax_s > cmax_t) throw std::exception();  // must be replaced by amanziException

  Epetra_Import importer(target_cmap, source_cmap);
  target.Import(source, importer, Insert);
#endif
}


/* *******************************************************************
 * Routine imports Darcy flux to internal Transport state                 
 ****************************************************************** */
void Transport_State::copymemory_vector(Epetra_Vector& source, Epetra_Vector& target)
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
  if (fmax_s > fmax_t) throw std::exception(); 

  Epetra_Import importer(target_fmap, source_fmap);
  target.Import(source, importer, Insert);
#endif
}


/* *******************************************************************
 * DEBUG: create constant analytical Darcy velocity fieldx u     
 ****************************************************************** */
void Transport_State::analytic_darcy_flux(const AmanziGeometry::Point& u)
{
  const Epetra_BlockMap& fmap = (*darcy_flux).Map();

  for (int f=fmap.MinLID(); f<=fmap.MaxLID(); f++) { 
    const AmanziGeometry::Point& normal = mesh_maps->face_normal(f);    
    (*darcy_flux)[f] = u * normal;
  }
}
void Transport_State::analytic_darcy_flux(
    AmanziGeometry::Point f_vel(const AmanziGeometry::Point&, double), double t)
{
  const Epetra_BlockMap& fmap = (*darcy_flux).Map();

  for (int f=fmap.MinLID(); f<=fmap.MaxLID(); f++) { 
    const AmanziGeometry::Point& normal = mesh_maps->face_normal(f);
    const AmanziGeometry::Point& fc = mesh_maps->face_centroid(f);
    (*darcy_flux)[f] = f_vel(fc, t) * normal;
  }
}


/* *******************************************************************
 * DEBUG: create analytical concentration C = f(x, t)       
 ****************************************************************** */
void Transport_State::analytic_total_component_concentration(double f(const AmanziGeometry::Point&, double), double t)
{
  const Epetra_BlockMap& cmap = (*total_component_concentration).Map();

  for (int c=cmap.MinLID(); c<=cmap.MaxLID(); c++) { 
    const AmanziGeometry::Point& xc = mesh_maps->cell_centroid(c);    
    (*total_component_concentration)[0][c] = f(xc, t);
  }
}
void Transport_State::analytic_total_component_concentration(double tcc)
{
  const Epetra_BlockMap& cmap = (*total_component_concentration).Map();

  for (int c=cmap.MinLID(); c<=cmap.MaxLID(); c++) { 
    (*total_component_concentration)[0][c] = tcc;
  }
}


/* **************************************************************** */
void Transport_State::error_total_component_concentration(
    double f(const AmanziGeometry::Point&, double), double t, double* L1, double* L2)
{
  int i, j, c;
  double d;
  const Epetra_BlockMap& cmap = (*total_component_concentration).Map();

  *L1 = *L2 = 0.0;
  for (c=cmap.MinLID(); c<=cmap.MaxLID(); c++ ) { 
    const AmanziGeometry::Point& xc = mesh_maps->cell_centroid(c);
    d = (*total_component_concentration)[0][c] - f(xc, t); 

    double volume = mesh_maps->cell_volume(c);
    *L1 += fabs(d) * volume;
    *L2 += d * d * volume;
  }

  *L2 = sqrt( *L2 );
}


/* *******************************************************************
 * DEBUG: create constant analytical porosity                    
 ****************************************************************** */
void Transport_State::analytic_porosity(double phi)
{
  const Epetra_BlockMap& cmap = (*porosity).Map();

  for (int c=cmap.MinLID(); c<=cmap.MaxLID(); c++) { 
    (*porosity)[c] = phi;  // default is 0.2
  }
}


/* ******************************************************************
 * DEBUG: create constant analytical water saturation            
 ***************************************************************** */
void Transport_State::analytic_water_saturation(double ws)
{
  const Epetra_BlockMap &  cmap = (*water_saturation).Map();

  for (int c=cmap.MinLID(); c<=cmap.MaxLID(); c++) { 
    (*water_saturation)[c] = ws;  // default is 1.0 
  }
}


/* *****************************************************************
 * DEBUG: create constant analytical water density               
 **************************************************************** */
void Transport_State::analytic_water_density(double wd)
{
  const Epetra_BlockMap &  cmap = (*water_density).Map();

  for (int c=cmap.MinLID(); c<=cmap.MaxLID(); c++) { 
    (*water_density)[c] = wd;  // default is 1000.0
  }
}

}  // namespace AmanziTransport
}  // namespace Amanzi
