#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"

#include "State.hpp"
#include "cell_geometry.hh"
#include "Transport_State.hpp"


using namespace std;
using namespace Teuchos;


/* ************************************************************* */
/* at the moment a transport state is a copy of the global state */  
/* ************************************************************* */
Transport_State::Transport_State( State & S )
{
  total_component_concentration = S.get_total_component_concentration();
  porosity                      = S.get_porosity();
  darcy_flux                    = S.get_darcy_flux();
  water_saturation              = S.get_water_saturation();
  water_density                 = S.get_water_density();
  mesh_maps                     = S.get_mesh_maps();
}




/* ************************************************************* */
/* mode = CopyPointers (default0 means a trivial copy of the     */
/*                     given transport state                     */
/* mode = ViewMemory   creates the transport from internal one   */ 
/*                     as the MPC expected                       */
/* mode = CopyMemory   creates internal transport state based on */
/*                     ovelapped mesh maps                       */
/* ************************************************************* */
Transport_State::Transport_State( Transport_State & S, Transport_CreateMode mode )
{
  if ( mode == CopyPointers ) {
     total_component_concentration = S.get_total_component_concentration();
     porosity                      = S.get_porosity();
     darcy_flux                    = S.get_darcy_flux();
     water_saturation              = S.get_water_saturation();
     water_density                 = S.get_water_density();
     mesh_maps                     = S.get_mesh_maps();
  }

  else if ( mode == CopyMemory ) { 
     porosity         = S.get_porosity(); 
     water_saturation = S.get_water_saturation(); 
     water_density    = S.get_water_density();
     mesh_maps        = S.get_mesh_maps();

     /* allocate memory for internal state */
     const Epetra_Map &  cmap = mesh_maps->cell_map( true );
     const Epetra_Map &  fmap = mesh_maps->face_map( true );

     int number_vectors = S.get_total_component_concentration()->NumVectors();

     total_component_concentration = rcp( new Epetra_MultiVector( cmap, number_vectors ) );
     darcy_flux = rcp( new Epetra_Vector( fmap ) );

     copymemory_multivector( S.ref_total_component_concentration(), *total_component_concentration );
     copymemory_vector( S.ref_darcy_flux(), *darcy_flux );
  }

  else if ( mode == ViewMemory ) {
     porosity         = S.get_porosity(); 
     water_saturation = S.get_water_saturation(); 
     water_density    = S.get_water_density();
     mesh_maps        = S.get_mesh_maps();

     double*  data_df;
     double**  data_tcc;
     const Epetra_Map &  cmap = mesh_maps->cell_map( false );
     const Epetra_Map &  fmap = mesh_maps->face_map( false );

     Epetra_Vector & df = S.ref_darcy_flux();
     df.ExtractView( &data_df );     
     darcy_flux = rcp( new Epetra_Vector( View, fmap, data_df ) );

     Epetra_MultiVector & tcc = S.ref_total_component_concentration();
     tcc.ExtractView( &data_tcc );     
     total_component_concentration = rcp( new Epetra_MultiVector( View, cmap, data_tcc, tcc.NumVectors() ) );
  }
}




/* ************************************************************* */
/* import concentrations to internal Transport state             */
/* ************************************************************* */
void Transport_State::copymemory_multivector( Epetra_MultiVector & source, Epetra_MultiVector & target )
{
  int  i, c, cmin, cmax, cmax_s, cmax_t, number_vectors;

  const Epetra_BlockMap &  source_cmap = source.Map();
  const Epetra_BlockMap &  target_cmap = target.Map();

  cmin   = source_cmap.MinLID();
  cmax_s = source_cmap.MaxLID();
  cmax_t = target_cmap.MaxLID();
  cmax   = min( cmax_s, cmax_t );

  number_vectors = source.NumVectors();
  for( c=cmin; c<=cmax; c++ ) 
     for( i=0; i<number_vectors; i++ ) target[i][c] = source[i][c];

#ifdef HAVE_MPI
  if ( cmax_s > cmax_t ) throw exception();

  Epetra_Import  importer( target_cmap, source_cmap );
  target.Import( source, importer, Insert );
#endif
}




/* ************************************************************* */
/* import darcy flux to internal Transport state                 */
/* ************************************************************* */
void Transport_State::copymemory_vector( Epetra_Vector & source, Epetra_Vector & target )
{
  int  f, fmin, fmax, fmax_s, fmax_t;

  const Epetra_BlockMap &  source_fmap = source.Map();
  const Epetra_BlockMap &  target_fmap = target.Map();

  fmin   = source_fmap.MinLID();
  fmax_s = source_fmap.MaxLID();
  fmax_t = target_fmap.MaxLID();
  fmax   = min( fmax_s, fmax_t );

  for( f=fmin; f<=fmax; f++ ) target[f] = source[f];

#ifdef HAVE_MPI
  if ( fmax_s > fmax_t ) throw exception(); 

  Epetra_Import  Importer( target_fmap, source_fmap );
  target.Import( source, Importer, Insert );
#endif
}




/* ************************************************************* */
/* DEBUG: create constant analytical Darcy velocity fieldx u     */
/* ************************************************************* */
void Transport_State::analytic_darcy_flux( double* u )
{
  int  i, f;
  double x[4][3], normal[3];

  const Epetra_BlockMap &  fmap = (*darcy_flux).Map();

  for( f=fmap.MinLID(); f<=fmap.MaxLID(); f++ ) { 
     mesh_maps->face_to_coordinates( f, (double*) x, (double*) x+12 );

     quad_face_normal(x[0], x[1], x[2], x[3], normal);

     (*darcy_flux)[f] = u[0] * normal[0] + u[1] * normal[1] + u[2] * normal[2];
  }
}




/* ************************************************************* */
/* DEBUG: create constant analytical concentration C_0 = x       */
/* ************************************************************* */
void Transport_State::analytic_total_component_concentration( double f(double*, double), double t )
{
  int  i, j, c;
  double x[8][3], center[3];

  const Epetra_BlockMap &  cmap = (*total_component_concentration).Map();

  for( c=cmap.MinLID(); c<=cmap.MaxLID(); c++ ) { 
     mesh_maps->cell_to_coordinates( c, (double*) x, (double*) x+24);

     for( i=0; i<3; i++ ) { 
        center[i] = 0;
        for( j=0; j<8; j++ ) center[i] += x[j][i];
        center[i] /= 8;
     }

     (*total_component_concentration)[0][c] = f( center, t );
  }
}




void Transport_State::error_total_component_concentration( double f(double*, double), double t, vector<double> & cell_volume, double* L1, double* L2 )
{
  int  i, j, c;
  double  d, x[8][3], center[3];

  const Epetra_BlockMap &  cmap = (*total_component_concentration).Map();

  *L1 = *L2 = 0.0;
  for( c=cmap.MinLID(); c<=cmap.MaxLID(); c++ ) { 
     mesh_maps->cell_to_coordinates( c, (double*) x, (double*) x+24);

     for( i=0; i<3; i++ ) { 
        center[i] = 0;
        for( j=0; j<8; j++ ) center[i] += x[j][i];
        center[i] /= 8;
     }

     d = (*total_component_concentration)[0][c] - f( center, t ); 

     *L1 += fabs( d ) * cell_volume[c];
     *L2 += d * d * cell_volume[c];
  }

  *L2 = sqrt( *L2 );
}




/* ************************************************************* */
/* DEBUG: create constant analytical porosity                    */
/* ************************************************************* */
void Transport_State::analytic_porosity( double phi )
{
  int  c;
  const Epetra_BlockMap &  cmap = (*porosity).Map();

  for( c=cmap.MinLID(); c<=cmap.MaxLID(); c++ ) { 
     (*porosity)[c] = phi;  /* default is 0.2 */
  }
}




/* ************************************************************* */
/* DEBUG: create constant analytical water saturation            */
/* ************************************************************* */
void Transport_State::analytic_water_saturation( double ws )
{
  int  c;
  const Epetra_BlockMap &  cmap = (*water_saturation).Map();

  for( c=cmap.MinLID(); c<=cmap.MaxLID(); c++ ) { 
     (*water_saturation)[c] = ws;  /* default is 1.0 */
  }
}




/* ************************************************************* */
/* DEBUG: create constant analytical water density               */
/* ************************************************************* */
void Transport_State::analytic_water_density( double wd )
{
  int  c;
  const Epetra_BlockMap &  cmap = (*water_density).Map();

  for( c=cmap.MinLID(); c<=cmap.MaxLID(); c++ ) { 
     (*water_density)[c] = wd;  /* default is 1000.0 */
  }
}
