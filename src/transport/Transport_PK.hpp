#ifndef __Transport_PK_hpp__
#define __Transport_PK_hpp__

#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Teuchos_RCP.hpp"

#include "Mesh_maps_simple.hh"
#include "State.hpp"
#include "Transport_State.hpp"
#include "Transport_BCs.hpp"


/* 
   Amanzi: Demo 1: Transport Process Kernel Interface:

   The transport PK receives a reduced (optional) copy of 
   a physical state at time n and returns a different state 
   at time n+1. 

   Unmodified physical quantaties in the returned state are
   the smart pointers to the original variables.
*/


namespace Amanzi_Transport{
  const int TRANSPORT_NULL           = 0;
  const int TRANSPORT_STATE_BEGIN    = 1;
  const int TRANSPORT_STATE_COMPLETE = 2;

  const double TRANSPORT_LARGE_TIME_STEP = 1e+99;
  const double TRANSPORT_SMALL_TIME_STEP = 1e-12;
}



using namespace std;
using namespace Teuchos;
using namespace Amanzi_Transport;



class Transport_PK {

public:
  /* three constructors */
  Transport_PK( ParameterList &parameter_list_MPC,
		RCP<Transport_State> TS_MPC );
  Transport_PK();

  ~Transport_PK() {};

  /* primary members */
  double calculate_transport_dT();
  void advance();
  void commit_state( RCP<Transport_State> TS );

  void process_parameter_list();
  void identify_upwind_cells();
  void extract_darcy_flux();

  vector<double> calculate_accumulated_influx();
  vector<double> calculate_accumulated_outflux();

  void geometry_package();

  /* access members */ 
  RCP<Transport_State> get_transport_state()      const { return TS; }
  RCP<Transport_State> get_transport_state_next()       { return TS_next; }
  RCP<Transport_State> get_transport_state_next() const { return TS_next; }

  double get_transport_dT()      { return dT; }
  int    get_transport_status()  { return status; }
 

public:
  /* member for debugging only */
  double get_face_area( int f )   { return face_area[f]; }
  double get_cell_volume( int c ) { return cell_volume[c]; }


private:
  /* original and proposed transport states */
  RCP<Transport_State>  TS;
  RCP<Transport_State>  TS_next;
  
  /* darcy_flux */
  vector<double>  darcy_flux;

  /* parameter list with Transport specific parameters */
  ParameterList  parameter_list;

  /* part of the future geometry package */
  vector<double>  face_area;
  vector<double>  cell_volume;

  vector<double>  upwind_cell;
  vector<double>  downwind_cell;

  /* transport time step, CFL, and status */
  double  cfl, dT;
  int     number_components;
  int     status;

  /* boundary conditions for each components and each side set */
  /* it will be converted to a separate class                  */
  vector<Transport_BCs>  bcs;

  /* accumulated influx and outflux for each side */
  vector<double>   influx;
  vector<double>  outflux;

  /* frequently used data */
  int  cmin, cmax, number_cells;
  int  fmin, fmax, number_faces;
};

#endif
