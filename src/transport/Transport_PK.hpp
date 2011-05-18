#ifndef __Transport_PK_hpp__
#define __Transport_PK_hpp__

#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"

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

namespace Amanzi
{

namespace Amanzi_Transport{
  const int TRANSPORT_NULL           = 0;
  const int TRANSPORT_FLOW_AVAILABLE = 1;
  const int TRANSPORT_STATE_BEGIN    = 2;
  const int TRANSPORT_STATE_COMPLETE = 3;

  const double TRANSPORT_LARGE_TIME_STEP = 1e+99;
  const double TRANSPORT_SMALL_TIME_STEP = 1e-12;

  const int TRANSPORT_BC_CONSTANT_INFLUX = 1;
  const int TRANSPORT_BC_NULL = 2;

  const double TRANSPORT_CONCENTRATION_OVERSHOOT = 1e-6;

  const int TRANSPORT_AMANZI_VERSION = 1;
}



class Transport_PK {

public:
  /* three constructors */
  Transport_PK( Teuchos::ParameterList &parameter_list_MPC,
		Teuchos::RCP<Transport_State> TS_MPC );
  Transport_PK();

  ~Transport_PK() {};

  /* primary members */
  double calculate_transport_dT();
  void advance( double dT );
  void commit_state( Teuchos::RCP<Transport_State> TS );

  void process_parameter_list();
  void identify_upwind_cells();

  void check_divergence_property();
  void check_GEDproperty( Epetra_MultiVector & tracer ); 
  void print_statistics();

  std::vector<double>  calculate_accumulated_influx();
  std::vector<double>  calculate_accumulated_outflux();

  void geometry_package();

  /* access members */ 
  Teuchos::RCP<Transport_State>  get_transport_state()      const { return TS; }
  Teuchos::RCP<Transport_State>  get_transport_state_next() const { return TS_nextMPC; }

  double get_transport_dT()      { return dT; }
  double get_cfl()               { return cfl; }
  int    get_transport_status()  { return status; }
 
  Transport_State & ref_transport_state_next()  { return *TS_nextBIG; }


public:
  /* parallel information: will be moved to private */
  int  MyPID;

  /* output information */
  int  verbosity_level, internal_tests;
  double  tests_tolerance;

  /* member for debugging only */
  double get_face_area( int f )   { return face_area[f]; }
  double get_cell_volume( int c ) { return cell_volume[c]; }
  std::vector<double> & get_cell_volume() { return cell_volume; }


private:
  /* original and proposed (MPC and BIG) transport states */
  Teuchos::RCP<Transport_State>  TS;
  Teuchos::RCP<Transport_State>  TS_nextMPC;   /* uses memory of BIG */
  Teuchos::RCP<Transport_State>  TS_nextBIG; 
  
  /* parameter list with Transport specific parameters */
  Teuchos::ParameterList  parameter_list;

  /* part of the future geometry package */
  std::vector<double>  face_area;
  std::vector<double>  cell_volume;

  /* internal data */
  Teuchos::RCP<Epetra_IntVector>  upwind_cell;
  Teuchos::RCP<Epetra_IntVector>  downwind_cell;

  /* communication patterns */
  Teuchos::RCP<Epetra_Import>  cell_importer;
  Teuchos::RCP<Epetra_Import>  face_importer;

  /* transport time step, CFL, and status */
  double  cfl, dT, dT_debug;
  int     number_components;
  int     status;

  /* boundary conditions for each components and each side set */
  /* it will be converted to a separate class                  */
  std::vector<Transport_BCs>  bcs;

  /* frequently used data */
  int  cmin, cmax_owned, cmax, number_owned_cells, number_wghost_cells;
  int  fmin, fmax_owned, fmax, number_owned_faces, number_wghost_faces;
};

} // close namespace Amanzi

#endif
