#ifndef __Transport_PK_hpp__
#define __Transport_PK_hpp__

#include "Teuchos_RCP.hpp"
#include "State.hpp"
#include "Transport_State.hpp"


/* 
   Amanzi: Demo 1: Transport Process Kernel Interface:

   The transport PK receives a reduced (optional) copy of 
   a physical state at time n and returns a different state 
   at time n+1. 

   Unmodified physical quantaties in the returned state are
   the smart pointers to the original variables.
*/


class Transport_PK {

public:
  Transport_PK ( Teuchos::RCP<Transport_State> TS_MPC );

  ~Transport_PK ();

  /* primary members */
  void calculate_transport_dT ();
  void advance_transport_state ();

  /* access members */ 
  Teuchos::RCP<Transport_State> get_transport_state ()      const { return TS; }
  Teuchos::RCP<Transport_State> get_transport_state_next () const { return TS_next; }

  double get_transport_dT ()     { return dT; }
  int    get_transport_status () { return status; }


private:
  /* smart pointer to the transport state for process kernel */
  Teuchos::RCP<Transport_State> TS;

  /* proposed new transport state */ 
  Teuchos::RCP<Transport_State> TS_next;

  /* part of the future geometry package */
  Epetra_Vector face_area();

  /* transport time step and status */
  double dT;
  int    status;
};

#endif
