#ifndef __Transport_PK_hpp__
#define __Transport_PK_hpp__

#include "Teuchos_RCP.hpp"
#include "State.hpp"
#include "Transport_State.hpp"


// Transport Process Kernel Interface

class Transport_PK {

public:
  Transport_PK (Teuchos::RCP<Transport_State> TS_);

  ~Transport_PK ();

  void advance();
  void commit_state ( Teuchos::RCP<Transport_State> );
  double time_step_estimate();

private:

  // auxilary state for process kernel
  const Teuchos::RCP<Transport_State> TS_mirror;
  Teuchos::RCP<Transport_State> TS_copy;

};

#endif
