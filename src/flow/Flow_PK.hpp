#ifndef __Flow_PK_hpp__
#define __Flow_PK_hpp__

#include "Teuchos_RCP.hpp"
#include "State.hpp"
#include "Flow_State.hpp"


// Flow Process Kernel Interface

class Flow_PK {

public:
  Flow_PK (Teuchos::RCP<Flow_State> CS_);

  ~Flow_PK ();

  void advance( Teuchos::RCP<Epetra_MultiVector> );
  void commit_state ( Teuchos::RCP<Flow_State> );

private:

  // auxilary state for process kernel
  Teuchos::RCP<Flow_State> CS;

};

#endif
