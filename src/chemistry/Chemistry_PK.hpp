#ifndef __Chemistry_PK_hpp__
#define __Chemistry_PK_hpp__

#include "Teuchos_RCP.hpp"
#include "State.hpp"
#include "Chemistry_State.hpp"


// Chemistry Process Kernel Interface

class Chemistry_PK {

public:
  Chemistry_PK (Teuchos::ParameterList &parameter_list_, 
		Teuchos::RCP<Chemistry_State> CS_);

  ~Chemistry_PK ();

  void advance( );
  void commit_state ( Teuchos::RCP<Chemistry_State> );

private:

  // auxilary state for process kernel
  Teuchos::RCP<Chemistry_State> CS;

  // parameter list
  Teuchos::ParameterList parameter_list;

};

#endif
