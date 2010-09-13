#ifndef __Geochemistry_hpp__
#define __Geochemistry_hpp__

#include "Teuchos_RCP.hpp"
#include "Chemistry_PK.hpp"

#include "Geochemistry.hpp"



// Chemistry Process Kernel Interface

class Chemistry_PK_OperatorSplitting : public Chemistry_PK {

public:
  Chemistry_PK_OperatorSplitting(Teuchos::RCP<Chemistry_State> CS_);
  ~Chemistry_PK_OperatorSplitting();



private:

  Geochemistry geochem; // evalates geochemistry at a single node

};

#endif
