#ifndef __Chemistry_PK_OperatorSplitting_hpp__
#define __Chemistry_PK_OperatorSplitting_hpp__

#include "Teuchos_RCP.hpp"
#include "Chemistry_PK.hpp"

#include "Geochemistry.hpp"



// Chemistry Process Kernel Interface for Operator Splitting
// 
// manages looping through each set of primary species provided by the
// Chemistry_State object and handing off the relevant data to the
// Geochemistry object which actually performs the calculations.

class Chemistry_PK_OperatorSplitting : public Chemistry_PK {

public:
  Chemistry_PK_OperatorSplitting(Teuchos::RCP<Chemistry_State> CS_);
  ~Chemistry_PK_OperatorSplitting();



private:

  Geochemistry geochem; // evalates geochemistry at a single node

};

#endif
