#ifndef __FlowBCs_hpp__
#define __FlowBCs_hpp__

#include "Teuchos_ParameterList.hpp"

// we assume that the BCs have already read into
// the internal Teuchos parameter list database
// 
// we also assume that the flow BCs are named
// FlowBC01 ... FlowBCXX, where XX is the number
// of flow BCs (really, this is the number of 
// side sets

class Flow_BCs {

public:
  Flow_BCs(Teuchos::ParameterList &parameterList):
    paramList(parameterList) {};
  ~Flow_BCs() {};

  void read_Flow_BCs();
  
private:
  Teuchos::ParameterList paramList; 

  
};

#endif
