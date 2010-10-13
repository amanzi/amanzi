
#include "Flow_BCs.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"
#include "limits.h"
#include "float.h"

void Flow_BCs::read_Flow_BCs() {
  
  num_BCs = paramList.get<int>("number of BCs", INT_MAX);
  // we cannot handle more than 100 BC paramerter lists

  if ( num_BCs > 100) throw std::exception(); 

  BCs.resize(num_BCs);

  for (int i=0; i<num_BCs; i++) {
    
    // by convention, boundary conditions have the name "BCXX"
    // where XX is a number from 00 to 99.
    std::string bc_name = "BC  ";
    
    int i1 = i%10;
    int i10 = i/10;

    bc_name[3] = '0' + i1;
    bc_name[2] = '0' + i10;

    if ( ! paramList.isSublist(bc_name) ) throw std::exception();
    Teuchos::ParameterList bc_list = paramList.sublist(bc_name);
    
    BCs[i].side_set = bc_list.get<int>("Side set ID", INT_MAX);
    if ( BCs[i].side_set == INT_MAX ) throw std::exception();
    
    BCs[i].value = bc_list.get<double>("BC value", DBL_MAX);
    if ( BCs[i].value == DBL_MAX ) throw std::exception();
    
    BCs[i].bc_type = bc_list.get<std::string>("Type","NOT DEFINED");    
    if ( BCs[i].bc_type == "NOT DEFINED") throw std::exception();
    
  }
}
