
#include "Flow_BCs.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

void Flow_BCs::read_Flow_BCs() {
  
  int num_BCs = paramList.get<int>("number of BCs", 0);

  BCs.resize(num_BCs);
  
  for (int i=0; i<num_BCs; i++) {
    
    std::string bc_name = "BC0 ";
    bc_name[3] = '1' + i;

    Teuchos::ParameterList bc_list = paramList.sublist(bc_name);
    
    BCs[i].side_set = bc_list.get<int>("Side set ID", 0);
    BCs[i].value = bc_list.get<double>("Dirichlet constant", 0.0);
    BCs[i].bc_type = bc_list.get<std::string>("Type","NOT DEFINED");    
  }
}
