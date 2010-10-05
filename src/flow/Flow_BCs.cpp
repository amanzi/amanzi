
#include "Flow_BCs.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

void Flow_BCs::read_Flow_BCs() {
  
  int num_BCs = paramList.get<int>("number of BCs", 0);


  Teuchos::ParameterList bc_list = paramList.sublist("BC01");
  bc_list.print(std::cout,2,true,true);



  std::cout << num_BCs << std::endl;

  // Teuchos::Array<int> SSIDs;
  // SSIDs = Teuchos::getArrayFromStringParameter(bc_list, "Side Set IDs");
  

}
