#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "UnitTest++.h"
#include "../Flow_BCs.hpp"

#include "Teuchos_Version.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"


TEST(Flow_BCs) {


  std::string    xmlInFileName = "test/test_Flow_BCs.xml";

  Teuchos::ParameterList flow_BCs_list;
  Teuchos::updateParametersFromXmlFile(xmlInFileName,&flow_BCs_list);



  Flow_BCs  fbcs ( flow_BCs_list );

  fbcs.read_Flow_BCs();
  
  std::vector<flow_bc> bc = fbcs.get_BCs();

  for (int i=0; i<bc.size(); i++) {
    std::cout << "type...     " << bc[i].bc_type << std::endl;
    std::cout << "value...    " << bc[i].value << std::endl;
    std::cout << "side set... " << bc[i].side_set << std::endl;
    std::cout << std::endl;
  }
  
  CHECK_EQUAL(bc.size(),fbcs.get_num_BCs());
}


