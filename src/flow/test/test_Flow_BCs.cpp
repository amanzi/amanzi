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



  Flow_BCs  fbcs ( flow_BCs_list);

  fbcs.read_Flow_BCs();
  
}
