#include <iostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_ParameterListAcceptorHelpers.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"


#include "simple_mesh_input.hh"
#include "transport_input.hh"
#include "test_input.hh"

typedef Teuchos::ParameterList::PrintOptions PLPrintOptions;

void print_line() { std::cout << "----------------------------------------" << std::endl; }

int main(int argc, char *argv[])
{

  using std::cout;
  using std::endl;
  using std::string;

  /* Create a simple test input object */
  TestInput testme;
  print_line();
  cout << "Test Input Object (Default)" << endl;
  Teuchos::printValidParameters(testme,cout,true);
  print_line();

  /* Add a simple list to testme */
  testme.add_sublist(5,1.23456e-9,true);
  testme.add_sublist(23,1.0045,false);
  Teuchos::RCP<const Teuchos::ParameterList> testme_list =
      testme.getParameterList();
  print_line();
  cout << "Dump Test Input List information to XML format" << endl;
  Teuchos::writeParameterListToXmlOStream(*testme_list,std::cout); 
  print_line();





  /* Create a input object to print out the defaults */
  SimpleMeshInput simple_mesh_info;
  print_line();
  cout << "Valid Parameters for Simple Mesh" << endl;
  print_line();
  Teuchos::printValidParameters(simple_mesh_info,cout,true);
  print_line();

  /* Example of creating Simple input object */
  print_line();
  cout << "Sample simple mesh input setup" << endl;
  print_line();
  SimpleMeshInput simple_mesh(0.0,0.0,0.0,1.0,1.0,1.0,30,30,30);
  simple_mesh.add_block(0.0,0.1);
  simple_mesh.add_block(0.1,0.5);
  simple_mesh.add_block(0.5,1.0);

  Teuchos::RCP<const Teuchos::ParameterList> simple_list =
      simple_mesh.getParameterList();

  simple_list->print(std::cout,PLPrintOptions().showTypes(true).showDoc(true));
  print_line();
  cout << "Dump information to XML format" << endl;
  Teuchos::writeParameterListToXmlOStream(*simple_list,std::cout); 
  print_line();

  /* Transport example */
  TransportInput transport;

  print_line();
  cout << "Transport Default Parameters (XML)" << endl;
  print_line();
  string transport_xml = transport.getXmlString(); 
  cout << transport_xml << endl;
  print_line();

  print_line();
  cout << "Change transport parameters CFL and max dT" << endl;
  print_line();
  transport.set_cfl(0.5);
  transport.set_max_dt(1.0e-03);
  transport.set_verbosity(1);
  transport.getParameterList()->print(cout,PLPrintOptions().showTypes(true).showDoc(true));
  print_line();


  return 0;
}    
