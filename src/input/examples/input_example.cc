#include <iostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_ParameterListAcceptorHelpers.hpp"

#include "simple_mesh_input.hh"

typedef Teuchos::ParameterList::PrintOptions PLPrintOptions;

void print_line() { std::cout << "----------------------------------------" << std::endl; }

int main(int argc, char *argv[])
{

  using std::cout;
  using std::endl;
  using std::string;


  /* Create a input object to print out the defaults */
  SimpleMeshInput simple_mesh_info;
  print_line();
  cout << "Valid Parameters for Simple Mesh" << endl;
  Teuchos::printValidParameters(simple_mesh_info,cout,true);
  print_line();

  Teuchos::RCP<const Teuchos::ParameterList> info_list =
     simple_mesh_info.getParameterList(); 

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

  
  return 0;
}    
