
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_getConst.hpp"
#include "Teuchos_Version.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

#ifdef HAVE_TEUCHOS_EXTENDED
#include "Teuchos_XMLParameterListHelpers.hpp"
#endif

#ifdef USE_MPI
#include "mpi.h"
#endif

using Teuchos::CommandLineProcessor;
using Teuchos::ParameterList;
using Teuchos::getParameter;
typedef ParameterList::PrintOptions PLPrintOptions;
using Teuchos::ParameterEntry;
using Teuchos::OSTab;
using Teuchos::rcp;

void print_break() { std::cout << "---------------------------------------------------" << std::endl; }
double Plus ( double a, double b ) { return a+b; }


int main(int argc, char *argv[])
{

  using std::cout;
  using std::endl;

  bool verbose = true;
  bool result;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  const int procRank = Teuchos::GlobalMPISession::getRank();

  bool success = true;

  //-----------------------------------------------------------
  // Create Main Parameter List / Sublist Structure
  //-----------------------------------------------------------

  /* This creates the general structure of the input list */
  ParameterList MainPList("Main");
  ParameterList& MPCPList = MainPList.sublist("MPC");
  ParameterList& MeshPList = MainPList.sublist("Mesh");
  ParameterList& StatePList = MainPList.sublist("State");
  ParameterList& ChemPList = MainPList.sublist("Chemistry");
  ParameterList& FlowPList = MainPList.sublist("Flow");
  ParameterList& TransportPList = MainPList.sublist("Transport");


  /* Mesh */

  /* Each Mesh list MUST define a mesh class */
  MeshPList.get("Mesh Class", "MOAB");
 Teuchos::setStringToIntegralParameter<int>(
    "Mesh Class", "MOAB",
    "Defines the underlying mesh framework for the mesh class",
    Teuchos::tuple<std::string>("Simple","MOAB", "STK"),
    &MeshPList
    );

 
  ParameterList& MOABPList = MeshPList.sublist("MOAB Mesh Parameters");
  MOABPList.set("Exodus file name", "fbasin_unstr_025_V02_128P.h5m");

  /* Create a validator that accepts double or std::string input */
  typedef Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes AcceptedTypes;
  Teuchos::RCP<Teuchos::AnyNumberParameterEntryValidator>
    DoubleStringValidator = rcp(
      new Teuchos::AnyNumberParameterEntryValidator(
        Teuchos::AnyNumberParameterEntryValidator::PREFER_DOUBLE,
        AcceptedTypes(false).allowDouble(true).allowString(true)
        )
      );

  /* Create a validator that accepts int or std::string input */
  typedef Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes AcceptedTypes;
  Teuchos::RCP<Teuchos::AnyNumberParameterEntryValidator>
    IntStringValidator = rcp(
      new Teuchos::AnyNumberParameterEntryValidator(
        Teuchos::AnyNumberParameterEntryValidator::PREFER_INT,
        AcceptedTypes(false).allowInt(true).allowString(true)
        )
      );


  /* MPC */
  MPCPList.set("Start Time", 0.0, "Initial simulation time", DoubleStringValidator);
  MPCPList.set("End Time", 0.1, "Final simulation time", DoubleStringValidator);
  MPCPList.set("End Cycle", 1000, "Final cycle count", IntStringValidator);
  ParameterList& PKList = MPCPList.sublist("Enabled PK");
  PKList.get("Flow",true);
  PKList.get("Transport",true);
  PKList.get("Chemistry",false);

  MPCPList.get("Flow Model","Darcy");
  Teuchos::setStringToIntegralParameter<int>(
    "Flow Model", "Darcy",
    "Define the flow model",
    Teuchos::tuple<std::string>("Darcy","Richards"),
    &MPCPList
    );

  ParameterList& VizPList = MPCPList.sublist("Viz");
  ParameterList& CGNSPList = VizPList.sublist("CGNS");
  CGNSPList.set("Viz dump cycle frequency", 10);
  CGNSPList.set("Viz dump time frequency", 0.05);
  CGNSPList.set("File name", "test1.cgns");

  /* State */
  const Teuchos::Array<double> 
    gravity = Teuchos::tuple<double>(0.0,0.0,0.0);
  StatePList.set("Gravity",gravity);

  ParameterList& H2OPList = StatePList.sublist("Water");
  H2OPList.set("saturation",1.0);
  H2OPList.set("density",1000.0);
  H2OPList.set("viscosity",0.001308);


  /* Dump List to the screen */
  MainPList.print(cout,PLPrintOptions().showTypes(true).showDoc(true));

  return 0;
}    
