#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"


#include "Chem_MPC.hpp"
#include "DataLayout.hpp"
#include "MeshWrapper.hpp"
#include "STKMesh1D.hpp"
#include "STK1DDisc.hpp"


#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
typedef int MPI_Comm;
#define MPI_COMM_WORLD 1
#include "Epetra_SerialComm.h"
#endif

int main ( int argc, char * argv[] )
{

  // initialize MPI
#ifdef HAVE_MPI
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  Teuchos::RCP<Epetra_Comm> Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  Teuchos::RCP<Epetra_Comm> Comm = Teuchos::rcp(new Epetra_SerialComm);
#endif
  
  // Command-line argument for input file
  char * xmlfilename=0;
  char defaultfile[10]={"input.xml"};
  if(argc>1){
    if(!strcmp(argv[1],"--help")){
      printf("cxx_main.exe [inputfile.xml]\n");
      exit(1);
    }
    else
      xmlfilename=argv[1];
  }
  else
    xmlfilename=defaultfile;
  

  // read the parameters from the input file
  Teuchos::RCP<Teuchos::ParameterList> Params = Teuchos::rcp(new Teuchos::ParameterList("Parameters"));
  Teuchos::updateParametersFromXmlFile(xmlfilename, Params.get());
  

  // create the mesh objects
  // this will change when we move to 3D

  // create the 1D mesh
  Teuchos::RCP<STKMesh1D> mesh = Teuchos::rcp(new STKMesh1D(Comm, Params));

  // create maps for the data that we will use in the discretization
  Teuchos::RCP<STK1DDisc> disc = Teuchos::rcp(new STK1DDisc(mesh, Comm));


  // initialize the Data Layout
  Teuchos::RCP<DataLayout> data_layout_1D 
    = Teuchos::rcp(new DataLayout( disc->getMap(), disc->getOverlapMap(),
   				   disc->getMap(), disc->getOverlapMap(), // in 1D face=vertex
  				   disc->getElementMap(), Teuchos::null));


  Teuchos::RCP<MeshWrapper> mesh_wrapper = Teuchos::rcp(new MeshWrapper(mesh, disc, data_layout_1D));

  // create the MPC object
  Teuchos::RCP<Chem_MPC> MPC = Teuchos::rcp(new Chem_MPC(Params, data_layout_1D, mesh_wrapper) );

  MPC->cycle_driver();

}
