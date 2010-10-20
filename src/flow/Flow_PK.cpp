
#include "Flow_PK.hpp"
#include "Mesh_maps_simple.hh"

#include "Epetra_MultiVector.h"
#include "Epetra_MpiComm.h"
#include "NOX.H"
#include "NOX_Utils.H"
#include "NOX_GlobalData.H"
#include "NOX_Epetra_Vector.H"
#include "NOX_Epetra_Group.H"
#include <ml_MultiLevelPreconditioner.h>
#include "NOX_Epetra_LinearSystem_AztecOO.H"

#include "Darcy_problem.hpp"
#include "Problem_interface.hpp"




Flow_PK::Flow_PK( Teuchos::RCP<Flow_State> FS_):
  FS(FS_)
{ 
  // use the constructor to initialize the flow process kernel


  Darcy_problem DP(FS);

  Epetra_Vector InitialGuess(*DP.get_NL_map());
  InitialGuess.PutScalar(0.0);

  // Set up the problem interface
  Teuchos::RCP<Problem_interface> interface = 
    Teuchos::rcp(new Problem_interface(&DP) );

  // Create the top level parameter list
  Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr =
    Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& nlParams = *(nlParamsPtr.get());

  // Set the nonlinear solver method
  nlParams.set("Nonlinear Solver", "Line Search Based");

  // Set the printing parameters in the "Printing" sublist
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  // printParams.set("MyPID", Comm.MyPID()); 
  printParams.set("Output Precision", 3);
  printParams.set("Output Processor", 0);
  printParams.set("Output Information", 
                        NOX::Utils::OuterIteration + 
                        NOX::Utils::OuterIterationStatusTest + 
                        NOX::Utils::InnerIteration +
                        NOX::Utils::Parameters + 
                        NOX::Utils::Details + 
                        NOX::Utils::Warning);


  // start definition of nonlinear solver parameters
  // Sublist for line search 
  Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
  searchParams.set("Method","Full Step");
  
  // Sublist for direction
  // here we define our own direction method
  Teuchos::ParameterList &dirParams = nlParams.sublist("Direction");
  dirParams.set("Method", "User Defined");

  // this is probably wrong... FIXME
  Teuchos::ParameterList& lsParams = nlParams.sublist("Linear Solver");
  lsParams.set("Preconditioner","User Defined");


  // we need a GlobalData variable fo initialize the 
  // Direction factory
  Teuchos::RCP<NOX::GlobalData> gd = Teuchos::rcp(new NOX::GlobalData(nlParamsPtr));

  Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = interface;
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;
  Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = interface;

  // Need a NOX::Epetra::Vector for constructor
  NOX::Epetra::Vector noxInitGuess(InitialGuess, NOX::DeepCopy);
  
  Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> MLPrec = 
    Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*DP.get_PrecMat(), true) );

  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys = 
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
						      iReq, iPrec, MLPrec, InitialGuess)); 


  Teuchos::RCP<NOX::Epetra::Group> grpPtr = 
    Teuchos::rcp(new NOX::Epetra::Group(printParams, 
                                        iReq, 
                                        noxInitGuess, 
                                        linSys)); 
  // Set up the status tests
  Teuchos::RCP<NOX::StatusTest::NormF> testNormF = 
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-6));
  Teuchos::RCP<NOX::StatusTest::MaxIters> testMaxIters = 
    Teuchos::rcp(new NOX::StatusTest::MaxIters(100));
  // this will be the convergence test to be used
  Teuchos::RCP<NOX::StatusTest::Combo> combo = 
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR, 
                                            testNormF, testMaxIters));
  
  // Create the solver
  solver = NOX::Solver::buildSolver(grpPtr, combo, nlParamsPtr);



};

Flow_PK::~Flow_PK()
{ 

};


void Flow_PK::advance( Teuchos::RCP<Epetra_MultiVector> tcc_star )
{
  cout << "advancing the state of the flow process model here" << endl;

  // this is how to get the element volumes...
  // mesh_wrapper->get_element_volumes() 

  // this is how to get the total component concentration
  // CS->get_total_component_concentration()

  // please update the argument to this function called tcc_star
  // with the result of your flow computation which is
  // the total component concentration ^star

  // see the Flow_State for the other available 
  // data in the flow specific state



  // Solve the nonlinear system
  NOX::StatusTest::StatusType status = solver->solve();

  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group & finalGroup = 
    dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
  const Epetra_Vector & finalSolution = 
    (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

  

};


void Flow_PK::commit_state( Teuchos::RCP<Flow_State> )
{
  cout << "committing the internal state of the flow process model" << endl;

  // use this function to commit the internal state
  // of the flow process kernel


  
  

};


