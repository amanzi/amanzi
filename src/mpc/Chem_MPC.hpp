#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "State.hpp"
#include "Chemistry_State.hpp"
#include "Chemistry_PK.hpp"

class Chem_MPC {

public:
  Chem_MPC (Teuchos::RCP<Teuchos::ParameterList> Parameters_,
	    Teuchos::RCP<STK_mesh::Mesh_maps_stk> mesh_maps_);
  ~Chem_MPC () {};

  void cycle_driver ();
  

private:
  
  // states
  Teuchos::RCP<State> S;
  Teuchos::RCP<Chemistry_State> CS;

  // misc setup information
  Teuchos::RCP<Teuchos::ParameterList> Parameters;
  Teuchos::RCP<STK_mesh::Mesh_maps_stk> mesh_maps;

  // storage for chemistry's return value
  Teuchos::RCP<Epetra_MultiVector> total_component_concentration_star;

  // chemistry process kernel
  Teuchos::RCP<Chemistry_PK> CPK;

};


