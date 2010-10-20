#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "State.hpp"
#include "Chemistry_State.hpp"
#include "Chemistry_PK.hpp"

class MPC {

public:
  MPC (Teuchos::RCP<Teuchos::ParameterList> Parameters_,
	    Teuchos::RCP<Mesh_maps_simple> mesh_maps_);
  ~MPC () {};

  void cycle_driver ();
  

private:
  
  // states
  Teuchos::RCP<State> S;
  Teuchos::RCP<Chemistry_State> CS;

  // misc setup information
  Teuchos::RCP<Teuchos::ParameterList> Parameters;
  Teuchos::RCP<Mesh_maps_simple> mesh_maps;

  // storage for chemistry's return value
  Teuchos::RCP<Epetra_MultiVector> total_component_concentration_star;

  // chemistry process kernel
  Teuchos::RCP<Chemistry_PK> CPK;

};


