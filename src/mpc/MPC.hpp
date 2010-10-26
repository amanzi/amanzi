#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "State.hpp"
#include "Chemistry_State.hpp"
#include "Chemistry_PK.hpp"
#include "Transport_State.hpp"
#include "Transport_PK.hpp"
#include "Flow_State.hpp"
#include "Flow_PK.hpp"

class MPC {

public:
  MPC (Teuchos::ParameterList parameter_list_,
       Teuchos::RCP<Mesh_maps_base> mesh_maps_);
  ~MPC () {};

  void cycle_driver ();
  void write_mesh();

private:
  
  // states
  Teuchos::RCP<State> S;
  Teuchos::RCP<Chemistry_State> CS;
  Teuchos::RCP<Transport_State> TS; 
  Teuchos::RCP<Flow_State> FS;

  // misc setup information
  Teuchos::ParameterList parameter_list;
  Teuchos::RCP<Mesh_maps_base> mesh_maps;

  // storage for chemistry's return value
  Teuchos::RCP<Epetra_MultiVector> total_component_concentration_star;

  // process kernels
  Teuchos::RCP<Chemistry_PK> CPK;
  Teuchos::RCP<Transport_PK> TPK;
  Teuchos::RCP<Flow_PK> FPK; 


};


