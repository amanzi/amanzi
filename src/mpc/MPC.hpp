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
  //void write_mesh();
  void write_mesh_data(std::string gmv_meshfile, std::string gmv_datafile, 
		       const int iter, const int digits);
#ifdef ENABLE_CGNS
  void write_cgns_data(std::string filename, int iter);
#endif
  void read_parameter_list();

private:
  
  // states
  Teuchos::RCP<State> S;
  Teuchos::RCP<Chemistry_State> CS;
  Teuchos::RCP<Transport_State> TS; 
  Teuchos::RCP<Flow_State> FS;

  // misc setup information
  Teuchos::ParameterList parameter_list;
  Teuchos::RCP<Mesh_maps_base> mesh_maps;

  // storage for the component concentration intermediate values
  Teuchos::RCP<Epetra_MultiVector> total_component_concentration_star;

  // process kernels
  Teuchos::RCP<Chemistry_PK> CPK;
  Teuchos::RCP<Transport_PK> TPK;
  Teuchos::RCP<Flow_PK> FPK; 

  Teuchos::ParameterList mpc_parameter_list;

  double T0, T1;
  int end_cycle;

  bool flow_enabled, transport_enabled, chemistry_enabled;


};


