#ifndef _MPC_HPP_
#define _MPC_HPP_


#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Epetra_MpiComm.h"
#include "State.hh"
#include "Chemistry_State.hh"
#include "chemistry_pk_base.hh"
#include "Transport_PK.hh"
#include "Flow_PK.hh"
#include "ObservationData.hh"
#include "Unstructured_observations.hh"
#include "visualization.hh"
#include "checkpoint.hh"
#include "chemistry_data.hh"

namespace Amanzi {

  enum time_integration_mode { STEADY, TRANSIENT, INIT_TO_STEADY, TRANSIENT_STATIC_FLOW };

class MPC : public Teuchos::VerboseObject<MPC> {
 public:
  MPC(Teuchos::ParameterList parameter_list_,
      Teuchos::RCP<AmanziMesh::Mesh> mesh_maps_,
      Epetra_MpiComm* comm_,
      Amanzi::ObservationData& output_observations_); 
  ~MPC () {};
    
  void cycle_driver ();

  // special function for walkabout
  void populate_walkabout_data();

 private:
  void mpc_init();
  void read_parameter_list();
  double time_step_limiter (double T, double dT, double T_end);

  // states
  Teuchos::RCP<State> S;
  Teuchos::RCP<AmanziChemistry::Chemistry_State> CS;
    
  // misc setup information
  Teuchos::ParameterList parameter_list;
  Teuchos::RCP<AmanziMesh::Mesh> mesh_maps;
    
  // storage for the component concentration intermediate values
  Teuchos::RCP<Epetra_MultiVector> total_component_concentration_star;
    
  // process kernels
  Teuchos::RCP<AmanziChemistry::Chemistry_PK_Base> CPK;
  Teuchos::RCP<AmanziTransport::Transport_PK> TPK;
  Teuchos::RCP<AmanziFlow::Flow_PK> FPK; 
    
  Teuchos::ParameterList mpc_parameter_list;
    
  double T0, T1, dTsteady, dTtransient, Tswitch;
  int end_cycle;
  time_integration_mode ti_mode;
  
  bool steady, init_to_steady, static_steady;
  bool flow_enabled, transport_enabled, chemistry_enabled;

  int transport_subcycling; 
  double chem_trans_dt_ratio;

  double ti_rescue_factor_;

  std::string flow_model;
    
  // restart from checkpoint file
  bool restart_requested;
  std::string restart_from_filename;

  // Epetra communicator
  Epetra_MpiComm* comm;
    
  // observations
  Amanzi::ObservationData&  output_observations;
  Amanzi::Unstructured_observations* observations;
    
  // visualization
  Teuchos::Ptr<Amanzi::Visualization> visualization;
  std::vector<std::string> auxnames;
    
  // checkpoint/restart 
  Teuchos::Ptr<Amanzi::Checkpoint> restart;
 
  // walkabout
  Teuchos::Ptr<Amanzi::Checkpoint> walkabout;  

  // time period control
  std::vector<std::pair<double,double> > reset_info_;

  // picard flag
  bool do_picard_;

  // stor for chemistry data to allow repeat of chemistry step
  Teuchos::RCP<chemistry_data> chem_data_;
};
    
} // namespace Amanzi

#endif
