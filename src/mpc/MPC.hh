#ifndef _MPC_HPP_
#define _MPC_HPP_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Epetra_MpiComm.h"
#include "State.hh"
#include "PK.hh"
#include "PK_Factory.hh"
#include "ObservationData.H"
#include "Unstructured_observations.hpp"
#include "Vis.hh"

namespace Amanzi
{

  class MPC : public Teuchos::VerboseObject<MPC> {

  public:
    // types
    typedef std::vector<Teuchos::RCP<PK> > PKs;

    MPC(Teuchos::ParameterList parameter_list_,
        Teuchos::RCP<AmanziMesh::Mesh> mesh_maps_,
        Epetra_MpiComm* comm_,
        Amanzi::ObservationData& output_observations_);

    ~MPC() {};

    void initialize_state(Teuchos::RCP<State>);
    double get_dT();
    bool advance_transient(double, Teuchos::RCP<State>, Teuchos::RCP<State>);
    void commit_state(double, Teuchos::RCP<State>);
    void cycle_driver ();

  private:
    void mpc_init();
    void read_parameter_list();

    // PK container and factory
    PK_Factory pk_factory;
    PKs pks;

    // states
    Teuchos::RCP<State> S;

    // misc setup information
    Teuchos::ParameterList parameter_list;
    Teuchos::RCP<AmanziMesh::Mesh> mesh_maps;

    Teuchos::ParameterList mpc_parameter_list;

    double T0, T1;
    int end_cycle;

    // Epetra communicator
    Epetra_MpiComm* comm;

    // observations
    Amanzi::ObservationData& output_observations;
    Amanzi::Unstructured_observations* observations;

    // visualization
    Amanzi::Vis *visualization;

  };

} // close namespace Amanzi

#endif
