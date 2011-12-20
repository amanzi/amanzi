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
#include "Vis.hpp"

namespace Amanzi
{

  class WeakMPC : public Teuchos::VerboseObject<MPC>, public PK {

  public:
    WeakMPC(Teuchos::ParameterList &parameter_list,
        Teuchos::RCP<State> &S);

    ~WeakMPC() {};

    // PK methods
    void initialize(Teuchos::RCP<State>& S);
    double get_dT();
    bool advance(double dt, const Teuchos::RCP<State> &S0,
                 Teuchos::RCP<State> &S1, Teuchos::RCP<Vector> &solution);
    void compute_f(const double t, const Vector& u, const Vector& udot,
                   Vector& f);
    void commit_state(double dt, Teuchos::RCP<State> &S);

  private:
    void mpc_init();
    void read_parameter_list();

    // PK container and factory
    PK_Factory pk_factory_;
    void cycle_driver ();
    std::vector< Teuchos::RCP<PK> > sub_pks_

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
