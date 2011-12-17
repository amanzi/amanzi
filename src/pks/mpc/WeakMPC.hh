#ifndef _MPC_HPP_
#define _MPC_HPP_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Epetra_MpiComm.h"

#include "State.hh"
#include "Vector.hh"
#include "TreeVector.hh"

#include "PK.hh"
#include "PK_Factory.hh"

namespace Amanzi
{

  class WeakMPC : public Teuchos::VerboseObject<MPC>, public PK {

  public:
    WeakMPC(Teuchos::ParameterList &mpc_plist,
        Teuchos::RCP<State> &S);

    ~WeakMPC() {};

    // PK methods
    void initialize();

    double get_dT();

    bool advance(double dt, const Teuchos::RCP<State> &S0,
                 Teuchos::RCP<State> &S1, Teuchos::RCP<TreeVector> &solution);
    void compute_f(const double t, const Vector& u, const Vector& udot,
                   Vector& f);
    void commit_state(double dt, Teuchos::RCP<State> &S);

    void solution_to_state(const TreeVector& u, const TreeVector& udot);
  private:
    
    // PK container and factory
    PK_Factory pk_factory_;
    std::vector< Teuchos::RCP<PK> > sub_pks_

    // states
    Teuchos::RCP<State> S_;

    // misc setup information
    Teuchos::ParameterList mpc_plist_;
  };

} // close namespace Amanzi

#endif
