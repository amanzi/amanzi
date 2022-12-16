/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy
*/

/*
  This is the mpc_pk component of the Amanzi code.

  Process kernel for coupling of Transport_PK and Chemistry_PK.
*/


#ifndef AMANZI_REACTIVETRANSPORT_PK_HH_
#define AMANZI_REACTIVETRANSPORT_PK_HH_

#include "Teuchos_RCP.hpp"

#include "PK.hh"
#include "Transport_PK.hh"
#include "Chemistry_PK.hh"
#include "PK_Factory.hh"
#include "PK_MPCAdditive.hh"

namespace Amanzi {

class ReactiveTransport_PK : public PK_MPCAdditive<PK> {
 public:
  ReactiveTransport_PK(Teuchos::ParameterList& pk_tree,
                       const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                       const Teuchos::RCP<State>& S,
                       const Teuchos::RCP<TreeVector>& soln);

  ~ReactiveTransport_PK(){};

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt();
  virtual void set_dt(double dt);

  // -- advance each sub pk dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false);

  virtual void Initialize();

  virtual void CommitStep(double t_old, double t_new, const Tag& tag);

  std::string name() { return "reactive transport"; }

 protected:
  bool chem_step_succeeded;
  bool storage_created;
  double dTtran_, dTchem_;

 private:
  // storage for the component concentration intermediate values
  Teuchos::RCP<Epetra_MultiVector> total_component_concentration_stor;

  Teuchos::RCP<Transport::Transport_PK> transport_pk_;
  Teuchos::RCP<AmanziChemistry::Chemistry_PK> chemistry_pk_;
  // int master_, slave_;


  // factory registration
  static RegisteredPKFactory<ReactiveTransport_PK> reg_;
};

} // namespace Amanzi
#endif
