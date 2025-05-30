/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy
*/

/*!

Process kernel for coupling of Transport_PK and Chemistry_PK.
Reactive transport can be setup using a steady-state flow.
The two PKs are executed consequitively.
The input spec requires new keyword *reactive transport*.

.. code-block:: xml

  <ParameterList name="PK tree">  <!-- parent list -->
  <ParameterList name="_REACTIVE TRANSPORT">
    <Parameter name="PK type" type="string" value="reactive transport"/>
    <ParameterList name="_TRANSPORT">
      <Parameter name="PK type" type="string" value="transport"/>
    </ParameterList>
    <ParameterList name="_CHEMISTRY">
      <Parameter name="PK type" type="string" value="chemistry amanzi"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>

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

  // PK methods
  // -- dt is the minimum of the sub pks
  double get_dt() override;
  void set_dt(double dt) override;
  bool AdvanceStep(double t_old, double t_new, bool reinit) override;
  void CommitStep(double t_old, double t_new, const Tag& tag) override;

  std::string name() override { return "reactive transport"; }

 protected:
  bool chem_step_succeeded;
  double dTtran_, dTchem_;

 private:
  Teuchos::RCP<Transport::Transport_PK> transport_pk_;
  Teuchos::RCP<AmanziChemistry::Chemistry_PK> chemistry_pk_;

  // factory registration
  static RegisteredPKFactory<ReactiveTransport_PK> reg_;
};

} // namespace Amanzi
#endif
