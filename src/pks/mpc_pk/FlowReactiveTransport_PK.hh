/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy
*/

/*!

PK for coupling of Flow PK with Transport_PK and Chemistry_PK
Amanzi uses operator splitting approach for coupled physical kernels.
The coupling of PKs is described as a tree where flow and reactive
transport are executed consequitively.
The input spec requires new keyword *flow reactive transport*.

.. code-block:: xml

  <ParameterList name="PK tree">  <!-- parent list -->
  <ParameterList name="_FLOW and REACTIVE TRANSPORT">
    <Parameter name="PK type" type="string" value="flow reactive transport"/>
    <ParameterList name="_FLOW">
      <Parameter name="PK type" type="string" value="darcy"/>
    </ParameterList>
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
  </ParameterList>

This example describe four PKs identified by keywords *darcy*, *reactive transport*,
*transport*, and *chemistry amanzi*.
The flow is fully saturated.
The transport of reactive chemicals is based on the native chemistry package *chemistry amanzi*.

Details of PKs are organized as a plain list of ParameterLists.
Note that *reactive transport* is MPC-PK and hence its description is short.

.. code-block:: xml

  <ParameterList name="PKs">
  <ParameterList name="_FLOW and REACTIVE TRANSPORT">
    <Parameter name="PK type" type="string" value="flow reactive transport"/>
    <Parameter name="PKs order" type="Array(string)" value="{_FLOW, _REACTIVE TRANSPORT}"/>
    <Parameter name="master PK index" type="int" value="0"/>
  </ParameterList>

  <ParameterList name="_REACTIVE TRANSPORT">
    <Parameter name="PK type" type="string" value="reactive transport"/>
    <Parameter name="PKs order" type="Array(string)" value="{_CHEMISTRY, _TRANSPORT}"/>
  </ParameterList>

  <ParameterList name="_FLOW">
    ...
  </ParameterList>

  <ParameterList name="_TRANSPORT">
    ...
  </ParameterList>

  <ParameterList name="_CHEMISTRY">
    ...
  </ParameterList>
  </ParameterList>

*/

#ifndef AMANZI_FLOWREACTIVETRANSPORT_PK_HH_
#define AMANZI_FLOWREACTIVETRANSPORT_PK_HH_

#include "Teuchos_RCP.hpp"

#include "PK.hh"
#include "PK_Factory.hh"
#include "PK_MPCSubcycled.hh"

namespace Amanzi {

class FlowReactiveTransport_PK : public PK_MPCSubcycled {
 public:
  FlowReactiveTransport_PK(Teuchos::ParameterList& pk_tree,
                           const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& soln);

  ~FlowReactiveTransport_PK(){};

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt();
  virtual void set_dt(double dt);

  // -- advance each sub pk from t_old to t_new.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false);

  virtual void CommitStep(double t_old, double t_new, const Tag& tag);

  std::string name() { return "flow reactive transport"; }

 private:
  // factory registration
  static RegisteredPKFactory<FlowReactiveTransport_PK> reg_;
};

} // namespace Amanzi
#endif
