/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Multi process coupler for sequential coupling across subdomains.

/*!

Noniterative sequential coupling simply calls each PK's AdvanceStep() method in
order.  Each child PK is a subdomain of the parent PK.

.. _mpc-weak-domain-decomposition-spec:
.. admonition:: mpc-weak-domain-decomposition-spec

    * `"domain`" ``[string]`` Name of the domain that is decomposed, the parent
      mesh.

    INCLUDES:

    - ``[weak-mpc-spec]`` *Is a* WeakMPC_.

*/

#pragma once

#include "mpc.hh"

namespace Amanzi {

class MPCWeakDomainDecomposition : public MPC<PK> {

public:
  MPCWeakDomainDecomposition(Teuchos::ParameterList& pk_tree,
                             const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                             const Teuchos::RCP<State>& S,
                             const Teuchos::RCP<TreeVector>& solution)
    : MPC<PK>(pk_tree, global_list, S, solution),
      PK(pk_tree, global_list, S, solution)
    {
      MPCWeakDomainDecomposition::init_(S);
    }

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt() override;
  virtual void set_dt(double dt) override;

  // -- advance each sub pk dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit) override;

protected:
  void init_(const Teuchos::RCP<State>& S);

protected:
  Comm_ptr_type comm_;

private:
  // factory registration
  static RegisteredPKFactory<MPCWeakDomainDecomposition> reg_;

};

} // namespace Amanzi
