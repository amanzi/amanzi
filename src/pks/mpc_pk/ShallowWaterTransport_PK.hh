/*
  MPC PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov

  Sequential coupling of shallow water and solute transport.
*/

#ifndef AMANZI_SHALLOW_WATER_TRANSPORT_PK_HH_
#define AMANZI_SHALLOW_WATER_TRANSPORT_PK_HH_

#include "PK_MPCWeak.hh"

namespace Amanzi {

class ShallowWaterTransport_PK : public PK_MPCWeak {
 public:
  ShallowWaterTransport_PK(
      Teuchos::ParameterList& pk_tree,
      const Teuchos::RCP<Teuchos::ParameterList>& global_list,
      const Teuchos::RCP<State>& S,
      const Teuchos::RCP<TreeVector>& soln)
    : PK_MPCWeak(pk_tree, global_list, S, soln) {};

  // PK methods
  // -- setup coupling information
  virtual void Setup(const Teuchos::Ptr<State>& S);

 private:
  // factory registration
  static RegisteredPKFactory<ShallowWaterTransport_PK> reg_;
};


// -----------------------------------------------------------------------------
// Setup of PK
// -----------------------------------------------------------------------------
void ShallowWaterTransport_PK::Setup(const Teuchos::Ptr<State>& S)
{
  PK_MPCWeak::Setup(S);

  // tell transport to use Riemann velocity
  std::string domain;
  Teuchos::ParameterList plist;
  plist.sublist("physical models and assumptions")
       .set<std::string>("darcy flux key", Keys::getKey(domain, "riemann_flux"));
  std::cout << name_ << std::endl;
}

}  // namespace Amanzi

#endif
