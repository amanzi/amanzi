/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  MPC PK

  Sequential coupling of shallow water and solute transport.
*/

#ifndef AMANZI_SHALLOW_WATER_TRANSPORT_PK_HH_
#define AMANZI_SHALLOW_WATER_TRANSPORT_PK_HH_

#include "PK_MPCWeak.hh"

namespace Amanzi {

class ShallowWaterTransport_PK : public PK_MPCWeak {
 public:
  ShallowWaterTransport_PK(Teuchos::ParameterList& pk_tree,
                           const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& soln);

  // PK methods
  virtual double get_dt() override;
  virtual void Setup() override;
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit) override;

 private:
  double cfl_;
  int failed_steps_;

 private:
  // factory registration
  static RegisteredPKFactory<ShallowWaterTransport_PK> reg_;
};

} // namespace Amanzi

#endif
