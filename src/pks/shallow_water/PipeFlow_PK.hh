/*
  Pipe Flow PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Pipe flow model inherited from shallow water model. 
*/

#ifndef AMANZI_PIPE_FLOW_PK_HH_
#define AMANZI_PIPE_FLOW_PK_HH_

#include "ShallowWater_PK.hh"

namespace Amanzi {
namespace ShallowWater {

class PipeFlow_PK : public ShallowWater_PK {

 public:
  PipeFlow_PK(Teuchos::ParameterList& pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList>& glist,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& soln);

  ~PipeFlow_PK() {};

  virtual double NumericalSourceFriction(double h, double qx, double WettedAngle) override;

  virtual void Initialize() override;

  virtual void UpdateWettedAngle() override;

 private:
  static RegisteredPKFactory<PipeFlow_PK> reg_;

  double pipe_diameter_;

  double Manning_coeff_;

  double celerity_;

};


}  // namespace ShallowWater
}  // namespace Amanzi

#endif
