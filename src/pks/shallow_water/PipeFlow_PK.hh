/*
  Pipe Flow PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Pipe flow model inherited from shallow water model.

  Author: Giacomo Capodaglio (gcapodaglio@lanl.gov)
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

  virtual double ComputePondedDepth(double WettedAngle) override;

  virtual double ComputeWettedAngle(double PondedDepth) override;

  virtual double ComputeWettedArea(double WettedAngle) override;

  virtual double ComputeWettedAngleNewton(double WettedArea) override;

  virtual std::vector<double> ComputeWettedQuantitiesEdge(int c, int e, double WettedAreaCell, double WettedAngleCell,
                                                          double Bc, double Bmax, const Epetra_MultiVector& B_n) override;

 private:
  static RegisteredPKFactory<PipeFlow_PK> reg_;

  double pipe_cross_section_;

  double Manning_coeff_;

};


}  // namespace ShallowWater
}  // namespace Amanzi

#endif
