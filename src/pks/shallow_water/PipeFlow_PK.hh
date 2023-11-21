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

  virtual void Setup() override;

  virtual double NumericalSourceFriction(double h, double qx, double qy, double WettedAngle, int component) override;

  virtual std::vector<double> NumericalSourceBedSlope(int c, double htc, double Bc,
                                                      double Bmax, const Epetra_MultiVector& B_n,
                                                      std::vector<int> bc_model, std::vector<double> bc_value_h) override;

  virtual std::vector<double> NumericalSourceBedSlope(int c, double htc, double Bc,
                                              double Bmax, const Epetra_MultiVector& B_n) override;

  virtual void UpdateSecondaryFields() override;

  virtual double ComputeTotalDepth(double PrimaryVar, double Bathymetry, double WettedAngle) override;

  virtual double ComputeWaterDepth(double WettedAngle) override;

  virtual double ComputeWettedAngle(double WaterDepth) override;

  virtual double ComputeWettedArea(double WettedAngle) override;

  virtual double ComputeWettedAngleNewton(double WettedArea) override;

  virtual double ComputeHydrostaticPressureForce(std::vector<double> SolArray) override;

  virtual double ComputePressureHead(double WettedAngle) override;

  virtual void SetupPrimaryVariableKeys() override;

  virtual void SetupExtraEvaluatorsKeys() override;

  virtual void ScatterMasterToGhostedExtraEvaluators() override;

  virtual void UpdateExtraEvaluators() override;

  virtual void SetPrimaryVariableBC(Teuchos::RCP<Teuchos::ParameterList> &bc_list) override;

  virtual void InitializeFields() override;

  virtual void ComputeExternalForcingOnCells(std::vector<double> &forcing) override;

  void GetDx(const int & cell, double & dx);

  virtual void ComputeCellArrays() override;

  virtual void ProjectNormalOntoMeshDirection(int c, AmanziGeometry::Point &normal) override;

  bool IsJunction(const int & cell);

  virtual std::vector<double> ComputeFieldsOnEdge(int c, int e, double htc, double Bc, double Bmax, const Epetra_MultiVector& B_n) override;

 private:
  static RegisteredPKFactory<PipeFlow_PK> reg_;

  double pipe_cross_section_;

  double Manning_coeff_;

 protected:
 Key water_depth_key_, pressure_head_key_;
 // unit vector that defines the pipe direction
 // (both components are zero for junction cell)
 Key direction_key_;

};


}  // namespace ShallowWater
}  // namespace Amanzi

#endif
