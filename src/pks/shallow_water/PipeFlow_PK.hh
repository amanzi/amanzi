/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Giacomo Capodaglio (gcapodaglio@lanl.gov)
*/

/*
  Pipe Flow PK

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

  ~PipeFlow_PK(){};

  virtual void Setup() override;

  virtual double get_dt() override;

  virtual void FunctionalTimeDerivative(double t, const TreeVector& A, TreeVector& f) override;

  double NumericalSourceFriction(double h,
                                 double qx,
                                 double qy,
                                 double WettedAngle,
                                 int component,
                                 double PipeD);

  std::vector<double> NumericalSourceBedSlope(int c,
                                              double htc,
                                              double Bc,
                                              double Bmax,
                                              const Epetra_MultiVector& B_n,
                                              double PipeD,
                                              std::vector<int> bc_model,
                                              std::vector<double> bc_value_h);

  std::vector<double> NumericalSourceBedSlope(int c,
                                              double htc,
                                              double Bc,
                                              double Bmax,
                                              const Epetra_MultiVector& B_n,
                                              double PipeD);

  virtual void UpdateSecondaryFields() override;

  double ComputeTotalDepth(double PrimaryVar, double Bathymetry, double WettedAngle, double PipeD);

  double ComputeWaterDepth(double WettedAngle, double PipeD);

  double ComputeWettedAngle(double WaterDepth, double PipeD);

  double ComputeWettedArea(double WettedAngle, double PipeD);

  double ComputeWettedAngleNewton(double WettedArea, double PipeD);

  virtual double ComputeHydrostaticPressureForce(std::vector<double> Data) override;

  virtual void ScatterMasterToGhostedExtraEvaluators() override;

  virtual void UpdateExtraEvaluators() override;

  virtual void SetPrimaryVariableBC(Teuchos::RCP<Teuchos::ParameterList>& bc_list) override;

  virtual void InitializeFields() override;

  virtual void ComputeExternalForcingOnCells(std::vector<double>& forcing) override;

  void SkipFace(AmanziGeometry::Point normal, bool& skipFace);

  void GetDx(const int& cell, double& dx);

  void ProjectNormalOntoMeshDirection(int c, AmanziGeometry::Point& normal);

  std::vector<double> ComputeFieldsOnEdge(int c,
                                          int e,
                                          double htc,
                                          double Bc,
                                          double Bmax,
                                          const Epetra_MultiVector& B_n,
                                          double PipeD);

 private:
  void ComputeCellArrays_();
  double ComputePressureHead_(double WettedAngle, double PipeD);
  bool IsJunction_(int c);

 private:
  static RegisteredPKFactory<PipeFlow_PK> reg_;

  double Manning_coeff_;
  double celerity_;

 protected:
  Key pressure_head_key_;
  // unit vector that defines the pipe direction
  // (both components are zero for junction cell)
  Key direction_key_;
  Key wetted_angle_key_;
  Key diameter_key_;

  std::vector<int> model_cells_owned_;
  std::vector<int> junction_cells_owned_;
  std::vector<int> model_cells_wghost_;
  std::vector<int> junction_cells_wghost_;
  bool cellArraysInitDone_ = false;
};


} // namespace ShallowWater
} // namespace Amanzi

#endif
