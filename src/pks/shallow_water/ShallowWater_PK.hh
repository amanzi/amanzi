/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.


 Authors: Svetlana Tokareva (tokareva@lanl.gov)
          Giacomo Capodaglio (gcapodaglio@lanl.gov)
          Naren Vohra (vohra@lanl.gov)

The mathematical model describing two-dimensional shallow water flow is

.. math::
  \begin{align*}
  & h_t + (hu)_x + (hv)_y = 0, \\
  & (hu)_t + (hu^2 + \frac{1}{2} gh^2)_x + (huv)_y = -ghB_x \\
  & (hv)_t + (huv)_x + (hv^2 + \frac{1}{2} gh^2)_y = -ghB_y
  \end{align*}

Here
:math:`h` [m] is water depth,
:math:`g` [m/s^2] is gravity acceleration,
:math:`u` [m/s] is depth averaged velocity in x direction,
:math:`v` [m/s] is depth averaged velocity in y direction,
:math:`B` [m] is bottom elevation (bathymetry),
:math:`H = h + B` [m] is water surface elevation.


Global parameters
.................
Global parameters are placed in the sublist `"shallow water`".
The list of global parameters include:

.. admonition:: shallow_water_params-spec

  * `"domain name`" ``[string]`` specifies mesh name that defined domain of this PK.
    Default is `"domain`".

  * `"cfl`" ``[double]`` is a safety factor (less than 1) applied to a stable
    timestep estimate. Default value is 1.

  * `"use limiter`" ``[bool]`` turns on/off limiters on all linear constructions.
    Default value is *false*.

*/

#ifndef AMANZI_SHALLOW_WATER_PK_HH_
#define AMANZI_SHALLOW_WATER_PK_HH_

#include <memory>

// TPLs
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "BCs.hh"
#include "CompositeVector.hh"
#include "DenseMatrix.hh"
#include "DenseVector.hh"
#include "Explicit_TI_RK.hh"
#include "Key.hh"
#include "LimiterCell.hh"
#include "NumericalFlux.hh"
#include "PK.hh"
#include "PK_Utils.hh"
#include "PK_Explicit.hh"
#include "PK_Factory.hh"
#include "PK_Physical.hh"
#include "PK_DomainFunction.hh"
#include "ReconstructionCellLinear.hh"
#include "State.hh"
#include "Tensor.hh"
#include "Units.hh"
#include "VerboseObject.hh"
#include "WhetStoneMeshUtils.hh"

#include "ShallowWaterBoundaryFunction.hh"

namespace Amanzi {
namespace ShallowWater {

class ShallowWater_PK : public PK_Physical, public PK_Explicit<TreeVector> {
 public:
  ShallowWater_PK(Teuchos::ParameterList& pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList>& glist,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& soln);
  ~ShallowWater_PK() {};

  virtual void parseParameterList() override {};
  virtual void Setup() override;
  virtual void Initialize() override;

  virtual double get_dt() override;
  virtual void set_dt(double dt) override {};

  // Advance PK by step size dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false) override;

  virtual void FunctionalTimeDerivative(double t, const TreeVector& A, TreeVector& f) override;

  virtual void ModifySolution(double t, TreeVector& A) override { VerifySolution_(A); }

  virtual void ScatterMasterToGhostedExtraEvaluators() {};

  virtual void UpdateExtraEvaluators() {};

  virtual void SetPrimaryVariableBC(Teuchos::RCP<Teuchos::ParameterList>& bc_list);

  virtual void InitializeFields();

  virtual void ComputeExternalForcingOnCells(std::vector<double>& forcing);

  // Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override;

  // Calculate any diagnostics prior to doing vis
  virtual void CalculateDiagnostics(const Tag& tag) override {};

  virtual std::string name() override { return "shallow water"; }

  // Bathymetry reconstruction on cell edge midpoints
  double BathymetryEdgeValue(int e, const Epetra_MultiVector& Bn);

  // Recalculate total depth for positivity of ponded depth
  double TotalDepthEdgeValue(int c,
                             int e,
                             double htc,
                             double Bc,
                             double Bmax,
                             const Epetra_MultiVector& B_n);

  double ComputeFieldOnEdge(int c,
                            int e,
                            double htc,
                            double Bc,
                            double Bmax,
                            const Epetra_MultiVector& B_n);


  std::vector<double>
  NumericalSourceBedSlope(int c, double htc, double Bc, double Bmax, const Epetra_MultiVector& B_n);

  virtual void UpdateSecondaryFields();

  virtual double ComputeHydrostaticPressureForce(std::vector<double> Data)
  {
    return g_ * 0.5 * Data[0] * Data[0];
  };

  double inverse_with_tolerance(double h, double tol);

  // access
  double get_total_source() const { return total_source_; }

 protected:
  void InitializeFieldFromField_(const std::string& field0,
                                 const std::string& field1,
                                 bool call_evaluator);

  void VerifySolution_(TreeVector& A);
  int ErrorDiagnostics_(double t, int c, double h);

 protected:
  Teuchos::RCP<Teuchos::ParameterList> sw_list_;
  Teuchos::RCP<TreeVector> soln_;

  Key domain_;

  // numerical flux
  std::shared_ptr<NumericalFlux> numerical_flux_;

  // names of state fields
  Key velocity_key_, discharge_key_;
  // the primary variable is:
  // ponded depth for shallow water
  // wetted area for pipe flow
  Key primary_variable_key_, prev_primary_variable_key_;
  Key water_depth_key_, total_depth_key_, bathymetry_key_;
  Key hydrostatic_pressure_key_;
  Key riemann_flux_key_;
  Key source_key_;

  std::string passwd_;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int dim_;
  int ncells_owned_;

  // source terms
  std::vector<Teuchos::RCP<PK_DomainFunction>> srcs_;
  double total_source_;

  // gravity magnitude
  double g_;

  double velocity_desingularization_eps_;

  // boundary conditions
  std::vector<Teuchos::RCP<ShallowWaterBoundaryFunction>> bcs_;
  std::vector<Teuchos::RCP<Operators::BCs>> op_bcs_;

  // limited reconstruction
  bool use_limiter_;
  Teuchos::RCP<Operators::ReconstructionCellLinear> total_depth_grad_;
  Teuchos::RCP<Operators::ReconstructionCellLinear> discharge_x_grad_, discharge_y_grad_;
  Teuchos::RCP<Operators::LimiterCell> limiter_;

  // advanced cfl control
  double cfl_, cfl_positivity_;
  int iters_, max_iters_;

  // control of PK
  int temporal_disc_order_;
  double cell_area2_max_;

 private:
  // factory registration
  static RegisteredPKFactory<ShallowWater_PK> reg_;
};

} // namespace ShallowWater
} // namespace Amanzi

#endif
