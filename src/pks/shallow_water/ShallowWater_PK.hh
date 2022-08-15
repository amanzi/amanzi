/*
 Shallow Water PK

 Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
 Amanzi is released under the three-clause BSD License.
 The terms of use and "as is" disclaimer for this license are
 provided in the top-level COPYRIGHT file.

 Author: Svetlana Tokareva (tokareva@lanl.gov)
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

// inversion operation protected for small values
double
inverse_with_tolerance(double h, double tol);

class ShallowWater_PK : public PK_Physical, public PK_Explicit<TreeVector> {
 public:
  ShallowWater_PK(Teuchos::ParameterList& pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList>& glist,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& soln);

  ShallowWater_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
                  Teuchos::RCP<State> S, const std::string& pk_list_name,
                  std::vector<std::string>& component_names);

  ~ShallowWater_PK(){};

  virtual void Setup() override;
  virtual void Initialize() override;

  virtual double get_dt() override;
  virtual void set_dt(double dt) override{};

  // Advance PK by step size dt.
  virtual bool
  AdvanceStep(double t_old, double t_new, bool reinit = false) override;

  virtual void FunctionalTimeDerivative(double t, const TreeVector& A, TreeVector& f) override;

  virtual void ModifySolution(double t, TreeVector& A) override;

  // Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override;

  // Calculate any diagnostics prior to doing vis
  virtual void CalculateDiagnostics(const Tag& tag) override{};

  virtual std::string name() override { return "shallow water"; }

  // Bathymetry reconstruction on cell edge midpoints
  double BathymetryRectangularCellValue(int c, const AmanziGeometry::Point& xp,
                                        const Epetra_MultiVector& Bn);
  double BathymetryEdgeValue(int e, const Epetra_MultiVector& Bn);

  // Recalculate total depth for positivity of ponded depth
  double TotalDepthEdgeValue(int c, int e);

  // due to rotational invariance of SW equations, we need flux in the
  // x-direction only.
  std::vector<double> PhysicalFlux_x(const std::vector<double>&);

  std::vector<double>
  NumericalFlux_x(std::vector<double>&, std::vector<double>&);
  std::vector<double> NumericalFlux_x_Rusanov(const std::vector<double>&,
                                              const std::vector<double>&);
  std::vector<double> NumericalFlux_x_CentralUpwind(const std::vector<double>&,
                                                    const std::vector<double>&);

  std::vector<double> NumericalSource(const std::vector<double>&, int);

  // access
  double get_total_source() const { return total_source_; }

  // temporal discretization order
  int temporal_disc_order;
  // maximum cell area
  double cell_area_max_;
  
 private:
  void
  InitializeFieldFromField_(const std::string& field0,
                            const std::string& field1, bool call_evaluator);
  bool ErrorDiagnostics_(double t, int c, double h, double B, double ht);

 protected:
  Teuchos::RCP<Teuchos::ParameterList> glist_;
  Teuchos::RCP<Teuchos::ParameterList> sw_list_;
  Teuchos::RCP<TreeVector> soln_;
  Teuchos::RCP<State> S_;

  Key domain_;

  // numerical flux
  std::shared_ptr<NumericalFlux> numerical_flux_;

  // names of state fields
  Key velocity_key_, discharge_key_;
  Key ponded_depth_key_, prev_ponded_depth_key_;
  Key total_depth_key_, bathymetry_key_;
  Key hydrostatic_pressure_key_;
  Key riemann_flux_key_;

  std::string passwd_;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int dim_;
  
  // source terms
  std::vector<Teuchos::RCP<PK_DomainFunction>> srcs_;
  double total_source_;

 private:
  // boundary conditions
  std::vector<Teuchos::RCP<ShallowWaterBoundaryFunction>> bcs_;
  std::vector<Teuchos::RCP<Operators::BCs>> op_bcs_;

  // gravity magnitude
  double g_;

  // limited reconstruction
  bool use_limiter_;
  Teuchos::RCP<Operators::ReconstructionCellLinear> total_depth_grad_,
    bathymetry_grad_;
  Teuchos::RCP<Operators::ReconstructionCellLinear> discharge_x_grad_,
    discharge_y_grad_;
  Teuchos::RCP<Operators::LimiterCell> limiter_;

  // advanced cfl control
  double cfl_;
  int iters_, max_iters_;

 private:
  // factory registration
  static RegisteredPKFactory<ShallowWater_PK> reg_;
};

} // namespace ShallowWater
} // namespace Amanzi

#endif
