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
#include "Key.hh"
#include "LimiterCell.hh"
#include "NumericalFlux.hh"
#include "PK.hh"
#include "PK_Explicit.hh"
#include "PK_Factory.hh"
#include "PK_Physical.hh"
#include "PK_DomainFunction.hh"
#include "ReconstructionCell.hh"
#include "State.hh"
#include "Tensor.hh"
#include "Units.hh"
#include "VerboseObject.hh"
#include "WhetStoneMeshUtils.hh"

#include "ShallowWaterBoundaryFunction.hh"

namespace Amanzi {
namespace ShallowWater {
    
class ShallowWater_PK : public PK_Physical,
                        public PK_Explicit<Epetra_Vector> {
 public:
  ShallowWater_PK(Teuchos::ParameterList& pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList>& glist,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& soln);

  ShallowWater_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
                  Teuchos::RCP<State> S,
                  const std::string& pk_list_name,
                  std::vector<std::string>& component_names);

  ~ShallowWater_PK() {};

  virtual void Setup(const Teuchos::Ptr<State>& S) override;
  virtual void Initialize(const Teuchos::Ptr<State>& S) override;

  virtual double get_dt() override;
  virtual void set_dt(double dt) override {};

  // Advance PK by step size dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit=false) override;

  virtual void FunctionalTimeDerivative(double t, const Epetra_Vector& component,
                                        Epetra_Vector& f_component) override;

  // Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) override;

  // Calculate any diagnostics prior to doing vis
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) override {};

  virtual std::string name() override { return "Shallow water PK"; }
                            
  // Bathymetry reconstruction on cell edge midpoints
  double BathymetryRectangularCellValue(int c, const AmanziGeometry::Point& xp, const Epetra_MultiVector& Bn);
  double BathymetryEdgeValue(int e, const Epetra_MultiVector& Bn);

  // due to rotational invariance of SW equations, we need flux in the x-direction only.
  std::vector<double> PhysicalFlux_x(const std::vector<double>&);

  std::vector<double> NumericalFlux_x(std::vector<double>&, std::vector<double>&);
  std::vector<double> NumericalFlux_x_Rusanov(const std::vector<double>&, const std::vector<double>&);
  std::vector<double> NumericalFlux_x_CentralUpwind(const std::vector<double>&, const std::vector<double>&);

  std::vector<double> PhysicalSource(const std::vector<double>&);
  std::vector<double> NumericalSource(const std::vector<double>&, int);
  std::vector<double> PhysFlux_x(std::vector<double> U);
  std::vector<double> PhysFlux_y(std::vector<double> U);
    
  std::vector<double> ComputePhiTotal(int K, std::vector<std::vector<double> >& U);
  std::vector<double> ResidualsLF(int K, int j, std::vector<std::vector<double> > U);
  std::vector<double> ResidualsTimeSpace(int c, int i, std::vector<std::vector<double> > U, std::vector<std::vector<double> > U_pr, double dt);
                                                                          
  std::vector<double> EvalSol_vol(std::vector<std::vector<double>> U, int qp, int c);
  std::vector<double> EvalSol_x_vol(std::vector<std::vector<double>> U, int qp, int c);
  std::vector<double> EvalSol_y_vol(std::vector<std::vector<double>> U, int qp, int c);
  std::vector<double> EvalSol_face(std::vector<std::vector<double>> U, int qpf, int f);
  std::vector<double> EvalPhySource_vol(std::vector<std::vector<double>> U, int qp, int c);
  double basis_value(int i, int c, AmanziGeometry::Point x);
  double basis_value_quad(int i, int c, AmanziGeometry::Point x);
  std::vector<double> basis_grad(int i, int c, AmanziGeometry::Point x);
  std::vector<double> basis_grad_quad(int i, int c, AmanziGeometry::Point x);
  std::vector<double> get_barycentric(std::vector<AmanziGeometry::Point> vertices, AmanziGeometry::Point x);

  // access
  double get_total_source() const { return total_source_; }

 private:
  bool ErrorDiagnostics_(int c, double h, double B, double ht);

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
  Key ponded_depth_key_;
  Key total_depth_key_;
  Key bathymetry_key_;
  Key hydrostatic_pressure_key_;

  std::string passwd_;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int dim_;
     
  // P1 basis, basis gradient
  std::vector<std::vector<std::vector<double>>> phi_; // P1 basis evaluated at volume quadrature points
  std::vector<std::vector<std::vector<double>>> phi_x_; // x partial derivative
  std::vector<std::vector<std::vector<double>>> phi_y_; // y partial derivative
  std::vector<std::vector<std::vector<double>>> phi_face_; // P1 basis evaluated at face quadrature points
                          
  // Volume, face quadrature weights
  std::vector<std::vector<double>> weights_vol_;
  std::vector<std::vector<double>> weights_face_;
                          
  // source terms
  std::vector<Teuchos::RCP<PK_DomainFunction> > srcs_;
  double total_source_;

 private:
  // boundary conditions
  std::vector<Teuchos::RCP<ShallowWaterBoundaryFunction> > bcs_;
  std::vector<Teuchos::RCP<Operators::BCs> > op_bcs_;

  // gravity magnitude
  double g_;

  // limited reconstruction
  Teuchos::RCP<Operators::ReconstructionCell> total_depth_grad_, bathymetry_grad_;
  Teuchos::RCP<Operators::ReconstructionCell> discharge_x_grad_, discharge_y_grad_;
  Teuchos::RCP<Operators::LimiterCell> limiter_;

  // advanced cfl control
  double cfl_;
  int iters_, max_iters_;

 private:
  // factory registration
  static RegisteredPKFactory<ShallowWater_PK> reg_;
};
    
}  // namespace ShallowWater
}  // namespace Amanzi

#endif
