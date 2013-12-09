/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_FLOW_MATRIX_BASE_HH_
#define AMANZI_FLOW_MATRIX_BASE_HH_

#include <strings.h>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "boundary_function.hh"

#include "FlowTypeDefs.hh"
#include "RelativePermeability.hh"


namespace Amanzi {
namespace AmanziFlow {

/*
  This class is used in the factory of linear solvers and requires that 
  all exposed members be defined. At the moment, I do not have a better design.
*/
template<class Vector, class VectorSpace>
class Matrix {
 public:
  Matrix() {};
  Matrix(Teuchos::RCP<State> S, Teuchos::RCP<RelativePermeability> rel_perm)
      : S_(S), rel_perm_(rel_perm) { mesh_ = S_->GetMesh(); }
  ~Matrix() {};

  // main methods
  virtual void Init() {};
  virtual void CreateMassMatrices(int method, std::vector<WhetStone::Tensor>& K) {};
  virtual void CreateStiffnessMatricesDarcy(int method, std::vector<WhetStone::Tensor>& K) {};
  virtual void CreateStiffnessMatricesRichards() {};

  virtual void SymbolicAssemble() {};
  virtual void Assemble() {};
  virtual void AssembleSchurComplement(std::vector<int>& bc_model, std::vector<bc_tuple>& bc_values) {};

  virtual void AddGravityFluxesDarcy(double rho, const AmanziGeometry::Point& gravity,
                                     std::vector<WhetStone::Tensor>& K) {};
  virtual void AddGravityFluxesRichards(double rho, const AmanziGeometry::Point& gravity,
                                        std::vector<WhetStone::Tensor>& K) {};

  virtual void AddTimeDerivative(
      const Epetra_MultiVector& p, const Epetra_MultiVector& phi, double rho, double dT) {};
  virtual void AddTimeDerivativeSpecificStorage(
      const Epetra_MultiVector& p, const Epetra_MultiVector& ss, double g, double dT) {};
  virtual void AddTimeDerivativeSpecificYield(
      const Epetra_MultiVector& p, const Epetra_MultiVector& sy, double g, double dT) {};

  virtual void ApplyBoundaryConditions(std::vector<int>& bc_model, std::vector<bc_tuple>& bc_values) {};

  virtual int Apply(const Vector& v, Vector& av) const = 0;
  virtual int ApplyInverse(const Vector& v, Vector& hv) const = 0;

  virtual const VectorSpace& DomainMap() const = 0;
  virtual const VectorSpace& RangeMap() const = 0;

  virtual void InitPreconditioner(const std::string& prec_name, const Teuchos::ParameterList& prec_list) {};
  virtual void UpdatePreconditioner() {};
  virtual void DestroyPreconditioner() {};

  virtual void DeriveMassFlux(const Vector& p, Vector& flux,
                              std::vector<int>& bc_model, std::vector<bc_tuple>& bc_values) {};

  virtual void CreateRHSVectors() {};

  // common members
  void AddActionProperty(int action) { actions_ |= action; }
  void SetSymmetryProperty(bool flag_symmetry) { flag_symmetry_ = flag_symmetry; }

  double ComputeResidual(const Vector& v, Vector& r);
  double ComputeNegativeResidual(const Vector& v, Vector& r);

  Teuchos::RCP<Vector> rhs() { return rhs_; }

  // members collecting statistics
  int nokay() { return nokay_; }
  int npassed() { return npassed_; }

 protected:
  Teuchos::RCP<State> S_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<RelativePermeability> rel_perm_;

  Teuchos::RCP<Vector> rhs_;

  int actions_;  // aplly, apply inverse, or both
  bool flag_symmetry_;
  int nokay_, npassed_;  // performance of algorithms generating matrices 
};


typedef Matrix<CompositeVector, CompositeVectorSpace> FlowMatrix;


/* ******************************************************************
* Linear algebra operations with matrices: r = f - A * u
****************************************************************** */
template<class Vector, class VectorSpace>
double Matrix<Vector, VectorSpace>::ComputeResidual(const Vector& u, Vector& r)
{
  Apply(u, r);
  r.Update(1.0, *rhs_, -1.0);

  double rnorm;
  r.Norm2(&rnorm);
  return rnorm;
}


/* ******************************************************************
* Linear algebra operations with matrices: r = A * u - f                                                 
****************************************************************** */
template<class Vector, class VectorSpace>
double Matrix<Vector, VectorSpace>::ComputeNegativeResidual(const Vector& u, Vector& r)
{
  Apply(u, r);
  r.Update(-1.0, *rhs_, 1.0);

  double rnorm;
  r.Norm2(&rnorm);
  return rnorm;
}


}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
