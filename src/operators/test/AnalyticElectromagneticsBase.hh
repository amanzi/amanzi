/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Base class for testing electromagnetics problems.
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_MAXWELL_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_MAXWELL_BASE_HH_

#include "Mesh.hh"
#include "Tensor.hh"
#include "WhetStoneFunction.hh"

class AnalyticElectromagneticsBase : public Amanzi::WhetStone::WhetStoneFunction {
 public:
  AnalyticElectromagneticsBase(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) : mesh_(mesh) {};
  ~AnalyticElectromagneticsBase() {};

  // analytic solution for Maxwell's equations
  // -- resitivity tensor T
  virtual Amanzi::WhetStone::Tensor Tensor(const Amanzi::AmanziGeometry::Point& p, double t) = 0;
  // -- analytic solution E
  virtual Amanzi::AmanziGeometry::Point electric_exact(const Amanzi::AmanziGeometry::Point& p, double t) = 0;
  virtual Amanzi::AmanziGeometry::Point magnetic_exact(const Amanzi::AmanziGeometry::Point& p, double t) = 0;
  // -- source term
  virtual Amanzi::AmanziGeometry::Point source_exact(const Amanzi::AmanziGeometry::Point& p, double t) = 0;

  // interface to function (dummy implementation)
  virtual double Value(const Amanzi::AmanziGeometry::Point& p, double t) const override { return 0.0; }

  // error calculation
  void ComputeFaceError(Epetra_MultiVector& u, double t, double& unorm, double& l2_err, double& inf_err);
  void ComputeEdgeError(Epetra_MultiVector& u, double t, double& unorm, double& l2_err, double& inf_err);
  void ComputeNodeError(Epetra_MultiVector& u, double t, double& unorm, double& l2_err, double& inf_err);

  // communications
  void GlobalOp(std::string op, double* val, int n);

 protected:
  Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh_;
};

#endif

