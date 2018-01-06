/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Base class for testing DG schemes.
  List of problems:

  AnalyticDG0n: polynomial solution of order n, where n=0,1,2, and 3
  AnalyticDG04: sin(3x) six(6y)
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_DG_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_DG_BASE_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"
#include "Polynomial.hh"

class AnalyticDGBase {
 public:
  AnalyticDGBase(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, int order)
    : mesh_(mesh),
      order_(order),
      d_(mesh_->space_dimension()) {};
  ~AnalyticDGBase() {};

  // analytic solution in conventional Taylor basis
  virtual void TaylorCoefficients(const Amanzi::AmanziGeometry::Point& p, double t,
                                  Amanzi::WhetStone::Polynomial& coefs) = 0;

  double function_exact(const Amanzi::AmanziGeometry::Point& p, double t) {
    Amanzi::WhetStone::Polynomial coefs;
    TaylorCoefficients(p, t, coefs);
    return coefs(0, 0);
  }

 protected:
  Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh_;
  int order_, d_;
};

#endif

