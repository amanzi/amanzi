/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)

*/

//! Simple test harness for a prescribed matrix with preconditioner.

#pragma once

#include "AmanziTypes.hh"
#include "AmanziComm.hh"

#include "Preconditioner.hh"
#include "PreconditionerFactory.hh"

using namespace Amanzi;

struct Matrix {

  explicit Matrix(const Map_ptr_type& map)
      : map_(map) {

    A_ = Teuchos::rcp(new Matrix_type(map_, map_, 3));

    double v0[2] = {1.0, -1.0};
    int inds0[2] = {0,1};
    A_->insertLocalValues(0, 2, v0, inds0);

    for (int i = 1; i < map_->getNodeNumElements()-1; i++) {
      int indices[3];
      double values[3] = { double(-i), double(2 * i + 1), double(-i - 1) };
      for (int k = 0; k < 3; k++) indices[k] = i + k - 1;
      A_->insertLocalValues(i, 3, values, indices);
    }

    int i = map_->getNodeNumElements()-1;
    double vN[2] = {double(-i), double(2*i+1)};
    int indsN[2] = {i-1,i};
    A_->insertLocalValues(i, 2, vN, indsN);

    A_->fillComplete(map_, map_);
  };

  explicit Matrix(const int N)
      : Matrix(Teuchos::rcp(new Map_type(N,0,Amanzi::getDefaultComm()))) {}
  
  void Init(const std::string& name,
            const ParameterList_ptr_type& plist)
  {
    plist->set("preconditioner type", name);
    AmanziPreconditioners::PreconditionerFactory<Matrix, Vector_type> fac;
    pc = fac.Create(plist);
    pc->Update(Teuchos::rcpFromRef(*this));
  }

  Matrix_ptr_type A() { return A_; }
  void AssembleMatrix() {}
  Map_ptr_type getDomainMap() const { return map_; }
  Map_ptr_type getRangeMap() const { return map_; }
  Map_ptr_type getRowMap() const { return map_; }
  Map_ptr_type getColumnMap() const { return map_; }

  void getLocalDiagCopy(Vector_type& diag) {
    A_->getLocalDiagCopy(diag);
  }
  
  int apply(const Vector_type& x, Vector_type& y) const {
    A_->apply(x,y);
    return 0;
  }

  int applyInverse(const Vector_type& y, Vector_type& x) const {
    pc->applyInverse(y,x);
    return 0;
  }

 protected:
  Map_ptr_type map_;
  Teuchos::RCP<AmanziPreconditioners::Preconditioner<Matrix, Vector_type>> pc;
  Matrix_ptr_type A_;

};
  
