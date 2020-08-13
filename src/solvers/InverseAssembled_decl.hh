/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/
//! Provides ApplyInverse() using assembled methods.

#pragma once

#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"

#include "SuperMap.hh"
#include "InverseHelpers.hh"
#include "Inverse.hh"

namespace Amanzi {
namespace AmanziSolvers {


//
// Class for assembled inverse methods.
//
template<class Operator,
         class Assembler=Operator,
         class Vector=typename Operator::Vector_t,
         class VectorSpace=typename Vector::VectorSpace_t>
class InverseAssembled :
      public Inverse<Operator,Assembler,Vector,VectorSpace> {
 public:
  InverseAssembled(const std::string& method_name) :
      method_name_(method_name),
      updated_(false),
      computed_once_(false),
      Inverse<Operator,Assembler,Vector,VectorSpace>()
  {}

  virtual void InitializeInverse(Teuchos::ParameterList& plist) override final;

  virtual void UpdateInverse() override final;
  virtual void ComputeInverse() override final;
  virtual int ApplyInverse(const Vector& X, Vector& Y) const override final;

  virtual double residual() const override final {
    return solver_->residual();
  }

  virtual int num_itrs() const override final {
    return solver_->num_itrs();
  }

  virtual void add_criteria(int criteria) override final {
    return solver_->add_criteria(criteria);
  }

  virtual int returned_code() const override final {
    return solver_->returned_code();
  }

  virtual std::string returned_code_string() const override final {
    return solver_->returned_code_string();
  }

 protected:
  bool updated_, computed_once_;
  std::string method_name_;

  using Inverse<Operator,Assembler,Vector,VectorSpace>::m_;
  using Inverse<Operator,Assembler,Vector,VectorSpace>::h_;

  Teuchos::RCP<Inverse<Epetra_CrsMatrix,Epetra_CrsMatrix,Epetra_Vector,Epetra_Map>> solver_;
  Teuchos::RCP<Operators::SuperMap> smap_;
  mutable Teuchos::RCP<Epetra_Vector> x_, y_;
};
  
    
} // namespace AmanziSolvers
} // namespace Amanzi
