/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/
//! Base class for providing Apply() and ApplyInverse() methods.

#pragma once

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "dbc.hh"
#include "errors.hh"

namespace Amanzi {

template<class Vector,
         class VectorSpace=typename Vector::VectorSpace_t>
class Matrix {
 public:
  using Vector_t = Vector;
  using VectorSpace_t = VectorSpace;
  
  virtual ~Matrix() = default;

  virtual int Apply(const Vector& x, Vector& y) const = 0;

  virtual void InitializeInverse(Teuchos::ParameterList& inv_list) = 0;
  virtual void UpdateInverse() = 0;
  virtual void ComputeInverse() = 0;
  virtual int ApplyInverse(const Vector& y, Vector& x) const = 0;
  
  virtual const VectorSpace& DomainMap() const = 0;
  virtual const VectorSpace& RangeMap() const = 0;

  //  virtual double TrueResidual(const Vector& x, const Vector& y) const = 0;

  // control and statistics -- must be valid for both iterative and
  // non-iterative methods, approximate and exact methods.
  virtual double residual() const = 0;
  virtual int num_itrs() const = 0;
  //  virtual void add_criteria(int criteria) = 0;

  virtual int returned_code() const = 0;
  virtual std::string returned_code_string() const = 0;

  virtual std::string name() const = 0;

};

}  // namespace Amanzi
 

               

