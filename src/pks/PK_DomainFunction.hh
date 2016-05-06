/*
  Process Kernels

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Base class for computing source terms. PK's classes for source
  terms should inherit from this class.
*/

#ifndef AMANZI_PK_DOMAIN_FUNCTION_HH_
#define AMANZI_PK_DOMAIN_FUNCTION_HH_

#include <string>
#include <vector>

#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"

namespace Amanzi {

class PK_DomainFunction {
 public:
  PK_DomainFunction() {};
  virtual ~PK_DomainFunction() {};

  // source term on time interval (t0, t1]
  virtual void Compute(double t0, double t1) { ASSERT(false); }

  // model name
  virtual std::string name() { return "undefined"; } 

  // iterator methods
  typedef std::map<int, double>::const_iterator Iterator;
  Iterator begin() const { return value_.begin(); }
  Iterator end() const { return value_.end(); }
  std::map<int, double>::size_type size() { return value_.size(); }

 protected:
  std::map<int, double> value_;
};

}  // namespace Amanzi

#endif

