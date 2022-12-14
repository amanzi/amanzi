/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Process Kernels

  Prototype class for computing boundary and source terms. PK's classes
  for source and boundary terms may inherit from this class or mimic
  its functionality.
*/

#ifndef AMANZI_PK_DOMAIN_FUNCTION_HH_
#define AMANZI_PK_DOMAIN_FUNCTION_HH_

#include <string>
#include <vector>

#include "Epetra_Vector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

namespace Amanzi {

class PK_DomainFunction {
 public:
  PK_DomainFunction() : domain_volume_(-1.0){};

  PK_DomainFunction(const Teuchos::ParameterList& plist) : domain_volume_(-1.0){};

  virtual ~PK_DomainFunction() = default;

  // source term on time interval (t0, t1]
  virtual void Compute(double t0, double t1) { AMANZI_ASSERT(false); }

  // model name
  virtual std::string name() const { return "undefined"; }

  // access
  // -- volume of the regions
  double domain_volume() { return domain_volume_; }
  // -- nick-name of the function
  std::string keyword() { return keyword_; }

  // iterator methods
  typedef std::map<int, std::vector<double>>::iterator Iterator;
  typename std::map<int, std::vector<double>>::iterator begin() { return value_.begin(); }
  typename std::map<int, std::vector<double>>::iterator end() { return value_.end(); }
  typename std::map<int, std::vector<double>>::size_type size() { return value_.size(); }

 protected:
  std::map<int, std::vector<double>> value_;
  std::map<int, double> linear_term_;
  double domain_volume_;
  std::string keyword_;
};

} // namespace Amanzi

#endif
