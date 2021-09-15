/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson
           Ethan Coon
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_TRANSPORT_DOMAIN_FUNCTION_HH_
#define AMANZI_TRANSPORT_DOMAIN_FUNCTION_HH_

#include <vector>
#include <map>
#include <string>

#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"
#include "State.hh"

#include "CommonDefs.hh"
#include "Mesh.hh"

namespace Amanzi {
namespace Transport {

class TransportDomainFunction {
 public:
  TransportDomainFunction() : domain_volume_(-1.0) {};
  TransportDomainFunction(const Teuchos::ParameterList& plist) : domain_volume_(-1.0) {};
  ~TransportDomainFunction() {};

  // source term on time interval (t0, t1]
  virtual void Compute(double t0, double t1) { AMANZI_ASSERT(false); }
  virtual void ComputeSubmodel(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                               Teuchos::RCP<CompositeVector> tcc) {};

  // model name
  virtual std::string name() const { return "undefined"; } 

  // access
  // -- volume of the regions
  double domain_volume() { return domain_volume_; }
  // -- nick-name of the function
  std::string keyword() { return keyword_; }

  std::vector<std::string>& tcc_names() { return tcc_names_; }
  void set_tcc_names(const std::vector<std::string>& tcc_names){ tcc_names_ = tcc_names; }
  std::vector<int>& tcc_index() { return tcc_index_; }
  virtual void set_state(const Teuchos::RCP<State>& S) { S_ = S; }

  // iterator methods
  typedef std::map<int, std::vector<double> >::iterator Iterator;
  Iterator begin() { return value_.begin(); }
  Iterator end() { return value_.end(); }
  std::map<int, std::vector<double> >::size_type size() { return value_.size(); }

  // derivatives
  const std::map<int, double>& linear_term() const { return linear_term_; }

 protected:
  double domain_volume_;
  std::map<int, std::vector<double> > value_;  // tcc values on boundary faces or 
                                               // src values in domain cells
  std::map<int, double> linear_term_;  // linearized term, e.g. [mol / s] for sources

  std::string keyword_;
  Teuchos::RCP<const State> S_;

  std::vector<std::string> tcc_names_;  // list of component names
  std::vector<int> tcc_index_;  // index of component in the global list
};

}  // namespace Transport
}  // namespace Amanzi

#endif
