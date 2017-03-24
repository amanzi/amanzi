/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Derived class for transport sources and sinks.
*/

#ifndef AMANZI_TRANSPORT_SOURCE_FUNCTION_HH_
#define AMANZI_TRANSPORT_SOURCE_FUNCTION_HH_

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"

#include "CommonDefs.hh"
#include "Mesh.hh"
#include "PK_DomainFunction.hh"

namespace Amanzi {
namespace Transport {

class TransportSourceFunction : public PK_DomainFunction {
 public:
  TransportSourceFunction() {};
  TransportSourceFunction(const Teuchos::ParameterList& plist) {};

  // access
  const std::string& tcc_name() { return tcc_name_; }
  int tcc_index() { return tcc_index_; }
  void set_tcc_name(const std::string& name) { tcc_name_ = name; }
  void set_tcc_index(int i) { tcc_index_ = i; }
  virtual void set_state(const Teuchos::RCP<State>& S) {S_=S;}
  virtual std::string name() { return "TransportSourceFunction\n";}

 private:
  std::string tcc_name_;
  int tcc_index_; // index the global list of components
  Teuchos::RCP<const State> S_;
};

}  // namespace Transport
}  // namespace Amanzi

#endif
