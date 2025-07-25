/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon
           Konstantin Lipnikov
           Daniil Svyatsky
*/

/*
  Multi-Process Coordinator

  Virtual base class.
*/

#pragma once

#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"

#include "IOEvent.hh"
#include "GeometryDefs.hh"
#include "State.hh"
#include "Units.hh"
#include "Key.hh"

namespace Amanzi {

class ObservableAmanzi : public Utils::IOEvent {
 public:
  ObservableAmanzi(std::string variable,
                   std::string region,
                   std::string functional,
                   Teuchos::ParameterList& plist,
                   Teuchos::ParameterList& units_plist,
                   Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : Utils::IOEvent(plist),
      variable_(variable),
      functional_(functional),
      sum_(0.0),
      mesh_(mesh),
      region_(region)
  {
    ReadParameters_();

    units_.Init(units_plist);

    domain_ = plist.get<std::string>("domain name", "domain");

    // for now we can only observe Integrals and Values
    if (functional_ != "observation data: integral" && functional_ != "observation data: point") {
      Errors::Message msg;
      msg << "FlexibleObservations: can only handle Functional == observation data:"
          << " integral, or functional == observation data: point";
      Exceptions::amanzi_throw(msg);
    }
  }

  virtual ~ObservableAmanzi() = default;

  virtual void ComputeObservation(State& S,
                                  double* value,
                                  double* volume,
                                  std::string& unit,
                                  double dt) = 0;
  virtual int ComputeRegionSize() { return region_size_; }

 public:
  std::string variable_;
  std::string functional_;
  double sum_;

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  AmanziMesh::cEntity_ID_View entity_ids_;
  int region_size_;

  std::string region_;
  double volume_;
  Key domain_;
};

} // namespace Amanzi
