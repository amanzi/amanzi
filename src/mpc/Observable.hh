

#ifndef AMANZI_OBSERVABLE_HH
#define AMANZI_OBSERVABLE_HH

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "IOEvent.hh"
#include "GeometryDefs.hh"
#include "State.hh"
#include "Units.hh"

namespace Amanzi{

class Observable : public IOEvent{
 public:
  Observable(std::string variable,
	     std::string region,
	     std::string functional,
	     Teuchos::ParameterList& plist,
	     Teuchos::ParameterList& units_plist,
	     Teuchos::RCP<AmanziMesh::Mesh> mesh):
    variable_(variable), region_(region),
    functional_(functional), plist_(plist),
    mesh_(mesh), IOEvent(plist)
  {
    ReadParameters_();
    
    Errors::Message msg;

    units_.Init(units_plist);

    // for now we can only observe Integrals and Values
    if (functional_ != "observation data: integral"  &&
        functional_ != "observation data: point" )  {
      msg << "FlexibleObservations: can only handle Functional == observation data:"
          << " integral, or functional == observation data: point";
      Exceptions::amanzi_throw(msg);
    }
  }


  virtual void ComputeObservation(State& S, double* value, double* volume) { assert(false); }

  virtual int ComputeRegionSize(){return region_size_;};

  // protected:    

  std::string variable_;
  std::string region_;
  std::string functional_;
  const Teuchos::ParameterList& plist_;
  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
  double volume_;
  int region_size_;
  AmanziMesh::Entity_ID_List entity_ids_;
  Utils::Units units_;
};

} // namespace Amanzi

#endif
