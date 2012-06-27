/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

------------------------------------------------------------------------- */

#include "function-factory.hh"
#include "vector_function.hh"
#include "composite_function.hh"
#include "errors.hh"

#include "energy_bc_factory.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Process Dirichet BC (temperature), step 1.
****************************************************************** */
Teuchos::RCP<BoundaryFunction> EnergyBCFactory::CreateTemperature() const {
  Teuchos::RCP<BoundaryFunction> bc = Teuchos::rcp(new BoundaryFunction(mesh_));
  try {
    ProcessTemperatureList_(plist_.sublist("temperature"), bc);
  } catch (Errors::Message& msg) {
    Errors::Message message;
    message << "EnergyBCFactory: \"temperature\" sublist error: " << msg.what();
    Exceptions::amanzi_throw(message);
  } catch (Teuchos::Exceptions::InvalidParameterType& msg) {
    Errors::Message message;
    message << "EnergyBCFactory: \"temperature\" sublist error: not a sublist: "
            << msg.what();
    Exceptions::amanzi_throw(message);
  }
  return bc;
};


/* ******************************************************************
* Process Neumann BC (enthalpy flux), step 1.
****************************************************************** */
Teuchos::RCP<BoundaryFunction> EnergyBCFactory::CreateEnthalpyFlux() const {
  Teuchos::RCP<BoundaryFunction> bc = Teuchos::rcp(new BoundaryFunction(mesh_));
  try {
    ProcessEnthalpyFluxList_(plist_.sublist("enthalpy flux"), bc);
  } catch (Errors::Message& msg) {
    Errors::Message message;
    message << "EnergyBCFactory: \"enthalpy flux\" sublist error: " << msg.what();
    Exceptions::amanzi_throw(message);
  } catch (Teuchos::Exceptions::InvalidParameterType& msg) {
    Errors::Message message;
    message << "EnergyBCFactory: \"enthalpy flux\" sublist error: not a sublist: " << msg.what();
    Exceptions::amanzi_throw(message);
  }
  return bc;
};


/* ******************************************************************
* Process Dirichet BC (temperature), step 2.
****************************************************************** */
void EnergyBCFactory::ProcessTemperatureList_(const Teuchos::ParameterList& list,
                                        const Teuchos::RCP<BoundaryFunction>& bc) const {
  // Iterate through the BC specification sublists in the list.
  // All are expected to be sublists of identical structure.
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i) {
    std::string name = i->first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList spec = list.sublist(name);
      try {
        ProcessTemperatureSpec_(spec, bc);
      } catch (Errors::Message& msg) {
        Errors::Message message;
        message << "in sublist \"" << spec.name().c_str() << "\": " << msg.what();
        Exceptions::amanzi_throw(message);
      }
    } else { // ERROR -- parameter is not a sublist
      Errors::Message message;
      message << "parameter \"" << name.c_str() << "\" is not a sublist";
      Exceptions::amanzi_throw(message);
    }
  }
};


/* ******************************************************************
* Process Dirichet BC (temperature), step 3.
****************************************************************** */
void EnergyBCFactory::ProcessTemperatureSpec_(const Teuchos::ParameterList& list,
        const Teuchos::RCP<BoundaryFunction>& bc) const {
  Errors::Message message;
  std::vector<std::string> regions;

  if (list.isParameter("regions")) {
    if (list.isType<Teuchos::Array<std::string> >("regions")) {
      regions = list.get<Teuchos::Array<std::string> >("regions").toVector();
    } else {
      message << "parameter \"regions\" is not of type \"Array string\"";
      Exceptions::amanzi_throw(message);
    }
  } else {  // Parameter "regions" is missing.
    message << "parameter \"regions\" is missing";
    Exceptions::amanzi_throw(message);
  }

  Teuchos::ParameterList f_list;
  if (list.isParameter("boundary temperature")) {
    if (list.isSublist("boundary temperature")) {
      // validate function-factory sublist
      f_list = list.sublist("boundary temperature");
    } else {  // "boundary temperature" is not a sublist
      message << "parameter \"boundary temperature\" is not a sublist";
      Exceptions::amanzi_throw(message);
    }
  } else {  // "boundary temperature" sublist is missing.
    message << "sublist \"boundary temperature\" is missing";
    Exceptions::amanzi_throw(message);
  }

  // Make the boundary temperature function.
  Teuchos::RCP<Function> f;
  FunctionFactory f_fact;
  try {
    f = Teuchos::rcp(f_fact.Create(f_list));
  } catch (Errors::Message& msg) {
    message << "error in sublist \"boundary temperature\": " << msg.what();
    Exceptions::amanzi_throw(message);
  }

  // A bit hacky -- this entire code needs to be revisited in light of the new
  // options for mesh_functions.
  Teuchos::RCP<VectorFunction> func = Teuchos::rcp(new CompositeFunction(f));

  // Add this BC specification to the boundary function.
  bc->Define(regions, func);
};


/* ******************************************************************
* Process Neumann BC (enthalpy flux), step 2.
****************************************************************** */
void EnergyBCFactory::ProcessEnthalpyFluxList_(const Teuchos::ParameterList& list,
        const Teuchos::RCP<BoundaryFunction>& bc) const {
  // Iterate through the BC specification sublists in the list.
  // All are expected to be sublists of identical structure.
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i) {
    std::string name = i->first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList spec = list.sublist(name);
      try {
        ProcessEnthalpyFluxSpec_(spec, bc);
      } catch (Errors::Message& msg) {
        Errors::Message message;
        message << "in sublist \"" << spec.name().c_str() << "\": " << msg.what();
        Exceptions::amanzi_throw(message);
      }
    } else {  // Parameter is not a sublist
      Errors::Message message;
      message << "parameter \"" << name.c_str() << "\" is not a sublist";
      Exceptions::amanzi_throw(message);
    }
  }
}


/* ******************************************************************
* Process Neumann BC (enthalpy flux), step 3.
****************************************************************** */
void EnergyBCFactory::ProcessEnthalpyFluxSpec_(const Teuchos::ParameterList& list,
        const Teuchos::RCP<BoundaryFunction>& bc) const {
  Errors::Message message;
  // Get the regions parameter value.
  std::vector<std::string> regions;
  if (list.isParameter("regions")) {
    if (list.isType<Teuchos::Array<std::string> >("regions")) {
      regions = list.get<Teuchos::Array<std::string> >("regions").toVector();
    } else {
      message << "parameter \"regions\" is not of type \"Array string\"";
      Exceptions::amanzi_throw(message);
    }
  } else {  // Parameter "regions" is missing.
    message << "parameter \"regions\" is missing";
    Exceptions::amanzi_throw(message);
  }

  Teuchos::ParameterList f_list;
  if (list.isParameter("outward enthalpy flux")) {
    if (list.isSublist("outward enthalpy flux")) {
      // validate function-factory sublist
      f_list = list.sublist("outward enthalpy flux");
    } else {  // "outward enthalpy flux" is not a sublist.
      message << "parameter \"outward enthalpy flux\" is not a sublist";
      Exceptions::amanzi_throw(message);
    }
  } else {  // "outward enthalpy flux" sublist is missing
    message << "sublist \"outward enthalpy flux\" is missing";
    Exceptions::amanzi_throw(message);
  }

  // Make the outward enthalpy flux function.
  Teuchos::RCP<Function> f;
  FunctionFactory f_fact;
  try {
    f = Teuchos::rcp(f_fact.Create(f_list));
  } catch (Errors::Message& msg) {
    message << "error in sublist \"outward enthalpy flux\": " << msg.what();
    Exceptions::amanzi_throw(message);
  }

  // A bit hacky -- this entire code needs to be revisited in light of the new
  // options for mesh_functions.
  Teuchos::RCP<VectorFunction> func = Teuchos::rcp(new CompositeFunction(f));

  // Add this BC specification to the boundary function.
  bc->Define(regions, func);
};

}  // namespace Energy
}  // namespace Amanzi
