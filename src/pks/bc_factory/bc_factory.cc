/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (versions 1 & 2)  (nnc@lanl.gov)
         Ethan Coon (ATS version)
*/

#include "constant-function.hh"
#include "function-factory.hh"
#include "vector_function.hh"
#include "composite_function.hh"
#include "errors.hh"

#include "bc_factory.hh"

namespace Amanzi {

/* ******************************************************************
* Process Dirichet BC (pressure), step 1.
****************************************************************** */
Teuchos::RCP<Functions::BoundaryFunction>
BCFactory::CreateWithFunction(std::string list_name, std::string function_name) const {
  Teuchos::RCP<Functions::BoundaryFunction> bc =
      Teuchos::rcp(new Functions::BoundaryFunction(mesh_));

  if (plist_.isParameter(list_name)) {
    if (plist_.isSublist(list_name)) {
      ProcessListWithFunction_(plist_.sublist(list_name), function_name, bc);
    } else {
      std::stringstream mstream;
      mstream << "BCFactory: " << list_name << " sublist error: not a sublist.";
      Errors::Message m(mstream.str());
      Exceptions::amanzi_throw(m);
    }
  }
  return bc;
}

/* ******************************************************************
* Process Dirichet BC (pressure), step 1.
****************************************************************** */
Teuchos::RCP<Functions::BoundaryFunction>
BCFactory::CreateWithoutFunction(std::string list_name) const {
  Teuchos::RCP<Functions::BoundaryFunction> bc =
      Teuchos::rcp(new Functions::BoundaryFunction(mesh_));

  if (plist_.isParameter(list_name)) {
    if (plist_.isSublist(list_name)) {
      ProcessListWithoutFunction_(plist_.sublist(list_name), bc);
    } else {
      std::stringstream mstream;
      mstream << "BCFactory: " << list_name << " sublist error: not a sublist.";
      Errors::Message m(mstream.str());
      Exceptions::amanzi_throw(m);
    }
  }
  return bc;
}


/* ******************************************************************
* Process Dirichet BC (pressure), step 2.
****************************************************************** */
void BCFactory::ProcessListWithFunction_(const Teuchos::ParameterList& list,
        std::string function_name,
        const Teuchos::RCP<Functions::BoundaryFunction>& bc) const {
  // Iterate through the BC specification sublists in the list.
  // All are expected to be sublists of identical structure.
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i) {
    std::string name = i->first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList spec = list.sublist(name);
      try {
        ProcessSpecWithFunction_(spec, function_name, bc);
      } catch (Errors::Message& msg) {
        Errors::Message m;
        m << "in sublist " << spec.name().c_str() << ": " << msg.what();
        Exceptions::amanzi_throw(m);
      }
    } else { // ERROR -- parameter is not a sublist
      Errors::Message m;
      m << "parameter " << name.c_str() << " is not a sublist";
      Exceptions::amanzi_throw(m);
    }
  }
}


/* ******************************************************************
* Process Dirichet BC (pressure), step 3.
****************************************************************** */
void BCFactory::ProcessSpecWithFunction_(const Teuchos::ParameterList& list,
        std::string function_name,
        const Teuchos::RCP<Functions::BoundaryFunction>& bc) const {
  std::stringstream m;
  std::vector<std::string> regions;

  if (list.isParameter("regions")) {
    if (list.isType<Teuchos::Array<std::string> >("regions")) {
      regions = list.get<Teuchos::Array<std::string> >("regions").toVector();
    } else {
      m << "parameter \"regions\" is not of type Array string";
      Errors::Message msg(m.str());
      Exceptions::amanzi_throw(msg);
    }
  } else {  // Parameter "regions" is missing.
    m << "parameter \"regions\" is missing";
    Errors::Message msg(m.str());
    Exceptions::amanzi_throw(msg);
  }

  Teuchos::ParameterList f_list;
  if (list.isParameter(function_name)) {
    if (list.isSublist(function_name)) {
      // validate function-factory sublist
      f_list = list.sublist(function_name);
    } else {  // function_name is not a sublist
      m << "parameter " << function_name << " is not a sublist";
      Errors::Message msg(m.str());
      Exceptions::amanzi_throw(msg);
    }
  } else {  // function_name sublist is missing.
    m << "sublist " << function_name << " is missing";
    Errors::Message msg(m.str());
    Exceptions::amanzi_throw(msg);
  }

  // Make the boundary pressure function.
  Teuchos::RCP<Function> f;
  FunctionFactory f_fact;
  try {
    f = Teuchos::rcp(f_fact.Create(f_list));
  } catch (Errors::Message& msg) {
    m << "error in sublist " << function_name << ": " << msg.what();
    Errors::Message msg(m.str());
    Exceptions::amanzi_throw(msg);
  }

  // A bit hacky -- this entire code needs to be revisited in light of the new
  // options for mesh_functions.
  Teuchos::RCP<VectorFunction> func = Teuchos::rcp(new CompositeFunction(f));

  // Add this BC specification to the boundary function.
  bc->Define(regions, func);
}


/* ******************************************************************
* Process Dirichet BC (pressure), step 2.
****************************************************************** */
void BCFactory::ProcessListWithoutFunction_(const Teuchos::ParameterList& list,
        const Teuchos::RCP<Functions::BoundaryFunction>& bc) const {
  // Iterate through the BC specification sublists in the list.
  // All are expected to be sublists of identical structure.
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i) {
    std::string name = i->first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList spec = list.sublist(name);
      try {
        ProcessSpecWithoutFunction_(spec, bc);
      } catch (Errors::Message& msg) {
        std::stringstream m;
        m << "in sublist " << spec.name().c_str() << ": " << msg.what();
        Errors::Message msg(m.str());
        Exceptions::amanzi_throw(msg);
      }
    } else { // ERROR -- parameter is not a sublist
      std::stringstream m;
      m << "parameter " << name.c_str() << " is not a sublist";
      Errors::Message msg(m.str());
      Exceptions::amanzi_throw(msg);
    }
  }
}


/* ******************************************************************
* Process Dirichet BC (pressure), step 3.
****************************************************************** */
void BCFactory::ProcessSpecWithoutFunction_(const Teuchos::ParameterList& list,
        const Teuchos::RCP<Functions::BoundaryFunction>& bc) const {
  Errors::Message m;
  std::vector<std::string> regions;

  if (list.isParameter("regions")) {
    if (list.isType<Teuchos::Array<std::string> >("regions")) {
      regions = list.get<Teuchos::Array<std::string> >("regions").toVector();
    } else {
      m << "parameter \"regions\" is not of type \"Array string\"";
      Exceptions::amanzi_throw(m);
    }
  } else {  // Parameter "regions" is missing.
    m << "parameter \"regions\" is missing";
    Exceptions::amanzi_throw(m);
  }

  // Make the boundary pressure function.
  Teuchos::RCP<Function> f = Teuchos::rcp(new ConstantFunction(0.));

  // A bit hacky -- this entire code needs to be revisited in light of the new
  // options for mesh_functions.
  Teuchos::RCP<VectorFunction> func = Teuchos::rcp(new CompositeFunction(f));

  // Add this BC specification to the boundary function.
  bc->Define(regions, func);
}

}  // namespace Amanzi
