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
namespace TestPKs {

/* ******************************************************************
* Process Dirichet BC (dirichlet), step 1.
****************************************************************** */
Teuchos::RCP<Functions::BoundaryFunction> BCFactory::CreateDirichlet() const {
  Teuchos::RCP<Functions::BoundaryFunction> bc = Teuchos::rcp(new Functions::BoundaryFunction(mesh_));
  if (plist_.isParameter("dirichlet")) {
    if (plist_.isSublist("dirichlet")) {
      ProcessDirichletList(plist_.sublist("dirichlet"), bc);
    } else {
      Errors::Message m;
      m << "BCFactory: \"dirichlet\" sublist error: not a sublist.";
      Exceptions::amanzi_throw(m);
    }
  }
  return bc;
}

/* ******************************************************************
* Process Neumann BC (flux), step 1.
****************************************************************** */
Teuchos::RCP<Functions::BoundaryFunction> BCFactory::CreateNeumann() const {
  Teuchos::RCP<Functions::BoundaryFunction> bc = Teuchos::rcp(new Functions::BoundaryFunction(mesh_));
  if (plist_.isParameter("mass flux")) {
    if (plist_.isSublist("mass flux")) {
      ProcessNeumannList(plist_.sublist("mass flux"), bc);
    } else {
      Errors::Message m;
      m << "BCFactory: \"mass flux\" sublist error: not a sublist.";
      Exceptions::amanzi_throw(m);
    }
  }
  return bc;
}

/* ******************************************************************
* Process Dirichet BC (dirichlet), step 2.
****************************************************************** */
void BCFactory::ProcessDirichletList(const Teuchos::ParameterList& list,
        const Teuchos::RCP<Functions::BoundaryFunction>& bc) const {
  // Iterate through the BC specification sublists in the list.
  // All are expected to be sublists of identical structure.
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i) {
    std::string name = i->first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList spec = list.sublist(name);
      try {
        ProcessDirichletSpec(spec, bc);
      } catch (Errors::Message& msg) {
        Errors::Message m;
        m << "in sublist \"" << spec.name().c_str() << "\": " << msg.what();
        Exceptions::amanzi_throw(m);
      }
    } else { // ERROR -- parameter is not a sublist
      Errors::Message m;
      m << "parameter \"" << name.c_str() << "\" is not a sublist";
      Exceptions::amanzi_throw(m);
    }
  }
}


/* ******************************************************************
* Process Dirichet BC (dirichlet), step 3.
****************************************************************** */
void BCFactory::ProcessDirichletSpec(const Teuchos::ParameterList& list,
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

  Teuchos::ParameterList f_list;
  if (list.isParameter("boundary value")) {
    if (list.isSublist("boundary value")) {
      // validate function-factory sublist
      f_list = list.sublist("boundary value");
    } else {  // "boundary dirichlet" is not a sublist
      m << "parameter \"boundary value\" is not a sublist";
      Exceptions::amanzi_throw(m);
    }
  } else {  // "boundary value" sublist is missing.
    m << "sublist \"boundary value\" is missing";
    Exceptions::amanzi_throw(m);
  }

  // Make the boundary value function.
  Teuchos::RCP<Function> f;
  FunctionFactory f_fact;
  try {
    f = Teuchos::rcp(f_fact.Create(f_list));
  } catch (Errors::Message& msg) {
    m << "error in sublist \"boundary value\": " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  // A bit hacky -- this entire code needs to be revisited in light of the new
  // options for mesh_functions.
  Teuchos::RCP<VectorFunction> func = Teuchos::rcp(new CompositeFunction(f));

  // Add this BC specification to the boundary function.
  bc->Define(regions, func);
}


/* ******************************************************************
* Process Neumann BC (mass flux), step 2.
****************************************************************** */
void BCFactory::ProcessNeumannList(const Teuchos::ParameterList& list,
                                        const Teuchos::RCP<Functions::BoundaryFunction>& bc) const {
  // Iterate through the BC specification sublists in the list.
  // All are expected to be sublists of identical structure.
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i) {
    std::string name = i->first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList spec = list.sublist(name);
      try {
        ProcessNeumannSpec(spec, bc);
      } catch (Errors::Message& msg) {
        Errors::Message m;
        m << "in sublist \"" << spec.name().c_str() << "\": " << msg.what();
        Exceptions::amanzi_throw(m);
      }
    } else {  // Parameter is not a sublist
      Errors::Message m;
      m << "parameter \"" << name.c_str() << "\" is not a sublist";
      Exceptions::amanzi_throw(m);
    }
  }
}


/* ******************************************************************
* Process Neumann BC (mass flux), step 3.
****************************************************************** */
void BCFactory::ProcessNeumannSpec(const Teuchos::ParameterList& list,
                                        const Teuchos::RCP<Functions::BoundaryFunction>& bc) const {
  Errors::Message m;
  // Get the regions parameter value.
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

  Teuchos::ParameterList f_list;
  if (list.isParameter("outward flux")) {
    if (list.isSublist("outward flux")) {
      // validate function-factory sublist
      f_list = list.sublist("outward flux");
    } else {  // "outward flux" is not a sublist.
      m << "parameter \"outward flux\" is not a sublist";
      Exceptions::amanzi_throw(m);
    }
  } else {  // "outward flux" sublist is missing
    m << "sublist \"outward flux\" is missing";
    Exceptions::amanzi_throw(m);
  }

  // Make the outward flux function.
  Teuchos::RCP<Function> f;
  FunctionFactory f_fact;
  try {
    f = Teuchos::rcp(f_fact.Create(f_list));
  } catch (Errors::Message& msg) {
    m << "error in sublist \"outward flux\": " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  // A bit hacky -- this entire code needs to be revisited in light of the new
  // options for mesh_functions.
  Teuchos::RCP<VectorFunction> func = Teuchos::rcp(new CompositeFunction(f));

  // Add this BC specification to the boundary function.
  bc->Define(regions, func);
}


}  // namespace
}  // namespace Amanzi
