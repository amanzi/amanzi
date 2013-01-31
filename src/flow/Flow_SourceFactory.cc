/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Konstantin Lipnikov (version 2)  (lipnikov@lanl.gov)
*/

#include <string>
#include <vector>

#include "FunctionFactory.hh"
#include "errors.hh"

#include "Flow_SourceFactory.hh"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Process source, step 1.
****************************************************************** */
DomainFunction* FlowSourceFactory::createSource() const
{
  DomainFunction* src = new DomainFunction(mesh_);

  // Iterate through the source specification sublists in the params_.
  // All are expected to be sublists of identical structure.
  for (Teuchos::ParameterList::ConstIterator i = params_->begin(); i != params_->end(); ++i) {
    std::string name = i->first;
    if (params_->isSublist(name)) {
      Teuchos::ParameterList& spec = params_->sublist(name);
      try {
        ProcessSourceSpec(spec, src);
      } catch (Errors::Message& msg) {
        Errors::Message m;
        m << "in sublist \"" << spec.name().c_str() << "\": " << msg.what();
        Exceptions::amanzi_throw(m);
      }
    } else {
      Errors::Message m;
      m << "parameter \"" << name.c_str() << "\" is not a sublist";
      Exceptions::amanzi_throw(m);
    }
  }
  return src;
}


/* ******************************************************************
* Process source, step 2.
****************************************************************** */
void FlowSourceFactory::ProcessSourceSpec(Teuchos::ParameterList& list, DomainFunction* src) const
{
  Errors::Message m;
  std::vector<std::string> regions;

  if (list.isParameter("regions")) {
    if (list.isType<Teuchos::Array<std::string> >("regions")) {
      regions = list.get<Teuchos::Array<std::string> >("regions").toVector();
    } else {
      m << "parameter \"regions\" is not of type \"Array string\"";
      Exceptions::amanzi_throw(m);
    }
  } else {
    m << "parameter \"regions\" is missing";
    Exceptions::amanzi_throw(m);
  }

  Teuchos::ParameterList* f_list;
  if (list.isSublist("sink")) {
    f_list = &list.sublist("sink");
  } else {
    m << "parameter \"sink\" is not a sublist";
    Exceptions::amanzi_throw(m);
  }

  // Make the source function.
  Teuchos::RCP<Function> f;
  FunctionFactory f_fact;
  try {
    f = Teuchos::rcp(f_fact.Create(*f_list));
  } catch (Errors::Message& msg) {
    m << "error in source sublist : " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  // Add this source specification to the domain function.
  int method;
  std::string action_name = list.get<std::string>("spatial distribution method", "none");
  ProcessStringActions(action_name, &method);
  
  src->Define(regions, f, method);
}


/* ****************************************************************
* Process string for a source specipic action.
**************************************************************** */
void FlowSourceFactory::ProcessStringActions(const std::string& name, int* method) const
{
  Errors::Message msg;
  if (name == "none") {
    *method = Amanzi::DOMAIN_FUNCTION_ACTION_NONE;
  } else if (name == "volume") {
    *method = Amanzi::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_VOLUME;
  } else if (name == "permeability") {
    *method = Amanzi::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY;
  } else {
    msg << "Flow PK: unknown source distribution method has been specified.";
    Exceptions::amanzi_throw(msg);
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi
