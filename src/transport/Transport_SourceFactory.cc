/*
This is the transport component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <string>
#include <vector>

#include "function-factory.hh"
#include "errors.hh"

#include "Transport_SourceFactory.hh"

namespace Amanzi {
namespace AmanziTransport {

/* ******************************************************************
* Process source, step 1.
****************************************************************** */
DomainFunction* TransportSourceFactory::CreateSource()
{
  DomainFunction* src = new DomainFunction(mesh_);

  // Iterate through the source specification sublists in the params_.
  // All are expected to be sublists of identical structure.
  for (Teuchos::ParameterList::ConstIterator i = params_->begin(); i != params_->end(); ++i) {
    std::string name = i->first;
    if (params_->isSublist(name)) {
      Teuchos::ParameterList& spec = params_->sublist(name);
      try {
        ProcessSourceSpec(spec, name, src);
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
void TransportSourceFactory::ProcessSourceSpec(
    Teuchos::ParameterList& list, const std::string& name, UniqueMeshFunction* src) const
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
  
  src->DefineMultiValue(regions, f, method, name);
}


/* ****************************************************************
* Process string for a source specipic action.
**************************************************************** */
void TransportSourceFactory::ProcessStringActions(const std::string& name, int* method) const
{
  Errors::Message msg;
  if (name == "none") {
    *method = Amanzi::DOMAIN_FUNCTION_ACTION_NONE;
  } else if (name == "volume") {
    *method = Amanzi::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_VOLUME;
  } else if (name == "permeability") {
    *method = Amanzi::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY;
  } else {
    msg << "Transport PK: unknown source distribution method has been specified.";
    Exceptions::amanzi_throw(msg);
  }
}

}  // namespace AmanziTransport
}  // namespace Amanzi
