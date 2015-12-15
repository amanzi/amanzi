/*
  Flow PK
 
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <string>
#include <vector>

#include "errors.hh"
#include "MultiFunction.hh"

#include "Flow_SourceFactory.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Process source, step 1.
****************************************************************** */
void FlowSourceFactory::Create(std::vector<FlowDomainFunction*>& srcs) const
{
  // Iterate through the source specification sublists in the plist_.
  // All are expected to be sublists of identical structure.
  for (Teuchos::ParameterList::ConstIterator i = plist_->begin(); i != plist_->end(); ++i) {
    std::string name = i->first;
    if (plist_->isSublist(name)) {
      Teuchos::ParameterList& spec = plist_->sublist(name);
      try {
        FlowDomainFunction* src = new FlowDomainFunction(mesh_);
        ProcessSourceSpec(spec, src);
        srcs.push_back(src);
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
}


/* ******************************************************************
* Process source, step 2.
****************************************************************** */
void FlowSourceFactory::ProcessSourceSpec(Teuchos::ParameterList& list, FlowDomainFunction* src) const
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
  Teuchos::RCP<Amanzi::MultiFunction> f;
  try {
    f = Teuchos::rcp(new MultiFunction(*f_list));
  } catch (Errors::Message& msg) {
    m << "error in source sublist : " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  // Add this source specification to the domain function.
  int method;
  std::string action_name = list.get<std::string>("spatial distribution method", "none");
  ProcessStringActions(action_name, &method);

  int submodel(CommonDefs::DOMAIN_FUNCTION_SUBMODEL_RATE);
  std::string submodel_name = list.get<std::string>("submodel", "rate");
  if (submodel_name == "integrated source") submodel = CommonDefs::DOMAIN_FUNCTION_SUBMODEL_INTEGRAL;

  src->Define(regions, f, method, submodel);
}


/* ****************************************************************
* Process string for a source specipic action.
**************************************************************** */
void FlowSourceFactory::ProcessStringActions(const std::string& name, int* method) const
{
  Errors::Message msg;
  if (name == "none") {
    *method = CommonDefs::DOMAIN_FUNCTION_ACTION_NONE;
  } else if (name == "volume") {
    *method = CommonDefs::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_VOLUME;
  } else if (name == "permeability") {
    *method = CommonDefs::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY;
  } else {
    msg << "Flow PK: unknown source distribution method has been specified.";
    Exceptions::amanzi_throw(msg);
  }
}

}  // namespace Flow
}  // namespace Amanzi
