/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <string>
#include <vector>

#include "errors.hh"

#include "TransportSourceFactory.hh"

namespace Amanzi {
namespace Transport {

/* ******************************************************************
* Process source, step 1.
****************************************************************** */
void TransportSourceFactory::Create(std::vector<TransportDomainFunction*>& srcs) 
{
  Errors::Message msg;

  if (plist_->isSublist("concentration")) {
    Teuchos::ParameterList& clist = plist_->get<Teuchos::ParameterList>("concentration");

    // Iterate through the source specification sublists in the clist.
    // All are expected to be sublists of identical structure.
    for (Teuchos::ParameterList::ConstIterator it = clist.begin(); it != clist.end(); ++it) {
      std::string name = it->first;
      if (clist.isSublist(name)) {
        Teuchos::ParameterList& srclist = clist.sublist(name);
	for (Teuchos::ParameterList::ConstIterator it1 = srclist.begin(); it1 != srclist.end(); ++it1) {
	  std::string specname = it1->first;

	  if (srclist.isSublist(specname)) {
	    Teuchos::ParameterList& spec = srclist.sublist(specname);
            try {
              TransportDomainFunction* src = new TransportDomainFunction(mesh_);
              ProcessSourceSpec(spec, name, src);
              srcs.push_back(src);
            } catch (Errors::Message& m) {
              msg << "in sublist \"" << specname.c_str() << "\": " << m.what();
              Exceptions::amanzi_throw(msg);
            }
          } else {
            msg << "parameter \"" << specname.c_str() << "\" is not a sublist";
            Exceptions::amanzi_throw(msg);
          }
        }
      } else {
        msg << "parameter \"" << name.c_str() << "\" is not a sublist";
        Exceptions::amanzi_throw(msg);
      }
    }
  } else {
    msg << "Transport PK: \"source terms\" has no sublist \"concentration\".\n";
    Exceptions::amanzi_throw(msg);  
  }
}


/* ******************************************************************
* Process source, step 2.
****************************************************************** */
void TransportSourceFactory::ProcessSourceSpec(
  Teuchos::ParameterList& list, const std::string& name, TransportDomainFunction* src) const 
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
  Teuchos::RCP<MultiFunction> f;
  try {
    f = Teuchos::rcp(new Amanzi::MultiFunction(*f_list));
  } catch (Errors::Message& msg) {
    m << "error in source sublist : " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  // Add this source specification to the domain function.
  int method;
  std::string action_name = list.get<std::string>("spatial distribution method", "none");
  ProcessStringActions(action_name, &method);

  std::string submodel_name = list.get<std::string>("submodel", "rate");
  int submodel(CommonDefs::DOMAIN_FUNCTION_SUBMODEL_RATE);
  if (submodel_name == "integrated source") submodel = CommonDefs::DOMAIN_FUNCTION_SUBMODEL_INTEGRAL;

  src->Define(regions, f, method, submodel, name);
}


/* ****************************************************************
* Process string for a source specipic action.
**************************************************************** */
void TransportSourceFactory::ProcessStringActions(const std::string& name, int* method) const 
{
  Errors::Message msg;
  if (name == "none") {
    *method = CommonDefs::DOMAIN_FUNCTION_ACTION_NONE;
  } else if (name == "volume") {
    *method = CommonDefs::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_VOLUME;
  } else if (name == "permeability") {
    *method = CommonDefs::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY;
  } else {
    msg << "Transport PK: unknown source distribution method has been specified.";
    Exceptions::amanzi_throw(msg);
  }
}

}  // namespace Transport
}  // namespace Amanzi
