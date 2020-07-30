/*
  Input Converter

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <algorithm>
#include <sstream>
#include <string>

//TPLs
#include <xercesc/dom/DOM.hpp>

// Amanzi's
#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"

#include "InputConverterU.hh"
#include "InputConverterU_Defs.hh"

namespace Amanzi {
namespace AmanziInput {

XERCES_CPP_NAMESPACE_USE

/* ******************************************************************
* Create energy list.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateShallowWater_(const std::string& domain)
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating shallow water, domain=" << domain << std::endl;

  MemoryManager mm;
  DOMNode* node;

  out_list.set<std::string>("domain name", (domain == "matrix") ? "domain" : domain);

  out_list.sublist("reconstruction")
      .set<int>("polynomial order", 0)
      .set<std::string>("limiter", "Barth-Jespersen")
      .set<std::string>("limiter stencil", "cell to closest cells")
      .set<std::string>("limiter location", "node")
      .set<double>("limiter cfl", 0.1);

  // boundary conditions
  out_list.sublist("boundary conditions");

  out_list.sublist("verbose object") = verb_list_.sublist("verbose object");

  return out_list;
}

}  // namespace AmanziInpit
}  // namespace Amanzi
