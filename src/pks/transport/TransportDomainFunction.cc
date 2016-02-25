/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson
           Ethan Coon
           Markus Berndt
           Konstantin Lipnikov
*/


#include <algorithm>
#include "errors.hh"

#include "TransportDomainFunction.hh"

namespace Amanzi {
namespace Transport {

/* ******************************************************************
* Calculate pairs <list of cells, function>
****************************************************************** */
void TransportDomainFunction::Define(
    const std::vector<std::string>& regions,
    const Teuchos::RCP<const MultiFunction>& f,
    int action, int submodel, const std::string& name)
{
  PK_DomainFunction::Define(regions, f, action, submodel);
  tcc_name_ = name;
}


/* ******************************************************************
* Calculate pairs <list of cells, function>
****************************************************************** */
void TransportDomainFunction::Define(
    const std::string& region,
    const Teuchos::RCP<const MultiFunction>& f,
    int action, int submodel, const std::string& name)
{
  PK_DomainFunction::Define(region, f, action, submodel);
  tcc_name_ = name;
}

}  // namespace Transport
}  // namespace Amanzi
