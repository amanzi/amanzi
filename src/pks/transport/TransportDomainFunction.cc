/*
This is the transport component of the Amanzi code. 

License: see $AMANZI_DIR/COPYRIGHT
Author (v1): Neil Carlson
       (v2): Ethan Coon, Markus Berndt, Konstantin Lipnikov
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
