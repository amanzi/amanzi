/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
           Ethan Coon (ecoon@lanl.gov)
*/

/*
  Nonmember function for creating regions.

*/

#ifndef AMANZI_REGION_FACTORY_HH_
#define AMANZI_REGION_FACTORY_HH_

#include "AmanziTypes.hh"
#include "Teuchos_ParameterList.hpp"

#include "GeometryDefs.hh"

namespace Amanzi {
namespace AmanziGeometry {

class Region;

Teuchos::RCP<Region>
createRegion(const std::string& reg_name,
             const std::string& reg_type,
             int reg_id,
             Teuchos::ParameterList& reg_spec,
             const Comm_type& comm);

// deprecate this, this old style
inline
Teuchos::RCP<Region>
createRegion(const std::string& reg_name,
             int reg_id,
             Teuchos::ParameterList& plist,
             const Comm_type& comm) {
  std::string reg_type = plist.begin()->first;
  Teuchos::ParameterList& reg_list = plist.sublist(reg_type);
  // strip the "region: "
  reg_type = reg_type.substr(8, std::string::npos);
  return createRegion(reg_name, reg_type, reg_id, reg_list, comm);
}


} // namespace AmanziGeometry
} // namespace Amanzi

#endif
