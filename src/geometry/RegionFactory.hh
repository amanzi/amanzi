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

} // namespace AmanziGeometry
} // namespace Amanzi

#endif
