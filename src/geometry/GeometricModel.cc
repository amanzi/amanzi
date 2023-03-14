/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: William Perkins
           Ethan Coon (ecoon@lanl.gov)
*/

/*
  Collection of Regions which decompose the domain into subdomains.

*/

#include "dbc.hh"
#include "errors.hh"

#include "RegionFactory.hh"
#include "GeometricModel.hh"

namespace Amanzi {
namespace AmanziGeometry {

// Constructor with no regions
GeometricModel::GeometricModel(unsigned int dim) : dim_(dim) {}

GeometricModel::GeometricModel(unsigned int dim,
                               Teuchos::ParameterList& gm_params,
                               const Comm_type& comm)
  : dim_(dim)
{
  // Go through the parameter list and populate the geometric model with regions
  for (const auto& reg_item : gm_params) {
    const std::string& reg_name = reg_item.first;
    if (gm_params.isSublist(reg_name)) {
      // Extract sublist specifying region
      Teuchos::ParameterList& reg_spec = gm_params.sublist(reg_name);

      // Create the region
      // Two types of specs -- legacy format is a sublist, new format is a typeed list
      Teuchos::RCP<Region> reg;
      if (reg_spec.isParameter("region type")) {
        // new style, with type parameter
        std::string region_type = reg_spec.get<std::string>("region type");
        reg = createRegion(reg_name, region_type, -1, reg_spec, comm);
      } else {
        // deprecate this -- old style e.g. "region: box" sublist, see #181
        if (reg_spec.numParams() != 1) {
          Errors::Message msg;
          msg << "Region spec \"" << reg_name << "\" should have exactly one shape sublist.";
          Exceptions::amanzi_throw(msg);
        }
        reg = createRegion(reg_name, -1, reg_spec, comm);
      }

      // Add it to the geometric model
      AddRegion(reg);

    } else {
      Errors::Message mesg("Error: Improper region specification");
      Exceptions::amanzi_throw(mesg);
    }
  }
}


// Add a Region
void
GeometricModel::AddRegion(const Teuchos::RCP<Region>& reg)
{
  if (dim_ < reg->get_manifold_dimension()) {
    Errors::Message mesg;
    mesg << "Dimension of geometric model less than that of region \"" << reg->get_name()
         << "\" topological dimension";
    Exceptions::amanzi_throw(mesg);
  }

  reg->set_id(size());

  Teuchos::RCP<const Region> rc = reg;
  regions_name_[reg->get_name()] = rc;


  Teuchos::RCP<const Region> rc2 = reg;
  regions_id_[reg->get_id()] = rc2;

  regions_.push_back(reg);
}


} // namespace AmanziGeometry
} // namespace Amanzi
