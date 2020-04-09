/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      William Perkins
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

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
  for (auto i = gm_params.begin(); i != gm_params.end(); ++i) {
    std::string region_name = gm_params.name(i);
    if (gm_params.isSublist(region_name)) {
      // Extract sublist specifying region
      Teuchos::ParameterList& reg_spec = gm_params.sublist(region_name);

      // Create the region
      Teuchos::RCP<Region> reg = createRegion(region_name, -1, reg_spec, comm);

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
  if (dim_ < reg->space_dimension()) {
    Errors::Message mesg;
    mesg << "Dimension of geometric model less than that of region \""
         << reg->name() << "\" space dimension";
    Exceptions::amanzi_throw(mesg);
  }

  if (dim_ < reg->manifold_dimension()) {
    Errors::Message mesg;
    mesg << "Dimension of geometric model less than that of region \""
         << reg->name() << "\" topological dimension";
    Exceptions::amanzi_throw(mesg);
  }

  reg->set_id(RegionSize());
  regions_.push_back(reg);
  regions_name_[reg->name()] = reg;
  regions_id_[reg->id()] = reg;
}


} // namespace AmanziGeometry
} // namespace Amanzi
