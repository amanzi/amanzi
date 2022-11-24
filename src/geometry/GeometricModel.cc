/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
  Collection of Regions which decompose the domain into subdomains.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: William Perkins
           Ethan Coon (ecoon@lanl.gov)
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
    const std::string& region_name = reg_item.first;
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
  // NOTE: extracted regions may break this -- for instance a 3D region on a 2D
  // surface mesh is the restriction of the parent entities of the 3D region on
  // the parent mesh.
  //
  // if (dim_ < reg->get_space_dimension()) {
  //   Errors::Message mesg;
  //   mesg << "Dimension of geometric model less than that of region \""
  //        << reg->get_name() << "\" space dimension";
  //   Exceptions::amanzi_throw(mesg);
  // }

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
