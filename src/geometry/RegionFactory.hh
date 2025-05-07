/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
           Ethan Coon (ecoon@lanl.gov)
*/

/*!

Regions are geometrical constructs used to define subsets of the computational
domain.  Regions are used for many things -- defining piecewise initial
conditions, boundary conditions, source terms, observations, and more.  Regions
may represents one-, two- or three-dimensional subsets of physical
space.

Together, a region + mesh + entity kind (cell, face, etc) define a "resolved"
region.  This is a discrete list of entities.  Most regions exist outside of
the mesh, and may be resolved on any mesh whose spatial dimension matches the
region's spatial dimension.  The exception to this are `"enumerated`" and
`"labeled set`" regions, which may only be resolved on the mesh for which they
were created, or for meshes derived from that mesh (see below).

<<<INSERT INFO ON THE RESOLUTION PROCESS AND "INSIDE" FOR GEOMETRIC REGIONS>>>

.. _region-typed-spec:
.. admonition:: region-typed-spec

   * `"region type`" ``[string]`` One of:

      - `"all`" An All_ region.
      - `"boundary`"  A Boundary_
      - `"point`"
      - `"line segment`"
      - `"plane`"
      - `"polygon`"
      - `"box`"
      - `"box volume fractions`"
      - `"cylinder`"
      - `"halfspace`"
      - `"labeled set`"
      - `"enumerated set`"
      - `"enumerated set from file`"
      - `"color function`"
      - `"level set`"
      - `"logical`"

   * `"_region_type_ parameters`" ``[_region_type_-spec]`` See below.

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
inline Teuchos::RCP<Region>
createRegion(const std::string& reg_name,
             int reg_id,
             Teuchos::ParameterList& plist,
             const Comm_type& comm)
{
  std::string reg_type = plist.begin()->first;
  Teuchos::ParameterList& reg_list = plist.sublist(reg_type);
  // strip the "region: "
  reg_type = reg_type.substr(8, std::string::npos);
  return createRegion(reg_name, reg_type, reg_id, reg_list, comm);
}


} // namespace AmanziGeometry
} // namespace Amanzi

#endif
