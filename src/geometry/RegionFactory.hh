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
region.  This is a discrete list of entities.

Most regions exist outside of the mesh, and may be resolved on any mesh whose
spatial dimension matches the region's spatial dimension.  For these
**geometric** regions, the resolution process is based on **centroids**.  An
entity on a given mesh of a given type is "inside" the region if the centroid
of that entity is within a tolerance (by default, 1.e-8, but this may be
changed using the `"expert parameters->tolerance`" parameter) of the region.
So, a cell, face, or edge need not be entirely within the region; only its
centroid, and if only a small portion of the cell is inside the region, it may
not be included in the set.

The exception to this rule is for `"point`", `"line segment`", and `"box volume
fractions`" regions.  In the case of resolving CELL entities against these
regions, we try to find the intersection of the region and the entity, rather
than simply checking centroids.  For `"line segment`" and `"box volume
fractions`" regions, the intersection is checked properly for FACE entities as
well (this should probably be done for `"point`" regions, but currently is
not).

Note that this may result in multiple cells being in the `"point`" region.  For
instance, if the point is within tolerance of a node, then all cells adjacent
to that node will be in the region.

Unlike geometric regions, **discrete** regions are regions that consider a
discrete set of entities.  `"labeled set,`" `"enumerated set`", and `"color
function`" regions are ways of identifying specific discrete entities, and may
be resolved on the mesh upon which they were labeled.

Special consideration is taken for resolving discrete regions on extracted or
derived meshes.  For instance, one may label FACEs of a 3D volumetric mesh,
then extract a surface mesh, then resolve the labeled set region against CELL
entities of the extracted mesh.  In this case, if surface cell's corresponding
face is labeled in the volumetric mesh, it is in the region.  Similar
strategies are made to resolve entities on all extracted meshes; wherever
possible, we try to do the logical thing to provide as much flexibility to the
user as possible.  Note that this is not done for geometric regions.  So, for
instance, rather than resolving a 3D box region on a 2D extracted surface mesh,
one must define a flattened 2D box region instead by dropping the z coordinate
of the 3D box.

Lastly, `"logical`" regions define logical relationships between other regions.


.. _region-typed-spec:
.. admonition:: region-typed-spec

   * `"region type`" ``[string]`` One of:

      - `"all`" See All_
      - `"boundary`"  See Boundary_
      - `"point`" See Point_
      - `"line segment`" See `Line Segment`_
      - `"plane`" See Plane_
      - `"polygon`" See Polygon_
      - `"box`" See Box_
      - `"box volume fractions`" See `Box Volume Fractions`_
      - `"cylinder`" See Cylinder_
      - `"halfspace`" See `Half Space`_
      - `"level set`" See `Level Set`_.
      - `"labeled set`" ` See `Labeled Set`_
      - `"enumerated set`" See `Enumerated Set`_
      - `"enumerated set from file`" See `Enumerated Set`_
      - `"color function`" See `Color Function`_.
      - `"logical`" See Logical_.

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

Teuchos::RCP<Region> createRegion(const std::string& reg_name,
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
