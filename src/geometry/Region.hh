/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      William Perkins
      Ethan Coon (coonet@ornl.gov)
*/

//! Region: a geometric or discrete subdomain (abstract)

/*!

Regions are geometrical constructs used in Amanzi to define subsets of
the computational domain in order to specify the problem to be solved, and the
output desired. Regions may represents zero-, one-, two- or three-dimensional
subsets of physical space.  for a three-dimensional problem, the simulation
domain will be a three-dimensional region bounded by a set of two-dimensional
regions.  If the simulation domain is N-dimensional, the boundary conditions
must be specified over a set of regions are (N-1)-dimensional.

Amanzi automatically defines the special region labeled *All*, which is the
entire simulation domain. Currently, the unstructured framework does
not support the *All* region, but it is expected to do so in the
near future.

Amanzi supports parameterized forms for a number of analytic shapes, as well
as more complex definitions based on triangulated surface files.


ONE OF:
* `"region: box`" ``[region-box-spec]``
OR:
* `"region: plane`" ``[region-plane-spec]``
OR:
* `"region: labeled set`" ``[region-labeled-set-spec]``
OR:
* `"region: color function`" ``[region-color-function-spec]``
OR:
* `"region: point`" ``[region-point-spec]``
OR:
* `"region: logical`" ``[region-logical-spec]``
OR:
* `"region: polygon`" ``[region-polygon-spec]``
OR:
* `"region: enumerated`" ``[region-enumerated-spec]``
OR:
* `"region: all`" ``[list]``
OR:
* `"region: boundary`" ``[region-boundary-spec]``
OR:
* `"region: box volume fractions`" ``[region-box-volume-fractions-spec]``
OR:
* `"region: line segment`" ``[region-line-segment-spec]``
END


Notes:

- Surface files contain labeled triangulated face sets.  The user is
  responsible for ensuring that the intersections with other surfaces
  in the problem, including the boundaries, are *exact* (*i.e.* that
  surface intersections are *watertight* where applicable), and that
  the surfaces are contained within the computational domain.  If
  nodes in the surface fall outside the domain, the elements they
  define are ignored.

  Examples of surface files are given in the *Exodus II* file
  format here.

- Region names must NOT be repeated.

Example:

.. code-block:: xml

   <ParameterList>  <!-- parent list -->
     <ParameterList name="regions">
       <ParameterList name="TOP SECTION">
         <ParameterList name="region: box">
           <Parameter name="low coordinate" type="Array(double)" value="{2, 3,
5}"/> <Parameter name="high coordinate" type="Array(double)" value="{4, 5, 8}"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="MIDDLE SECTION">
         <ParameterList name="region: box">
           <Parameter name="low coordinate" type="Array(double)" value="{2, 3,
3}"/> <Parameter name="high coordinate" type="Array(double)" value="{4, 5, 5}"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="BOTTOM SECTION">
         <ParameterList name="region: box">
           <Parameter name="low coordinate" type="Array(double)" value="{2, 3,
0}"/> <Parameter name="high coordinate" type="Array(double)" value="{4, 5, 3}"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="INFLOW SURFACE">
         <ParameterList name="region: labeled set">
           <Parameter name="label"  type="string" value="sideset_2"/>
           <Parameter name="file"   type="string" value="F_area_mesh.exo"/>
           <Parameter name="format" type="string" value="Exodus II"/>
           <Parameter name="entity" type="string" value="face"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="OUTFLOW PLANE">
         <ParameterList name="region: plane">
           <Parameter name="point" type="Array(double)" value="{0.5, 0.5,
0.5}"/> <Parameter name="normal" type="Array(double)" value="{0, 0, 1}"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="BLOODY SAND">
         <ParameterList name="region: color function">
           <Parameter name="file" type="string" value="F_area_col.txt"/>
           <Parameter name="value" type="int" value="25"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="FLUX PLANE">
         <ParameterList name="region: polygon">
           <Parameter name="number of points" type="int" value="5"/>
           <Parameter name="points" type="Array(double)" value="{-0.5, -0.5,
-0.5, 0.5, -0.5, -0.5, 0.8, 0.0, 0.0, 0.5,  0.5, 0.5, -0.5, 0.5, 0.5}"/>
          </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>

In this example, *TOP SESCTION*, *MIDDLE SECTION* and *BOTTOM SECTION*
are three box-shaped volumetric regions. *INFLOW SURFACE* is a
surface region defined in an Exodus II-formatted labeled set
file and *OUTFLOW PLANE* is a planar region. *BLOODY SAND* is a volumetric
region defined by the value 25 in color function file.

*/

#ifndef AMANZI_REGION_HH_
#define AMANZI_REGION_HH_

#include <string>

#include "VerboseObject.hh"
#include "GeometryDefs.hh"

namespace Amanzi {
namespace AmanziGeometry {

class Point;

class Region {
 public:
  // Destructor
  virtual ~Region(){};

  // Dimension of the subdomain
  unsigned int manifold_dimension() const { return manifold_dimension_; }
  void set_manifold_dimension(unsigned int dimension)
  {
    manifold_dimension_ = dimension;
  }

  // Dimension of points in the subdomain
  unsigned int space_dimension() const { return space_dimension_; }
  void set_space_dimension(unsigned int dimension)
  {
    space_dimension_ = dimension;
  }

  // Name of the region -- no setter (set by constructor)
  std::string name() const { return name_; }

  // Integer identifier of the region
  int id() const { return id_; }
  void set_id(int id) { id_ = id; }

  // Geometric/enumerated
  bool is_geometric() const { return geometric_; }

  // Type of the region
  RegionType type() const { return type_; }

  // Get the Lifecycle of this region - Do mesh entity sets derived from
  // it have to be kept around or are they temporary and can be destroyed
  // as soon as they are used?
  LifeCycleType lifecycle() const { return lifecycle_; }

  // Is the specified point inside the closure of the Region
  virtual bool inside(const Point& p) const = 0;

  // Calculate intersection measure of object with the Region. The intersection
  // is defined when the object and Region have same dimensionality.
  //
  // -- counter clockwise ordered polygon does not require faces
  double intersect(const Kokkos::View<Point*>& polytope) const
  {
    std::vector<std::vector<int>> faces;
    return intersect(polytope, faces);
  }

  // Polyhedron with counter clockwise ordered faces (wrt normals)
  virtual double intersect(const Kokkos::View<Point*>& polytope,
                           const std::vector<std::vector<int>>& faces) const
  {
    return -1.0;
  }

 protected:
  // Constructor -- protected as it should never be called directly
  Region(const std::string& name, int id, bool geometric, RegionType type,
         unsigned int dim, unsigned int geom_dim,
         LifeCycleType lifecycle = PERMANENT)
    : name_(name),
      id_(id),
      geometric_(geometric),
      type_(type),
      manifold_dimension_(dim),
      space_dimension_(geom_dim),
      lifecycle_(lifecycle){};

 protected:
  // Lifecycle (Temporary or Permanent)
  LifeCycleType lifecycle_;

  // Topological dimension of region (0, 1, 2, 3)
  unsigned int manifold_dimension_;
  unsigned int space_dimension_;

  // Name of identifier
  std::string name_;

  // Identifier of region
  int id_;

  // Region type
  RegionType type_;

  // Geometric or enumerated?
  bool geometric_;

 private:
  Region(const Region& other);            // prevent copy constructor
  Region& operator=(const Region& other); // prevent operator=
};


typedef std::vector<Teuchos::RCP<Region>> RegionVector;

} // namespace AmanziGeometry
} // namespace Amanzi

#endif
