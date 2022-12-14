/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
      Ethan Coon (ecoon@lanl.gov)
*/

//! RegionLabeledSet: A region defined by a set of mesh entities in a mesh file
/*!

The list *region: labeled set* defines a named set of mesh entities
existing in an input mesh file. This is the same file that contains
the computational mesh. The name of the entity set is given
by *label*.  For example, a mesh file in the Exodus II
format can be processed to tag cells, faces and/or nodes with
specific labels, using a variety of external tools. Regions based
on such sets are assigned a user-defined label for Amanzi, which may
or may not correspond to the original label in the exodus file.
Note that the file used to express this labeled set may be in any
Amanzi-supported mesh format (the mesh format is specified in the
parameters for this option).  The *entity* parameter may be
necessary to specify a unique set.  For example, an Exodus file
requires *cell*, *face* or *node* as well as a label (which is
an integer).  The resulting region will have the dimensionality
associated with the entities in the indicated set.

.. _region-labeled-set-spec:
.. admonition:: region-labeled-set-spec

    * `"label`" ``[string]`` Set per label defined in the mesh file.
    * `"file`" ``[string]`` File name.
    * `"format`" ``[string]`` Currently, we only support mesh files in the "Exodus II" format.
    * `"entity`" ``[string]`` Type of the mesh object (cell, face, etc).

Example:

.. code-block:: xml

   <ParameterList name="AQUIFER">
     <ParameterList name="region: labeled set">
       <Parameter name="entity" type="string" value="cell"/>
       <Parameter name="file" type="string" value="porflow4_4.exo"/>
       <Parameter name="format" type="string" value="Exodus II"/>
       <Parameter name="label" type="string" value="1"/>
     </ParameterList>
   </ParameterList>

*/


/*
  Strictly speaking, we should tie this region class to a particular
  mesh or mesh file but that cause a circular dependency of meshes
  on regions and of labeled set regions on meshes. We will rely on the
  fact that when a mesh is created specifying a geometric model, it
  will create mesh entity sets based on the labeled sets in that
  geometric model.

  If we need to change this behavior, then we can make a forward
  declaration of AmanziMesh::Mesh, make the Mesh class a friend, add
  a mesh variable to this class and have a protected method to set
  the mesh
*/

#ifndef AMANZI_REGION_LABELED_SET_HH_
#define AMANZI_REGION_LABELED_SET_HH_

#include "Region.hh"

namespace Amanzi {
namespace AmanziGeometry {

class RegionLabeledSet : public Region {
 public:
  // constructor
  RegionLabeledSet(const std::string& name,
                   const int id,
                   const std::string& entity_str,
                   const std::string& file,
                   const std::string& format,
                   const std::string& label,
                   const LifeCycleType lifecycle = LifeCycleType::PERMANENT);

  // Label in the file
  const std::string& label() const { return label_; }

  // Is the the specified point inside this region
  bool inside(const Point& p) const;

  const std::string& entity_str() const { return entity_str_; }

 protected:
  const std::string entity_str_; // what kind of entities make up this set
  const std::string file_;       // which file are we supposed to read it from
  const std::string format_;     // format of the file
  const std::string label_; // Label used to identify set in the file (may be different from name)
};

} // namespace AmanziGeometry
} // namespace Amanzi


#endif
