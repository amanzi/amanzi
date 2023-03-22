/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! A factory for creating control-volume logical meshes from segments.
/*!

Builds up a logical mesh from a series of segments.  Note that logical meshes
include cell centroids, but only to make plotting easier -- no geometric
information from these are used.

Logical meshes need the following information:

- cell volume
- face area (and an orientation/normal for problems with gravity)
- length of the cell-to-face connector

This factory looks to build a logical mesh from a series of 1D segments which
may connect to other segments.  From this concept, the needed information can
be provided in a few ways, or inferred.  The primary distinction is between
cell volumes that are provided in the input and cell volumes that are inferred
by the equation:

V_c = sum_{each face in c} length of cell-to-face connector * face area.

``[mesh-logical-spec]``

* `"cell volumes provided`" ``[bool]`` **false** If false, calculate via the
    above formula.  Otherwise each segment must provide the `"cell volume`" or
    `"cell volumes`" parameter described below.

* `"segments`" ``[mesh-logical-segment-spec-list]``  List of segment specs below:

Each segment consists of a collection of cells, and it is assumed that
face-to-cell connectors within the segment are half of the cell length.  For
each segment, cell volumes, face areas, an orientation, and the length of
cell-to-face connectors must be specified:

Cell volumes are either determined by the above formula if the `"cell volumes
provided`" option is specified, or through the following option:

``[mesh-logical-segment-spec]``

ONE OF:
* `"cell volume [m^3]`" ``[double]`` uniform volume of all cells
OR:
* `"cell volumes [m^3]`" ``[Array(double)]`` list of volumes
END

Cell lengths may be specified either as a list of cell lengths or as a segment
length and a number of cells (describing uniform cell length).

ONE OF:
* `"cell lengths [m]`" ``[Array(double)]`` List of cell lengths.
OR
* `"segment length [m]`" ``[double]`` total length, [m]
* `"number of cells`" ``[int]`` number of cells
OR
* `"segment begin point`" ``[Array(double)]`` The segment start point.
* `"segment end point`" ``[Array(double)]`` The segment end point.
* `"number of cells`" ``[int]`` number of cells
END

Note the order of checking is as above -- if both begin/end points and segment
lengths are provided, segment lengths will be used as specified, while the
begin/end points will only be used for visualization.

Face areas are specified through:

ONE OF:
* `"face area [m^2]`" ``[double]`` Face area of all faces in the segment.
OR
* `"face areas [m^2]`" ``[Array(double)]`` List of face areas.
END

Note that the number of faces, and hence the length of the list of face areas,
depends upon how exactly tips are to be handled.  See below.

Orientation of the faces is given by:

* `"orientation`" ``[Array(double)]`` Vector of length dimension describing
    the face normals.  This could be exteneded to be an Array of vectors, but
    currently that is not supported.  Note this is only used if gravity is
    involved.

If segment begin/end points are provided, the orientation is given by
end-begin.

Finally, often tips of segments may be boundaries (faces with boundary
conditions), no type (i.e. there is no face), junction (at least one, likely
more, faces will be connected to the end cell eventually), or branch (the tip
cell will add a face to a junction tip of another segment).

* `"first tip type`" ``[string]`` Type of the first tip.  One of `"none`",
   `"boundary`", `"junction`", or `"branch`".
* `"last tip type`" ``[string]`` Type of the last tip, see above.

If either first or (respectively) last is a branch, then the following should
also be provided:

* `"first tip branch segment`" ``[string]`` Name of the segment to connect to.
* `"first tip branch segment tip`" ``[string]`` One of `"first`" or `"last`",
    the tip to connect this branch face to.

Note that segments with a branch tip MUST appear AFTER the junction tip to
which they will be connected.

This information then lets us know how many face areas to expect.  The face
area array must be at least as long as the number of cells - 1 (interior
faces).  On `"none`" and `"junction`" tips, face areas are not included, as
that face either doesn't exist (for `"none`") or will have its area specified
by any corresponding `"branch`" segments.  On `"boundary`" or `"branch`" tips,
face areas of the end face must be is included.

Finally, a note on centroids.  Centroids can be calculated by providing a
`"segment begin point`" and `"segment end point`".  Note that there is no
expectation that these coordinates differ by the length of the segment.
Instead, cell centroids are located proportionally to the full segment legnth.
This allows meanding streams, for instance, to be visualized as straight lines
while still have the correct solution.

Alternatively, cell centroids can be manually specified (note this is checked
first):

* `"cell centroids`" ``[Array(double)]`` **optional** A list of dimension * n_cells length,
    where the first d values specify cell 0's centroid, the next d specify cell
    1's centroid, etc.

Finally, there are a few factory-global options:

Cell centroids can be inferred by specifying the FIRST segment's first tip,
and orientations and lengths for each segment that do not use first/last tip
coordinates.  Centroids can be spaced out, assuming the first tip of each
branch is a branch tip and the junction it connects to already has a last tip
coordinate.  In this case, specify:

* `"infer cell centroids`" ``[bool]`` **false**.  When
    true, uses the above algorithm.

This option is provided in the `"logical from segments parameters`"
list.

As an example, consider the confluence of the Allegheny and the Monongahela,
which come together in Pittsburgh to form the Ohio.  This could be simply
modeled by three segments:

.. code-block:: xml

  <ParameterList name="mesh">
    <ParameterList name="my mesh">
      <Parameter name="mesh type" type="string" value="logical from segments" />
      <ParameterList name="logical from segments parameters">
        <Parameter name="cell volumes provided" type="bool" value="false" />
        <Parameter name="infer segment begin and end coordaintes" type="bool" value="true" />
        <ParameterList name="segments">
          <ParameterList name="ohio">
            <Parameter name="number of cells" type="int" value="10" />
            <Parameter name="segment length [m]" type="double" value="10000." />
            <Parameter name="cross sectional area [m^2]" type="double" value="300.0" />
            <Parameter name="first tip" type="Array(double)" value="{0.0,0.0}" />
            <Parameter name="first tip type" type="string" value="boundary" />
            <Parameter name="last tip type" type="string" value="junction" />
          </ParameterList>
          <ParameterList name="allegheny">
            <Parameter name="number of cells" type="int" value="10" />
            <Parameter name="segment length [m]" type="double" value="7000." />
            <Parameter name="cross sectional area [m^2]" type="double" value="120.0" />
            <Parameter name="orientation" type="Array(double)" value="{1,1}" />
            <Parameter name="first tip type" type="string" value="branch" />
            <Parameter name="first tip branch segment" type="string" value="ohio" />
            <Parameter name="first tip branch segment tip" type="string" value="last" />
            <Parameter name="last tip type" type="string" value="boundary" />
          </ParameterList>
          <ParameterList name="monongahela">
            <Parameter name="number of cells" type="int" value="10" />
            <Parameter name="segment length [m]" type="double" value="8000." />
            <Parameter name="cross sectional area [m^2]" type="double" value="180.0" />
            <Parameter name="orientation" type="Array(double)" value="{-1,1}" />
            <Parameter name="first tip type" type="string" value="branch" />
            <Parameter name="first tip branch segment" type="string" value="ohio" />
            <Parameter name="first tip branch segment tip" type="string" value="last" />
            <Parameter name="last tip type" type="string" value="boundary" />
          </ParameterList>
        </ParameterList>
      </ParameterList>
    </ParameterList>
  </ParameterList>

 */

#ifndef AMANZI_LOGICAL_MESH_FACTORY2_HH_
#define AMANZI_LOGICAL_MESH_FACTORY2_HH_

#include "Teuchos_ParameterList.hpp"
#include "Epetra_Map.h"
#include "AmanziComm.hh"

#include "MeshLogical.hh"
#include "GeometricModel.hh"

#include "VerboseObject.hh"
#include "dbc.hh"
#include "errors.hh"

namespace Amanzi {
namespace AmanziMesh {

class MeshLogicalFactory {
 public:

  enum class LogicalTip_t {
    NONE,     // do nothing
    BOUNDARY, // tip is a boundary, include the face
    JUNCTION, // tip is a junction -- face will be added later
    BRANCH    // tip branches from a junction.  Add the face.
  };

  
  MeshLogicalFactory(const Comm_ptr_type& comm,
                     const Teuchos::RCP<AmanziGeometry::GeometricModel>& gm,
                     bool calculated_volume = true) :
    comm_(comm),
    gm_(gm),
    calculated_volume_(calculated_volume),
    tracking_centroids_(false)    
  {}

  // Create, assuming AddSegment, etc have all been called previously.
  Teuchos::RCP<MeshLogical>
  Create() const;

  // Create from parameterlist
  Teuchos::RCP<MeshLogical>
  Create(Teuchos::ParameterList& plist);

  //
  // Convenience function to add a segment, mostly for testing
  //
  void
  AddSegment(
      int n_cells,
      const AmanziGeometry::Point& begin,
      const AmanziGeometry::Point& end,
      double face_area,
      MeshLogicalFactory::LogicalTip_t first_tip_type,
      MeshLogicalFactory::LogicalTip_t last_tip_type,
      std::string const& name,
      Entity_ID_List *const cells,
      Entity_ID_List *const faces);
  

  // Add a segment
  //
  // Centroids and cell volumes are optional arguments.
  //
  // Cells and faces are optional return values containing the list of
  // entities in the new segment.
  void
  AddSegment(
      Point_List const * const cell_centroids,
      Double_List const *const cell_volumes,
      Double_List const& cell_lengths,
      Double_List const& face_areas,
      AmanziGeometry::Point const& orientation,
      MeshLogicalFactory::LogicalTip_t first_tip_type,
      MeshLogicalFactory::LogicalTip_t last_tip_type,
      std::string const& name,
      Entity_ID_List *const cells,
      Entity_ID_List *const faces);

  // Add segment from sublist
  void
  AddSegment(Teuchos::ParameterList& sublist);
  
  // Reserve a slot for a face (likely to be added via AddFace!)
  int
  ReserveFace();
  // Manually add a connection, returning the face id.
  int
  AddFace(int f,
          const Entity_ID_List& cells,
          const Point_List& bisectors,
          const std::vector<int>& dirs,
          double area);

  // Add a set from a list of entities
  int
  AddSet(const std::string& set_name,
         const std::string& ent,
         const Entity_ID_List& ents);
  
 protected:
  // Reads a tip type from plist
  MeshLogicalFactory::LogicalTip_t
  GetTipType_(Teuchos::ParameterList& plist, const std::string& pname);

  // Reads a point from plist
  AmanziGeometry::Point
  GetPoint_(Teuchos::ParameterList& plist, const std::string& pname);
  
 protected:
  Double_List cell_volumes_;
  Double_List cell_lengths_;
  Double_List face_areas_;
  std::vector<Entity_ID_List> face_cell_list_;
  std::vector<Point_List> face_cell_bisectors_;
  std::vector<std::vector<int> > face_cell_dirs_;

  Point_List cell_centroids_;

  Comm_ptr_type comm_;
  Teuchos::RCP<AmanziGeometry::GeometricModel> gm_;

  std::map<std::string, Entity_ID_List> seg_cells_;
  std::map<std::string, Entity_ID_List> seg_faces_;
  std::map<std::string, AmanziGeometry::Point> seg_orientations_;

  bool calculated_volume_;

  bool tracking_centroids_;
  std::map<std::string, AmanziGeometry::Point> tracking_end_points_;
};
  

} // namespace AmanziMesh
} // namespace Amanzi


#endif
