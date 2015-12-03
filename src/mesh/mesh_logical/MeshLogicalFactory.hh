/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_LOGICAL_MESH_FACTORY_H_
#define AMANZI_LOGICAL_MESH_FACTORY_H_

#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"

#include "MeshLogical.hh"
#include "GeometricModel.hh"

#include "VerboseObject.hh"
#include "dbc.hh"
#include "errors.hh"

namespace Amanzi {
namespace AmanziMesh {

class MeshLogicalFactory {
 public:
  MeshLogicalFactory(const Epetra_MpiComm* incomm,
		     AmanziGeometry::GeometricModelPtr& gm) :
    comm_(incomm),
    gm_(gm),
    centroids_good_(true)
  {}

  // Add a vertical or horizontal segment of uniform size.
  //
  // Adds a segment of with a specified length and number of cells,
  // subdividing uniformly.  include_last_face and include_first_face
  // set whether to add a boundary face -- if this segment is to be
  // connected to another segment, these boundary faces should be
  // manually added to connect to both cells.
  //
  // Also makes a region of this segment of name set_name.
  void
  AddSegment(int n_cells,
	     double length,
	     bool vertical,
	     double cross_section_area,
	     bool include_first_face,
	     bool include_last_face,
	     const std::string& set_name,
	     std::vector<Entity_ID>* cells,
	     std::vector<Entity_ID>* faces);

  // Add a segment of uniform size.
  //
  // Adds a segment of with a specified length and number of cells,
  // subdividing uniformly.  include_last_face and include_first_face
  // set whether to add a boundary face -- if this segment is to be
  // connected to another segment, these boundary faces should be
  // manually added to connect to both cells.
  //
  // Also makes a region of this segment of name set_name.
  void
  AddSegment(int n_cells,
	     const AmanziGeometry::Point& begin,
	     const AmanziGeometry::Point& end,
	     double cross_section_area,
	     bool include_first_face,
	     bool include_last_face,
	     const std::string& set_name,
	     std::vector<Entity_ID>* cells,
	     std::vector<Entity_ID>* faces);


  // Manually add a connection, returning the face id.
  int
  AddConnection(const Entity_ID_List& cells,
		const AmanziGeometry::Point& normal,
		const std::vector<double>& lengths,
		double area);		

  int
  AddSet(const std::string& set_name,
         const Entity_ID_List* cells,
         const Entity_ID_List* faces);
  
  Teuchos::RCP<MeshLogical>
  Create();

  // Teuchos::RCP<MeshLogical>
  // Create(Teuchos::ParameterList& plist) {}

  
 protected:
  std::vector<double> cell_volume_;
  std::vector<Entity_ID_List> face_cell_list_;
  std::vector<std::vector<double> > face_cell_lengths_;
  std::vector<AmanziGeometry::Point> face_cell_normals_;
  std::vector<AmanziGeometry::Point> cell_centroids_;
  bool centroids_good_;

  const Epetra_MpiComm* comm_;
  AmanziGeometry::GeometricModelPtr gm_;
  

};
  

} // namespace AmanziMesh
} // namespace Amanzi


#endif
