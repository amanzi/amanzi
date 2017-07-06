/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_LOGICAL_MESH_FACTORY_H_
#define AMANZI_LOGICAL_MESH_FACTORY_H_

#include "Teuchos_ParameterList.hpp"
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

  enum LogicalTip_t {
    TIP_NULL,
    TIP_BOUNDARY, // tip is a root boundary, i.e. face.  add cell and face.
    TIP_DEFERRED, // tip included, but face will be added later
    TIP_JUNCTION, // tip is a root junction cell.  add cell.
    TIP_BRANCH // tip branches from a junction.  add neither cell nor face.
  };
  
  MeshLogicalFactory(const Epetra_MpiComm* incomm_,
                     const Teuchos::RCP<AmanziGeometry::GeometricModel>& gm) :
    comm_(incomm_),
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
	     MeshLogicalFactory::LogicalTip_t first_tip,
	     MeshLogicalFactory::LogicalTip_t last_tip,
	     const std::string& set_name,
	     std::vector<Entity_ID>* cells,
	     std::vector<Entity_ID>* faces,
             double* cell_length);

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
	     MeshLogicalFactory::LogicalTip_t first_tip,
	     MeshLogicalFactory::LogicalTip_t last_tip,
	     const std::string& set_name,
	     std::vector<Entity_ID>* cells,
	     std::vector<Entity_ID>* faces,
             double* cell_length);


  // Add a segment totally generically
  void
  AddSegment(const AmanziGeometry::Point& begin,
             const AmanziGeometry::Point& end,
             std::vector<double> lengths,
             std::vector<double> areas,
             std::vector<double> vols,
             MeshLogicalFactory::LogicalTip_t first_tip,
             MeshLogicalFactory::LogicalTip_t last_tip,
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
         const std::string& ent,
         const Entity_ID_List& ents);
  
  Teuchos::RCP<MeshLogical>
  Create();

  Teuchos::RCP<MeshLogical>
  Create(Teuchos::ParameterList& plist);

  
 protected:
  std::vector<double> cell_volume_;
  std::vector<Entity_ID_List> face_cell_list_;
  std::vector<std::vector<double> > face_cell_lengths_;
  std::vector<AmanziGeometry::Point> face_cell_normals_;
  std::vector<AmanziGeometry::Point> cell_centroids_;
  bool centroids_good_;

  const Epetra_MpiComm* comm_;
  Teuchos::RCP<AmanziGeometry::GeometricModel> gm_;
  

};
  

} // namespace AmanziMesh
} // namespace Amanzi


#endif
