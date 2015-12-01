/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "EnumeratedSetRegion.hh"
#include "MeshLogicalFactory.hh"

namespace Amanzi {
namespace AmanziMesh {

// Add a vertical or horizontal segment of uniform size.
void
MeshLogicalFactory::AddSegment(int n_cells,
			       double length,
			       bool vertical,
			       double cross_section_area,
			       bool include_first_face,
			       bool include_last_face,
			       const std::string& set_name,
			       std::vector<Entity_ID>* cells,
			       std::vector<Entity_ID>* faces) {
  if (faces) faces->clear();
  if (cells) cells->clear();

  int cell_first = cell_volume_.size();
  int cell_last = cell_first + n_cells;
  double ds = length / n_cells;
  std::vector<double> ds_halflengths(2, ds/2.0);
  std::vector<double> ds_halflengths_boundary(1, ds/2.0);

  AmanziGeometry::Point normal(0.,0.,0.);
  if (vertical) {
    normal[2] = cross_section_area;
  } else {
    normal[0] = cross_section_area;
  }

  // set cell ids
  std::vector<Entity_ID> new_cells(n_cells);
  for (int i=0; i!=n_cells; ++i) {
    new_cells[i] = cell_first + i;
  }
  if (cells) *cells = new_cells;
  
  // extend cell volumes
  std::vector<double> vols(n_cells, ds * cross_section_area);
  cell_volume_.insert(cell_volume_.end(), vols.begin(), vols.end());

  // extend face_cell_list
  // - first face
  if (include_first_face) {
    if (faces) faces->push_back(face_cell_list_.size());
    Entity_ID_List my_cells(1,cell_first);
    face_cell_list_.push_back(my_cells);
    face_cell_lengths_.push_back(ds_halflengths_boundary);
    face_cell_normals_.push_back(normal);
  }

  // - internal faces
  for (int i=cell_first; i!=cell_last-1; ++i) {
    if (faces) faces->push_back(face_cell_list_.size());
    Entity_ID_List my_cells(2);
    my_cells[0] = i; my_cells[1] = i+1;
    face_cell_list_.push_back(my_cells);
    face_cell_lengths_.push_back(ds_halflengths);
    face_cell_normals_.push_back(normal);
  }

  // -last face
  if (include_last_face) {
    if (faces) faces->push_back(face_cell_list_.size());
    Entity_ID_List my_cells(1,cell_last-1);
    face_cell_list_.push_back(my_cells);
    face_cell_lengths_.push_back(ds_halflengths_boundary);
    face_cell_normals_.push_back(normal);
  }

  // create the region
  // - these are destroyed when the gm is destroyed
  AmanziGeometry::EnumeratedSetRegionPtr enum_rgn =
    new AmanziGeometry::EnumeratedSetRegion(set_name, gm_->Num_Regions(),
			    "CELL", new_cells);
  gm_->Add_Region(enum_rgn);

  centroids_good_ = false;
  return;
}


// Add a segment of uniform size.
void
MeshLogicalFactory::AddSegment(int n_cells,
			       const AmanziGeometry::Point& begin,
			       const AmanziGeometry::Point& end,
			       double cross_section_area,
			       bool include_first_face,
			       bool include_last_face,
			       const std::string& set_name,
			       std::vector<Entity_ID>* cells,
			       std::vector<Entity_ID>* faces) {
  if (faces) faces->clear();
  if (cells) cells->clear();

  int cell_first = cell_volume_.size();
  int cell_last = cell_first + n_cells;

  double ds = AmanziGeometry::norm(end-begin) / n_cells;
  std::vector<double> ds_halflengths(2, ds/2.0);
  std::vector<double> ds_halflengths_boundary(1, ds/2.0);

  AmanziGeometry::Point normal(end-begin);
  normal *= cross_section_area / AmanziGeometry::norm(normal);  

  // set cell ids
  std::vector<Entity_ID> new_cells(n_cells);
  for (int i=0; i!=n_cells; ++i) {
    new_cells[i] = cell_first + i;
  }
  if (cells) *cells = new_cells;
  
  // extend cell volumes
  std::vector<double> vols(n_cells, ds * cross_section_area);
  cell_volume_.insert(cell_volume_.end(), vols.begin(), vols.end());

  // extend face_cell_list
  // - first face
  if (include_first_face) {
    if (faces) faces->push_back(face_cell_list_.size());
    Entity_ID_List my_cells(1,cell_first);
    face_cell_list_.push_back(my_cells);
    face_cell_lengths_.push_back(ds_halflengths_boundary);
    face_cell_normals_.push_back(-normal); // negate as it must be outward normal
  }

  // - internal faces
  for (int i=cell_first; i!=cell_last-1; ++i) {
    if (faces) faces->push_back(face_cell_list_.size());
    Entity_ID_List my_cells(2);
    my_cells[0] = i; my_cells[1] = i+1;
    face_cell_list_.push_back(my_cells);
    face_cell_lengths_.push_back(ds_halflengths);
    face_cell_normals_.push_back(normal);
  }

  // -last face
  if (include_last_face) {
    if (faces) faces->push_back(face_cell_list_.size());
    Entity_ID_List my_cells(1,cell_last-1);
    face_cell_list_.push_back(my_cells);
    face_cell_lengths_.push_back(ds_halflengths_boundary);
    face_cell_normals_.push_back(normal);
  }

  // cell centroids
  normal /= cross_section_area;
  for (int i=0; i!=n_cells; ++i) {
    cell_centroids_.push_back(begin + (i+0.5)*ds*normal);
  }

  // create the region
  // - these are destroyed when the gm is destroyed
  AmanziGeometry::EnumeratedSetRegionPtr enum_rgn =
    new AmanziGeometry::EnumeratedSetRegion(set_name, gm_->Num_Regions(),
			    "CELL", new_cells);
  gm_->Add_Region(enum_rgn);

  return;
}


Teuchos::RCP<MeshLogical>
MeshLogicalFactory::Create() {
  Teuchos::RCP<AmanziMesh::MeshLogical> mesh;
  if (centroids_good_) {
    mesh = Teuchos::rcp(new MeshLogical(comm_,
					cell_volume_,
					face_cell_list_,
					face_cell_lengths_,
					face_cell_normals_,
					&cell_centroids_));
    
  } else {
    mesh = Teuchos::rcp(new MeshLogical(comm_,
					cell_volume_,
					face_cell_list_,
					face_cell_lengths_,
					face_cell_normals_));
  }
  mesh->set_geometric_model(gm_);
  return mesh;
}


// Manually add a connection, returning the face id.
int
MeshLogicalFactory::AddConnection(const Entity_ID_List& cells,
				  const AmanziGeometry::Point& normal,
				  const std::vector<double>& lengths,
				  double area) {

  int f = face_cell_list_.size();
  face_cell_list_.push_back(cells);
  face_cell_normals_.push_back(area/AmanziGeometry::norm(normal)*normal);
  face_cell_lengths_.push_back(lengths);
  return f;
}

  
} // namespace AmanziMesh
} // namespace Amanzi
