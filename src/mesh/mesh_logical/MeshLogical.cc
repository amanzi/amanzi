/*
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!
 Logical mesh that can be modified and constructed on the fly.

 Logical mesh is a topologically defined mesh with no real coordinate
 geometry.  By definition it is perfectly parallel with no ghost entities,
 as it is intended to be used along with a normal mesh as a subgrid model.
 As it is not a geomtric mesh, it cannot work with all (many) spatial
 discretizations -- currently only Finite Volume.

 In particular:
  1. nodes do not exist
*/

#include <set>
#include "RegionEnumerated.hh"
#include "MeshLogical.hh"

namespace Amanzi {
namespace AmanziMesh {

// lumped things for more efficient calculation
std::pair<double, AmanziGeometry::Point>
MeshLogicalAlgorithms::computeCellGeometry(const MeshFramework& mesh, const Entity_ID c) const
{
  double volume = mesh.getCellVolume(c);
  AmanziGeometry::Point centroid = mesh.getCellCentroid(c);
  return std::make_pair(volume, centroid);
}

std::tuple<double, AmanziGeometry::Point, cPoint_View>
MeshLogicalAlgorithms::computeFaceGeometry(const MeshFramework& mesh, const Entity_ID f) const
{
  double area = mesh.getFaceArea(f);
  AmanziGeometry::Point centroid = mesh.getFaceCentroid(f);

  cEntity_ID_View fcells;
  mesh.getFaceCells(f, Parallel_kind::ALL, fcells);
  Point_View normals("normals", fcells.size());
  for (int i = 0; i != fcells.size(); ++i) {
    normals[i] = mesh.getFaceNormal(f, fcells[i], nullptr);
  }
  return std::make_tuple(area, centroid, normals);
}

std::pair<AmanziGeometry::Point, AmanziGeometry::Point>
MeshLogicalAlgorithms::computeEdgeGeometry(const MeshFramework& mesh, const Entity_ID e) const
{
  Errors::Message msg("There are no edges in a MeshLogical.");
  Exceptions::amanzi_throw(msg);
  return std::make_pair(AmanziGeometry::Point(), AmanziGeometry::Point());
}

//
// Topological constructor of a MeshLogical splits topology from geometry.
//
MeshLogical::MeshLogical(const Comm_ptr_type& comm,
                         const std::vector<Entity_ID_List>& face_cell_ids,
                         const std::vector<std::vector<int>>& face_cell_dirs,
                         const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : MeshFramework(comm, Teuchos::null, plist), face_cell_ids_(face_cell_ids)
{
  setSpaceDimension(3);
  setManifoldDimension(1);
  setAlgorithms(Teuchos::rcp(new MeshLogicalAlgorithms()));

  // Count number of cells referenced, and check that the number of cells
  // referenced is equal to the largest id referenced (+1 for 0)
  int c_max = -1;
  std::set<int> cells;
  for (auto c : view<MemSpace_kind::HOST>(face_cell_ids_.entries)) {
    cells.insert(c);
    c_max = std::max(c, c_max);
  }
  int num_cells = cells.size();
  if (num_cells != c_max + 1) {
    Errors::Message mesg("MeshLogical expects cells to be 0-indexed and contiguous.");
    Exceptions::amanzi_throw(mesg);
  }

  // geometric info that must get set later
  Kokkos::resize(cell_volumes_, num_cells);
  initView(cell_volumes_, -1.);
  Kokkos::resize(face_areas_, face_cell_ids_.size<MemSpace_kind::HOST>());
  initView(face_areas_, -1.);

  std::vector<Entity_ID_List> cell_face_ids_v(num_cells);
  std::vector<std::vector<int>> cell_face_dirs_v(num_cells);
  std::vector<Point_List> cell_face_bisectors_v(num_cells);

  int f_id = 0;
  for (auto& f : face_cell_ids) {
    cell_face_ids_v[f[0]].push_back(f_id);
    cell_face_dirs_v[f[0]].push_back(face_cell_dirs[f_id][0]);
    cell_face_bisectors_v[f[0]].emplace_back(AmanziGeometry::Point());

    if (f.size() > 1 && f[1] >= 0) {
      cell_face_ids_v[f[1]].push_back(f_id);
      cell_face_dirs_v[f[1]].push_back(face_cell_dirs[f_id][1]);
      cell_face_bisectors_v[f[1]].emplace_back(AmanziGeometry::Point());
    }
    f_id++;
  }

  cell_face_ids_ = RaggedArray_DualView<Entity_ID>(cell_face_ids_v);
  cell_face_dirs_ = RaggedArray_DualView<int>(cell_face_dirs_v);
  cell_face_bisectors_ = RaggedArray_DualView<AmanziGeometry::Point>(cell_face_bisectors_v);
}

//
// MeshLogical constructor includes geometry.
//  Includes gravity.
//
//  - cell_volume             : length ncells array of cell volumes
//  - face_cell_list          : length nfaces array of length 2 arrays
//                              defining the topology
//  - face_cell_lengths       : length of the cell-to-face connection
//  - face_area_normals       : length nfaces array of normals of the
//                              face, points from cell 1 to 2 in
//                              face_cell_list topology, magnitude
//                              is area
//  - cell_centroids          : (optional, for plotting) length ncell
//                              array of centroids
//
// Breaks standards following the rest of the mesh infrastructure.
//
MeshLogical::MeshLogical(const Comm_ptr_type& comm,
                         const std::vector<Entity_ID_List>& face_cell_ids,
                         const std::vector<std::vector<int>>& face_cell_dirs,
                         const Double_List& cell_volumes,
                         const Double_List& face_areas,
                         const std::vector<Point_List>& face_cell_bisectors,
                         const Point_List* cell_centroids,
                         const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : MeshLogical(comm, face_cell_ids, face_cell_dirs, plist)
{
  // set geometry
  setLogicalGeometry(&cell_volumes, &face_areas, &face_cell_bisectors, cell_centroids);
}


void
MeshLogical::getLogicalGeometry(Double_List* const cell_volumes,
                                Double_List* const face_areas,
                                std::vector<Point_List>* const face_cell_bisectors,
                                Point_List* const cell_centroids) const
{
  if (cell_volumes) *cell_volumes = asVector(cell_volumes_);
  if (face_areas) *face_areas = asVector(face_areas_);
  if (cell_centroids) *cell_centroids = asVector(cell_centroids_);

  // not sure what would be useful for user here... wait til we have a user and fix it!
  if (face_cell_bisectors) {
    face_cell_bisectors->resize(getNumEntities(Entity_kind::FACE, Parallel_kind::ALL));
  }
}


void
MeshLogical::setLogicalGeometry(Double_List const* const cell_volumes,
                                Double_List const* const face_areas,
                                std::vector<Point_List> const* const face_cell_bisectors,
                                Point_List const* const cell_centroids)
{
  auto n_cells = getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  auto n_faces = getNumEntities(Entity_kind::FACE, Parallel_kind::OWNED);

  if (cell_volumes && n_cells != cell_volumes->size()) {
    Errors::Message mesg("MeshLogical::setLogicalGeometry() called with bad data");
    Exceptions::amanzi_throw(mesg);
  }
  if (face_areas && n_faces != face_areas->size()) {
    Errors::Message mesg("MeshLogical::setLogicalGeometry() called with bad data");
    Exceptions::amanzi_throw(mesg);
  }
  if (face_cell_bisectors && n_faces != face_cell_bisectors->size()) {
    Errors::Message mesg("MeshLogical::setLogicalGeometry() called with bad data");
    Exceptions::amanzi_throw(mesg);
  }
  if (cell_centroids && n_cells != cell_centroids->size()) {
    Errors::Message mesg("MeshLogical::setLogicalGeometry() called with bad data");
    Exceptions::amanzi_throw(mesg);
  }


  if (cell_volumes) vectorToView(cell_volumes_, *cell_volumes);
  if (face_areas) vectorToView(face_areas_, *face_areas);
  if (cell_centroids) vectorToView(cell_centroids_, *cell_centroids);

  if (face_cell_bisectors) {
    std::vector<Point_List> cell_face_bisectors_v(n_cells);

    for (Entity_ID f = 0; f != n_faces; ++f) {
      cEntity_ID_View f_cells;
      getFaceCells(f, Parallel_kind::ALL, f_cells);
      AMANZI_ASSERT((*face_cell_bisectors)[f].size() == f_cells.size());

      for (int c_index = 0; c_index != f_cells.size(); ++c_index) {
        Entity_ID c = f_cells[c_index];
        cell_face_bisectors_v[c].push_back((*face_cell_bisectors)[f][c_index]);
      }
    }
    cell_face_bisectors_ = RaggedArray_DualView<AmanziGeometry::Point>(cell_face_bisectors_v);
  }
}


// testing purposes -- checks if the caches match
bool
MeshLogical::operator==(const MeshLogical& other)
{
  double _eps = 1.e-10;

  if (&other == this) return true;
  if (cell_face_ids_ != other.cell_face_ids_) return false;
  if (face_cell_ids_ != other.face_cell_ids_) return false;

  if (cell_volumes_.size() != other.cell_volumes_.size()) return false;
  for (size_t i = 0; i != cell_volumes_.size(); ++i) {
    if (std::abs(cell_volumes_[i] - other.cell_volumes_[i]) > _eps) return false;
  }

  if (cell_centroids_.size() != other.cell_centroids_.size()) return false;
  for (size_t i = 0; i != cell_centroids_.size(); ++i) {
    if (AmanziGeometry::norm(cell_centroids_[i] - other.cell_centroids_[i]) > _eps) return false;
  }

  if (cell_face_bisectors_.size<MemSpace_kind::HOST>() !=
      other.cell_face_bisectors_.size<MemSpace_kind::HOST>())
    return false;
  for (size_t i = 0; i != cell_face_bisectors_.size<MemSpace_kind::HOST>(); ++i) {
    if (cell_face_bisectors_.getRow<MemSpace_kind::HOST>(i).size() !=
        other.cell_face_bisectors_.getRow<MemSpace_kind::HOST>(i).size())
      return false;
    for (size_t j = 0; j != cell_face_bisectors_.getRow<MemSpace_kind::HOST>(i).size(); ++j)
      if (AmanziGeometry::norm(cell_face_bisectors_.get<MemSpace_kind::HOST>(i, j) -
                               other.cell_face_bisectors_.get<MemSpace_kind::HOST>(i, j)) > _eps)
        return false;
  }

  return true;
}


std::size_t
MeshLogical::getNumEntities(const Entity_kind kind, const Parallel_kind ptype) const
{
  if (kind == Entity_kind::CELL) {
    return cell_face_ids_.size<MemSpace_kind::HOST>();
  } else if (kind == Entity_kind::FACE) {
    return face_cell_ids_.size<MemSpace_kind::HOST>();
  } else {
    return 0;
  }
}

//
// Nodal methods
//
AmanziGeometry::Point
MeshLogical::getNodeCoordinate(const Entity_ID node) const
{
  Errors::Message mesg("There are no nodes in a MeshLogical.");
  Exceptions::amanzi_throw(mesg);
  return AmanziGeometry::Point();
}

void
MeshLogical::getFaceNodes(const Entity_ID f, cEntity_ID_View& nodes) const
{
  Errors::Message mesg("There are no nodes in a MeshLogical.");
  Exceptions::amanzi_throw(mesg);
}

void
MeshLogical::getNodeFaces(const Entity_ID nodeid,
                          const Parallel_kind ptype,
                          cEntity_ID_View& faceids) const
{
  Errors::Message mesg("There are no nodes in a MeshLogical.");
  Exceptions::amanzi_throw(mesg);
}


//
// Cell methods
//
void
MeshLogical::getCellFacesAndBisectors(const Entity_ID cellid,
                                      cEntity_ID_View& faceids,
                                      cPoint_View* const bisectors) const
{
  faceids = cell_face_ids_.getRow<MemSpace_kind::HOST>(cellid);
  if (bisectors) *bisectors = cell_face_bisectors_.getRow<MemSpace_kind::HOST>(cellid);
}

void
MeshLogical::getCellFacesAndDirs(const Entity_ID c,
                                 cEntity_ID_View& faces,
                                 cEntity_Direction_View* const dirs) const
{
  faces = cell_face_ids_.getRow<MemSpace_kind::HOST>(c);
  if (dirs) (*dirs) = cell_face_dirs_.getRow<MemSpace_kind::HOST>(c);
}

double
MeshLogical::getCellVolume(const Entity_ID c) const
{
  return cell_volumes_[c];
}

AmanziGeometry::Point
MeshLogical::getCellCentroid(const Entity_ID c) const
{
  if (cell_centroids_.size() > 0) {
    return cell_centroids_[c];
  } else {
    return AmanziGeometry::Point();
  }
}

//
// Face methods
//
void
MeshLogical::getFaceCells(const Entity_ID f,
                          const Parallel_kind ptype,
                          cEntity_ID_View& cells) const
{
  cells = face_cell_ids_.getRow<MemSpace_kind::HOST>(f);
}

double
MeshLogical::getFaceArea(const Entity_ID f) const
{
  return face_areas_[f];
}

AmanziGeometry::Point
MeshLogical::getFaceCentroid(const Entity_ID f) const
{
  if (cell_centroids_.size() > 0) {
    cEntity_ID_View fcells;
    getFaceCells(f, Parallel_kind::ALL, fcells);
    auto c0 = fcells[0];

    cEntity_ID_View cfaces;
    cPoint_View bisectors;
    getCellFacesAndBisectors(c0, cfaces, &bisectors);

    int i = 0;
    for (; i != cfaces.size(); ++i)
      if (cfaces[i] == f) break;
    AMANZI_ASSERT(i != cfaces.size());
    auto face_centroid_left = getCellCentroid(c0) + bisectors[i];

    if (fcells.size() == 2) {
      auto c1 = fcells[1];
      getCellFacesAndBisectors(c1, cfaces, &bisectors);

      i = 0;
      for (; i != cfaces.size(); ++i)
        if (cfaces[i] == f) break;
      AMANZI_ASSERT(i != cfaces.size());
      auto face_centroid_right = getCellCentroid(c1) + bisectors[i];

      face_centroid_left = (face_centroid_left + face_centroid_right) / 2.;
    }
    return face_centroid_left;
  } else {
    return AmanziGeometry::Point();
  }
}

AmanziGeometry::Point
MeshLogical::getFaceNormal(const Entity_ID f, const Entity_ID c, int* const orientation) const
{
  Entity_ID cc = c;
  if (c == -1) {
    cEntity_ID_View fcells;
    getFaceCells(f, Parallel_kind::ALL, fcells);
    cc = fcells[0];
  }

  cEntity_ID_View cfaces;
  cPoint_View bisectors;
  getCellFacesAndBisectors(cc, cfaces, &bisectors);
  int i = 0;
  for (; i != cfaces.size(); ++i)
    if (cfaces[i] == f) break;

  double bisector_norm = AmanziGeometry::norm(bisectors[i]);
  auto normal = getFaceArea(f) / bisector_norm * bisectors[i];
  if (orientation) *orientation = (cc == 0) ? 1 : -1;
  return normal;
}


//
// Note this works on Mesh, but is less useful for a general mesh
// --------------------------------------------------------------------------------
bool
viewMeshLogical(const Mesh& m, std::ostream& os)
{
  // if (m.getComm()->NumProc() != 1) {
  //   return true;
  // }

  // os << "cell_centroids, volumes =" << std::endl;
  // for (int c=0; c!=m.getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED); ++c) {
  //   os << m.cell_centroid(c) << " " << m.cell_volume(c) << std::endl;
  // }
  // os << "face_connections, areas =" << std::endl;
  // for (int f=0; f!=m.getNumEntities(Entity_kind::FACE, Parallel_kind::OWNED); ++f) {
  //   AmanziMesh::Entity_ID_View fcells;
  //   m.face_get_cells(f, Parallel_kind::ALL, &fcells);
  //   for (auto c : fcells) os << c << " ";
  //   os << m.face_area(f) << std::endl;
  // }

  // os << "cell_sets =" << std::endl;
  // for (auto& r : *m.getGeometricModel()) {
  //   if (r->type() == AmanziGeometry::ENUMERATED) {
  //     AmanziMesh::Entity_ID_View set;
  //     m.get_set_entities(r->name(), AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED, &set);
  //     os << r->name() << " ";
  //     for (auto e : set) os << e << " ";
  //     os << std::endl;
  //   }
  // }
  return false;
}


} // namespace AmanziMesh
} // namespace Amanzi
