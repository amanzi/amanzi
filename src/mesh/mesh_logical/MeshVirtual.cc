/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.
*/

#include <set>

#include "MeshVirtual.hh"

namespace Amanzi {
namespace AmanziMesh {

/* ******************************************************************
* Algorithms
****************************************************************** */
std::pair<double, AmanziGeometry::Point>
MeshVirtualAlgorithms::computeCellGeometry(const Mesh& mesh, const Entity_ID c) const
{
  double volume = mesh.getCellVolume(c);
  AmanziGeometry::Point centroid = mesh.getCellCentroid(c);
  return std::make_pair(volume, centroid);
}


std::tuple<double, AmanziGeometry::Point, cPoint_View>
MeshVirtualAlgorithms::computeFaceGeometry(const Mesh& mesh, const Entity_ID f) const
{
  double area = mesh.getFaceArea(f);
  AmanziGeometry::Point centroid = mesh.getFaceCentroid(f);

  MeshHost::cEntity_ID_View fcells;
  mesh.getFaceCells(f, fcells);
  Point_View normals("normals", fcells.size());
  for (int i = 0; i != fcells.size(); ++i) {
    normals[i] = mesh.getFaceNormal(f, fcells[i], nullptr);
  }

  return std::make_tuple(area, centroid, normals);
}


std::pair<AmanziGeometry::Point, AmanziGeometry::Point>
MeshVirtualAlgorithms::computeEdgeGeometry(const Mesh& mesh, const Entity_ID e) const
{
  Errors::Message msg("There are no edges in a MeshVirtual.");
  Exceptions::amanzi_throw(msg);
  return std::make_pair(AmanziGeometry::Point(), AmanziGeometry::Point());
}


double
MeshVirtualAlgorithms::computeCellVolume(const Mesh& mesh, const Entity_ID c) const
{
  return static_cast<const MeshVirtual*>(mesh.getMeshFramework().get())->getCellVolume(c);
}


AmanziGeometry::Point
MeshVirtualAlgorithms::computeCellCentroid(const Mesh& mesh, const Entity_ID c) const
{
  return static_cast<const MeshVirtual*>(mesh.getMeshFramework().get())->getCellCentroid(c);
}


double
MeshVirtualAlgorithms::computeFaceArea(const Mesh& mesh, const Entity_ID f) const
{
  return static_cast<const MeshVirtual*>(mesh.getMeshFramework().get())->getFaceArea(f);
}

AmanziGeometry::Point
MeshVirtualAlgorithms::computeFaceCentroid(const Mesh& mesh, const Entity_ID f) const
{
  return static_cast<const MeshVirtual*>(mesh.getMeshFramework().get())->getFaceCentroid(f);
}

AmanziGeometry::Point
MeshVirtualAlgorithms::computeFaceNormal(const Mesh& mesh,
                                         const Entity_ID f,
                                         const Entity_ID c,
                                         int* const orientation) const
{
  return static_cast<const MeshVirtual*>(mesh.getMeshFramework().get())
    ->getFaceNormal(f, c, orientation);
}


/* ******************************************************************
* Constructor
****************************************************************** */
MeshVirtual::MeshVirtual(const Comm_ptr_type& comm,
                         const Teuchos::RCP<Teuchos::ParameterList>& plist,
                         const std::vector<Entity_ID_List>& face_cells,
                         const std::vector<AmanziGeometry::Point>& cell_centroids,
                         const std::vector<double>& cell_volumes,
                         const std::vector<AmanziGeometry::Point>& face_centroids,
                         const std::vector<AmanziGeometry::Point>& face_normals)


  : MeshFramework(comm, Teuchos::null, plist),
    cell_centroids_(cell_centroids),
    cell_volumes_(cell_volumes),
    face_centroids_(face_centroids),
    face_normals_(face_normals)
{
  int d = cell_centroids[0].dim();
  setSpaceDimension(d);
  setManifoldDimension(d);

  int ncells = cell_centroids.size();
  int nfaces = face_cells.size();

  num_entities_[CELL] = ncells;
  num_entities_[FACE] = nfaces;
  num_entities_[NODE] = 0;

  // face -> cells
  std::vector<Entity_ID_List> face_cells_v(face_cells);
  for (int f = 0; f < nfaces; ++f) {
    if (face_cells_v.size() == 2) {
      if (face_cells_v[f][0] > face_cells_v[f][1]) {
        std::swap(face_cells_v[f][0], face_cells_v[f][1]);
      }
    }
  }
  face_cell_ids_ = RaggedArray_DualView<Entity_ID>(face_cells_v);

  // cell -> faces 
  // -- internal faces go first
  std::vector<Entity_ID_List> cell_face_ids_v(ncells);
  std::vector<std::vector<int>> cell_face_dirs_v(ncells);

  int fid = 0;
  for (int n = 0; n < face_cells.size(); ++n) {
    int c1 = face_cells[n][0];

    if (face_cells[n].size() == 1) {
      cell_face_ids_v[c1].push_back(fid);
      cell_face_dirs_v[c1].push_back(1);

    } else {
      int c2 = face_cells[n][1];
      cell_face_ids_v[c1].push_back(fid);
      cell_face_ids_v[c2].push_back(fid);

      cell_face_dirs_v[c1].push_back(1);
      cell_face_dirs_v[c2].push_back(-1);
    }
    fid++;
  }

  cell_face_ids_ = RaggedArray_DualView<Entity_ID>(cell_face_ids_v);
  cell_face_dirs_ = RaggedArray_DualView<int>(cell_face_dirs_v);

  if (vo_.get() && vo_->getVerbLevel() >= Teuchos::VERB_LOW) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "\nvirtual mesh with " << ncells << " cells and " << nfaces << " faces" << std::endl;
  }
}


/* ******************************************************************
* Get nodes of a cell
****************************************************************** */
void
MeshVirtual::getCellNodes(const Entity_ID c, cEntity_ID_View& nodes) const
{
  Errors::Message mesg("No nodes in MeshVirtual.");
  Exceptions::amanzi_throw(mesg);
}


/* ******************************************************************
* Get face normal
****************************************************************** */
AmanziGeometry::Point
MeshVirtual::getFaceNormal(const Entity_ID f, const Entity_ID c, int* orientation) const
{
  auto cells = face_cell_ids_.getRowUnmanaged<MemSpace_kind::HOST>(f);

  if (cells.size() == 1 || c < 0) {
    if (orientation) *orientation = 1;
    return face_normals_[f];
  }

  if (c == cells[0]) {
    if (orientation) *orientation = 1;
    return face_normals_[f];
  }
   
  if (c == cells[1]) {
    if (orientation) *orientation = -1;
    return -face_normals_[f];
  }

  AMANZI_ASSERT(false);
  return AmanziGeometry::Point(2);
}


/* ******************************************************************
* face -> cells map
****************************************************************** */
void
MeshVirtual::getFaceCells(const Entity_ID f, cEntity_ID_View& cells) const
{
  cells = face_cell_ids_.getRowUnmanaged<MemSpace_kind::HOST>(f);
}


/* ******************************************************************
* Get nodes of face
****************************************************************** */
void
MeshVirtual::getFaceNodes(const Entity_ID f, cEntity_ID_View& nodes) const
{
  Errors::Message mesg("No nodes in MeshVirtual.");
  Exceptions::amanzi_throw(mesg);
}


/* ******************************************************************
* Faces of type 'ptype' connected to a node.
****************************************************************** */
void
MeshVirtual::getNodeFaces(const Entity_ID n, cEntity_ID_View& faces) const
{
  Errors::Message msg("No nodes in MeshVirtual.");
  Exceptions::amanzi_throw(msg);
}


/* ******************************************************************
* Node coordinates
****************************************************************** */
AmanziGeometry::Point 
MeshVirtual::getNodeCoordinate(const Entity_ID n) const
{
  Errors::Message msg("No nodes in MeshVirtual.");
  Exceptions::amanzi_throw(msg);
  return AmanziGeometry::Point(2);
}


/* ******************************************************************
* Get faces of a cell and directions in which it is used. 
* The result must be cached in the base class to reduce cost.
****************************************************************** */
void
MeshVirtual::getCellFacesAndDirs(const Entity_ID c,
                                 cEntity_ID_View& faces,
                                 cDirection_View* const dirs) const
{
  faces = cell_face_ids_.getRowUnmanaged<MemSpace_kind::HOST>(c);
  if (dirs) *dirs = cell_face_dirs_.getRowUnmanaged<MemSpace_kind::HOST>(c);
}

} // namespace AmanziMesh
} // namespace Amanzi
