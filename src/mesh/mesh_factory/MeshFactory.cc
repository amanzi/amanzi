/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      William Perkins, Ethan Coon
*/

//! <MISSING_ONELINE_DOCSTRING>

#include <boost/format.hpp>

#include "MeshException.hh"
#include "MeshFactory.hh"
#include "FileFormat.hh"

#include "Mesh_simple.hh"

#ifdef HAVE_MESH_MSTK
#  include "Mesh_MSTK.hh"
#endif
#ifdef HAVE_MESH_MOAB
#  include "Mesh_MOAB.hh"
#endif

namespace Amanzi {
namespace AmanziMesh {


// -------------------------------------------------------------
//  class MeshFactory
// -------------------------------------------------------------

// -------------------------------------------------------------
// MeshFactory:: constructors / destructor
// -------------------------------------------------------------
MeshFactory::MeshFactory(
  const Comm_ptr_type& comm,
  const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm,
  const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : comm_(comm), gm_(gm), plist_(plist), preference_()
{
  if (comm_ == Teuchos::null) { comm_ = Amanzi::getDefaultComm(); }

  set_preference(default_preference());

  if (plist_ == Teuchos::null) {
    plist_ = Teuchos::rcp(new Teuchos::ParameterList());
  }
}


// -------------------------------------------------------------
// MeshFactory::preference
// -------------------------------------------------------------
void
MeshFactory::set_preference(const Preference& pref)
{
  preference_ = filter_preference(pref);
  if (preference_.size() == 0) {
    Exceptions::amanzi_throw(FrameworkMessage(
      "Preference requested includes no frameworks that are enabled."));
  }
}

// -------------------------------------------------------------
// MeshFactory::create
// -------------------------------------------------------------
Teuchos::RCP<Mesh>
MeshFactory::create(const std::string& filename, const bool request_faces,
                    const bool request_edges)
{
  // check the file format
  FileFormat fmt = fileFormatFromFilename(*comm_, filename);

  if (fmt == FileFormat::UNKNOWN) {
    FileMessage e(
      std::string(
        boost::str(boost::format("%s: unknown file format") % filename))
        .c_str());
    Exceptions::amanzi_throw(e);
  }

  Teuchos::RCP<Mesh> mesh = Teuchos::null;
  for (auto p : preference_) {
    int nproc = comm_->getSize();

#ifdef HAVE_MESH_MSTK
    if (p == Framework::MSTK) {
      if ((nproc == 1 && fmt == FileFormat::EXODUS_II) ||
          (nproc > 1 &&
           (fmt == FileFormat::EXODUS_II || fmt == FileFormat::NEMESIS))) {
        mesh = Teuchos::rcp(new Mesh_MSTK(
          filename, comm_, gm_, plist_, request_faces, request_edges));
        break;
      }
    }
#endif

#ifdef HAVE_MESH_MOAB
    if (p == Framework::MOAB) {
      if (fmt == FileFormat::MOAB_HDF5 ||
          (nproc == 1 && fmt == FileFormat::EXODUS_II)) {
        mesh = Teuchos::rcp(new Mesh_MOAB(
          filename, comm_, gm_, plist_, request_faces, request_edges));
        break;
      }
    }
#endif
  }

  if (mesh != Teuchos::null) {
    // Initialize the cache after constructing the mesh
    return mesh;
  }

  Message m("No construct was found in preferences that is available and can "
            "read this file/file format.");
  Exceptions::amanzi_throw(m);
  return Teuchos::null;
}

/**
 * Coellective
 *
 * This creates a mesh by generating a block of hexahedral cells.
 *
 * Hopefully, if any one process has an error, all processes will
 * throw an Mesh::Message exception.
 *
 * @param x0 origin x-coordinate
 * @param y0 origin y-coordinate
 * @param z0 origin z-coordinate
 * @param x1 maximum x-coordinate
 * @param y1 maximum y-coordinate
 * @param z1 maximum z-coordinate
 * @param nx number of cells in the x-direction
 * @param ny number of cells in the y-direction
 * @param nz number of cells in the z-direction
 *
 * @return mesh instance
 */
Teuchos::RCP<Mesh>
MeshFactory::create(const double x0, const double y0, const double z0,
                    const double x1, const double y1, const double z1,
                    const int nx, const int ny, const int nz,
                    const bool request_faces, const bool request_edges)
{
  int nproc = comm_->getSize();

  Teuchos::RCP<Mesh> mesh = Teuchos::null;

  for (auto p : preference_) {
    if (p == Framework::SIMPLE && nproc == 1) {
      mesh = Teuchos::rcp(new Mesh_simple(x0,
                                          y0,
                                          z0,
                                          x1,
                                          y1,
                                          z1,
                                          nx,
                                          ny,
                                          nz,
                                          comm_,
                                          gm_,
                                          plist_,
                                          request_faces,
                                          request_edges));
      break;
    }

#ifdef HAVE_MESH_MSTK
    if (p == Framework::MSTK) {
      mesh = Teuchos::rcp(new Mesh_MSTK(x0,
                                        y0,
                                        z0,
                                        x1,
                                        y1,
                                        z1,
                                        nx,
                                        ny,
                                        nz,
                                        comm_,
                                        gm_,
                                        plist_,
                                        request_faces,
                                        request_edges));
      break;
    }
#endif
  }

  if (mesh != Teuchos::null) {
    // Initialize the cache after constructing the mesh
    return mesh;
  }

  Exceptions::amanzi_throw(Message("No construct was found in preferences that "
                                   "is available and can generate in 3D."));
  return Teuchos::null;
}


/**
 * Coellective
 *
 * This creates a mesh by generating a block of quadrilateral cells.
 *
 * Hopefully, if any one process has an error, all processes will
 * throw an Mesh::Message exception.
 *
 * @param x0 origin x-coordinate
 * @param y0 origin y-coordinate
 * @param x1 maximum x-coordinate
 * @param y1 maximum y-coordinate
 * @param nx number of cells in the x-direction
 * @param ny number of cells in the y-direction
 *
 * @return mesh instance
 */
Teuchos::RCP<Mesh>
MeshFactory::create(const double x0, const double y0, const double x1,
                    const double y1, const int nx, const int ny,
                    const bool request_faces, const bool request_edges)
{
  int nproc = comm_->getSize();

  Teuchos::RCP<Mesh> mesh = Teuchos::null;

  for (auto p : preference_) {
#ifdef HAVE_MESH_MSTK
    if (p == Framework::MSTK) {
      mesh = Teuchos::rcp(new Mesh_MSTK(x0,
                                        y0,
                                        x1,
                                        y1,
                                        nx,
                                        ny,
                                        comm_,
                                        gm_,
                                        plist_,
                                        request_faces,
                                        request_edges));
      break;
    }
#endif
  }

  if (mesh != Teuchos::null) {
    // Initialize the cache after constructing the mesh
    return mesh;
  }

  Exceptions::amanzi_throw(Message("No construct was found in preferences that "
                                   "is available and can generate in 2D."));
  return Teuchos::null;
}

/**
 * This creates a mesh by generating a block of
 * quadrilateral/hexahedral cells, but using a parameter list with the
 * limits and cell counts.
 *
 * @param parameter_list
 *
 * @return
 */
Teuchos::RCP<Mesh>
MeshFactory::create(const Teuchos::ParameterList& parameter_list,
                    const bool request_faces, const bool request_edges)
{
  Message e("MeshFactory::create: error: ");
  if (!parameter_list.isParameter("number of cells")) {
    e.add_data("missing \"number of cells\"");
    Exceptions::amanzi_throw(e);
  }
  auto ncells = parameter_list.get<Teuchos::Array<int>>("number of cells");

  if (!parameter_list.isParameter("domain low coordinate")) {
    e.add_data("missing \"domain low coordinate\"");
    Exceptions::amanzi_throw(e);
  }
  auto low =
    parameter_list.get<Teuchos::Array<double>>("domain low coordinate");

  if (!parameter_list.isParameter("domain high coordinate")) {
    e.add_data("missing \"domain high coordinate\"");
    Exceptions::amanzi_throw(e);
  }
  auto high =
    parameter_list.get<Teuchos::Array<double>>("domain high coordinate");

  if (!((ncells.size() == 2) || (ncells.size() == 3))) {
    e.add_data("invalid size \"number of cells\", must be 2 or 3");
    Exceptions::amanzi_throw(e);
  }

  if (low.size() != ncells.size()) {
    e.add_data(
      "invalid size \"domain low coordinate\", must match \"number of cells\"");
    Exceptions::amanzi_throw(e);
  }
  if (high.size() != ncells.size()) {
    e.add_data("invalid size \"domain high coordinate\", must match \"number "
               "of cells\"");
    Exceptions::amanzi_throw(e);
  }

  if (ncells.size() == 2) {
    return create(low[0],
                  low[1],
                  high[0],
                  high[1],
                  ncells[0],
                  ncells[1],
                  request_faces,
                  request_edges);
  } else {
    return create(low[0],
                  low[1],
                  low[2],
                  high[0],
                  high[1],
                  high[2],
                  ncells[0],
                  ncells[1],
                  ncells[2],
                  request_faces,
                  request_edges);
  }
}


/**
 * This creates a mesh by extracting subsets of entities from an existing
 * mesh possibly flattening it by removing the last dimension or (in the
 * future) extruding it when it makes sense
 *
 * @param inmesh
 * @param setids
 * @param setkind
 * @param flatten
 * @param extrude
 *
 * @return
 */
Teuchos::RCP<Mesh>
MeshFactory::create(const Teuchos::RCP<const Mesh>& inmesh,
                    const Entity_ID_List& setids,
                    const Entity_kind setkind, const bool flatten,
                    const bool request_faces, const bool request_edges)
{
  // we have sane defaults from the parent mesh for some things
  auto gm = Teuchos::RCP<const AmanziGeometry::GeometricModel>(gm_);
  if (gm == Teuchos::null) { gm = inmesh->geometric_model(); }

  auto comm = Comm_ptr_type(comm_);
  if (comm == Teuchos::null) { comm = inmesh->get_comm(); }

  Teuchos::RCP<Mesh> mesh = Teuchos::null;
  // extract
  for (auto p : preference_) {
#ifdef HAVE_MESH_MSTK
    if (p == Framework::MSTK) {
      return Teuchos::rcp(new Mesh_MSTK(inmesh,
                                        setids,
                                        setkind,
                                        flatten,
                                        comm,
                                        gm,
                                        plist_,
                                        request_faces,
                                        request_edges));
    }
#endif
  }

  if (mesh != Teuchos::null) {
    // Initialize the cache after constructing the mesh
    return mesh;
  }

  Exceptions::amanzi_throw(Message("No construct was found in preferences that "
                                   "is available and can extract."));
  return Teuchos::null;
}

/**
 * This creates a mesh by extracting subsets of entities from an existing
 * mesh possibly flattening it by removing the last dimension or (in the
 * future) extruding it when it makes sense
 *
 * @param inmesh
 * @param setnames
 * @param setkind
 * @param flatten
 * @param extrude
 *
 * @return
 */
Teuchos::RCP<Mesh>
MeshFactory::create(const Teuchos::RCP<const Mesh>& inmesh,
                    const std::vector<std::string>& setnames,
                    const Entity_kind setkind, const bool flatten,
                    const bool request_faces, const bool request_edges)
{
  Entity_ID_List ids;
  for (auto name : setnames) {
    Entity_ID_List ids_l;
    inmesh->get_set_entities(name, setkind, Parallel_type::OWNED, ids_l);
    ids.insert(ids.end(),ids_l.begin(),ids_l.end()); 
  }
  return create(inmesh, ids, setkind, flatten, request_faces, request_edges);
}
} // namespace AmanziMesh
} // namespace Amanzi
