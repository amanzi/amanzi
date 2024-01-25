/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: William Perkins, Ethan Coon
*/

#include <boost/format.hpp>

#include "RegionLogical.hh"
#include "MeshException.hh"
#include "MeshFrameworkFactory.hh"
#include "FileFormat.hh"

#include "MeshExtractedManifold.hh"
// #include "MeshColumn.hh"
#include "Mesh_simple.hh"
#include "RegionLabeledSet.hh"

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
// MeshFrameworkFactory:: constructors / destructor
// -------------------------------------------------------------
MeshFrameworkFactory::MeshFrameworkFactory(const Comm_ptr_type& comm,
                                           const Teuchos::RCP<AmanziGeometry::GeometricModel>& gm,
                                           const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : comm_(comm), preference_(), plist_(plist), gm_(gm)
{
  if (comm_ == Teuchos::null) { comm_ = Amanzi::getDefaultComm(); }

  set_preference(default_preference());

  if (plist_ == Teuchos::null) { plist_ = Teuchos::rcp(new Teuchos::ParameterList()); }

  vo_ = Teuchos::rcp(new VerboseObject(comm_, "Amanzi::MeshFactory", *plist_));

  // submesh parameter
  extraction_method_ = "none";
  if (plist != Teuchos::null) {
    const auto& umesh = plist->sublist("unstructured");
    if (umesh.isSublist("submesh")) {
      const auto& submesh = umesh.sublist("submesh");
      if (submesh.isParameter("extraction method")) {
        extraction_method_ = submesh.get<std::string>("extraction method");
      }
    }
  }
}


// -------------------------------------------------------------
// MeshFrameworkFactory::preference
// -------------------------------------------------------------
void
MeshFrameworkFactory::set_preference(const Preference& pref)
{
  preference_ = filter_preference(pref);
  if (preference_.size() == 0) {
    Exceptions::amanzi_throw(
      FrameworkMessage("Preference requested includes no frameworks that are enabled."));
  }
}


// -------------------------------------------------------------
// MeshFrameworkFactory::create
// -------------------------------------------------------------
Teuchos::RCP<MeshFramework>
MeshFrameworkFactory::create(const std::string& filename)
{
  auto tab = vo_->getOSTab();

  // check the file format
  FileFormat fmt = fileFormatFromFilename(*comm_, filename);

  if (fmt == FileFormat::UNKNOWN) {
    FileMessage e(
      std::string(boost::str(boost::format("%s: unknown file format") % filename)).c_str());
    Exceptions::amanzi_throw(e);
  }

  for (auto p : preference_) {
    int nproc = comm_->getSize();

#ifdef HAVE_MESH_MSTK
    if (p == Framework::MSTK) {
      if (vo_->os_OK(Teuchos::VERB_HIGH)) {
        *vo_->os() << "Creating mesh from file \"" << filename << "\" using format \"MSTK\""
                   << std::endl;
      }
      if ((nproc == 1 && fmt == FileFormat::EXODUS_II) ||
          (nproc > 1 && (fmt == FileFormat::EXODUS_II || fmt == FileFormat::NEMESIS))) {
        auto mesh = Teuchos::rcp(
          new Mesh_MSTK(filename, comm_, gm_, Teuchos::rcp(new Teuchos::ParameterList(*plist_))));
        return mesh;
      }
    }
#endif

#ifdef HAVE_MESH_MOAB
    if (p == Framework::MOAB) {
      if (vo_->os_OK(Teuchos::VERB_HIGH)) {
        *vo_->os() << "Creating mesh from file \"" << filename << "\" using format \"MOAB\""
                   << std::endl;
      }
      if (fmt == FileFormat::MOAB_HDF5 || (nproc == 1 && fmt == FileFormat::EXODUS_II)) {
        auto mesh = Teuchos::rcp(
          new Mesh_MOAB(filename, comm_, gm_, Teuchos::rcp(new Teuchos::ParameterList(*plist_))));
        return mesh;
      }
    }
#endif
  }

  Message m(
    "No construct was found in preferences that is available and can read this file/file format.");
  Exceptions::amanzi_throw(m);
  return Teuchos::null;
}


// -------------------------------------------------------------
// This creates a mesh by generating a block of hexahedral cells.
//
// Hopefully, if any one process has an error, all processes will
// throw an Mesh::Message exception.
//
// Input:
//   (x0, y0, z0) left-bottom-front corner
//   (x1, y1, z1) right-top-back corner
//   nx x ny x nz number of cells in the box mesh
//
// Returns mesh instance
// -------------------------------------------------------------
Teuchos::RCP<MeshFramework>
MeshFrameworkFactory::create(const double x0,
                             const double y0,
                             const double z0,
                             const double x1,
                             const double y1,
                             const double z1,
                             const int nx,
                             const int ny,
                             const int nz)
{
  int nproc = comm_->getSize();

  for (auto p : preference_) {
    if (p == Framework::SIMPLE && nproc == 1) {
      auto mesh = Teuchos::rcp(new Mesh_simple(x0,
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
                                               Teuchos::rcp(new Teuchos::ParameterList(*plist_))));
      return mesh;
    }

#ifdef HAVE_MESH_MSTK
    if (p == Framework::MSTK) {
      if (vo_->os_OK(Teuchos::VERB_HIGH)) {
        *vo_->os() << "Creating 3D block mesh using format \"MSTK\"" << std::endl;
      }
      auto mesh = Teuchos::rcp(new Mesh_MSTK(x0,
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
                                             Teuchos::rcp(new Teuchos::ParameterList(*plist_))));
      return mesh;
    }
#endif
  }

  Exceptions::amanzi_throw(
    Message("No construct was found in preferences that is available and can generate in 3D."));
  return Teuchos::null;
}


// -------------------------------------------------------------
// This creates a mesh by generating a block of quadrilateral cells.
//
// Hopefully, if any one process has an error, all processes will
// throw an Mesh::Message exception.
//
// Input:
//   (x0, y0) left-bottom-front corner
//   (x1, y1) right-top-back corner
//   nx x ny number of cells in the box mesh
//
// Returns mesh instance
// -------------------------------------------------------------
Teuchos::RCP<MeshFramework>
MeshFrameworkFactory::create(const double x0,
                             const double y0,
                             const double x1,
                             const double y1,
                             const int nx,
                             const int ny)
{
  for (auto p : preference_) {
#ifdef HAVE_MESH_MSTK
    if (p == Framework::MSTK) {
      if (vo_->os_OK(Teuchos::VERB_HIGH)) {
        *vo_->os() << "Creating 2D block mesh using format \"MSTK\"" << std::endl;
      }
      auto mesh = Teuchos::rcp(new Mesh_MSTK(
        x0, y0, x1, y1, nx, ny, comm_, gm_, Teuchos::rcp(new Teuchos::ParameterList(*plist_))));
      return mesh;
    }
#endif
  }

  Exceptions::amanzi_throw(
    Message("No construct was found in preferences that is available and can generate in 2D."));
  return Teuchos::null;
}


// -------------------------------------------------------------
// This creates a mesh by generating a block of quadrilateral
// or hexahedral cells, using a parameter list with the limits
// and cell counts.
// -------------------------------------------------------------
Teuchos::RCP<MeshFramework>
MeshFrameworkFactory::create(const Teuchos::ParameterList& parameter_list)
{
  Message e("MeshFrameworkFactory::create: error: ");
  if (!parameter_list.isParameter("number of cells")) {
    e.add_data("missing \"number of cells\"");
    Exceptions::amanzi_throw(e);
  }
  auto ncells = parameter_list.get<Teuchos::Array<int>>("number of cells");

  if (!parameter_list.isParameter("domain low coordinate")) {
    e.add_data("missing \"domain low coordinate\"");
    Exceptions::amanzi_throw(e);
  }
  auto low = parameter_list.get<Teuchos::Array<double>>("domain low coordinate");

  if (!parameter_list.isParameter("domain high coordinate")) {
    e.add_data("missing \"domain high coordinate\"");
    Exceptions::amanzi_throw(e);
  }
  auto high = parameter_list.get<Teuchos::Array<double>>("domain high coordinate");

  if (!((ncells.size() == 2) || (ncells.size() == 3))) {
    e.add_data("invalid size \"number of cells\", must be 2 or 3");
    Exceptions::amanzi_throw(e);
  }

  if (low.size() != ncells.size()) {
    e.add_data("invalid size \"domain low coordinate\", must match \"number of cells\"");
    Exceptions::amanzi_throw(e);
  }
  if (high.size() != ncells.size()) {
    e.add_data("invalid size \"domain high coordinate\", must match \"number of cells\"");
    Exceptions::amanzi_throw(e);
  }

  if (ncells.size() == 2) {
    return create(low[0], low[1], high[0], high[1], ncells[0], ncells[1]);
  } else {
    return create(
      low[0], low[1], low[2], high[0], high[1], high[2], ncells[0], ncells[1], ncells[2]);
  }
}


// -------------------------------------------------------------
// This creates a mesh by extracting subsets of entities from an
// existing mesh possibly flattening it by removing the last
// dimension or (in the future) extruding it when it makes sense
// -------------------------------------------------------------
Teuchos::RCP<MeshFramework>
MeshFrameworkFactory::create(const Teuchos::RCP<const Mesh>& inmesh,
                             const MeshFramework::cEntity_ID_View& setids,
                             const Entity_kind setkind,
                             const bool flatten)
{
  // we have sane defaults from the parent mesh for some things
  auto gm = Teuchos::RCP<const AmanziGeometry::GeometricModel>(gm_);
  if (gm == Teuchos::null) { gm = inmesh->getGeometricModel(); }

  // create the comm via split
  Comm_ptr_type comm = Teuchos::null;
  auto parent_comm = Teuchos::rcp_dynamic_cast<const MpiComm_type>(inmesh->getComm());
  if (parent_comm != Teuchos::null && plist_->get<bool>("create subcommunicator", false)) {
    // mpi comm split
    MPI_Comm child_comm;
    bool empty = setids.size() == 0;
    int ierr = MPI_Comm_split(*parent_comm->getRawMpiComm(), empty ? MPI_UNDEFINED : 1, 0, &child_comm);
    if (ierr) {
      Errors::Message msg("Error in MPI_Comm_split in creating extracted mesh.");
      Exceptions::amanzi_throw(msg);
    }
    if (!empty) comm = Teuchos::rcp(new MpiComm_type(child_comm));
  } else {
    // note, this could be parent_comm or it could be COMM_SELF or any other
    // comm that the factory was made with.
    comm = comm_;
  }

  // check that ids are unique
  for (auto it1 = setids.begin(); it1 != setids.end(); ++it1) {
    for (auto it2 = it1 + 1; it2 != setids.end(); ++it2) {
      if (*it2 == *it1) {
        Exceptions::amanzi_throw(Message("Extrated mesh has geometrically identical elements."));
      }
    }
  }

  // extract
  for (auto p : preference_) {
#ifdef HAVE_MESH_MSTK
    if (p == Framework::MSTK) {
      *vo_->os() << "Pref = MSTK\n";
      if (vo_->os_OK(Teuchos::VERB_HIGH)) {
        *vo_->os() << "Creating extracted mesh using format \"MSTK\" with " << setids.size()
                   << " local entities." << std::endl;
      }

      if (comm != Teuchos::null) {
        auto mesh = Teuchos::rcp(
          new Mesh_MSTK(inmesh->getMeshFramework(), setids, setkind, flatten, comm, gm, plist_));
        return mesh;
      } else {
        return Teuchos::null;
      }
    }
#endif
  }

  Exceptions::amanzi_throw(
    Message("No construct was found in preferences that is available and can extract."));
  return Teuchos::null;
}

// -------------------------------------------------------------
// This creates a mesh by extracting subsets of entities from an existing
//  mesh possibly flattening it by removing the last dimension or (in the
// future) extruding it when it makes sense
// -------------------------------------------------------------
Teuchos::RCP<MeshFramework>
MeshFrameworkFactory::create(const Teuchos::RCP<const Mesh>& inmesh,
                             const std::vector<std::string>& setnames,
                             const Entity_kind setkind,
                             const bool flatten)
{
  if (extraction_method_ == "manifold mesh") {
    const auto& comm = inmesh->getComm();
    const auto& gm = inmesh->getGeometricModel();

    // greedy solution for more than one sets. A new region is forced into GM
    std::string setname(setnames[0]);
    if (setnames.size() > 1) {
      for (int i = 1; i < setnames.size(); ++i) setname += "_" + setnames[i];
      auto rgn = Teuchos::rcp(new AmanziGeometry::RegionLogical(setname, -1, "union", setnames));
      Teuchos::rcp_const_cast<AmanziGeometry::GeometricModel>(gm)->AddRegion(rgn);
    }

    auto mesh =
      Teuchos::rcp(new MeshExtractedManifold(inmesh, setname, setkind, comm, gm, plist_, flatten));

    return mesh;

  } else {
    MeshFramework::Entity_ID_View ids;
    auto inmesh_host = onMemSpace<MemSpace_kind::HOST>(inmesh);
    for (auto name : setnames) {
      auto ids_l = inmesh->getSetEntities(name, setkind, Parallel_kind::OWNED);
      ids.insert(ids.end(), ids_l.begin(), ids_l.end());
    }
    return create(inmesh, ids, setkind, flatten);
  }
}


} // namespace AmanziMesh
} // namespace Amanzi
