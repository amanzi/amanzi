/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  Output

  Silo implementation of an Output object.
*/

#include <locale>
#include <iomanip>

#include "errors.hh"
#include "Point.hh"
#include "Mesh.hh"
#include "OutputFactory.hh"
#include "OutputSilo.hh"

namespace Amanzi {

OutputSilo::OutputSilo(Teuchos::ParameterList& plist,
                       const Teuchos::RCP<const AmanziMesh::MeshHost>& mesh)
  : count_(0),
    mesh_(mesh),
    fid_(nullptr)
{
  if (mesh_->getVisMesh().getSpaceDimension() != 3 ||
      mesh_->getVisMesh().getManifoldDimension() != 3) {
    Errors::Message msg("OutputSilo is untested on non-3D meshes.");
    Exceptions::amanzi_throw(msg);
  }

  formatter_ = OutputFactory::createDirectoryFormatter(plist);
}


// destructor must release file resource on non-finalized
OutputSilo::~OutputSilo()
{
  closeFile_();
}


// open and close files
void
OutputSilo::createTimestep(double time, int cycle)
{
  // check not open
  AMANZI_ASSERT(fid_ == NULL);
  closeFile_();
  int comm_size = mesh_->getComm()->getSize();
  int comm_rank = mesh_->getComm()->getRank();

  // open file
  auto fname = getFilename(count_);

  if (comm_rank == 0) {
    std::cout << "Creating " << fname << " at " << cycle << " time " << time << std::endl;
    fid_ = DBCreate(fname.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_HDF5);
    closeFile_();
  }

  for (int rank = 0; rank != comm_size; ++rank) {
    if (rank == comm_rank) {
      fid_ = DBOpen(fname.c_str(), DB_HDF5, DB_APPEND);

      // directory
      std::string dirname = "/domain_";
      dirname = dirname + std::to_string(rank);
      int ierr = DBMkDir(fid_, dirname.c_str());
      AMANZI_ASSERT(!ierr);
      ierr |= DBSetDir(fid_, dirname.c_str());
      AMANZI_ASSERT(!ierr);

      // -- DBOptions
      DBoptlist* optlist = DBMakeOptlist(3);
      std::string zonelist = "zonelist";
      ierr |= DBAddOption(optlist, DBOPT_PHZONELIST, (void*)zonelist.c_str());
      ierr |= DBAddOption(optlist, DBOPT_CYCLE, &cycle);
      ierr |= DBAddOption(optlist, DBOPT_DTIME, &time);
      AMANZI_ASSERT(!ierr);

      // write mesh
      // -- coordinate names
      // This is optional for now, but we'll give it anyway.
      char* coordnames[3];
      coordnames[0] = (char*)"x-coords";
      coordnames[1] = (char*)"y-coords";
      coordnames[2] = (char*)"z-coords";

      // -- nodal coordinates
      int nnodes = mesh_->getVisMesh().getNumEntities(AmanziMesh::NODE, AmanziMesh::Parallel_kind::ALL);
      std::vector<double> x(nnodes), y(nnodes), z(nnodes);

      for (int i = 0; i != nnodes; ++i) {
        AmanziGeometry::Point xyz = mesh_->getVisMesh().getNodeCoordinate(i);
        x[i] = xyz[0];
        y[i] = xyz[1];
        z[i] = xyz.dim() > 2 ? xyz[2] : 0.0;
      }
      double* coords[3];
      coords[0] = &x[0];
      coords[1] = &y[0];
      coords[2] = &z[0];

      // -- write the base mesh UCD object
      int ncells =
        mesh_->getVisMesh().getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
      ierr |= DBPutUcdmesh(fid_,
                           "mesh",
                           mesh_->getVisMesh().getSpaceDimension(),
                           (char const* const*)coordnames,
                           coords,
                           nnodes,
                           ncells,
                           0,
                           0,
                           DB_DOUBLE,
                           optlist);
      AMANZI_ASSERT(!ierr);

      // -- Construct the silo face-node info.
      // We rely on the mesh having the faces nodes arranged counter-clockwise
      // around the face.  This should be satisfied by AmanziMesh.
      int nfaces = mesh_->getVisMesh().getNumEntities(AmanziMesh::FACE, AmanziMesh::Parallel_kind::ALL);
      std::vector<int> face_node_counts(nfaces);
      std::vector<char> ext_faces(nfaces, 0x0);
      std::vector<int> face_node_list;

      // const auto& nmap = mesh_->getVisMesh().node_map(true);
      for (int f = 0; f != nfaces; ++f) {
        auto fnodes = mesh_->getVisMesh().getFaceNodes(f);
        face_node_counts[f] = fnodes.size();
        face_node_list.insert(face_node_list.end(), begin(fnodes), end(fnodes));

        auto fcells = mesh_->getVisMesh().getFaceCells(f);
        if (fcells.size() == 1) { ext_faces[f] = 0x1; }
      }

      // -- Construct the silo cell-face info
      std::vector<int> cell_face_counts(ncells);
      std::vector<int> cell_face_list;

      for (int c = 0; c != ncells; ++c) {
        auto cfaces = mesh_->getVisMesh().getCellFaces(c);
        cell_face_counts[c] = cfaces.size();
        cell_face_list.insert(cell_face_list.end(), begin(cfaces), end(cfaces));
      }

      ierr |= DBPutPHZonelist(fid_,
                              "zonelist",
                              nfaces,
                              &face_node_counts[0],
                              face_node_list.size(),
                              &face_node_list[0],
                              &ext_faces[0],
                              ncells,
                              &cell_face_counts[0],
                              cell_face_list.size(),
                              &cell_face_list[0],
                              0,
                              0,
                              ncells - 1,
                              optlist);
      AMANZI_ASSERT(!ierr);

      // -- clean up (could be done on finalize if needed? --etc)
      DBFreeOptlist(optlist);
      closeFile_();
    }
    mesh_->getComm()->barrier();
  }

  if (comm_rank == 0) {
    fid_ = DBOpen(fname.c_str(), DB_HDF5, DB_APPEND);
    DBSetDir(fid_, "/");
    std::string meshname("mesh");
    std::vector<std::string> meshnames_str(comm_size);
    std::vector<char*> meshnames(comm_size, nullptr);
    std::vector<int> meshtypes(comm_size);
    for (int i = 0; i != comm_size; ++i) {
      meshnames_str[i] = std::string("/domain_") + std::to_string(i) + "/mesh";
      meshnames[i] = const_cast<char*>(meshnames_str[i].c_str());
      meshtypes[i] = DB_UCDMESH;
    }

    int ierrl = DBPutMultimesh(
      fid_, meshname.c_str(), comm_size, (char**)meshnames.data(), (int*)meshtypes.data(), nullptr);
    AMANZI_ASSERT(!ierrl);
    closeFile_();
  }
}

void
OutputSilo::finalizeTimestep()
{
  count_++;
}


void
OutputSilo::closeFile_() const
{
  if (fid_) {
    DBClose(fid_);
    fid_ = nullptr;
  }
}

std::string
OutputSilo::fixName_(const std::string& s) const
{
  int n = s.size(), wp = 0;
  std::vector<char> result(n);
  for (int i = 0; i < n; ++i) {
    if (!std::isalnum(s[i])) {
      result[wp++] = '_';
    } else {
      result[wp++] = s[i];
    }
  }
  return std::string(&result[0], &result[wp]);
}

int
OutputSilo::getSiloLocation_(AmanziMesh::Entity_kind entity_kind) const {
  switch(entity_kind) {
    case AmanziMesh::Entity_kind::CELL :
      return DB_ZONECENT; break;
    case AmanziMesh::Entity_kind::NODE :
      return DB_NODECENT; break;
    default:
      Errors::Message msg("OutputSilo only knows how to write NODEs and CELLs.");
      Exceptions::amanzi_throw(msg);
      return -1;
  }
}

} // namespace Amanzi
