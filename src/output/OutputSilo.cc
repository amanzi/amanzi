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
#include "OutputSilo.hh"

namespace Amanzi {

OutputSilo::OutputSilo(Teuchos::ParameterList& plist,
                       const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                       bool is_vis,
                       bool is_dynamic)
  : count_(0), mesh_(mesh), fid_(nullptr)
{
  if (mesh_->vis_mesh().space_dimension() != 3 || mesh_->vis_mesh().manifold_dimension() != 3) {
    Errors::Message msg("OutputSilo is untested on non-3D meshes.");
    Exceptions::amanzi_throw(msg);
  }
  Init_(plist);
}


// destructor must release file resource on non-finalized
OutputSilo::~OutputSilo()
{
  CloseFile_();
}


// open and close files
void
OutputSilo::InitializeCycle(double time, int cycle, const std::string& tag)
{
  // check not open
  AMANZI_ASSERT(fid_ == NULL);
  CloseFile_();
  int comm_size = mesh_->get_comm()->NumProc();
  int comm_rank = mesh_->get_comm()->MyPID();

  // open file
  std::stringstream filename;
  filename << filenamebase_ << "_" << std::setfill('0') << std::setw(sigfigs_) << count_;
  std::string fname = filename.str() + ".silo";

  if (comm_rank == 0) {
    std::cout << "Creating " << fname << " at " << cycle << " time " << time << std::endl;
    fid_ = DBCreate(fname.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_HDF5);
    CloseFile_();
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
      int nnodes = mesh_->vis_mesh().num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_kind::ALL);
      AmanziMesh::Double_List x(nnodes), y(nnodes), z(nnodes);

      AmanziGeometry::Point xyz;
      for (int i = 0; i != nnodes; ++i) {
        mesh_->vis_mesh().node_get_coordinates(i, &xyz);
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
        mesh_->vis_mesh().num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
      ierr |= DBPutUcdmesh(fid_,
                           "mesh",
                           mesh_->vis_mesh().space_dimension(),
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
      int nfaces = mesh_->vis_mesh().num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_kind::ALL);
      std::vector<int> face_node_counts(nfaces);
      std::vector<char> ext_faces(nfaces, 0x0);
      std::vector<int> face_node_list;

      // const auto& nmap = mesh_->vis_mesh().node_map(true);
      for (int f = 0; f != nfaces; ++f) {
        AmanziMesh::Entity_ID_View fnodes;
        mesh_->vis_mesh().face_get_nodes(f, &fnodes);
        //for (int i=0; i!=fnodes.size(); ++i) fnodes[i] = nmap.GID(fnodes[i]);
        face_node_counts[f] = fnodes.size();
        face_node_list.insert(face_node_list.end(), fnodes.begin(), fnodes.end());

        AmanziMesh::Entity_ID_View fcells;
        // mesh_->vis_mesh().face_get_cells(f, AmanziMesh::Parallel_kind::ALL, &fcells);
        mesh_->vis_mesh().face_get_cells(f, AmanziMesh::Parallel_kind::OWNED, &fcells);
        if (fcells.size() == 1) { ext_faces[f] = 0x1; }
      }

      // -- Construct the silo cell-face info
      std::vector<int> cell_face_counts(ncells);
      std::vector<int> cell_face_list;
      // const auto& fmap = mesh_->vis_mesh().face_map(true);

      for (int c = 0; c != ncells; ++c) {
        AmanziMesh::Entity_ID_View cfaces;
        std::vector<int> dirs;
        mesh_->vis_mesh().cell_get_faces_and_dirs(c, &cfaces, &dirs, false);
        for (int i = 0; i != cfaces.size(); ++i) {
          if (dirs[i] < 0) cfaces[i] = ~cfaces[i];
          // cfaces[i] = fmap.GID(cfaces[i]);
        }

        cell_face_counts[c] = cfaces.size();
        cell_face_list.insert(cell_face_list.end(), cfaces.begin(), cfaces.end());
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
      CloseFile_();
    }
    mesh_->get_comm()->Barrier();
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
    CloseFile_();
  }
}

void
OutputSilo::FinalizeCycle()
{
  count_++;
}


// write data to file
void
OutputSilo::WriteVector(const Epetra_Vector& vec,
                        const std::string& name,
                        const AmanziMesh::Entity_kind& kind) const
{
  int comm_size = mesh_->get_comm()->NumProc();
  int comm_rank = mesh_->get_comm()->MyPID();
  auto varname = FixName_(name);

  // open file
  std::stringstream filename;
  filename << filenamebase_ << "_" << std::setfill('0') << std::setw(sigfigs_) << count_;
  std::string fname = filename.str() + ".silo";

  for (int rank = 0; rank != comm_size; ++rank) {
    if (rank == comm_rank) {
      fid_ = DBOpen(fname.c_str(), DB_HDF5, DB_APPEND);

      // directory
      std::string dirname = "/domain_";
      dirname = dirname + std::to_string(rank);
      int ierr = DBSetDir(fid_, dirname.c_str());
      AMANZI_ASSERT(!ierr);

      if (kind == AmanziMesh::CELL) {
        int ierr = DBPutUcdvar1(fid_,
                                varname.c_str(),
                                "mesh",
                                (void*)&vec[0],
                                vec.MyLength(),
                                NULL,
                                0,
                                DB_DOUBLE,
                                DB_ZONECENT,
                                NULL);
        AMANZI_ASSERT(!ierr);
      } else if (kind == AmanziMesh::NODE) {
        int ierr = DBPutUcdvar1(fid_,
                                varname.c_str(),
                                "mesh",
                                (void*)&vec[0],
                                vec.MyLength(),
                                NULL,
                                0,
                                DB_DOUBLE,
                                DB_NODECENT,
                                NULL);
        AMANZI_ASSERT(!ierr);
      } else {
        Errors::Message msg("OutputSilo only knows how to write CELL and NODE based quantities.");
        Exceptions::amanzi_throw(msg);
      }
      CloseFile_();
    }
    mesh_->get_comm()->Barrier();
  }

  if (comm_rank == 0) {
    fid_ = DBOpen(fname.c_str(), DB_HDF5, DB_APPEND);
    DBSetDir(fid_, "/");
    std::string meshname("mesh");
    std::vector<std::string> varnames_str(comm_size);
    std::vector<char*> varnames(comm_size, nullptr);
    std::vector<int> vartypes(comm_size);
    for (int i = 0; i != comm_size; ++i) {
      varnames_str[i] = std::string("/domain_") + std::to_string(i) + "/" + varname;
      varnames[i] = const_cast<char*>(varnames_str[i].c_str());
      vartypes[i] = DB_UCDVAR;
    }

    int ierrl = DBPutMultivar(fid_,
                              FixName_(name).c_str(),
                              comm_size,
                              (char**)varnames.data(),
                              (int*)vartypes.data(),
                              nullptr);
    AMANZI_ASSERT(!ierrl);
    CloseFile_();
  }
}


void
OutputSilo::WriteMultiVector(const Epetra_MultiVector& vec,
                             const std::vector<std::string>& names,
                             const AmanziMesh::Entity_kind& kind) const
{
  AMANZI_ASSERT(vec.NumVectors() == names.size());
  for (int i = 0; i != vec.NumVectors(); ++i) { WriteVector(*vec(i), names[i], kind); }
}

// can we template this?
void
OutputSilo::WriteAttribute(const double& val, const std::string& name) const
{}

void
OutputSilo::WriteAttribute(const int& val, const std::string& name) const
{}


void
OutputSilo::WriteAttribute(const std::string& val, const std::string& name) const
{}


void
OutputSilo::Init_(Teuchos::ParameterList& plist)
{
  filenamebase_ = plist.get<std::string>("file name base", "amanzi_vis");
  sigfigs_ = plist.get<int>("file name counter digits", 5);
}


void
OutputSilo::CloseFile_() const
{
  if (fid_) {
    DBClose(fid_);
    fid_ = nullptr;
  }
}

std::string
OutputSilo::FixName_(const std::string& s) const
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

void
OutputSilo::ReadThrowsError_() const
{
  Errors::Message msg("OutputSilo does not yet support read -- only write.");
  Exceptions::amanzi_throw(msg);
}

} // namespace Amanzi
