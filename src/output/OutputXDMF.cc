/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  Output

  XDMF implementation of an Output object.
*/

#include "HDF5_MPI.hh"
#include "OutputXDMF.hh"

namespace Amanzi {

OutputXDMF::OutputXDMF(Teuchos::ParameterList& plist,
                       const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                       bool is_vis,
                       bool is_dynamic,
                       bool include_io_set)
  : is_vis_(is_vis),
    is_dynamic_(is_dynamic),
    include_io_set_(include_io_set),
    mesh_written_(false),
    mesh_(mesh)
{
  Init_(plist);
}


// open and close files
void
OutputXDMF::InitializeCycle(double time, int cycle, const std::string& tag)
{
  if (is_dynamic_ || (!mesh_written_)) {
    io_->writeMesh(time, cycle);
    mesh_written_ = true;
  }

  io_->createTimestep(time, cycle, tag);
  io_->open_h5file();
}


void
OutputXDMF::FinalizeCycle()
{
  io_->close_h5file();
  io_->endTimestep();
}


// write data to file
void
OutputXDMF::WriteVector(const Epetra_Vector& vec,
                        const std::string& name,
                        const AmanziMesh::Entity_kind& kind) const
{
  if (mesh_->isLogical()) {
    io_->writeNodeDataReal(vec, name);
  } else if (kind == AmanziMesh::Entity_kind::CELL) {
    io_->writeCellDataReal(vec, name);
  } else if (kind == AmanziMesh::Entity_kind::NODE) {
    const auto& map = vec.Map();
    if (map.MaxElementSize() > 1) {
      Epetra_Vector vec_node(mesh_->getMap(AmanziMesh::Entity_kind::NODE, false));
      int nnodes =
        mesh_->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);

      for (int n = 0; n < nnodes; ++n) {
        int g = map.FirstPointInElement(n);
        vec_node[n] = vec[g];
      }
      io_->writeNodeDataReal(vec_node, name);
    } else {
      io_->writeNodeDataReal(vec, name);
    }
  }
}


void
OutputXDMF::WriteMultiVector(const Epetra_MultiVector& vec,
                             const std::vector<std::string>& names,
                             const AmanziMesh::Entity_kind& kind) const
{
  AMANZI_ASSERT(names.size() == vec.NumVectors());
  for (int i = 0; i != vec.NumVectors(); ++i) {
    WriteVector(*vec(i), names[i], kind);
  }
}


void
OutputXDMF::WriteAttribute(const double& val, const std::string& name) const
{
  io_->writeAttrReal(val, name);
}


void
OutputXDMF::WriteAttribute(const int& val, const std::string& name) const
{
  io_->writeAttrInt(val, name);
}


void
OutputXDMF::WriteAttribute(const std::string& val, const std::string& name) const
{
  io_->writeAttrString(val, name);
}


// read data from file
void
OutputXDMF::ReadVector(Epetra_Vector& vec, const std::string& name) const
{
  io_->readData(vec, name);
}


void
OutputXDMF::ReadMultiVector(Epetra_MultiVector& vec, const std::vector<std::string>& names) const
{
  AMANZI_ASSERT(names.size() == vec.NumVectors());
  for (int i = 0; i != vec.NumVectors(); ++i) {
    io_->readData(*vec(i), names[i]);
  }
}


void
OutputXDMF::ReadAttribute(double& val, const std::string& name) const
{
  io_->readAttrReal(val, name);
}


void
OutputXDMF::ReadAttribute(int& val, const std::string& name) const
{
  io_->readAttrInt(val, name);
}


void
OutputXDMF::ReadAttribute(std::string& val, const std::string& name) const
{
  io_->readAttrString(val, name);
}


void
OutputXDMF::Init_(Teuchos::ParameterList& plist)
{
  // create and set up the HDF5_MPI object
  io_ = Teuchos::rcp(new HDF5_MPI(mesh_->getComm(), include_io_set_));
  io_->setTrackXdmf(is_vis_);
  io_->setDynMesh(is_dynamic_);

  std::string filenamebase = plist.get<std::string>("file name base", "amanzi_vis");
  io_->createMeshFile(mesh_, filenamebase + "_mesh");
  io_->createDataFile(filenamebase + "_data");
}


} // namespace Amanzi
