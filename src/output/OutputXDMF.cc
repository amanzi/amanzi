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
  if (mesh_->is_logical()) {
    io_->writeNodeDataReal(vec, name);
  } else if (kind == AmanziMesh::CELL) {
    io_->writeCellDataReal(vec, name);
  } else if (kind == AmanziMesh::NODE) {
    io_->writeNodeDataReal(vec, name);
  }
}


void
OutputXDMF::WriteMultiVector(const Epetra_MultiVector& vec,
                             const std::vector<std::string>& names,
                             const AmanziMesh::Entity_kind& kind) const
{
  AMANZI_ASSERT(names.size() == vec.NumVectors());
  for (int i = 0; i != vec.NumVectors(); ++i) {
    if (mesh_->is_logical()) {
      io_->writeNodeDataReal(*vec(i), names[i]);
    } else if (kind == AmanziMesh::CELL) {
      io_->writeCellDataReal(*vec(i), names[i]);
    } else if (kind == AmanziMesh::NODE) {
      io_->writeNodeDataReal(*vec(i), names[i]);
    }
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
  for (int i = 0; i != vec.NumVectors(); ++i) { io_->readData(*vec(i), names[i]); }
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
  io_ = Teuchos::rcp(new HDF5_MPI(mesh_->get_comm(), include_io_set_));
  io_->setTrackXdmf(is_vis_);
  io_->setDynMesh(is_dynamic_);

  std::string filenamebase = plist.get<std::string>("file name base", "amanzi_vis");
  io_->createMeshFile(mesh_, filenamebase + "_mesh");
  io_->createDataFile(filenamebase + "_data");
}


} // namespace Amanzi
