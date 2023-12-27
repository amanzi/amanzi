/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Input/Output via HDF5 file for checkpoint/restart.
/*
  Input/Output via HDF5 file for checkpoint/restart.
*/

#include "Key.hh"
#include "CompositeVector.hh"
#include "OutputFactory.hh"
#include "InputOutputHDF5.hh"

namespace Amanzi {

InputHDF5::InputHDF5(const Comm_ptr_type& comm, const std::string& filename)
  : file_(std::make_unique<FileHDF5>(comm, filename, FILE_READONLY))
{}

void
InputHDF5::read(const Teuchos::ParameterList& attrs, double& val) const
{
  val = file_->readAttribute<double>(attrs.name(), "/");
}

void
InputHDF5::read(const Teuchos::ParameterList& attrs, int& val) const
{
  val = file_->readAttribute<int>(attrs.name(), "/");
}

void
InputHDF5::read(const Teuchos::ParameterList& attrs, std::string& val) const
{
  val = file_->readAttribute<std::string>(attrs.name(), "/");
}

void
InputHDF5::read(const Teuchos::ParameterList& attrs, Vector_type& vec) const
{
  file_->readVector(attrs.name(), vec);
}

void
InputHDF5::read(const Teuchos::ParameterList& attrs, IntVector_type& vec) const
{
  file_->readVector(attrs.name(), vec);
}

void
InputHDF5::read(const Teuchos::ParameterList& attrs, MultiVector_type& vec) const
{
  // read as a block
  file_->readMultiVector(attrs.name(), vec);
}

void
InputHDF5::read(const Teuchos::ParameterList& attrs, IntMultiVector_type& vec) const
{
  file_->readMultiVector(attrs.name(), vec);
}


OutputHDF5::OutputHDF5(Teuchos::ParameterList& plist, const Comm_ptr_type& comm) : comm_(comm)
{
  single_file_ = plist.get<bool>("single file", true);
  domain_name_ = Keys::cleanName(plist.name());
  formatter_ = OutputFactory::createDirectoryFormatter(plist);
}


void
OutputHDF5::createTimestep(double time, int cycle)
{
  cycle_ = cycle;
  file_ = std::make_unique<FileHDF5>(comm_, getFilename_(cycle), FILE_CREATE);
}

void
OutputHDF5::finalizeTimestep()
{
  file_.reset(); // destroy the file object, which closes files
}

std::string
OutputHDF5::getFilename(int cycle) const
{
  if (single_file_) {
    return formatter_(cycle) + ".h5";
  } else {
    return formatter_(cycle);
  }
}

std::string
OutputHDF5::getFilename_(int cycle) const
{
  std::string fname = getFilename(cycle);
  if (!single_file_) { fname = fname + "/" + domain_name_ + ".h5"; }
  return fname;
}

void
OutputHDF5::write(const Teuchos::ParameterList& attrs, const int& val) const
{
  file_->writeAttribute(attrs.name(), "/", val);
}

void
OutputHDF5::write(const Teuchos::ParameterList& attrs, const double& val) const
{
  file_->writeAttribute(attrs.name(), "/", val);
}

void
OutputHDF5::write(const Teuchos::ParameterList& attrs, const std::string& val) const
{
  file_->writeAttribute(attrs.name(), "/", val);
}

void
OutputHDF5::write(const Teuchos::ParameterList& attrs, const Vector_type& vec) const
{
  file_->writeVector(attrs.name(), vec);
}

void
OutputHDF5::write(const Teuchos::ParameterList& attrs, const IntVector_type& vec) const
{
  file_->writeVector(attrs.name(), vec);
}

// void
// OutputHDF5::write(const Teuchos::ParameterList& attrs, const MultiVector_type& vec) const
// {
//   file_->writeMultiVector(attrs.name(), vec);
// }

// void
// OutputHDF5::write(const Teuchos::ParameterList& attrs, const IntMultiVector_type& vec) const
// {
//   file_->writeMultiVector(attrs.name(), vec);
// }


} // namespace Amanzi
