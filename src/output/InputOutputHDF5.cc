/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Input/Output via HDF5 file for checkpoint/restart.

/*
  Input/Output via HDF5 file for checkpoint/restart.
*/

#include "InputOutputHDF5.hh"

namespace Amanzi {

InputHDF5::InputHDF5(const Comm_ptr_type& comm, const std::string& filename)
  : file_(std::make_unique<FileHDF5>(comm, filename, FILE_READONLY))
{}

void
InputHDF5::Read(const Teuchos::ParameterList& attrs, double& val) const
{
  val = file_->ReadAttribute<double>(attrs.name(), "/");
}

void
InputHDF5::Read(const Teuchos::ParameterList& attrs, int& val) const
{
  val = file_->ReadAttribute<int>(attrs.name(), "/");
}

void
InputHDF5::Read(const Teuchos::ParameterList& attrs, std::string& val) const
{
  val = file_->ReadAttribute<std::string>(attrs.name(), "/");
}

void
InputHDF5::Read(const Teuchos::ParameterList& attrs, Vector_type& vec) const
{
  file_->ReadVector(attrs.name(), vec);
}

void
InputHDF5::Read(const Teuchos::ParameterList& attrs, IntVector_type& vec) const
{
  file_->ReadVector(attrs.name(), vec);
}

void
InputHDF5::Read(const Teuchos::ParameterList& attrs,
                MultiVector_type& vec) const
{
  std::vector<std::string> fieldnames;
  for (std::size_t i = 0; i != vec.getNumVectors(); ++i)
    fieldnames.push_back(attrs.name() + "." + std::to_string(i));
  file_->ReadMultiVector(fieldnames, vec);
}

void
InputHDF5::Read(const Teuchos::ParameterList& attrs,
                IntMultiVector_type& vec) const
{
  std::vector<std::string> fieldnames;
  for (std::size_t i = 0; i != vec.getNumVectors(); ++i)
    fieldnames.push_back(attrs.name() + "." + std::to_string(i));
  file_->ReadMultiVector(fieldnames, vec);
}

void
InputHDF5::Read(const Teuchos::ParameterList& attrs,
                CompositeVector_<double>& vec) const
{
  for (const auto& compname : vec) {
    std::string vname = attrs.name() + "." + compname;
    Teuchos::ParameterList attrs_comp(attrs);
    attrs_comp.setName(vname);
    Read(attrs_comp, *vec.GetComponent(compname, false));
  }
}

void
InputHDF5::Read(const Teuchos::ParameterList& attrs,
                CompositeVector_<int>& vec) const
{
  for (const auto& compname : vec) {
    std::string vname = attrs.name() + "." + compname;
    Teuchos::ParameterList attrs_comp(attrs);
    attrs_comp.setName(vname);
    Read(attrs_comp, *vec.GetComponent(compname, false));
  }
}


OutputHDF5::OutputHDF5(Teuchos::ParameterList& plist, const Comm_ptr_type& comm)
  : comm_(comm),
    filename_base_(
      plist.template get<std::string>("file name base", "checkpoint")),
    filename_digits_(plist.template get<int>("file name digits", 5))
{}


void
OutputHDF5::CreateFile(double time, int cycle)
{
  cycle_ = cycle;
  file_ = std::make_unique<FileHDF5>(comm_, Filename(), FILE_CREATE);
}

void
OutputHDF5::FinalizeFile()
{
  file_.reset(); // destroy the file object, which closes files
}

std::string
OutputHDF5::Filename() const
{
  std::stringstream oss;
  oss << filename_base_;
  oss.fill('0');
  oss.width(filename_digits_);
  oss << std::right << cycle_ << ".h5";
  return oss.str();
}

void
OutputHDF5::Write(const Teuchos::ParameterList& attrs, const int& val) const
{
  file_->WriteAttribute(attrs.name(), "/", val);
}

void
OutputHDF5::Write(const Teuchos::ParameterList& attrs, const double& val) const
{
  file_->WriteAttribute(attrs.name(), "/", val);
}

void
OutputHDF5::Write(const Teuchos::ParameterList& attrs,
                  const std::string& val) const
{
  file_->WriteAttribute(attrs.name(), "/", val);
}

void
OutputHDF5::Write(const Teuchos::ParameterList& attrs,
                  const Vector_type& vec) const
{
  file_->WriteVector(attrs.name(), vec);
}

void
OutputHDF5::Write(const Teuchos::ParameterList& attrs,
                  const IntVector_type& vec) const
{
  file_->WriteVector(attrs.name(), vec);
}

void
OutputHDF5::Write(const Teuchos::ParameterList& attrs,
                  const MultiVector_type& vec) const
{
  std::vector<std::string> fieldnames;
  for (std::size_t i = 0; i != vec.getNumVectors(); ++i)
    fieldnames.push_back(attrs.name() + "." + std::to_string(i));
  file_->WriteMultiVector(fieldnames, vec);
}

void
OutputHDF5::Write(const Teuchos::ParameterList& attrs,
                  const IntMultiVector_type& vec) const
{
  std::vector<std::string> fieldnames;
  for (std::size_t i = 0; i != vec.getNumVectors(); ++i)
    fieldnames.push_back(attrs.name() + "." + std::to_string(i));
  file_->WriteMultiVector(fieldnames, vec);
}

void
OutputHDF5::Write(const Teuchos::ParameterList& attrs,
                  const CompositeVector_<int>& vec) const
{
  for (const auto& compname : vec) {
    std::string vname = attrs.name() + "." + compname;
    Teuchos::ParameterList attrs_comp(attrs);
    attrs_comp.setName(vname);
    Write(attrs_comp, *vec.GetComponent(compname, false));
  }
}

void
OutputHDF5::Write(const Teuchos::ParameterList& attrs,
                  const CompositeVector_<double>& vec) const
{
  for (const auto& compname : vec) {
    std::string vname = attrs.name() + "." + compname;
    Teuchos::ParameterList attrs_comp(attrs);
    attrs_comp.setName(vname);
    Write(attrs_comp, *vec.GetComponent(compname, false));
  }
}

} // namespace Amanzi
