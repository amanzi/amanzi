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

#ifndef AMANZI_INPUT_OUTPUT_HDF5_HH_
#define AMANZI_INPUT_OUTPUT_HDF5_HH_

#include <string>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "AmanziTypes.hh"
#include "Output.hh"
#include "Input.hh"
#include "FileHDF5.hh"


namespace Amanzi {

class InputHDF5 : public Input {
 public:
  InputHDF5(const Comm_ptr_type& comm, const std::string& filename);
  virtual ~InputHDF5() override = default;

  // read data from file
  virtual void
  Read(const Teuchos::ParameterList& attrs, double& val) const override;
  virtual void
  Read(const Teuchos::ParameterList& attrs, int& val) const override;
  virtual void
  Read(const Teuchos::ParameterList& attrs, std::string& val) const override;
  virtual void
  Read(const Teuchos::ParameterList& attrs, Vector_type& vec) const override;
  virtual void
  Read(const Teuchos::ParameterList& attrs, IntVector_type& vec) const override;
  virtual void Read(const Teuchos::ParameterList& attrs,
                    MultiVector_type& vec) const override;
  virtual void Read(const Teuchos::ParameterList& attrs,
                    IntMultiVector_type& vec) const override;
  virtual void Read(const Teuchos::ParameterList& attrs,
                    CompositeVector_<double>& vec) const override;
  virtual void Read(const Teuchos::ParameterList& attrs,
                    CompositeVector_<int>& vec) const override;

 protected:
  std::unique_ptr<FileHDF5> file_;
};


class OutputHDF5 : public Output {
 public:
  OutputHDF5(Teuchos::ParameterList& plist, const Comm_ptr_type& comm);
  virtual ~OutputHDF5() override = default;

  virtual void CreateFile(double time, int cycle) override;
  virtual void FinalizeFile() override;
  virtual std::string Filename() const override;

  virtual void
  Write(const Teuchos::ParameterList& attrs, const int& val) const override;
  virtual void
  Write(const Teuchos::ParameterList& attrs, const double& val) const override;
  virtual void Write(const Teuchos::ParameterList& attrs,
                     const std::string& val) const override;

  virtual void Write(const Teuchos::ParameterList& attrs,
                     const Vector_type& vec) const override;
  virtual void Write(const Teuchos::ParameterList& attrs,
                     const IntVector_type& vec) const override;

  virtual void Write(const Teuchos::ParameterList& attrs,
                     const MultiVector_type& vec) const override;
  virtual void Write(const Teuchos::ParameterList& attrs,
                     const IntMultiVector_type& vec) const override;

  virtual void Write(const Teuchos::ParameterList& attrs,
                     const CompositeVector_<int>& vec) const override;
  virtual void Write(const Teuchos::ParameterList& attrs,
                     const CompositeVector_<double>& vec) const override;

 protected:
  std::string filename_base_;
  int filename_digits_;
  Comm_ptr_type comm_;
  int cycle_;

  std::unique_ptr<FileHDF5> file_;
};


} // namespace Amanzi

#endif
