/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
      Markus Berndt
*/

//! Checkpoint writes ALL data, but no meshes to files.

/*
  Reads/writes to/from file using a generic Input/Output object.
*/

#ifndef AMANZI_STATE_CHECKPOINT_HH_
#define AMANZI_STATE_CHECKPOINT_HH_

#include "Teuchos_ParameterList.hpp"

#include "IOEvent.hh"
#include "Output.hh"
#include "Input.hh"
#include "Data_Initializers.hh"

namespace Amanzi {

class Checkpoint : public IOEvent {
 public:
  Checkpoint(Teuchos::ParameterList& plist, const Comm_ptr_type& comm,
             bool read = false);

  // start/finish checkpoint writing
  void CreateFile(double time, int cycle);
  void FinalizeFile(bool final = false);

  // user-provided writing
  template <typename T>
  typename std::enable_if<!Output::writes<T>::value>::type
  Write(const Teuchos::ParameterList& attrs, const T& t) const
  {
    Data_Initializers::UserWriteCheckpoint(*this, attrs, t);
  }

  // output-provided writing
  template <typename T>
  typename std::enable_if<Output::writes<T>::value>::type
  Write(const Teuchos::ParameterList& attrs, const T& t) const
  {
    output_->Write(attrs, t);
  }

  // user-provided reading
  template <typename T>
  typename std::enable_if<!Input::reads<T>::value>::type
  Read(const Teuchos::ParameterList& attrs, T& t) const
  {
    Data_Initializers::UserReadCheckpoint(*this, attrs, t);
  }

  // input-provided reading
  template <typename T>
  typename std::enable_if<Input::reads<T>::value>::type
  Read(const Teuchos::ParameterList& attrs, T& t) const
  {
    input_->Read(attrs, t);
  }

 protected:
  Comm_ptr_type comm_;

  std::unique_ptr<Output> output_;
  std::unique_ptr<Input> input_;

  std::string filename_base_;
};


} // namespace Amanzi

#endif
