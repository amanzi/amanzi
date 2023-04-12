/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)
*/

//! Manages checkpoint/restart capability.
/*!

A user may request periodic dumps of ATS Checkpoint Data in the
`"checkpoint`" sublist.  The user has no explicit control over the
content of these files, but has the guarantee that the ATS run will be
reproducible (with accuracies determined by machine round errors and
randomness due to execution in a parallel computing environment).
Therefore, output controls for Checkpoint Data are limited to file
name generation and writing frequency, by numerical cycle number.
Unlike `"visualization`", there is only one `"checkpoint`" list for
all domains/meshes.

.. _checkpoint-spec:
.. admonition:: checkpoint-spec

    * `"file name base`" ``[string]`` **"checkpoint"**
    * `"file name digits`" ``[int]`` **5**
    * `"single file checkpoint`" ``[bool]`` **true** If true, writes all
      checkpoint to one file.  If false, uses a subdirectory with one file per
      mesh.  false is required if meshes exist on other communicators than
      MPI_COMM_WORLD, but this is toggled if the code detects that this is
      necessary.

    INCLUDES:
    - ``[io-event-spec]`` An IOEvent_ spec

Example:

.. code-block:: xml

  <ParameterList name="checkpoint">
    <Parameter name="cycles start period stop" type="Array(int)" value="{0, 100, -1}" />
    <Parameter name="cycles" type="Array(int)" value="{999, 1001}" />
    <Parameter name="times start period stop 0" type="Array(double)" value="{0.0, 10.0, 100.0}"/>
    <Parameter name="times start period stop 1" type="Array(double)" value="{100.0, 25.0, -1.0}"/>
    <Parameter name="times" type="Array(double)" value="{101.0, 303.0, 422.0}"/>
  </ParameterList>

In this example, checkpoint files are written when the cycle number is
a multiple of 100, every 10 seconds for the first 100 seconds, and
every 25 seconds thereafter, along with times 101, 303, and 422.  Files will be written in the form: `"checkpoint00000.h5`".

*/

/*

Developer notes:

As per the above description, Checkpoint must handle ALL domains, which may
exist or not exist on any given rank of MPI_COMM_WORLD.  If all domains are
defined on MPI_COMM_WORLD, everything is typically fine, and we can use a
"single file" checkpoint.  If domains are on a subset of MPI_COMM_WORLD, then
writes happen on that COMM, and so we use separate files for each domain.

*/


#ifndef AMANZI_STATE_CHECKPOINT_HH_
#define AMANZI_STATE_CHECKPOINT_HH_

#include <fstream>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"

#include "IOEvent.hh"
#include "Input.hh"
#include "Output.hh"
#include "ObservationData.hh"
#include "Key.hh"

namespace Amanzi {
class State;

class Checkpoint : public IOEvent {
 public:
  enum class WriteType { STANDARD = 0, FINAL, POST_MORTEM };

  // constructor for do-it-yourself
  Checkpoint(bool single_file = true);

  // constructor for writing
  Checkpoint(Teuchos::ParameterList& plist, const State& S);

  // constructor for reading, potentially on multiple communicators from
  // multiple files in this directory.
  Checkpoint(const std::string& file_or_dirname, const State& S);

  // constructor for reading a single field from a specific file
  Checkpoint(const std::string& filename,
             const Comm_ptr_type& comm,
             const std::string& domain = "domain");

  void write(const State& S,
             WriteType write_type = WriteType::STANDARD,
             Amanzi::ObservationData* obs_data = nullptr);

  void read(State& S);

  // NOTE: everything below is part of the public API, but should probably be
  // made protected.
  //
  // write each data type -- typically not called by the user, but by the Data object
  template <typename T>
  void write(const Teuchos::ParameterList& attrs, const T& t) const {
    auto domain = single_file_ ? std::string("domain") : Keys::getDomain(attrs.name());
    if (!output_.count(domain)) domain = "domain";
    return output_.at(domain)->write(attrs, t);
  }

  // write each data type -- typically not called by the user, but by the Data object
  template <typename T>
  bool read(const Teuchos::ParameterList& attrs, T& t) const
  {
    auto domain = single_file_ ? std::string("domain") : Keys::getDomain(attrs.name());
    if (!input_.count(domain)) domain = "domain";
    input_.at(domain)->read(attrs, t);
    return true; // legacy return value ierr, now will simply fail
  }

 protected:
  void readParameters_();
  void createFile_(const int cycle);
  void createFinalFile_(const int cycle);
  void finalizeFile_();
  void writeObservations_(Amanzi::ObservationData* obs_data);

  bool single_file_;
  std::string filenamebase_;

  std::map<std::string, std::unique_ptr<Output>> output_;
  std::map<std::string, std::unique_ptr<Input>> input_;
};


} // namespace Amanzi
#endif
