/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)
*/
/*!

A user may request periodic dumps of checkpoint data in the `"checkpoint`"
sublist.  The user has no explicit control over the content of these files, but
has the guarantee that the run can be continued perfectly from that file, in
the sense that a simulation started from a checkpoint file will have the
bitwise same solution as if the run was continued without stopping.  Unlike
`"visualization`", there is only one `"checkpoint`" list for all
domains/meshes.

.. _checkpoint-spec:
.. admonition:: checkpoint-spec

   * `"file name base`" ``[string]`` **"checkpoint"**
   * `"file name digits`" ``[int]`` **5**
   * `"single file checkpoint`" ``[bool]`` **true** If true, writes all
     checkpoint to one file.  If false, uses a subdirectory with one file per
     mesh.  false is required if meshes exist on other communicators than
     MPI_COMM_WORLD, but this is toggled automatically if the code detects that
     this is necessary.

   INCLUDES:

   - ``[io-event-spec]`` An :ref:`IOEvent` spec.

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


#ifndef AMANZI_STATE_CHECKPOINT_HH_
#define AMANZI_STATE_CHECKPOINT_HH_

#include <fstream>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"


#include "HDF5_MPI.hh"
#include "IOEvent.hh"
#include "ObservationData.hh"
#include "Key.hh"

#define DEBUG_COMM_HANGS 0

namespace Amanzi {
class State;

class Checkpoint : public Utils::IOEvent {
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

  // public interface for coordinator clients
  template<typename T>
  void Write(const std::string& name, const T& t) const;

  template<typename T>
  bool Read(const std::string& name, T& t) const
  {
    return true;
  }

  void Write(const State& S,
             WriteType write_type = WriteType::STANDARD,
             Amanzi::ObservationData* obs_data = nullptr);

  void CreateFile(int cycle);
  void CreateFinalFile(int cycle);
  void Finalize();

  // i/o
  void WriteVector(const Epetra_MultiVector& vec, const std::vector<std::string>& names) const;
  void WriteObservations(ObservationData* obs_data);

  void ReadAttributes(State& S);

  // access
  void set_filebasename(std::string base) { filebasename_ = base; }

 protected:
  void ReadParameters_();

  std::string filebasename_;
  int filenamedigits_;
  int restart_cycle_;
  bool single_file_;

  std::map<std::string, Teuchos::RCP<Amanzi::HDF5_MPI>> output_;
};


template<>
void Checkpoint::Write<Epetra_Vector>(const std::string& name, const Epetra_Vector& t) const;

template<>
inline void
Checkpoint::Write<double>(const std::string& name, const double& t) const
{
  auto domain = single_file_ ? std::string("domain") : Keys::getDomain(name);
  if (!output_.count(domain) ) domain = "domain";
  output_.at(domain)->writeAttrReal(t, name);
}

template<>
inline void
Checkpoint::Write<int>(const std::string& name, const int& t) const
{
  auto domain = single_file_ ? std::string("domain") : Keys::getDomain(name);
  if (!output_.count(domain) ) domain = "domain";
  output_.at(domain)->writeAttrInt(t, name);
}

template<>
inline bool
Checkpoint::Read<Epetra_Vector>(const std::string& name, Epetra_Vector& t) const
{
  auto domain = single_file_ ? std::string("domain") : Keys::getDomain(name);
  if (!output_.count(domain) ) domain = "domain";

#if DEBUG_COMM_HANGS
  {
    auto comm_world = getDefaultComm();
    std::string fname = std::string("./err.") + std::to_string(comm_world->MyPID());
    std::fstream str(fname, std::fstream::out | std::fstream::app);
    str << "rank (" << comm_world->MyPID() << "/" << comm_world->NumProc() << ") on comm ("
        << output_.at(domain)->Comm()->MyPID() << "/" << output_.at(domain)->Comm()->NumProc()
        << ") reading \"" << name << "\" of type Epetra_Vector from domain \"" << domain
        << "\" in file \"" << output_.at(domain)->H5DataFilename() << std::endl;
  }
#endif

  output_.at(domain)->readData(t, name);
  return true;
}

template<>
inline bool
Checkpoint::Read<double>(const std::string& name, double& t) const
{
  auto domain = single_file_ ? std::string("domain") : Keys::getDomain(name);
  if (!output_.count(domain) ) domain = "domain";

#if DEBUG_COMM_HANGS
  {
    auto comm_world = getDefaultComm();
    std::string fname = std::string("./err.") + std::to_string(comm_world->MyPID());
    std::fstream str(fname, std::fstream::out | std::fstream::app);
    str << "rank (" << comm_world->MyPID() << "/" << comm_world->NumProc() << ") on comm ("
        << output_.at(domain)->Comm()->MyPID() << "/" << output_.at(domain)->Comm()->NumProc()
        << ") reading \"" << name << "\" of type double from domain \"" << domain << "\" in file \""
        << output_.at(domain)->H5DataFilename() << std::endl;
  }
#endif

  output_.at(domain)->readAttrReal(t, name);
  return true; // FIXME
}

template<>
inline bool
Checkpoint::Read<int>(const std::string& name, int& t) const
{
  auto domain = single_file_ ? std::string("domain") : Keys::getDomain(name);
  if (!output_.count(domain) ) domain = "domain";

#if DEBUG_COMM_HANGS
  {
    auto comm_world = getDefaultComm();
    std::string fname = std::string("./err.") + std::to_string(comm_world->MyPID());
    std::fstream str(fname, std::fstream::out | std::fstream::app);
    str << "rank (" << comm_world->MyPID() << "/" << comm_world->NumProc() << ") on comm ("
        << output_.at(domain)->Comm()->MyPID() << "/" << output_.at(domain)->Comm()->NumProc()
        << ") reading \"" << name << "\" of type int from domain \"" << domain << "\" in file \""
        << output_.at(domain)->H5DataFilename() << std::endl;
  }
#endif

  output_.at(domain)->readAttrInt(t, name);
  return true; // FIXME
}

} // namespace Amanzi
#endif
