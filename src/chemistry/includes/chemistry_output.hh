/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_CHEMISTRY_IO_HH_
#define AMANZI_CHEMISTRY_CHEMISTRY_IO_HH_

/*
**
**  Base class for all chemisty io. Handles common code to deal with
**  output to a file and/or stdout and verbosity.
**
**  TODO(bandre): We are creating a global variable, similar to
**  std::cout: ChemistryOutput chem_out; Inorder to be in compliance with
**  the standard, we need to make the global variable a pointer,
**  new/delete the pointer from main. Need to make sure that this is
**  thread safe!
**
**  TODO(bandre): need a way to silence the output for parallel
**  jobs. Something like if the silent flag is set, the rest of the
**  flags are ignored or overridden to zero.
**
**  TODO(bandre): need to add a function that will add and remove
**  output flags, so output for a certain level can be enable only in
**  one function, etc.
**
**  Use:
// create a global ChemistryOutput object in the amanzi::chemisry
// namespace that can be used by an other chemistry object
namespace amanzi {
namespace chemistry {
ChemistryOutput chem_out;
}  // end namespace chemistry
}  // end namespace amanzi

int main ()
{
  namespace ac = amanzi::chemistry;

  ac::VerbosityFlags verbosity_flags;
  ac::VerbosityMap verbosity_map = ac::CreateVerbosityMap();

  std::string verbosity_level = ac::strings::kVerbosityVerbose;
  verbosity_flags.set(verbosity_map.at(verbosity_level));

  verbosity_level = ac::strings::kVerbosityDebugChemistryCoordinator;
  verbosity_flags.set(verbosity_map.at(verbosity_level));

  bool use_stdout = true;
  std::string output_file_name = "test.txt";
  ac::chem_out.Initialize(verbosity_flags, output_file_name, use_stdout);
}

In other header files:
namespace amanzi {
namespace chemistry {
extern ChemistryOutput chem_out;
}  // end namespace chemistry
}  // end namespace amanzi

**
*/

#include <string>
#include <sstream>
#include <vector>
#include <ostream>

#include "chemistry_containers.hh"
#include "chemistry_verbosity.hh"

namespace amanzi {
namespace chemistry {

class ChemistryOutput {
 public:
  ChemistryOutput();
  virtual ~ChemistryOutput();

  void Initialize(const OutputOptions& options);

  void OpenFileStream(const std::string& file_name);
  void CloseFileStream(void);
  void Write(const Verbosity level, const std::stringstream& data);
  void Write(const Verbosity level, const std::string& data);

  void set_verbosity_flags(const VerbosityFlags& verbosity_flags) {
    this->verbosity_flags_ = verbosity_flags;
  }

  const VerbosityFlags& verbosity_flags(void) const {
    return this->verbosity_flags_;
  }

  void set_use_stdout(const bool use_stdout) {
    this->use_stdout_ = use_stdout;
  }

  bool use_stdout(void) const {
    return this->use_stdout_;
  }

 protected:
  VerbosityFlags verbosity_flags_;
  bool use_stdout_;
  std::ofstream* file_stream_;

 private:


};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_CHEMISTRY_OUT_HH_
