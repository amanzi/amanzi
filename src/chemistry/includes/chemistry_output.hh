/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_CHEMISTRY_IO_HH_
#define AMANZI_CHEMISTRY_CHEMISTRY_IO_HH_

/*
**
**  Base class for all chemisty io. Handles common code to deal with
**  output to a file and/or stdout and verbosity.
**
**  NOTE(bandre): We are creating a global variable, similar to
**  std::cout: ChemistryOutput chem_out; Inorder to be in compliance
**  with the standard, we need to make the global variable a pointer,
**  ChemistryOutput* chem_out, then new/delete the pointer from
**  main. 
**
**  TODO(bandre): Is this thread safe?!
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
ChemistryOutput* chem_out;
}  // end namespace chemistry
}  // end namespace amanzi

int main ()
{
  namespace ac = amanzi::chemistry;

  ac::OutputOptions output_options;
  output_options.use_stdout = true;
  output_options.file_name = "chemistry-unit-test-results.txt";
  output_options.verbosity_levels.push_back(ac::strings::kVerbosityError);
  output_options.verbosity_levels.push_back(ac::strings::kVerbosityWarning);
  output_options.verbosity_levels.push_back(ac::strings::kVerbosityVerbose);
  ac::chem_out = new ac::ChemistryOutput();
  ac::chem_out->Initialize(output_options);

  .... do something interesting....

  delete ac::chem_out;
}

In other header files:
namespace amanzi {
namespace chemistry {
extern ChemistryOutput* chem_out;
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

void SetupDefaultChemistryOutput(void);

class ChemistryOutput {
 public:
  ChemistryOutput();
  virtual ~ChemistryOutput();

  void Initialize(const OutputOptions& options);
  void AddLevel(const std::string& level);
  void RemoveLevel(const std::string& level);

  void AddLevel(const Verbosity& level);
  void RemoveLevel(const Verbosity& level);

  void Write(const Verbosity level, const std::stringstream& data);
  void Write(const Verbosity level, const std::string& data);

  void set_use_stdout(const bool use_stdout) {
    this->use_stdout_ = use_stdout;
  }

  bool use_stdout(void) const {
    return this->use_stdout_;
  }

 protected:
 private:
  void OpenFileStream(const std::string& file_name);
  void CloseFileStream(void);

  const VerbosityMap& verbosity_map(void) const {
    return this->verbosity_map_;
  }

  const VerbosityFlags& verbosity_flags(void) const {
    return this->verbosity_flags_;
  }

  VerbosityMap verbosity_map_;
  VerbosityFlags verbosity_flags_;
  bool use_stdout_;
  std::ofstream* file_stream_;

};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_CHEMISTRY_OUT_HH_
