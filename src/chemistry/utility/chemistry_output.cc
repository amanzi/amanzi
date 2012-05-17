/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "chemistry_output.hh"

#include <iostream>
#include <fstream>

#include "chemistry_exception.hh"
#include "chemistry_verbosity.hh"

namespace amanzi {
namespace chemistry {

// global chem_out object, new/delete should be called from the driver!
ChemistryOutput* chem_out = NULL;

// default method to create and initialize the chemistry output object
// with some sane defaults. Can be replaced by an driver if desired.
void SetupDefaultChemistryOutput(void) {
  OutputOptions output_options;
  output_options.use_stdout = true;
  output_options.file_name = "chemistry-output.txt";
  output_options.verbosity_levels.push_back(strings::kVerbosityError);
  output_options.verbosity_levels.push_back(strings::kVerbosityWarning);
  output_options.verbosity_levels.push_back(strings::kVerbosityVerbose);

  if (chem_out == NULL) {
    chem_out = new ChemistryOutput();
  }
  chem_out->Initialize(output_options);
}  // end SetupDefaultChemistryOutput()


ChemistryOutput::ChemistryOutput(void)
    : verbosity_map_(),
      verbosity_flags_(),
      use_stdout_(false),
      file_stream_(NULL) {
  verbosity_map_ = CreateVerbosityMap();
  verbosity_flags_.reset();
}  // end ChemistryOutput constructor

ChemistryOutput::~ChemistryOutput() {
  CloseFileStream();
}  // end ChemistryOutput destructor

void ChemistryOutput::Initialize(const OutputOptions& options) {
  verbosity_flags_.reset();
  std::vector<std::string>::const_iterator level;
  for (level = options.verbosity_levels.begin();
       level != options.verbosity_levels.end(); ++level) {
    AddLevel(*level);
  }
  set_use_stdout(options.use_stdout);
  OpenFileStream(options.file_name);
}  // end Initialize

void ChemistryOutput::AddLevel(const std::string& level) {
  // if a valid level was provided by the user, set the flag
  if (verbosity_map().count(level) > 0) {
    verbosity_flags_.set(verbosity_map().at(level), true);
  }
}  // end AddLevel()

void ChemistryOutput::RemoveLevel(const std::string& level) {
  // valid level provided by the user, unset the flag
  if (verbosity_map().count(level) > 0) {
    verbosity_flags_.set(verbosity_map().at(level), false);
  }
}  // end RemoveLevel()

void ChemistryOutput::AddLevel(const Verbosity& level) {
  verbosity_flags_.set(level, true);
}  // end AddLevel()

void ChemistryOutput::RemoveLevel(const Verbosity& level) {
  verbosity_flags_.set(level, false);
}  // end RemoveLevel()

void ChemistryOutput::OpenFileStream(const std::string& file_name) {
  // close the current file if it exists
  CloseFileStream();
  if (file_name.size()) {
    file_stream_ = new std::ofstream;
    file_stream_->open(file_name.c_str());
    if (!(*file_stream_)) {
      std::ostringstream ost;
      ost << ChemistryException::kChemistryError 
          << "ChemistryOutput::Initialize(): failed to open output file: " 
          << file_name << std::endl;
      throw ChemistryInvalidInput(ost.str());
    }
  }
}  // end OpenFileStream()

void ChemistryOutput::CloseFileStream(void) {
  if (file_stream_) {
    if (file_stream_->is_open()) {
      file_stream_->close();
    }
  }
  //delete file_stream_;
}  // end CloseFileStream()

void ChemistryOutput::Write(const Verbosity level, const std::stringstream& data) {
  //Write(level, data.str().c_str());
  Write(level, data.str());
}  // end Write()

void ChemistryOutput::Write(const Verbosity level, const std::string& data) {
  if (!verbosity_flags().test(kSilent)) {
    if (verbosity_flags().test(level)) {
      if (file_stream_) {
        *file_stream_ << data;
      }
      if (use_stdout()) {
        std::cout << data;
      }
    }
  }
}  // end Write()

}  // namespace chemistry
}  // namespace amanzi
