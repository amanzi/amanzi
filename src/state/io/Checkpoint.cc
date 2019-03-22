/* -------------------------------------------------------------------------
Amanzi

License:
Author: Markus Berndt
        Ethan Coon (ecoon@lanl.gov)

Checkpointing for state.
------------------------------------------------------------------------- */
#include <iomanip>
#include <iostream>

#include "boost/filesystem.hpp"

#include "Checkpoint.hh"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

namespace Amanzi {

Checkpoint::Checkpoint(Teuchos::ParameterList& plist,
                       const Comm_ptr_type& comm)
    : IOEvent(plist),
      comm_(comm)
{
  // ReadParameters_();

  // // set the line prefix for output
  // this->setLinePrefix("Amanzi::Checkpoint     ");
  // // make sure that the line prefix is printed
  // this->getOStream()->setShowLinePrefix(true);

  // // Read the sublist for verbosity settings.
  // Teuchos::readVerboseObjectSublist(&plist_, this);

  // // Set up the HDF5
  // //checkpoint_output_ = Teuchos::rcp(new HDF5_MPI(comm_));
  // //checkpoint_output_->setTrackXdmf(false);
  // final_ = false;
}

Checkpoint::Checkpoint(const std::string& filename, const Comm_ptr_type& comm)
    : IOEvent(),
      comm_(comm) {
  //checkpoint_output_ = Teuchos::rcp(new HDF5_MPI(comm_, filename));
  //checkpoint_output_->open_h5file();
}

// this constructor makes an object that will not create any output
Checkpoint::Checkpoint()
    : IOEvent(),
      comm_(Teuchos::null) {}

// -----------------------------------------------------------------------------
// Set up control from parameter list.
// -----------------------------------------------------------------------------
void Checkpoint::ReadParameters_() {
  filebasename_ = plist_.get<std::string>("file name base", "checkpoint");
  filenamedigits_ = plist_.get<int>("file name digits", 5);
};

void Checkpoint::CreateFile(const int cycle) {
  // // create the restart file
  // std::stringstream oss;
  // oss.flush();
  // oss << filebasename_;
  // oss.fill('0');
  // oss.width(filenamedigits_);
  // oss << std::right << cycle;
  // //checkpoint_output_->createDataFile(oss.str());
  // //checkpoint_output_->open_h5file();

  // if (final_) {
  //   std::string ch_file = oss.str() + ".h5";

  //   std::stringstream oss_final;
  //   oss_final << filebasename_ << "_final.h5";
  //   std::string ch_final = oss_final.str();

  //   boost::filesystem::remove(ch_final.data());
  //   boost::filesystem::create_hard_link(ch_file.data(), ch_final.data());
  // }
};

void Checkpoint::Finalize() { //checkpoint_output_->close_h5file(); 
}

} // namespace Amanzi
