/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_GEOCHEMISTRY_DRIVER_HH_
#define AMANZI_CHEMISTRY_GEOCHEMISTRY_DRIVER_HH_

//
// Main driver class for all geochemistry computations. This is the
// class called by the MPC or stand alone executables. It handles most
// of the I/O, accepts setup options, reaction network / database
// information, REV info, and thermodynamic state info, then sets up
// and calls the appropriate evaluator: MixingCellEvaluator for batch
// chemistry, GlobalImplicitEvaluator for the GIA appreach.
//

#include <vector>
#include <ostream>


namespace amanzi {
namespace chemistry {

class GeochemistryDriver {
 public:
  GeochemistryDriver();
  virtual ~GeochemistryDriver();

  virtual void Display(std::ostream& output) const;

 protected:
 private:

};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_GEOCHEMISTRY_DRIVER_HH_
