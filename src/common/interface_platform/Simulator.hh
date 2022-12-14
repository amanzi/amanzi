/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Jeffrey Johnson (jnjohnson@lbl.gov)
      Markus Berndt (berndt@lanl.gov)
*/

/*
  Interface Platform

*/

#ifndef AMANZI_SIMULATOR_HH_
#define AMANZI_SIMULATOR_HH_

#include <mpi.h>
#include <map>
#include <vector>

// TPLs
#include "xercesc/dom/DOM.hpp"

#include "AmanziComm.hh"
#include "ObservationData.hh"

namespace Amanzi {
class Simulator {
 public:
  enum ReturnType { SUCCESS, FAIL, NUM_RETURN_TYPES };

  Simulator(){};
  explicit Simulator(xercesc::DOMDocument* input){};
  virtual ~Simulator(){};

  virtual Amanzi::Simulator::ReturnType
  Run(const Comm_ptr_type& mpi_comm, Amanzi::ObservationData& output_observations) = 0;

 private:
  // Disallowed operations.
  Simulator(const Simulator&);
  Simulator& operator=(const Simulator&);
};

} // namespace Amanzi

#endif
