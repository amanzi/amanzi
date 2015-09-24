#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <map>
#include <vector>

// TPLs
#include "xercesc/dom/DOM.hpp"

#include "ObservationData.hh"

namespace Amanzi
{
  class Simulator
  {
    enum ReturnType {SUCCESS, FAIL, NUM_RETURN_TYPES};

    // Constructor accepts an XML DOM.
    explicit Simulator(xercesc::DOMDocument* input) {}

    virtual ~Simulator() {}
    
    virtual Amanzi::Simulator::ReturnType Run(const MPI_Comm&               mpi_comm,
                                              Amanzi::ObservationData&      output_observations) = 0;

   private:

    // Disallowed operations.
    Simulator();
    Simulator(const Simulator&);
    Simulator& operator=(const Simulator&);
  };

} // end namespace amanzi

#endif
