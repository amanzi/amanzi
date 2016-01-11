#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <mpi.h>
#include <map>
#include <vector>

// TPLs
#include "xercesc/dom/DOM.hpp"

#include "ObservationData.hh"

namespace Amanzi
{
  class Simulator
  {
   public: 

    enum ReturnType {SUCCESS, FAIL, NUM_RETURN_TYPES};

    // Legacy constructor -- move to private section when we get rid of v1.2 spec.
    Simulator() {}

    // Constructor accepts an XML DOM.
    explicit Simulator(xercesc::DOMDocument* input) {}

    virtual ~Simulator() {}
    
    virtual Amanzi::Simulator::ReturnType Run(const MPI_Comm&               mpi_comm,
                                              Amanzi::ObservationData&      output_observations) = 0;

   private:

    // Disallowed operations.
    Simulator(const Simulator&);
    Simulator& operator=(const Simulator&);
  };

} // end namespace amanzi

#endif
