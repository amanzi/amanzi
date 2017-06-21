/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $AMANZI_DIR/COPYRIGHT
Author: ??

Effectively stolen from Amanzi, with few modifications.
------------------------------------------------------------------------- */

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"

struct SimulationDriver
  : public Teuchos::VerboseObject<SimulationDriver>
{
  virtual int Run (const MPI_Comm&               mpi_comm,
                   Teuchos::ParameterList&       input_parameter_list);

};
