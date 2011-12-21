/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $AMANZI_DIR/COPYRIGHT
Author: ??

Effectively stolen from Amanzi, with few modifications.
------------------------------------------------------------------------- */

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "ObservationData.H"
#include "Simulator.H"

struct AmanziUnstructuredGridSimulationDriver
  : Amanzi::Simulator, public Teuchos::VerboseObject<AmanziUnstructuredGridSimulationDriver>
{
  virtual ReturnType Run (const MPI_Comm&               mpi_comm,
                          Teuchos::ParameterList&       input_parameter_list,
                          Amanzi::ObservationData&      output_observations);
};
