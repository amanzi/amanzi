/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include <mpi.h>

#include "Teuchos_GlobalMPISession.hpp"
#include <TestReporterStdout.h>
#include <UnitTest++.h>

#include "VerboseObject_objs.hh"
#include "state_evaluators_registration.hh"

#include "test/factory_models.hh"

template<>
Utils::RegisteredFactory<Evaluator, EvaluatorModel_CompositeVector<ModelCompressiblePorosity>> ModelCompressiblePorosity<cVectorView_type<AmanziDefaultDevice>, VectorView_type<AmanziDefaultDevice>>::fac1_("compressible porosity");

template<>
Utils::RegisteredFactory<Evaluator, EvaluatorModelByMaterial<ModelCompressiblePorosity>> ModelCompressiblePorosity<cVectorView_type<AmanziDefaultDevice>, VectorView_type<AmanziDefaultDevice>>::fac2_("compressible porosity by material");

template<>
Utils::RegisteredFactory<Evaluator, EvaluatorModel_CompositeVector<ModelWaterRetentionVanGenuchten>> ModelWaterRetentionVanGenuchten<cVectorView_type<AmanziDefaultDevice>, VectorView_type<AmanziDefaultDevice>>::fac1_("water retention van Genuchten");

template<>
Utils::RegisteredFactory<Evaluator, EvaluatorModelByMaterial<ModelWaterRetentionVanGenuchten>> ModelWaterRetentionVanGenuchten<cVectorView_type<AmanziDefaultDevice>, VectorView_type<AmanziDefaultDevice>>::fac2_("water retention van Genuchten by material");


int
main(int argc, char* argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::initialize();
  auto status = UnitTest::RunAllTests();
  Kokkos::finalize();
  return status;
}
