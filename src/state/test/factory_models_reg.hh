/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! This provides an example of what a reg file looks like for a user-provided
//! model.

#include "factory_models.hh"

// Each instantiable Evaluator must be given a unique name to allow it to be
// included in the evaluator factory.  These names are, like the model,
// user-provided, but must be static.  Therefore we define them in a
// "registration" file.  These must be done separately from the main header
// file because, while that can be included more than once, this definition
// must be included exactly once.  Furthermore, these models may never be
// referenced in the compiled executable (since they are often only ever
// created by the factory).  So this (and all other) registration file must be
// included in the final executable.
template<>
Utils::RegisteredFactory<Evaluator, EvaluatorModel_CompositeVector<ModelCompressiblePorosity>>
EvaluatorModel_CompositeVector<ModelCompressiblePorosity>::fac_("compressible porosity");

template<>
Utils::RegisteredFactory<Evaluator, EvaluatorModelByMaterial<ModelCompressiblePorosity>>
EvaluatorModelByMaterial<ModelCompressiblePorosity>::fac_("compressible porosity by material");

template<>
Utils::RegisteredFactory<Evaluator, EvaluatorModel_CompositeVector<ModelWaterRetentionVanGenuchten>>
EvaluatorModel_CompositeVector<ModelWaterRetentionVanGenuchten>::fac_("water retention van Genuchten");

template<>
Utils::RegisteredFactory<Evaluator, EvaluatorModelByMaterial<ModelWaterRetentionVanGenuchten>>
EvaluatorModelByMaterial<ModelWaterRetentionVanGenuchten>::fac_("water retention van Genuchten by material");

