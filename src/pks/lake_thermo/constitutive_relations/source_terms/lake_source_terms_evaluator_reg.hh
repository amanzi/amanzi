/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Source term evaluator for enthalpy of mass source.

  Authors: Svetlana Tokareva (tokareva@lanl.gov)
*/

#include "lake_source_terms_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

Utils::RegisteredFactory<FieldEvaluator,LakeThermoSourceEvaluator> LakeThermoSourceEvaluator::factory_("lake thermo source");

} //namespace
} //namespace

