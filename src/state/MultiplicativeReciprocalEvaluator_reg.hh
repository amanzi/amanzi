/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  A generic evaluator for multiplying and diving a collection of fields.
*/

#include "MultiplicativeReciprocalEvaluator.hh"

namespace Amanzi {

Utils::RegisteredFactory<FieldEvaluator, MultiplicativeReciprocalEvaluator>
  MultiplicativeReciprocalEvaluator::factory_("multiplicative reciprocal");

} // namespace
