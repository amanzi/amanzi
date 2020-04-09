/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Evaluates a function to provide Dirichlet data on faces.

#ifndef AMANZI_EVALUATOR_BCS_DIRICHLET_HH_
#  define AMANZI_EVALUATOR_BCS_DIRICHLET_HH_

#  include "EvaluatorIndependent.hh"
#  include "BoundaryFunctionFactory.hh"

namespace Amanzi {

class
