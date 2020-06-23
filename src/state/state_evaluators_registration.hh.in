/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator for an unchanging cell volume.

------------------------------------------------------------------------- */

#include "cell_volume_evaluator.hh"

namespace Amanzi {

Utils::RegisteredFactory<FieldEvaluator,CellVolumeEvaluator>
CellVolumeEvaluator::fac_("cell volume");

} // namespace
/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator for an unchanging MPI Comm rank.

------------------------------------------------------------------------- */

#include "rank_evaluator.hh"

namespace Amanzi {

Utils::RegisteredFactory<FieldEvaluator,RankEvaluator> RankEvaluator::fac_("mpi_comm_rank");

} // namespace
/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator with no dependencies solved for by a PK.

------------------------------------------------------------------------- */

#include "primary_variable_field_evaluator.hh"

namespace Amanzi {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,PrimaryVariableFieldEvaluator> PrimaryVariableFieldEvaluator::fac_("primary variable");

} // namespace
#include "independent_variable_field_evaluator_fromfunction.hh"
namespace Amanzi {

Utils::RegisteredFactory<FieldEvaluator,IndependentVariableFieldEvaluatorFromFunction> IndependentVariableFieldEvaluatorFromFunction::fac_("independent variable");

} // namespace

#include "secondary_variable_field_evaluator_fromfunction.hh"
namespace Amanzi {

Utils::RegisteredFactory<FieldEvaluator,SecondaryVariableFieldEvaluatorFromFunction> SecondaryVariableFieldEvaluatorFromFunction::fac_("secondary variable from function");

} // namespace

#include "independent_variable_field_evaluator_fromfile.hh"
namespace Amanzi {

Utils::RegisteredFactory<FieldEvaluator,IndependentVariableFieldEvaluatorFromFile> IndependentVariableFieldEvaluatorFromFile::fac_("independent variable from file");

} // namespace

#include "constant_variable_field_evaluator.hh"

namespace Amanzi {

Utils::RegisteredFactory<FieldEvaluator, ConstantVariableFieldEvaluator>
    ConstantVariableFieldEvaluator::fac_("constant variable");

} // namespace Amanzi
/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator for an changing cell volume.

------------------------------------------------------------------------- */

#include "deforming_cell_volume_evaluator.hh"

namespace Amanzi {

Utils::RegisteredFactory<FieldEvaluator,DeformingCellVolumeEvaluator>
DeformingCellVolumeEvaluator::fac_("deforming cell volume");

} // namespace
