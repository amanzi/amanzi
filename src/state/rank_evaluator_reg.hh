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
