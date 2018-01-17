/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator for an unchanging MPI Comm rank.

------------------------------------------------------------------------- */

#include "rank_evaluator.hh"

namespace Amanzi {

Utils::RegisteredFactory<Evaluator, RankEvaluator>
    RankEvaluator::fac_("mpi_comm_rank");

} // namespace Amanzi
