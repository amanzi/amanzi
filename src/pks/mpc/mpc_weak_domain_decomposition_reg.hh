/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "mpc_weak_domain_decomposition.hh"

namespace Amanzi {

RegisteredPKFactory<MPCWeakDomainDecomposition> MPCWeakDomainDecomposition::reg_("weakly coupled subdomain PKs");

} // namespace
