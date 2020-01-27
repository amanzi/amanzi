/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A DomainSet coupler, couples a bunch of domains of the same structure.

------------------------------------------------------------------------- */

#include "DomainSetMPC.hh"

namespace Amanzi {

// must communicate dts since columns are serial
double DomainSetMPC::get_dt() {
  double dt = 1.0e99;
  for (const auto& pk : sub_pks_) {
    dt = std::min<double>(dt,pk->get_dt());
  }
  
  double dt_local = dt;
  solution_->Comm()->MinAll(&dt_local, &dt, 1);
  
  return dt;
}

// -----------------------------------------------------------------------------
// Set timestep for sub PKs 
// -----------------------------------------------------------------------------
void DomainSetMPC::set_dt( double dt) {
  for (const auto& pk : sub_pks_) {
    pk->set_dt(dt);
  }

};


//-------------------------------------------------------------------------------------
// Semi coupled thermal hydrology
bool 
DomainSetMPC::AdvanceStep(double t_old, double t_new, bool reinit) {
  int nfailed = 0;
  for (const auto& pk : sub_pks_) {
    bool fail = pk->AdvanceStep(t_old, t_new, reinit);
    if (fail) {
      nfailed++;
      break;
    }
  }

  int nfailed_global(0);
  solution_->Comm()->SumAll(&nfailed, &nfailed_global, 1);
  if (nfailed_global) return true;
  return false;
}
  

} // namespace Amanzi
