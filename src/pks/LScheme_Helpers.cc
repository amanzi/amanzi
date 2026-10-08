/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

#include "LScheme_Helpers.hh"

namespace Amanzi {

void
LSchemeDataPK::reset()
{
  num_itrs = 0;
  patience = 0;
  for (int i = 0; i < 3; ++i) {
    ns_itrs[i] = 0;
    seq_error[i] = 0.0;
  }
}


void
LSchemeDataPK::shift()
{
  for (int i = 1; i >= 0; --i) {
    seq_error[i + 1] = seq_error[i];
    ns_itrs[i + 1] = ns_itrs[i];
  }
  seq_error[0] = 0.0;
  ns_itrs[0] = 0;
  num_itrs++;
}


double
LSchemeDataPK::update()
{
  double r(0.0);
  if (num_itrs > 0) {
    r = seq_error[0] / seq_error[1];
    if (r < 0.6) patience--;
    else if (r > 0.95) patience++;

    if (patience == -2 || ns_itrs[0] < 4) {
      safety_factor *= 0.92;
      patience = 0;
    } else if (patience == 2 || ns_itrs[0] > 10) {
      safety_factor *= 1.1;
      patience = 0;
    } 

    if (num_itrs > 5) safety_factor *= 1.1;

    safety_factor = std::min(safety_factor, safety_factor_max);
    safety_factor = std::max(safety_factor, safety_factor_min);
  }
  return r;
}
 

void
LSchemeDataPK::print(std::ostream& os) const
{
  os << std::scientific 
     << "res: " << seq_error[0] << " " << seq_error[1] << " " << seq_error[2]
     << "  inc: " << last_step_increment
     << "  sf: " << safety_factor 
     << "  patience=" << patience 
     << "  nka: " << ns_itrs[0] << " " << ns_itrs[1] << " " << ns_itrs[2] << std::endl;
}

} // namespace Amanzi

