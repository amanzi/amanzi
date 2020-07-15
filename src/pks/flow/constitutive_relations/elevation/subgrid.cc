/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  License: see $ATS_DIR/COPYRIGHT
  Authors: Ahmad Jan (jana@ornl.gov)
           Ethan Coon (ecoon@ornl.gov)
*/

#include <cmath>
#include "subgrid.hh"

namespace Amanzi {
namespace Flow {

double
subgrid_VolumetricDepth(double depth, double del_max, double del_ex)
{
  double vol_pd = depth >= del_max ? depth - del_ex :
      std::pow(depth/del_max, 2) * (2*del_max - 3*del_ex)
      + std::pow(depth/del_max,3) * (2*del_ex - del_max);
  
  if (vol_pd < 0) {
    vol_pd = std::pow(depth/del_max, 2) * (0.01*del_max + 0.005*del_ex)
      + std::pow(depth/del_max,3) * (0.15*del_ex + 0.075*del_max);
  }
  
  return vol_pd;
}

double
subgrid_DVolumetricDepth_DDepth(double depth, double del_max, double del_ex)
{
  double dvol_pd = depth >= del_max ? 1 :
      2 * depth/del_max * (2*del_max - 3*del_ex) / del_max
      + 3 * std::pow(depth/del_max,2) * (2*del_ex - del_max) / del_max;

  if (dvol_pd < 0) {
    dvol_pd = 2 * depth/del_max * (0.01*del_max + 0.005*del_ex) / del_max
      + 3 * std::pow(depth/del_max,2) * (0.15*del_ex + 0.075*del_max) / del_max;
  }
  
  return dvol_pd;
}  

} // namespace
} // namespace
