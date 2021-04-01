/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ben Andre
*/

#include <cmath>
#include <iostream>

#include "virial_coefficient.hh"

namespace Amanzi {
namespace AmanziChemistry {

VirialCoefficient::VirialCoefficient()
    :  isp1(-1),
       isp2(-1),
       isp3(-1),
       npol(0),
       ifun1(-1),
       ifun2(-1),
       ifun3(-1),
       virial(0.0) {
}


/*!
    @brief UpdateVirial

    @class VirialCoefficient

    @details Update the virial coefficient as a function of the temperature
    and liquid pressure (last wasn't implemented yet)
*/
void VirialCoefficient::UpdateVirial(const double& temp, const double& pressure) {
  for (int i = 0; i < npol; i++) {
    if (i == 0) { virial = pol[i]; }
    else if (i == 1) { virial += pol[i] * temp; }
    else if (i == 2) { virial += pol[i] / temp; }
    else if (i == 3) { virial += pol[i] * std::log10(temp); }
    else if (i == 4) { virial += pol[i] / (temp*temp); }
  }
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
