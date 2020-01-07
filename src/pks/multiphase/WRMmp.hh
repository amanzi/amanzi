/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (version 1) 
           Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)

  The base class for water retention models.
*/

#ifndef AMANZI_MULTIPHASE_WATER_RETENTION_MODEL_HH_
#define AMANZI_MULTIPHASE_WATER_RETENTION_MODEL_HH_

#include <string>

namespace Amanzi {
namespace Multiphase {

class WRMmp {
 public:
  WRMmp() {};
  virtual ~WRMmp() {};

  virtual double k_relative(double Sw, std::string phase_name) = 0;
  virtual double capillaryPressure(double s) = 0;
  virtual double dPc_dS(double s) = 0;
  virtual double dKdPc(double pc) { return 0.0; }
  virtual double dKdS(double Sw, std::string phase_name) = 0;
};

}  // namespace Flow
}  // namespace Amanzi
  
#endif
  
