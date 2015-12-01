/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (version 1) 
           Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#ifndef AMANZI_WATER_RETENTION_MODEL_HH_
#define AMANZI_WATER_RETENTION_MODEL_HH_

#include <string>

namespace Amanzi {
namespace Flow {

class WRM {
 public:
  virtual ~WRM() {};

  virtual double k_relative(double pc) const = 0;
  virtual double saturation(double pc) const = 0;
  virtual double dSdPc(double pc) const = 0;  // derivative of saturation w.r.t. to capillary pressure
  virtual double capillaryPressure(double s) const = 0;
  virtual double residualSaturation() const = 0;
  virtual double dKdPc(double pc) const { return 0.0; }

  const std::string region() const { return region_; };

 protected:
  void set_region(std::string region) { region_ = region; };
  std::string region_;
};

typedef double(WRM::*KRelFn)(double pc) const; 

}  // namespace Flow
}  // namespace Amanzi
  
#endif
  
