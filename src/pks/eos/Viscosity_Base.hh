/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Basic interface of a Viscosity.
*/

#ifndef AMANZI_EOS_VISCOSITY_BASE_HH_
#define AMANZI_EOS_VISCOSITY_BASE_HH_

namespace Amanzi {
namespace AmanziEOS {

// Equation of State model
class Viscosity_Base {
 public:
  virtual ~Viscosity_Base() {};

  // Virtual methods that form the Viscosity
  virtual double Viscosity(double T) = 0;
  virtual double DViscosityDT(double T) = 0;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
