/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Basic interface of viscosity.
*/

#ifndef AMANZI_EOS_VISCOSITY_BASE_HH_
#define AMANZI_EOS_VISCOSITY_BASE_HH_

namespace Amanzi {
namespace AmanziEOS {

// Equation of State model
class ViscosityBase {
 public:
  ViscosityBase(Teuchos::ParameterList& eos_plist) : eos_plist_(eos_plist) {};
  virtual ~ViscosityBase() {};

  virtual double Viscosity(double T, double p) = 0;
  virtual double DViscosityDT(double T, double p) = 0;
  virtual double DViscosityDp(double T, double p) = 0;

 protected:
  Teuchos::ParameterList eos_plist_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
