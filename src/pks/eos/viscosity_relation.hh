/*
  This is the EOS component of the ATS and Amanzi codes.
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  Basic interface of a Viscosity.
*/

#ifndef AMANZI_EOS_VISCOSITY_MODEL_HH_
#define AMANZI_EOS_VISCOSITY_MODEL_HH_

namespace Amanzi {
namespace Relations {

// Equation of State model
class ViscosityRelation {
 public:
  virtual ~ViscosityRelation() {};

  // Virtual methods that form the Viscosity
  virtual double Viscosity(double T) = 0;
  virtual double DViscosityDT(double T) = 0;
};

}  // namespace Relations
}  // namespace Amanzi

#endif
