/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  EOS -- purely virtual base class for an EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_EOS_HH_
#define AMANZI_RELATIONS_EOS_HH_

namespace Amanzi {
namespace Relations {

class EOS {

 public:
  virtual ~EOS() {};

  // Virtual methods that form the EOS
  virtual double MassDensity(double T, double p) = 0;
  virtual double DMassDensityDT(double T, double p) = 0;
  virtual double DMassDensityDp(double T, double p) = 0;

  virtual double MolarDensity(double T, double p) = 0;
  virtual double DMolarDensityDT(double T, double p) = 0;
  virtual double DMolarDensityDp(double T, double p) = 0;

  // If molar mass is constant, we can take some shortcuts if we need both
  // molar and mass densities.  MolarMass() is undefined if
  // !IsConstantMolarMass()
  virtual bool IsConstantMolarMass() = 0;
  virtual double MolarMass() = 0;
};

} // namespace
} // namespace

#endif
