/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*!

The class for mineral reaction, should be written with the mineral as the reactant:

    Calcite = 1.0 Ca++ + 1.0 HCO3- -1.0 H+

The `"mineral kinetics`" section includes at the moment both thermodynamic and kinetic data.
The list of parameters for each reaction includes

* `"rate model`" [string] is the model name for reaction rate [mol/s]. Available option is *TST*.

* `"equilibrium constant`" [double] is logarithm of the equilibrium constant.

* `"rate constant`" [double] is log10 of the reaction rate constant.

* `"modifiers`" [string] is the list of pairs: species name and value of modyfying exponent,
  so that the string has always even number of words.

* `"gram molecular weight`" [double] is amount of a molecular substance whose weight, 
  in grams, is numerically equal to the molecular weight of that substance. 

* `"reaction`" [string] is the mineral reaction equation.

* `"molar volume`" [double] is the molar volume, [m^3 / mol].

* `"specific surface area`" [double] is the specific mineral surface area.

The reaction rate is the dissolution rate for the mineral, so it is positive 
for dissolution and negative for precipitation.
We assume that the sublist name is the mineral name.
We also assume that the mineral reaction includes only primary species and that
the modifying species are only primary species.
 
.. code-block:: xml

  <ParameterList name="mineral kinetics">
    <ParameterList name="Kaolinite">
      <Parameter name="rate model" type="string" value="TST"/>
      <Parameter name="rate constant" type="double" value="-8.967"/>
      <Parameter name="modifiers" type="string" value="H+  7.77000E-01"/>
      <Parameter name="gram molecular weight" type="double" value="258.16"/>
      <Parameter name="reaction" type="string" value="5.0 H2O  -6.0 H+  2.0 Al+++  2.0 SiO2(aq)"/>
      <Parameter name="equilibrium constant" type="double" value="7.57"/>
      <Parameter name="molar volume" type="double" value="9.952e-05"/>
      <Parameter name="specific surface area" type="double" value="100.0"/>
    </ParameterList>
  </ParameterList>

A few examples of mineral reactions is beflow.
Each line in has five columns: mineral reaction, logarithm of equilibrium constant, 
gram molecular weight [g/mol], molar volume [cm^3/mol],
and specific surface area [cm^2 mineral / cm^3 bulk].

.. code-block:: txt

   Halite       = 1.0 Na+   1.0 Cl-                                         1.58550   58.4425   27.0150  1.0
   Gypsum       = 2.0 H2O   1.0 SO4-2   1.0 Ca++                           -4.581    172.1722   74.21216 1.0
   Calcite      =-1.0 H+    1.0 HCO3-   1.0 Ca++                            1.8487   100.087    36.934   1.0
   Gibbsite     = 3.0 H2O   1.0 Al+++  -3.0 H+                              7.756     78.0036   31.956   1.0
   Schoepite    = 3.0 H2O  -2.0 H+      1.0 UO2++                           4.8333   322.058    66.08    1.0
   Basaluminite =15.0 H2O -10.0 H+      4.0 Al+++    1.0 SO4--             22.2511   464.140   218.934   1.0
   Ferrihydrite = 2.5 H2O   1.0 Fe++   -2.0 H+       0.25 O2(aq)           -3.594    106.869    23.99    1.0
   Jurbanite    = 6.0 H2O   1.0 Al+++  -1.0 H+       1.0 SO4--             -3.23     230.064   126.0     1.0
   Kaolinite    = 5.0 H2O  -6.0 H+      2.0 Al+++    2.0 SiO2(aq)           7.570    258.160    99.520   1.0
   Soddyite     = 4.0 H2O  -4.0 H+      1.0 SiO2(aq) 2.0 UO2++              0.392    668.169   131.27    1.0
   K-Feldspar   = 2.0 H2O   1.0 K+      1.0 Al+++   -4.0 H+  3.0 SiO2(aq)  -0.2753   278.332   108.87    1.0
   Polyhalite   = 2.0 H2O   1.0 Mg++    2.0 Ca++     2.0 K+  4.0 SO4-2    -13.7440   218.1     100.9722  1.0

*/

#ifndef AMANZI_CHEMISTRY_MINERAL_HH_
#define AMANZI_CHEMISTRY_MINERAL_HH_

#include <cmath>
#include <string>
#include <vector>

#include "SecondarySpecies.hh"
#include "Species.hh"

namespace Amanzi {
namespace AmanziChemistry {

// forward declarations from chemistry
class MatrixBlock;

class Mineral : public SecondarySpecies {
 public:
  Mineral();
  Mineral(int id,
          const std::string& name,
          const Teuchos::ParameterList& plist,
          const std::vector<Species>& primary_species);
  ~Mineral(){};

  // update molalities
  void Update(const std::vector<Species>& primary_species, const Species& water_species);

  // add stoichiometric contribution of complex to total
  void AddContributionToTotal(std::vector<double>* total);

  // add derivative of total with respect to free-ion to dtotal
  void AddContributionToDTotal(const std::vector<Species>& primary_species, MatrixBlock* dtotal);

  void Display(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const;

  double Q_over_K() const { return std::exp(lnQK_); };
  double saturation_index() const { return std::log10(Q_over_K()); }; // SI = log10(Q/Keq)

  double specific_surface_area() const { return specific_surface_area_; }
  void set_specific_surface_area(double d) { specific_surface_area_ = d; }

  double molar_volume() const { return molar_volume_; }

  // not supported yet
  void UpdateSpecificSurfaceArea(){};

  double volume_fraction() const { return volume_fraction_; }
  void set_volume_fraction(double d) { volume_fraction_ = d; }

  void UpdateVolumeFraction(double rate, double dt);

 private:
  double molar_volume_;          // [m^3 / moles]
  double specific_surface_area_; // [m^2 mineral / m^3 bulk]
  double volume_fraction_;       // [m^3 mineral / m^3 bulk]
};

} // namespace AmanziChemistry
} // namespace Amanzi
#endif
