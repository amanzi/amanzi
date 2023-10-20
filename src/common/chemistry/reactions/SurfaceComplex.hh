/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*!

The `"surface complexes`" is the list of surface complexation reactions.
Each reaction is defined by the following parameters:

* `"reaction`" [string] is a surface complexation reaction involing the complex site
  and primary species.

* `"charge`" [int] is the charge of the complex.

* `"equilibrium constant`" [double] is the logarithm of the reaction equilibrium coefficeint.

.. code-block:: xml

  <ParameterList name="surface complexes">
    <ParameterList name="(>davis_O)UO2+">
      <Parameter name="charge" type="int" value="1"/>
      <Parameter name="reaction" type="string" value="1.0 >davis_OH  -1.0 H+  1.0 UO2++"/>
      <Parameter name="equilibrium constant" type="double" value="-0.444"/>
    </ParameterList>
  </ParameterList>

A few additional examples.
Each line has three fields: reaction, logarithm of equailibrium coefficient, and charge.

.. code-block:: txt

   >FeOH2+_s     = 1.0 >FeOH_s   1.0 H+                        -7.18   1.0
   >FeO-_w       = 1.0 >FeOH_w  -1.0 H+                         8.82  -1.0
   >FeOHUO2++    = 1.0 >FeOH     1.0 UO2++                     -6.63   2.0
   >SiO-         =-1.0 H+        1.0 >SiOH                      0.0
   >SiOH2+       = 1.0 H+        1.0 >SiOH                      0.0
   >AlOUO2+      = 1.0 >AlOH    -1.0 H+   1.0 UO2++            -3.13   1.0
   >FeOHZn+_s    = 1.0 >FeOH_s  -1.0 H+   1.0 Zn++             -0.66   1.0
   >SiOUO3H3++   = 1.0 >SiOH     1.0 H2O  1.0 UO2++             5.18   2.0
   >UO2++        = 1.0 UO2++    -1.0 Ca++ 1.0 >Ca++            -5.12   0.0
   >SiOUO3H2+    = 1.0 >SiOH     1.0 H2O  -1.0 H+  1.0 UO2++    5.18   1.0
   >SiOUO3H      = 1.0 >SiOH     1.0 H2O  -2.0 H+  1.0 UO2++    5.18   0.0
   >SiOUO3-      = 1.0 >SiOH     1.0 H2O  -3.0 H+  1.0 UO2++   12.35  -1.0
   >SiOUO2(OH)2- = 1.0 >SiOH     2.0 H2O  -3.0 H+  1.0 UO2++   12.35  -1.0
   >FeOHUO3      = 1.0 >FeOH     1.0 H2O  -2.0 H+  1.0 UO2++    3.05   0.0

*/

#ifndef AMANZI_CHEMISTRY_SURFACECOMPLEX_HH_
#define AMANZI_CHEMISTRY_SURFACECOMPLEX_HH_

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"

#include "FunctionTabular.hh"

#include "Species.hh"
#include "SurfaceSite.hh"

namespace Amanzi {
namespace AmanziChemistry {

// forward declarations from chemistry
class MatrixBlock;

class SurfaceComplex {
 public:
  SurfaceComplex(){};
  SurfaceComplex(const std::string& name,
                 int id,
                 const std::vector<Species>& primary_species,
                 const std::vector<SurfaceSite>& surface_sites,
                 const Teuchos::ParameterList& plist);
  ~SurfaceComplex(){};

  // update molalities
  void Update(const std::vector<Species>& primary_species, const SurfaceSite& surface_site);

  // update temperature dependent quantities
  void UpdateTemperatureDependentCoefs(double T);

  // add stoichiometric contribution of complex to total
  void AddContributionToTotal(std::vector<double>* total);

  // add derivative of total with respect to free-ion to dtotal
  void AddContributionToDTotal(const std::vector<Species>& primary_species, MatrixBlock* dtotal);

  void Display(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const;

  std::string name() const { return name_; }

  int free_site_id() const { return free_site_id_; }
  double free_site_stoichiometry() const { return free_site_stoichiometry_; }
  double stoichiometry(int i) const { return stoichiometry_[i]; }

  int ncomp() const { return ncomp_; };
  int species_id(int i) const { return species_ids_[i]; };
  double surface_concentration() const { return surface_concentration_; };

 private:
  std::string name_;
  int id_;
  double charge_;

  double surface_concentration_; // units? ?[mol/m^3 bulk]?

  int ncomp_; // numebr components in reaction
  std::vector<std::string> species_names_;
  std::vector<int> species_ids_;      // ids of primary species in rxn
  std::vector<double> stoichiometry_; // stoich of primary species in rxn

  std::string free_site_name_;
  double free_site_stoichiometry_; // stoichiometry of free site in rxn
  int free_site_id_;

  double h2o_stoichiometry_;

  double lnK_;  // log value of equlibrium constant
  double lnQK_; // store lnQK for derivatives later
  double logK_;
  Teuchos::RCP<FunctionTabular> func_;
};

} // namespace AmanziChemistry
} // namespace Amanzi
#endif
