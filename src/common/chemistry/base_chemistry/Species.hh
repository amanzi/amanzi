/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ben Andre
*/

/*!

The `"primary species`" section is a list of primary species, one sublist for a species.
Each sublist is named after the species and contains the following parameters:

* `"ion size parameter`" [double] is an empirical parameter that provides agreement
  between measured activity coefficients and ionic strength. In theory, it is the diameter
  of the hydrated ion.

* `"charge`" [int] is the ion charge. The net charge of an ion is non-zero since the
  total number of electrons is unequal to the total number of protons.

* `"gram molecular weight`" [double] is amount of a molecular substance whose weight,
  in grams, is numerically equal to the molecular weight of that substance.

.. code-block:: xml

  <ParameterList name="thermodynamic database">
    <ParameterList name="primary species">
      <ParameterList name="H+">
        <Parameter name="ion size parameter" type="double" value="9.0"/>
        <Parameter name="charge" type="int" value="1"/>
        <Parameter name="gram molecular weight" type="double" value="1.0079"/>
      </ParameterList>
      <ParameterList name="Ca++">
        <Parameter name="ion size parameter" type="double" value="6.0"/>
        <Parameter name="charge" type="int" value="2"/>
        <Parameter name="gram molecular weight" type="double" value="40.078"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>

Here is the short list of species that could be used for models.
Each line has four data fields: name of a species, ion size parameter, charge, and atomic mass [u].

.. code-block:: txt

   Al+++      9.0     3.0    26.9815
   HCO3-      4.0    -1.0    61.0171
   HPO4--     4.0    -2.0    95.9793
   Cl-        3.0    -1.0    35.4527
   CO2(aq)    3.0     0.0    44.01
   Cs137      2.5     1.0   132.9054
   F-         3.5    -1.0    18.9984
   Fe++       6.0     2.0    55.847
   K+         3.0     1.0    39.0983
   Mg++       8.0     2.0    24.30
   Na+        4.0     1.0    22.9898
   N2(aq)     3.0     0.0    28.0135
   NO3-       3.0    -1.0    62.0049
   O2(aq)     3.0     0.0    31.9988
   Pb_210     1.0     0.0   210.00
   Pu_238     1.0     0.0   238.00
   Ra_226     1.0     0.0   226.00
   SiO2(aq)   3.0     0.0    60.0843
   SO4--      4.0    -2.0    96.0636
   Sr90       5.0     2.0    87.6200
   Tc_99      1.0     0.0    99.00
   Th_230     1.0     0.0   230.00
   Tritium    9.0     0.0     1.01
   U_234      1.0     0.0   234.00
   UO2++      4.5     2.0   270.028
   Zn++       6.0     2.0    65.39

*/

#ifndef AMANZI_CHEMISTRY_SPECIES_HH_
#define AMANZI_CHEMISTRY_SPECIES_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "VerboseObject.hh"

namespace Amanzi {
namespace AmanziChemistry {

class Species {
 public:
  Species();
  Species(int id, const std::string& name, const Teuchos::ParameterList& plist);
  virtual ~Species(){};

  // update(): calculate the new activity coefficient, set the molarity,
  // activity and associated log values. Need to look at different
  // ActivityCoefficient models and determine what the most generic
  // version of this interface will require.
  virtual void update();
  virtual void update(double molality);

  // accessor methods for calculated values
  double molality() const { return molality_; }
  double activity() const { return activity_; }
  double act_coef() const { return act_coef_; }

  // natural log versions
  double ln_molality() const { return ln_molality_; }
  double ln_activity() const { return ln_activity_; }
  double ln_act_coef() const { return ln_act_coef_; }

  // access invariant data
  int identifier() const { return identifier_; }
  double charge() const { return charge_; }
  double gram_molecular_weight() const { return gram_molecular_weight_; }
  double ion_size_parameter() const { return ion_size_parameter_; }
  std::string name() const { return name_; }

  // these should only be used by the activity coefficient model
  void set_act_coef(double d) { act_coef_ = d; }

  void Display(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const;

 protected:
  // these are dangerous, should only be used internally. Use the
  // public update() to ensure that all related data gets updated!
  void molality(double d) { molality_ = d; }
  void activity(double d) { activity_ = d; }
  void ln_molality(double d) { ln_molality_ = d; }
  void ln_activity(double d) { ln_activity_ = d; }

  double molality_, activity_, act_coef_;
  double ln_molality_, ln_activity_, ln_act_coef_;

 private:
  int identifier_;
  double charge_; // why is this a double rather than int...?
  double gram_molecular_weight_;
  double ion_size_parameter_; // angstroms
  std::string name_;
};

typedef std::vector<Species> SpeciesArray;

} // namespace AmanziChemistry
} // namespace Amanzi

#endif
