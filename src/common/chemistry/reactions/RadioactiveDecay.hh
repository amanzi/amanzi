/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*!

The `"radiaoctive decay`" section is the list of decay reactions for aqueous and 
sorbed species. It not deal with decay of solid phase.
Currently, it does not deal with decay of solid phase.

* `"reactant`" [string] is the species name.

* `"product`" [string] is the species name.

* `"half life`" [double] is the half-life time of decay.

.. code-block:: xml

  <ParameterList name="radioactive decay">
    <ParameterList name="complex_0">
      <Parameter name="reactant" type="string" value="A"/>
      <Parameter name="product" type="string" value="B"/>
      <Parameter name="half life" type="double" value="3.9447e+08"/>
    </ParameterList>
    <ParameterList name="complex_1">
      <Parameter name="reactant" type="string" value="B"/>
      <Parameter name="product" type="string" value="C"/>
      <Parameter name="half life" type="double" value="7.8894e+08"/>
    </ParameterList>
    <ParameterList name="complex_2">
      <Parameter name="reactant" type="string" value="C"/>
      <Parameter name="product" type="string" value=""/>
      <Parameter name="half life" type="double" value="1.57788e+08"/>
    </ParameterList>
  </ParameterList>

A few additionale examples:

.. code-block:: txt

   Cs137  -->  1.0 Cs137    half_life 30.2 years
   Pb_210 -->               half_life 22.2 years
   Pu_238 -->  1.0 U_234    half_life 87.7 years
   Ra_226 -->  1.0 Pb_210   half_life 1.6e3 years
   Th_230 -->  1.0 Ra_226   half_life 7.54e4 years
   U_234  -->  1.0 Th_230   half_life 2.45e5 years
   Tc_99  -->               half_life 2.111e5 years
   Sr90   -->  1.0 Sr90     half_life 28.8 years

*/

#ifndef AMANZI_CHEMISTRY_RADIOACTIVE_DECAY_HH_
#define AMANZI_CHEMISTRY_RADIOACTIVE_DECAY_HH_

#include <string>
#include <vector>

#include "Species.hh"

namespace Amanzi {
namespace AmanziChemistry {

// forward declarations from chemistry
class MatrixBlock;

class RadioactiveDecay {
 public:
  RadioactiveDecay();
  RadioactiveDecay(const Teuchos::ParameterList& plist,
                   const std::map<std::string, int>& name_to_id);
  ~RadioactiveDecay(){};

  // update forward and reverse effective reaction rates
  void UpdateRate(const std::vector<double>& total,
                  const std::vector<double>& total_sorbed,
                  const double porosity,
                  const double saturation,
                  const double bulk_volume);
  void AddContributionToResidual(std::vector<double>* residual);
  void AddContributionToJacobian(const MatrixBlock& dtotal,
                                 const MatrixBlock& dtotal_sorbed,
                                 const double porosity,
                                 const double saturation,
                                 const double bulk_volume,
                                 MatrixBlock* J);
  void Display(const Teuchos::Ptr<VerboseObject> vo) const;

  int parent_id() const { return species_ids_.at(0); }
  std::string parent_name() const { return species_names_.at(0); }
  double rate() const { return rate_; }
  double rate_constant() const { return rate_constant_; }

 private:
  void ConvertHalfLifeToRateConstant();
  std::vector<std::string> species_names_;
  std::vector<int> species_ids_;      // ids of primary species in rxn
  std::vector<double> stoichiometry_; // stoich of primary species in rxn

  double rate_constant_; // rate constant [1/sec]
  double half_life_;
  double rate_;
};

} // namespace AmanziChemistry
} // namespace Amanzi

#endif
