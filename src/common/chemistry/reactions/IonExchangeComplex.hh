/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*!

The class for ion exchange complexation reaction

  NaX <===> Na+ + X-

The `"ion exchange complexes`" is the list of ion exchange complexation reactions.
We assume that the reactions are always written as one complex equals one primary 
species plus one exchange site.
Each complex is defined by the following parameters:

* `"reaction`" [string] is the exchange reaction.

* `"equilibrium constant`" [double] is the logarithm of the equilibrium constant.

.. code-block:: xml

  <ParameterList name="ion exchange complexes">
    <ParameterList name="Na+X">
      <Parameter name="reaction" type="string" value="1.0 Na+  1.0 X-"/>
      <Parameter name="equilibrium constant" type="double" value="1.0"/>
    </ParameterList>
    <ParameterList name="Ca++X">
      <Parameter name="reaction" type="string" value="1.0 Ca++  2.0 X-"/>
      <Parameter name="equilibrium constant" type="double" value="0.316228"/>
    </ParameterList>
  </ParameterList>

A few additional examples, reaction equation and the equilibrium coefficient:

.. code-block:: txt

   Al+++X = 1.0 Al+++ 3.0 X-    1.71133
   Ca0.5X = 0.5 Ca++  1.0 X-   -0.99
   H+X    = 1.0 H+    1.0 X-    0.0251189
   Mg++X  = 1.0 Mg++  2.0 X-    0.1666

*/

#ifndef AMANZI_CHEMISTRY_IONEXCHANGECOMPLEX_HH_
#define AMANZI_CHEMISTRY_IONEXCHANGECOMPLEX_HH_

#include <vector>

#include "Species.hh"

namespace Amanzi {
namespace AmanziChemistry {

class IonExchangeComplex {
 public:
  IonExchangeComplex(){};
  IonExchangeComplex(const std::string& complex_name,
                     int complex_id,
                     const Teuchos::ParameterList& plist,
                     const std::vector<Species>& primary_species);
  virtual ~IonExchangeComplex(){};

  void Display(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayReaction(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const;

  // setters and getters
  std::string name() const { return name_; }
  std::string primary_name() const { return primary_name_; }
  std::string site_name() const { return site_name_; }

  int primary_id() const { return primary_id_; };
  double K() const { return K_; };
  double X() const { return X_; };
  double concentration() const { return concentration_; };

  void set_X(double d) { X_ = d; };
  void set_concentration(double d) { concentration_ = d; };

 private:
  std::string name_, primary_name_, site_name_;
  int id_, primary_id_;
  double concentration_, K_, X_;
};

} // namespace AmanziChemistry
} // namespace Amanzi

#endif
