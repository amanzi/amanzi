/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Chemistry

  Base class for sorption isotherms
*/

#ifndef AMANZI_CHEMISTRY_SORPTION_ISOTHERM_HH_
#define AMANZI_CHEMISTRY_SORPTION_ISOTHERM_HH_

#include <string>

#include "Species.hh"

namespace Amanzi {
namespace AmanziChemistry {

class SorptionIsotherm {
 public:
  enum SorptionIsothermType { FREUNDLICH, LANGMUIR, LINEAR };

  SorptionIsotherm(const std::string& name, const SorptionIsothermType type);
  virtual ~SorptionIsotherm(){};

  // returns sorbed concentration
  virtual double Evaluate(const Species& primary_species) = 0;
  virtual double EvaluateDerivative(const Species& primary_species) = 0;

  virtual void Display(const Teuchos::Ptr<VerboseObject> vo) const {};

  virtual const std::vector<double>& GetParameters() = 0;
  virtual void SetParameters(const std::vector<double>& params) = 0;

  std::string name() const { return name_; }
  SorptionIsothermType isotherm_type() const { return isotherm_type_; }

 private:
  std::string name_;
  SorptionIsothermType isotherm_type_;
};

} // namespace AmanziChemistry
} // namespace Amanzi
#endif
