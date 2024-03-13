/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ben Andre
*/

/*
  Chemistry

  Class for activity calculations based on the Debye-Huckel B-dot equation.
*/

#ifndef AMANZI_CHEMISTRY_ACTIVITY_MODEL_DEBYE_HUCKEL_HH_
#define AMANZI_CHEMISTRY_ACTIVITY_MODEL_DEBYE_HUCKEL_HH_

#include "ActivityModel.hh"

namespace Amanzi {
namespace AmanziChemistry {

class Species;

class ActivityModelDebyeHuckel : public ActivityModel {
 public:
  ActivityModelDebyeHuckel() : ActivityModel(), max_log_gamma_(7.0){};
  ~ActivityModelDebyeHuckel(){};

  virtual double Evaluate(const Species& species) final;

  virtual void EvaluateVector(const std::vector<Species>& primary_species,
                              const std::vector<AqueousEquilibriumComplex>& secondary_species,
                              std::vector<double>* gamma,
                              double* actw) final;

  virtual void Display() const override;

 private:
  static const double debyeA;
  static const double debyeB;
  static const double debyeBdot;

  double max_log_gamma_;
};

} // namespace AmanziChemistry
} // namespace Amanzi

#endif
