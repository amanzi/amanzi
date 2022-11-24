/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ben Andre
*/

#ifndef AMANZI_CHEMISTRY_ACTIVITY_MODEL_UNIT_HH_
#define AMANZI_CHEMISTRY_ACTIVITY_MODEL_UNIT_HH_

#include "ActivityModel.hh"

namespace Amanzi {
namespace AmanziChemistry {

class Species;

class ActivityModelUnit : public ActivityModel {
 public:
  ActivityModelUnit() : ActivityModel(){};
  ~ActivityModelUnit(){};

  virtual double Evaluate(const Species& species) final;

  virtual void EvaluateVector(const std::vector<Species>& primary_species,
                              const std::vector<AqueousEquilibriumComplex>& secondary_species,
                              std::vector<double>* gamma,
                              double* actw) final;

  virtual void Display() const override;
};

} // namespace AmanziChemistry
} // namespace Amanzi

#endif
