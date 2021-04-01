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

#include "activity_model.hh"

namespace Amanzi {
namespace AmanziChemistry {

class Species;

class ActivityModelUnit : public ActivityModel {
 public:
  ActivityModelUnit() : ActivityModel() {};
  ~ActivityModelUnit() {};

  double Evaluate(const Species& species);

  void EvaluateVector(const std::vector<Species>& prim, 
                      const std::vector<AqueousEquilibriumComplex>& sec,
                      std::vector<double>* gamma,
                      double* actw);

  void Display() const;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi

#endif
