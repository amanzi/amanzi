/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ben Andre

  Base class for activity calculations
*/

#ifndef AMANZI_CHEMISTRY_ACTIVITY_MODEL_HH_
#define AMANZI_CHEMISTRY_ACTIVITY_MODEL_HH_

#include <vector>
#include <string>

// Amanzi
#include "VerboseObject.hh"

#include "aqueous_equilibrium_complex.hh"
#include "species.hh"

namespace Amanzi {
namespace AmanziChemistry {

class ActivityModel {
 public:
  ActivityModel();
  virtual ~ActivityModel() {};

  struct ActivityModelParameters {
    // if the activity model requires a database put the name here
    std::string database_filename;
    // This the name of the approach for the J's function in the Pitzer model
    std::string pitzer_jfunction;
  };

  virtual void Setup(const ActivityModelParameters& parameters,
                     const std::vector<Species>& primary_species,
                     const std::vector<AqueousEquilibriumComplex>& secondary_species);

  virtual double Evaluate(const Species& species) = 0;

  virtual void EvaluateVector(
      const std::vector<Species>& primary_species,
      const std::vector<AqueousEquilibriumComplex>& secondary_species,
      std::vector<double>* gamma,
      double* actw) = 0;

  void CalculateIonicStrength(
      const std::vector<Species>& primary_species,
      const std::vector<AqueousEquilibriumComplex>& secondary_species);

  void CalculateSumAbsZ(
      const std::vector<Species>& primary_species,
      const std::vector<AqueousEquilibriumComplex>& secondary_species);

  void CalculateSumC(
      const std::vector<Species>& primary_species,
      const std::vector<AqueousEquilibriumComplex>& secondary_species);

  void CalculateActivityCoefficients(
      std::vector<Species>* primary_species,
      std::vector<AqueousEquilibriumComplex>* secondary_species,
      Species* water);

  // access
  double ionic_strength() const { return I_; }

  void name(const std::string& name) { name_ = name; }
  std::string name() { return name_; }

  // i/o
  virtual void Display() const = 0;
  void set_verbosity(Teuchos::Ptr<VerboseObject> vo) { vo_ = vo; } 

 protected:
  void ResizeGamma(int size);

 protected:
  double I_;  // ionic strength
  double Z_;  // sum ( m_i * abs(z_i) )
  double M_;  // sum ( m_i )

  Teuchos::Ptr<VerboseObject> vo_;

 private:
  std::string name_;
  int num_species_;
  std::vector<double> gamma_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi

#endif
