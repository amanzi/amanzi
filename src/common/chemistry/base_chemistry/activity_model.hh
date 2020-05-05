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

  void CalculateIonicStrength(
      const std::vector<Species>& primarySpecies,
      const std::vector<AqueousEquilibriumComplex>& secondarySpecies);
  void CalculateSumAbsZ(
        const std::vector<Species>& primarySpecies,
        const std::vector<AqueousEquilibriumComplex>& secondarySpecies);
  void CalculateSumC(
          const std::vector<Species>& primarySpecies,
          const std::vector<AqueousEquilibriumComplex>& secondarySpecies);
  void CalculateActivityCoefficients(
      std::vector<Species>* primarySpecies,
      std::vector<AqueousEquilibriumComplex>* secondarySpecies,
      Species* water);
  virtual double Evaluate(const Species& species) = 0;
  virtual void EvaluateVector(
      const std::vector<Species>& primarySpecies,
      const std::vector<AqueousEquilibriumComplex>& secondarySpecies,
      std::vector<double>* gamma,
      double* actw) = 0;

  double ionic_strength(void) const { return this->I_; }
  int num_species(void) const { return num_species_; }

  virtual void Display(void) const = 0;

  void name(const std::string name) { this->name_ = name; }
  std::string name(void) { return this->name_; }

  void set_verbosity(Teuchos::Ptr<VerboseObject> vo) { vo_ = vo; } 
  Teuchos::Ptr<VerboseObject> verbosity() { return vo_; } 

 protected:
  void ResizeGamma(const int size);

  double log_to_ln(double d) { return d * 2.30258509299; }
  double ln_to_log(double d) { return d * 0.434294481904; }

  void ionic_strength(double d) { this->I_ = d; }

  double I_;  // ionic strength
  double Z_;  // sum ( m_i * abs(z_i) )
  double M_;  // sum ( m_i )

  Teuchos::Ptr<VerboseObject> vo_;

 private:
  void set_num_species(const int value) { this->num_species_ = value; }

  std::string name_;
  int num_species_;
  std::vector<double> gamma_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi

#endif
