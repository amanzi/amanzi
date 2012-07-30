/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_ACTIVITY_MODEL_HH_
#define AMANZI_CHEMISTRY_ACTIVITY_MODEL_HH_

// Base class for activity calculations

#include <vector>
#include <string>

#include "species.hh"
#include "aqueous_equilibrium_complex.hh"
#include "chemistry_verbosity.hh"

namespace amanzi {
namespace chemistry {

class ActivityModel {
 public:
  struct ActivityModelParameters {
    Verbosity verbosity;
    // if the activity model requires a database put the name here
    std::string database_filename;
    // This the name of the approach for the J's function in the Pitzer model
    std::string pitzer_jfunction;
  };

  ActivityModel();
  virtual ~ActivityModel();

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

  double ionic_strength(void) const {
    return this->I_;
  }

  virtual void Display(void) const = 0;

  void name(const std::string name) {
    this->name_ = name;
  }
  std::string name(void) {
    return this->name_;
  }

  virtual void set_verbosity(const Verbosity verbosity) {
    this->verbosity_ = verbosity;
  };
  virtual Verbosity verbosity(void) const {
    return this->verbosity_;
  };

 protected:
  double log_to_ln(double d) {
    return d * 2.30258509299;
  }
  double ln_to_log(double d) {
    return d * 0.434294481904;
  }

  void ionic_strength(double d) {
    this->I_ = d;
  }

  double I_;  // ionic strength

  double Z_;  // sum ( m_i * abs(z_i) )

  double M_;  // sum ( m_i )

 private:
  Verbosity verbosity_;
  std::string name_;
};

}  // namespace chemistry
}  // namespace amanzi

#endif  // AMANZI_CHEMISTRY_ACTIVITY_MODEL_HH_
