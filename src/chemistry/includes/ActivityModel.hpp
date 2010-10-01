#ifndef __Activity_hpp__
#define __Activity_hpp__

// Base class for activity calculations

#include <string>
#include <vector>

#include "Species.hpp"
#include "AqueousEquilibriumComplex.hpp"

class Activity {

 public:
  Activity();
  ~Activity();

  void calculateIonicStrength(std::vector<Species> primarySpecies,
                     std::vector<AqueousEquilibriumComplex> secondarySpecies);
  void calculateActivityCoefficients(
                     std::vector<Species> &primarySpecies,
                     std::vector<AqueousEquilibriumComplex> &secondarySpecies);
  double calculateActivityCoefficient(double charge, 
                                      double ion_size_parameter);


  double ionic_strength(void) const { return this->I_; }

  void ionic_strength(double d) { this->I_ = d; }

  void display(void) const;

 protected:
 private:

  double log_to_ln(double d) { return d*2.30258509299; }
  double ln_to_log(double d) { return d*0.434294481904; }

  double I_; // ionic strength

  const double debyeA;
  const double debyeB;
  const double debyeBdot;

};

#endif // __Activity_hpp__

