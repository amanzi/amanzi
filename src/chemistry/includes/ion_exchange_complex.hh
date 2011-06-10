/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_IONEXCHANGECOMPLEX_HH_
#define AMANZI_CHEMISTRY_IONEXCHANGECOMPLEX_HH_

/*
**  Class for ion exchange complexation reaction
**
**  NaX <===> Na+ + X-
**
*/


#include <vector>

#include "species.hh"
#include "ion_exchange_site.hh"

// forward declarations
class Block;

class IonExchangeComplex : Species {
 public:
  IonExchangeComplex();
  IonExchangeComplex(const SpeciesName name,
                     const SpeciesId complex_id,
                     const SpeciesName primary_name,
                     const double primary_stoichiometry,
                     const SpeciesId primary_id,
                     const SpeciesName exchange_site_name,
                     const double exchange_site_stoichiometry,
                     const SpeciesId exchange_site_id,
                     const double h2o_stoich,
                     const double log_Keq);
  virtual ~IonExchangeComplex();

  // over ride Species.update()!
  virtual void update(void);
  virtual void update(const double molality);

  // update molalities
  virtual void Update(const std::vector<Species>primary_species,
                      const std::vector<IonExchangeSite>exchange_sites);
  // add stoichiometric contribution of complex to total
  virtual void AddContributionToTotal(std::vector<double>* total);
  // add derivative of total with respect to free-ion to dtotal
  virtual void AddContributionToDTotal(const std::vector<Species>& primary_species,
                                       Block* dtotal);

  void display(void) const;
  void Display(void) const;
  void DisplayReaction(void) const;
  void DisplayResultsHeader(void) const;
  void DisplayResults(void) const;

  double complex_stoich_coeff(void) const {
    return complex_stoich_coeff_;
  };
  SpeciesName primary_name(void) const {
    return primary_name_;
  };
  double primary_stoichiometry(void) const {
    return primary_stoichiometry_;
  };
  SpeciesId primary_id(void) const {
    return primary_id_;
  };

  SpeciesName exchange_site_name(void) const {
    return exchange_site_name_;
  };
  double exchange_site_stoichiometry(void) const {
    return exchange_site_stoichiometry_;
  };
  SpeciesId exchange_site_id(void) const {
    return exchange_site_id_;
  };

  double log_Keq(void) const {
    return this->log_Keq_;
  };
  double ln_Keq(void) const {
    return this->ln_Keq_;
  };

 protected:
  void set_ln_QKeq(double d) {
    this->ln_QKeq_ = d;
  };

 private:
  double complex_stoich_coeff_;  // stoichiometric coefficient of this species!
  SpeciesName primary_name_;
  double primary_stoichiometry_;
  SpeciesId primary_id_;

  SpeciesName exchange_site_name_;
  double exchange_site_stoichiometry_;
  SpeciesId exchange_site_id_;

  double h2o_stoich_;               // stoichiometry of water in equation
  double log_Keq_;                  // log10 value of equlibrium constant
  double ln_Keq_;                   // natural log value of equlibrium constant
  double ln_QKeq_;                  // store lnQK for derivatives later
};

#endif  // AMANZI_CHEMISTRY_IONEXCHANGECOMPLEX_HH_
