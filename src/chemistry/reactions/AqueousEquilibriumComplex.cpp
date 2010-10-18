#include <iomanip>

#include "AqueousEquilibriumComplex.hpp"

AqueousEquilibriumComplex::AqueousEquilibriumComplex() 
                          : SecondarySpecies()
{
  species_names_.clear();
  species_ids_.clear();
  stoichiometry_.clear();
  logK_array_.clear();

} // end AqueousEquilibriumComplex() constructor

AqueousEquilibriumComplex::AqueousEquilibriumComplex(std::string s) 
                          : SecondarySpecies()
{
  static_cast<void>(s);
  // string = "name ncomp stoich1 comp1 stoich2 comp2 ... stoichN compN
  //           logK1 logK2 ... logKN a0 charge mol_wt"
} // end AqueousEquilibriumComplex() constructor

AqueousEquilibriumComplex::AqueousEquilibriumComplex(const SpeciesName name, 
                            const SpeciesId id,
                            std::vector<SpeciesName>species,
                            std::vector<double>stoichiometries,
                            std::vector<int>species_ids,
                            const double h2o_stoich, 
                            const double charge, 
                            const double mol_wt,
                            const double size, 
                            const double logK) 
                            : SecondarySpecies(name, id, species, 
                                               stoichiometries, species_ids,
                                               h2o_stoich, charge, mol_wt, 
                                               size, logK)
{

} // end AqueousEquilibriumComplex() constructor

AqueousEquilibriumComplex::~AqueousEquilibriumComplex() 
{
} // end AqueousEquilibriumComplex() destructor

// temporary location for member functions
// ask Ben!!!
void AqueousEquilibriumComplex::Update_kludge(const std::vector<Species> primarySpecies) 
{
  double lnQK = -lnK_;
  for (int i = 0; i < ncomp_; i++) {
    lnQK += stoichiometry_[i] * primarySpecies[species_ids_[i]].ln_activity();
  }
  lnQK_ = lnQK;
//  molality_ = std::exp(lnQK) / act_coef_;
  update(std::exp(lnQK) / act_coef_);
  
} // end update()

void AqueousEquilibriumComplex::AddContributionToTotal(std::vector<double> &total) 
{
  for (int i = 0; i < ncomp_; i++) {
    total[species_ids_[i]] += stoichiometry_[i] * molality(); 
  }
} // end addContributionToTotal()

void AqueousEquilibriumComplex::AddContributionToDTotal(
                                   const std::vector<Species> primarySpecies,
                                   Block *dtotal) 
{
  
  // taking derivative of contribution to residual in row i with respect
  // to species in column j

  // column loop
  for (int j = 0; j < ncomp_; j++) {
    int jcomp = species_ids_[j];
    double tempd = stoichiometry_[j]*                              
      std::exp(lnQK_ - primarySpecies[jcomp].ln_molality()) / 
      act_coef_; // here act_coef is from complex
    // row loop
    for (int i = 0; i < ncomp_; i++)
      dtotal->addValue(species_ids_[i], jcomp, stoichiometry_[i]*tempd);
  }
} // end addContributionToDTotal()

void AqueousEquilibriumComplex::display(void) const
{
  std::cout << "    " << name() << " = ";
  for (int i = 0; i < (int)species_names_.size(); i++) {
    std::cout << stoichiometry_[i] << " " << species_names_[i];
    if (i < (int)species_names_.size() - 1) {
      std::cout << " + ";
    }
  }
  std::cout << std::endl;
  std::cout << "        logK = " << logK_ << std::endl;
  std::cout << "        charge = " << charge() << std::endl;
  std::cout << "        mol wt = " << gram_molecular_weight() << std::endl;
} // end display()

void AqueousEquilibriumComplex::Display(void) const
{
  std::cout << "    " << name() << " = ";
  for (int i = 0; i < (int)species_names_.size(); i++) {
    std::cout << stoichiometry_[i] << " " << species_names_[i];
    if (i < (int)species_names_.size() - 1) {
      std::cout << " + ";
    }
  }
  std::cout << std::endl;
  std::cout << std::setw(40) << " " 
            << std::setw(10) << logK_
            << std::setw(10) << charge()
            << std::setw(10) << gram_molecular_weight()
            << std::endl;
} // end Display()
