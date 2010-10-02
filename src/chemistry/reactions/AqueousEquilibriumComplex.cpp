#include "AqueousEquilibriumComplex.hpp"

AqueousEquilibriumComplex::AqueousEquilibriumComplex() 
                          : Species()
{
  
  ncomp_ = 0; // # components in reaction
  species_names_.clear();
  species_ids_.clear();
  stoichiometry_.clear();
  logK_array_.clear();
  h2o_stoich_ = 0.;
  lnK_ = 0.;
  lnQK_ = 0.0;
  logK_ = 0.0;

} // end AqueousEquilibriumComplex() constructor

AqueousEquilibriumComplex::AqueousEquilibriumComplex(std::string s) 
                          : Species()
{
  static_cast<void>(s);
  // string = "name ncomp stoich1 comp1 stoich2 comp2 ... stoichN compN
  //           logK1 logK2 ... logKN a0 charge mol_wt"
} // end AqueousEquilibriumComplex() constructor

AqueousEquilibriumComplex::AqueousEquilibriumComplex(SpeciesName name, 
                            std::vector<SpeciesName>species,
                            std::vector<double>stoichiometries,
                            std::vector<int>species_ids,
                            double h2o_stoich, double charge, double mol_wt,
                            double size, double logK) 
                            : Species()
{

  ncomp_ = (int)species.size();
  // geh - my compiler needs the specific reference to this object
  this->name(name);
  
  for (std::vector<SpeciesName>::const_iterator i = species.begin(); 
       i != species.end(); i++)
    species_names_.push_back(*i); 
  for (std::vector<double>::const_iterator i = stoichiometries.begin();
       i != stoichiometries.end(); i++)
    stoichiometry_.push_back(*i);
  for (std::vector<int>::const_iterator i = species_ids.begin();
       i != species_ids.end(); i++)
    species_ids_.push_back(*i);

  h2o_stoich_ = h2o_stoich;
  // geh - my compiler needs the specific reference to this object
  this->charge(charge);
  gram_molecular_weight(mol_wt);
  ion_size_parameter(size);
  logK_ = logK;
  lnK_ = log_to_ln(logK);

} // end AqueousEquilibriumComplex() constructor

AqueousEquilibriumComplex::~AqueousEquilibriumComplex() 
{
} // end AqueousEquilibriumComplex() destructor

// temporary location for member functions
// ask Ben!!!
void AqueousEquilibriumComplex::update_kludge(const std::vector<Species> primarySpecies) 
{
  double lnQK = -lnK_;
  for (int i = 0; i < ncomp_; i++) {
    lnQK += stoichiometry_[i] * primarySpecies[species_ids_[i]].ln_activity();
  }
  lnQK_ = lnQK;
//  molality_ = std::exp(lnQK) / act_coef_;
  update(std::exp(lnQK) / act_coef_);
  
} // end update()

void AqueousEquilibriumComplex::addContributionToTotal(std::vector<double> &total) 
{
  for (int i = 0; i < ncomp_; i++) {
    total[species_ids_[i]] += stoichiometry_[i] * molality(); 
  }
} // end addContributionToTotal()

void AqueousEquilibriumComplex::addContributionToDTotal(
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
