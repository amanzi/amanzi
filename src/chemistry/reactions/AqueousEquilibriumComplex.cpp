#include "AqueousEquilibriumComplex.hpp"

AqueousEquilibriumComplex::AqueousEquilibriumComplex() 
{
  
  ncomp_ = 0; // # components in reaction
//  species_names_;
//  species_ids;
//  stoichiometry;
//  logK_array;
  lnK_ = 0.;
  logK_ = 0.0;
  
}

AqueousEquilibriumComplex::AqueousEquilibriumComplex(std::string s) 
{
  // string = "name ncomp stoich1 comp1 stoich2 comp2 ... stoichN compN
  //           logK1 logK2 ... logKN a0 charge mol_wt"

}

AqueousEquilibriumComplex::AqueousEquilibriumComplex(SpeciesName name, 
                            std::vector<SpeciesName>species,
                            std::vector<double>stoichiometries,
                            std::vector<int>species_ids,
                            double h2o_stoich, double charge, double mol_wt,
                            double size, double logK) 
{

  ncomp_ = (int)species.size();
  set_name(name);
  
  for (std::vector<SpeciesName>::const_iterator i=species.begin(); 
       i!=species.end(); i++)
    species_names_.push_back(*i); 
  for (std::vector<double>::const_iterator i=stoichiometries.begin();
       i!=stoichiometries.end(); i++)
    stoichiometry_.push_back(*i);
  for (std::vector<int>::const_iterator i=species_ids.begin();
       i!=species_ids.end(); i++)
    species_ids_.push_back(*i);

  h2o_stoich_ = h2o_stoich;
  set_charge(charge);
  set_gram_molecular_weight(mol_wt);
  set_ion_size_parameter(size);
  logK_ = logK;
  lnK_ = log_to_ln(logK);

}

AqueousEquilibriumComplex::~AqueousEquilibriumComplex() 
{
  
}

// temporary location for member functions
void AqueousEquilibriumComplex::update(const std::vector<Species> primarySpecies) 
{

  double lnQK = -lnK_;
  for (int i=0; i<ncomp_; i++) {
    lnQK += stoichiometry_[i]*primarySpecies[species_ids_[i]].get_ln_activity();
  }
  lnQK_ = lnQK;
  molality_ = exp(lnQK)/act_coef_;

}

void AqueousEquilibriumComplex::addContributionToTotal(std::vector<double> &total) 
{

  for (int i=0; i<ncomp_; i++) {
    total[species_ids_[i]] += stoichiometry_[i]*molality_; 
  }

}

void AqueousEquilibriumComplex::addContributionToDTotal(const std::vector<Species> primarySpecies,
                                                        Block *dtotal) 
{
  
  // taking derivative of contribution to residual in row i with respect
  // to species in column j

  // column loop
  for (int j=0; j<ncomp_; j++) {
    int jcomp = species_ids_[j];
    double tempd = stoichiometry_[j]*                              
      exp(lnQK_-primarySpecies[jcomp].get_ln_molality())/
      act_coef_; // here act_coef is from complex
    // row loop
    for (int i=0; i<ncomp_; i++)
      dtotal->addValue(species_ids_[i],jcomp,stoichiometry_[i]*tempd);
  }

}

void AqueousEquilibriumComplex::display(void) const
{
  std::cout << "    " << get_name() << " = ";
  for (int i = 0; i < (int)species_names_.size(); i++) {
    std::cout << stoichiometry_[i] << " " << species_names_[i];
    if (i < (int)species_names_.size() - 1) {
      std::cout << " + ";
    }
  }
  std::cout << std::endl;
  std::cout << "        logK = " << logK_ << std::endl;
  std::cout << "        charge = " << get_charge() << std::endl;
  std::cout << "        mol wt = " << get_gram_molecular_weight() << std::endl;
}
