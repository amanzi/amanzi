/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include "AqueousEquilibriumComplex.hpp"
#include "Mineral.hpp"
#include "ThermoDatabase.hpp"
#include "Beaker.hpp"
#include "Species.hpp"
#include "StringTokenizer.hpp"

ThermoDatabase::ThermoDatabase(void)
    : Beaker(),
      primary_id_(0),
      mineral_id_(0)
{
}  // end ThermoDatabase constructor

ThermoDatabase::~ThermoDatabase(void)
{
}  // end ThermoDatabase destructor

void ThermoDatabase::setup(std::vector<double> &total, 
                           const Beaker::BeakerParameters parameters)
{
  SetParameters(parameters);
  ReadFile(parameters.thermo_database_file);
  this->SetupActivityModel(parameters.activity_model_name);
  this->resize(this->primary_species().size());

  if (static_cast<unsigned int>(this->ncomp()) != total.size()) {
    // initial conditions and database input don't match. Print a
    // helpful message and exit gracefully.
  }

  this->SetupMineralKinetics(parameters.mineral_kinetics_file);
}  // end setup()

/*******************************************************************************
 **
 **  Thermodynamic database file format
 **
 **  all species information is contained on a single semicolon delimited line.
 **
 **  any line which begins with a # or a space is considered a
 **  comment. blank lines are ignored
 **
 **  File sections are indicated by a < character.
 **
 **  Valid section names are: 
 **
 **    "Primary Species" "Aqueous Equilibrium Complexes" "Minerals" 
 **
 **    "Ion Exchange"
 **
 **  For example <Primary Species  <Aqueous Equilibrium Complexes
 **
 *******************************************************************************/
void ThermoDatabase::ReadFile(const std::string file_name)
{
  if (verbosity() > kDebugInputFile) {
    std::cout << "ThermoDatabase::ReadFile()...." << std::endl;
  }

  std::ifstream input(file_name.c_str());
  if (!input) {
    // should be some type of helpful error message and graceful exit here....
  }

  enum LineType { kCommentLine, kPrimarySpeciesLine, 
                  kAqueousEquilibriumComplexLine, kMineralLine,
                  kUnknownLine };
  enum SectionType { kPrimarySpeciesSection, kAqueousEquilibriumComplexSection, 
                     kMineralSection, kUnknownSection };

  std::string kSectionPrimary("<Primary Species");
  std::string kSectionAqueousEquilibriumComplex("<Aqueous Equilibrium Complexes");
  std::string kSectionMineral("<Minerals");

  LineType line_type;
  SectionType current_section;
  bool parsed_primaries = false;
  int count = 0;
  while (!input.eof() && count < 100) {
    count++;
    std::string line;
    getline(input, line);
    if ((line.size() > 0) && (line[line.size() - 1] == '\r')) {
      // getline only searches for \n line ends. windows files use \r\n
      // check for a hanging \r and remove it if it is there
      line.resize(line.size() - 1);
    }
    char first = line[0];
    if (first == '#' || first == '\0' || first == ' ') {
      line_type = kCommentLine;
    } else if (first == '<') {
      if (line == kSectionPrimary) {
        line_type = kPrimarySpeciesLine;
        current_section = kPrimarySpeciesSection;
      } else if (line == kSectionAqueousEquilibriumComplex) {
        line_type = kAqueousEquilibriumComplexLine;
        current_section = kAqueousEquilibriumComplexSection;
      } else if (line == kSectionMineral) {
        line_type = kMineralLine;
        current_section = kMineralSection;
      } else {
        std::cout << "ThermoDatabase::ReadFile(): unknown section string \'"
                  << line << "\'" << std::endl;
        current_section = kUnknownSection;
      }
    } else {
      // assume we are on a data line within current_section....
      if (current_section == kPrimarySpeciesSection) {
        ParsePrimarySpecies(line);
        parsed_primaries = true;
      } else if (current_section == kAqueousEquilibriumComplexSection) {
        if (parsed_primaries) {
          ParseAqueousEquilibriumComplex(line);
        } else {
          // print a helpful message and exit gracefully
          std::cout << "ERROR: ThermoDatabase::ReadFile() : "
                    << "Attempting to parse aqueous equilibrium complexes before "
                    << "primary species have been specified. Please check for "
                    << "additional error messages and verify database file is "
                    << "correct." << std::endl;
        }
      } else if (current_section == kMineralSection) {
        if (parsed_primaries) {
          ParseMineral(line);
        } else {
          // print a helpful message and exit gracefully
          std::cout << "ERROR: ThermoDatabase::ReadFile() : "
                    << "Attempting to parse minerals before "
                    << "primary species have been specified. Please check for "
                    << "additional error messages and verify database file is "
                    << "correct." << std::endl;
        }
      } else {
        // should never be here....
      }
    }

  } // end while
}  // end ReadFile()

/*******************************************************************************
 **
 **  Thermodynamic database file format
 **
 **  <Primary Species  <Aqueous Equilibrium Complexes
 **
 **  Primary Species Fields:
 **
 **  Name ; size parameter ; charge ; gram molecular weight
 **
 *******************************************************************************/
void ThermoDatabase::ParsePrimarySpecies(const std::string data)
{
  if (verbosity() > kDebugInputFile) {
    std::cout << "ThermoDatabase::ParsePrimarySpecies()...." << std::endl;
  }

  std::string semicolon(";");
  std::string space(" ");
  StringTokenizer primary_data(data, semicolon);
  StringTokenizer no_spaces;
  
  // get name
  no_spaces.tokenize(primary_data.at(0), space);
  std::string name(no_spaces.at(0));
  //std::cout << "name: " << name << "  id: " << primary_id_ << std::endl;

  // get size parameter
  no_spaces.tokenize(primary_data.at(1), space);
  double size_parameter(std::atof(no_spaces.at(0).c_str()));
  //std::cout << "size parameter: " << size_parameter << std::endl;
                  
  // get charge
  no_spaces.tokenize(primary_data.at(2), space);
  double charge(std::atof(no_spaces.at(0).c_str()));
  //std::cout << "change: " << charge << std::endl;

  // get gram molecular weight
  no_spaces.tokenize(primary_data.at(3), space);
  double gram_molecular_weight(std::atof(no_spaces.at(0).c_str()));
  //std::cout << "gmw: " << gram_molecular_weight << std::endl;
  
  Species primary(primary_id_++, name, charge, gram_molecular_weight, size_parameter);
  this->addPrimarySpecies(primary);
  if (verbosity() > kDebugInputFile) {
    primary.display();
  }
                  
}  // end ParsePrimarySpecies()

/*******************************************************************************
 **
 **  Thermodynamic database file format
 **
 **  <Aqueous Equilibrium Complexes
 **
 **  Secondary Species Fields:
 **
 **  Name = coeff reactant ... ; log Keq ; size parameter ; charge ; gram molecular weight
 **
 *******************************************************************************/
void ThermoDatabase::ParseAqueousEquilibriumComplex(const std::string data)
{
  if (verbosity() > kDebugInputFile) {
    std::cout << "ThermoDatabase::ParseAqueousEquilibriumComplex()...." << std::endl;
  }
  std::string semicolon(";");
  std::string space(" ");
  StringTokenizer no_spaces;

  StringTokenizer aqueous_eq(data, semicolon);

  std::string name;
  std::vector<SpeciesName> species;
  std::vector<double> stoichiometries;
  std::vector<int> species_ids;
  double h2o_stoich;
  std::string reaction(aqueous_eq.at(0));
  ParseReaction(reaction, &name, &species, &stoichiometries, &species_ids, &h2o_stoich);

  no_spaces.tokenize(aqueous_eq.at(1), space);
  double logKeq(std::atof(no_spaces.at(0).c_str()));

  no_spaces.tokenize(aqueous_eq.at(2), space);
  double size_parameter(std::atof(no_spaces.at(0).c_str()));

  no_spaces.tokenize(aqueous_eq.at(3), space);
  double charge(std::atof(no_spaces.at(0).c_str()));

  no_spaces.tokenize(aqueous_eq.at(4), space);
  double gram_molecular_weight(std::atof(no_spaces.at(0).c_str()));

  AqueousEquilibriumComplex secondary(name,
                                      species,
                                      stoichiometries,
                                      species_ids,
                                      h2o_stoich,
                                      charge, gram_molecular_weight, size_parameter, logKeq);
  this->addAqueousEquilibriumComplex(secondary);
  if (verbosity() > kDebugInputFile) {
    secondary.display();
  }
         
}  // end ParseAqueousEquilibriumComplex()

/*******************************************************************************
 **
 **  Thermodynamic database file format
 **
 **  <Mineral
 **
 **  Secondary Species Fields:
 **
 **  Name = coeff reactant ... ; log Keq ; gram molecular weight ; molar density
 **
 *******************************************************************************/
void ThermoDatabase::ParseMineral(const std::string data)
{
  if (verbosity() > kDebugInputFile) {
    std::cout << "ThermoDatabase::ParseMineral()...." << std::endl;
  }
  std::string semicolon(";");
  std::string space(" ");
  StringTokenizer no_spaces;

  StringTokenizer aqueous_eq(data, semicolon);

  std::string name;
  std::vector<SpeciesName> species;
  std::vector<double> stoichiometries;
  std::vector<int> species_ids;
  double h2o_stoich;
  std::string reaction(aqueous_eq.at(0));
  ParseReaction(reaction, &name, &species, &stoichiometries, &species_ids, &h2o_stoich);

  no_spaces.tokenize(aqueous_eq.at(1), space);
  double logKeq(std::atof(no_spaces.at(0).c_str()));

  no_spaces.tokenize(aqueous_eq.at(2), space);
  double gram_molecular_weight(std::atof(no_spaces.at(0).c_str()));

  no_spaces.tokenize(aqueous_eq.at(3), space);
  double molar_density(std::atof(no_spaces.at(0).c_str()));

  Mineral mineral(name, mineral_id_++,
                  species,
                  stoichiometries,
                  species_ids,
                  h2o_stoich,
                  gram_molecular_weight, logKeq, molar_density);
  this->addMineral(mineral);
  if (verbosity() > kDebugInputFile) {
    mineral.display();
  }
         
}  // end ParseMineral()

/*******************************************************************************
 **
 **  Thermodynamic database file format
 **
 **  reaction
 **
 **  Fields:
 **
 **  SpeciesName = coeff PrimaryName coeff PrimaryName ... 
 **
 *******************************************************************************/
void ThermoDatabase::ParseReaction(const std::string reaction, 
                                   std::string* name,
                                   std::vector<SpeciesName>* species, 
                                   std::vector<double>* stoichiometries, 
                                   std::vector<int>* species_ids, 
                                   double* h2o_stoich)
{
  if (verbosity() > kDebugInputFile) {
    std::cout << "ThermoDatabase::ParseReaction()...." << std::endl;
  }
  std::string equal("=");
  std::string space(" ");
  StringTokenizer no_spaces;
  std::vector<Species> primary_species = this->primary_species();

  StringTokenizer list(reaction, equal);
  no_spaces.tokenize(list.at(0), space);
  *name = no_spaces.at(0);
  //std::cout << "  name: " << *name << std::endl;

  StringTokenizer products(list.at(1), space);

  for (StringTokenizer::iterator s = products.begin();
       s != products.end(); s++) {
    double coeff = std::atof(s->c_str());
    s++;
    std::string primary_name(*s);

    //std::cout << "  " << coeff << " " << primary_name << std::endl;

    if (primary_name == "H2O") {
      *h2o_stoich = coeff;
      //std::cout << "  " << *h2o_stoich << " H2O" << std::endl;
    } else {
      species->push_back(primary_name);
      stoichiometries->push_back(coeff);
      int id = -1;
      for (SpeciesArray::const_iterator primary = primary_species.begin();
           primary != primary_species.end(); primary++) {
        //std::cout << "    " << primary->name() << "  ??  " << primary_name << std::endl;
        if (primary->name() == primary_name) {
          id = primary->identifier();
          primary = primary_species.end() - 1;
        }
      }
      if (id < 0) {
        std::cout << "ThermoDatabase::ParseReaction(): reaction primary species \'" 
                  << primary_name << "\' was not found in the primary species list...."
                  << std::endl;
      } else {
        species_ids->push_back(id);
      }
    }
  }  // end for(s)
}  // end ParseReaction()

