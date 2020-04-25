/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
*/
 
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include "aqueous_equilibrium_complex.hh"
#include "general_rxn.hh"
#include "radioactive_decay.hh"
#include "mineral_kinetics_factory.hh"
#include "mineral.hh"
#include "surface_site.hh"
#include "surface_complex.hh"
#include "ion_exchange_rxn.hh"
#include "ion_exchange_site.hh"
#include "ion_exchange_complex.hh"
#include "sorption_isotherm_factory.hh"
#include "sorption_isotherm.hh"
#include "sorption_isotherm_linear.hh"
#include "beaker.hh"
#include "species.hh"
#include "string_tokenizer.hh"
#include "chemistry_exception.hh"

#include "exceptions.hh"
#include "simple_thermo_database.hh"

namespace Amanzi {
namespace AmanziChemistry {

/*
** TODO(bandre): need a lot more error checking and helpfull error messages from here.
*/
SimpleThermoDatabase::SimpleThermoDatabase(Teuchos::RCP<VerboseObject> vo)
    : Beaker(),
      primary_id_(0),
      aqueous_equilibrium_complex_id_(0),
      mineral_id_(0),
      ion_exchange_complex_id_(0),
      surface_site_id_(0),
      surface_complex_id_(0),
      surface_complexation_rxn_id_(0) {
  surface_sites_.clear();
  surface_complexation_reactions_.clear();
  vo_ = vo;
}


void SimpleThermoDatabase::Setup(const Beaker::BeakerComponents& components,
                                 const Beaker::BeakerParameters& parameters) {
  ReadFile(parameters.thermo_database_file);
  SetParameters(parameters);
  SetupActivityModel(parameters.activity_model_name, 
                     parameters.pitzer_database, parameters.jfunction_pitzer);
  ResizeInternalMemory(primary_species().size());
  VerifyComponentSizes(components);
  CopyComponentsToBeaker(components);
}


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
void SimpleThermoDatabase::ReadFile(const std::string& file_name) {
  if (verbosity() == kDebugInputFile) {
    std::cout << "SimpleThermoDatabase::ReadFile()...." << std::endl;
  }

  std::ifstream input(file_name.c_str());
  if (!input) {
    std::stringstream msg;
    msg << "SimpleThermoDatabase cannot open file \"" << file_name << "\"\n";
    vo_->WriteWarning(Teuchos::VERB_NONE, msg);
    Exceptions::amanzi_throw(ChemistryInvalidInput(msg.str()));
  }

  enum LineType { kCommentLine, kPrimarySpeciesLine,
                  kAqueousEquilibriumComplexLine,
                  kMineralLine, kMineralKineticsLine,
                  kGeneralKineticsLine, kRadioactiveDecayLine,
                  kIonExchangeSiteLine, kIonExchangeComplexLine,
                  kSurfaceComplexSiteLine, kSurfaceComplexLine,
                  kIsothermLine,
                  kUnknownLine
  };
  enum SectionType { kPrimarySpeciesSection, kAqueousEquilibriumComplexSection,
                     kMineralSection, kMineralKineticsSection,
                     kGeneralKineticsSection, kRadioactiveDecaySection,
                     kIonExchangeSiteSection, kIonExchangeComplexSection,
                     kSurfaceComplexSiteSection, kSurfaceComplexSection,
                     kIsothermSection,
                     kUnknownSection
  };

  std::string kSectionPrimary("<Primary Species");
  std::string kSectionAqueousEquilibriumComplex("<Aqueous Equilibrium Complexes");
  std::string kSectionMineral("<Minerals");
  std::string kSectionMineralKinetics("<Mineral Kinetics");
  std::string kSectionGeneralKinetics("<General Kinetics");
  std::string kSectionRadioactiveDecay("<Radioactive Decay");
  std::string kSectionIonExchangeSites("<Ion Exchange Sites");
  std::string kSectionIonExchangeComplexes("<Ion Exchange Complexes");
  std::string kSectionSurfaceComplexSites("<Surface Complex Sites");
  std::string kSectionSurfaceComplexes("<Surface Complexes");
  std::string kSectionIsotherms("<Isotherms");

  LineType line_type;
  SectionType current_section;
  bool parsed_primaries = false;
  bool parsed_minerals = false;
  bool parsed_ion_exchange_sites = false;
  bool parsed_surface_complex_sites = false;
  int count = 0;
  while (!input.eof() && count < 5000) {
    count++;
    int data_order_error = 0;
    std::string error_section("");
    std::string line;
    getline(input, line);

    if ((line.size() > 0) && (line.at(line.size() - 1) == '\r')) {
      // getline only searches for \n line ends. windows files use \r\n
      // check for a hanging \r and remove it if it is there
      line.resize(line.size() - 1);
    }
    char first = '\0';
    if (line.length() > 0) first = line[0];
    if (first == '#' || first == '\0') {
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
      } else if (line == kSectionMineralKinetics) {
        line_type = kMineralKineticsLine;
        current_section = kMineralKineticsSection;
      } else if (line == kSectionGeneralKinetics) {
        line_type = kGeneralKineticsLine;
        current_section = kGeneralKineticsSection;
      } else if (line == kSectionRadioactiveDecay) {
        line_type = kRadioactiveDecayLine;
        current_section = kRadioactiveDecaySection;
      } else if (line == kSectionIonExchangeSites) {
        line_type = kIonExchangeSiteLine;
        current_section = kIonExchangeSiteSection;
      } else if (line == kSectionIonExchangeComplexes) {
        line_type = kIonExchangeComplexLine;
        current_section = kIonExchangeComplexSection;

      } else if (line == kSectionSurfaceComplexSites) {
        line_type = kSurfaceComplexSiteLine;
        current_section = kSurfaceComplexSiteSection;
      } else if (line == kSectionSurfaceComplexes) {
        line_type = kSurfaceComplexLine;
        current_section = kSurfaceComplexSection;

      } else if (line == kSectionIsotherms) {
        line_type = kIsothermLine;
        current_section = kIsothermSection;
      } else {
        std::cout << "SimpleThermoDatabase::ReadFile(): unknown section string \'"
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
          data_order_error = 1;
          error_section = "aqueous equilibrium complexes";
        }
      } else if (current_section == kGeneralKineticsSection) {
        if (parsed_primaries) {
          ParseGeneralKinetics(line);
        } else {
          data_order_error = 1;
          error_section = "general kinetics";
        }
      } else if (current_section == kRadioactiveDecaySection) {
        if (parsed_primaries) {
          ParseRadioactiveDecay(line);
        } else {
          data_order_error = 1;
          error_section = "radioactive decay";
        }
      } else if (current_section == kMineralSection) {
        if (parsed_primaries) {
          ParseMineral(line);
          parsed_minerals = true;
        } else {
          data_order_error = 1;
          error_section = "minerals";
        }
      } else if (current_section == kMineralKineticsSection) {
        if (parsed_primaries && parsed_minerals) {
          ParseMineralKinetics(line);
        } else {
          data_order_error = 2;
          error_section = "mineral kinetics";
        }
      } else if (current_section == kIonExchangeSiteSection) {
        ParseIonExchangeSite(line);
        parsed_ion_exchange_sites = true;
      } else if (current_section == kIonExchangeComplexSection) {
        if (parsed_primaries && parsed_ion_exchange_sites) {
          ParseIonExchangeComplex(line);
        } else {
          error_section = "ion exchange complex";
          if (parsed_primaries) {
            if (!parsed_ion_exchange_sites) {
              data_order_error = 3;
            }
          } else if (parsed_ion_exchange_sites) {
            if (!parsed_primaries) {
              data_order_error = 1;
            }
          } else {
            data_order_error = 4;
          }
        }
      } else if (current_section == kSurfaceComplexSiteSection) {
        ParseSurfaceComplexSite(line);
        parsed_surface_complex_sites = true;
      } else if (current_section == kSurfaceComplexSection) {
        if (parsed_primaries && parsed_surface_complex_sites) {
          ParseSurfaceComplex(line);
        } else {
          error_section = "surface complex";
          if (parsed_primaries) {
            if (!parsed_surface_complex_sites) {
              data_order_error = 5;
            }
          } else if (parsed_surface_complex_sites) {
            if (!parsed_primaries) {
              data_order_error = 1;
            }
          } else {
            data_order_error = 4;
          }
        }
      } else if (current_section == kIsothermSection) {
        if (parsed_primaries) {
          ParseSorptionIsotherm(line);
        } else {
          error_section = "isotherms";
          data_order_error = 1;
        }
      } else {
        // should never be here....
      }
      if (data_order_error) {
        // print a helpful message and exit gracefully
        std::ostringstream error_stream;
        error_stream << "SimpleThermoDatabase::ReadFile() : "
                     << "Attempting to parse " << error_section << " before ";
        std::string temp = "";
        if (data_order_error == 5) {
          temp += "surface complex sites ";
        } else if (data_order_error == 3) {
          temp += "ion exchange sites ";
        } else if (data_order_error == 2) {
          temp += "minerals ";
        } else if (data_order_error == 1) {
          temp += "primary species ";
        }
        error_stream << temp << "have been specified. Please check for "
                     << "additional error messages and verify database file is "
                     << "correct." << std::endl;
        error_stream << "SimpleThermoDatabase::ReadFile(): \n";
        error_stream << "a data order error has occured reading file " << file_name
                     << "       please see output for details.\n";
        Exceptions::amanzi_throw(ChemistryInvalidInput(error_stream.str()));
      }
    }
  }  // end while

  FinishSurfaceComplexation();
}


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
void SimpleThermoDatabase::ParsePrimarySpecies(const std::string& data) {
  if (verbosity() == kDebugInputFile) {
    std::cout << "SimpleThermoDatabase::ParsePrimarySpecies()...." << std::endl;
    std::cout << "  data: " << data << std::endl;
  }

  std::string semicolon(";");
  std::string space(" \t");
  StringTokenizer primary_data(data, semicolon);
  StringTokenizer no_spaces;

  // get name
  no_spaces.tokenize(primary_data.at(0), space);
  std::string name(no_spaces.at(0));
  // std::cout << "name: " << name << "  id: " << primary_id_ << std::endl;

  // get size parameter
  no_spaces.tokenize(primary_data.at(1), space);
  double size_parameter(std::atof(no_spaces.at(0).c_str()));
  // std::cout << "size parameter: " << size_parameter << std::endl;

  // get charge
  no_spaces.tokenize(primary_data.at(2), space);
  double charge(std::atof(no_spaces.at(0).c_str()));
  // std::cout << "change: " << charge << std::endl;

  // get gram molecular weight
  no_spaces.tokenize(primary_data.at(3), space);
  double gram_molecular_weight(std::atof(no_spaces.at(0).c_str()));
  // std::cout << "gmw: " << gram_molecular_weight << std::endl;

  Species primary(primary_id_++, name, charge, gram_molecular_weight, size_parameter);
  this->AddPrimarySpecies(primary);
  if (verbosity() == kDebugInputFile) {
    primary.display();
  }
}


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
void SimpleThermoDatabase::ParseAqueousEquilibriumComplex(const std::string& data) {
  if (verbosity() == kDebugInputFile) {
    std::cout << "SimpleThermoDatabase::ParseAqueousEquilibriumComplex()...." << std::endl;
    std::cout << "  data: " << data << std::endl;
  }
  std::string semicolon(";");
  std::string space(" \t");
  StringTokenizer no_spaces;

  StringTokenizer aqueous_eq(data, semicolon);

  std::string name;
  std::vector<SpeciesName> species;
  std::vector<double> stoichiometries;
  std::vector<SpeciesId> species_ids;
  double h2o_stoich = 0;
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
                                      aqueous_equilibrium_complex_id_++,
                                      species,
                                      stoichiometries,
                                      species_ids,
                                      h2o_stoich,
                                      charge, gram_molecular_weight, size_parameter, logKeq);
  this->AddAqueousEquilibriumComplex(secondary);
  if (verbosity() == kDebugInputFile) {
    secondary.display(vo_);
  }
}


/*******************************************************************************
 **
 **  Thermodynamic database file format
 **
 **  <Isotherms
 **  primary species name ; type ; parameters
 **
 **  type is one of: linear, langmuir, freundlich
 **
 **  parameters is a space delimited list of numbers. The number of
 **  parameters and their meaning depends on the isotherm type.
 **
 *******************************************************************************/
void SimpleThermoDatabase::ParseSorptionIsotherm(const std::string& data) {
  std::stringstream message;
  message << "SimpleThermoDatabase::ParseSorptionIsotherm()...." << std::endl
          << "  data: " << data << std::endl;
  vo_->Write(Teuchos::VERB_EXTREME, message);

  std::string semicolon(";");
  std::string space(" ");
  StringTokenizer no_spaces;
  // reaction_string; forward stoichs species; forward rate;
  //   backward stoichs species; backward rate

  StringTokenizer substrings(data, semicolon);

  SorptionIsothermFactory sif;

  // parse main reaction string
  no_spaces.tokenize(substrings.at(0),space);
  std::string species_name = no_spaces.at(0);
  no_spaces.tokenize(substrings.at(1),space);
  std::string isotherm_type = no_spaces.at(0);
  StringTokenizer parameters(substrings.at(2),space);

  SorptionIsotherm *sorption_isotherm = sif.Create(isotherm_type,
                                                   parameters);

  SorptionIsothermRxn rxn(species_name,SpeciesNameToID(species_name),
                          sorption_isotherm);
  this->AddSorptionIsothermRxn(rxn);
}


/*******************************************************************************
 **
 **  Thermodynamic database file format
 **
 **  <General Kinetics
 **
 *******************************************************************************/
void SimpleThermoDatabase::ParseGeneralKinetics(const std::string& data) {
  if (verbosity() == kDebugInputFile) {
    std::cout << "SimpleThermoDatabase::ParseGeneralKinetics()...." << std::endl;
    std::cout << "  data: " << data << std::endl;
  }
  std::string semicolon(";");
  std::string space(" ");
  std::string lr_arrow("<->");

  // reaction_string; forward stoichs species; forward rate;
  //   backward stoichs species; backward rate

  StringTokenizer substrings(data, semicolon);
  StringTokenizer substring;

  // parse main reaction string
  std::vector<SpeciesName> species;
  std::vector<double> stoichiometries;
  ParseReactionString(substrings.at(0), lr_arrow, &species, &stoichiometries);

  std::vector<int> species_ids;
  for (std::vector<SpeciesName>::iterator s = species.begin();
       s != species.end(); s++) {
    species_ids.push_back(SpeciesNameToID(*s));
  }

  // parse forward rates
  std::vector<double> forward_stoichiometries;
  std::vector<int> forward_species_ids;
  forward_stoichiometries.clear();
  forward_species_ids.clear();
  StringTokenizer forward_list(substrings.at(1),space);
  for (StringTokenizer::iterator i = forward_list.begin();
       i != forward_list.end(); ++i) {
    std::stringstream tempstring(std::stringstream::in |
                                 std::stringstream::out);
    tempstring << *i; // load string into buffer
    double d;
    tempstring >> d;
    if (tempstring.fail())
      d = 1.;  // no stoich specified
    else
      ++i; // increment substring
    std::string s = *i;  // read species name
    RemoveLeadingAndTrailingSpaces(&s);
    forward_species_ids.push_back(SpeciesNameToID(s));
    forward_stoichiometries.push_back(d);
  }

  // parse backward rates
  std::vector<double> backward_stoichiometries;
  std::vector<int> backward_species_ids;
  backward_stoichiometries.clear();
  backward_species_ids.clear();
  StringTokenizer backward_list(substrings.at(3),space);
  for (StringTokenizer::iterator i = backward_list.begin();
       i != backward_list.end(); ++i) {
    std::stringstream tempstring(std::stringstream::in |
                                 std::stringstream::out);
    tempstring << *i; // load string into buffer
    double d;
    tempstring >> d;
    if (tempstring.fail())
      d = 1.;  // no stoich specified
    else
      ++i; // increment substring
    std::string s = *i;  // read species name
    RemoveLeadingAndTrailingSpaces(&s);
    backward_species_ids.push_back(SpeciesNameToID(s));
    backward_stoichiometries.push_back(d);
  }

  double forward_rate_constant(0.);
  double backward_rate_constant(0.);
  std::stringstream tempstream(std::stringstream::in |
                               std::stringstream::out);
  if (substrings.at(2).size() > 0 &&
      substrings.at(2).find_first_not_of(space) != std::string::npos) {
    tempstream << substrings.at(2);
    tempstream >> forward_rate_constant;
    if (tempstream.fail()) std::cout << "Error reading forward rate constant\n";
  }
  if (substrings.at(4).size() > 0 &&
      substrings.at(4).find_first_not_of(space) != std::string::npos) {
    tempstream << substrings.at(4);
    tempstream >> backward_rate_constant;
    if (tempstream.fail()) std::cout << "Error reading backward rate constant\n";
  }

  GeneralRxn general("",species,stoichiometries,species_ids,
                     forward_stoichiometries,forward_species_ids,
                     backward_stoichiometries,backward_species_ids,
                     forward_rate_constant,backward_rate_constant);
  this->AddGeneralRxn(general);
}


/* ************************************************************************** /
   Name: ParseReactionString
   Purpose: Reads in a reaction string of format:
   30 A(aq) + 2 B(aq) <-> C(aq) + .3 D(aq) + -4 E(aq)
   and returns a list of species and ids:
   species = [5]("A(aq)","B(aq)","C(aq)","D(aq)","E(aq)")
   stoichiometries = [5](-30.,-2.,1.,0.3,-4.)
   Reactants and products have negative and positive stoichiometires,
   respectively.
   Author: Glenn Hammond
   Date: 07/28/11
   / ************************************************************************** */
void SimpleThermoDatabase::ParseReactionString(const std::string reaction,
                                               const std::string arrow,
                                               std::vector<std::string>* species,
                                               std::vector<double>* stoichiometries) {
  species->clear();
  stoichiometries->clear();
  std::string space(" ");
  std::string::size_type offset = reaction.find(arrow);
  std::string reactants = reaction.substr(0, offset);
  std::string products = reaction.substr(offset+arrow.size(),
                                         reaction.size());
  StringTokenizer list;
  // reactants
  list.tokenize_leave_delimiters(reactants,"+-");
  bool negate = true;
  for (StringTokenizer::iterator i = list.begin();
       i != list.end(); ++i) {
    if (*i == "-" || *i == "+" || *i == " ") {
      if (*i == "-") negate = !negate;
    }
    else {
      std::string s = *i;
      RemoveLeadingAndTrailingSpaces(&s);
      std::stringstream tempstring(std::stringstream::in |
                                   std::stringstream::out);
      if (s.size() > 0) {
        tempstring << s; // load string into buffer
        double d;
        tempstring >> d;  // read stoich
        if (tempstring.fail())
          d = 1.;  // no stoich specified
        else
          tempstring >> s; // read species name
        if (negate) d *= -1.; // swap sign on stoich
        RemoveLeadingAndTrailingSpaces(&s);
        species->push_back(s);
        stoichiometries->push_back(d);
        negate = true;
      }
    }
  }
  // products
  list.tokenize_leave_delimiters(products,"+-");
  negate = false;
  for (StringTokenizer::iterator i = list.begin();
       i != list.end(); ++i) {
    if (*i == "-" || *i == "+" || *i == " ") {
      if (*i == "-") negate = !negate;
    }
    else {
      std::string s = *i;
      RemoveLeadingAndTrailingSpaces(&s);
      std::stringstream tempstring(std::stringstream::in |
                                   std::stringstream::out);
      if (s.size() > 0) {
        tempstring << s; // load string into buffer
        double d;
        tempstring >> d;  // read stoich
        if (tempstring.fail())
          d = 1.;  // no stoich specified
        else
          tempstring >> s; // read species name
        if (negate) d *= -1.; // swap sign on stoich
        RemoveLeadingAndTrailingSpaces(&s);
        species->push_back(s);
        stoichiometries->push_back(d);
        negate = false;
      }
    }
  }
}


/* ************************************************************************** /
   Name: SpeciesNameToID
   Purpose: Maps a primary species name to an id
   Author: Glenn Hammond
   Date: 07/28/11
/ ************************************************************************** */
int SimpleThermoDatabase::SpeciesNameToID(const SpeciesName species_name) {
  for (SpeciesArray::const_iterator primary_species =
           this->primary_species().begin();
       primary_species != this->primary_species().end(); primary_species++) {
    if (primary_species->name() == species_name) {
      return primary_species->identifier();
    }
  }
  return -1;
}


/* ************************************************************************** /
Name: RemoveLeadingAndTrailingSpaces
Purpose: Removes leading and trailing spaces in a string
Author: Glenn Hammond
Date: 07/28/11
/ ************************************************************************** */
void SimpleThermoDatabase::RemoveLeadingAndTrailingSpaces(std::string* s) {
  std::size_t offset = s->find_first_not_of(" "); // leading spaces
  if (offset > 0) s->erase(0, offset);  // remove them
  offset = s->find_last_not_of(" ");   // trailing spaces
  if (s->size() > offset) s->erase(offset+1, s->size());  // remove them
}


/*******************************************************************************
**
**  Thermodynamic database file format
**
**  <Radioactive Decay
**  parent --> (coeff species (+ coeff species...)); half_life XXXX units
**
**  The stoichiometric coefficient of the parent should always be one!
** 
**  where units is one of: years, days, hours, minutes, seconds
**
*******************************************************************************/
void SimpleThermoDatabase::ParseRadioactiveDecay(const std::string& data) {

  std::stringstream message;
  message << "SimpleThermoDatabase::ParseRadioactiveDecay()....\n"
          << "  data: " << data << std::endl;
  vo_->Write(Teuchos::VERB_EXTREME, message);
  
  std::string semicolon(";");
  std::string space(" ");
  std::string right_arrow("-->");

  StringTokenizer substrings(data, semicolon);
  if (substrings.size() != 2) {
    message.str("");
    message << "ERROR: SimpleThermoDatabase::ParseRadioactiveDecay():\n"
            << "  Radioactive decay data should be of the form: "
            << "  'parent --> (coeff species (+ coeff species...)); half_life XXXX units'"
            << "  where the reaction is seperated from the half life by a semicolon.\n";
    Exceptions::amanzi_throw(ChemistryInvalidInput(message.str()));
  }

  std::vector<SpeciesName> species;
  std::vector<double> stoichiometry;
  std::vector<int> species_ids;

  // parse main reaction string
  StringTokenizer split_rxn(substrings.at(0), right_arrow);
  std::string parent = split_rxn.at(0);
  utilities::RemoveLeadingAndTrailingWhitespace(&parent);
  int parent_id = SpeciesNameToID(parent);
  if (parent_id < 0) {
    message.str("");
    message << "ERROR: SimpleThermoDatabase::ParseRadioactiveDecay():\n"
            << "  Unknown parent species '" << parent << "'.\n"
            << "  Parent species must be in the primary species list.\n";
    Exceptions::amanzi_throw(ChemistryInvalidInput(message.str()));
  }
  species.push_back(parent);
  species_ids.push_back(parent_id);
  stoichiometry.push_back(-1.0);

  std::string progeny_seperator(" + ");
  std::string progeny_string(split_rxn.at(1));
  utilities::RemoveLeadingAndTrailingWhitespace(&progeny_string);

  // NOTE: we allow zero progeny
  if (progeny_string.size() > 0) {
    size_t begin = 0;
    size_t end = 0;
    do {
      // grab the next progeny group
      begin = progeny_string.find_first_not_of(progeny_seperator);
      end = progeny_string.find(progeny_seperator, begin);
      std::string data_tmp = progeny_string.substr(begin, end);
      // remove the data we just grabbed
      progeny_string.erase(begin, end);

      // process the progeny group
      utilities::RemoveLeadingAndTrailingWhitespace(&data_tmp);
      StringTokenizer progeny_group(data_tmp, " ");
      std::string name;
      double coeff;
      if (progeny_group.size() == 1) {
        // no coeff provided, assume 1.0
        coeff = 1.0;
        name = progeny_group.at(0);
      } else if (progeny_group.size() == 2) {
        coeff = std::atof(progeny_group.at(0).c_str());
        name = progeny_group.at(1);
      } else {
        message.str("");
        message << "ERROR: SimpleThermoDatabase::ParseRadioactiveDecay():\n"
                << "  Progeny group '" << name << "' is in an unrecognized format.\n"
                << "  Progeny group format must be a name or 'coeff name'.\n";
        Exceptions::amanzi_throw(ChemistryInvalidInput(message.str()));
      }
      utilities::RemoveLeadingAndTrailingWhitespace(&name);  // needed?
      int id = SpeciesNameToID(name);
      if (id < 0) {
        message.str("");
        message << "ERROR: SimpleThermoDatabase::ParseRadioactiveDecay():\n"
                << "  Unknown progeny species '" << name << "'.\n"
                << "  Progeny species must be in the primary species list.\n";
        Exceptions::amanzi_throw(ChemistryInvalidInput(message.str()));
      }
      species.push_back(name);
      species_ids.push_back(id);
      stoichiometry.push_back(coeff);
    } while (end != std::string::npos);
  }
  double half_life(0.0);
  std::string half_life_units;
  std::string half_life_data = substrings.at(1);

  utilities::RemoveLeadingAndTrailingWhitespace(&half_life_data);
  StringTokenizer hl(half_life_data, space);
  assert(hl.size() == 3);
  half_life = std::atof(hl.at(1).c_str());
  half_life_units = hl.at(2);
  
  RadioactiveDecay rxn(species, species_ids, stoichiometry,
                       half_life, half_life_units);
  this->AddRadioactiveDecayRxn(rxn);
}


/*******************************************************************************
 **
 **  Thermodynamic database file format
 **
 **  <Mineral
 **
 **  Secondary Species Fields:
 **
 **  Name = coeff reactant ... ; log Keq ; gram molecular weight [g/mole] ; 
 **  molar volume [cm^3/mole] ; specific surface area [cm^2 mineral / cm^3 bulk]
 **
 *******************************************************************************/
void SimpleThermoDatabase::ParseMineral(const std::string& data) {
  if (verbosity() == kDebugInputFile) {
    std::cout << "SimpleThermoDatabase::ParseMineral()...." << std::endl;
    std::cout << "  data: " << data << std::endl;
  }
  std::string semicolon(";");
  std::string space(" \t");
  StringTokenizer no_spaces;

  StringTokenizer mineral_eq(data, semicolon);

  std::string name;
  std::vector<SpeciesName> species;
  std::vector<double> stoichiometries;
  std::vector<int> species_ids;
  double h2o_stoich = 0;
  std::string reaction(mineral_eq.at(0));
  ParseReaction(reaction, &name, &species, &stoichiometries, &species_ids, &h2o_stoich);

  no_spaces.tokenize(mineral_eq.at(1), space);
  double logKeq(std::atof(no_spaces.at(0).c_str()));

  no_spaces.tokenize(mineral_eq.at(2), space);
  double gram_molecular_weight(std::atof(no_spaces.at(0).c_str()));

  no_spaces.tokenize(mineral_eq.at(3), space);
  double molar_volume(std::atof(no_spaces.at(0).c_str()));
  // convert: [cm^3/mole] --> [m^3/mole]
  molar_volume /= 1000000.0;

  no_spaces.tokenize(mineral_eq.at(4), space);
  double specific_surface_area(std::atof(no_spaces.at(0).c_str()));
  // convert: [cm^2 mineral / cm^3 bulk] --> [m^2 mineral / m^3 bulk]
  specific_surface_area *= 100.0;

  Mineral mineral(name, mineral_id_++,
                  species,
                  stoichiometries,
                  species_ids,
                  h2o_stoich,
                  gram_molecular_weight,
                  logKeq,
                  molar_volume,
                  specific_surface_area);
  mineral.set_verbosity(verbosity());
  this->AddMineral(mineral);
  if (verbosity() == kDebugInputFile) {
    mineral.display();
  }
}


/*******************************************************************************
 **
 **  Thermodynamic database file format
 **
 **  <Mineral Kinetics
 **
 **  all rate information is contained on a single semicolon delimited line.
 **
 **  Field 0 : MineralName : assumed to have the same stoichiometry as the mineral definition
 **
 **  Field 1 : rate_name ("TST", ....)
 **
 **  remaining fields are processed by the rate class
 **
 *******************************************************************************/
void SimpleThermoDatabase::ParseMineralKinetics(const std::string& data) {
  if (verbosity() == kDebugInputFile) {
    std::cout << "SimpleThermoDatabase::ParseMineralKinetics()...." << std::endl;
    std::cout << "  data: " << data << std::endl;
  }
  std::string semicolon(";");
  std::string space(" \t");
  StringTokenizer no_spaces;

  StringTokenizer rate_data(data, semicolon);
  no_spaces.tokenize(rate_data.at(0), space);
  std::string mineral_name = no_spaces.at(0);
  std::string rate_type = rate_data.at(1);

  rate_data.erase(rate_data.begin());  // erase mineral name
  rate_data.erase(rate_data.begin());  // erase reaction type string

  MineralKineticsFactory mkf;
  mkf.set_debug(false);
  SpeciesId mineral_id = mkf.VerifyMineralName(mineral_name, minerals());
  Mineral mineral = minerals().at(mineral_id);
  KineticRate* kinetic_rate = mkf.Create(rate_type, rate_data, mineral, primary_species());

  this->AddMineralKineticRate(kinetic_rate);
  if (verbosity() == kDebugInputFile || verbosity() == kDebugMineralKinetics) {
    kinetic_rate->Display(vo_);
  }
}


/*******************************************************************************
 **
 **  Thermodynamic database file format
 **
 **  <Ion Exchange Site
 **
 **  Name ; change ; location
 **
 **  where location is the mineral where the exchanger is located, i.e. kaolinite
 **
 **  TODO(bandre): eventually something like "bulk" will be used as a dummy
 **  mineral for bulk soil rather than a specific mineral. need to
 **  coordinate this with surface complexation.
 **
 *******************************************************************************/
void SimpleThermoDatabase::ParseIonExchangeSite(const std::string& data) {
  if (verbosity() == kDebugInputFile) {
    std::cout << "SimpleThermoDatabase::ParseIonExchangeSite()...." << std::endl;
    std::cout << "  data: " << data << std::endl;
  }

  // double mol_wt = 0.0;  // not used in ion exchange sites
  // double size = 0.0;  // not used in ion exchange sites

  std::string semicolon(";");
  std::string space(" \t");
  StringTokenizer no_spaces;

  StringTokenizer exchanger_data(data, semicolon);
  no_spaces.tokenize(exchanger_data.at(0), space);
  std::string exchanger_name(no_spaces.at(0));

  no_spaces.tokenize(exchanger_data.at(1), space);
  double exchanger_charge(std::atof(no_spaces.at(0).c_str()));

  no_spaces.tokenize(exchanger_data.at(2), space);
  std::string exchanger_location(no_spaces.at(0));

  IonExchangeSite exchanger(exchanger_name, exchanger_charge, exchanger_location);
  IonExchangeRxn ionx_rxn(exchanger);

  this->AddIonExchangeRxn(ionx_rxn);
  if (verbosity() == kDebugInputFile) {
    exchanger.Display(vo_);
  }
}


/*******************************************************************************
 **
 **  Thermodynamic database file format
 **
 **  <Ion Exchange Complexes
 **
 **  Secondary Species Fields:
 **
 **  Name = coeff primary coeff exchanger ; Keq
 **
 **  Assume:
 **
 **    - that the coefficient of the ion exchange complex is one.
 **
 **    - each complexation reaction is written between a single
 **    primary species and a single exchange site
 **
 *******************************************************************************/
void SimpleThermoDatabase::ParseIonExchangeComplex(const std::string& data) {
  if (verbosity() == kDebugInputFile) {
    std::cout << "SimpleThermoDatabase::ParseIonExchangeComplex()...." << std::endl;
    std::cout << "  data: " << data << std::endl;
  }
  std::string semicolon(";");
  std::string space(" \t");
  StringTokenizer tokenizer;

  StringTokenizer complex_data(data, semicolon);

  tokenizer.tokenize(complex_data.at(0),"=");
  std::string complex_name = tokenizer.at(0);
  std::string products = tokenizer.at(1);
  tokenizer.tokenize(products,space);
  SpeciesName primary_name = tokenizer.at(1);
  IonxSiteName site_name = tokenizer.at(3);
  tokenizer.tokenize(complex_data.at(1), space);
  double K(std::atof(tokenizer.at(0).c_str()));

  int primary_id(-999);
  for (unsigned int i = 0; i < primary_species().size(); ++i) {
    if (primary_name == primary_species().at(i).name()) {
      primary_id = i;
      break;
    }
  }

  IonExchangeComplex exchange_complex(complex_name,
                                      ion_exchange_complex_id_++,
                                      primary_name,
                                      primary_id,
                                      K);

  for (unsigned int i = 0; i < ion_exchange_rxns().size(); i++) {
    if (site_name == ion_exchange_rxns().at(i).site().name()) {
      AddIonExchangeComplex(i,exchange_complex);
      break;
    }
  }

  if (verbosity() == kDebugInputFile) {
    exchange_complex.display(vo_);
  }
}


/*******************************************************************************
 **
 **  Thermodynamic database file format
 **
 **  <Surface Complex Site
 **
 **  Name ; density 
 **
 **  where location is the mineral where the exchanger is located, i.e. kaolinite
 **
 **  TODO(bandre): eventually something like "bulk" will be used as a dummy
 **  mineral for bulk soil rather than a specific mineral. need to
 **  coordinate this with surface complexation.
 **
 *******************************************************************************/
void SimpleThermoDatabase::ParseSurfaceComplexSite(const std::string& data) {
  if (verbosity() == kDebugInputFile) {
    std::cout << "SimpleThermoDatabase::ParseSurfaceComplexSite()...." << std::endl;
    std::cout << "  data: " << data << std::endl;
  }

  // double mol_wt = 0.0;  // not used in ion exchange sites
  // double size = 0.0;  // not used in ion exchange sites

  std::string semicolon(";");
  std::string space(" \t");
  StringTokenizer no_spaces;

  StringTokenizer site_data(data, semicolon);
  no_spaces.tokenize(site_data.at(0), space);
  std::string site_name(no_spaces.at(0));

  no_spaces.tokenize(site_data.at(1), space);
  double site_density(std::atof(no_spaces.at(0).c_str()));

  SurfaceSite site(site_name, surface_site_id_++,
                   site_density);
  surface_sites_.push_back(site);  // local storage to make parsing reactions easier...

  SurfaceComplexationRxn rxn(site);
  surface_complexation_reactions_.push_back(rxn);
  surface_complexation_rxn_id_++;

  // this->AddSurfaceComplexationRxn(rxn);
  if (verbosity() == kDebugInputFile) {
    site.display();
    // rxn.Display();
  }
}


/*******************************************************************************
 **
 **  Thermodynamic database file format
 **
 **  <Surface Complexes
 **
 **  Secondary Species Fields:
 **
 **  Name = coeff primary coeff exchanger ; Keq ; charge
 **
 **  Assume:
 **
 **    - ?
 **
 *******************************************************************************/
void SimpleThermoDatabase::ParseSurfaceComplex(const std::string& data) {
  if (verbosity() == kDebugInputFile) {
    std::cout << "SimpleThermoDatabase::ParseSurfaceComplex()...." << std::endl;
    std::cout << "  data: " << data << std::endl;
  }
  std::string semicolon(";");
  std::string space(" \t");
  StringTokenizer no_spaces;

  StringTokenizer complex_data(data, semicolon);

  std::string name;
  std::vector<SpeciesName> primary_name;
  std::vector<double> primary_stoichiometry;
  std::vector<SpeciesId> primary_id;
  SpeciesName surface_site_name;
  double surface_site_stoichiometry;
  SpeciesId surface_site_id;
  double h2o_stoich = 0;
  std::string reaction(complex_data.at(0));
  ParseSurfaceComplexReaction(reaction, &name,
                              &primary_name, &primary_stoichiometry, &primary_id,
                              &surface_site_name, &surface_site_stoichiometry, &surface_site_id,
                              &h2o_stoich);

  no_spaces.tokenize(complex_data.at(1), space);
  double logKeq(std::atof(no_spaces.at(0).c_str()));

  no_spaces.tokenize(complex_data.at(2), space);
  double charge(std::atof(no_spaces.at(0).c_str()));

  SurfaceComplex surface_complex(name,
                                 surface_complex_id_++,
                                 primary_name,
                                 primary_stoichiometry,
                                 primary_id,
                                 h2o_stoich,
                                 surface_site_name,
                                 surface_site_stoichiometry,
                                 surface_site_id,
                                 charge,
                                 logKeq);
  surface_complexation_reactions_[surface_site_id].AddSurfaceComplex(surface_complex);

  if (verbosity() == kDebugInputFile) {
    surface_complex.display(vo_);
  }
}  // end ParseSurfaceComplex()

void SimpleThermoDatabase::FinishSurfaceComplexation(void) {
  if (verbosity() == kDebugInputFile) {
    std::cout << "SimpleThermoDatabase::FinishSurfaceComplexation() :" << std::endl;
  }
  std::vector<SurfaceComplexationRxn>::iterator rxn;
  for (rxn = surface_complexation_reactions_.begin();
       rxn != surface_complexation_reactions_.end(); rxn++) {
    this->AddSurfaceComplexationRxn(*rxn);
  }
  // just to maintain our sanity, clear our temporary working copies
  // to make sure they aren't used anywhere else!
  surface_sites_.clear();
  surface_complexation_reactions_.clear();
}


/*******************************************************************************
 **
 **  pares a reaction whose products are all primary species
 **
 **  Fields:
 **
 **  SpeciesName = coeff PrimaryName coeff PrimaryName ...
 **
 *******************************************************************************/
void SimpleThermoDatabase::ParseReaction(const std::string& reaction,
                                         std::string* name,
                                         std::vector<SpeciesName>* species,
                                         std::vector<double>* stoichiometries,
                                         std::vector<int>* species_ids,
                                         double* h2o_stoich) {
  if (verbosity() == kDebugInputFile) {
    std::cout << "    SimpleThermoDatabase::ParseReaction()...." << std::endl;
    std::cout << "      data: " << reaction << std::endl;
  }
  std::string equal("=");
  std::string space(" \t");
  StringTokenizer no_spaces;
  std::vector<Species> primary_species = this->primary_species();

  StringTokenizer list(reaction, equal);
  no_spaces.tokenize(list.at(0), space);
  *name = no_spaces.at(0);
  // std::cout << "  name: " << *name << std::endl;

  StringTokenizer products(list.at(1), space);

  for (StringTokenizer::iterator s = products.begin();
       s != products.end(); s++) {
    double coeff = std::atof(s->c_str());
    s++;
    std::string primary_name(*s);

    // std::cout << "  " << coeff << " " << primary_name << std::endl;

    if (primary_name == "H2O") {
      *h2o_stoich = coeff;
      // std::cout << "  " << *h2o_stoich << " H2O" << std::endl;
    } else {
      species->push_back(primary_name);
      stoichiometries->push_back(coeff);
      int id = -1;
      for (SpeciesArray::const_iterator primary = primary_species.begin();
           primary != primary_species.end(); primary++) {
        // std::cout << "    " << primary->name() << "  ??  " << primary_name << std::endl;
        if (primary->name() == primary_name) {
          id = primary->identifier();
          primary = primary_species.end() - 1;
        }
      }
      if (id < 0) {
        std::cout << "SimpleThermoDatabase::ParseReaction(): reaction primary species \'"
                  << primary_name << "\' was not found in the primary species list...."
                  << std::endl;
      } else {
        species_ids->push_back(id);
      }
    }
  }  // end for(s)
}


/*******************************************************************************
 **
 **  parse ion exchange reaction, reaction products are single primary
 **  species and a single exchange site. The order of primary species
 **  and exchange species does not matter.
 **
 **  Fields:
 **
 **  SpeciesName = coeff PrimaryName coeff IonExchangeSite
 **
 *******************************************************************************/
void SimpleThermoDatabase::ParseIonExchangeReaction(const std::string& reaction,
                                                    std::string* name,
                                                    SpeciesName* primary_name,
                                                    SpeciesId* primary_id) {
  if (verbosity() == kDebugInputFile) {
    std::cout << "    SimpleThermoDatabase::ParseReaction()...." << std::endl;
    std::cout << "      data: " << reaction << std::endl;
  }
}


/*******************************************************************************
 **
 **  parse surface complex reaction, reaction products are an array of
 **  primary species and a single surface site. The order of
 **  primary species and surface sites does not matter.
 **
 **  Fields:
 **
 **  SpeciesName = coeff PrimaryName coeff SurfaceSite ...
 **
 *******************************************************************************/
void SimpleThermoDatabase::ParseSurfaceComplexReaction(const std::string& reaction,
                                                       std::string* name,
                                                       std::vector<SpeciesName>* primaries,
                                                       std::vector<double>* primary_stoichiometries,
                                                       std::vector<SpeciesId>* primary_ids,
                                                       SpeciesName* surface_site_name,
                                                       double* surface_site_stoichiometry,
                                                       SpeciesId* surface_site_id,
                                                       double* h2o_stoich) {
  if (verbosity() == kDebugInputFile) {
    std::cout << "    SimpleThermoDatabase::ParseReaction()...." << std::endl;
    std::cout << "      data: " << reaction << std::endl;
  }
  std::string equal("=");
  std::string space(" \t");
  StringTokenizer no_spaces;
  std::vector<Species> primary_species = this->primary_species();

  StringTokenizer list(reaction, equal);
  no_spaces.tokenize(list.at(0), space);
  *name = no_spaces.at(0);
  // std::cout << "  name: " << *name << std::endl;

  StringTokenizer products(list.at(1), space);

  for (StringTokenizer::iterator search_species = products.begin();
       search_species != products.end(); search_species++) {
    double coeff = std::atof(search_species->c_str());
    search_species++;
    std::string search_name(*search_species);
    // std::cout << "  " << coeff << " " << search_name << std::endl;

    if (search_name == "H2O") {
      *h2o_stoich = coeff;
      // std::cout << "  " << *h2o_stoich << " H2O" << std::endl;
    } else {
      // check to see if we have a primary species
      int id = -1;
      for (SpeciesArray::const_iterator primary = primary_species.begin();
           primary != primary_species.end(); primary++) {
        // std::cout << "    " << primary->name() << "  ??  " << search_name << std::endl;
        if (primary->name() == search_name) {
          id = primary->identifier();
          primary = primary_species.end() - 1;
        }
      }
      if (id >= 0) {
        // matched a primary species
        primaries->push_back(search_name);
        primary_stoichiometries->push_back(coeff);
        primary_ids->push_back(id);
      } else if (id < 0) {
        // did not match a primary. check to see if it is an surface site
        for (std::vector<SurfaceSite>::const_iterator surface = surface_sites_.begin();
             surface != surface_sites_.end(); surface++) {
          // std::cout << "    " << surface->name() << "  ??  " << search_name << std::endl;
          if (surface->name() == search_name) {
            id = surface->identifier();
            surface = surface_sites_.end() - 1;
          }
        }
        if (id >= 0) {
          // matched an surface site
          *surface_site_name = search_name;
          *surface_site_stoichiometry = coeff;
          *surface_site_id = id;
        }
      } else {
        // did not match an surface site or primary species
        std::cout << "SimpleThermoDatabase::ParseSurfaceComplexReaction(): reaction species \'"
                  << search_name << "\' was not found in the primary species "
                  << "list or the surface site list...."
                  << std::endl;
      }  // else unknown species
    }  // else not water
  }  // end for(search_species)
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
