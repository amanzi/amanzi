/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "simple_thermo_database.hh"

#include <cstdlib>

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include "aqueous_equilibrium_complex.hh"
#include "mineral_kinetics_factory.hh"
#include "mineral.hh"
#include "surface_site.hh"
#include "surface_complex.hh"
#include "ion_exchange_site.hh"
#include "ion_exchange_complex.hh"
#include "beaker.hh"
#include "species.hh"
#include "string_tokenizer.hh"
#include "chemistry_exception.hh"

#include "exceptions.hh"

namespace amanzi {
namespace chemistry {

/*
**
** TODO(bandre): need a lot more error checking and helpfull error messages from here.
**
*/

SimpleThermoDatabase::SimpleThermoDatabase(void)
    : Beaker(),
      primary_id_(0),
      aqueous_equilibrium_complex_id_(0),
      mineral_id_(0),
      ion_exchange_site_id_(0),
      ion_exchange_complex_id_(0),
      surface_site_id_(0),
      surface_complex_id_(0),
      surface_complexation_rxn_id_(0) {
  surface_sites_.clear();
  surface_complexation_reactions_.clear();
}  // end SimpleThermoDatabase constructor

SimpleThermoDatabase::~SimpleThermoDatabase(void) {
}  // end SimpleThermoDatabase destructor

void SimpleThermoDatabase::Setup(const Beaker::BeakerComponents& components,
                                 const Beaker::BeakerParameters& parameters) {
  this->SetParameters(parameters);
  this->ReadFile(parameters.thermo_database_file);
  this->SetupActivityModel(parameters.activity_model_name);
  this->resize(this->primary_species().size());
  this->VerifyComponentSizes(components);
  this->SetComponents(components);
}  // end Setup()

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
    std::ostringstream error_stream;
    error_stream << "SimpleThermoDatabase::ReadFile(): \n";
    error_stream << "file could not be opened.... " << file_name << "\n";
    Exceptions::amanzi_throw(ChemistryInvalidInput(error_stream.str()));
  }

  enum LineType { kCommentLine, kPrimarySpeciesLine,
                  kAqueousEquilibriumComplexLine,
                  kMineralLine, kMineralKineticsLine, kGeneralKineticsLine,
                  kIonExchangeSiteLine, kIonExchangeComplexLine,
                  kSurfaceComplexSiteLine, kSurfaceComplexLine,
                  kUnknownLine
  };
  enum SectionType { kPrimarySpeciesSection, kAqueousEquilibriumComplexSection,
                     kMineralSection, kMineralKineticsSection, kGeneralKineticsSection,
                     kIonExchangeSiteSection, kIonExchangeComplexSection,
                     kSurfaceComplexSiteSection, kSurfaceComplexSection,
                     kUnknownSection
  };

  std::string kSectionPrimary("<Primary Species");
  std::string kSectionAqueousEquilibriumComplex("<Aqueous Equilibrium Complexes");
  std::string kSectionMineral("<Minerals");
  std::string kSectionMineralKinetics("<Mineral Kinetics");
  std::string kSectionGeneralKinetics("<General Kinetics");
  std::string kSectionIonExchangeSites("<Ion Exchange Sites");
  std::string kSectionIonExchangeComplexes("<Ion Exchange Complexes");
  std::string kSectionSurfaceComplexSites("<Surface Complex Sites");
  std::string kSectionSurfaceComplexes("<Surface Complexes");

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
      } else if (line == kSectionMineralKinetics) {
        line_type = kMineralKineticsLine;
        current_section = kMineralKineticsSection;
      } else if (line == kSectionGeneralKinetics) {
        line_type = kGeneralKineticsLine;
        current_section = kGeneralKineticsSection;
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
void SimpleThermoDatabase::ParsePrimarySpecies(const std::string& data) {
  if (verbosity() == kDebugInputFile) {
    std::cout << "SimpleThermoDatabase::ParsePrimarySpecies()...." << std::endl;
    std::cout << "  data: " << data << std::endl;
  }

  std::string semicolon(";");
  std::string space(" ");
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
  this->addPrimarySpecies(primary);
  if (verbosity() == kDebugInputFile) {
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
void SimpleThermoDatabase::ParseAqueousEquilibriumComplex(const std::string& data) {
  if (verbosity() == kDebugInputFile) {
    std::cout << "SimpleThermoDatabase::ParseAqueousEquilibriumComplex()...." << std::endl;
    std::cout << "  data: " << data << std::endl;
  }
  std::string semicolon(";");
  std::string space(" ");
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
  this->addAqueousEquilibriumComplex(secondary);
  if (verbosity() == kDebugInputFile) {
    secondary.display();
  }
}  // end ParseAqueousEquilibriumComplex()

/*******************************************************************************
 **
 **  Thermodynamic database file format
 **
 **  <General Kinetics
 **
 *******************************************************************************/
void SimpleThermoDatabase::ParseGeneralKinetics(const std::string& data) {
  std::ostringstream error_stream;
  error_stream << "SimpleThermoDatabase::ParseGeneralKinetics() : not implemented."
               << std::endl;
  Exceptions::amanzi_throw(ChemistryInvalidInput(error_stream.str()));
}  // end ParseGeneralKinetics()

/*******************************************************************************
 **
 **  Thermodynamic database file format
 **
 **  <Mineral
 **
 **  Secondary Species Fields:
 **
 **  Name = coeff reactant ... ; log Keq ; gram molecular weight [g/mole] ; molar volume [cm^3/mole] ; specific surface area [m^2]
 **
 *******************************************************************************/
void SimpleThermoDatabase::ParseMineral(const std::string& data) {
  if (verbosity() == kDebugInputFile) {
    std::cout << "SimpleThermoDatabase::ParseMineral()...." << std::endl;
    std::cout << "  data: " << data << std::endl;
  }
  std::string semicolon(";");
  std::string space(" ");
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

  no_spaces.tokenize(mineral_eq.at(4), space);
  double specific_surface_area(std::atof(no_spaces.at(0).c_str()));

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
  this->addMineral(mineral);
  if (verbosity() == kDebugInputFile) {
    mineral.display();
  }
}  // end ParseMineral()

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
  std::string space(" ");
  StringTokenizer no_spaces;

  StringTokenizer rate_data(data, semicolon);
  no_spaces.tokenize(rate_data.at(0), space);
  std::string mineral_name = no_spaces.at(0);
  std::string rate_type = rate_data.at(1);

  rate_data.erase(rate_data.begin());  // erase mineral name
  rate_data.erase(rate_data.begin());  // erase reaction type string

  MineralKineticsFactory mkf;
  mkf.set_verbosity(verbosity());
  SpeciesId mineral_id = mkf.VerifyMineralName(mineral_name, minerals());
  Mineral mineral = minerals().at(mineral_id);
  KineticRate* kinetic_rate = mkf.Create(rate_type, rate_data, mineral, primary_species());

  this->AddMineralKineticRate(kinetic_rate);
  if (verbosity() == kDebugInputFile || verbosity() == kDebugMineralKinetics) {
    kinetic_rate->Display();
  }
}  // end ParseMineralKinetics()


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

  double mol_wt = 0.0;  // not used in ion exchange sites
  double size = 0.0;  // not used in ion exchange sites

  std::string semicolon(";");
  std::string space(" ");
  StringTokenizer no_spaces;

  StringTokenizer exchanger_data(data, semicolon);
  no_spaces.tokenize(exchanger_data.at(0), space);
  std::string exchanger_name(no_spaces.at(0));

  no_spaces.tokenize(exchanger_data.at(1), space);
  double exchanger_charge(std::atof(no_spaces.at(0).c_str()));

  no_spaces.tokenize(exchanger_data.at(2), space);
  std::string exchanger_location(no_spaces.at(0));

  IonExchangeSite exchanger(exchanger_name, ion_exchange_site_id_++,
                            exchanger_charge, exchanger_location,
                            mol_wt, size);

  this->AddIonExchangeSite(exchanger);
  if (verbosity() == kDebugInputFile) {
    exchanger.display();
  }
}  // end ParseIonExchangeSite()

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
  std::string space(" ");
  StringTokenizer no_spaces;

  StringTokenizer complex_data(data, semicolon);

  std::string name;
  SpeciesName primary_name;
  double primary_stoichiometry;
  SpeciesId primary_id;
  SpeciesName exchange_site_name;
  double exchange_site_stoichiometry;
  SpeciesId exchange_site_id;
  double h2o_stoich = 0;
  std::string reaction(complex_data.at(0));
  ParseIonExchangeReaction(reaction, &name,
                           &primary_name, &primary_stoichiometry, &primary_id,
                           &exchange_site_name, &exchange_site_stoichiometry, &exchange_site_id,
                           &h2o_stoich);

  no_spaces.tokenize(complex_data.at(1), space);
  double logKeq(std::atof(no_spaces.at(0).c_str()));

  IonExchangeComplex exchange_complex(name,
                                      ion_exchange_complex_id_++,
                                      primary_name,
                                      primary_stoichiometry,
                                      primary_id,
                                      exchange_site_name,
                                      exchange_site_stoichiometry,
                                      exchange_site_id,
                                      h2o_stoich,
                                      logKeq);
  this->AddIonExchangeComplex(exchange_complex);

  if (verbosity() == kDebugInputFile) {
    exchange_complex.display();
  }
}  // end ParseIonExchangeComplex()

/*******************************************************************************
 **
 **  Thermodynamic database file format
 **
 **  <Surface Complex Site
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
void SimpleThermoDatabase::ParseSurfaceComplexSite(const std::string& data) {
  if (verbosity() == kDebugInputFile) {
    std::cout << "SimpleThermoDatabase::ParseSurfaceComplexSite()...." << std::endl;
    std::cout << "  data: " << data << std::endl;
  }

  double mol_wt = 0.0;  // not used in ion exchange sites
  double size = 0.0;  // not used in ion exchange sites

  std::string semicolon(";");
  std::string space(" ");
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

  // this->addSurfaceComplexationRxn(rxn);
  if (verbosity() == kDebugInputFile) {
    site.display();
    // rxn.Display();
  }
}  // end ParseSurfaceComplexSite()

/*******************************************************************************
 **
 **  Thermodynamic database file format
 **
 **  <Surface Complexes
 **
 **  Secondary Species Fields:
 **
 **  Name = coeff primary coeff exchanger ; Keq
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
  std::string space(" ");
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
    surface_complex.display();
  }
}  // end ParseSurfaceComplex()

void SimpleThermoDatabase::FinishSurfaceComplexation(void) {
  if (verbosity() == kDebugInputFile) {
    std::cout << "SimpleThermoDatabase::FinishSurfaceComplexation() :" << std::endl;
  }
  std::vector<SurfaceComplexationRxn>::iterator rxn;
  for (rxn = surface_complexation_reactions_.begin();
       rxn != surface_complexation_reactions_.end(); rxn++) {
    this->addSurfaceComplexationRxn(*rxn);
  }
}  // end FinishSurfaceComplexation()

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
  std::string space(" ");
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
}  // end ParseReaction()

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
                                                    double* primary_stoichiometry,
                                                    SpeciesId* primary_id,
                                                    SpeciesName* exchanger_name,
                                                    double* exchanger_stoichiometry,
                                                    SpeciesId* exchanger_id,
                                                    double* h2o_stoich) {
  if (verbosity() == kDebugInputFile) {
    std::cout << "    SimpleThermoDatabase::ParseReaction()...." << std::endl;
    std::cout << "      data: " << reaction << std::endl;
  }
  std::string equal("=");
  std::string space(" ");
  StringTokenizer no_spaces;
  std::vector<Species> primary_species = this->primary_species();
  std::vector<IonExchangeSite> ion_exchange_sites = this->ion_exchange_sites();

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
        *primary_name = search_name;
        *primary_stoichiometry = coeff;
        *primary_id = id;
      } else if (id < 0) {
        // did not match a primary. check to see if it is an ion exchange site
        for (std::vector<IonExchangeSite>::const_iterator exchanger = ion_exchange_sites.begin();
             exchanger != ion_exchange_sites.end(); exchanger++) {
          if (exchanger->name() == search_name) {
            id = exchanger->identifier();
            exchanger = ion_exchange_sites.end() - 1;
          }
        }
        if (id >= 0) {
          // matched an exchange site
          *exchanger_name = search_name;
          *exchanger_stoichiometry = coeff;
          *exchanger_id = id;
        }
      } else {
        // did not match an exchange site or primary species
        std::cout << "SimpleThermoDatabase::ParseIonExchangeReaction(): reaction species \'"
                  << search_name << "\' was not found in the primary species "
                  << "list or the ion exchange site list...."
                  << std::endl;
      }  // else unknown species
    }  // else not water
  }  // end for(search_species)
}  // end ParseIonExchangeReaction()

/*******************************************************************************
 **
 **  parse ion exchange reaction, reaction products are an array of
 **  primary species and an array of exchange sites. The order of
 **  primary species and exchange sites does not matter.
 **
 **  Fields:
 **
 **  SpeciesName = coeff PrimaryName coeff IonExchangeSite ...
 **
 *******************************************************************************/
void SimpleThermoDatabase::ParseIonExchangeReaction(const std::string& reaction,
                                                    std::string* name,
                                                    std::vector<SpeciesName>* primaries,
                                                    std::vector<double>* primary_stoichiometries,
                                                    std::vector<SpeciesId>* primary_ids,
                                                    std::vector<SpeciesName>* exchange_sites,
                                                    std::vector<double>* exchanger_stoichiometries,
                                                    std::vector<SpeciesId>* exchanger_ids,
                                                    double* h2o_stoich) {
  if (verbosity() == kDebugInputFile) {
    std::cout << "    SimpleThermoDatabase::ParseReaction()...." << std::endl;
    std::cout << "      data: " << reaction << std::endl;
  }
  std::string equal("=");
  std::string space(" ");
  StringTokenizer no_spaces;
  std::vector<Species> primary_species = this->primary_species();
  std::vector<IonExchangeSite> ion_exchange_sites = this->ion_exchange_sites();

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
        // did not match a primary. check to see if it is an ion exchange site
        for (std::vector<IonExchangeSite>::const_iterator exchanger = ion_exchange_sites.begin();
             exchanger != ion_exchange_sites.end(); exchanger++) {
          if (exchanger->name() == search_name) {
            id = exchanger->identifier();
            exchanger = ion_exchange_sites.end() - 1;
          }
        }
        if (id >= 0) {
          // matched an exchange site
          exchange_sites->push_back(search_name);
          exchanger_stoichiometries->push_back(coeff);
          exchanger_ids->push_back(id);
        }
      } else {
        // did not match an exchange site or primary species
        std::cout << "SimpleThermoDatabase::ParseIonExchangeReaction(): reaction species \'"
                  << search_name << "\' was not found in the primary species "
                  << "list or the ion exchange site list...."
                  << std::endl;
      }  // else unknown species
    }  // else not water
  }  // end for(search_species)
}  // end ParseIonExchangeReaction()

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
  std::string space(" ");
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
}  // end ParseSurfaceComplexReaction()

}  // namespace chemistry
}  // namespace amanzi
