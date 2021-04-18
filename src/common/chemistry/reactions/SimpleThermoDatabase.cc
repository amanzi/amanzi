/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Glenn Hammond
*/
 
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include "exceptions.hh"
#include "errors.hh"

#include "AqueousEquilibriumComplex.hh"
#include "Beaker.hh"
#include "GeneralRxn.hh"
#include "IonExchangeComplex.hh"
#include "IonExchangeRxn.hh"
#include "IonExchangeSite.hh"
#include "KineticRateFactory.hh"
#include "Mineral.hh"
#include "RadioactiveDecay.hh"
#include "Species.hh"
#include "SorptionIsotherm.hh"
#include "SorptionIsothermFactory.hh"
#include "SorptionIsothermLinear.hh"
#include "SurfaceSite.hh"
#include "SurfaceComplex.hh"

#include "SimpleThermoDatabase.hh"

namespace Amanzi {
namespace AmanziChemistry {

/* *******************************************************************
* TBW
******************************************************************* */
SimpleThermoDatabase::SimpleThermoDatabase(Teuchos::RCP<Teuchos::ParameterList> plist,
                                           Teuchos::RCP<VerboseObject> vo)
  : plist_(plist),
    Beaker(vo.ptr()) {
  surface_sites_.clear();
  surface_complexation_reactions_.clear();
}


/* *******************************************************************
* Setup
******************************************************************* */
void SimpleThermoDatabase::Initialize(const BeakerParameters& parameters)
{
  // primary species
  const auto& pslist = plist_->sublist("primary species");
  
  int id(0);
  primary_species_.clear();
  std::map<std::string, int> name_to_id;

  for (auto it = pslist.begin(); it != pslist.end(); ++it, ++id) {
    std::string name = it->first;
    if (pslist.isSublist(name)) {
      const auto& tmp = pslist.sublist(name);
      Species primary(id, name, tmp);
      primary_species_.push_back(primary);
      name_to_id[name] = id;
    }
  }

  // aqueous complexes
  const auto& aqlist = plist_->sublist("aqueous equilibrium complexes");

  id = 0;
  for (auto it = aqlist.begin(); it != aqlist.end(); ++it, ++id) {
    std::string name = it->first;
    if (aqlist.isSublist(name)) {
      const auto& tmp = aqlist.sublist(name);

      double h2o_stoich(0.0);

      AqueousEquilibriumComplex secondary(id, name, tmp, this->primary_species());
      AddAqueousEquilibriumComplex(secondary);
    }
  }

  // minerals
  const auto& mnlist = plist_->sublist("mineral kinetics");

  id = 0;
  for (auto it = mnlist.begin(); it != mnlist.end(); ++it, ++id) {
    std::string name = it->first;
    if (mnlist.isSublist(name)) {
      const auto& tmp = mnlist.sublist(name);

      Mineral mineral(id, name, tmp, this->primary_species());
      AddMineral(mineral);

      KineticRateFactory krf;
      KineticRate* kinetic_rate = krf.Create(tmp, mineral, primary_species());

      AddMineralKineticRate(kinetic_rate);
    }
  }

  // general kinetics
  const auto& gklist = plist_->sublist("general kinetics");

  for (auto it = gklist.begin(); it != gklist.end(); ++it) {
    std::string name = it->first;
    if (gklist.isSublist(name)) {
      const auto& tmp = gklist.sublist(name);

      // parse main reaction string
      std::string reactants = tmp.get<std::string>("reactants");
      std::string products = tmp.get<std::string>("products");

      std::vector<std::string> species;
      std::vector<double> stoichiometries;
      ParseReaction_(reactants, products, &species, &stoichiometries);

      std::vector<int> species_ids;
      for (auto s = species.begin(); s != species.end(); s++) {
        species_ids.push_back(name_to_id.at(*s));
      }

      // parse forward rates
      double coeff;
      std::vector<double> forward_stoichiometries;
      std::vector<int> forward_species_ids;

      std::istringstream iss1(reactants);
      while (iss1 >> coeff || !iss1.eof()) {
        iss1 >> name;
        forward_species_ids.push_back(name_to_id.at(name));
        forward_stoichiometries.push_back(coeff);
      }

      // parse backward rates
      std::vector<double> backward_stoichiometries;
      std::vector<int> backward_species_ids;

      std::istringstream iss2(reactants);
      while (iss2 >> coeff || !iss2.eof()) {
        iss2 >> name;
        backward_species_ids.push_back(name_to_id.at(name));
        backward_stoichiometries.push_back(coeff);
      }

      double forward_rate_constant = tmp.get<double>("forward rate");
      double backward_rate_constant = tmp.get<double>("backward rate");

      GeneralRxn general("", species, stoichiometries, species_ids,
                         forward_stoichiometries, forward_species_ids,
                         backward_stoichiometries, backward_species_ids,
                         forward_rate_constant, backward_rate_constant);
      AddGeneralRxn(general);
    }
  }

  // ion excange sites
  // location is the mineral where the exchanger is located, i.e. kaolinite
  // TODO(bandre): eventually something like "bulk" will be used as a dummy
  // mineral for bulk soil rather than a specific mineral. need to
  // coordinate this with surface complexation.
  const auto& islist = plist_->sublist("ion exchange sites");

  for (auto it = islist.begin(); it != islist.end(); ++it) {
    std::string name = it->first;
    if (islist.isSublist(name)) {
      const auto& tmp = islist.sublist(name);

      double charge = tmp.get<int>("charge");
      std::string location = tmp.get<std::string>("location");

      IonExchangeSite exchanger(name, charge, location);
      IonExchangeRxn ionx_rxn(exchanger);

      AddIonExchangeRxn(ionx_rxn);
    }
  }

  // ion exchange complexes
  // We assume that (a) the coefficient of the ion exchange complex is one
  // and (b) each complexation reaction is written between a single primary
  // species and a single exchange site.
  const auto& iclist = plist_->sublist("ion exchange complexes");

  id = 0;
  for (auto it = iclist.begin(); it != iclist.end(); ++it, ++id) {
    std::string name = it->first;
    if (iclist.isSublist(name)) {
      const auto& tmp = iclist.sublist(name);

      double coeff;
      std::string primary_name, site_name;

      double lnKeq = tmp.get<double>("equilibrium constant");
      std::istringstream iss(tmp.get<std::string>("reaction"));
      iss >> coeff;
      iss >> primary_name;
      iss >> coeff;
      iss >> site_name;

      int primary_id(-999);
      for (int i = 0; i < primary_species().size(); ++i) {
        if (primary_name == primary_species().at(i).name()) {
          primary_id = i;
          break;
        }
      }

      IonExchangeComplex ion_complex(name, id,
                                     primary_name, primary_id,
                                     lnKeq);

      for (int i = 0; i < ion_exchange_rxns().size(); i++) {
        if (site_name == ion_exchange_rxns().at(i).site().get_name()) {
          AddIonExchangeComplex(i, ion_complex);
          break;
        }
      }
    }
  }

  // surface complex sites
  const auto& sslist = plist_->sublist("surface complex sites");

  id = 0;
  for (auto it = sslist.begin(); it != sslist.end(); ++it, ++id) {
    std::string name = it->first;
    if (sslist.isSublist(name)) {
      const auto& tmp = sslist.sublist(name);

      double density = tmp.get<double>("density");

      SurfaceSite site(name, id, density);
      surface_sites_.push_back(site);  // local storage to make parsing reactions easier...

      SurfaceComplexationRxn rxn(site);
      surface_complexation_reactions_.push_back(rxn);
    }
  }

  // surface complexes
  const auto& sclist = plist_->sublist("surface complexes");

  id = 0;
  for (auto it = sclist.begin(); it != sclist.end(); ++it, ++id) {
    std::string name = it->first;
    if (sclist.isSublist(name)) {
      const auto& tmp = sclist.sublist(name);

      std::string reaction = tmp.get<std::string>("reaction");
      double lnKeq = tmp.get<double>("equilibrium constant");
      double charge = tmp.get<int>("charge");

      std::vector<std::string> primary_name;
      std::vector<double> primary_stoichiometry;
      std::vector<int> primary_id;
      std::string surface_site_name;
      double surface_site_stoichiometry, h2o_stoich(0.0);
      int surface_site_id;

      ParseReaction_(reaction,
                     &primary_name, &primary_stoichiometry, &primary_id,
                     &surface_site_name, &surface_site_stoichiometry, &surface_site_id,
                     &h2o_stoich);

      SurfaceComplex surf_complex(name, id,
                                  primary_name, primary_stoichiometry, primary_id,
                                  h2o_stoich,
                                  surface_site_name, surface_site_stoichiometry, surface_site_id,
                                  charge,
                                  lnKeq);

      surface_complexation_reactions_[surface_site_id].AddSurfaceComplex(surf_complex);
    }
  }

  // sortion isotherms
  const auto& silist = plist_->sublist("isotherms");

  for (auto it = silist.begin(); it != silist.end(); ++it) {
    std::string name = it->first;
    if (silist.isSublist(name)) {
      const auto& tmp = silist.sublist(name);

      std::string model = tmp.get<std::string>("model");
      std::vector<double> params = tmp.get<Teuchos::Array<double> >("parameters").toVector();
      
      SorptionIsothermFactory sif;
      auto sorption_isotherm = sif.Create(model, params);

      SorptionIsothermRxn rxn(name, name_to_id.at(name), sorption_isotherm);
      AddSorptionIsothermRxn(rxn);
    }
  }

  // parent --> (coeff species (+ coeff species...)); half_life XXXX units
  // The stoichiometric coefficient of the parent should always be one!
  // where units is one of: years, days, hours, minutes, seconds
  const auto& rdlist = plist_->sublist("radioactive decay");

  for (auto it = rdlist.begin(); it != rdlist.end(); ++it) {
    std::string name = it->first;
    if (rdlist.isSublist(name)) {
      const auto& tmp = rdlist.sublist(name);

      std::vector<std::string> species;
      std::vector<double> stoichiometry;
      std::vector<int> species_ids;

      std::string parent = tmp.get<std::string>("reactant");
      int parent_id = name_to_id.at(parent);
      if (parent_id < 0) {
        std::stringstream ss;
        ss << "Unknown parent species '" << parent << "'.\n"
           << "Parent species must be in the primary species list.\n";
        Exceptions::amanzi_throw(Errors::Message(ss.str()));
      }
      species.push_back(parent);
      species_ids.push_back(parent_id);
      stoichiometry.push_back(-1.0);

      std::string progeny = tmp.get<std::string>("product");

      // NOTE: we allow zero progeny
      if (progeny.size() > 0) {
        int id2 = name_to_id.at(progeny);
        if (id2 < 0) {
          std::stringstream ss;
          ss << "Unknown progeny species '" << name << "'.\n"
             << "Progeny species must be in the primary species list.\n";
          Exceptions::amanzi_throw(Errors::Message(ss.str()));
        }
        species.push_back(progeny);
        species_ids.push_back(id2);
        stoichiometry.push_back(1.0);
      }

      double half_life = tmp.get<double>("half life");
      std::string units("seconds");
  
      RadioactiveDecay rxn(species, species_ids, stoichiometry, half_life, units);
      AddRadioactiveDecayRxn(rxn);
    }
  }

  // cleaning
  for (auto rxn = surface_complexation_reactions_.begin();
       rxn != surface_complexation_reactions_.end(); rxn++) {
    this->AddSurfaceComplexationRxn(*rxn);
  }
  surface_sites_.clear();
  surface_complexation_reactions_.clear();

  // this will allocate internal memory
  Beaker::Initialize(parameters);
}


/* *******************************************************************
* Parse surface complex reaction, reaction products are an array of
* primary species and a single surface site. The order of
* primary species and surface sites does not matter.
*
*  SpeciesName = coeff PrimaryName coeff SurfaceSite ...
******************************************************************* */
void SimpleThermoDatabase::ParseReaction_(const std::string& reaction,
                                          std::vector<std::string>* primaries,
                                          std::vector<double>* primary_stoichiometries,
                                          std::vector<int>* primary_ids,
                                          std::string* surface_site_name,
                                          double* surface_site_stoichiometry,
                                          int* surface_site_id,
                                          double* h2o_stoich)
{
  double coeff;
  std::string search_name;
  std::vector<Species> primary_species = this->primary_species();

  std::istringstream iss(reaction);
  while (iss >> coeff || !iss.eof()) {
    iss >> search_name;

    if (search_name == "H2O") {
      *h2o_stoich = coeff;
    } else {
      // check to see if we have a primary species
      int id = -1;
      for (auto it = primary_species.begin(); it != primary_species.end(); ++it) {
        if (it->name() == search_name) {
          id = it->identifier();
          break;
        }
      }

      if (id >= 0) {
        primaries->push_back(search_name);
        primary_stoichiometries->push_back(coeff);
        primary_ids->push_back(id);

      } else if (id < 0) {
        // did not match a primary. check to see if it is an surface site
        for (auto it = surface_sites_.begin(); it != surface_sites_.end(); ++it) {
          if (it->name() == search_name) {
            id = it->identifier();
            break;
          }
        }

        if (id >= 0) {
          *surface_site_name = search_name;
          *surface_site_stoichiometry = coeff;
          *surface_site_id = id;
        }

      } else {
        // did not match an surface site or primary species
        std::cout << "Reaction species \'" << search_name << "\' was not found in the primary"
                  << " species list or the surface site list.\n";
      }
    }
  }
}


/* *******************************************************************
* Reads in a reaction string of format: reactants -> products
* Example:
*   30 A(aq) + 2 B(aq) <-> C(aq) + .3 D(aq) + -4 E(aq)
* returns 
*   species = [5]("A(aq)","B(aq)","C(aq)","D(aq)","E(aq)")
*   stoichiometries = [5](-30.,-2.,1.,0.3,-4.)
* NOTE: Reactants and products have negative and positive 
* stoichiometries, respectively.
******************************************************************* */
void SimpleThermoDatabase::ParseReaction_(const std::string& reactants,
                                          const std::string& products,
                                          std::vector<std::string>* species,
                                          std::vector<double>* stoichiometries)
{
  double coeff;
  std::string name;

  species->clear();
  stoichiometries->clear();

  // reactants
  std::istringstream iss1(reactants);
  while (iss1 >> coeff || !iss1.eof()) {
    iss1 >> name;

    species->push_back(name);
    stoichiometries->push_back(-coeff);
  }

  // products
  std::istringstream iss2(products);
  while (iss2 >> coeff || !iss2.eof()) {
    iss2 >> name;

    species->push_back(name);
    stoichiometries->push_back(coeff);
  }
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
