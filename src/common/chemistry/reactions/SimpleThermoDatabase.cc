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
#include "ReactionString.hh"
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
void SimpleThermoDatabase::Initialize(const BeakerState& state,
                                      const BeakerParameters& parameters)
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
      auto tmp = aqlist.sublist(name);
      tmp.set<double>("temperature", state.temperature);

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

      GeneralRxn general(tmp, name_to_id);
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

      IonExchangeSite exchanger(name, tmp);
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

      IonExchangeComplex ion_complex(name, id, tmp, primary_species());

      const std::string& site_name = ion_complex.site_name();
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

      SurfaceSite site(name, id, tmp);
      surface_sites_.push_back(site);  // local storage to make parsing reactions easier

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

      SurfaceComplex surf_complex(name, id, primary_species(), surface_sites_, tmp);
      surface_complexation_reactions_[surf_complex.free_site_id()].AddSurfaceComplex(surf_complex);
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
      RadioactiveDecay rxn(species, species_ids, stoichiometry, half_life);
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
  Beaker::Initialize(state, parameters);
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
