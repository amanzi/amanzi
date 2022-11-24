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
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include "dbc.hh"
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
  : Beaker(vo.ptr()), plist_(plist){};


/* *******************************************************************
* Setup
******************************************************************* */
void
SimpleThermoDatabase::Initialize(BeakerState& state, const BeakerParameters& parameters)
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
  const auto& aux = plist_->sublist("aqueous equilibrium complexes");
  auto aqlist = RebuildAqueousComplexes_(aux);

  id = 0;
  for (auto it = aqlist.begin(); it != aqlist.end(); ++it, ++id) {
    std::string name = it->first;
    if (aqlist.isSublist(name)) {
      auto tmp = aqlist.sublist(name);
      tmp.set<double>("temperature", state.temperature);

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
  std::vector<SurfaceSite> surface_sites;
  std::vector<SurfaceComplexationRxn> surface_complexation_reactions;

  const auto& sslist = plist_->sublist("surface complex sites");

  id = 0;
  for (auto it = sslist.begin(); it != sslist.end(); ++it, ++id) {
    std::string name = it->first;
    if (sslist.isSublist(name)) {
      const auto& tmp = sslist.sublist(name);

      SurfaceSite site(name, id, tmp);
      surface_sites.push_back(site); // local storage to make parsing reactions easier

      SurfaceComplexationRxn rxn(site);
      surface_complexation_reactions.push_back(rxn);
    }
  }

  // surface complexes
  const auto& sclist = plist_->sublist("surface complexes");

  id = 0;
  for (auto it = sclist.begin(); it != sclist.end(); ++it, ++id) {
    std::string name = it->first;
    if (sclist.isSublist(name)) {
      const auto& tmp = sclist.sublist(name);

      SurfaceComplex surf_complex(name, id, primary_species(), surface_sites, tmp);
      surface_complexation_reactions[surf_complex.free_site_id()].AddSurfaceComplex(surf_complex);
    }
  }

  // sortion isotherms
  const auto& silist = plist_->sublist("isotherms");

  for (auto it = silist.begin(); it != silist.end(); ++it) {
    std::string name = it->first;
    if (silist.isSublist(name)) {
      const auto& tmp = silist.sublist(name);

      SorptionIsothermFactory sif;
      auto sorption_isotherm = sif.Create(tmp);

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

      RadioactiveDecay rxn(tmp, name_to_id);
      AddRadioactiveDecayRxn(rxn);
    }
  }

  // cleaning
  for (auto rxn = surface_complexation_reactions.begin();
       rxn != surface_complexation_reactions.end();
       rxn++) {
    this->AddSurfaceComplexationRxn(*rxn);
  }

  // this will allocate internal memory
  Beaker::Initialize(state, parameters);
}


/* *******************************************************************
* Setup
******************************************************************* */
Teuchos::ParameterList
SimpleThermoDatabase::RebuildAqueousComplexes_(const Teuchos::ParameterList& plist)
{
  Teuchos::ParameterList aqlist = plist;
  if (plist.numParams() == 0) return aqlist;

  std::vector<std::string> primaries;
  std::vector<std::string> secondaries;

  for (int i = 0; i < this->primary_species().size(); ++i) {
    primaries.push_back(primary_species()[i].name());
  }
  primaries.push_back("H2O");

  for (auto it = plist.begin(); it != plist.end(); ++it) {
    std::string name = it->first;
    if (plist.isSublist(name)) secondaries.push_back(name);
  }

  std::string empty("");
  int npri = primaries.size();
  int nsec = secondaries.size();

  MatrixBlock A(nsec, nsec), B(nsec, npri), Bnew(nsec, npri);
  A.Zero();
  B.Zero();

  // two ways to provide the equilibrium constant are reduced to one
  int nT(1);
  for (int i = 0; i < secondaries.size(); ++i) {
    auto& tmp = plist.sublist(secondaries[i]);
    if (tmp.isSublist("equilibrium constant")) {
      int n =
        tmp.sublist("equilibrium constant").get<Teuchos::Array<double>>("Keq").toVector().size();
      nT = std::max(nT, n);
    }
  }
  MatrixBlock logK(nsec, nT), logKnew(nsec, nT);

  for (int i = 0; i < secondaries.size(); ++i) {
    auto& tmp = plist.sublist(secondaries[i]);

    std::vector<double> stoich;
    std::vector<std::string> names;

    std::string reaction = tmp.get<std::string>("reaction");
    ParseReaction(reaction, empty, &names, &stoich);

    // two ways to provide the equilibrium constant
    if (!tmp.isSublist("equilibrium constant")) {
      double Keq = tmp.get<double>("equilibrium constant");
      for (int k = 0; k < nT; ++k) logK(i, k) = Keq;
    } else {
      auto Keq = tmp.sublist("equilibrium constant").get<Teuchos::Array<double>>("Keq").toVector();
      if (Keq.size() != nT) { AMANZI_ASSERT(false); }
      for (int k = 0; k < nT; ++k) logK(i, k) = Keq[k];
    }

    A(i, i) = 1.0;
    for (int k = 0; k < names.size(); ++k) {
      auto pos = std::find(secondaries.begin(), secondaries.end(), names[k]);
      if (pos != secondaries.end()) {
        int n = std::distance(secondaries.begin(), pos);
        A(i, n) = -stoich[k];
      }
    }

    for (int k = 0; k < names.size(); ++k) {
      auto pos = std::find(primaries.begin(), primaries.end(), names[k]);
      if (pos != primaries.end()) {
        int n = std::distance(primaries.begin(), pos);
        B(i, n) = -stoich[k];
      }
    }
  }

  // calculate new aqueous reactions
  A.Inverse();
  Multiply(A, B, Bnew);
  Multiply(A, logK, logKnew);

  for (int i = 0; i < secondaries.size(); ++i) {
    auto& tmp = aqlist.sublist(secondaries[i]);

    std::stringstream ss;
    for (int k = 0; k < npri; ++k) {
      if (Bnew(i, k) != 0.0) { ss << Bnew(i, k) << " " << primaries[k] << " "; }
    }
    tmp.set<std::string>("reaction", ss.str());

    // two ways to provide the equilibrium constant
    if (tmp.isSublist("equilibrium constant")) {
      Teuchos::Array<double> Keq(nT);
      for (int k = 0; k < nT; ++k) Keq[k] = logKnew(i, k);
      tmp.sublist("equilibrium constant").set<Teuchos::Array<double>>("Keq", Keq);
    } else {
      tmp.set<double>("equilibrium constant", logKnew(i, 0));
    }
  }

  return aqlist;
}

} // namespace AmanziChemistry
} // namespace Amanzi
