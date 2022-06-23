/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Class for ion exchange reaction
*/

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "MatrixBlock.hh"

#include "exceptions.hh"
#include "IonExchangeRxn.hh"

namespace Amanzi {
namespace AmanziChemistry {

IonExchangeRxn::IonExchangeRxn()
  : uniform_z_(false),
    uniform_z_set_(false),
    ref_cation_sorbed_conc_(1.e-9)
{
  ionx_site_.clear();
  ionx_complexes_.clear();
}


IonExchangeRxn::IonExchangeRxn(
    IonExchangeSite* ionx_sites,
    const std::vector<IonExchangeComplex>& ionx_complexes) 
    : uniform_z_(false),
      uniform_z_set_(false),
      ref_cation_sorbed_conc_(1.e-9) {
  // surface site
  ionx_site_.push_back(*ionx_sites);

  // surface complexes
  for (std::vector<IonExchangeComplex>::const_iterator i = ionx_complexes.begin();
       i != ionx_complexes.end(); i++) {
    ionx_complexes_.push_back(*i);
  }
}


IonExchangeRxn::IonExchangeRxn(IonExchangeSite ionx_sites) 
    : uniform_z_(false),
      uniform_z_set_(false),
      ref_cation_sorbed_conc_(1.e-9)
{
  ionx_site_.push_back(ionx_sites);
  ionx_complexes_.clear();
}


void IonExchangeRxn::AddIonExchangeComplex(const IonExchangeComplex& ionx_complex) {
  ionx_complexes_.push_back(ionx_complex);
}


void IonExchangeRxn::AddIonExchangeSite(const IonExchangeSite& site) {
  ionx_site_.push_back(site);
}


void IonExchangeRxn::Update(const std::vector<Species>& primary_species)
{
  // pflotran: reaction.F90, function RTotalSorbEqIonX

  // Note: "X" is the fraction of sites occupied by a particular
  // cation (accounting for charge). X = sorbed_conc * charge / CEC
  // sum(X) = 1 (res == 0), should be a requirement for convergence of
  // the update loop.
  bool one_more;
  double tol = 1.e-12;

  if (!uniform_z_set_) {
    CheckUniformZ(primary_species);
  }

  for (auto ionx = ionx_complexes_.begin(); ionx != ionx_complexes_.end(); ionx++) {
    ionx->set_X(0.0);
  }

  double omega = ionx_site_[0].get_cation_exchange_capacity();
  if (!uniform_z_) { // Z_i /= Z_j for all i,j
    int interation_count = 0;
    int ref_cation = ionx_complexes_[0].primary_id();
    double ref_cation_act = primary_species[ref_cation].activity();
    double ref_cation_Z = primary_species[ref_cation].charge();
    double ref_cation_K = ionx_complexes_[0].K();
    double ref_cation_X = ref_cation_Z*ref_cation_sorbed_conc_/omega;

    one_more = false;
    while(true) {
      if (ref_cation_X < 0.0) ref_cation_X = 0.99;
      ionx_complexes_[0].set_X(ref_cation_X);
      double ref_cation_quotient = ref_cation_X*ref_cation_K/ref_cation_act;
      double total = ref_cation_X;
      for (int i = 1; i < ionx_complexes_.size(); i++) {
        int icomp = ionx_complexes_[i].primary_id();
        double value = primary_species[icomp].activity() / ionx_complexes_[i].K() *
                       pow(ref_cation_quotient, primary_species[icomp].charge() / ref_cation_Z);
        ionx_complexes_[i].set_X(value);
        total += value;
      }
      if (false) {
        std::cout << "-- ionx total: " << total << std::endl;
      }
      ++interation_count;
      if (one_more) break;

      double res = 1.0 - total;
      double dres_dref_cation_X = 1.0;
      for (int i = 1; i < ionx_complexes_.size(); i++) {
        int icomp = ionx_complexes_[i].primary_id();
        dres_dref_cation_X += (primary_species[icomp].charge() / ref_cation_Z) *
            (ionx_complexes_[i].X() / ref_cation_X);
      }
      double dref_cation_X = -res / dres_dref_cation_X;
      ref_cation_X -= dref_cation_X;
      if (std::fabs(dref_cation_X/ref_cation_X) < tol && std::fabs(res) < tol) {
        one_more = true;
      }
    }
    ref_cation_sorbed_conc_ = ref_cation_X*omega/ref_cation_Z;
  }
  else { // Z_i == Z_j for all i,j
    double sumkm = 0.0;
    for (auto ionx = ionx_complexes_.begin(); ionx != ionx_complexes_.end(); ionx++) {
      int icomp = ionx->primary_id();
      double value = primary_species[icomp].activity() * ionx->K();
      ionx->set_X(value);
      sumkm += value;
    }
    for (auto ionx = ionx_complexes_.begin(); ionx != ionx_complexes_.end(); ionx++) {
      double temp = ionx->X() / sumkm;
      ionx->set_X(temp);
    }
  }

  for (auto ionx = ionx_complexes_.begin(); ionx != ionx_complexes_.end(); ionx++) {
    int icomp = ionx->primary_id();
    // NOTE: pflotran is doing a += here, but the array was zeroed out.
    ionx->set_concentration(ionx->X() * omega / primary_species[icomp].charge());
  }
}


void IonExchangeRxn::AddContributionToTotal(std::vector<double> *total)
{
  // pflotran: reaction.F90, function RTotalSorbEqIonX
  for (auto it = ionx_complexes_.begin(); it != ionx_complexes_.end(); ++it) {
    int icomp = it->primary_id();
    (*total)[icomp] += it->concentration();
  }
}


void IonExchangeRxn::AddContributionToDTotal(
    const std::vector<Species>& primary_species,
    MatrixBlock* dtotal)
{
  // pflotran: reaction.F90, function RTotalSorbEqIonX

  // sum up charges
  double sumZX = 0.0;
  for (auto ionx = ionx_complexes_.begin(); ionx != ionx_complexes_.end(); ionx++) {
    int icomp = ionx->primary_id();
    sumZX += primary_species[icomp].charge() * ionx->X();
  }

  // add contribution to derivatives
  for (auto ionx = ionx_complexes_.begin(); ionx != ionx_complexes_.end(); ionx++) {
    int icomp = ionx->primary_id();
    double temp = primary_species[icomp].charge() / sumZX;
    for (auto ionx2 = ionx_complexes_.begin(); ionx2 != ionx_complexes_.end(); ionx2++) {
      int jcomp = ionx2->primary_id();
      double value;
      if (ionx == ionx2) {
        value = ionx->concentration() * (1.0 - temp*ionx2->X()) / primary_species[jcomp].molality();
      } else {
        value = -ionx->concentration() * temp * ionx2->X() / primary_species[jcomp].molality();
      }
      dtotal->AddValue(icomp, jcomp, value);
    }
  }
}


void IonExchangeRxn::CheckUniformZ(const std::vector<Species>& primary_species)
{
  bool flag = true;
  for (int i = 0; i < ionx_complexes_.size(); i++) {
    for (int j = i+1; j < ionx_complexes_.size(); j++) {
      if (primary_species[ionx_complexes_[i].primary_id()].charge() != 
          primary_species[ionx_complexes_[j].primary_id()].charge()) {
        flag = false;
        i = ionx_complexes_.size();
        break;
      }
    }
  }
  uniform_z_ = flag;
  uniform_z_set_ = true;
}


void IonExchangeRxn::DisplaySite(const Teuchos::Ptr<VerboseObject> vo) const
{
  for (auto site = ionx_site_.begin(); site != ionx_site_.end(); site++) {
    site->Display(vo);
  }
}


void IonExchangeRxn::Display(const Teuchos::Ptr<VerboseObject> vo) const
{
  for (auto it = ionx_complexes_.begin(); it != ionx_complexes_.end(); ++it) {
    it->Display(vo);
  }
}


void IonExchangeRxn::DisplayResultsHeader() const {
  std::cout << std::setw(15) << "---" << std::endl;
}


void IonExchangeRxn::DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const
{
  ionx_site_[0].DisplayResultsHeader(vo);
  for (auto site = ionx_site_.begin(); site != ionx_site_.end(); site++) {
    site->DisplayResults(vo);
  }

  ionx_complexes_[0].DisplayResultsHeader(vo);
  for (auto it = ionx_complexes_.begin(); it != ionx_complexes_.end(); ++it) {
    it->DisplayResults(vo);
  }
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
