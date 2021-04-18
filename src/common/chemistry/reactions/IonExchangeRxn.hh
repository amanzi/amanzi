/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Class for ion exchange reaction
*/

#ifndef AMANZI_CHEMISTRY_IONEXCHANGERXN_HH_
#define AMANZI_CHEMISTRY_IONEXCHANGERXN_HH_

#include <vector>

#include "IonExchangeComplex.hh"
#include "IonExchangeSite.hh"
#include "Species.hh"

namespace Amanzi {
namespace AmanziChemistry {

// forward declarations from chemistry
class MatrixBlock;

class IonExchangeRxn {
 public:
  IonExchangeRxn();
  IonExchangeRxn(IonExchangeSite* ionx_sites,
                 const std::vector<IonExchangeComplex>& ionx_complexes);
  explicit IonExchangeRxn(IonExchangeSite ionx_sites);
  ~IonExchangeRxn() {};

  // add complexes to the reaction
  void AddIonExchangeComplex(const IonExchangeComplex& complex);
  void AddIonExchangeSite(const IonExchangeSite& site);

  // update sorbed concentrations
  void Update(const std::vector<Species>& primary_species);
  // add stoichiometric contribution of complex to sorbed total
  void AddContributionToTotal(std::vector<double> *total);
  // add derivative of total with respect to free-ion to sorbed dtotal
  void AddContributionToDTotal(const std::vector<Species>& primary_species,
                               MatrixBlock* dtotal);
  void CheckUniformZ(const std::vector<Species>& primary_species);
  void set_uniform_z(const bool flag) { uniform_z_ = flag; };
  bool uniform_z() const { return uniform_z_; };
  bool uniform_z_set() const { return uniform_z_set_; };

  IonExchangeSite site() const { return ionx_site_[0]; };
  std::vector<IonExchangeSite> ionx_sites() const { return ionx_site_; };
  std::vector<IonExchangeComplex>  ionx_complexes() const { return ionx_complexes_; };

  void set_cation_exchange_capacity(const double d) { ionx_site_[0].set_cation_exchange_capacity(d); };

  void set_ref_cation_sorbed_conc(double value) { ref_cation_sorbed_conc_ = value; }
  double ref_cation_sorbed_conc() const { return ref_cation_sorbed_conc_; }

  void Display(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplaySite(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayComplexes(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResultsHeader() const;
  void DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const;

 private:
  std::vector<IonExchangeComplex> ionx_complexes_;
  std::vector<IonExchangeSite> ionx_site_;
  bool uniform_z_;
  bool uniform_z_set_;
  double ref_cation_sorbed_conc_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
