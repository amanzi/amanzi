/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_IONEXCHANGERXN_HH_
#define AMANZI_CHEMISTRY_IONEXCHANGERXN_HH_

// Class for ion exchange reaction

#include <vector>

#include "ion_exchange_complex.hh"
#include "ion_exchange_site.hh"
#include "species.hh"

namespace amanzi {
namespace chemistry {

// forward declarations from chemistry
class Block;

class IonExchangeRxn {
 public:
  IonExchangeRxn();
  IonExchangeRxn(IonExchangeSite* ionx_sites,
                 const std::vector<IonExchangeComplex>& ionx_complexes);
  explicit IonExchangeRxn(IonExchangeSite ionx_sites);
  ~IonExchangeRxn();

  // add complexes to the reaction
  void AddIonExchangeComplex(const IonExchangeComplex& complex);
  void AddIonExchangeSite(const IonExchangeSite& site);

  // update sorbed concentrations
  void Update(const std::vector<Species>& primarySpecies);
  // add stoichiometric contribution of complex to sorbed total
  void AddContributionToTotal(std::vector<double> *total);
  // add derivative of total with respect to free-ion to sorbed dtotal
  void AddContributionToDTotal(const std::vector<Species>& primarySpecies,
                               Block* dtotal);
  void CheckUniformZ(const std::vector<Species>& primarySpecies);
  void set_uniform_z(const bool flag) { uniform_z_ = flag; };
  bool uniform_z(void) const { return uniform_z_; };
  bool uniform_z_set(void) const { return uniform_z_set_; };
  IonExchangeSite site(void) const { return ionx_site_[0]; };
  std::vector<IonExchangeSite>  ionx_sites(void) const { return ionx_site_; };
  std::vector<IonExchangeComplex>  ionx_complexes(void) const { return ionx_complexes_; };
  void set_cation_exchange_capacity(const double d) { ionx_site_[0].set_cation_exchange_capacity(d); };

  void display(void) const;
  void Display(void) const;
  void DisplaySite(void) const;
  void DisplayComplexes(void) const;
  void DisplayResultsHeader(void) const;
  void DisplayResults(void) const;

 protected:

 private:
  std::vector<IonExchangeComplex> ionx_complexes_;
  std::vector<IonExchangeSite> ionx_site_;
  bool uniform_z_;
  bool uniform_z_set_;
  double ref_cation_sorbed_conc_;

};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_IONEXCHANGERXN_HH_
