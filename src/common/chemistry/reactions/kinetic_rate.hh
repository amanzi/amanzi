/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Abstract base class for all kinetic rates
*/

#ifndef AMANZI_CHEMISTRY_KINETIC_RATE_HH_
#define AMANZI_CHEMISTRY_KINETIC_RATE_HH_

#include <vector>
#include <string>

#include "VerboseObject.hh"

#include "species.hh"
#include "secondary_species.hh"
#include "mineral.hh"
#include "string_tokenizer.hh"

namespace Amanzi {
namespace AmanziChemistry {

class MatrixBlock;

class KineticRate {
 public:
  virtual ~KineticRate(void) {};

  virtual void Setup(const SecondarySpecies& reaction,
                     const StringTokenizer& reaction_data,
                     const SpeciesArray& primary_species) = 0;
  virtual void Update(const SpeciesArray& primary_species,
                      const std::vector<Mineral>& minerals) = 0;
  virtual void AddContributionToResidual(const std::vector<Mineral>& minerals,
                                         const double bulk_volume,
                                         std::vector<double> *residual) = 0;
  virtual void AddContributionToJacobian(const SpeciesArray& primary_species,
                                         const std::vector<Mineral>& minerals,
                                         const double bulk_volume,
                                         MatrixBlock* J) = 0;
  virtual void Display(const Teuchos::Ptr<VerboseObject> vo) const = 0;

  virtual void ParseParameters(const StringTokenizer& rate_parameters) = 0;

  void SetSpeciesIds(const SpeciesArray& species,
                     const std::string& species_type,
                     const std::vector<SpeciesName>& in_names,
                     const std::vector<double>& in_stoichiometry,
                     std::vector<SpeciesId>* out_ids,
                     std::vector<double>* out_stoichiometry);

  void DisplayReaction(const Teuchos::Ptr<VerboseObject> vo) const;

  void set_debug(const bool value) {
    this->debug_ = value;
  };

  bool debug(void) const {
    return this->debug_;
  };
  std::string name(void) const {
    return this->name_;
  };
  SpeciesId identifier(void) const {
    return this->identifier_;
  };

  double reaction_rate(void) const {
    return this->reaction_rate_;
  }

 protected:
  KineticRate(void);

  void set_name(const std::string in_name) {
    this->name_ = in_name;
  };
  void set_identifier(const SpeciesId in_id) {
    this->identifier_ = in_id;
  };

  void set_reaction_rate(const double rate) {
    this->reaction_rate_ = rate;
  };

  std::vector<SpeciesName> reactant_names;
  std::vector<double> reactant_stoichiometry;
  std::vector<SpeciesId> reactant_ids;

 private:
  bool debug_;
  std::string name_;
  SpeciesId identifier_;  // the index identifier of the associated mineral!
  double reaction_rate_;  // volumetric rate: [moles/sec/m^3 bulk]
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif

