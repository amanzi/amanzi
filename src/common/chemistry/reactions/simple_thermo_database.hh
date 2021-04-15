/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
*/
 
#ifndef AMANZI_CHEMISTRY_SIMPLETHERMODATABASE_HH_
#define AMANZI_CHEMISTRY_SIMPLETHERMODATABASE_HH_

#include <vector>
#include <string>

#include "species.hh"
#include "beaker.hh"

namespace Amanzi {
namespace AmanziChemistry {

class SimpleThermoDatabase : public Beaker {
 public:
  SimpleThermoDatabase(Teuchos::RCP<Teuchos::ParameterList> plist,
                       Teuchos::RCP<VerboseObject> vo);
  virtual ~SimpleThermoDatabase() {};

  virtual void Initialize(const BeakerParameters& parameters);

 private:
  void ParseReaction_(const std::string& reactants,
                      const std::string& products,
                      std::vector<std::string>* species,
                      std::vector<double>* stoichiometries);

  void ParseReaction_(const std::string& reaction,
                      std::vector<std::string>* primary_name,
                      std::vector<double>* primary_stoichiometry,
                      std::vector<int>* primary_id,
                      std::string* surface_name,
                      double* surface_stoichiometry,
                      int* surface_id,
                      double* h2o_stoich);

 private:
  Teuchos::RCP<Teuchos::ParameterList> plist_;

  std::vector<SurfaceSite> surface_sites_;
  std::vector<SurfaceComplexationRxn> surface_complexation_reactions_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi

#endif
