/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Class for ion exchange complexation reaction

  NaX <===> Na+ + X-
*/

#ifndef AMANZI_CHEMISTRY_IONEXCHANGECOMPLEX_HH_
#define AMANZI_CHEMISTRY_IONEXCHANGECOMPLEX_HH_

#include <vector>

#include "species.hh"

namespace Amanzi {
namespace AmanziChemistry {

typedef std::string IonxComplexName;
typedef int IonxComplexId; 

class IonExchangeComplex {
 public:
  IonExchangeComplex() {};
  IonExchangeComplex(const IonxComplexName,
                     const IonxComplexId complex_id,
                     const SpeciesName primary_name,
                     const SpeciesId primary_id,
                     const double K);
  virtual ~IonExchangeComplex() {};

  void display(const Teuchos::Ptr<VerboseObject> vo) const;
  void Display(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayReaction(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const;

  std::string name(void) const { return name_; };
  std::string primary_name(void) const { return primary_name_; };

  int primary_id(void) const { return primary_id_; };
  double K(void) const { return K_; };
  double X(void) const { return X_; };
  double concentration(void) const { return concentration_; };

  void set_X(const double d) { X_ = d; };
  void set_concentration(const double d) { concentration_ = d; };

 private:
  IonxComplexName name_;
  IonxComplexId id_;
  SpeciesName primary_name_;
  SpeciesId primary_id_;

  double concentration_;
  double K_;
  double X_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
