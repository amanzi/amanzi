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

class IonExchangeComplex {
 public:
  IonExchangeComplex() {};
  IonExchangeComplex(const std::string& complex_name,
                     int complex_id,
                     const std::string& primary_name,
                     int primary_id,
                     double K);
  virtual ~IonExchangeComplex() {};

  void display(const Teuchos::Ptr<VerboseObject> vo) const;
  void Display(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayReaction(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const;

  // setters and getters
  std::string name() const { return name_; };
  std::string primary_name() const { return primary_name_; };

  int primary_id() const { return primary_id_; };
  double K() const { return K_; };
  double X() const { return X_; };
  double concentration() const { return concentration_; };

  void set_X(double d) { X_ = d; };
  void set_concentration(double d) { concentration_ = d; };

 private:
  std::string name_, primary_name_;
  int id_, primary_id_;

  double concentration_, K_, X_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
