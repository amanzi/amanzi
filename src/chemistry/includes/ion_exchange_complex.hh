/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_IONEXCHANGECOMPLEX_HH_
#define AMANZI_CHEMISTRY_IONEXCHANGECOMPLEX_HH_

/*
**  Class for ion exchange complexation reaction
**
**  NaX <===> Na+ + X-
**
*/

#include <vector>

#include "species.hh"

namespace amanzi {
namespace chemistry {

typedef std::string IonxComplexName;
typedef int IonxComplexId; 

class IonExchangeComplex {
 public:
  IonExchangeComplex();
  IonExchangeComplex(const IonxComplexName,
                     const IonxComplexId complex_id,
                     const SpeciesName primary_name,
                     const SpeciesId primary_id,
                     const double K);
  virtual ~IonExchangeComplex();

  void display(void) const;
  void Display(void) const;
  void DisplayReaction(void) const;
  void DisplayResultsHeader(void) const;
  void DisplayResults(void) const;

  std::string name(void) const {
    return name_;
  };
  std::string primary_name(void) const {
    return primary_name_;
  };
  int primary_id(void) const {
    return primary_id_;
  };
  double K(void) const {
    return K_;
  };
  double X(void) const {
    return X_;
  };
  double concentration(void) const {
    return concentration_;
  };

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

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_IONEXCHANGECOMPLEX_HH_
