/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Amanzi

License:
Author: Markus Berndt
        Ethan Coon (ecoon@lanl.gov)

Observable data object

------------------------------------------------------------------------- */

#ifndef AMANZI_OBSERVABLE_HH_
#define AMANZI_OBSERVABLE_HH_

#include "AmanziTypes.hh"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "IOEvent.hh"
#include "MeshDefs.hh"
#include "ObservationData.hh"

namespace Amanzi {

class State;

double ObservableSum(double a, double b);
double ObservableMin(double a, double b);
double ObservableMax(double a, double b);

class Observable : public IOEvent {

public:
  Observable(Teuchos::ParameterList &plist, Epetra_MpiComm *comm);

  std::string name() { return name_; }
  std::string variable() { return variable_; }

  // DO NOT OVERRIDE -- instead, use the virtual, protected version.
  void Update(const State &S,
              Amanzi::ObservationData::DataTriple &data_triplet);

  void Flush();

protected:
  virtual void Update_(const State &S,
                       Amanzi::ObservationData::DataTriple &data_triplet);

  virtual void WriteHeader_();

protected:
  bool write_;
  int interval_;
  int count_;

  std::string filenamebase_;
  Teuchos::RCP<std::ofstream> out_;

  bool flux_normalize_;
  std::string name_;
  std::string variable_;
  std::string region_;
  std::string functional_;
  std::string location_;
  std::string delimiter_;

  double (*function_)(double a, double b, double vol);
};

} // namespace Amanzi

#endif
