/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

/* -------------------------------------------------------------------------
Amanzi

License:
Author: Markus Berndt
        Ethan Coon (coonet@ornl.gov)

Observable data object

------------------------------------------------------------------------- */

#ifndef AMANZI_OBSERVABLE_HH_
#define AMANZI_OBSERVABLE_HH_

#include "Epetra_MpiComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "IOEvent.hh"
#include "MeshDefs.hh"
#include "ObservationData.hh"

namespace Amanzi {

class State;

double
ObservableSum(double a, double b);
double
ObservableMin(double a, double b);
double
ObservableMax(double a, double b);

class Observable : public IOEvent {
 public:
  Observable(const Teuchos::RCP<Teuchos::ParameterList>& plist,
             const Comm_ptr_type& comm);

  std::string name() { return name_; }
  std::string variable() { return variable_; }

  // DO NOT OVERRIDE -- instead, use the virtual, protected version.
  void
  Update(const State& S, Amanzi::ObservationData::DataQuadruple& data_triplet);

  void Flush();

 protected:
  virtual void
  Update_(const State& S, Amanzi::ObservationData::DataQuadruple& data_triplet);

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
