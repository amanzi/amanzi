/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)
*/

//!  Collection of Observations on an unstructured mesh to be written to a common file.

/*!


.. _observation-spec:
.. admonition:: observation-spec

    * `"observation output filename`" ``[string]`` user-defined name for the file
      that the observations are written to.
    * `"delimiter`" ``[string]`` **\, ** Delimiter to split columns of the file
    *  `"write interval`" ``[int]`` **1** Interval of observations on which to flush IO files.
    * `"time units`" ``[string]`` **s** Controls the unit of the time column in the observation file.
    * `"domain`" ``[string]`` **""** Can be used to set the communicator which writes (usually not needed).
    * `"observed quantities`" ``[observable-spec-list]`` A list of Observable_
       objects that are all put in the same file.

    INCLUDES:
    - ``[io-event-spec]`` An IOEvent_ spec

Note, for backwards compatibility, an ``observable-spec`` may be directly
included within the `observation-spec` if it is the only variable to be
observed in this file.

*/

#ifndef AMANZI_UNSTRUCTURED_OBSERVATIONS_HH_
#define AMANZI_UNSTRUCTURED_OBSERVATIONS_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "IOEvent.hh"
#include "ObservationData.hh"
#include "Observable.hh"

#include "TimeStepManager.hh"

namespace Amanzi {

class UnstructuredObservations : public IOEvent {

 public:
  UnstructuredObservations(Teuchos::ParameterList& obs_list);

  void Setup(const Teuchos::Ptr<State>& S);
  void MakeObservations(const Teuchos::Ptr<State>& S);
  void Flush();

 private:
  // NOTE: this should only be called on one rank
  void InitFile_();
  void Write_(double time, const std::vector<double>& obs);

 private:
  std::vector<Teuchos::RCP<Observable>> observables_;

  std::string writing_domain_;
  bool write_;
  int num_total_;
  std::string filename_;
  std::string delimiter_;
  int interval_;
  int count_;
  bool time_integrated_;
  std::vector<double> integrated_observation_;
  bool observed_once_;

  double time_unit_factor_;
  std::string time_unit_;

  std::unique_ptr<std::ofstream> fid_;
};

}


#endif
