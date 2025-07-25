/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)
*/

/*!
Observations are a localized-in-space but frequent-in-time view of simulation
output, designed to get at useful diagnostic quantities such as hydrographs,
total water content, quantities at a point, etc.  These allow frequent
collection in time without saving huge numbers of visualization files to do
postprocessing.  In fact, these should be though of as orthogonal data queries
to visualization -- visualization is pointwise in time but complete in space,
while observations are pointwise/finite in space but complete in time.

.. _observation-spec:
.. admonition:: observation-spec

   * `"observation output filename`" ``[string]`` user-defined name for the file
     that the observations are written to.
   * `"delimiter`" ``[string]`` **COMMA** Delimiter to split columns of the file
   *  `"write interval`" ``[int]`` **1** Interval of observations on which to flush IO files.
   * `"time units`" ``[string]`` **s** Controls the unit of the time column in the observation file.
   * `"domain`" ``[string]`` **"domain"** Can be used to set the communicator which writes (defaults to the standard subsurface domain).
   * `"observed quantities`" ``[observable-spec-list]`` A list of Observable_
     objects that are all put in the same file.

   INCLUDES:

   - ``[io-event-spec]`` An :ref:`IOEvent` spec.

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

class UnstructuredObservations : public Utils::IOEvent {
 public:
  UnstructuredObservations(Teuchos::ParameterList& obs_list);

  void Setup(const Teuchos::Ptr<State>& S);
  void MakeObservations(const Teuchos::Ptr<State>& S);
  void Flush();

 private:
  // NOTE: this should only be called on one rank
  void InitFile_(const Teuchos::Ptr<const State>& S);
  void Write_(double time, const std::vector<double>& obs);

 private:
  std::string writing_domain_;
  Comm_ptr_type comm_;
  bool write_;

  std::vector<Teuchos::RCP<Observable>> observables_;

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


// non-member functions
bool portable_name(const std::string& name);
bool portable_filename(const std::string& filename);


} // namespace Amanzi


#endif
