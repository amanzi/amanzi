/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! Collects, reduces, and writes observations during a simulation.

/*
  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.
*/

/*!

Observations are a localized-in-space but frequent-in-time view of simulation
output, designed to get at useful diagnostic quantities such as hydrographs,
total water content, quantities at a point, etc.  These allow frequent
collection in time without saving huge numbers of visualization files to do
postprocessing.  In fact, these should be though of as orthogonal data queries
to visualization -- vis is pointwise in time but complete in space, while
observations are pointwise/finite in space but complete in time.

Observables describe what is observed -- the variable, region,
integration/averaging scheme, etc.

A user may request any number of specific observables.  Each observable spec
involves a field quantity, a functional reduction operator, a region from which
it will extract its source data, and a list of discrete times for its
evaluation.  The observations are evaluated during the simulation and written
to disk by the UnstructuredObservation_ object.

.. _observable-spec:
.. admonition:: observable-spec

    * `"variable`" ``[string]`` any ATS variable used by any PK, e.g. `"pressure`"
      or `"surface-water_content`"

    * `"region`" ``[string]`` the label of a user-defined region

    * `"location name`" ``[string]`` **cell** The mesh location of the thing to
      be measured, i.e. `"cell`", `"face`", or `"node`"

    * `"number of vectors`" ``[int]`` **1** Number of vector components to write.

    * `"functional`" ``[string]`` the type of function to apply to the variable
      on the region.  One of:

      - `"point`" returns the value of the field quantity at a
        point.  The region and location name should result in a single entity being
        selected.

      - `"average`" returns the volume-weighted average of
        the field across all entities in the region selected.  This is likely
        what you want for intensive state variables such as `"temperature`" or
        `"saturation_liquid`".

      - `"integral`" returns the volume-weighted sum of a
        variable over the region.  This should be used for example on intensive
        sources, for instance `"surface-precipitation`", to get the total
        source/sink.

      - `"extensive integral`" returns the sum of an variable
        over the region.  This should be used for extensive quantities such as
        `"water_content`" or `"energy`" which already include the volume in
        their value.

      - `"minimum`" returns the min value over the region

      - `"maximum`" returns the max value over the region

    * `"direction normalized flux`" ``[bool]`` **false** For flux observations,
      dots the face-normal flux with a vector to ensure fluxes are integrated
      pointing the same direction.

    * `"direction normalized flux direction`" ``[Array(double)]`` **optional** For
      flux observations, provides the vector to dot the face normal with.  If this
      is not provided, then it is assumed that the faces integrated over are all
      boundary faces and that the default vector is the outward normal direction
      for each face.

    * `"time integrated`" ``[bool]`` **false** If true, observe the
      time-integral, observing on all cycles and accumulating the
      backwards-Euler product of dt times the observable at the new time.

*/

#ifndef AMANZI_OBSERVABLE_HH_
#define AMANZI_OBSERVABLE_HH_

#include <limits>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "Point.hh"
#include "MeshDefs.hh"
#include "IOEvent.hh"

namespace Amanzi {

class State;

class Observable {
 public:
  static const double nan;

  Observable(Teuchos::ParameterList& plist);

  const std::string& get_name() { return name_; }
  const std::string& get_variable() { return variable_; }
  const std::string& get_region() { return region_; }
  const std::string& get_location() { return location_; }
  const std::string& get_functional() { return functional_; }
  bool is_time_integrated() { return time_integrated_; }
  int get_num_vectors() { return num_vectors_; }

  void Setup(const Teuchos::Ptr<State>& S);
  void FinalizeStructure(const Teuchos::Ptr<State>& S);
  void Update(const Teuchos::Ptr<State>& S, std::vector<double>& data, int start_loc);

 protected:
  bool flux_normalize_;
  Teuchos::RCP<AmanziGeometry::Point> direction_;

  std::string name_;
  std::string variable_;
  std::string region_;
  std::string functional_;
  std::string location_;
  int num_vectors_;
  bool time_integrated_;
  double old_time_;
  bool has_eval_;

  double (*function_)(double a, double b, double vol);
};


} // namespace

#endif
