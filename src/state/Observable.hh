/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! Observable: Collects, reduces, and writes observations during a simulation.

/*
  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
*/

/*!  
Observations are a localized-in-space but frequent in time view of
data, designed to get at useful diagnostic quantities such as
hydrographs, total water content, quantities at a point, etc.  These
are designed to allow frequent collection in time without saving huge
numbers of visualization files to do postprocessing.  In fact, these
should be though of as orthogonal data queries to visualization -- vis
is pointwise in time but complete in space, while observations are
pointwise/finite in space but complete in time.

A user may request any number of specific observations from ATS.  Each
observation spec involves a field quantity, a functional reduction
operator, a region from which it will extract its source data, and a
list of discrete times for its evaluation.  The observations are
evaluated during the simulation and written to disk.

* `"observations`" [list] can accept multiple ``[observation-spec]`` entries.

An ``[observation-spec]`` consists of the following quantities:

* `"observation output filename`" [string] user-defined name for the file that the observations are written to.

* `"variable`" [string] any ATS variable used by any PK, e.g. `"pressure`" or `"surface-water_content`"

* `"region`" [string] the label of a user-defined region

* `"location name`" [string] the mesh location of the thing to be measured, i.e. `"cell`", `"face`", or `"node`"

* `"functional`" [string] the label of a function to apply to the variable across the region.  Valid functionals include:
 * `"observation data: point`" returns the value of the field quantity at a point.  The region and location name must result in a single entity being selected.
 * `"observation data: extensive integral`" returns the sum of an (extensive) variable over the region.  This should be used for extensive quantities such as `"water_content`" or `"energy`".
 * `"observation data: intensive integral`" returns the volume-weighted average of an (intensive) variable over the region.  This should be used for intensive quantities such as `"temperature`" or `"saturation_liquid`".

* Additionally, each ``[observation-spec]`` contains all parameters as in a IOEvent_ spec, which are used to specify at which times/cycles the observation is collected.

For flux observations, and additional option is available:

* `"direction normalized flux`" [bool] *false* Normalize the flux to point in the outward-normal direction.  This is important when looking at fluxes across a boundary, for instance to plot a hydrograph.


Example:

.. code-block:: xml
  
  <ParameterList name="observations" type="ParameterList">
    <!-- This measures the hydrograph out the "east" face of the surface domain -->
    <ParameterList name="surface outlet flux" type="ParameterList">
      <Parameter name="variable" type="string" value="surface-mass_flux" />
      <Parameter name="direction normalized flux" type="bool" value="true" />
      <Parameter name="region" type="string" value="east" />
      <Parameter name="functional" type="string" value="observation data: extensive integral" />
      <Parameter name="delimiter" type="string" value=" " />
      <Parameter name="location name" type="string" value="face" />
      <Parameter name="observation output filename" type="string" value="surface_outlet_flux.dat" />
      <Parameter name="times start period stop" type="Array(double)" value="{0.0,86400.0,-1.0}" />
    </ParameterList>
    <!-- This measures the total water, in mols, in the entire subsurface domain -->
    <ParameterList name="subsurface water content" type="ParameterList">
      <Parameter name="variable" type="string" value="water_content" />
      <Parameter name="region" type="string" value="computational domain" />
      <Parameter name="functional" type="string" value="observation data: extensive integral" />
      <Parameter name="delimiter" type="string" value=" " />
      <Parameter name="location name" type="string" value="cell" />
      <Parameter name="observation output filename" type="string" value="water_content.dat" />
      <Parameter name="times start period stop" type="Array(double)" value="{0.0,86400.0,-1.0}" />
    </ParameterList>
    <!-- This tracks the temperature at a point -->
    <ParameterList name="temperature_probeA" type="ParameterList">
      <Parameter name="variable" type="string" value="temperature" />
      <Parameter name="region" type="string" value="probeA" />
      <Parameter name="functional" type="string" value="observation data: point" />
      <Parameter name="delimiter" type="string" value=" " />
      <Parameter name="location name" type="string" value="cell" />
      <Parameter name="observation output filename" type="string" value="temperature_probeA.dat" />
      <Parameter name="times start period stop" type="Array(double)" value="{0.0,86400.0,-1.0}" />
    </ParameterList>
  </ParameterList>

*/

#ifndef AMANZI_OBSERVABLE_HH_
#define AMANZI_OBSERVABLE_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_MpiComm.h"

#include "ObservationData.hh"
#include "MeshDefs.hh"
#include "IOEvent.hh"

namespace Amanzi {

class State;

double ObservableSum(double a, double b);
double ObservableMin(double a, double b);
double ObservableMax(double a, double b);

class Observable : public IOEvent {

 public:

  Observable(Teuchos::ParameterList& plist, Epetra_MpiComm *comm);

  std::string name() { return name_; }
  std::string variable() { return variable_; }

  // DO NOT OVERRIDE -- instead, use the virtual, protected version.
  void Update(const State& S,
              Amanzi::ObservationData::DataQuadruple& data_quad);

  void Flush();

 protected:
  virtual void Update_(const State& S,
                       Amanzi::ObservationData::DataQuadruple& data_quad);

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


} // namespace

#endif
