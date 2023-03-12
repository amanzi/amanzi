/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*!

A user may request any number of specific observations from Amanzi.
Each labeled observation data quantity involves a field quantity, a model, a region from which it will extract its source data, and a list of discrete times 
for its evaluation.  The observations are evaluated during the simulation and returned to the calling process through one of Amanzi arguments.

* `"observation data`" [list] can accept multiple lists for named observations.

  * `"observation output filename`" [string] user-defined name for the file that the observations are written to.
    The file name can contain relative or absolute path to an *existing* directory only. 

  * `"time unit`" [string] defines time unit for output data.
    Available options are `"s`", `"h`", `"d`", and `"y`". Default is `"s`".

  * `"mass unit`" [string] defines mass unit for output data. 
    Available options are `"g`", `"lb`", and `"lb`". Default is `"kg`".

  * `"length unit`" [string] defines length unit for output data.
     Available options are `"cm`", `"in`", `"ft`", `"yd`" , `"m`", and `"km`". Default is `"m`".

  * `"concentration unit`" [string] defines concentration unit for output data.
     Available options are `"molar`", and `"SI`". Default is `"molar`".

  * `"precision`" [int] defines the number of significant digits. Default is 16.

  * OBSERVATION [list] user-defined label, can accept values for `"variables`", `"functional`",
    `"region`", `"times`", and TSPS (see below).

    * `"domain name`" [string] name of the domain. Typically, it is either `"domain`" for
      the matrix/subsurface or `"fracture`" for the fracture network.

    * `"variable`" [string] a list of field quantities taken from the list of 
      available field quantities:

      * volumetric water content [-] (volume water / bulk volume)
      * aqueous saturation [-] (volume water / volume pore space)
      * aqueous pressure [Pa]
      * hydraulic head [m] 
      * permeability-weighted hydraulic head [m] 
      * drawdown [m] 
      * permeability-weighted drawdown [m] 
      * volumetric water content [-]
      * gravimetric water content [-]
      * water table [m]
      * SOLUTE aqueous concentration [mol/m^3]
      * SOLUTE gaseous concentration [mol/m^3]
      * SOLUTE sorbed concentration [mol/kg] 
      * SOLUTE free ion concentration
      * x-, y-, z- aqueous volumetric flux [m/s]
      * material id [-]
      * aqueous mass flow rate [kg/s] (when funtional="integral")
      * aqueous volumetric flow rate [m^3/s] (when functional="integral")
      * fractures aqueous volumetric flow rate [m^3/s] (when functional="integral")
      * SOLUTE volumetric flow rate [mol/s] (when functional="integral")
      * SOLUTE breakthrough curve [mol] (when functional="integral")
      * pH [-] 
      * centroid x [m]

    Observations *drawdown* and *permeability-weighted* are calculated with respect to the value 
    registered at the first time it was requested.

    The following observations are point-type obervations: "water table", "drawdown".

    The following observations are integrated continuously in time but saved only at specified
    times: "SOLUTE breakthrough curve". 

    * `"functional`" [string] the label of a function to apply to each of the variables
      in the variable list (Function options detailed below)

    * `"region`" [string] the label of a user-defined region

    * `"cycles start period stop`" [Array(int)] the first entry is the start cycle,
      the second is the cycle period, and the third is the stop cycle or -1 in which case 
      there is no stop cycle. A visualization dump shall be written at such cycles that 
      satisfy cycle = start + n*period, for n=0,1,2,... and cycle < stop if stop != -1.0.

    * `"cycles start period stop n`" [Array(int)] if multiple cycles start-period-stop
      parameters are needed, then use these parameters with n=0,1,2,..., and not the single 
      `"cycles start period stop`" parameter.

    * `"cycles`" [Array(int)] an array of discrete cycles that at which a visualization dump shall be written. 

    * `"times start period stop`" [Array(double)] the first entry is the start time,
      the second is the time period, and the third is the stop time or -1 in which case 
      there is no stop time. A visualization dump shall be written at such times that 
      satisfy time = start + n*period, for n=0,1,2,... and time < stop if stop != -1.0.

    * `"times start period stop n`" [Array(double) if multiple start-period-stop parameters
      are needed, then use this these parameters with n=0,1,2,..., and not the 
      single  `"times start period stop`" parameter.

    * `"times`" [Array(double)] an array of discrete times that at which a visualization
      dump shall be written.

    * `"delimiter`" [string] the string used to delimit columns in the observation file
      output, default is `",`".

    * `"interpolation`" [string] the string which defines
      the interpolation method to compute observation. Works ONLY with
      Line Segment region at the moment. Available options `"linear`"
      and `"constant`". Default is `"linear`" 

    * `"weighting`"  [string] the string defined the weighting
      function applied to compute observation. Works ONLY with
      Line Segment region at the moment. Available options `"flux norm`"
      and `"none`". Default is `"none`". `"flux norm`" is the absolute
      value of the Darcy flux in a cell.

The following observation functionals are currently supported.
All of them operate on the variables identified.

* `"observation data: point`" returns the value of the field quantity at a point.

* `"observation data: integral`" returns the integral of the field quantity over the region specified.

* `"observation data: extensive integral`" returns the integral of an extensive variable
  over the region specified.  Note that this should be used over the above Integral when 
  the variable to be integrated is an extensive quantity, i.e. water content or flux.

* `"observation data: minimum`" and `"observation data: maximum`" returns the minimum 
  (respectively maximum) of the field quantity over the region specified.

.. code-block:: xml

  <ParameterList>  <!-- parent list -->
  <ParameterList name="observation data">
    <Parameter name="observation output filename" type="string" value="_OUTPUT.out"/>
    <Parameter name="precision" type="int" value="10"/>
    <ParameterList name="_ANY OBSERVATION NAME">
      <Parameter name="region" type="string" value="_REGION"/>
      <Parameter name="functional" type="string" value="observation data: point"/>
      <Parameter name="variable" type="string" value="volumetric water content"/>

      <Parameter name="cycles" type="Array(int)" value="{100000, 200000, 400000, 500000}"/>
      <Parameter name="cycles start period stop" type="Array(int)" value="{0, 100, -1}"/>

      <Parameter name="times" type="Array(double)" value="{101.0, 303.0, 422.0}"/>
      <Parameter name="times start period stop 0" type="Array(double)" value="{0.0, 10.0, 100.0}"/>
      <Parameter name="times start period stop 1" type="Array(double)" value="{100.0, 25.0, -1.0}"/>
    </ParameterList>

    <ParameterList name="_ANY OBSERVATION NAME B">  <!-- another observation -->
      ...
    </ParameterList>
  </ParameterList>
  </ParameterList>

In this example, we collect `"volumetric water content`" on four selected cycles,
every 100 cycles, three selected times, every 10 seconds from 0 to 100, and every 25 seconds after that.

*/

#ifndef AMANZI_FLEXIBLE_OBSERVATIONS_HH_
#define AMANZI_FLEXIBLE_OBSERVATIONS_HH_

#include <map>

// TPLs
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

// Amanzi
#include "IOEvent.hh"
#include "ObservationData.hh"
#include "State.hh"
#include "TimeStepManager.hh"
#include "Units.hh"
#include "VerboseObject.hh"
#include "ObservableAmanzi.hh"
#include "Key.hh"

namespace Amanzi {

class FlexibleObservations {
 public:
  FlexibleObservations(Teuchos::RCP<Teuchos::ParameterList> coordinator_list,
                       Teuchos::RCP<Teuchos::ParameterList> obs_list,
                       Teuchos::RCP<Teuchos::ParameterList> units_list,
                       Amanzi::ObservationData& observation_data,
                       Teuchos::RCP<const State> S);

  ~FlexibleObservations()
  {
    if (vo_ != NULL) delete vo_;
  }

  void RegisterComponentNames(std::vector<std::string>& comp_names,
                              std::vector<double>& comp_mol_masses,
                              int num_liquid)
  {
    comp_names_ = comp_names;
    comp_mol_masses_ = comp_mol_masses;
    num_liquid_ = num_liquid;
  }

  int MakeObservations(State& S);
  int MakeContinuousObservations(State& S);

  bool DumpRequested(const int);
  bool DumpRequested(const double);
  bool DumpRequested(const int, const double);

  void RegisterWithTimeStepManager(const Teuchos::Ptr<TimeStepManager>& tsm);

  void FlushObservations();

 private:
  double CalculateWaterTable_(State& S, AmanziMesh::Entity_ID_List& ids);

 protected:
  VerboseObject* vo_;

 private:
  int rank_;
  Teuchos::RCP<Teuchos::ParameterList> obs_list_;
  Teuchos::RCP<Teuchos::ParameterList> coordinator_list_;
  Amanzi::ObservationData& observation_data_;
  std::map<std::string, Teuchos::RCP<Observable>> observations;

  std::vector<std::string> comp_names_;
  std::vector<double> comp_mol_masses_;
  int num_liquid_;

  Utils::Units units_;
};

} // namespace Amanzi

#endif
