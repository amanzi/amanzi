/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*!

A user may request periodic dumps of Walkabout data. Output controls for Walkabout data are 
limited to file name generation and writing frequency, by numerical cycle number or time.

* `"walkabout data`" [list] can accept a file name base [string] and cycle data [list] 
  used to generate the file base name or directory base name that is used in writing Walkabout data. 

  * `"file name base`" [string] The file name can contain relative or absolute path to an *existing* 
    directory only.  Default is `"walkabout`".
  
  * `"file name digits`" [int] specify the number of digits that should be appended to the file 
    name for the cycle number. Default is 5.

  * `"cycles start period stop`" [Array(int)] the first entry is the start cycle, 
    the second is the cycle period, and the third is the stop cycle or -1 in which case 
    there is no stop cycle. A visualization dump shall be written at such cycles that 
    satisfy cycle = start + n*period, for n=0,1,2,... and cycle < stop if stop != -1.0.

  * `"cycles start period stop n`" [Array(int)] if multiple cycles start-period-stop parameters 
    are needed, then use these parameters with n=0,1,2,..., and not the single
    `"cycles start period stop`" parameter.

  * `"cycles`" [Array(int)] an array of discrete cycles that at which a visualization dump shall be written. 

  * `"times start period stop`" [Array(double)] the first entry is the start time, 
    the second is the time period, and the third is the stop time or -1 in which case 
    there is no stop time. A visualization dump shall be written at such times that 
    satisfy time = start + n*period, for n=0,1,2,... and time < stop if stop != -1.0.

  * `"times start period stop n`" [Array(double) if multiple start-period-stop parameters 
    are needed, then use this these parameters with n=0,1,2,..., and not the single
    `"times start period stop`" parameter.

  * `"times`" [Array(double)] an array of discrete times that at which a visualization dump shall be written.

  * `"write regions`" [list] contains three lists of equal size with region names,
    material names, and material ids to write into the output file.

    * `"region names`" [Array(string)] specifies names of regions.
    * `"material names`" [Array(int)] specifies names of materials. 
    * `"material ids`" [Array(int)] specifies material ids. 

.. code-block:: xml

  <ParameterList>  <!-- parent list -->
  <ParameterList name="walkabout data">
    <Parameter name="file name base" type="string" value="_WALKABOUT"/>
    <Parameter name="file name digits" type="int" value="5"/>

    <Parameter name="cycles start period stop" type="Array(int)" value="{0, 100, -1}"/>
    <Parameter name="cycles" type="Array(int)" value="{999, 1001}"/>

    <Parameter name="times start period stop 0" type="Array(double)" value="{0.0, 10.0, 100.0}"/>
    <Parameter name="times start period stop 1" type="Array(double)" value="{100.0, 25.0, -1.0}"/>
    <Parameter name="times" type="Array(double)" value="{101.0, 303.0, 422.0}"/>

    <ParameterList name="write regions">
      <Parameter name="region names" type="Array(string)" value="{_REGION1, _REGION2}"/>
      <Parameter name="material names" type="Array(string)" value="{_MAT1, _MAT2}"/>
      <Parameter name="material ids" type="Array(int)" value="{1000, 2000}"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>

In this example, walkabout data files are written when the cycle number is 
a multiple of 100.

*/

#ifndef AMANZI_WALKABOUT_CHECKPOINT_HH_
#define AMANZI_WALKABOUT_CHECKPOINT_HH_

#include <map>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "PK.hh"
#include "State.hh"
#include "Checkpoint.hh"

#include "TimeStepManager.hh"
#include "VerboseObject.hh"


namespace Amanzi {

class WalkaboutCheckpoint : public Checkpoint {
 public:
  WalkaboutCheckpoint(Teuchos::ParameterList& plist, const State& S) : Checkpoint(plist, S){};
  WalkaboutCheckpoint() : Checkpoint(true){};

  // output of fields
  void WriteDataFile(Teuchos::RCP<State>& S, Teuchos::RCP<PK> pk);

  // recontruct vector velocity at mesh nodes
  void CalculateDarcyVelocity(Teuchos::RCP<State>& S,
                              std::vector<AmanziGeometry::Point>& xyz,
                              std::vector<AmanziGeometry::Point>& velocity) const;

  // interpolate various fileds to mesh nodes
  void CalculateData(Teuchos::RCP<State>& S,
                     std::vector<AmanziGeometry::Point>& xyz,
                     std::vector<AmanziGeometry::Point>& velocity,
                     std::vector<double>& porosity,
                     std::vector<double>& saturation,
                     std::vector<double>& pressure,
                     std::vector<double>& isotherm_kd,
                     std::vector<int>& material_ids);

 private:
  Teuchos::RCP<PK> pk_;
};

} // namespace Amanzi

#endif
