/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
           Daniil Svyatskiy
*/

/*!

CycleDriver is a class to hold the cycle driver, which runs the overall, top level timestep
loop. It instantiates states, ensures they are initialized, and runs the timestep loop
including Vis and restart/checkpoint dumps.  It contains one and only one PK
-- most likely this PK is an MPC of some type -- to do the actual work.

The new multi-processor cycle driver provides more flexibility
to handle multiphysics process kernels (PKs) and multiple time periods.

* `"component names`" [Array(string)] provides the list of species names.
  It is required for reactive transport.

* `"component molar masses`" [Array(string)] provides the list of
  molar masses of species. It is required for proper conversion to and from
  dimensionless units. Default is 1.

* `"number of liquid components`" [int] is the number of liquid components.

* `"time periods`" [list] contains the list of time periods involved in the simulation.
  The number of time periods is not limited.

  * `"TP #`" [list] defines a particular time period. The numbering
    should be sequential starting with 0.

    * `"PK tree`" [list] describes a hierarchical structure of the process kernels
      that reflect their weak and strong coupling.

      * `"PKNAME`"  [list] name of PK which is used in the
        simulation. Name can be arbitrary but the sublist with the same name
        should exist in the list of PKs (see below).

      * `"PK type`" [string] specifies the type of PK supported by Amanzi. At the moment
        available options are (`"darcy`", `"richards`", `"transport`", `"one-phase energy`",
        `"two-phase energy`", `"reactive transport`", `"flow reactive transport`",
        `"thermal richards`", `"chemistry`", `"transport implicit`", `"transport matrix fracture`",
        `"transport matrix fracture implicit`", `"flow`", and `"darcy matrix fracture`").

      * `"start period time`" [double] is the start time of the current time period.

      * `"end period time`" [double] is the end time of the current time period.

      * `"maximum cycle number`" [int] is the maximum allowed number of cycles in
        the current time period. Special value -1 means unlimited number of cycles.

      * `"initial time step`" is the initial time step for the current time period.

.. code-block:: xml

  <ParameterList>  <!-- parent list -->
  <ParameterList name="cycle driver">
    <Parameter name="component names" type="Array(string)" value="{H+, Na+, NO3-, Zn++}"/>
    <Parameter name="component molar masses" type="Array(double)" value="{1.0e-3, 23.0e-3, 62.0e-3, 65.4e-3}"/>
    <Parameter name="number of liquid components" type="int" value="4"/>
    <ParameterList name="time periods">
      <ParameterList name="TP 0">
        <ParameterList name="PK tree">
          <ParameterList name="_FLOW and REACTIVE TRANSPORT">
            <Parameter name="PK type" type="string" value="flow reactive transport"/>
            <ParameterList name="_REACTIVE TRANSPORT">
              <Parameter name="PK type" type="string" value="reactive transport"/>
              <ParameterList name="_TRANSPORT">
                <Parameter name="PK type" type="string" value="transport"/>
              </ParameterList>
              <ParameterList name="_CHEMISTRY">
                <Parameter name="PK type" type="string" value="chemistry"/>
              </ParameterList>
            </ParameterList>
            <ParameterList name="_FLOW">
              <Parameter name="PK type" type="string" value="darcy"/>
            </ParameterList>
          </ParameterList>
        </ParameterList>
        <Parameter name="start period time" type="double" value="0.0"/>
        <Parameter name="end period time" type="double" value="1.5778463e+09"/>
        <Parameter name="maximum cycle number" type="int" value="-1"/>
        <Parameter name="initial time step" type="double" value="1.57680e+05"/>
      </ParameterList>

      <ParameterList name="TP 1">
      ...
      </ParameterList>
    </ParameterList>
  </ParameterList>
  </ParameterList>

In this simulation, we use the PK labeled as *flow reactive transport*. It is
defined internally as sequential application of two PKs, *flow* and *reactive transport*.
The latter is defined as sequential application of two PKs, *transport* and *chemistry*.
Process kernel *reactive transport* can susbcycle with respect to *flow*.
Process kernel *chemistry* can susbcycle with respect to *transport*.

*/

#ifndef AMANZI_CYCLE_DRIVER_HH_
#define AMANZI_CYCLE_DRIVER_HH_

#include "Epetra_MpiComm.h"
#include "Teuchos_Time.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "FlexibleObservations.hh"
#include "ObservationData.hh"
#include "WalkaboutCheckpoint.hh"
#include "VerboseObject.hh"

namespace Amanzi {

class TimeStepManager;
class Visualization;
class Checkpoint;
class State;
class TreeVector;
class PK;

class CycleDriver {
 public:
  CycleDriver(Teuchos::RCP<Teuchos::ParameterList> glist,
              Teuchos::RCP<Amanzi::State>& S,
              const Comm_ptr_type& comm,
              Amanzi::ObservationData& observations_data);

  // PK methods
  void Setup();
  void Initialize();
  void Init_PK(int);
  void Reset_PK();
  void Finalize();
  void ReportMemory();
  double Advance(double dt);
  void Visualize(bool force = false, const Tag& tag = Tags::DEFAULT);
  void Observations(bool force = false, bool integrate = false);
  void WriteCheckpoint(double dt, bool force = false);
  void WriteWalkabout(bool force);
  //  void RegisterOutput();
  double get_dt(bool after_failuer = false);
  void set_dt(double dt);
  void ResetDriver(int time_period_id);
  // one stop shopping
  Teuchos::RCP<State> Go();

  // access (for unit tests only)
  Teuchos::RCP<const Amanzi::WalkaboutCheckpoint> walkabout() const { return walkabout_; }

 private:
  void CoordinatorInit_();
  void ReadParameterList_();

 private:
  // PK container and factory
  Teuchos::RCP<PK> pk_;

  // states
  Teuchos::RCP<State> S_, S_old_;
  Teuchos::RCP<TreeVector> soln_;

  // time step manager
  Teuchos::RCP<TimeStepManager> tsm_;

  // misc setup information
  Teuchos::RCP<Teuchos::ParameterList> glist_;
  Teuchos::RCP<Teuchos::ParameterList> coordinator_list_;

  std::vector<double> t_, tp_start_, tp_end_, tp_dt_, tp_max_cycle_, tp_max_dt_;
  double max_dt_, min_dt_;
  int cycle0_, cycle1_;
  int num_time_periods_;
  int time_period_id_;

  // Amanzi communicator
  Comm_ptr_type comm_;

  // observations
  Amanzi::ObservationData& observations_data_;
  Teuchos::RCP<FlexibleObservations> observations_;

  // vis and checkpointing
  std::vector<Teuchos::RCP<Visualization>> visualization_;
  std::vector<Teuchos::RCP<Visualization>> failed_visualization_;
  Teuchos::RCP<Checkpoint> checkpoint_;
  bool restart_requested_;
  //  bool output_registered_;
  std::string restart_filename_;

  // time period control
  std::vector<std::pair<double, double>> reset_info_;
  std::vector<std::pair<double, double>> reset_max_;

  // //  checkpoint/restart
  // Teuchos::RCP<Amanzi::Checkpoint> restart_;

  // walkabout
  Teuchos::RCP<Amanzi::WalkaboutCheckpoint> walkabout_;

  // fancy OS
  Teuchos::RCP<VerboseObject> vo_;
};


// non-meber function
// names of all fields that go into a vis file
inline std::set<std::string>
StateVisFields(const State& S)
{
  std::set<std::string> fields;
  for (auto it = S.data_begin(); it != S.data_end(); ++it) {
    for (auto& e : *it->second) {
      if (e.second->io_vis()) fields.insert(it->first);
    }
  }
  return fields;
}

} // namespace Amanzi

#endif
