.. _Amanzi XML Schema :

============================================================
Input File XML Schema 
============================================================

Overview
++++++++

The present document describes how to specify the data required to execute Amanzi and perform a simulation.  This specification should be regarded as a companion to the mathematical requirements document entitled *Mathematical Formulation Requirements and Specifications for the Process Models* (see :ref:`ASCEM Overview <ASCEM Overview>`), and parameterizations of the individual submodels are consistent between Amanzi, the mathematical requirements document and this document.

The open-source, platform independent Akuna_ user environment can generate *Amanzi* models and generate corresponding valid, human-readable XML input files that can then be executed by *Amanzi*.  Example input files are available in the Amanzi source repository.

XML Schema 2.3.0
++++++++++++++++

Amanzi solves a set of parameterized models for multiphase flow in porous media. An Amanzi simulation is specified by providing:

* values for a parameterized PDE-based transport model, including boundary and initial conditions, constitutive laws, and parameterized/phenomenological models for fluid and chemical sources and characterizations of the porous medium,
* parameters controlling the selection of key algorithmic options and output,
* a description of the (discrete) state of the computational system, including a list of the independent variables and instructions for obtaining or generating the discrete mesh, and a characterization of the (parallel) computing environment.

The primary input to *Amanzi* is through an XML file. The Amanzi input XML format is defined in terms of the XML schema that can be found in the Amanzi source code repository.  Users can construct models and generate compliant XML input files using the Akuna_ tool suite.  Users can also choose to generate compliant file using a text editor or other method.

.. In practice, Amanzi is called by a "simulation coordinator" which
.. manages the simulation instructions and orchestrates the flow of
.. data. A basic simulation coordinator is provided with the Amanzi
.. source code distribution. This simple stand-alone coordinator can be
.. used to drive a simple sequence of Amanzi runs, or can serve as a
.. template for user-generated extensions supporting more intricate
.. workflows.

The following is a description of each of the sections with the XML input schema.  Each section includes a short description of what is defined within the section followed by the XML elements available for the section.  Not all elements described are required to generated a valid input file.  The limited set of required elements is noted in each section.  Note that tags for all sections must be present for an input file to be valid even if no elements within that section are required.

Please note, many attributes within the XML list a limited set of specified values.  During validation of the input file or initialization of Amanzi the values in the user provided input file will be compared against the limited set provided in the XML Schema document.  Errors will occur is the values do not match exactly.  These values are CASE SENSITIVE.  The Amanzi schema has been designed will all LOWER CASE values.  Please note this when writing input file.  In particular, `"Exodus II`" will be evaluated as `"exodus ii`".

Amanzi Input
------------

The entire input file is encapsulated under the root tag ``amanzi_input``.  At this level the version of the input schema and the type of simulation to be conducted are specified through required attributes.  The version is specified with the attribute ``version``.  The current schema can be found in the repository and install location amanzi.xsd.  The simulation type is specified using the attribute ``type`` as either *structured* or *unstructured*.  

For an unstructured simulations the root tag would looks as the following.

.. code-block:: xml

    <amanzi_input version="2.3.0" type="unstructured">
    </amanzi_input>

Model Description
-------------------

This section allows the user to provide information about the model being developed and how and when it was developed.  Default units for the model are also stored in this section.  This entire section is optional but encouraged for documentation purposes.

The opening tag ``model_description`` accepts an attribute ``name`` in which the user may give the current model a name.  The available elements within ``model_description`` include: ``comments``, ``author``, ``created``, ``modified``, ``model_name``, ``description``, ``purpose``, and ``units``.  Under the ``units`` element, the user may define the default units to used throughout the rest of the input file.  The options available for the units element are as follows:

+----------------------------+---------------+
| Model_Description Elements | Value Type    |
+============================+===============+
| comments                   | string        |
+----------------------------+---------------+
| author                     | string        |
+----------------------------+---------------+
| created                    | string        |
+----------------------------+---------------+
| modified                   | string        |
+----------------------------+---------------+
| model_id                   | string        |
+----------------------------+---------------+
| description                | string        |
+----------------------------+---------------+
| purpose                    | string        |
+----------------------------+---------------+
| units                      | element block |
+----------------------------+---------------+

The units block specifies the default units to be used through out the input file unless otherwise specified.

+----------------+------------+--------------------+
| Units Elements | Value Type | Value Options      |
+================+============+====================+
| length_unit    | string     | ``m, cm``          |
+----------------+------------+--------------------+
| time_unit      | string     | ``y, d, h, s``     |
+----------------+------------+--------------------+
| mass_unit      | string     | ``kg``             |
+----------------+------------+--------------------+
| conc_unit      | string     | ``molar, mol/m^3`` |
+----------------+------------+--------------------+

Here is an overall example for the model description element.

.. code-block:: xml

  <model_description name="BC Cribs">
    <comments>Some comment here</comments>
    <model_id>DVZ-212a-3217</model_id>
    <author>David Moulton</author>
    <units>
      <length_unit>m</length_unit>
      <time_unit>s</time_unit>
      <mass_unit>kg</mass_unit>
      <conc_unit>molar</conc_unit>
    </units>
  </model_description>


Definitions
-----------

This section allows the user to provide useful definitions to be used throughout the other sections.  Definitions are grouped as element blocks: `constants`_ and `macros`_.

Constants
_________

The user may specify as many constants as desired.  The available constants fall into the following types and descriptions:

+-------------------------+------------------------------------------------------------------+
| Constants Elements      | Description                                                      |
+=========================+==================================================================+
| constant                | general constant definition, can be of any of the following types|
+-------------------------+------------------------------------------------------------------+
| time_constant           | define a constant with a time value                              |
+-------------------------+------------------------------------------------------------------+
| numerical_constant      | define a constant with a numerical value (no units specified)    |
+-------------------------+------------------------------------------------------------------+
| area_mass_flux_constant | define a constant with an area mass flux value                   |
+-------------------------+------------------------------------------------------------------+

Each element has the following format:

+-------------------------+-----------------+----------------+----------------------------------+
| Constants Elements      | Attribute Names | Attribute Type | Attribute Values                 |
+=========================+=================+================+==================================+
| constant                | name            | string         | (user specified name)            |
|                         | value           | string         | (value of constant)              |
|                         | type            | string         | ``none, time, area_mass_flux``   |
+-------------------------+-----------------+----------------+----------------------------------+
| time_constant           | name            | string         | (user specified name)            |
|                         | value           | time(,char)    | (time value with optional units) |
+-------------------------+-----------------+----------------+----------------------------------+
| numerical_constant      | name            | string         | (user specified name)            |
|                         | value           | exponential    | (numerical constant value)       |
+-------------------------+-----------------+----------------+----------------------------------+
| area_mass_flux_constant | name            | string         | (user specified name)            |
|                         | value           | exponential    | (flux value)                     |
+-------------------------+-----------------+----------------+----------------------------------+

Here is an overall structure for the constants element.

.. code-block:: xml

  <constants>
    <constant name="Name of Constant" type="none | time | area_mass_flux" value="constant_value"/>
    <time_constant  name="Name of Time"  value="value,y|d|h|s"/>
    <numerical_constant name="Name of Numerical Constant" value="value_constant"/>
    <area_mass_flux_constant name="Name of Flux Constant" value="value_of_flux"/>
  </constants>

Macros
______

The ``macros`` section defines time and cycle macros.  These specify a series of times or cycles for writing out visualization, checkpoint, walkabout, or observation files.  Each macro type is described in the following table.  The macro can contain a list of specific time values at which to perform an action or a time/cycle interval at which to perform an action.

+--------------------------+-----------------+----------------+-------------------+------------------------------------------------+
| Macros Elements          | Attribute Names | Attribute Type | Sub-Elements      | Sub-Element Type/Value                         |
+==========================+=================+================+===================+================================================+
| time_macro (time series) | name            | string         | time              | time(,unit) / value of time with optional unit |
+--------------------------+-----------------+----------------+-------------------+------------------------------------------------+
| time_macro (interval)    | name            | string         | start             | time(,unit) / value of start time              |
|                          |                 |                | timestep_interval | time(,unit) / time interval between actions    |
|                          |                 |                | stop              | time(,unit) / final time value                 |
|                          |                 |                |                   | ( -1 specifies final time )                    | 
+--------------------------+-----------------+----------------+-------------------+------------------------------------------------+
| cycle_macro (interval)   | name            | string         | start             | integer / cycle number to start action         |
|                          |                 |                | timestep_interval | integer / number of cycles between actions     |
|                          |                 |                | stop              | integer / cycle number to stop action          | 
|                          |                 |                |                   | ( -1 specifies final step )                    | 
+--------------------------+-----------------+----------------+-------------------+------------------------------------------------+


Here are examples of the macros:

.. code-block:: xml

  <time_macro name="Name of Macro">
    <time>Value</time>
  </time_macro>

.. code-block:: xml

  <cycle_macro name="Name of Macro">
    <start>Value</start>
    <timestep_interval>Value</timestep_interval>
    <stop>Value|-1</stop>
  </cycle_macro>


Here is an overall example for the ``definition`` element.

.. code-block:: xml

   <definitions>
     <constants>
       <constant name="zero" type="none" value="0.000"/>
       <constant name="start" type="time" value="1956.0,y"/>
       <constant name="future_recharge" type="area_mass_flux" value="1.48666E-6"/>
       <time_constant name="start_time" value="1956.0,y"/>
       <numerical_constant name="zero" value="0.000"/>
     </constants>
     <macros>
       <time_macro name="Macro 1">
         <time>6.17266656E10</time>
         <time>6.3372710016E10</time>
         <time>6.33834396E10</time>
       </time_macro>
       <cycle_macro name = "Every_1000_timesteps">
         <start>0</start>
         <timestep_interval>1000</timestep_interval>
         <stop>-1 </stop>
       </cycle_macro>
     </macros>
   </definitions>


Execution Control
-----------------

The ``execution_controls`` section defines the general execution of the Amanzi simulation.  Amanzi can execute in four modes: steady state, transient, transient with static flow, or initialize to a steady state and then continue to transient.  The transient with static flow mode does not compute the flow solution at each time step.  During initialization the flow field is set in one of two ways: (1) A constant Darcy velocity is specified in the initial condition; (2) Boundary conditions for the flow (e.g., pressure), along with the initial condition for the pressure field are used to solve for the Darcy velocity. At present this mode only supports the "Single Phase" flow model.

Default values for execution are defined in the ``execution_control_defaults`` element.  These values are used for any time period during the simulation for which the controls were not specified.  Individual time periods of the simulation are defined using ``execution_control`` elements.  For a steady state simulation, only one ``execution_control`` element will be defined.  However, for a transient simulation a series of controls may be defined during which different control values will be used.  For a valid ``execution_controls`` section the ``execution_control_defaults`` element and at least one ``execution_control`` element must appear.

The ``execution_controls`` element has the following subelements:

+-----------------------------+--------------+--------------------------------------------------------------------------+
| Execution_controls Elements | Element Type | Description                                                              |
+=============================+==============+==========================================================================+
| comments                    | string       | user specified comments                                                  |
+-----------------------------+--------------+--------------------------------------------------------------------------+
| verbosity                   | string       | verbosity level ``extreme, high, medium, low, none``                     |
+-----------------------------+--------------+--------------------------------------------------------------------------+
| restart                     | string       | name of Amanzi checkpoint file to be used for restart                    |
+-----------------------------+--------------+--------------------------------------------------------------------------+
| execution_control_defaults  | see below    | default values to be used if not specified in execution_control elements |
+-----------------------------+--------------+--------------------------------------------------------------------------+
| execution_control           | see below    | execution control values for a given time period                         |
+-----------------------------+--------------+--------------------------------------------------------------------------+

Execution_control_defaults
__________________________

The ``execution_control_defaults`` element has the following attributes.

+------------------+----------------+----------------------------------+
| Attribute Names  | Attribute Type | Attribute Values                 |
+==================+================+==================================+
| init_dt          | time           | time value(,unit)                |
+------------------+----------------+----------------------------------+
| max_dt           | time           | time value(,unit)                |
+------------------+----------------+----------------------------------+
| reduction_factor | exponential    | factor for reducing time step    |
+------------------+----------------+----------------------------------+
| increase_factor  | exponential    | factor for increasing time step  |
+------------------+----------------+----------------------------------+
| mode             | string         | ``steady, transient``            |
+------------------+----------------+----------------------------------+
| method           | string         | ``bdf1, picard``                 |
+------------------+----------------+----------------------------------+

Execution_control
_________________

The ``execution_control`` element has the following attributes. 

+------------------+----------------+-----------------------------------------------------+
| Attribute Names  | Attribute Type | Attribute Values                                    |
+==================+================+=====================================================+
| start            | time           | time value(,unit) (start time for this time period) |
+------------------+----------------+-----------------------------------------------------+
| end              | time           | time value(,unit) (stop time for this time period)  |
+------------------+----------------+-----------------------------------------------------+
| max_cycles       | integer        | max cycles to use for structured                    |
+------------------+----------------+-----------------------------------------------------+
| init_dt          | time           | time value(,unit)                                   |
+------------------+----------------+-----------------------------------------------------+
| max_dt           | time           | time value(,unit)                                   |
+------------------+----------------+-----------------------------------------------------+
| reduction_factor | exponential    | factor for reducing time step                       |
+------------------+----------------+-----------------------------------------------------+
| increase_factor  | exponential    | factor for increasing time step                     |
+------------------+----------------+-----------------------------------------------------+
| mode             | string         | ``steady, transient``                               |
+------------------+----------------+-----------------------------------------------------+
| method           | string         | ``bdf1, picard``                                    |
+------------------+----------------+-----------------------------------------------------+

Each ``execution_control`` is required to define a ``start`` time.  The final control period must define an ``end`` time.  It is assumed that the start time of the next control period is the end time of the previous period.  Therefore, it is not required that each ``execution_control`` element have an ``end`` time defined.

The attribute ``max_cycles`` is only valid for transient and transient with static flow execution modes.

The ``execution_control`` section also provides the elements ``comments`` and ``verbosity``.  Users may provide any text within the ``comment`` element to annotate this section.  ``verbosity`` takes the attribute level=`` extreme | high | medium | low | none``.  This triggers increasing levels of reporting from inside Amanzi.  For debugging purposes use the level extreme.

Restarting a simulation is available using the ``restart`` element.  The text given for the ``restart`` element is the name of the Amanzi checkpoint file to be read in and initialized from.

Here is an overall example for the ``execution_control`` element.

.. code-block:: xml

  <execution_controls>
    <execution_control_defaults init_dt="3.168E-08" max_dt="0.01" reduction_factor="0.8" 
                                increase_factor="1.25" mode="transient" method="bdf1"/>
    <execution_control start="0.0 y" end="1956.0 y" init_dt="1 d" max_dt="500.0 y" 
                       reduction_factor="0.8" mode="steady" />
    <execution_control start="1956.0 y" end="3000.0 y" init_dt="1 s" max_dt="10.0 y" 
                       reduction_factor="0.8" mode="transient" />
  </execution_controls>


Numerical Controls
------------------

This section allows the user to define control parameters associated with the underlying numerical implementation.  The list of available options is lengthy.  However, none are required for a valid input file.  The ``numerical_controls`` section is divided up into the subsections: `common_controls`_, `unstructured_controls`_, and `structured_controls`_.  

Common_controls
_______________

The section is currently empty.  However, in future versions controls that are common between the unstructured and structured executions will be moved to this section and given common terminology.

Unstructured_controls
_____________________


The ``unstructured_controls`` sections is divided in the subsections: ``unstr_steady-state_controls``, ``unstr_transient_controls``, ``unstr_linear_solver``, ``unstr_nonlinear_solver``, ``unstr_flow_controls``, ``unstr_transport_controls``, and ``unstr_chemistry_controls``.  The list of available options is as follows:

.. code-block:: xml

  <unstructured_controls>

    <unstr_flow_controls>
      <discretization_method> fv-default | fv-monotone | fv-multi_point_flux_approximation |
                              fv-extended_to_boundary_edges | mfd-default | mfd-optimized_for_sparsity |
                              mfd-support_operator | mfd-optimized_for_monotonicity | 
                              mfd-two_point_flux_approximation </discretization_method>
      <rel_perm_method> upwind-darcy_velocity (default) | upwind-gravity | upwind-amanzi | 
                        other-arithmetic_average | other-harmonic_average </rel_perm_method>
      <preconditioning_strategy> diffusion_operator | linearized_operator (default) </preconditioning_strategy>
      <atmospheric_pressure> exp </atmospheric_pressure>
    </unstr_flow_controls>

    <unstr_transport_controls>
      <algorithm> explicit first-order (default) | explicit second-order | implicit upwind </algorithm> 
      <sub_cycling> on (defulat) | off </sub_cycling> 
      <cfl> exp </cfl>
    </unstr_transport_controls>

    <unstr_chemistry_controls>
      <process_model> implicit operator split | none </process_model>

      <!-- Amanzi native chemistry -->
      <activity_model> unit (default) | debye-huckel </activity_model> 
      <tolerance> exp </tolerance> <!-- default: 100 -->
      <maximum_newton_iterations> int </maximum_newton_iterations> <!-- default: 1e-12 -->
      <auxiliary_data> pH </auxiliary_data> 

      <!-- Pflotran chemistry -->
      <activity_coefficients> timestep (default) | off </activity_coefficients>
      <max_relative_change_tolerance> exp </max_relative_change_tolerance> <!-- suggested 1.0e-16 -->
      <max_residual_tolerance> exp </max_residual_tolerance> <!-- suggested 1.0e-16 -->
      <log_formulation> on (default) | off </log_formulation>
    </unstr_chemistry_controls>

    <unstr_steady-state_controls>
      <min_iterations> int </min_iterations> 
      <max_iterations> int </max_iterations>
      <limit_iterations> int </limit_iterations>
      <nonlinear_tolerance> exp </nonlinear_tolerance> 
      <error_control_options> pressure (default) | residual </error_control_options> 
      <nonlinear_iteration_damping_factor> exp </nonlinear_iteration_damping_factor>
      <max_preconditioner_lag_iterations> int </max_preconditioner_lag_iterations> 
      <max_divergent_iterations> int </max_divergent_iterations> 
      <nonlinear_iteration_divergence_factor> exp </nonlinear_iteration_divergence_factor> 
      <restart_tolerance_relaxation_factor> exp </restart_tolerance_relaxation_factor> 
      <restart_tolerance_relaxation_factor_damping> exp </restart_tolerance_relaxation_factor_damping> 
      <preconditioner> hypre_amg (default) | trilinos_ml | block_ilu </preconditioner> 
      <nonlinear_iteration_initial_guess_extrapolation_order> int </nonlinear_iteration_initial_guess_extrapolation_order> 
      <unstr_initialization>
	<!-- NOTE: including an empty section here turns intialization on with default values
	     To deactive intialization, remove section completely -->
        <clipping_saturation> exp </clipping_saturation> 
        <clipping_pressure> exp </clipping_pressure> 
        <method> picard (default) | darcy_solver </method> 
        <preconditioner> hypre_amg (default) | trilinos_ml | block_ilu </preconditioner> 
        <linear_solver>aztec00 | aztecoo | AztecOO</linear_solver>
        <error_control_options> pressure (default) | residual </error_control_options>
        <convergence_tolerance> exp </convergence_tolerance>
        <max_iterations> int </max_iterations>
      </unstr_initialization>
    </unstr_steady-state_controls>

    <unstr_transient_controls>
      <min_iterations> int </min_iterations>
      <max_iterations> int </max_iterations>
      <limit_iterations> int </limit_iterations>
      <nonlinear_tolerance> exp </nonlinear_tolerance>
      <nonlinear_iteration_damping_factor> exp </nonlinear_iteration_damping_factor>
      <max_preconditioner_lag_iterations> int </max_preconditioner_lag_iterations>
      <max_divergent_iterations> int </max_divergent_iterations>
      <nonlinear_iteration_divergence_factor> exp </nonlinear_iteration_divergence_factor>
      <restart_tolerance_relaxation_factor> exp </restart_tolerance_relaxation_factor> 
      <restart_tolerance_relaxation_factor_damping> exp </restart_tolerance_relaxation_factor_damping>
      <error_control_options> pressure,residual (default) </error_control_options>
      <preconditioner> hypre_amg (default) | trilinos_ml | block_ilu </preconditioner> 
      <initialize_with_darcy> true | false (default) </initialize_with_darcy>
      <nonlinear_iteration_initial_guess_extrapolation_order>int</nonlinear_iteration_initial_guess_extrapolation_order> 
    </unstr_transient_controls>

    <unstr_linear_solver>
      <method> gmres (default) | pcg </method>
      <max_iterations>int </max_iterations> 
      <tolerance> exp </tolerance> 
      <preconditioner> hypre_amg (default) | trilinos_ml | block_ilu </preconditioner> 
    </unstr_linear_solver>

    <unstr_nonlinear_solver name="nka | newton | jfnk | newton_picard" >
      <modify_correction> true | false (default) </modify_correction>
      <update_upwind_frequency> every_timestep (default) | every_nonlinear_iteration </update_upwind_frequency> 
    </unstr_nonlinear_solver>

    <unstr_preconditioners>
      <hypre_amg>
        <hypre_cycle_applications> int </hypre_cycle_applications> <!-- default: 5 suggested range: 1-5 -->
        <hypre_smoother_sweeps> int </hypre_smoother_sweeps> <!-- default: 3 suggested range: 1-5 -->
        <hypre_tolerance> exp </hypre_tolerance> <!-- default: 0.0 suggested range: 0.0-0.1 -->
        <hypre_strong_threshold> exp </hypre_strong_threshold> <!-- default: 0.5 suggested range: 0.2-0.8 -->
      </hypre_amg>
      <trilinos_ml>
        <trilinos_cycle_applications> int </trilinos_cycle_applications> 
        <trilinos_smoother_sweeps> int </trilinos_smoother_sweeps> 
        <trilinos_threshold> exp </trilinos_threshold>  
        <trilinos_smoother_type> jacobi (default) | gauss_seidel | ilu </trilinos_smoother_type> 
      </trilinos_ml>
      <block_ilu>
        <ilu_overlap> int </ilu_overlap> 
        <ilu_relax> exp </ilu_relax> 
        <ilu_rel_threshold> exp </ilu_rel_threshold> 
        <ilu_abs_threshold> exp </ilu_abs_threshold> 
        <ilu_level_of_fill> int </ilu_level_of_fill> 
      </block_ilu>
    </unstr_preconditioners>
  </unstructured_controls>

Here is an overall example for the ``unstructured_controls`` element.

.. code-block:: xml

   <unstructured_controls>
     <comments>Numerical controls comments here</comments>
     <unstr_steady-state_controls>
       <comments>Note that this section contained data on timesteps, which was moved into the execution control section.</comments>
       <min_iterations>10</min_iterations>
       <max_iterations>15</max_iterations>
       <max_preconditioner_lag_iterations>30</max_preconditioner_lag_iterations>
       <nonlinear_tolerance>1.0e-5</nonlinear_tolerance>
     </unstr_steady-state_controls>
     <unstr_transient_controls>
       <comments>Proposed comments section.</comments>
       <bdf1_integration_method min_iterations="10" max_iterations="15" max_preconditioner_lag_iterations="5" />
     </unstr_transient_controls>
     <unstr_linear_solver>
       <comments>Proposed comment section.</comments>
       <method>gmres</method>
       <max_iterations>20</max_iterations>
       <tolerance>1.0e-18</tolerance>
       <preconditioner> trilinos_ml | hypre_amg | block_ilu </preconditioner>
     </unstr_linear_solver>
   </unstructured_controls>


Structured_controls
___________________

The ``structured_controls`` sections is divided in the subsections: ``str_steady-state_controls``, ``str_transient_controls``, ``str_amr_controls``, ``<petsc_options_file>``, and ``max_n_subcycle_transport``.  The list of available options is as follows:

.. code-block:: xml

  <structured_controls>

    <comments>Numerical controls comments here</comments>

    <petsc_options_file> String </petsc_options_file>
    <str_steady-state_controls>
      <max_pseudo_time> Exponential </max_pseudo_time>
      <limit_iterations> Integer </limit_iterations>
      <min_iterations> Integer </min_iterations>
      <min_iterations_2> Integer </min_iterations_2>
      <time_step_increase_factor> Exponential </time_step_increase_factor>
      <time_step_increase_factor_2> Exponential </time_step_increase_factor_2>
      <max_consecutive_failures_1> Integer </max_consecutive_failures_1>
      <time_step_retry_factor_1> Exponential </time_step_retry_factor_1>
      <max_consecutive_failures_2> Integer </max_consecutive_failures_2>
      <time_step_retry_factor_2> Exponential </time_step_retry_factor_2>
      <time_step_retry_factor_f> Exponential </time_step_retry_factor_f>
      <max_num_consecutive_success> Integer </max_num_consecutive_success>
      <extra_time_step_increase_factor> Exponential </extra_time_step_increase_factor>
      <abort_on_psuedo_timestep_failure> true | false </abort_on_psuedo_timestep_failure>
      <limit_function_evals> Integer </limit_function_evals>
      <do_grid_sequence> true | false </do_grid_sequence>
      <grid_sequence_new_level_dt_factor>
        <dt_factor> Exponential </dt_factor> <!-- one element for each AMR level -->
      </grid_sequence_new_level_dt_factor>
    </str_steady-state_controls>


    <str_transient_controls>
      <max_ls_iterations> Integer </max_ls_iterations>
      <ls_reduction_factor> Exponential </ls_reduction_factor>
      <min_ls_factor> Exponential </min_ls_factor>
      <ls_acceptance_factor> Exponential </ls_acceptance_factor>
      <monitor_line_search> Integer </monitor_line_search>
      <monitor_linear_solve> Integer </monitor_linear_solve>
      <perturbation_scale_for_J> Exponential </perturbation_scale_for_J>
      <use_dense_Jacobian> true | false </use_dense_Jacobian>
      <upwind_krel> true | false </upwind_krel>
      <pressure_maxorder> Integer </pressure_maxorder>
      <scale_solution_before_solve> true | false </scale_solution_before_solve>
      <semi_analytic_J> true | false </semi_analytic_J>
      <cfl> Exponential </cfl>
    </str_transient_controls>

    <str_transient_controls>
      <amr_levels> Integer </amr_levels>
      <refinement_ratio>Integer Integer</refinement_ratio> <!-- amr_levels-1 number of integers should be listed-->
      <do_amr_cubcycling> true | false </do_amr_cubcycling>
      <regrid_interval>Integer Integer</regrid_interval> <!-- amr_levels number of integers should be listed-->
      <blocking_factor>Integer Integer</blocking_factor> <!-- amr_levels number of integers should be listed-->
      <number_error_buffer_cells>Integer Integer</number_error_buffer_cells> <!-- amr_levels-1 number of integers should be listed-->
      <max_grid_size>Integer Integer</max_grid_size> <!-- amr_levels number of integers should be listed-->
      <refinement_indicators> 
        <field_name> String </field_name>
        <regions> String </regions>
        <max_refinement_level> Integer </max_refinement_level>
        <start_time> Exponential </start_time>
        <end_time> Exponential </end_time>
        <!-- user may specify exactly 1 of the following -->
        <value_greater> Exponential </value_greater>
        <inside_region> true | false </inside_region>
        <value_less> Exponential </value_less>
        <adjacent_difference_greater> Exponential </adjacent_difference_greater>
        <inside_region> true | false </inside_region>
      </refinement_indicators>
    </str_transient_controls>

    <max_n_subcycle_transport> Integer </max_n_subcycle_transport>
  </structured_controls>

Mesh
----

A mesh must be defined for the simulation to be conducted on.  The mesh can be structured or unstructured.  Structured meshes are always internally generated while unstructured meshes may be generated internally or imported from an existing `Exodus II <http://sourceforge.net/projects/exodusii/>`_ file. Generated meshes in both frameworks are always regular uniformly spaced meshes.

Mesh - Generate (Structured)
____________________________


The ``mesh`` section takes a ``dimension`` element which indicates if the mesh is 2D or 3D. A 2D mesh can be given in 3D space with a third coordinate of 0. If a 2D mesh is specified this impacts other aspects of the input file.  It is up to the user to ensure consistency within the input file.  Other effected parts of the input file include region definitions and initial conditions which use coordinates, the material property permeability which must be specified using the correct subset of x, y, and z coordinates, and the initial condition velocity which also requires the correct subset of x, y, and z coordinates.

This section also takes an element indicating how the mesh is to be internally generated. The ``generate`` element specifies the details about the number of cells in each direction and the low and high coordinates of the bounding box.  It should be noted that in order to accommodate mesh refinement, the number of cells in each direction must be even.

Finally, as in other sections, a ``comments`` element is provide to include any comments or documentation the user wishes.

Here is an example specification for internally generated ``mesh`` element for structured.

.. code-block:: xml

  <mesh> 
    <comments>3D block</comments>
    <dimension>3</dimension>
    <generate>
      <number_of_cells nx="400"  ny="200"  nz="10"/>
      <box  low_coordinates="0.0,0.0,0.0" high_coordinates="200.0,200.0,1.0"/>
    </generate>
  </mesh>

Mesh - Generate (Unstructured)
______________________________

The unstructured portion of Amanzi can utilize different mesh frameworks.  Therefore the framework is specified as an attribute to the ``mesh`` element.  
The available options are: ``mstk``, ``moab``, and ``simple``.  
If no framework is specified, the default ``mstk`` is used.

The ``mesh`` section takes a ``dimension`` element which indicates if the mesh is 2D or 3D. A 2D mesh can be given in 3D space with a third coordinate of 0. If a 2D mesh is specified this impacts other aspects of the input file.  It is up to the user to ensure consistency within the input file.  Other effected parts of the input file include region definitions and initial conditions which use coordinates, the material property permeability which must be specified using the correct subset of x, y, and z coordinates, and the initial condition velocity which also requires the correct subset of x, y, and z coordinates.

This section also takes an element indicating how the mesh is to be internally generated. The ``generate`` element specifies the details about the number of cells in each direction and the low and high coordinates of the bounding box.  

Finally, as in other sections, a ``comments`` element is provide to include any comments or documentation the user wishes.

The following is an example specification for a generated unstructured mesh.

.. code-block:: xml

  <mesh framework="mstk"> 
    <comments>Pseudo 2D</comments>
    <dimension>3</dimension>
    <generate>
      <number_of_cells nx="432"  ny="1"  nz="256"/>
      <box low_coordinates="0.0,0.0,0.0" high_coordinates="216.0,1.0,107.52"/>
    </generate>
  </mesh>

Mesh - Read (Unstructured)
__________________________

The unstructured mode of Amanzi can utilize different mesh frameworks.  Therefore the framework is specified as an attribute to the ``mesh`` element.
The available options are: ``mstk``, ``moab``, and ``simple``.  
If no framework is specified, the default ``mstk`` is used.

The ``mesh`` section takes a ``dimension`` element which indicates if the mesh is 2D or 3D. A 2D mesh can be given in 3D space with a third coordinate of 0. If a 2D mesh is specified this impacts other aspects of the input file.  It is up to the user to ensure consistency within the input file.  Other effected parts of the input file include region definitions and initial conditions which use coordinates, the material property permeability which must be specified using the correct subset of x, y, and z coordinates, and the initial condition velocity which also requires the correct subset of x, y, and z coordinates.

The unstructured mode of Amanzi can read meshes in the Exodus II format.  The ``read`` element contains the subelements ``file``, ``format``, and ``verify`` for specifying the mesh file format (currently only Exodus II) and the mesh file name (relative path allowed).  Any regions or attributes specified in the mesh file will also be read.  These names of the regions and attributes can be utilized in appropriate sections of the input file.  The subelement ``verify`` turns on or off checks performed on the mesh when it is read it.  The mesh verification takes time and is only recommended for debugging meshes on first use.

Finally, as in other sections, a ``comments`` element is provide to include any comments or documentation the user wishes.

Finally, an example of reading an unstructured mesh from a file is given below.

.. code-block:: xml

  <mesh framework="mstk"> 
    <comments>Read from Exodus II</comments>
    <dimension>3</dimension>
    <read>
      <file>dvz.exo</file>
      <format>exodus ii</format>
      <verify>true</verify>
    </read>
  </mesh>

Regions
-------

Regions are geometrical constructs used in Amanzi to define subsets of the computational domain in order to specify the problem to be solved, and the output desired. Regions are commonly used to specify material properties, boundary conditions and observation domains. Regions may represent zero-, one-, two- or three-dimensional subsets of physical space. For a three-dimensional problem, the simulation domain will be a three-dimensional region bounded by a set of two-dimensional regions. If the simulation domain is N-dimensional, the boundary conditions must be specified over a set of regions are (N-1)-dimensional.

Amanzi automatically defines the special region labeled "All", which is the entire simulation domain. Under the "Structured" option, Amanzi also automatically defines regions for the coordinate-aligned planes that bound the domain, using the following labels: "XLOBC", "XHIBC", "YLOBC", "YHIBC", "ZLOBC", "ZHIBC"

The ``regions`` block is required.  Within the region block no regions are required to be defined.  The optional elements valid for both structured and unstructured include ``region``, ``box``, ``point``, ``plane``, and ``logical``.  As in other sections there is also an options ``comments`` element.

The elements ``box``, ``point``, and ``plane`` allow for in-line description of regions.  The ``region`` element uses a subelement to either define a ``box`` or ``plane`` region or specify a region file.  Below are further descriptions of these elements.

Additional regions valid only for unstructured are ``polygonal_surface``.  Additional regions valid only for structured include ``polygon`` and ``ellipse`` in 2D and ``rotated_polygon`` and ``swept_polygon`` in 3D.

Each region definition requires a ``name`` attribute.  
These names must be unique to avoid confusion when other sections refer to the regions.

Box
___

A box region region is defined by a low corner coordinates and high corner coordinates. Box regions can be degenerate in one or more directions.

.. code-block:: xml

  <box name="MyBox" low_coordinates="x_low,y_low,z_low"
                    high_coordinates="x_high,y_high,z_high"/>


Point
_____

A point region is defined by a point coordinates.

.. code-block:: xml

  <point name="Well" coordinate="x,y,z" />

Plane
_____

A plane region is defined by a point on the plane and the normal direction of the plane.

.. code-block:: xml

  <plane name="plane name" location="x,y,z" normal="nx,ny,nz" tolerance="optional exp"/>

The attribute ``tolerance`` is optional.
This value prescribes an absolute tolerance for determining the cell face centroids that lie on the defined plane.

Labeled Set
___________

A labeled set region is a predefined set of mesh entities defined in the Exodus II mesh file. This type of region is useful when applying boundary conditions on an irregular surface that has been tagged in the external mesh generator.  Please note that both the format and entity attribute values are case sensitive. Also not that the attribute ``label`` refers to the name of the region used in the mesh file.  Currently the label/name needs to be an integer value.  Also the region names in the mesh file should be unique to avoid errors and confusion as to which region is being referred to.

.. code-block:: xml

  <region name="region name">
    <region_file label="integer label" name="filename" type="labeled set"
                 format="exodus ii" entity=["cell"|"face"] />
  </region>

Color function
______________

A color function region defines a region based on a specified integer color in a structured color function file. The color values may be specified at the nodes or cells of the color function grid. A computational cell is assigned the color of the data grid cell containing its cell centroid or the data grid nearest its cell-centroid. Computational cell sets are then build from all cells with the specified color value. In order to avoid gaps and overlaps in specifying materials, it is strongly recommended that regions be defined using a single color function file.  At this time, Exodus II is the only file format available.   Please note that both the format and entity attribute values are case sensitive.

.. code-block:: xml

  <region name="region name">
    <region_file label="integer label" name="filename" type="color" 
                 format="exodus ii"  entity=["cell"|"face"]/>
  </region>

Logical
_______

Logical regions are compound regions formed from other primitive type regions using boolean operations. Supported operators are union, intersection, subtraction and complement.  This region type is only valid for the unstructured algorithm.

.. code-block:: xml

    <logical name="logical name">
      <operation>union|intersection|subtraction|complement</operation>
      <region_list>region1, region2, region3<region_list/>
    </logical>


Polygonal_Surface (unstructured only)
_____________________________________

A polygonal_surface region is used to define a bounded planar region and is specified by the number of points and a list of points.  The points must be listed in order and this ordering is maintained during input translation.  This region type is only valid for the unstructured algorithm.

.. code-block:: xml

    <polygonal_surface name="polygon name" num_points="3" tolerance="optional exp">
      <point>X1, Y1, Z1</point>
      <point>X2, Y2, Z2</point>
      <point>X3, Y3, Z3</point>
      <point>X4, Y4, Z4</point>
    </polygonal_surface>

The attribute ``tolerance`` is optional.  
This value prescribes an absolute tolerance for determining the cell face centroids that lie on the defined plane.

Polygon (structured 2D only)
____________________________

A polygon region is used to define a bounded planar region and is specified by the number of points and a list of points.  The points must be listed in order and this ordering is maintained during input translation.  This region type is only valid for the structured algorithm in 2D.

.. code-block:: xml

    <polygon name="polygon name" num_points="3">
      <point> (X1, Y1) </point>
      <point> (X2, Y2) </point>
      <point> (X3, Y3) </point>
    </polygon>

Ellipse (structured 2D only)
____________________________

An ellipse region is used to define a bounded planar region and is specified by a center and X and Y radii.  This region type is only valid for the structured algorithm in 2D.

.. code-block:: xml

    <ellipse name="polygon name" num_points="3">
      <center> (X, Y) </center>
      <radius> (radiusX, radiusY) </radius>
    </ellipse>

Rotated Polygon (structured 3D only)
____________________________________

A rotated_polygon region is defined by a list of points defining the polygon, the plane in which the points exist, the axis about which to rotate the polygon, and a reference point for the rotation axis.  The points listed for the polygon must be in order and the ordering will be maintained during input translation. This region type is only valid for the structured algorithm in 3D.

.. code-block:: xml

    <rotated_polygon name="rotated_polygon name">
      <vertex> (X, Y, Z) </vertex>
      <vertex> (X, Y, Z) </vertex>
      <vertex> (X, Y, Z) </vertex>
      <xyz_plane> (XY | YZ | XZ) </xyz_plane>
      <axis> (X | Y | Z) </axis>
      <reference_point> (X, Y) </reference_point>
    </rotated_polygon>

Swept Polygon (structured 3D only)
__________________________________

A swept_polygon region is defined by a list of points defining the polygon, the plane in which the points exist, the extents (min,max) to sweep the polygon normal to the plane.  The points listed for the polygon must be in order and the ordering will be maintained during input translation. This region type is only valid for the structured algorithm in 3D.

.. code-block:: xml

    <swept_polygon name="swept_polygon name">
      <vertex> (X, Y, Z) </vertex>
      <vertex> (X, Y, Z) </vertex>
      <vertex> (X, Y, Z) </vertex>
      <xyz_plane> (XY | YZ | XZ) </xyz_plane>
      <extent_min> exponential </extent_min>
      <extent_max> exponential </extent_max>
    </swept_polygon>

Geochemistry
------------

The ``geochemistry`` section allows for geochemical constraints to be named and defined.  These named constraints are referred to in the ``initial_conditions``, ``boundary_conditions``, and ``sources`` sections.

Note, at this time the geochemistry constraints must be named and defined in the external chemistry engine input file.  This chemistry input file must be provided by the user.  The capability to create the chemistry input file based on information provided within the XML input file is currently under development.

The ``geochemistry`` section has the subelements ``comments``, ``verbosity``, and an unbounded number of ``constraint`` elements.  The ``verbosity`` element specifies the verbosity level to be used by the chemistry engine.  The current options are: silent, terse, verbose, warnings, and errors.  

Currently, each named constraint must correspond to a constraint with matching name in the external chemistry engine input file.  With the XML input file, the ``constraint`` element takes an attribute ``name``.  This name must match the name of a constraint within the chemistry input file.  In addition, the use may include the constraint information as subelements to the ``constraint`` element for information purposes.  The ``constraint`` element accepts three version of a ``primary`` subelement.  The ``primary`` element describes the constraint on a primary specie, mineral, or gas.  

For each type of primary, the ``primary`` element takes the attributes ``name``, ``type``, and ``value``.  The ``name`` should match a primary listed in the ``phases`` section.  The  ``type`` can be: free_ion, pH, total, mineral, gas, total+sorbed, or charge.  Note, for a non-reacting primary (i.e. tracer or solute), only the type total is valid.  If the type selected is mineral, an additional attribute ``mineral`` giving the mineral name is expected.   If the type selected is gas, an additional attribute ``gas`` giving the gas name is expected.  Again, the names given should match names specified in the ``phases`` section.

.. code-block:: xml

    <geochemistry>
      <verbosity>silent | terse | verbose | warnings | errors</verbosity>
      <constraints> <!-- REQUIRED -->
        <constraint name="string"> <!-- REQUIRED -->
          <primary name="primary_name_string"  initial_guess="exp"   type="free_ion | pH | total | total+sorbed | charge"/>
          <primary name="nonreactive_primary"  initial_guess="exp"   type="total"/>
          <primary name="primary_name_string"  initial_guess="exp"   type="mineral" mineral="mineral_name"/>
          <primary name="primary_name_string"  initial_guess="exp"   type="gas" gas = "gas_name"/>
          <mineral name="mineral_name_string"  volume_fraction="exp" surface_area ="exp"/>
        </constraint>
      </constraints>
    </geochemistry>

Materials
---------

The ``materials`` section allows for 1 or more material to be defined.  Each material element requires an attribute ``name`` to distinguish the material definitions. 

The ``material`` in this context is meant to represent the media through which fluid phases are transported. In the literature, this is also referred to as the "soil", "rock", "matrix", etc. Properties of the material must be specified over the entire simulation domain, and is carried out using the Region constructs defined above. For example, a single material may be defined over the "All" region (see above), or a set of materials can be defined over subsets of the domain via user-defined regions. If multiple regions are used for this purpose, they should be disjoint, but should collectively tile the entire domain. Each material requires a label and the following set of physical properties using the supported models described below.

While many material properties are available for the user to define, the minimum requirements for a valid material definition are specifying the ``assigned_regions``, either ``permeability`` or ``hydraulic_conductivity``, and the ``porosity``.  However, if a capillary pressure model or relative permeability model is chosen (other than ``none``), the associated parameters must also be provided.  Likewise, all model specific parameters must be provided for the chosen dispersion tensor model.

Assigned_regions
________________

The ``assigned_regions`` element list the regions to which the following material properties are to be assigned.  If only 1 material exists, the ``All`` region should be used.  If the material properties are to be assigned to multiple regions, provide a comma separated list of the region names. Leading and trailing white space will be trimmed.  Also, spaces within the region names will be preserved.

.. code-block:: xml

   <assigned_regions>Comma seperated list of Regions</assigned_regions>

Mechanical_properties
_____________________

This element collects a series of mechanical properties as subelements.  As mentioned above, the only required subelement is ``porosity``.

For ``dispersion_tensor`` several models are available.  The model is specified using the ``type`` attribute and additional attributes are used to specify the properties of the given model.  The available dispersion models are described in the following table.


+--------------------------+-----------------+-------------------+
| Dispersion Model         | Attribute Names | Attribute Values  |
+==========================+=================+===================+
| uniform_isotropic        | alpha_l         | Exponential value |
|                          | alpha_t         | Exponential value |
+--------------------------+-----------------+-------------------+
| burnett_frind            | alpha_l         | Exponential value |
|                          | alpha_th        | Exponential value |
|                          | alpha_tv        | Exponential value |
+--------------------------+-----------------+-------------------+
| lichtner_kelkar_robinson | alpha_lh        | Exponential value |
|                          | alpha_tv        | Exponential value |
|                          | alpha_th        | Exponential value |
|                          | alpha_tv        | Exponential value |
+--------------------------+-----------------+-------------------+

.. code-block:: xml

   <mechanical_properties>
     <porosity value="Exponential"/>
     <particle_density value="Exponential"/>
     <specific_storage value="Exponential"/>
     <specific_yield value="Exponential"/>
     <dispersion_tensor type="uniform_isotropic" alpha_l="Exponential" alpha_t="Exponential"/>
     <tortuosity value="Exponential"/>
   </mechanical_properties>

Permeability
____________

For each material either the ``permeability`` or the ``hydraulic_conductivity`` must be specified, but not both.  If specifying the ``permeability`` either a single value for each direction can be given for the entire material or the ``permeability`` values can be read from a file by specifying the attribute ``type`` as "file".  The attribute ``attribute`` gives the name of the attribute with an Exodus II mesh file to be read.  Finally the file name is specified using the attribute ``filename``.

.. code-block:: xml

    <permeability x="Exponential" y="Exponential" z="Exponential"/>
    <permeability filename="file name" type="file" attribute="attribute name"/>

Hydraulic_conductivity
______________________

As noted above, either the ``permeability`` or ``hydraulic_conductivity`` must be specified for each material, but not both.  The values for each direction are listed as attributes.  Currently file read has not yet been implemented.

.. code-block:: xml

    <hydraulic_conductivity x="Exponential" y="Exponential" z="Exponential"/>

Cap_pressure
____________

Capillary pressure can be specified using the element ``cap_pressure``.  The attribute ``model`` specifies whether the van Genuchten or Brooks-Corey model is to be used.  The subelement ``parameters`` lists the model specific parameters.  Note, for both models a smoothing interval can be specified but is optional.  Also, if not specifying capillary pressure this element can be skipped or the value of ``model`` set to "none".

.. code-block:: xml

   <cap_pressure model="van_genuchten | brooks_corey | none">
     <!-- for van_Genuchten -->
     <parameters m="Exponential" alpha="Exponential" sr="Exponential" 
                 optional_krel_smoothing_interval="Exponential"/>
     <!-- for Brooks_Corey -->
     <parameters lambda="Exponential" alpha="Exponential" sr="Exponential" 
                 optional_krel_smoothing_interval="Exponential"/>
   </cap_pressure>

Rel_perm
________

Relative permeability can be specified using the element ``rel_perm``.  The attribute ``model`` specifies whether the Mualam or Burdine model is to be used.  The subelement ``exp`` lists the model specific parameters for Burdine.  Also, if not specifying capillary pressure this element can be skipped or the value of ``model`` set to "none".

.. code-block:: xml

   <rel_perm model="mualem | burdine | none">
     <!-- Burdine only -->
     <exp value="Exponential"/>
   </rel_perm>

Sorption_isotherms
__________________

Kd models can be specified for 1 or more primaries using the ``sorption_isotherms`` element. Note for the Kd model to be active the chemistry state under ``process_kernels`` must be "on" and an engine must be specified.  Also, all primaries must be listed for each material.  Three Kd models are available.  All model parameters for the given model must be present.

.. code-block:: xml

   <sorption_isotherms>
     <primary name="Name of Primary" >
       <kd_model model="linear" kd = "Exponential" />
     </primary>
     <primary name="Name of Primary" >
       <kd_model model="langmuir" kd="Exponential" b="Exponential"/>
     </primary>
     <primary name="Name of Primary" >
       <kd_model model="freundlich" kd="Exponential" n="Exponential" />
     </primary>
   </sorption_isotherms>

Minerals
________

Mineral concentrations are specified using the volume fraction and specific surface area attributes ``volume_fraction`` and ``specific_surface_area`` respectively in the ``minerals`` block

.. code-block:: xml

   <minerals>
     <mineral name="Calcite" volume_fraction="0.1" specific_surface_area="1.0"/>
   </minerals>

Ion_exchange
____________

Ion exchange reactions are specified in the ``ion_exchange`` block.  Cations active in the reaction are grouped under the subelement ``cations``.  The attribute ``cec`` specifies the cation exchange capacity for the reaction.  Each cation is listed in a ``cation`` subelement with the attributes ``name`` and ``value`` to specify the cation name and the associated selectivity coefficient.

.. code-block:: xml

   <ion_exchange>
     <cations cec="750.0">
       <cation name="Ca++" value="0.2953"/>
       <cation name="Mg++" value="0.1666"/>
       <cation name="Na+" value="1.0"/>
     </cations>
   </ion_exchange>

Surface_complexation
____________________

Surface complexation reactions are specified in the ``surface_complexation`` block.  Individual reactions are specified using the ``site`` block.  It has the attributes ``density`` and ``name`` to specify the site density and the name of the site.  Note, the site name must match a surface complexation site in the database file without any leading characters, such as `>`.  The subelement ``complexes`` provides a comma separated list of complexes.  Again, the names of the complexes must match names within the datafile without any leading characters.

.. code-block:: xml

   <surface_complexation>
     <site density="1.908e-3" name="FeOH_s">
       <complexes>FeOHZn+_s, FeOH2+_s, FeO-_s</complexes>
     </site>
     <site density="7.6355e-2" name="FeOH_w">
       <complexes>FeOHZn+_w, FeO-_w, FeOH2+_w</complexes>
     </site>
   </surface_complexation>

An example materials element would look like

.. code-block:: xml

  <materials>
    <material name="Facies_1">
      <comments>Material corresponds to region facies1</comments>
      <assigned_regions>Between_Planes_1_and_2</assigned_regions>
      <mechanical_properties>
        <porosity value="0.4082"/>
        <particle_density value="2720.0"/>
      </mechanical_properties>
      <permeability x="1.9976E-12" y="1.9976E-12" z="1.9976E-13"/>
      <cap_pressure model="van_genuchten">
        <parameters m="0.2294" alpha="1.9467E-04" sr="0.0"/>
      </cap_pressure>
      <rel_perm model="mualem"/>
    </material>
  </materials>

Process Kernels
---------------

Amanzi current employees three process kernels that need to be defined in the input file (1) flow, (2) transport, and (3) chemistry.  The ``process_kernels`` section allows the user to define which kernels are to be used during the section and select high level features of those kernels.  

Flow
____

Currently three scenarios are available for calculated the flow field.  `"richards`" is a single phase, variably saturated flow assuming constant gas pressure.  `"saturated`" is a single phase, fully saturated flow.  `"constant`" is equivalent to the flow model of single phase (saturated) with the time integration mode of transient with static flow in the version 1.2.1 input specification.  This flow model indicates that the flow field is static so no flow solver is called during time stepping. During initialization the flow field is set in one of two ways: (1) A constant Darcy velocity is specified in the initial condition; (2) Boundary conditions for the flow (e.g., pressure), along with the initial condition for the pressure field are used to solve for the Darcy velocity.

.. code-block:: xml

    <flow state = "on | off" 
          model = "richards | saturated | constant" />

Transport
_________

For `"transport`" a `"state`" must be specified.  

.. code-block:: xml

    <transport state = "on | off" />

Chemistry
_________

For `"chemistry`" a combination of `"state`", `"engine`", `"input_filename`", and `"database`" must be specified.  If `"state`" is `"off`" then `"engine`" is set to `"none`".  Otherwise the `"engine`" model must be specified.  If PFloTran is specified as the chemistry engine, the user must specify a PFloTran database filename in the `"database`" attribute.  Also, if PFloTran is the chemistry engine, the user may provide a PFloTran input file using the attribute ``input_filename`` or omit the attribute.  If the attribute is omitted, Amanzi will automatically generate the PFloTran input file using information in the XML input file.  It should be noted the automatically generated file is written to a file using the same name as the XML input file, but with the extension .in.  If such a file exists in the run directory, Amanzi will not overwrite the file and just use the existing file.

.. code-block:: xml

    <chemistry state = "on | off" 
               engine = "amanzi | pflotran | none" 
               input_filename = "string"
               database = "string" />

An example ``process_kernels`` is as follows:

.. code-block:: xml

  <process_kernels>
    <comments>This is a proposed comment field for process_kernels</comments>
    <flow state = "on" model = "richards"/>
    <transport state = "on" algorithm = "explicit first-order" sub_cycling = "on"/>
    <chemistry state = "off" engine="none"/>
  </process_kernels>

Phases
------

The ``phases`` section is used to specify components of each of the phases that are mobile, and species that are contained within them. For each phase, the list identifies the set of all independent variables that are to be stored on each discrete mesh cell.

The terminology for flow in porous media can be somewhat ambiguous between the multiphase and groundwater communities, particularly in regards to "components", "solutes", "primaries/secondaries", and "chemicals". Since Amanzi is designed to handle a wide variety of problems, we must settle on a nomenclature for our use here. In the general problem, multiple "phases" may coexist in the domain (e.g. gaseous, aqueous/liquid, etc), and each is comprised of a number of "components" (section 2.2). In turn, each component may carry a number of "primaries" and "secondaries" and some of these may participate in chemical reactions. As a result of reactions, a chemical source or sink term may appear for the primaries involved in the reaction, including primaries in other mobile phases or in the material matrix. Additionally, certain reactions such as precipitation may affect the flow properties of the material itself during the simulation, and some might affect the properties of the fluid (e.g. brines affect the liquid density). While Amanzi does not currently support chemical reactions and thermal processes, the specification here allows for the existence of the necessary data structures and input data framework. Note that if primary concentrations are significant, the system may be better modeled with that primary treated as a separate component. Clearly, these definitions are highly problem-dependent, so Amanzi provide a generalized interface to accommodate a variety of scenarios.

Currently in Amanzi, primaries are transported in the various phase components, and are treated in "complexes". Each complex is typically in chemical equilibrium with itself and does not undergo phase change. Under these conditions, knowledge of the local concentration of the "basis" or "primary" species (the terms are used here interchangeably) in a chemical complex is sufficient to determine the concentrations of all related secondary species in the phase. Each basis species has a total component concentration and a free ion concentration. The total component concentration for each basis species is a sum of the free ion concentrations in the phase components and its stoichiometric contribution to all secondary species. Amanzi splits the total component concentration into a set of totals for each of the transported phase components, and a total sorbed concentration. Given the free ion concentration of each basis species (and if there is more than one phase, a specification of the thermodynamic relationships that determine the partitioning between phase components (if mass transfer is allowed - not in current Amanzi), we can reconstruct the concentration of the primary and secondary species in each phase. As a result only the basis species are maintained in the state data structures for each phases component.

In addition to primaries in the transported phases, there may be various immobile chemical constituents within the porous media (material) matrix, such as "minerals" and "surface complexes". Bookkeeping for these constituents is managed in Amanzi data structures by generalizing the "primary" concept - a slot in the state is allocated for each of these immobile species, but their concentrations are not included in the transport/flow components of the numerical integration. To allow selective transport of the various primaries, Amanzi uses the concept of primary groups. The aqueous primary concentrations are typically treated together as a group, for example, and often represent the only chemical constituents that are mobile. Thus, the current Amanzi will assume that any other groups specified in an Aqueous phase are immobile.

This section specifies the phases present and specific properties about those phases.  The first grouping is by `liquid_phase`_, `solid_phase`_, and `gas_phase`_.  

Liquid_phase
____________

The ``liquid_phase`` element is required to produce a valid input file.  If primaries are being transported properties of the primaries in the current phase can be specified.  If primaries are being reacted a list of secondary chemical is also specified.

.. code-block:: xml

   <liquid_phase name = "water">
     <viscosity> Exponential </viscosity>
     <density> Exponential </density>
     <dissolved_components> 
       <primaries>
         <primary coefficient_of_diffusion="Exponential" first_order_decay_constant="Exponential"> PrimaryName </primary>
       </primaries> 
       <secondaries>
         <secondary>SecondaryName</secondary>
       </secondaries>
     </dissolved_components>
   </liquid_phase>

Solid_phase
___________

The ``solid_phase`` element allows the user to define a ``minerals`` element under which a series of ``mineral`` elements can be listed to specify any minerals present in the solid phase.  The ``mineral`` elements contain the name of the mineral.

.. code-block:: xml

   <solid_phase>
     <minerals>
       <mineral> MineralName </mineral>
     </minerals> 
   </solid_phase>

Gas_phase
_________

The ``gas_phase`` element allows the user to define a ``gases`` element under which a series of ``gas`` elements can be listed to specify any gases present in the gas phase.  The ``gas`` elements contain the name of the gas.

.. code-block:: xml

   <gas_phase>
     <gases>
       <gas> GasName </gas>
     </gases> 
   </gas_phase>

An example ``phases`` element looks like the following.

.. code-block:: xml

  <phases>
    <liquid_phase name = "water">
      <viscosity>1.002E-03</viscosity>
      <density>998.2</density>
      <dissolved_components> 
        <primaries>
          <primary coefficient_of_diffusion="1.0e-9" first_order_decay_constant="1.0e-9">Tc-99</primary>
        </primaries>
      </dissolved_components>
    </liquid_phase>
  </phases>


Initial Conditions
------------------

The `"initial_conditions`" section contains at least 1 and up to an unbounded number of `"initial_condition`" elements.  Each `"initial_condition`" element defines a single initial condition that is applied to one or more region specified in the ``assigned_regions`` element.  The initial condition can be applied to a liquid phase or solid phase using the appropriate subelement.

To specify a liquid phase the ``liquid_phase`` element is used.  At least one ``liquid_component`` must be specified.  In addition a ``geochemistry_component`` element can be specified.  Under the ``liquid_component`` element an initial condition can be defined.  Under the ``geochemistry_component`` element a geochemistry constraint is listed.

The initial conditions are defined using a specific elements.  The element name indicates the type of condition and the attributes define the necessary information.  Below is a table of the conditions available for the liquid phase and the attributes required to define them.

+-----------------------+------------------+-------------------------------+
| Initial Condition Type| Attributes       | Value Type                    |
+=======================+==================+===============================+
| uniform_pressure      | name             | string                        |
|                       | value            | double/time_constant/constant |
+-----------------------+------------------+-------------------------------+
| linear_pressure       | name             | string                        |
|                       | value            | double/time_constant/constant |
|                       | reference_coord  | coordinate                    |
|                       | gradient         | coordinate                    |
+-----------------------+------------------+-------------------------------+
| velocity              | name             | string                        |
|                       | x                | double/constant               |
|                       | y                | double/constant               |
|                       | (z)              | double/constant               |
+-----------------------+------------------+-------------------------------+
| uniform_saturation    | name             | string                        |
|                       | value            | double/time_constant/constant |
+-----------------------+------------------+-------------------------------+
| linear_saturation     | name             | string                        |
|                       | value            | double/time_constant/constant |
|                       | reference_coord  | coordinate                    |
|                       | gradient         | coordinate                    |
+-----------------------+------------------+-------------------------------+


For the geochemistry_component, an unbounded number of ``constraint`` elements can be listed.  Each ``constraint`` element has the attribute *name*.  The *name*  attribute refers to a geochemical constraint defined in the `"geochemisty`" section.

If in the ``process_kernels`` section the flow model is set to *constant* then the flow field is set in one of the following ways: (1) A constant Darcy velocity is specified in the initial condition (as above); (2) Boundary conditions for the flow (e.g., pressure), along with the initial condition for the pressure field are used to solve for the Darcy velocity.

An example ``initial_conditions`` element looks like the following.

.. code-block:: xml

   <initial_conditions>
     <initial_condition name="Pressure and concentration for entire domain">
       <comments>Initial Conditions Comments</comments>
       <assigned_regions>All</assigned_regions>
       <liquid_phase name = "water">
         <liquid_component name = "water">
           <linear_pressure value = "101325" reference_coord ="(0.0,0.0,0.5)" gradient="(0.0,0.0,-9793.5192)"/>
         </liquid_component>
         <geochemistry_component>
           <constraint name = "initial"/>
         </geochemistry_component>
       </liquid_phase>
     </initial_condition>
   </initial_conditions>

Boundary Conditions
-------------------

Boundary conditions are defined in a similar manor to the initial conditions.  Under the tag ``boundary_conditions`` and series of individual ``boundary_condition`` elements can be defined.  Within each ``boundary_condition`` element the ``assigned_regions`` and ``liquid_phase`` elements must appear.  The boundary condition can be applied to one or more region using a comma separated list of region names.  Under the ``liquid_phase`` element the ``liquid_component`` element must be define.  A ``geochemistry_component`` element may optionally be defined.

Under the ``liquid_component`` and ``geochemistry_component`` elements a time series of boundary conditions is defined using the boundary condition elements available in the table below.  Each component element can only contain one type of boundary condition.  Both elements also accept a *name* attribute to indicate the phase associated with the boundary condition.

+---------------------------+--------------------------------+----------------------------------------+
|Boundary Condition Type    | Attributes                     | Value Type                             |
+===========================+================================+========================================+
| | inward_mass_flux        | | name                         | | string                               | 
| | inward_volumetric_flux  | | start                        | | double/time_constant/constant        |
| | outward_mass_flux       | | value                        | | double                               |
| | outward_volumetric_flux | | function                     | | ``linear | uniform | constant``      |
+---------------------------+--------------------------------+----------------------------------------+
|uniform_pressure           | | name                         | | string                               |
|                           | | start                        | | double/time_constant/constant        |
|                           | | value                        | | double                               |
|                           | | function                     | | ``uniform | constant``               |
+---------------------------+--------------------------------+----------------------------------------+
|linear_pressure            | | name                         | | string                               |
|                           | | gradient_value               | | coordinate                           |
|                           | | reference_point              | | coordinate                           |
|                           | | reference_value              | | double                               |
+---------------------------+--------------------------------+----------------------------------------+ 
|hydrostatic                | | name                         | | string                               |
|                           | | start                        | | double/time_constant/constant        |
|                           | | value                        | | double                               |
|                           | | function                     | | ``uniform | constant``               |
|                           | | coordinate_system            | | ``absolute | relative to mesh top``  |
|                           | | submodel                     | | ``no_flow_above_water_table | none`` |
+---------------------------+--------------------------------+----------------------------------------+ 
|linear_hydrostatic         | | name                         | | string                               |
|                           | | gradient_value               | | coordinate                           |
|                           | | reference_point              | | coordinate                           |
|                           | | reference_water_table_height | | double                               |
|                           | | submodel                     | | ``no_flow_above_water_table | none`` |
+---------------------------+--------------------------------+----------------------------------------+ 
|seepage_face               | | name                         | | string                               |
|(unstructured only)        | | start                        | | double/time_constant/constant        |
|                           | | inward_mass_flux             | | double/time_constant/constant        |
|                           | | function                     | | ``linear | uniform | constant``      |
+---------------------------+--------------------------------+----------------------------------------+
|no_flow                    | | name                         | | string                               |
|                           | | start                        | | double/time_constant/constant        |
|                           | | function                     | | ``linear | uniform | constant``      |
+---------------------------+--------------------------------+----------------------------------------+

For the geochemistry_component, an unbounded number of ``constraint`` elements may be listed.  Each constraint has the attributes *name*, *function*, and *start*.  The function option available is *constant*.  The name should match a constraint defined in the `"phases`" section

An example ``boundary_conditions`` element looks like the following.

.. code-block:: xml

  <boundary_conditions>
    <boundary_condition name = "Recharge at top of the domain">
      <assigned_regions>Recharge_Boundary_WestofCribs,Recharge_Boundary_btwnCribs,Recharge_Boundary_EastofCribs</assigned_regions>
      <liquid_phase name = "water">
        <liquid_component name = "water">
          <inward_mass_flux start="0.0"    function= "constant"  value="pre_1956_recharge"/>
          <inward_mass_flux start="1956.0,y" function= "constant"  value="post_1956_recharge"/>
          <inward_mass_flux start="2012.0,y" function= "constant"  value="future_recharge"/>
          <inward_mass_flux start="3000.0,y" function= "constant"  value="future_recharge"/>
        </liquid_component>
        <geochemistry_component>
          <constraint name = "west"  start="0.0"      function= "constant"/>
          <constraint name = "west2" start="1956.0,y" function= "constant"/>
        </geochemistry_component>
      </liquid_phase>
    </boundary_condition>
  </boundary_conditions>


Sources
-------

Sources are defined in a similar manner to the boundary conditions.  Under the tag ``sources`` and series of individual ``source`` elements can be defined.  Within each ``source`` element the ``assigned_regions`` and ``liquid_phase`` elements must appear.  Sources can be applied to one or more region using a comma separated list of region names.  Under the ``liquid_phase`` element the ``liquid_component`` element must be define.  An unbounded number of ``solute_component`` elements and one ``geochemistry_component`` element may optionally be defined.

Under the ``liquid_component`` and ``solute_component`` elements a time series of boundary conditions is defined using the boundary condition elements available in the table below.  Each component element can only contain one type of source.  Both elements also accept a *name* attribute to indicate the phase associated with the source.

+-------------------------+--------------------+-----------------------------------+
|Liquid Phase Source Type | Attributes         | Value Type                        |
+=========================+====================+===================================+
|volume_weighted          | | start            | | double/time_constant/constant   |
|perm_weighted            | | value            | | double                          |
|                         | | function         | | ``linear | uniform | constant`` |
+-------------------------+--------------------+-----------------------------------+

For the solute component, the source available is ``aqueous_conc`` which has the attributes *name*, *value*, *function*, and *start*.  The function options available are *uniform*, *linear*, and *constant*.


An example ``sources`` element looks like the following.

.. code-block:: xml

  <sources>
    <source name = "Pumping Well" >
      <assigned_regions>Well</assigned_regions>
      <liquid_phase name = "water">
        <liquid_component name="water">
          <volume_weighted start="0.0" function="constant" value="-4.0"/>
        </liquid_component>
      </liquid_phase>
    </source>
  </sources>


Outputs
-------

Output data from Amanzi is currently organized into four specific elements: `"Vis`", `"Checkpoint`", `"Observations`", and `"Walkabout Data`".  Each of these is controlled in different ways, reflecting their intended use.

* `"Visualization`" is intended to represent snapshots of the solution at defined instances during the simulation to be visualized.  The ``vis`` element defines the naming and frequency of saving the visualization files.  The visualization files may include only a fraction of the state data, and may contain auxiliary "derived" information.

* `"Checkpoint`" is intended to represent all that is necessary to repeat or continue an Amanzi run.  The specific data contained in a checkpoint dump is specific to the algorithm options and mesh framework selected.  Checkpoint is special in that no interpolation is performed prior to writing the data files; the raw binary state is necessary.  As a result, the user is allowed to only write checkpoint at the discrete intervals of the simulation. The ``checkpoint`` element defines the naming and frequency of saving the checkpoint files.

* `"Observations`" is intended to represent diagnostic values to be returned to the calling routine from Amanzi's simulation driver.  Observations are typically generated at arbitrary times, and frequently involve various point samplings and volumetric reductions that are interpolated in time to the desired instant.  Observations may involve derived quantities or state fields.  The ``observations`` element may define one or more specific observation.

* `"Walkabout Data`" is intended to be used as input to the particle tracking software Walkabout.

NOTE: Each output type allows the user to specify the base_filename or filename for the output to be written to.  The string format of the element allows the user to specify the relative path of the file.  It should be noted that the Amanzi I/O library does not create any new directories.  Therefore, if a relative path to a location other than the current directory is specified Amanzi assumes the user (or the Agni controller) has already created any new directories.  If the relative path does not exist the user will see error meesages from the HDF5 library indicating failure to create and open the output file.

Vis
___

The ``vis`` element defines the visualization file naming scheme and how often to write out the files.  The ``base_filename`` element contain the text component of the how the visualization files will be named.  The ``base_filename`` is appended with an index number to indicate the sequential order of the visualization files.  The ``num_digits`` elements indicates how many digits to use for the index.  See the about NOTE about specifying a file location other than the current working directory. Finally, the ``time_macros`` or ``cycle_macros`` element indicates previously defined time_macros or cycle_macros to be used to determine the frequency at which to write the visualization files.  One or more macro can be listed in a comma separated list.  Amanzi will converted the list of macros to a single list of times or cycles contained by all of the macros listed and output accordingly.

The ``vis`` element also includes an optional subelement ``write_regions``.  This was primarily implemented for debugging purposes but is also useful for visualizing fields only on specific regions.  The subelement accepts an arbitrary number of subelements named ``field``, with attibutes ``name`` (a string) and ``regions`` (a comma separated list of region names).  For each such subelement, a field will be created in the vis files using the name as a label.  The field will be initialized to 0, and then, for region list R1, R2, R3..., cells in R1 will be set to 1, cells in R2 will be set to 2, etc.  When regions in the list overlap, later ones in the list will take precedence. 

An example ``vis`` element looks like the following.

.. code-block:: xml

   <vis>
     <base_filename>plot</base_filename>
     <num_digits>5</num_digits>
     <time_macros>Macro 1</time_macros>
     <write_regions>
       <field name="fieldname" regions="region1, region2, region3" />
     </write_regions>
   </vis>

Checkpoint
__________

The ``checkpoint`` element defines the file naming scheme and frequency for writing out the checkpoint files.  The ``base_filename`` element contains the text component of the how the checkpoint files will be named.  The ``base_filename`` is appended with an index number to indicate the sequential order of the checkpoint files.  The ``num_digits`` elements indicates how many digits to use for the index.  See the about NOTE about specifying a file location other than the current working directory.  Finally, the ``cycle_macro`` element indicates the previously defined cycle_macro to be used to determine the frequency at which to write the checkpoint files.

An example ``checkpoint`` element looks like the following.

.. code-block:: xml

   <checkpoint>
     <base_filename>chk</base_filename>
     <num_digits>5</num_digits>
     <cycle_macro>Every_1000_steps</cycle_macro>
   </checkpoint>

Observations
____________

The ``observations`` element defines the file for writing observations to and specifies individual observations to be made.  At this time, all observations are written to a single file defined in the ``filename`` element.  See the about NOTE about specifying a file location other than the current working directory. Also, observations are only available for the liquid phases.  Therefore individual observations are defined in subelements under the ``liquid_phase`` tag.  The ``liquid_phase`` tag takes an attribute ``name`` to identify which phase the observations are associated with.

The element name of individual observations indicate the quantity being observed.  Below is a list of currently available observations.  Individual observations require the subelements ``assigned_regions``, ``functional``, and ``time_macros``.  ``aqueous_conc`` and ``primary_volumetric_flow_rate`` observations also require the name of the primary.  This is specified with an extra subelement ``primary``. 

Available Observations:

- integrated_mass
- volumetric_water_content
- gravimetric_water_content
- aqueous_pressure
- x_aqueous_volumetric_flux
- y_aqueous_volumetric_flux
- z_aqueous_volumetric_flux
- material_id
- hydraulic_head
- aqueous_mass_flow_rate
- aqueous_volumetric_flow_rate
- aqueous_conc
- drawdown
- primary_volumetric_flow_rate

An example ``observations`` element looks like the following.

.. code-block:: xml

    <observations>
      <filename>observation.out</filename>
      <liquid_phase name="water">
        <aqueous_pressure>
          <assigned_regions>Obs_r1</assigned_regions>
          <functional>point</functional>
          <time_macros>Observation Times</time_macros>
        </aqueous_pressure>
        <aqueous_pressure>
          <assigned_regions>Obs_r2</assigned_regions>
          <functional>point</functional>
          <time_macros>Observation Times</time_macros>
        </aqueous_pressure>
        <aqueous_pressure>
          <assigned_regions>Obs_r2</assigned_regions>
          <functional>point</functional>
          <time_macros>Observation Times</time_macros>
        </aqueous_pressure>
      </liquid_phase>
    </observations>

.. _Akuna : http://esd.lbl.gov/research/projects/ascem/thrusts/platform/
.. _Mathematical Formulation Requirements and Specifications for the Process Models: http://software.lanl.gov/ascem/trac/attachment/wiki/Documents/ASCEM-HPC-ProcessModels_2011-01-0a.pdf

Walkabout
_________

The ''walkabout'' element defines the file naming scheme and frequency for writing out the walkabout files.  As mentioned above, the user does not influence what is written to the walkabout files only the writing frequency and naming scheme.  Thus, the ''walkabout'' element requires the subelements ``base_filename``, ``num_digits``, and ``cycle_macro``.

The *base_filename* element contain the text component of the how the walkabout files will be named.  The *base_filename* is appended with an index number to indicate the seqential order of the walkabout files.  The *num_digits* elements indicates how many digits to use for the index.  Final the *cycle_macro* element indicates the previously defined cycle_macro to be used to determine the frequency at which to write the walkabout files.

Example:

.. code-block:: xml

  <walkabout>
    <base_filename>chk</base_filename>
    <num_digits>5</num_digits>
    <cycle_macro>Every_100_steps</cycle_macro>
  </walkabout>

