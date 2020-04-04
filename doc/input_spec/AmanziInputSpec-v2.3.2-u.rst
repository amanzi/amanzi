===========================================================
Unstructured Amanzi XML Input Specification (Version 2.3.2)
===========================================================

.. contents:: **Table of Contents**


Overview
========

The Amanzi simulator evolves a system of conservation equations for reacting flows in porous media, as detailed in the ASCEM 
report entitled `"Amanzi Theory Guide, Mathematical Modeling Requirement`" (hereafter referred to as the 'Amanzi Theory Guide (ATG)'). 
The purpose of the present document is to specify the data required to execute Amanzi.  This specification should be regarded as a companion to the ATG, and parameterizations of the individual submodels are consistent between Amanzi, the ATG and this document. Where applicable, the relevant sections of the ATG are indicated.

All data required to execute Amanzi is specified within an XML formated file laid out according to the Amanzi input schema.
The current version of the Amanzi schema is located with the Amanzi source code repository.
The following discusses each section of the schema, its purpose and provides examples.
Further details can be found in the schema document doc/input_spec/schema/amanzi.xsd.

Please note, many attributes within the XML list a limited set of specified values.  During validation of the input file or initialization of Amanzi the values in the user provided input file will be compared against the limited set provided in the XML Schema document.  Errors will occur is the values do not match exactly.  These values are CASE SENSITIVE.  The Amanzi schema has been designed will all LOWER CASE values.  Please note this when writing input file.  In particular, `"Exodus II`" will be evaluated as `"exodus ii`".

All user-defined names are capitalized to highlight that they are not a part of the input spec.


Amanzi Input
============

Here, the user specifies which version of the input the input file adheres to. The user also specifies the overall type of simulation being run.  Amanzi supports both structured and unstructured numerical solution approaches.  This flexibility has a direct impact on the selection and design of the underlying numerical algorithms, the style of the software implementations, and, ultimately, the complexity of the user-interface. The attribute *type* is used to selected between the following:

* ``Structured``: This instructs Amanzi to use BoxLib data structures and an associated paradigm to numerically represent the flow equations.  Data containers in the BoxLib software library, developed by CCSE at LBNL, are based on a hierarchical set of uniform Cartesian grid patches.  ``Structured`` requires that the simulation domain be a single coordinate-aligned rectangle, and that the "base mesh" consists of a logically rectangular set of uniform hexahedral cells.  This option supports a block-structured approach to dynamic mesh refinement, wherein successively refined subregions of the solution are constructed dynamically to track "interesting" features of the evolving solution.  The numerical solution approach implemented under the ``Structured`` framework is highly optimized to exploit regular data and access patterns on massively parallel computing architectures. 

* ``Unstructured``: This instructs Amanzi to use data structures provided in the Trilinos software framework.  To the extent possible, the discretization algorithms implemented under this option are largely independent of the shape and connectivity of the underlying cells.  As a result, this option supports an arbitrarily complex computational mesh structure that enables users to work with numerical meshes that can be aligned with geometrically complex man-made or geostatigraphical features.  Under this option, the user typically provides a mesh file that was generated with an external software package.  The following mesh file formats are currently supported: `"Exodus II`".  Amanzi also provides a rudimentary capability to generate regular meshes within the unstructured framework internally.

An example root tag of an input file would look like the following.

.. code-block:: xml

  <amanzi_input version="2.2.1" type="unstructured"/>


Model Description
=================

This allows the users to provide a name and general description of model being developed.  This is also the section in which the units for the problem are stored. This entire section is optional but encouraged as documentation.

.. code-block:: xml

  <model_description name="NAME of MODEL" >
      Required Elements: NONE
      Optional Elements: comment, author, created, modified, model_id, description, purpose, units
  </model_description>


Units
-----

The ``units`` element defines the default units to be assumed for the entire input file.
Amanzi's internal default units are SI units.
Conversion from the listed units to Amanzi's internal default units is done during conversion
of this spec to the internal (developers') spec.

``units`` has the optional elements of length, time, mass, and concentration.  Each of those in turn have their own structure.  The structures are as follows.

REMINDER - UNITS ARE NOT IMPLEMENTED YET

.. code-block:: xml

  <units>
      Required Elements: NONE
      Optional Elements: length_unit, time_unit, mass_unit, conc_unit
  </units>

Acceptable values for each unit are as follows:

+------------------+----------------------------+
| Units Elements   | Value Options              |
+==================+============================+
| length_unit      | km, m, yr, ft, in, or cm   |
+------------------+----------------------------+
| time_unit        | y, noleap, d, h, min, or s |
+------------------+----------------------------+
| mass_unit        | ton, kg, g, or lb          |
+------------------+----------------------------+
| volume_unit      | m3, gal, or L              |
+------------------+----------------------------+
| amount_unit      | mol                        |
+------------------+----------------------------+
| conc_unit        | molar, SI, ppm, or ppb     |
+------------------+----------------------------+
| temperature_unit | K, X, or F                 |
+------------------+----------------------------+
| derived units    | Pa, and J                  |
+------------------+----------------------------+

Here is an overall example for the model description element.

.. code-block:: xml

  <model_description name="DVZ 3layer 2D">
    <comments>This is a simplified 3-layer DVZ problem in 2D with two cribs (Flow+Transport)</comments>
    <model_name>DVZ 3layer</model_name>
    <author>d3k870</author>
    <units>
      <length_unit>m</length_unit>
      <time_unit>s</time_unit>
      <mass_unit>kg</mass_unit>
      <conc_unit>molar</conc_unit>
    </units>
  </model_description>


Definitions
===========

Definitions allows the user the define and name constants, times, and macros to be used in later sections of the input file.  This is to streamline the look and readability of the input file.  The user should take care not to reuse names within this section or other sections.  This may have unindented consequences.

.. code-block:: xml

  <definitions>
      Required Elements: NONE
      Optional Elements: constants, macros
  </definitions>

Constants
---------

Here the user can define and name constants to be used in other sections of the input file.  Note that if a name is repeated the last read value will be retained and all others will be overwritten.  See `Constants`_ for specifying time units other than seconds.

.. code-block:: xml

  <constants>
      Required Elements: NONE
      Optional Elements: constant, time_constant, numerical_constant, area_mass_flux_constant 
  </constants>

A ``constant`` has three attributes ``name``, ``type``, and ``value``.  The user can provide any name, but note it should not be repeated anywhere within the input to avoid confusion.  The available types include: `"none`", `"time`", `"numerical`", and `"area_mass_flux`".  Values assigned to constants of type `"time`" can include known units, otherwise seconds will be assumed as the default. See `Constants`_ for specifying time units other than seconds.

.. code-block:: xml

    <constant name="STRING" type="none | time | numerical | area_mass_flux" value="constant_value"/>

A ``time_constant`` is a specific form of a constant assuming the constant type is a time.  It takes the attributes ``name`` and ``value`` where the value is a time (time unit optional).

.. code-block:: xml

    <time_constant name="NAME of TIME" value="time,y|d|h|s"/>

A ``numerical_constant`` is a specific form of a constant.  It takes the attributes ``name`` and ``value``. 

.. code-block:: xml

    <numerical_constant name="NAME of NUMERICAL CONSTANT" value="value_constant"/>

A ``area_mass_flux_constant`` is a specific form of a constant.  It takes the attributes ``name`` and ``value`` where the value is an area mass flux. 

.. code-block:: xml

    <area_mass_flux_constant name="NAME of FLUX CONSTANT" value="value_of_flux"/>

Macros
------

The ``macros`` section defines time, cycle, and variable macros.  These specify a list or interval for triggering an action, particularly, writing out visualization, checkpoint, walkabout, or observation files.  

.. code-block:: xml

  <constants>
      Required Elements: NONE
      Optional Elements: time_macro, cycle_macro
  </constants>

Time_macro
__________

The ``time_macro`` requires an attribute ``name``.  The macro can then either take the form of one or more labeled time subelements or the subelements ``start``, ``timestep_interval``, and ``stop`` again containing labeled times.  A ``stop`` value of -1 will continue the cycle macro until the end of the simulation.  The labeled times can be time values assuming the default time unit of seconds or including a known time unit.

.. code-block:: xml

  <time_macro name="NAME of MACRO">
    <time>value</time>
  </time_macro>

or 

.. code-block:: xml

  <time_macro name="NAME of MACRO">
    <start> time_value </start>
    <timestep_interval> time_interval_value </timestep_interval>
    <stop> time_value | -1 </stop>
  </time_macro>


Cycle_macro
___________

The ``cycle_macro`` requires an attribute ``name`` and the subelements ``start``, ``timestep_interval``, and ``stop`` with integer values.  A ``stop`` value of -1 will continue the cycle macro until the end of the simulation.

.. code-block:: xml

  <cycle_macro name="NAME of MACRO">
    <start> cycle_value </start>
    <timestep_interval>value</timestep_interval>
    <stop>value|-1</stop>
  </cycle_macro>

An example ``definition`` section would look as the following:

.. code-block:: xml

  <definitions>
    <constants>
      <constant name="BEGIN"            type="none"           value="0.000"/>
      <constant name="START"            type="time"           value="1956.0,y"/>
      <constant name="B-18_RELEASE_END" type="time"           value ="1956.3288,y"/>
      <constant name="future_recharge"  type="area_mass_flux" value="1.48666e-6"/>
      <numerical_constant name="ZERO" value="0.000"/>
    </constants>
    <macros>
      <time_macro name="MACRO 1">
        <time>6.17266656E10</time>
        <time>6.172982136E10</time>
        <time>6.173297712E10</time>
        <time>6.3372710016E10</time>
        <time>6.33834396E10</time>
      </time_macro>
      <cycle_macro name="EVERY_1000_TIMESTEPS">
        <start>0</start>
        <timestep_interval>1000</timestep_interval>
        <stop>-1</stop>
      </cycle_macro>
    </macros>
  </definitions>


Execution Controls
==================

The ``execution_controls`` section defines the general execution of the Amanzi simulation.  Amanzi can execute in four modes: steady state, transient, transient with static flow, or initialize to a steady state and then continue to transient.  The transient with static flow mode does not compute the flow solution at each time step.  During initialization the flow field is set in one of two ways: (1) A constant Darcy velocity is specified in the initial condition; (2) Boundary conditions for the flow (e.g., pressure), along with the initial condition for the pressure field are used to solve for the Darcy velocity. At present this mode only supports the "Single Phase" flow model.

.. code-block:: xml
  
  <execution_controls>
      Required Elements: execution_control_defaults, execution_control (1 or more)
      Optional Elements: comments, verbosity, restart | initialize
  </execution_controls>

The ``execution_controls`` block is required.

Verbosity
---------

The ``verbosity`` element specifies the level of output messages provided by Amanzi.  If not present, the default value of `"medium`" will be used.

.. code-block:: xml
  
  <verbosity level="none | low | medium | high | extreme" />
 
A level of `"extreme`" is recommended for developers.  For users trying to debug input files or monitor solver performance and convergence `"high`" is recommended.

Restart, Initialize
-------------------

The ``restart`` and ``initialize`` elements specify the name of an Amanzi checkpoint file used to initialize a run.  Only one of these two may be present.  ``restart`` indicates that the run is to be continued from where it left off.  ``initialize`` indicates that a completely new run is desired, but that the state fields in the named checkpoint file should be used to initialize the state, rather than the initial conditions block in the input.

TODO: DEFINE RESTART VS INITIALIZE HERE

Execution_control_defaults
--------------------------

The ``execution_control_defaults`` element specifies default values to be utilized when not specified in individual ``execution_control`` elements.   For a valid ``execution_controls`` section the ``execution_control_defaults`` element is *required*.  The attributes available are:

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
| method           | string         | ``bdf1``                         |
+------------------+----------------+----------------------------------+
| max_cycles       | integer        | max number of cycles to use      |
+------------------+----------------+----------------------------------+

Execution_control
-----------------

Individual time periods of the simulation are defined using ``execution_control`` elements.  For a steady state simulation, only one ``execution_control`` element will be defined.  However, for a transient simulation a series of controls may be defined during which different control values will be used.  For a valid ``execution_controls`` section at least one ``execution_control`` element is *required*.  Any attributes not specified in the ``execution_control`` element will use the value defined in the above ``execution_control_defaults`` element.  The attributes available are:
  
+------------------+----------------+----------------------------------------------------------+
| Attribute Names  | Attribute Type | Attribute Values                                         |
+==================+================+==========================================================+
| start            | time           | | time value(,unit) (start time for this time period)    |
|                  |                | | (*required* for each ``execution_control`` element)    |
+------------------+----------------+----------------------------------------------------------+
| end              | time           | | time value(,unit) (stop time for this time period)     |
|                  |                | | (only *required* once in ``execution_controls`` block) |
+------------------+----------------+----------------------------------------------------------+
| init_dt          | time           | time value(,unit)                                        |
+------------------+----------------+----------------------------------------------------------+
| max_dt           | time           | time value(,unit)                                        |
+------------------+----------------+----------------------------------------------------------+
| reduction_factor | exponential    | factor for reducing time step                            |
+------------------+----------------+----------------------------------------------------------+
| increase_factor  | exponential    | factor for increasing time step                          |
+------------------+----------------+----------------------------------------------------------+
| mode             | string         | ``steady, transient``                                    |
+------------------+----------------+----------------------------------------------------------+
| method           | string         | ``bdf1``                                                 |
+------------------+----------------+----------------------------------------------------------+
| max_cycles       | integer        | max number of cycles to use                              |
+------------------+----------------+----------------------------------------------------------+

Each ``execution_control`` element *requires* a start time.  If multiple ``execution_control`` elements are defined ``end`` times are not required for each element.  The ``start`` time of the next execution section is used as the ``end`` of the previous section.  However, at least one ``end`` time *must* defined within the ``execution_controls`` block.

Under the structure algorithm, the attribute ``max_cycles`` is only valid for transient and transient with static flow execution modes.

Here is an overall example for the ``execution_control`` element.

.. code-block:: xml

  <execution_controls>
    <execution_control_defaults init_dt="3.168E-08" max_dt="0.01" reduction_factor="0.8" increase_factor="1.25" mode="transient" method="bdf1"/>
    <execution_control start="0.0;y" end="1956.0,y" init_dt="0.01" max_dt="500.0" reduction_factor="0.8" mode="steady"   />
    <execution_control start="B-17_release_begin" />
    <execution_control start="B-17_release_end" />
    <execution_control start="B-18_release_begin" />
    <execution_control start="B-18_release_end" end="3000.0,y" />
  </execution_controls>

Numerical Controls
==================

This section allows the user to define control parameters associated with the underlying numerical implementation.  The list of available options is lengthy.  However, none are required for a valid input file.  The ``numerical_controls`` section is divided up into the subsections: `common_controls`_, `unstructured_controls`_, and `structured_controls`_.  The ``common_controls`` section is currently empty.  However, in future versions controls that are common between the unstructured and structured executions will be moved to this section and given common terminology.

.. code-block:: xml

  <numerical_controls>
      Required Elements: unstructured_controls [U] *OR* structured_controls [S]
      Optional Elements: comments, common_controls [U/S], unstructured_controls [U], structured_controls [S]
  </numerical_controls>

Common_controls
---------------

The section is currently empty.  However, in future versions controls that are common between the unstructured and structured executions will be moved to this section and given common terminology.

Unstructured_controls
---------------------

The ``unstructured_controls`` sections is divided in the subsections specific to the process kernels and the numerical solver mode.  Either the ``unstructured_controls`` or the ``structured_controls`` section headings are required.  However, no options within the sections are required.  The list of available options is as follows:

.. code-block:: xml

  <unstructured_controls>
      Required Elements: NONE
      Optional Elements: unstr_flow_controls, unstr_transport_controls, unstr_chemistry_controls, unstr_steady-state_controls, unstr_transient_controls, unstr_linear_solver, unstr_nonlinear_solver, unstr_preconditioners
  </unstructured_controls>

Unstr_flow_controls
___________________

``unstr_flow_controls`` specifies numerical controls for the flow process kernel available under the unstructured algorithm.  It has the following subelements:

+--------------------------+--------------+-------------------------------------------------------------+
| Element Names            | Content Type | Content Value                                               |
+==========================+==============+=============================================================+
| discretization_method    | string       | | ``fv-default, fv-monotone,``                              |
|                          |              | | ``fv-multi_point_flux_approximation,``                    |
|                          |              | | ``fv-extended_to_boundary_edges,``                        |
|                          |              | | ``mfd-default, mfd-optimized_for_sparsity,``              | 
|                          |              | | ``mfd-support_operator, mfd-optimized_for_monotonicity,`` | 
|                          |              | | ``mfd-two_point_flux_approximation``                      |
+--------------------------+--------------+-------------------------------------------------------------+
| rel_perm_method          | string       | | ``upwind-darcy_velocity, upwind-gravity, upwind-amanzi,`` | 
|                          |              | | ``other-arithmetic_average, other-harmonic_average``      |
|                          |              | | *default = upwind-darcy_velocity*                         |
+--------------------------+--------------+-------------------------------------------------------------+
| preconditioning_strategy | string       | | ``diffusion_operator, linearized_operator``               |
|                          |              | | *default = linearized_operator*                           |
+--------------------------+--------------+-------------------------------------------------------------+
| atmospheric_pressure     | exponential  | value of atmospheric pressure                               |
+--------------------------+--------------+-------------------------------------------------------------+

Unstr_transport_controls
________________________

``unstr_transport_controls`` specifies numerical controls for the transport process kernel available under the unstructured algorithm.  It has the following subelements:

+---------------+--------------+----------------------------------------------------+
| Element Names | Content Type | Content Value                                      |
+===============+==============+====================================================+
| algorithm     | string       | | ``explicit first-order, explicit second-order,`` |
|               |              | | ``implicit upwind, none``                        |
|               |              | | *default = explicit first-order*                 |
+---------------+--------------+----------------------------------------------------+
| sub_cycling   | string       | | ``on, off``                                      | 
|               |              | | *default = on*                                   |
+---------------+--------------+----------------------------------------------------+
| cfl           | exponential  | cfl condition number                               |
+---------------+--------------+----------------------------------------------------+


Unstr_chemistry_controls
________________________

``unstr_chemistry_controls`` specifies numerical controls for the chemistry process kernel available under the unstructured algorithm. Currently two chemistry engines are available through Amanzi.  They are the Amanzi native chemistry engine or the PFloTran chemistry engine available through the Alquimia interface.  Options for both engines are specified here. 

The subelements pertaining to the Amanzi native chemistry engine are:

+----------------------------------------+--------------+-----------------------------------+
| Element Names                          | Content Type | Content Value                     |
+========================================+==============+===================================+
| process_model                          | string       | ``implicit operator split, none`` |
+----------------------------------------+--------------+-----------------------------------+
| activity_model                         | string       | ``unit, debye-huckel``            |
+----------------------------------------+--------------+-----------------------------------+
| maximum_newton_iterations              | integer      |                                   |
+----------------------------------------+--------------+-----------------------------------+
| tolerance                              | exponential  |                                   |
+----------------------------------------+--------------+-----------------------------------+
+----------------------------------------+--------------+-----------------------------------+
| auxiliary_data                         | string       | ``pH``                            |
+----------------------------------------+--------------+-----------------------------------+

The subelements pertaining to the pflotran chemistry engine are:

+----------------------------------------+--------------+-----------------------------------+
| Element Names                          | Content Type | Content Value                     |
+========================================+==============+===================================+
| activity_coefficients                  | string       | ``timestep, off``                 |
+----------------------------------------+--------------+-----------------------------------+
| max_relative_change_tolerance          | exponential  |                                   |
+----------------------------------------+--------------+-----------------------------------+
| max_residual_tolerance                 | exponential  |                                   |
+----------------------------------------+--------------+-----------------------------------+
| min_time_step                          | exponential  |                                   |
+----------------------------------------+--------------+-----------------------------------+
| max_time_step                          | exponential  |                                   |
+----------------------------------------+--------------+-----------------------------------+
| initial_time_step                      | exponential  |                                   |
+----------------------------------------+--------------+-----------------------------------+
| time_step_control_method               | string       | ``fixed, simple``                 |
+----------------------------------------+--------------+-----------------------------------+
| time_step_cut_threshold                | integer      | (use only if method = simple)     |
+----------------------------------------+--------------+-----------------------------------+
| time_step_cut_factor                   | exponential  | (use only if method = simple)     |
+----------------------------------------+--------------+-----------------------------------+
| time_step_increase_threshold           | integer      | (use only if method = simple)     |
+----------------------------------------+--------------+-----------------------------------+
| time_step_increase_factor              | exponential  | (use only if method = simple)     |
+----------------------------------------+--------------+-----------------------------------+
| log_formulation                        | string       | ``on, off``                       |
+----------------------------------------+--------------+-----------------------------------+
| generate_chemistry_engine_inputfile    | string       |                                   |
+----------------------------------------+--------------+-----------------------------------+
| read_chemistry_engine_inputfile        | string       |                                   |
+----------------------------------------+--------------+-----------------------------------+

Unstr_steady-state_controls
___________________________

+---------------------------------------------+---------------+---------------------------------------+
| Element Names                               | Content Type  | Content Value                         |
+=============================================+===============+=======================================+
| min_iterations                              | integer       |                                       |
+---------------------------------------------+---------------+---------------------------------------+
| max_iterations                              | integer       |                                       |
+---------------------------------------------+---------------+---------------------------------------+
| max_preconditioner_lag_iterations           | integer       |                                       |
+---------------------------------------------+---------------+---------------------------------------+
| nonlinear_tolerance                         | exponential   |                                       |
+---------------------------------------------+---------------+---------------------------------------+
| limit_iterations                            | integer       |                                       |
+---------------------------------------------+---------------+---------------------------------------+
| nonlinear_iteration_damping_factor          | exponential   |                                       |
+---------------------------------------------+---------------+---------------------------------------+
| nonlinear_iteration_divergence_factor       | exponential   |                                       |
+---------------------------------------------+---------------+---------------------------------------+
| max_divergent_iterations                    | integer       |                                       |
+---------------------------------------------+---------------+---------------------------------------+
| initialize_with_darcy                       | boolean       | ``true, false``                       |
+---------------------------------------------+---------------+---------------------------------------+
| restart_tolerance_factor                    | exponential   |                                       |
+---------------------------------------------+---------------+---------------------------------------+
| restart_tolerance_relaxation_factor         | exponential   |                                       |
+---------------------------------------------+---------------+---------------------------------------+
| restart_tolerance_relaxation_factor_damping | exponential   |                                       |
+---------------------------------------------+---------------+---------------------------------------+
| preconditioner                              | string        | ``trilinos_ml, hypre_amg, block_ilu`` |
+---------------------------------------------+---------------+---------------------------------------+
| unstr_initialization                        | element block |                                       |
+---------------------------------------------+---------------+---------------------------------------+

Specifics about each ``preconditioner`` is defined in the `Unstr_preconditioners`_ section.

The ``unstr_initialization`` IS USED FOR SOMETHING.  If the ``unstr_initialization`` element is present, even without any subelements, initialization is turned on and default values are used.  Users should take care to only include the ``unstr_initialization`` element when its use is intended.  The ``unstr_initialization`` has the following subelements:

+-----------------------+---------------+---------------------------------------+
| Element Names         | Content Type  | Content Value                         |
+=======================+===============+=======================================+
| clipping_saturation   | exponential   |                                       |
+-----------------------+---------------+---------------------------------------+
| clipping_pressure     | exponential   |                                       |
+-----------------------+---------------+---------------------------------------+
| method                | string        | ``picard, darcy_solver``              |
+-----------------------+---------------+---------------------------------------+
| preconditioner        | string        | ``trilinos_ml, hypre_amg, block_ilu`` |
+-----------------------+---------------+---------------------------------------+
| linear_solver         | string        | ``aztec00``                           |
+-----------------------+---------------+---------------------------------------+
| error_control_options | string        | ``pressure``                          |
+-----------------------+---------------+---------------------------------------+
| convergence_tolerance | exponential   |                                       |
+-----------------------+---------------+---------------------------------------+
| max_iterations        | integer       |                                       |
+-----------------------+---------------+---------------------------------------+

Unstr_transient_controls
________________________

+-------------------------------------------------------+---------------+---------------------------------------+
| Element Names                                         | Content Type  | Content Value                         |
+=======================================================+===============+=======================================+
| min_iterations                                        | integer       | *default = 10*                        |
+-------------------------------------------------------+---------------+---------------------------------------+
| max_iterations                                        | integer       | *default = 15*                        |
+-------------------------------------------------------+---------------+---------------------------------------+
| limit_iterations                                      | integer       | *default = 20*                        |
+-------------------------------------------------------+---------------+---------------------------------------+
| nonlinear_tolerance                                   | exponential   | *default = 1.0e-5*                    |
+-------------------------------------------------------+---------------+---------------------------------------+
| nonlinear_iteration_damping_factor                    | exponential   | *default = 1.0*                       |
+-------------------------------------------------------+---------------+---------------------------------------+
| max_preconditioner_lag_iterations                     | integer       | *default = 5*                         |
+-------------------------------------------------------+---------------+---------------------------------------+
| max_divergent_iterations                              | integer       | *default = 3*                         |
+-------------------------------------------------------+---------------+---------------------------------------+
| nonlinear_iteration_divergence_factor                 | exponential   | *default = 1000.0*                    |
+-------------------------------------------------------+---------------+---------------------------------------+
| restart_tolerance_relaxation_factor                   | exponential   |                                       |
+-------------------------------------------------------+---------------+---------------------------------------+
| restart_tolerance_relaxation_factor_damping           | exponential   |                                       |
+-------------------------------------------------------+---------------+---------------------------------------+
| error_control_options                                 | string        | ``pressure, residual``                |
+-------------------------------------------------------+---------------+---------------------------------------+
| nonlinear_iteration_initial_guess_extrapolation_order | integer       |                                       |
+-------------------------------------------------------+---------------+---------------------------------------+
| preconditioner                                        | string        | ``trilinos_ml, hypre_amg, block_ilu`` |
+-------------------------------------------------------+---------------+---------------------------------------+
| initialize_with_darcy                                 | boolean       | | ``true, false``                     |
|                                                       |               | | *default = false*                   |
+-------------------------------------------------------+---------------+---------------------------------------+


Unstr_linear_solver
___________________

+----------------+--------------+---------------------------------------+
| Element Names  | Content Type | Content Value                         |
+================+==============+=======================================+
| method         | string       | ``gmres, pcg``                        |
+----------------+--------------+---------------------------------------+
| max_iterations | integer      |                                       |
+----------------+--------------+---------------------------------------+
| tolerance      | exponential  |                                       |
+----------------+--------------+---------------------------------------+
| preconditioner | string       | ``trilinos_ml, hypre_amg, block_ilu`` |
+----------------+--------------+---------------------------------------+


Unstr_nonlinear_solver
______________________

The nonlinear solver of choice is listed as the attribute ``name`` to the ``unstr_nonlinear_solver`` element.  The available options are: `"nka`", `"newton`", `"jfnk`", or `"newton_picard`".  Additional subelements are as follows:

+-------------------------+--------------+-----------------------------------------------+
| Element Names           | Content Type | Content Value                                 |
+=========================+==============+===============================================+
| modify_correction       | boolean      | | ``true, false``                             |
|                         |              | | *default = false*                           |
+-------------------------+--------------+-----------------------------------------------+
| update_upwind_frequency | string       | ``every_timestep, every_nonlinear_iteration`` |
+-------------------------+--------------+-----------------------------------------------+


Unstr_preconditioners
_____________________

Options for each available precondition are set in the ``unstr_preconditioners`` section.  The preconditioners assigned to each numerical solver are specified in the appropriate sections above.  Note that only one set of options may be specified for each precondition.  There is multiple solvers are assigned the preconditioner they will all utilize the same set of options.  The ``unstr_preconditioners`` element is defined as follows:

.. code-block:: xml

  <unstr_preconditioners>
      Required Elements: NONE
      Optional Elements: hypre_amg, trilinos_ml, block_ilu
  </unstr_preconditioners>

The subelements for the Hyper AMG preconditioner are as follows:

+-----------------------------+--------------+------------------------------------------+
| Element Names               | Content Type | Content Value                            |
+=============================+==============+==========================================+
| hypre_cycle_applications    | integer      |                                          |
+-----------------------------+--------------+------------------------------------------+
| hypre_smoother_sweeps       | integer      |                                          |
+-----------------------------+--------------+------------------------------------------+
| hypre_tolerance             | exponential  |                                          |
+-----------------------------+--------------+------------------------------------------+
| hypre_strong_threshold      | exponential  |                                          |
+-----------------------------+--------------+------------------------------------------+

The subelements for the Trilinos ML preconditioner are as follows:

+-----------------------------+--------------+------------------------------------------+
| Element Names               | Content Type | Content Value                            |
+=============================+==============+==========================================+
| trilinos_smoother_type      | string       | ``jacobi, gauss_seidel, ilu``            |
+-----------------------------+--------------+------------------------------------------+
| trilinos_threshold          | exponential  |                                          |
+-----------------------------+--------------+------------------------------------------+
| trilinos_smoother_sweeps    | integer      |                                          |
+-----------------------------+--------------+------------------------------------------+
| trilinos_cycle_applications | integer      |                                          |
+-----------------------------+--------------+------------------------------------------+

The subelements for the Block ILU preconditioner are as follows:

+-----------------------------+--------------+------------------------------------------+
| Element Names               | Content Type | Content Value                            |
+=============================+==============+==========================================+
| ilu_overlap                 | integer      |                                          |
+-----------------------------+--------------+------------------------------------------+
| ilu_relax                   | exponential  |                                          |
+-----------------------------+--------------+------------------------------------------+
| ilu_rel_threshold           | exponential  |                                          |
+-----------------------------+--------------+------------------------------------------+
| ilu_abs_threshold           | exponential  |                                          |
+-----------------------------+--------------+------------------------------------------+
| ilu_level_of_fill           | integer      |                                          |
+-----------------------------+--------------+------------------------------------------+

An example ``unstructured_controls`` section would look as the following:

.. code-block:: xml

       <unstructured_controls>
            <unstr_flow_controls>
                <discretization_method>fv-default</discretization_method>
                <rel_perm_method>upwind-darcy_velocity</rel_perm_method>
                <preconditioning_strategy>diffusion_operator</preconditioning_strategy>
            </unstr_flow_controls>
            <unstr_transport_controls>
                <algorithm>explicit first-order</algorithm>
                <sub_cycling>on</sub_cycling>
                <cfl>1</cfl>
            </unstr_transport_controls>
            <unstr_steady-state_controls>
                <min_iterations>10</min_iterations>
                <max_iterations>15</max_iterations>
                <limit_iterations>20</limit_iterations>
                <max_preconditioner_lag_iterations>5</max_preconditioner_lag_iterations>
                <nonlinear_tolerance>1.0e-5</nonlinear_tolerance>
                <error_control_options>pressure</error_control_options>
                <nonlinear_iteration_damping_factor>1</nonlinear_iteration_damping_factor>
                <nonlinear_iteration_divergence_factor>1000</nonlinear_iteration_divergence_factor>
                <max_divergent_iterations>3</max_divergent_iterations>
                <initialize_with_darcy>true</initialize_with_darcy>
                <restart_tolerance_relaxation_factor>1</restart_tolerance_relaxation_factor>
                <preconditioner>hypre_amg</preconditioner>
            </unstr_steady-state_controls>
            <unstr_transient_controls>
                <min_iterations>10</min_iterations>
                <max_iterations>15</max_iterations>
                <limit_iterations>20</limit_iterations>
                <nonlinear_tolerance>1.0e-5</nonlinear_tolerance>
                <nonlinear_iteration_damping_factor>1.0</nonlinear_iteration_damping_factor>
                <max_preconditioner_lag_iterations>5</max_preconditioner_lag_iterations>
                <max_divergent_iterations>3</max_divergent_iterations>
                <nonlinear_iteration_divergence_factor>1000</nonlinear_iteration_divergence_factor>
                <restart_tolerance_relaxation_factor>1</restart_tolerance_relaxation_factor>
                <error_control_options>pressure,residual</error_control_options>
                <preconditioner>hypre_amg</preconditioner>
                <initialize_with_darcy>true</initialize_with_darcy>
            </unstr_transient_controls>
            <unstr_preconditioners>
                <hypre_amg>
                    <hypre_cycle_applications>5</hypre_cycle_applications>
                    <hypre_smoother_sweeps>3</hypre_smoother_sweeps>
                    <hypre_tolerance>0.0</hypre_tolerance>
                    <hypre_strong_threshold>0.5</hypre_strong_threshold>
                </hypre_amg>
            </unstr_preconditioners>
            <unstr_linear_solver>
                <method>gmres</method>
                <max_iterations>100</max_iterations>
                <tolerance>1.0e-16</tolerance>
                <preconditioner>hypre_amg</preconditioner>
            </unstr_linear_solver>
            <unstr_nonlinear_solver name="nka">
                <modify_correction>false</modify_correction>
                <update_upwind_frequency>every_timestep</update_upwind_frequency>
            </unstr_nonlinear_solver>
        </unstructured_controls>

Structured_controls
-------------------

The ``structured_controls`` sections specifies numerical control options for the structured solver.  Either the ``unstructured_controls`` or the ``structured_controls`` section headings are required.  However, no options within the sections are required.  The list of available options is as follows:

.. code-block:: xml

  <structured_controls>
      Required Elements: NONE
      Optional Elements: comments, str_time_step_controls, str_flow_controls, str_transport_controls, str_amr_controls
  </structured_controls>

The subsections ``str_flow_controls`` and  ``str_transient_controls`` specify options specific to those process kernals.  The ``str_time_step_controls`` specify options for controlling the time step based on performance of the nonlinear solvers.  The subsection ``str_amr_controls`` specify options for AMR, including those for gridding and distribution granularity of data in parallel.

Str_time_step_controls
______________________

``str_time_step_controls`` has the following elements

+-----------------------------------+---------------+------------------------------------------+
| Element Names                     | Content Type  | Content Value                            |
+===================================+===============+==========================================+
| comments                          | string        |                                          |
+-----------------------------------+---------------+------------------------------------------+
| min_iterations                    | integer       |  *default = 10*                          |
+-----------------------------------+---------------+------------------------------------------+
| max_iterations                    | integer       |  *default = 15*                          |
+-----------------------------------+---------------+------------------------------------------+
| limit_iterations                  | integer       |  *default = 20*                          |
+-----------------------------------+---------------+------------------------------------------+
| min_iterations_2                  | integer       |  *default = 2*                           |
+-----------------------------------+---------------+------------------------------------------+
| time_step_increase_factor         | exponential   |  *default = 1.6*                         |
+-----------------------------------+---------------+------------------------------------------+
| time_step_increase_factor_2       | exponential   |  *default = 10*                          |
+-----------------------------------+---------------+------------------------------------------+
| max_consecutive_failures_1        | integer       |  *default = 3*                           |
+-----------------------------------+---------------+------------------------------------------+
| time_step_retry_factor_1          | exponential   |  *default = 0.2*                         |
+-----------------------------------+---------------+------------------------------------------+
| max_consecutive_failures_2        | integer       |  *default = 4*                           |
+-----------------------------------+---------------+------------------------------------------+
| time_step_retry_factor_2          | exponential   |  *default = 0.01*                        |
+-----------------------------------+---------------+------------------------------------------+
| time_step_retry_factor_f          | exponential   |  *default = 0.001*                       |
+-----------------------------------+---------------+------------------------------------------+
| max_num_consecutive_success       | integer       |  *default = 0*                           |
+-----------------------------------+---------------+------------------------------------------+
| extra_time_step_increase_factor   | exponential   |  *default = 10*                          |
+-----------------------------------+---------------+------------------------------------------+
| limit_function_evals              | integer       |  *default = 1000000*                     |
+-----------------------------------+---------------+------------------------------------------+
| do_grid_sequence                  | boolean       | ``true, false`` (*default = true*)       |
+-----------------------------------+---------------+------------------------------------------+
| grid_sequence_new_level_dt_factor | element block |  *see below*                             |
+-----------------------------------+---------------+------------------------------------------+

The element ``grid_sequence_new_level_dt_factor`` is an element block listing a series of dt_factors, one for each level.

Str_flow_controls
_________________

``str_flow_controls`` has the following elements

+-----------------------------------+---------------+------------------------------------------+
| Element Names                     | Content Type  | Content Value                            |
+===================================+===============+==========================================+
| comments                          | string        |                                          |
+-----------------------------------+---------------+------------------------------------------+
| petsc_options_file                | string        | *default = .petsc*                       |
+-----------------------------------+---------------+------------------------------------------+
| max_ls_iterations                 | integer       | *default = 10*                           |
+-----------------------------------+---------------+------------------------------------------+
| ls_reduction_factor               | exponential   | *default = 0.1*                          |
+-----------------------------------+---------------+------------------------------------------+
| min_ls_factor                     | exponential   | *default = 1.e-8*                        |
+-----------------------------------+---------------+------------------------------------------+
| ls_acceptance_factor              | exponential   | *default = 1.4*                          |
+-----------------------------------+---------------+------------------------------------------+
| monitor_line_search               | integer       | *default = 0*                            |
+-----------------------------------+---------------+------------------------------------------+
| monitor_linear_solve              | integer       | *default = 0*                            |
+-----------------------------------+---------------+------------------------------------------+
| use_fd_jac                        | boolean       | ``true, false`` (*default = true*)       |
+-----------------------------------+---------------+------------------------------------------+
| perturbation_scale_for_J          | exponential   | *default = 1.e-8*                        |
+-----------------------------------+---------------+------------------------------------------+
| use_dense_Jacobian                | boolean       | ``true, false`` (*default = false*)      |
+-----------------------------------+---------------+------------------------------------------+
| upwind_krel                       | string        | | ``upwind-darcy_velocity``,             |
|                                   |               | | ``other-arithmetic_average``,          |
|                                   |               | | ``other-harmonic_average``             |
+-----------------------------------+---------------+------------------------------------------+
| pressure_maxorder                 | integer       | *default = 3*                            |
+-----------------------------------+---------------+------------------------------------------+
| scale_solution_before_solve       | boolean       | ``true, false`` (*default = true*)       |
+-----------------------------------+---------------+------------------------------------------+
| semi_analytic_J                   | boolean       | ``true, false`` (*default = false*)      |
+-----------------------------------+---------------+------------------------------------------+
| atmospheric_pressure              | exponential   | *default = 1011325 (Pa)*                 |
+-----------------------------------+---------------+------------------------------------------+

Str_transport_controls
______________________

``str_transport_controls`` has the following elements

+-----------------------------------+---------------+------------------------------------------+
| Element Names                     | Content Type  | Content Value                            |
+===================================+===============+==========================================+
| comments                          | string        |                                          |
+-----------------------------------+---------------+------------------------------------------+
| max_n_subcycle_transport          | integer       | *default = 20*                           |
+-----------------------------------+---------------+------------------------------------------+
| cfl                               | exponential   | *default = 1*                            |
+-----------------------------------+---------------+------------------------------------------+

Str_amr_controls
________________

``str_amr_controls`` has the following elements

+-----------------------------------+------------------+-----------------------------------------------+
| Element Names                     | Content Type     | Content Value                                 |
+===================================+==================+===============================================+
| comments                          | string           |                                               |
+-----------------------------------+------------------+-----------------------------------------------+
| amr_levels                        | integer          | *default = 1*                                 |
+-----------------------------------+------------------+-----------------------------------------------+
| refinement_ratio                  | list of integers | *default = 2*                                 |
+-----------------------------------+------------------+-----------------------------------------------+
| do_amr_subcycling                 | boolean          | ``true, false`` *(default = true)*            |
+-----------------------------------+------------------+-----------------------------------------------+
| regrid_interval                   | list of integers | *default = 2*                                 |
+-----------------------------------+------------------+-----------------------------------------------+
| blocking_factor                   | list of integers | *default = 2*                                 |
+-----------------------------------+------------------+-----------------------------------------------+
| number_error_buffer_cells         | list of integers | *default = 1*                                 |
+-----------------------------------+------------------+-----------------------------------------------+
| max_grid_size                     | list of integers | *default = 64*                                |
+-----------------------------------+------------------+-----------------------------------------------+
| refinement_indicator              | element block    | *(see below)*                                 |
+-----------------------------------+------------------+-----------------------------------------------+


The user may define 1 or more refinement indicators.  Each refinement indicator is specified using the element block ``refinement_indicator`` with an attribute ``name`` to name the indicator.  The ``refinement_indicator`` has the following elements

+-----------------------------------+------------------+-----------------------------------------------+
| Element Names                     | Content Type     | Content Value                                 |
+===================================+==================+===============================================+
| field_name                        | string           |                                               |
+-----------------------------------+------------------+-----------------------------------------------+
| regions                           | string           |                                               |
+-----------------------------------+------------------+-----------------------------------------------+
| max_refinement_level              | integer          |                                               |
+-----------------------------------+------------------+-----------------------------------------------+
| start_time                        | exponential      |                                               |
+-----------------------------------+------------------+-----------------------------------------------+
| end_time                          | exponential      |                                               |
+-----------------------------------+------------------+-----------------------------------------------+
| | choose 1 of the following       | |                | |                                             |
| | value_greater                   | | exponential    | |                                             |
| | value_less                      | | exponential    | |                                             |
| | adjacent_difference_greater     | | exponential    | |                                             |
| | inside_region                   | | boolean        | | ``true, false``                             |
+-----------------------------------+------------------+-----------------------------------------------+

Mesh
====

Amanzi supports both structured and unstructured numerical solution approaches.  This flexibility has a direct impact on the selection and design of the underlying numerical algorithms, the style of the software implementations, and, ultimately, the complexity of the user-interface. The type of simulation is specified in the root tag ``amanzi_input``.  The ``mesh`` element varies slightly depending on whether the simulation type is ``structured`` or ``unstructured`` but is required for both.  For `"unstructured`", the ``mesh`` element specifies the internal mesh framework to be utilized and whether the mesh is to be internal generated or read in from an Exodus II file.  The default mesh framework is MSTK.  The other available frameworks are stk::mesh and simple (in serial). For `"structured`", the ``mesh`` element, specifies how the mesh is to be internally generated.

To internally generate a mesh the ``mesh`` element takes the following form.  The mesh framework attribute only applies to the `"unstructured`" and therefore is skipped for `"structured`" simulations.

Also, for parallel unstructured meshes, it is possible to choose a Partitioner from the available options, `"metis"`, `"zoltan_graph"` and `"zoltan_rcb"`. `"metis"` and `"zoltan_graph"` perform a graph partitioning of the mesh with no regard to the geometry of the mesh. `"zoltan_rcb"` partitions meshes using Recursive Coordinate Bisection which can lead to better partitioning in meshes that are thin in a particular direction. Additionally, the use of `"zoltan_rcb"` with the MSTK framework triggers an option to detect columns of elements in a mesh and adjust the partitioning such that no column is split over multiple partitions. If no partitioner is specified, a default method is used (`"metis"`).

.. code-block:: xml

   <mesh framework=["mstk"|"stk::mesh"|"simple"]>
      <comments> May be included in the Mesh element </comments>
      <dimension>3</dimension>
      <partitioner>"some partitioner keyword"</partitioner>
      <generate>
         <number_of_cells nx = "integer value"  ny = "integer value"  nz = "integer value"/>
         <box  low_coordinates = "x_low,y_low,z_low" high_coordinates = "x_high,y_high,z_high"/>
      </generate>
   </mesh>

For example:

.. code-block:: xml

  <mesh framework="mstk">
    <dimension>2</dimension>
    <partitioner>"zoltan_rcb"</partitioner>
    <generate>
      <number_of_cells nx="54" nz="60" />
      <box high_coordinates="216.0,120.0" low_coordinates="0.0, 0.0" />
    </generate>
  </mesh>

Currently Amanzi only read Exodus II mesh files for `"unstructured`" simulations.  An example ``mesh`` element would look as the following.

.. code-block:: xml

  <mesh framework="mstk"> 
    <comments> May be included in the Mesh element </comments>
    <dimension>3</dimension>
    <read>
      <file>mesh.exo</file>
      <format>exodus ii</format>
    </read>
  </mesh>

Note that the ``format`` content is case-sensitive and compared against a set of known and acceptable formats.  That set is [`"exodus ii`",`"exodus II`",`"Exodus II`",`"Exodus ii`"].  The set of all such limited options can always be verified by checking the Amanzi schema file.

Regions
=======

Regions are geometrical constructs used in Amanzi to define subsets of the computational domain in order to specify the problem to be solved, and the output desired. Regions are commonly used to specify material properties, boundary conditions and observation domains. Regions may represent zero-, one-, two- or three-dimensional subsets of physical space. For a three-dimensional problem, the simulation domain will be a three-dimensional region bounded by a set of two-dimensional regions. If the simulation domain is N-dimensional, the boundary conditions must be specified over a set of regions are (N-1)-dimensional.

Amanzi automatically defines the special region labeled "All", which is the entire simulation domain. Under the "Structured" option, Amanzi also automatically defines regions for the coordinate-aligned planes that bound the domain, using the following labels: `"XLOBC`", `"XHIBC`", `"YLOBC`", `"YHIBC`", `"ZLOBC`", `"ZHIBC`".

The ``regions`` block is required.  Within the region block at least one regions is required to be defined.  Most users define at least one region the encompasses the entire domain.  The optional elements valid for both structured and unstructured include `"region`", `"box`", `"point`", and `"plane`".  As in other sections there is also an options ``comments`` element.

The elements ``box``, ``point``, and ``plane`` allow for in-line description of regions.  The ``region`` element uses a subelement to either define a `"box`" or `"plane`" region or specify a region file.  Below are further descriptions of these elements.

Additional regions valid only for unstructured are ``polygonal_surface`` and ``logical``.  Additional regions valid only for structured include ``polygon`` and ``ellipse`` in 2D and ``rotated_polygon`` and ``swept_polygon`` in 3D.

.. code-block:: xml

  <regions>
      Required Elements: NONE
      Optional Elements: comments, box, point, region, (unstructured only - polygonal_surface, logical), (structured 2D only - polygon, ellipse), (structured 3D only - rotated_polygon, swept_polygon)
  </regions>

The elements box and point allow for in-line description of regions.  The region element uses a subelement to either define a box region or specify a region file.  

Box
---

A box region region is defined by a low corner coordinates and high corner coordinates.

.. code-block:: xml

  <box  name="box name" low_coordinates = "x_low,y_low,z_low" high_coordinates = "x_high,y_high,z_high"/>

Point
-----

A point region region is defined by a point coordinates.

.. code-block:: xml

  <point name="point name" coordinate = "x,y,z" />

Plane
-----

A plane region is defined by a point on the plane and the normal direction of the plane

.. code-block:: xml

  <plane name="plane name" location="x,y,z" normal="dx,dy,dz" tolerance="optional exp"/> 

The attribute ``tolerance`` is optional.  This value prescribes a tolerance for determining the cell face centroids that lie on the defined plane.

Region
------

A region allows for a box region, a point region, or a region file to be defined.

.. code-block:: xml

  <region name="Name of Region">
      Required Elements: 1 of the following - region_file, box, point  
      Optional Elements: comments
  </region>

A region is define as describe above.  A file is define as follows.


.. code-block:: xml

  <region_file name="filename" type=["color"|"labeled set"] format=["exodus ii"] entity=["cell"|"face"] label="integer"/>

Currently color functions and labeled sets can only be read from Exodus II files.  This will likely be the same file specified in the ``mesh`` element.  PLEASE NOTE the values listed within [] for attributes above are CASE SENSITIVE.  For many attributes within the Amanzi Input Schema the value is tested against a limited set of specific strings.  Therefore an user generated input file may generate errors due to a mismatch in cases.  Note that all specified names within this schema use lower case.

Polygonal_Surface [U]
---------------------

A polygonal_surface region is used to define a bounded planar region and is specified by the number of points and a list of points.  The points must be listed in order and this ordering is maintained during input translation.  This region type is only valid for the unstructured algorithm.

.. code-block:: xml

    <polygonal_surface name="polygon name" num_points="3" tolerance="optional exp">
      <point> X, Y, Z </point>
      <point> X, Y, Z </point>
      <point> X, Y, Z </point>
    </polygonal_surface>

The attribute ``tolerance`` is optional.  This value prescribes a tolerance for determining the cell face centroids that lie on the defined plane.

Logical
-------

Logical regions are compound regions formed from other primitive type regions using boolean operations. Supported operators are union, intersection, subtraction and complement.  This region type is only valid for the unstructured algorithm.


.. code-block:: xml

    <logical  name="logical name" operation = "union | intersection | subtraction | complement" region_list = "region1, region2, region3"/>


Polygon [S]
-----------

A polygon region is used to define a bounded planar region and is specified by the number of points and a list of points.  The points must be listed in order and this ordering is maintained during input translation.  This region type is only valid for the structured algorithm in 2D.

.. code-block:: xml

    <polygon name="polygon name" num_points="3">
      <point> X, Y </point>
      <point> X, Y </point>
      <point> X, Y </point>
    </polygon>

Ellipse [S]
-----------

An ellipse region is used to define a bounded planar region and is specified by a center and X and Y radii.  This region type is only valid for the structured algorithm in 2D.

.. code-block:: xml

    <ellipse name="polygon name" num_points="3">
      <center> X, Y </center>
      <radius> radiusX, radiusY </radius>
    </ellipse>

Rotated Polygon [S]
-------------------

A rotated_polygon region is defined by a list of points defining the polygon, the plane in which the points exist, the axis about which to rotate the polygon, and a reference point for the rotation axis.  The points listed for the polygon must be in order and the ordering will be maintained during input translation. This region type is only valid for the structured algorithm in 3D.

.. code-block:: xml

    <rotated_polygon name="rotated_polygon name">
        <vertex> X, Y, Z </vertex>
        <vertex> X, Y, Z </vertex>
        <vertex> X, Y, Z </vertex>
        <xyz_plane> XY | YZ | XZ </xyz_plane>
        <axis> X | Y | Z </axis>
        <reference_point> X, Y </reference_point>
    </rotated_polygon>

Swept Polygon [S]
-----------------

A swept_polygon region is defined by a list of points defining the polygon, the plane in which the points exist, the extents (min,max) to sweep the polygon normal to the plane.  The points listed for the polygon must be in order and the ordering will be maintained during input translation. This region type is only valid for the structured algorithm in 3D.

.. code-block:: xml

    <swept_polygon name="swept_polygon name">
        <vertex> X, Y, Z </vertex>
        <vertex> X, Y, Z </vertex>
        <vertex> X, Y, Z </vertex>
        <xyz_plane> XY | YZ | XZ </xyz_plane>
        <extent_min> exponential </extent_min>
        <extent_max> exponential </extent_max>
    </swept_polygon>

Geochemistry
============

Geochemistry allows users to define a reaction network and constraints to be associated with species defined under the ``dissolved_components`` section of the ``phases`` block.  Amanzi provides access to an internal geochemical engine as well as the Alquimia interface.  The Alquimia interface provides access to third-party geochemistry engines.  Currently available through Alquimia is the PFloTran engine. The user may specify engine specific information using the appropriate subelement.

.. code-block:: xml

  <geochemistry>
      Required Elements: NONE
      Optional Elements: verbosity, constraints
  </geochemistry>

Verbosity
---------

The ``verbosity`` element sets the verbosity for the geochemistry engine.  Available options are silent, terse, verbose, warnings, and errors.

Constraints
-----------

The ``constraints`` block is a list of ``constraint`` subelements identifying geochemical constraints and any relevant minerals for the reaction network.  Currently utilized by the PFloTran engine only.  If the attribute ``input_filename`` is missing from the ``process_kernels`` subelement ``chemistry``, Amanzi will automatically generating the PFloTran engine inputfile including the constraints defined here.  The constraints named and/or defined here can be referenced in the ``initial_conditions`` and ``boundary_conditions`` blocks.

* Each ``constraint`` has a ``name`` attribute.  If the user is providing the PFloTran input file, the name must match a constraint defined in the file.  Otherwise, the subelements defining the constraint must be provided and Amanzi will generate a constraint using this name. 

Individual constraints can have an unbounded number of chemical constraints defined under it.  The possible constraints are as follows.

  * Primary constraints are specified using the element ``primary``.  Attributes include ``name`` the name of the primary species, ``type`` the constraint type, and ``value`` the initial value to be used. For constraints based on equilibrium with a specific mineral or gas, an additional attribute specifying the mineral or gas is expected, ``mineral`` or ``gas`` respectively.  The table below lists the constraint types, which attributes are requires, and the corresponding value of the attribute ``type``.  Note, for non-reactive species/solutes, use the type "total".

  * Mineral constraints are specified using the element ``mineral``.  Attributes include ``name`` the name of the mineral, ``volume_fraction`` the volume fraction, and ``surface_area`` the specific surface area.


+------------------+---------------------+----------------+
| Constraint Type  | Required Attributes | ``type`` Value |
+==================+=====================+================+
| | Free ion       | | name              | free_ion       |
| | concentration  | | value             |                |
|                  | | type              |                |
+------------------+---------------------+----------------+
| | pH             | | name              | pH             |
|                  | | value             |                |
|                  | | type              |                |
+------------------+---------------------+----------------+
| | Total aquesous | | name              | total          |
| | concentration  | | value             |                |
|                  | | type              |                |
+------------------+---------------------+----------------+
| | Total aquesous | | name              | total+sorbed   |
| | + sorbed       | | value             |                |
| | concentration  | | type              |                |
+------------------+---------------------+----------------+
| | Charge balance | | name              | charge         |
|                  | | value             |                |
|                  | | type              |                |
+------------------+---------------------+----------------+
| | Concentration  | | name              | mineral        |
| | based on       | | value             |                |
| | mineral        | | type              |                |
|                  | | mineral           |                |
+------------------+---------------------+----------------+
| | Concentration  | | name              | gas            |
| | based on       | | value             |                |
| | mineral        | | type              |                |
|                  | | gas               |                |
+------------------+---------------------+----------------+

An example of a fully specified constraint is as follows.

.. code-block:: xml

  <constraints>
    <constraint name="initial">
        <primary name="Tc-99"   value="1e-3" type="total"/>
        <primary name="H2O"     value="1e-9"   type="mineral" mineral="Calcite"/>
        <primary name="CO2(aq)" value="1e-9"   type="gas" gas="CO2"/>
        <mineral name="Calcite" volume_fraction="1e-3" surface_area ="1e-5"/>
    </constraint>
  </constraints>

Note, if the user has provided a PFloTran input file, all that is required is the following,

.. code-block:: xml

  <constraints>
    <constraint name="initial"/>
  </constraints>

Any additional information provided is for the user's reference and will be ignored by Amanzi.

Materials
=========

The ``material`` in this context is meant to represent the media through with fluid phases are transported. In the literature, this is also referred to as the "soil", "rock", "matrix", etc. Properties of the material must be specified over the entire simulation domain, and is carried out using the Region constructs defined above. For example, a single material may be defined over the "All" region (see above), or a set of materials can be defined over subsets of the domain via user-defined regions. If multiple regions are used for this purpose, they should be disjoint, but should collectively tile the entire domain. The ``materials`` block is required.

Material
--------

Within the Materials block an unbounded number of ``material`` elements can be defined.  Each material requires a label and has the following requirements.

.. code-block:: xml

  <material>
      Required Elements: mechanical_properties, permeability or hydraulic_conductivity, assigned_regions
      Optional Elements: comments, cap_pressure, rel_perm, sorption_isotherms, minerals, ion_exchange, surface_complexation 
  </material>
 
Mechanical_properties
---------------------

.. code-block:: xml

  <mechanical_properties>
      Required Elements: porosity (FILE OPTION NOT IMPLEMENTED) 
      Optional Elements: particle_density, specific_storage, specific_yield, dispersion_tensor, tortuosity
  </mechanical_properties>

* ``mechanical_properties`` has six elements that can be either values or specified as files.  It has the following requirements.

    * ``porosity`` is defined in-line using attributes.  It is specified in one of three ways: as a value between 0 and 1 using value="<value>", through a file using type="file" and filename="<filename>", or as a gslib file using type="gslib", parameter_file="<filename>", value="<value>" and (optionally) data_file="<filename>" (defaults to ``porosity_data``.  NOTE - FILE OPTION NOT IMPLEMENTED YET.

    * ``particle_density`` is defined in-line using attributes.  Either it is specified as a value greater than 0 using ``value`` or it specified through a file using ``filename`` and ``type``.  NOTE - FILE OPTION NOT IMPLEMENTED YET.

    * ``specific_storage`` is defined in-line using attributes.  Either it is specified as a value greater than 0 using ``value`` or it specified through a file using ``filename`` and ``type``.  NOTE - FILE OPTION NOT IMPLEMENTED YET.

    * ``specific_yield`` is defined in-line using attributes.  Either it is specified as a value using ``value`` or it specified through a file using ``filename`` and ``type``.  NOTE - FILE OPTION NOT IMPLEMENTED YET.

    * ``dispersion_tensor`` is defined in-line using attributes.  The attribute ``type`` is used to specify either the model to utilize of that a file is to be read.  The ``type`` options are: uniform_isotropic, burnett_frind, lichtner_kelkar_robinson, or file.  For ``uniform_isotropic`` values are specified using the attributes ``alpha_l`` and ``alpha_t``.  For ``burnett_frind`` values are specified using the attributes ``alpha_l``, ``alpha_th``, and ``alpha_tv``. For ``lichtner_kelkar_robinson`` values are specified using the attributes ``alpha_l`h", ``alpha_lv``, ``alpha_th``, and ``alpha_tv``.  For ``file`` the file name is specified using ``filename``.  NOTE - FILE OPTION NOT IMPLEMENTED YET.

    * ``tortuosity`` is defined in-line using attributes.  Either it is specified as a value using ``value`` or it specified through a file using ``filename`` and ``type``.  NOTE - FILE OPTION NOT IMPLEMENTED YET.


.. code-block:: xml

  <mechanical_properties>
      <porosity value="exponential"/>
      <particle_density value="exponential"/>
      <specific_storage value="exponential"/>
      <specific_yield value="exponential"/>
      <dispersion_tensor type="uniform_isotropic" "alpha_l="exponential" alpha_t="exponential"/>
      <tortuosity value="exponential"/>
  </mechanical_properties>

Assigned_regions
----------------

* ``assigned_regions`` is a comma separated list of region names for which this material is to be assigned.  Region names must be from the regions defined in the ``regions`` sections.  Region names can contain spaces.

.. code-block:: xml

    <assigned_regions>Region1, Region_2, Region 3</assigned_regions>

Permeability
------------

Permeability or hydraulic_conductivity must be specified but not both. If specified as constant values, permeability has the attributes ``x``, ``y``, and ``z``.  Permeability may also be extracted from the attributes of an Exodus II file, or generated as a gslib file.

.. code-block:: xml

  <permeability x="exponential" y="exponential" z="exponential" />
  or
  <permeability type="file" filename="file name" attribute="attribute name"/>
  or
  <permeability type="gslib" parameter_file="file name" value="exponential" data_file="file name"/>

Hydraulic_conductivity
----------------------

* ``hydraulic_conductivity`` is the hydraulic conductivity and has the attributes ``x``, ``y``, and ``z``. Permeability or hydraulic_conductivity must be specified but not both.

.. code-block:: xml

  <hydraulic_conductivity x="exponential" y="exponential" z="exponential" />
  or
  <hydraulic_conductivity type="gslib" parameter_file="file name" value="exponential" data_file="file name"/>

Cap_pressure
------------

*  ``cap_pressure`` is an optional element.  The available models are ``van_genuchten``, ``brooks_corey``, and ``none``.  The model name is specified in an attribute and parameters are specified in a subelement.  Model parameters are listed as attributes to the parameter element.

* ``van_genuchten`` parameters include ``alpha``, ``sr``, ``m``, and ``optional_krel_smoothing_interval``.  ``brooks_corey`` parameters include ``alpha``, ``sr``, ``m``, and ``optional_krel_smoothing_interval``.

.. code-block:: xml

  <cap_pressure model="van_genuchten | brooks_corey | none" >
      Required Elements: alpha, Sr, m (van_genuchten and brooks_corey only)
      Optional Elements: optional_krel_smoothing_interval (van_genuchten and brooks_corey only)
  </cap_pressure>

Rel_perm
--------

*  ``rel_perm`` is an optional element.  The available models are ``mualem``, ``burdine``, and ``none``.  The model name is specified in an attribute and parameters are specified in a subelement.  Model parameters are listed as attributes to the parameter element.

* ``mualem`` has no parameters.  ``burdine`` parameters include ``exp``.

.. code-block:: xml

  <rel_perm model="mualem | burdine | none )" >
      Required Elements: none 
      Optional Elements: exp (burdine only)
  </rel_perm>

Sorption_isotherms
------------------

*  ``sorption_isotherms`` is an optional element for providing Kd models and molecular diffusion values for individual solutes.  All non-reactive primaries or solutes should be listed under each material.  Values of 0 indicate that the primary is not present/active in the current material.  The available Kd models are `"linear`", `"langmuir`", and `"freundlich`".  Different models and parameters are assigned per solute in sub-elements through attributes. The Kd and molecular diffusion parameters are specified in subelements.

.. code-block:: xml

    <sorption_isotherms>
	<solute name="string" />
            Required Elements: none
            Optional Elements: kd_model
    </sorption_isotherms>

.
    * ``kd_model`` takes the following form:

.. code-block:: xml
 
    <sorption_isotherms>
	<primary name="string" />
            <kd_model model="linear|langmuir|freundlich" kd="Value" b="Value (langmuir only)" n="Value (freundlich only)" />
	</primary>
    </sorption_isotherms>
  
Minerals
--------

* For each mineral, the concentrations are specified using the volume fraction and specific surface area using the attributes ``volume_fraction`` and ``specific_surface_area`` respectively.  

.. code-block:: xml

       <minerals>
           <mineral name="Calcite" volume_fraction="0.1" specific_surface_area"1.0"/>
       </minerals>

Ion_exchange
------------

* The ``ion_exhange`` block, specified parameters for an ion exchange reaction.  Cations active in the reaction are grouped under the element ``cations``.  The attribute ``cec`` specifies the cation exchange capacity for the reaction.  Each cation is listed in a ``cation`` subelement with the attributes ``name`` and ``value`` to specify the cation name and the associated selectivity coefficient.

.. code-block:: xml

        <ion_exchange>
            <cations cec="750.0">
                <cation name="Ca++" value="0.2953"/>
                <cation name="Mg++" value="0.1666"/>
                <cation name="Na+" value="1.0"/>
            </cations>
        </ion_exchange>

Surface_complexation
--------------------

* The ``surface_complexation`` block specifies parameters for surface complexation reactions.  Individual reactions are specified using the ``site`` block.  It has the attributes ``density`` and ``name`` to specify the site density and the name of the site.  Note, the site name must match a surface complexation site in the database file without any leading characters, such as `>`.  The subelement ``complexes`` provides a comma seperated list of complexes.  Again, the names of the complexes must match names within the datafile without any leading characters.

.. code-block:: xml

        <surface_complexation>
            <site density="1.908e-3" name="FeOH_s">
                <complexes>FeOHZn+_s, FeOH2+_s, FeO-_s</complexes>
            </site>
            <site density="7.6355e-2" name="FeOH_w">
                <complexes>FeOHZn+_w, FeO-_w, FeOH2+_w</complexes>
            </site>
        </surface_complexation>
    
Process Kernels
===============

The ``process_kernels`` block specifies which PKs are active.  This block is required for a valid input file.

.. code-block:: xml

  <process_kernels>
      Required Elements: flow, transport, chemistry
      Optional Elements: comments
  </process_kernels>

For each process kernel the element ``state`` indicates whether the solution is being calculated or not.  

Flow
----

* ``flow`` has the following attributes, 
      
      * ``state`` = "on | off"

      *  ``model`` = " richards | saturated | constant" 

Currently three scenarios are available for calculated the flow field.  ``richards`` is a single phase, variably saturated flow assuming constant gas pressure.  ``saturated`` is a single phase, fully saturated flow.  ``constant`` is equivalent to a flow model of single phase (saturated) with the time integration mode of transient with static flow in the version 1.2.1 input specification.  This flow model indicates that the flow field is static so no flow solver is called during time stepping. During initialization the flow field is set in one of two ways: (1) A constant Darcy velocity is specified in the initial condition; (2) Boundary conditions for the flow (e.g., pressure), along with the initial condition for the pressure field are used to solve for the Darcy velocity.

Note:  Unstructured options ``discretization_method``,  ``rel_perm_method``, and ``preconditioning_strategy`` have been moved to the ``unstr_flow_controls`` section under ``numerical_controls``/

Transport
---------

* ``transport`` has the following attributes,
      
      * ``state`` = "on | off"

For ``transport`` the ``state`` must be specified.  

Note:  Unstructured options ``algorithm`` and ``sub_cycling`` have been moved to the ``unstr_transport_controls`` section under ``numerical_controls``/

Chemistry
---------

* ``chemistry`` has the following attributes,
      
      * ``state`` = "on | off"
      
      * ``engine`` = "amanzi | pflotran | crunchflow | none"

      * ``input_filename`` is the name of the chemistry engine input file (filename.in).  If this is omitted Amanzi will automatically generate this file.

      * ``database`` is the name of the chemistry reaction database file (filename.dat).   

For ``chemistry`` a combination of ``state`` and ``engine`` must be specified.  If ``state`` is `"off`" then ``engine`` is set to `"none`".  Otherwise the ``engine`` must be specified. 

Phases
======

Some general discussion of the ``Phases`` section goes here.

.. code-block:: xml

  <Phases>
      Required Elements: liquid_phase 
      Optional Elements: solid_phase
      Optional Elements: gas_phase [U]
  </Phases>

Liquid_phase
------------

* ``liquid_phase`` has the following elements

.. code-block:: xml

  <liquid_phase>
      Required Elements: viscosity, density
      Optional Elements: dissolved_components, eos [S]
  </liquid_phase>

Here is more info on the ``liquid_phase`` elements:

    * ``eos`` = "string" 

    * ``viscosity`` = "exponential"

    * ``density`` = "exponential"

    * ``dissolved_components`` has the elements

        * ``primaries`` 
          
        * ``secondaries``

The subelement ``primaries`` is used for specifying reactive and non-reactive primary species.  An unbounded number of subelements ``primary`` can be specified.  The text body of the element lists the name of the primary.  Note, the name of the primary must match a species in the database file.  The ``primary`` element has the following attributes:

    * ``coefficient_of_diffusion`` = "exponential", this is an optional attribute

    * ``first_order_decay_constant`` = "exponential", this is an optional attribute

    * ``forward_rate`` = "exponential", this is a required attribute when being used with non-reactive primaries/solutes and automatically generating the chemistry engine input file

    * ``backward_rate`` = "exponential", this is a required attribute when being used with non-reactive primaries/solutes and automatically generating the chemistry engine input file

The subelement ``secondaries`` is used for specifying secondaries species for reactive chemistry.  An unbounded number of sublements ``secondary`` can be specified.  The body of the element lists the name of the secondary species.  Note, the name of the secondary must match a species in the database file.


Solid_phase
-----------

* ``solid_phase`` has the following elements

.. code-block:: xml

  <solid_phase>
      Required Elements: minerals
      Optional Elements: NONE
  </solid_phase>

Here is more info on the ``solid_phase`` elements:

    * ``minerals`` has the element 

        * ``mineral`` which contains the name of the mineral. Note, the name of the mineral must match a species in the database file.

Initial Conditions
==================

Some general discussion of the ``initial_condition`` section goes here.

The ``initial_conditions`` section requires at least 1 and up to an unbounded number of ``initial_condition`` elements.  Each ``initial_condition`` element defines a single initial condition that is applied to one or more region.  The following is a description of the ``initial_condition`` element.

.. code-block:: xml

  <initial_condition>
      Required Elements: assigned_regions
      Optional Elements: liquid_phase (, comments, solid_phase - SKIPPED)
  </initial_condition>

Assigned_regions
----------------

* ``assigned_regions`` is a comma separated list of regions to apply the initial condition to.

Liquid_phase
------------

* ``liquid_phase`` has the following elements

.. code-block:: xml

  <liquid_phase>
      Required Elements: liquid_component
      Optional Elements: geochemistry_component
  </liquid_phase>

*  Here is more info on the ``liquid_component`` block:

    * ``uniform_pressure`` is defined in-line using attributes.  Uniform specifies that the initial condition is uniform in space.  Value specifies the value of the pressure.  
      
    * ``linear_pressure`` is defined in-line using attributes.  Linear specifies that the initial condition is linear in space.  Gradient specifies the gradient value in each direction in the form of a coordinate (grad_x, grad_y, grad_z).  Reference_coord specifies a reference location as a coordinate.  Value specifies the value of the pressure.
      
    * ``uniform_saturation`` is defined in-line using attributes.  See ``uniform_pressure`` for details.
      
    * ``linear_saturation`` is defined in-line using attributes. See ``linear_pressure`` for details.
      
    * ``velocity`` is defined in-line using attributes.  Specify the velocity is each direction using the appropriate attributes x, y, and z.

.. code-block:: xml

    <uniform_pressure name="some name" value="exponential" />
    <linear_pressure name="some name" value="exponential" reference_coord="coordinate" gradient="coordinate"/>
    <uniform_saturation name="some name" value="exponential" />
    <linear_saturation name="some name" value="exponential" reference_coord="coordinate" gradient="coordinate"/>
    <velocity name="some name" x="exponential" y="exponential" z="exponential"/>

*  Here is more info on the ``geochemistry_component`` block:

    * ``geochemistry_component`` appears once.  An unbounded number of subelements ``constraint`` are used specify geochemical constraints to be applied at the beginning of the simulation.  Each ``constraint`` has an attribute ``name``.  The specified constraint must be defined in the external geochemistry file and the name must match.

.. code-block:: xml

     <geochemistry>
         <constraint name = "initial"/>
     </geochemistry>


Solid_phase
-----------

* ``solid_phase`` has the following elements - Reminder this element has been SKIPPED

.. code-block:: xml

  <solid_phase>
      Required Elements: geochemistry - SKIPPED
      Optional Elements: mineral, geochemistry - BOTH SKIPPED 
  </solid_phase>

Here is more info on the ``solid_phase`` elements: - NOT IMPLEMENTED YET

    * ``mineral`` has the element - SKIPPED 

        * ``mineral`` which contains the name of the mineral

    * ``geochemistry`` is an element with the following subelement: NOT IMPLEMENTED YET

        * ``constraint`` is an element with the following attributes: ONLY UNIFORM, for now

Boundary Conditions
===================

Some general discussion of the ``boundary_condition`` section goes here.

The ``boundary_conditions`` section contains an unbounded number of ``boundary_condition`` elements.  Each ``boundary_condition`` element defines a single initial condition that is applied to one or more region.  The following is a description of the ``boundary_condition`` element.

.. code-block:: xml

  <boundary_condition>
      Required Elements: assigned_regions, liquid_phase
      Optional Elements: comments - SKIPPED
  </boundary_condition>

Assigned_regions
----------------

* ``assigned_regions`` is a comma separated list of regions to apply the initial condition to.

Liquid_phase
------------

* ``liquid_phase`` has the following elements

.. code-block:: xml

  <liquid_phase>
      Required Elements: liquid_component
      Optional Elements: geochemistry_component
  </liquid_phase>

*  Here is more info on the ``liquid_component`` elements:

    * ``inward_mass_flux`` is defined in-line using attributes.  The attributes include "function", "start", and "value". Function specifies linear or constant temporal functional form during each time interval.  Start is a series of time values at which time intervals start.  Value is the value of the ``inward_mass_flux`` during the time interval. 

    * ``outward_mass_flux`` is defined in-line using attributes.  See ``inward_mass_flux`` for details.

    * ``inward_volumetric_flux`` is defined in-line using attributes.  See ``inward_mass_flux`` for details.

    * ``outward_volumetric_flux`` is defined in-line using attributes.  See ``inward_mass_flux`` for details.

    * ``uniform_pressure`` is defined in-line using attributes.  Uniform refers to uniform in spatial dimension.  See ``inward_mass_flux`` for details.

    * ``linear_pressure`` is defined in-line using attributes.  Linear refers to linear in spatial dimension. Gradient_value specifies the gradient value in each direction in the form of a coordinate (grad_x, grad_y, grad_z).  Reference_point specifies a reference location as a coordinate.  Reference_value specifies a reference value for the boundary condition. 

    * ``seepage_face`` is defined in-line using attributes.  The attributes include "function", "start", and "value". Function specifies linear or constant temporal functional form during each time interval.  Start is a series of time values at which time intervals start.  inward_mass_flux is the value of the inward_mass_flux during the time interval.
 
    * ``hydrostatic`` is an element with the attributes below.  By default the coordinate_system is set to "absolute".  Not specifying the attribute will result in the default value being used.  The attribute submodel is optional.  If not specified the submodel options will not be utilized.

    * ``linear_hydrostatic`` is defined in-line using attributes.  Linear refers to linear in spatial dimension. Gradient_value specifies the gradient value in each direction in the form of a coordinate (grad_x, grad_y, grad_z).  Reference_point specifies a reference location as a coordinate.  Reference_water_table_height specifies a reference value for the water table.  Optionally, the attribute "submodel" can be used to specify no flow above the water table height.

    * ``no_flow`` is defined in-line using attributes.  The attributes include "function" and "start". Function specifies linear or constant temporal functional form during each time interval.  Start is a series of time values at which time intervals start.  

.. code-block:: xml

     <inward_mass_flux value="exponential" function="linear | constant" start="time" />
     <outward_mass_flux value="exponential" function="linear | constant" start="time" />
     <inward_volumetric_flux value="exponential" function="linear | constant" start="time" />
     <outward_volumetric_flux value="exponential" function="linear | constant" start="time" />
     <uniform_pressure name="some name" value="exponential" function="uniform | constant" start="time" />
     <linear_pressure name="some name" gradient_value="coordinate" reference_point="coordinate" reference_value="exponential" />
     <seepage_face name="some name" inward_mass_flux="exponential" function="linear | constant" start="time" />
     <hydrostatic name="some name" value="exponential" function="uniform | constant" start="time" coordinate_system="absolute | relative to mesh top" submodel="no_flow_above_water_table | none"/>
     <linear_hydrostatic name="some name" gradient_value="exponential" reference_point="coordinate" reference_water_table_height="exponential" submodel="no_flow_above_water_table | none"/>
     <no_flow function="linear | constant" start="time" />

*  Here is more info on the ``geochemistry_component`` elements:

    * ``constraint`` is an element with the following attributes: ONLY UNIFORM, for now
    * If function is not specified and there is a geochemical constraint of the given name in the 
      ``geochemistry`` top-level element, information for that constraint will be taken from the 
      geochemical engine.

.. code-block:: xml

     <constraint name="some name" start="time" function="constant"/>

Sources
=======

Sources are defined in a similar manner to the boundary conditions.  Under the tag ``sources`` an unbounded number of individual ``source`` elements can be defined.  Within each ``source`` element the ``assigned_regions`` and ``liquid_phase`` elements must appear.  Sources can be applied to one or more region using a comma separated list of region names.  Under the ``liquid_phase`` element the ``liquid_component`` element must be define.  An unbounded number of ``solute_component`` elements and one ``geochemistry`` element may optionally be defined.

Under the ``liquid_component`` and ``solute_component`` elements a time series of boundary conditions is defined using the boundary condition elements available in the table below.  Each component element can only contain one type of source.  Both elements also accept a *name* attribute to indicate the phase associated with the source.

.. code-block:: xml

  <sources>
      Required Elements: assigned_regions, liquid_phase
      Optional Elements: comments - SKIPPED
  </sources>

Assigned_regions
----------------

* ``assigned_regions`` is a comma separated list of regions to apply the source to.

Liquid_phase
------------

* ``liquid_phase`` has the following elements

.. code-block:: xml

  <liquid_phase>
      Required Elements: liquid_component
      Optional Elements: solute_component (, geochemistry - SKIPPED)
  </liquid_phase>

*  Here is more info on the ``liquid_component`` elements:

    * ``volume_weighted`` is defined in-line using attributes.  The attributes include "function", "start", and "value". Function specifies linear or constant temporal functional form during each time interval.  Start is a series of time values at which time intervals start.  Value is the value of the ``volume_weighted`` during the time interval. 

    * ``perm_weighted`` is defined in-line using attributes.  See ``volume_weighted`` for details.

*  Here is more info on the ``solute_component`` elements:

    * ``uniform_conc`` is defined in-line using attributes.  The attributes include "name", "function", "start", and "value". Name is the name of a previously defined solute. Function specifies linear or constant temporal functional form during each time interval.  Start is a series of time values at which time intervals start.  Value is the value of the ``uniform_conc`` during the time interval. 

    * ``flow_weighted_conc`` is defined in-line using attributes.  See ``uniform_conc`` for details.

    * ``diffusion_dominated_release`` is defined in-line using attributes.  The attributes include "name", "start", "total_inventory", "mixing_length", and "effective_diffusion_coefficient". Name is the name of a previously defined solute. Start is a series of time values at which time intervals start.  Value is the value of the ``diffusion_dominated_release`` during the time interval. 

Output
======

Output data from Amanzi is currently organized into four specific elements: ``Vis``, ``Checkpoint``, ``Observations``, and ``Walkabout Data``.  Each of these is controlled in different ways, reflecting their intended use.

* ``Vis`` is intended to represent snapshots of the solution at defined instances during the simulation to be visualized.  The ''vis'' element defines the naming and frequencies of saving the visualization files.  The visualization files may include only a fraction of the state data, and may contain auxiliary "derived" information (see *elsewhere* for more discussion).

* ``Checkpoint`` is intended to represent all that is necessary to repeat or continue an Amanzi run.  The specific data contained in a Checkpoint Data dump is specific to the algorithm options and mesh framework selected.  Checkpoint is special in that no interpolation is performed prior to writing the data files; the raw binary state is necessary.  As a result, the user is allowed to only write Checkpoint at the discrete intervals of the simulation. The ''checkpoint'' element defines the naming and frequencies of saving the checkpoint files.

* ``Observations`` is intended to represent diagnostic values to be returned to the calling routine from Amanzi's simulation driver.  Observations are typically generated at arbitrary times, and frequently involve various point samplings and volumetric reductions that are interpolated in time to the desired instant.  Observations may involve derived quantities (see discussion below) or state fields.  The ''observations'' element may define one or more specific ''observation''.

* ``Walkabout Data`` is intended to be used as input to the particle tracking software Walkabout.

NOTE: Each output type allows the user to specify the base_filename or filename for the output to be written to.  The string format of the element allows the user to specify the relative path of the file.  It should be noted that the Amanzi I/O library does not create any new directories.  Therefore, if a relative path to a location other than the current directory is specified Amanzi assumes the user (or the Agni controller) has already created any new directories.  If the relative path does not exist the user will see error messages from the HDF5 library indicating failure to create and open the output file.

Vis
---

The ''vis'' element defines the visualization file naming scheme and how often to write out the files.  Thus, the ''vis'' element has the following requirements

.. code-block:: xml

  <vis>
      Required Elements: base_filename, num_digits 
      Optional Elements: time_macros, cycle_macros
  </vis>

The *base_filename* element contains the text component of the how the visualization files will be named.  The *base_filename* is appended with an index number to indicate the sequential order of the visualization files.  The *num_digits* elements indicates how many digits to use for the index. See the about NOTE about specifying a file location other than the current working directory.

The presence of the ''vis'' element means that visualization files will be written out after cycle 0 and the final cycle of the simulation.  The optional elements *time_macros* or *cycle_macros* indicate additional points during the simulation at which visualization files are to be written out.  Both elements allow one or more of the appropriate type of macro to be listed.  These macros will be determine the appropriate times or cycles to write out visualization files.  See the `Definitions`_ section for defining individual macros.

The ``vis`` element also includes an optional subelement ``write_regions``.  This was primarily implemented for debugging purposes but is also useful for visualizing fields only on specific regions.  The subelement accepts an arbitrary number of subelements named ``field``, with attributes ``name`` (a string) and ``regions`` (a comma separated list of region names).  For each such subelement, a field will be created in the vis files using the name as a label.  The field will be initialized to 0, and then, for region list R1, R2, R3..., cells in R1 will be set to 1, cells in R2 will be set to 2, etc.  When regions in the list overlap, later ones in the list will take precedence.

(*EIB NOTE* - there should be a comment here about how the output is controlled, i.e. for each PK where do you go to turn on and off fields.  This will probably get filled in as the other sections fill out.)

Example:

.. code-block:: xml

  <vis>
     <base_filename>plot</base_filename>
     <num_digits>5</num_digits>
     <time_macros>Macro 1</time_macros>
     <write_regions>
       <field name="Region List 1" regions="R1, R2, R3" />
       <field name="Region List 2" regions="All" />
     </write_regions>
  </vis>


Checkpoint
----------

The ''checkpoint'' element defines the file naming scheme and frequency for writing out the checkpoint files.  As mentioned above, the user does not influence what is written to the checkpoint files.  Thus, the ''checkpoint'' element has the following requirements

.. code-block:: xml

  <checkpoint>
      Required Elements: base_filename, num_digits, cycle_macros
      Optional Elements: NONE
  </checkpoint>

The *base_filename* element contain the text component of the how the checkpoint files will be named.  The *base_filename* is appended with an index number to indicate the sequential order of the checkpoint files.  The *num_digits* elements indicates how many digits to use for the index. (*EIB NOTE* - verify if this is sequence index or iteration id)  Final the *cycle_macros* element indicates the previously defined cycle_macro to be used to determine the frequency at which to write the checkpoint files. Multiple cycle macros may be specified in a comma separated list. See the about NOTE about specifying a file location other than the current working directory.

NOTE: Previously the ''walkabout'' element had the subelement ''cycle_macro''.  All output is moving away from only allowing a single macro to be specified to allowing multiple macros as a comma separated list.  To ease the transition for users both singular and plural are currently accepted.  However, the singular option will go away in the future.  Please update existing input files to use ''cycle_macros''.

Example:

.. code-block:: xml

  <checkpoint>
     <base_filename>chk</base_filename>
     <num_digits>5</num_digits>
     <cycle_macros>Every_100_steps</cycle_macros>
  </checkpoint>


Observations
------------

The Observations element holds all the observations that the user is requesting from Amanzi, as well as meta data, such as the name of the file that Amanzi will write observations to.  The observations are collected by their phase. Thus, the ''observations'' element has the following requirements

.. code-block:: xml

   <observations>
     Required Elements: filename, liquid_phase
     Optional Elements: NONE
   </observations>

The *filename* element contains the filename for the observation output, and may include the full path.  Currently, all observations are written to the same file.  See the about NOTE about specifying a file location other than the current working directory.

The *liquid_phase* element requires that the name of the phase be specified as an attribute and at least one observation.  The observation element is named according to what is being observed.  The observations elements available are as follows:

.. code-block:: xml

     <liquid_phase name="Name of Phase (Required)">
       Required Elements: NONE 
       Optional Elements: integrated_mass [S], volumetric_water_content, gravimetric_water_content, aqueous_pressure, 
                          x_aqueous_volumetric_flux, y_aqueous_volumetric_flux, z_aqueous_volumetric_flux, material_id, 
                          hydraulic_head, aqueous_mass_flow_rate, aqueous_volumetric_flow_rate, aqueous_conc, drawdown,
                          water_table, solute_volumetric_flow_rate
     </liquid_phase>

The observation element identifies the field quantity to be observed.  Subelements identify the elements for a region, a model (functional) with which it will extract its source data, and a list of discrete times for its evaluation.  The observations are evaluated during the simulation and returned to the calling process through one of Amanzi arguments. The elements for each observation type are as follows:

.. code-block:: xml

   <observation_type>
     Required Elements: assigned_region, functional, time_macros or cycle_macros 
     Optional Elements: NONE
   </observation_type>

The only exceptions are aqueous_conc and solute_volumetric_flow_rate which both require a solute to be specified.  An additional subelement "solute" gives the name of the solute to calculate the aqueous concentration or volumetric flow rate for.  Be sure the name of given for the solute matches a defined solute elsewhere in the input file.  

NOTE: Previously individual observation elements had the subelement ''cycle_macro'' or ''time_macro''.  All output is moving away from only allowing a single macro to be specified to allowing multiple macros as a comma separated list.  To ease the transition for users both singular and plural are currently accepted.  However, the singular option will go away in the future.  Please update existing input files to use ''cycle_macros'' or ''time_macros''.


NOTE: Observation "water_table" calculates maximum position of the water table (using a piecewise linear interpolation of cell-based pressures) in a given volume region. If the region is saturated, the code returns *1.0e+99*. If the region is dry, the code returns *-1.0e+99*.

Example:

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

Walkabout [U]
-------------

The ''walkabout'' element defines the file naming scheme and frequency for writing out the walkabout files.  As mentioned above, the user does not influence what is written to the walkabout files only the writing frequency and naming scheme.  Thus, the ''walkabout'' element has the following requirements

.. code-block:: xml

  <walkabout>
      Required Elements: base_filename, num_digits, cycle_macros
      Optional Elements: NONE
  </walkabout>

The *base_filename* element contain the text component of the how the walkabout files will be named.  The *base_filename* is appended with an index number to indicate the sequential order of the walkabout files.  The *num_digits* elements indicates how many digits to use for the index.  Final the *cycle_macros* element indicates the previously defined cycle_macro to be used to determine the frequency at which to write the walkabout files. See the about NOTE about specifying a file location other than the current working directory.

NOTE: Previously the ''walkabout'' element had the subelement ''cycle_macro''.  All output is moving away from only allowing a single macro to be specified to allowing multiple macros as a comma separated list.  To ease the transition for users both singular and plural are currently accepted.  However, the singular option will go away in the future.  Please update existing input files to use ''cycle_macros''.

Example:

.. code-block:: xml

  <walkabout>
     <base_filename>chk</base_filename>
     <num_digits>5</num_digits>
     <cycle_macros>Every_100_steps</cycle_macros>
  </walkabout>

Misc
====

This section includes a collection of miscellaneous global options, specified as root tags.  Each of these options has a default behavior that will occur if the parameter is omitted.  If the parameter appears with no attributes specified, the default values for the attributes will be assumed.

.. code-block:: xml

  <echo_translated_input file_name="some name"/>

* Write the input data after internal translation. If this parameter is missing, the default XML
  file `"XXX_native_v7.xml`" is written, where `"XXX.xml`" is the name of the original Amanzi input file.
  If this parameter is present but attribute ``file_name`` is either omitted of empty string, no 
  translated file is written.


Full Example
============

.. code-block:: xml

  <amanzi_input type="unstructured" version="2.2.1">
    <echo_translated_input format="unstructured_native" file_name="oldspec.xml"/>

    <model_description name="example of full unstructured schema">
      <comments>Example input file </comments>
      <units>
        <length_unit>m</length_unit>
        <time_unit>s</time_unit>
        <mass_unit>kg</mass_unit>
        <conc_unit>molar</conc_unit>
      </units>
    </model_description>

    <definitions>
      <macros>
        <time_macro name="Observation Times">
          <time>1.2096E+10</time>
        </time_macro>
        <time_macro name="EveryMonth">
          <start>1956,y</start>
          <timestep_interval>1,m</timestep_interval>
          <stop>1988,y</stop>
        </time_macro>
        <cycle_macro name="Every100Cycles">
          <start>0</start>
          <timestep_interval>100</timestep_interval>
        </cycle_macro>
      </macros>
    </definitions>

    <process_kernels>
      <comments>Variably saturated flow</comments>
      <flow model="richards" state="on"/>
      <transport state="on"/>
      <chemistry engine="none" state="off"/>
    </process_kernels>

    <phases>
      <liquid_phase name="water">
        <eos>false</eos>
        <viscosity>1.002E-03</viscosity>
        <density>998.2</density>
        <dissolved_components>
            <primaries>
                <primary coefficient_of_diffusion="1e-9">Tc-99</primary>
            </primaries>
        </dissolved_components>
      </liquid_phase>
    </phases>

    <execution_controls>
      <verbosity level="high"/>
      <execution_control_defaults init_dt="1.0" method="picard" mode="steady" />
      <execution_control end="1956,y" mode="steady" start="0.0" init_dt="1000.0"/>
      <execution_control end="3000,y" mode="transient" start="1956,y" />
    </execution_controls>

    <numerical_controls>
      <unstructured_controls>

        <unstr_flow_controls>
          <preconditioning_strategy>linearized_operator</preconditioning_strategy>
        </unstr_flow_controls>

        <unstr_transport_controls>
          <algorithm>explicit first-order</algorithm>
          <sub_cycling>on</sub_cycling>
          <cfl>1</cfl>
        </unstr_transport_controls>

        <unstr_steady-state_controls>
          <min_iterations>10</min_iterations>
          <max_iterations>15</max_iterations>
          <limit_iterations>20</limit_iterations>
          <max_preconditioner_lag_iterations>5</max_preconditioner_lag_iterations>
          <nonlinear_tolerance>1.0e-5</nonlinear_tolerance>
          <nonlinear_iteration_damping_factor>1</nonlinear_iteration_damping_factor>
          <nonlinear_iteration_divergence_factor>1000</nonlinear_iteration_divergence_factor>
          <max_divergent_iterations>3</max_divergent_iterations>
  
          <unstr_initialization>
            <method>darcy_solver</method>
            <linear_solver>aztecoo</linear_solver>
          </unstr_initialization>
        </unstr_steady-state_controls>
  
        <unstr_transient_controls>
          <min_iterations>10</min_iterations>
          <max_iterations>15</max_iterations>
          <limit_iterations>20</limit_iterations>
          <max_preconditioner_lag_iterations>5</max_preconditioner_lag_iterations>
          <nonlinear_tolerance>1.0e-5</nonlinear_tolerance>
          <nonlinear_iteration_damping_factor>1</nonlinear_iteration_damping_factor>
          <nonlinear_iteration_divergence_factor>1000</nonlinear_iteration_divergence_factor>
          <max_divergent_iterations>3</max_divergent_iterations>
        </unstr_transient_controls>

        <unstr_linear_solver>
          <max_iterations>100</max_iterations>
          <tolerance>1e-20</tolerance>
        </unstr_linear_solver>

        <unstr_preconditioners>
          <hypre_amg />
          <trilinos_ml />
          <block_ilu />
        </unstr_preconditioners>

      </unstructured_controls>
    </numerical_controls>

    <mesh framework="mstk">
      <dimension>2</dimension>
      <generate>
        <number_of_cells nx="54" nz="60"/>
        <box high_coordinates="216.0,120.0" low_coordinates="0.0, 0.0"/>
      </generate>
    </mesh>

    <regions>
      <region name="All">
        <box high_coordinates="216.0, 120.0" low_coordinates="0.0, 0.0" />
      </region>
      <region name="Bottom Surface">
        <box high_coordinates="216.0, 0.0" low_coordinates="0.0, 0.0" />
      </region>
      <region name="RegionBottom">
        <box high_coordinates="216.0, 40.0" low_coordinates="0.0, 0.0" />
      </region>
      <region name="RegionMiddle">
        <box high_coordinates="216.0, 80.0" low_coordinates="0.0, 40.0" />
      </region>
      <region name="RegionTop">
        <box high_coordinates="216.0, 120.0" low_coordinates="0.0, 80.0" />
      </region>
      <region name="Recharge_Boundary_WestOfCribs">
        <box high_coordinates="72.0, 120.0" low_coordinates="0.0, 120.0" />
      </region>
      <region name="Crib_216-B-17">
        <box high_coordinates="80.0, 120.0" low_coordinates="72.0, 120.0" />
      </region>
      <region name="Recharge_Boundary_btwnCribs">
        <box high_coordinates="136.0, 120.0" low_coordinates="80.0, 120.0" />
      </region>
      <region name="Crib_216-B-18">
        <box high_coordinates="148.0, 120.0" low_coordinates="136.0, 120.0" />
      </region>
      <region name="Recharge_Boundary_EastOfCribs">
        <box high_coordinates="216.0, 120.0" low_coordinates="148.0, 120.0" />
      </region>
      <region name="Well">
        <box high_coordinates="112.0, 60.0" low_coordinates="108.0, 40.0" />
      </region>
    </regions>

    <materials>
      <material name="Facies_1">
        <mechanical_properties>
          <porosity value="0.4082"/>
        </mechanical_properties>
        <permeability x="1.9976E-12" z="1.9976E-13" />
        <cap_pressure model="van_genuchten">
          <parameters alpha="1.9467E-04" m="0.2294" sr="0.0"/>
        </cap_pressure>
        <rel_perm model="mualem"/>
        <assigned_regions>RegionMiddle</assigned_regions>
      </material>
  
      <material name="Facies_2">
        <mechanical_properties>
          <porosity value="0.2206"/>
        </mechanical_properties>
        <permeability x="6.9365E-11" z="6.9365E-12" />
        <cap_pressure model="van_genuchten">
          <parameters alpha="2.0260E-03" m="0.2136" sr="0.0"/>
        </cap_pressure>
        <rel_perm model="mualem"/>
        <assigned_regions>RegionBottom</assigned_regions>
      </material>
  
      <material name="Facies_3">
        <mechanical_properties>
          <porosity value="0.2340"/>
        </mechanical_properties>
        <permeability x="2.0706E-09" z="2.0706E-10" />
        <cap_pressure model="van_genuchten">
          <parameters alpha="2.0674E-03" m="0.3006" sr="0.0"/>
        </cap_pressure>
        <rel_perm model="mualem"/>
        <assigned_regions>RegionTop</assigned_regions>
      </material>
    </materials>

     <geochemistry>
        <verbosity>silent</verbosity>
        <constraints>
            <constraint name="initial">
                <primary name="Tc-99" type="total" value="0.0"/>
            </constraint>
            <constraint name="Crib_216-B-17">
                <primary name="Tc-99" type="total" value="1.881389E-06"/>
            </constraint>
            <constraint name="Crib_216-B-18">
                <primary name="Tc-99" type="total" value="2.266885E-06"/>
            </constraint>
        </constraints>
    </geochemistry>
    <initial_conditions>
      <initial_condition name="All">
        <assigned_regions>All</assigned_regions>
        <liquid_phase name="water">
          <liquid_component name="water">
            <linear_pressure name="IC1" value="101325.0" reference_coord="0.0, 0.0" gradient="0,-9793.5192" />
          </liquid_component>
          <geochemistry_component>
            <constraint name="initial"/>
          </geochemistry_component>
        </liquid_phase>
      </initial_condition>
    </initial_conditions>

    <boundary_conditions>
      <boundary_condition name="BC For Bottom Surface">
        <assigned_regions>Bottom Surface</assigned_regions>
        <liquid_phase name="water">
          <liquid_component name="water">
            <hydrostatic function="uniform" start="0.0" value="0.0"/>
          </liquid_component>
          <geochemistry_component>
            <constraint function="constant" name="initial" start="0.0 y"/>
          </geochemistry_component>
        </liquid_phase>
      </boundary_condition>
  
      <boundary_condition name="BC For Crib_216-B-17">
        <assigned_regions>Crib_216-B-17</assigned_regions>
        <liquid_phase name="water">
          <liquid_component name="water">
            <inward_volumetric_flux value="1.1071e-10" function="constant" start="0.0" />
            <inward_volumetric_flux value="0.00254022e-3" function="constant" start="6.17266656e+10" />
            <inward_volumetric_flux value="1.48666E-9" function="constant" start="6.1729344E10" />
            <inward_volumetric_flux value="1.48666E-9" function="constant" start="9.4672798E10" />
          </liquid_component>
          <geochemistry_component>
            <constraint function="constant" name="initial" start="0.0"/>
            <constraint function="constant" name="Crib_216-B-17" start="6.17266656e+10"/>
            <constraint function="constant" name="initial" start="6.1729344E10"/>
          </geochemistry_component>
        </liquid_phase>
      </boundary_condition>
  
      <boundary_condition name="BC For Crib_216-B-18">
        <assigned_regions>Crib_216-B-18</assigned_regions>
        <liquid_phase name="water">
          <liquid_component name="water">
            <inward_volumetric_flux value="1.1071E-10" function="constant" start="0.0" />
            <inward_volumetric_flux value="1.48666E-9" function="constant" start="6.17266656e+10" />
            <inward_volumetric_flux value="0.00330423e-3" function="constant" start="6.173178481E10" />
            <inward_volumetric_flux value="1.48666E-9" function="constant" start="6.173705521E10" />
            <inward_volumetric_flux value="1.48666E-9" function="constant" start="9.4672798E10" />
          </liquid_component>
          <geochemistry_component>
            <constraint function="constant" name="initial" start="0.0"/>
            <constraint function="constant" name="Crib_216-B-17" start="6.173178481E10"/>
            <constraint function="constant" name="initial" start="6.173705521E10"/>
          </geochemistry_component>
        </liquid_phase>
      </boundary_condition>
  
      <boundary_condition name="BC Rest">
        <assigned_regions>Recharge_Boundary_WestOfCribs,
                          Recharge_Boundary_btwnCribs,
                          Recharge_Boundary_EastOfCribs</assigned_regions>
        <liquid_phase name="water">
          <liquid_component name="water">
            <inward_volumetric_flux value="1.1071E-10" function="constant" start="0.0" />
            <inward_volumetric_flux value="1.48666E-9" function="constant" start="6.17266656e+10" />
          </liquid_component>
          <geochemistry_component>
            <constraint function="constant" name="initial" start="0.0 y"/>
          </geochemistry_component>
        </liquid_phase>
      </boundary_condition>
    </boundary_conditions>

    <output>
       <vis>
        <base_filename>plot</base_filename>
        <num_digits>5</num_digits>
        <cycle_macros>Every100Cycles</cycle_macros>
      </vis>
      <checkpoint>
        <base_filename>chk</base_filename>
        <num_digits>5</num_digits>
        <cycle_macros>Every100Cycles</cycle_macros>
      </checkpoint>
    </output>
  </amanzi_input>

