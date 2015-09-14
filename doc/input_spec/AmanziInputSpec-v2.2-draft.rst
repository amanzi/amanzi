==============================================
Amanzi XML Input Specification (Version 2.1.1)
==============================================

.. contents:: **Table of Contents**

Overview
========

The Amanzi simulator evolves a system of conservation equations for reacting flows in porous media, as detailed in the ASCEM report entitled "Mathematical Formulation Requirements and Specifications for the Process Models`" (hereafter referred to as the 'Model Requirements Document (MRD)'). The purpose of the present document is to specify the data required to execute Amanzi.  This specification should be regarded as a companion to the MRD, and parameterizations of the individual submodels are consistent between Amanzi, the MRD and this document. Where applicable, the relevant sections of the MRD are indicated.

All data required to execute Amanzi is specified within an XML formated file layed out according to the Amanzi input schema.  The current version of the Amanzi schema is located with the Amanzi source code repository.  The following discusses each section of the schema, its purpose and provides examples.  Further details can be found in the schema document amanzi.xsd.

Please note, many attributes within the XML list a limited set of specified values.  During validation of the input file or initialization of Amanzi the values in the user provided input file will be compared against the limited set provided in the XML Schema document.  Errors will occur is the values do not match exactly.  These values are CASE SENSITIVE.  The Amanzi schema has been designed will all LOWER CASE values.  Please note this when writing input file.  In particular, `"Exodus II`" will be evaluated as `"exodus ii`".

Amanzi Input
============

Here, the user specifies which version of the input the input file adheres to. The user also specifies the overall type of simulation being run.  Amanzi supports both structured and unstructured numerical solution approaches.  This flexibility has a direct impact on the selection and design of the underlying numerical algorithms, the style of the software implementations, and, ultimately, the complexity of the user-interface. The attribute `"type`" is used to selected between the following:

* `"Structured`": This instructs Amanzi to use BoxLib data structures and an associated paradigm to numerically represent the flow equations.  Data containers in the BoxLib software library, developed by CCSE at LBNL, are based on a hierarchical set of uniform Cartesian grid patches.  `"Structured`" requires that the simulation domain be a single coordinate-aligned rectangle, and that the "base mesh" consists of a logically rectangular set of uniform hexahedral cells.  This option supports a block-structured approach to dynamic mesh refinement, wherein successively refined subregions of the solution are constructed dynamically to track "interesting" features of the evolving solution.  The numerical solution approach implemented under the `"Structured`" framework is highly optimized to exploit regular data and access patterns on massively parallel computing architectures. 

* `"Unstructured`": This instructs Amanzi to use data structures provided in the Trilinos software framework.  To the extent possible, the discretization algorithms implemented under this option are largely independent of the shape and connectivity of the underlying cells.  As a result, this option supports an arbitrarily complex computational mesh structure that enables users to work with numerical meshes that can be aligned with geometrically complex man-made or geostatigraphical features.  Under this option, the user typically provides a mesh file that was generated with an external software package.  The following mesh file formats are currently supported: `"Exodus II`".  Amanzi also provides a rudimentary capability to generate regular meshes within the unstructured framework internally.

An exmample root tag of an input file would look like the following.

.. code-block:: xml

  <amanzi_input version="2.1.0" type="unstructured"/>


Model Description
=================

This allows the users to provide a name and general description of model being developed.  This is also the section in which the units for the problem are stored. This entire section is optional but encouraged as documentation.

.. code-block:: xml

  <model_description name="Name of Model" >
      Required Elements: NONE
      Optional Elements: comment, author, created, modified, model_id, description, purpose (units - NOT IMPLEMENTED YET)
  </model_description>

Units
-----

The ``units`` element defines the default units to be assumed for the entire input file.  Amanzi's internal default units are SI units.  Conversion from the default units specified in the ``units`` element to SI units will be done by Amanzi's input translator.  A time unit can be specified with any time value in the input file.  

``units`` has the optional elements of length, time, mass, and concentration.  Each of those in turn have their own structure.  The structures are as follows.

REMINDER - UNITS ARE NOT IMPLEMENTED YET

.. code-block:: xml

  <units>
      Required Elements: NONE
      Optional Elements: length_unit, time_unit, mass_unit, conc_unit
  </units>

.. code-block:: xml

  <length_unit>
      Required Elements: m or cm
      Optional Elements: NONE
  </length_unit>

.. code-block:: xml

  <time_unit>
      Required Elements: y, d, h, or s
      Optional Elements: NONE
  </time_unit>

.. code-block:: xml

  <mass_unit>
      Required Elements: kg
      Optional Elements: NONE
  </mass_unit>

.. code-block:: xml

  <conc_unit>
      Required Elements: molar
      Optional Elements: NONE
  </conc_unit>


Here is an overall example for the model description element.

.. code-block:: xml

  <model_description name="BC Cribs">
    <comments>Added section on units definition</comments>
    <model_name>What should be in this field; originally TBD</model_name>
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
      Optional Elements: named_times, constants, macros
  </definitions>

Named Times
-----------

Here the user can specify and name times to be used in other sections of the input file.   Note that if a name is repeated the last read value will be retained and all others will be overwritten.

.. code-block:: xml

  <named_times>
      Required Elements: NONE
      Optional Elements: time [S]
  </named_times>

A *time* requires the attributes `"name`" and `"value`".  If a unit is not specified with the value seconds is taken as the default.

.. code-block:: xml

  <named_times>
    <time  name="String" value="time,y|d|h|s"/>
  </named_times>

Constants
---------

Here the user can define and name constants to be used in other sections of the input file.  Note that if a name is repeated the last read value will be retained and all others will be overwritten.

.. code-block:: xml

  <constants>
      Required Elements: NONE
      Optional Elements: constant, time_constant, numerical_constant, area_mass_flux_constant 
  </constants>

A *constant* has three attributes `"name`", `"type`", and `"value`".  The user can provide any name, but not it should not be repeated anywhere within the input to avoid confusion.  The available types include: `"none`", `"time`", `"numerical`", and `"area_mass_flux`".  Values assigned to constants of type `"time`" can include known units, otherwise seconds will be assumed as the default.

.. code-block:: xml

    <constant name="String" type="none | time | numerical | area_mass_flux" value="constant_value"/>

A *time_constant* is a specific form of a constant.  It takes the attributes `"name`" and `"value`" where the value is a time (time unit optional).

.. code-block:: xml

    <time_constant  name="Name of Time"  value="value,y|d|h|s"/>

A *numerical_constant* is a specific form of a constant.  It takes the attributes `"name`" and `"value`". 

.. code-block:: xml

    <numerical_constant name="Name of Numerical Constant" value="value_constant"/>

A *area_mass_flux_constant* is a specific form of a constant.  It takes the attributes `"name`" and `"value`" where the value is an area mass flux. 

.. code-block:: xml

    <area_mass_flux_constant name="Name of Flux Constant" value="value_of_flux"/>

Macros
------

The ``macros`` section defines time, cycle, and variable macros.  These specify a list or interval for triggering an action, particularly, writing out visualization, checkpoint, walkabout, or observation files.  

.. code-block:: xml

  <constants>
      Required Elements: NONE
      Optional Elements: time_macro, cycle_macro, variable_macro [S]
  </constants>

Time_macro
__________

The *time_macro* requires an attribute `"name`".  The macro can then either take the form of one or more labeled time subelements or the subelements `"start`", `"timestep_interval`", and `"stop`" again containing labeled times.  A `"stop`" value of -1 will continue the cycle macro until the end of the simulation.  The labeled times can be time values assuming the default time unit of seconds or including a known time unit.

.. code-block:: xml

  <time_macro name="Name of Macro">
    <time>Value</time>
  </time_macro>

or 

.. code-block:: xml

  <time_macro name="Name of Macro">
    <start> TimeValue </start>
    <timestep_interval> TimeIntervalValue </timestep_interval>
    <stop> TimeValue | -1 </stop>
  </time_macro>


Cycle_macro
___________

The *cycle_macro* requires an attribute `"name`" and the subelements `"start`", `"timestep_interval`", and `"stop`" with integer values.  A `"stop`" value of -1 will continue the cycle macro until the end of the simulation.

.. code-block:: xml

  <cycle_macro name="Name of Macro">
    <start>Value</start>
    <timestep_interval>Value</timestep_interval>
    <stop>Value|-1</stop>
  </cycle_macro>

Variable_macro
______________

The *variable_macro* requires an attribute `"name`"  and one or more subelements `"variable`" containing strings.

.. code-block:: xml

  <variable_macro name="Name of Macro">
    <variable> VariableString </variable>
  </variable_macro>


An example *definitions* section would look as the following:

.. code-block:: xml

  <definitions>

    <constants>
      <constant name="zero"              type="none"           value="0.000"/>
      <constant name ="start"            type="time"           value="1956.0;y"/>
      <constant name ="B-18_release_end" type="time"           value ="1956.3288;y"/>
      <constant name="future_recharge"   type="area_mass_flux" value="1.48666E-6"/>
      <numerical_constant name="zero" value="0.000"/>
    </constants>

    <macros>

      <time_macro name="Macro 1">
        <time>6.17266656E10</time>
        <time>6.172982136E10</time>
        <time>6.173297712E10</time>
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


Execution Controls
==================

The ``execution_controls`` section defines the general execution of the Amanzi simulation.  Amanzi can execute in four modes: steady state, transient, transient with static flow, or initialize to a steady state and then continue to transient.  The transient with static flow mode does not compute the flow solution at each time step.  During initialization the flow field is set in one of two ways: (1) A constant Darcy velocity is specified in the initial condition; (2) Boundary conditions for the flow (e.g., pressure), along with the initial condition for the pressure field are used to solve for the Darcy velocity. At present this mode only supports the "Single Phase" flow model.

.. code-block:: xml
  
  <execution_controls>
      Required Elements: execution_control_defaults, execution_control (1 or more)
      Optional Elements: comments, verbosity
  </execution_controls>

Some explanation of each element goes here.

Verbosity
---------

The ``verbosity`` element specifies the level of output messages provided by Amanzi.  If not present, the default value of *medium* will be set.

.. code-block:: xml
  
  <verbosity level="none | low | medium | high | extreme" />
 
Note, for debugging purposes use level="extreme". 

Execution_control_defaults
--------------------------

The ``execution_control_defaults`` element specifies default values to be utilized when not specified in individual ``execution_control`` elements.   For a valid ``execution_controls`` section the ``execution_control_defaults`` element is *required*.  The attributes available are:

    * init_dt = "labeled_time" 
      
    * max_dt = "labeled_time" 
      
    * reduction_factor = "exponential"
      
    * increase_factor = "exponential"
      
    * mode = "steady | transient" 
      
    * method = "bdf1 | picard" [S]

.. code-block:: xml

  <execution_control_defaults init_dt="labeled_time" max_dt="labeled_time" reduction_factor="exponential" increase_factor="exponential" mode="steady | transient" method=" bdf1 | picard" />

Execution_control
-----------------

Individual time periods of the simulation are defined using ``execution_control`` elements.  For a steady state simulation, only one ``execution_control`` element will be defined.  However, for a transient simulation a series of controls may be defined during which different control values will be used.  For a valid ``execution_controls`` section at least one ``execution_control`` element must appear.  The attributes available are:
  
    * start = "string", this attribute is required
      
    * end = "labeled_time", this attribute us required for the final execution_control element 
      
    * init_dt = "labeled_time" 
      
    * max_dt = "labeled_time" 
      
    * reduction_factor = "exponential"
      
    * increase_factor = "exponential"
      
    * mode = "steady | transient" 
      
    * method = "bdf1 | picard" [S]

    * restart = "string", this attribute specifies the name of the Amanzi checkpoint file previously created and to be used to restart the current simulation

    * initialize = "string" [U], this attribute specifies the name of the Amanzi checkpoint file previously created and to be used to initialize the current simulation
       
    * max_cycles = "integer" 

.. code-block:: xml

  <execution_control start="string" end="labeled_time" init_dt="labeled_time" max_dt="labeled_time" reduction_factor="exponential" increase_factor="exponential" mode="steady | transient" method=" bdf1 | picard" restart="string"/>


Numerical Controls
==================

This section allows the user to define control parameters associated with the underlying numerical implementation.  The list of available options is lengthy.  However, none are required for a valid input file.  The ``numerical_controls`` section is divided up into the subsections: ``common_controls``, ``unstructured_controls``, and ``structured_controls``.  The ``common_controls`` section is currently empty.  However, in future versions controls that are common between the unstructured and structured executions will be moved to this section and given common terminology.

.. code-block:: xml

  <numerical_controls>
      Required Elements: NONE
      Optional Elements: comments, common_controls [S], unstructured_controls [U], structured_controls [S]
  </numerical_controls>

Unstructured_controls
---------------------

The ``unstructured_controls`` sections is divided in the subsections: ``unstr_steady-state_controls``, ``unstr_transient_controls``, ``unstr_linear_solver``, ``unstr_nonlinear_solver``, and ``unstr_chemistry_controls``.  The list of available options is as follows:

.. code-block:: xml

  <unstructured_controls>
      Required Elements: NONE
      Optional Elements: unstr_steady-state_controls, unstr_transient_controls, unstr_linear_solver, unstr_nonlinear_solver, unstr_chemistry_controls
  </unstructured_controls>

`"unstructured_controls`" contains options specific to the unstructured modes.  It has the following structure and elements

  * `"unstr_flow_controls`" specifies numerical controls for the flow process kernel available under the unstructured algorithm.  It has the following elements

    * `"discretization_method`" specifies the spatial discretization method. Is has type "string" (options: fv-default, fv-monotone, fv-multi_point_flux_approximation, fv-extended_to_boundary_edges, mfd-default, mfd-optimized_for_sparsity, mfd-support_operator, mfd-optimized_for_monotonicity, mfd-two_point_flux_approximation)

    * `"rel_perm_method`" defines a method for calculating the upwinded relative permeability. It has type "string" (options: upwind-darcy_velocity(default), upwind-gravity, upwind-amanzi, other-arithmetic_average, other-harmonic_average)

    * `"preconditioning_strategy`" = "string" (options: linearized_operator(default), diffusion_operator)

  * `"unstr_transport_controls`" specifies numerical controls for the transport process kernel available under the unstructured algorithm.  It has the following elements

    * `"algorithm`" = "string" (options: explicit first-order(default), explicit second-order, implicit upwind)

    * `"sub_cycling`" = "string" (options: off(default), on)

  * `"unstr_steady-state_controls`"  has the following elements

    * `"comments`" = "string" - SKIPPED

    * `"min_iterations`" = "integer"

    * `"max_iterations`" = "integer"

    * `"max_preconditioner_lag_iterations`" = "integer"

    * `"nonlinear_tolerance`" = "exponential"

    * `"unstr_initialization`"  has the following elements

        * `"method`" = "string" (options: picard, darcy_solver)

        * `"preconditioner`" = "string" (options: trilinos_ml, hypre_amg, block_ilu)

        * `"linear_solver`" = "string" (options: aztec00)

        * `"control_options`" = "string"

        * `"max_iterations`" = "integer"

        * `"clipping_saturation`" = "exponential"

        * `"clipping_pressure`" = "exponential"

        * `"convergence_tolerance`" = "exponential"

    * `"limit_iterations`" = "integer"

    * `"nonlinear_iteration_damping_factor`" = "exponential"

    * `"nonlinear_iteration_divergence_factor`" = "exponential"

    * `"max_divergent_iterations`" = "integer"

    * `"initialize_with_darcy`" = "boolean"

    * `"restart_tolerance_factor`" = "exponential"
 
    * `"restart_tolerance_relaxation_factor`" = "exponential"

    * `"restart_tolerance_relaxation_factor_damping`" = "exponential"

    * `"preconditioner`" requires an attribute `"name`". (options: trilinos_ml, hypre_amg, block_ilu) See below for subelements based on preconditioner name.

  * `"unstr_transient_controls`"  has the same elements as `"unstr_steady-state_controls`"

  * `"unstr_linear_solver`"  has the following elements

    * `"comments`" = "string" - SKIPPED
 
    * `"method`" = "string" (options: aztec00)

    * `"max_iterations`" = "integer"

    * `"tolerance`" = "exponential"

    * `"preconditioner`" requires an attribute `"name`" (options: trilinos_ml, hypre_amg, block_ilu) See below for subelements based on preconditioner name.

  * `"unstr_nonlinear_solver`"  has an attribute `"name`" (options: nka, newton, inexact newton)

  * `"unstr_chemistry_controls`"  has the following elements

    * `"chem_tolerance`" = "exponential" 
 
    * `"chem_max_newton_iterations`" = "integer"

  * `"unstr_preconditioners`" has a list of named presonditioners. Available preconditioners 
    are Trilinos' ML, Hypre's AMG, and block ILU.  Below are the structures for each named preconditioner.

    * `"trilinos_ml'`" has the following optional elements

      * `"trilinos_smoother_type`" = "string" (options: jacobi, gauss_seidel, ilu)

      * `"trilinos_threshold`" = "exponential" 

      * `"trilinos_smoother_sweeps`" = "integer"

      * `"trilinos_cycle_applications`" = "integer"

    * `"hypre_amg'`" has the following optional elements

      * `"hypre_cycle_applications`" = "integer"

      * `"hypre_smoother_sweeps`" = "integer"

      * `"hypre_tolerance`" = "exponential" 

      * `"hypre_strong_threshold`" = "exponential" 

    * `"block_ilu'`" has the following optional elements

      * `"ilu_overlap`" = "integer"

      * `"ilu_relax`" = "exponential"

      * `"ilu_rel_threshold`" = "exponential" 

      * `"ilu_abs_threshold`" = "exponential" 

      * `"ilu_level_of_fill`" = "integer" 

Structured_controls
---------------------

.. code-block:: xml

  <unstructured_controls>
      Required Elements: NONE
      Optional Elements: str_steady-state_controls, str_transient_controls, str_amr_controls, max_n_subcycle_transport
  </unstructured_controls>

`"structured_controls`" contains options specific to the structured modes.  It has the following structure and elements

* `"structured_controls`" 

  * `"petsc_options_file`"  is an element that specifies the name of a petsc control options file.  By default, the filename is .petsc and will be read in automatically if it exists.  This options allows the user to specify a file with an alternative name.
  
  * `"str_steady-state_controls`"  has the following elements
  
    * `"max_pseudo_time`" = "exponential"

    * `"limit_iterations`" = "integer"

    * `"min_iterations`" = "integer"

    * `"min_iterations_2`" = "integer"
  
    * `"time_step_increase_factor_2`" = "exponential"
  
    * `"max_consecutive_failures_1`" = "integer"
  
    * `"time_step_retry_factor_1`" = "exponential"
  
    * `"max_consecutive_failures_2`" = "integer"
  
    * `"time_step_retry_factor_2`" = "exponential"
  
    * `"time_step_retry_factor_f`" = "exponential"
  
    * `"max_num_consecutive_success`" = "integer"
  
    * `"extra_time_step_increase_factor`" = "exponential"
  
    * `"abort_on_psuedo_timestep_failure`" = "integer"
  
    * `"use_PETSc_snes`" = "bool"
  
    * `"limit_function_evals`" = "exponential"
  
    * `"do_grid_sequence`" = "bool"
  
    * `"grid_sequence_new_level_dt_factor`" takes a sequence of exponential values as subelements

        * `"dt_factor`" = "exponential"

  * `"str_transient_controls`"  has the following elements
  
    * `"max_ls_iterations`" = "integer"
  
    * `"ls_reduction_factor`" = "exponential"
  
    * `"min_ls_factor`" = "exponential"
  
    * `"ls_acceptance_factor`" = "exponential"
  
    * `"monitor_line_search`" = "integer"
  
    * `"monitor_linear_solve`" = "integer"
  
    * `"use_fd_jac`" = "bool"
  
    * `"perturbation_scale_for_J`" = "exponential"
  
    * `"use_dense_Jacobian`" = "bool"
  
    * `"upwind_krel`" = "bool"
  
    * `"pressure_maxorder`" = "integer"
  
    * `"scale_solution_before_solve`" = "bool"
  
    * `"semi_analytic_J`" = "bool"

    * `"cfl`" = "exponential"

  * `"str_amr_controls`"  has the following elements
  
    * `"amr_levels`" = "integer"
  
    * `"refinement_ratio`" takes a space separated list of integer values
  
    * `"do_amr_cubcycling`" = "bool"
  
    * `"regrid_interval`" takes a space separated list of integer values
  
    * `"blocking_factor`" takes space separated list of integer values

    * `"number_error_buffer_cells`" takes space separated list of integer values

    * `"max_grid_size`" = "integer"
  
    * `"refinement_indicators`" takes the following subelements
    
      * `"field_name`" = "string"
    
      * `"regions`" = "string"
    
      * `"max_refinement_level`" = "string"
    
      * `"start_time`" = "exponential"
    
      * `"end_time`" = "exponential"
      
      * The user may also specify exactly 1 of the following
      
        * `"value_greater`" = "exponential"
      
        * `"value_less`" = "exponential"
      
        * `"adjacent_difference_greater`" = "exponential"
      
        * `"inside_region`" = "bool"

Mesh
====

Amanzi supports both structured and unstructured numerical solution approaches.  This flexibility has a direct impact on the selection and design of the underlying numerical algorithms, the style of the software implementations, and, ultimately, the complexity of the user-interface. The type of simulation is specified in the root tag `"amanzi_input`".  The `"mesh`" element specifies the internal mesh framework to be utilized and whether the mesh is to be internal generated or read in from an Exodus II file.  The default mesh framework is MSTK.  The other available frameworks are stk::mesh and simple (in serial).

To internally generate a mesh the `"mesh`" element takes the following form.


.. code-block:: xml

   <mesh framework=["mstk"|"stk::mesh"|"simple"]>
      <comments> May be included in the Mesh element </comments>
      <dimension>3</dimension>
      <generate>
         <number_of_cells nx = "integer value"  ny = "integer value"  nz = "integer value"/>
         <box  low_coordinates = "x_low,y_low,z_low" high_coordinates = "x_high,y_high,z_high"/>
      </generate>

   </mesh>

For example:

.. code-block:: xml

  <mesh framework="mstk"> 
   <generate>
     <number_of_cells nx = "64"  ny = "56"  nz = "107"/>
     <box  low_coordinates = "0.0,0.0,0.0" high_coordinates = "320.0,280.0,107.0"/>
   </generate>
  </mesh>

Currently Amanzi only read Exodus II mesh files.  An example `"mesh`" element would look as the following.

.. code-block:: xml

  <mesh framework="mstk"> 
    <comments> May be included in the Mesh element </comments>
    <dimension>3</dimension>
    <read>
      <file>mesh.exo</file>
      <format>exodus ii</format>
    </read>
  </mesh>

Note that the `"format`" content is case-sensitive and compared against a set of known and acceptable formats.  That set is ["exodus ii","exodus II","Exodus II","Exodus ii"].  The set of all such limited options can always be verified by checking the Amanzi schema file.

Regions
=======

Regions are geometrical constructs used in Amanzi to define subsets of the computational domain in order to specify the problem to be solved, and the output desired. Regions are commonly used to specify material properties, boundary conditions and observation domains. Regions may represent zero-, one-, two- or three-dimensional subsets of physical space. For a three-dimensional problem, the simulation domain will be a three-dimensional region bounded by a set of two-dimensional regions. If the simulation domain is N-dimensional, the boundary conditions must be specified over a set of regions are (N-1)-dimensional.

Amanzi automatically defines the special region labeled "All", which is the entire simulation domain. Under the "Structured" option, Amanzi also automatically defines regions for the coordinate-aligned planes that bound the domain, using the following labels: "XLOBC", "XHIBC", "YLOBC", "YHIBC", "ZLOBC", "ZHIBC"

The ``regions`` block is required.  Within the region block no regions are required to be defined.  The optional elements valid for both structured and unstructured include ``region``, ``box``, ``point``, and ``plane``.  As in other sections there is also an options ``comments`` element.

The elements ``box``, ``point``, and ``plane`` allow for in-line description of regions.  The ``region`` element uses a subelement to either define a ``box`` or ``plane`` region or specify a region file.  Below are further descriptions of these elements.

Additional regions valid only for unstructured are ``polygonal_surface`` and ``logical``.  Additional regions valid only for structured include ``polygon`` and ``ellipse`` in 2D and ``rotated_polygon`` and ``swept_polygon`` in 3D.

.. code-block:: xml

  <regions>
      Required Elements: NONE
      Optional Elements: comments, box, point, region, (unstructured only - polygonal_surface, logical), (structured 2D only - polygon, ellipse), (structured 3D only - rotated_polygon, swept_polygon)
  </regions>

The regions block is required.  Within the region block no regions are required to be defined.  

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

A region allows for a box region or a region file.

.. code-block:: xml

  <region name="Name of Region">
      Required Elements: region  
      Optional Elements: comments
  </region>

A region is define as describe above.  A file is define as follows.


.. code-block:: xml

  <region_file name="filename" type=["color"|"labeled set"] format=["exodus ii"] entity=["cell"|"face"] label="integer"/>

Currently color functions and labeled sets can only be read from Exodus II files.  This will likely be the same file specified in the `"mesh`" element.  PLEASE NOTE the values listed within [] for attributes above are CASE SENSITIVE.  For many attributes within the Amanzi Input Schema the value is tested against a limited set of specific strings.  Therefore an user generated input file may generate errors due to a mismatch in cases.  Note that all specified names within this schema use lower case.

Polygonal_Surface [U]
---------------------

A polygonal_surface region is used to define a bounded planar region and is specified by the number of points and a list of points.  The points must be listed in order and this ordering is maintained during input translation.  This region type is only valid for the unstructured algorithm.

.. code-block:: xml

    <polygonal_surface name="polygon name" num_points="3" tolerance="optional exp">
      <point> (X, Y, Z) </point>
      <point> (X, Y, Z) </point>
      <point> (X, Y, Z) </point>
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
      <point> (X, Y) </point>
      <point> (X, Y) </point>
      <point> (X, Y) </point>
    </polygon>

Ellipse [S]
-----------

An ellipse region is used to define a bounded planar region and is specified by a center and X and Y radii.  This region type is only valid for the structured algorithm in 2D.

.. code-block:: xml

    <ellipse name="polygon name" num_points="3">
      <center> (X, Y) </center>
      <radius> (radiusX, radiusY) </radius>
    </ellipse>

Rotated Polygon [S]
-------------------

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

Swept Polygon [S]
-----------------

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
============

Geochemistry allows users to define a reaction network and constraints to be associated with solutes defined under the `"dissolved_components`" section of the `"phases`" block.

.. code-block:: xml

  <geochemistry>
      Required Elements: reaction_network [S], constraint [S]
  </geochemistry>

PFLOTRAN Chemistry
------------------

For geochemistry simulated through PFLOTRAN, the user defines a reaction network and constraints.  These are defined within the same or separate text files through PFLOTRAN's input specification (see the CHEMISTRY and CONSTRAINT card definitions at https://bitbucket.org/pflotran/pflotran-dev/wiki/Documentation/QuickGuide).

`"reaction_network`" defines a file containing a PFLOTRAN CHEMISTRY block.

`"constraint`" defines a file containing a PFLOTRAN CONSTRAINT block.

.. code-block:: xml

  <geochemistry>
      <reaction_network file="calcite_flow_and_tran.in" format="simple"/>
      <constraint name="Initial" filename="calcite_flow_and_tran.in"/>
      <constraint name="Inlet" filename="calcite_flow_and_tran.in"/>
  </geochemistry>

Materials
=========

The "material" in this context is meant to represent the media through with fluid phases are transported. In the literature, this is also referred to as the "soil", "rock", "matrix", etc. Properties of the material must be specified over the entire simulation domain, and is carried out using the Region constructs defined above. For example, a single material may be defined over the "All" region (see above), or a set of materials can be defined over subsets of the domain via user-defined regions. If multiple regions are used for this purpose, they should be disjoint, but should collectively tile the entire domain. 

Material
--------

Within the Materials block an unbounded number of `"material`" elements can be defined.  Each material requires a label and has the following requirements.

.. code-block:: xml

  <material>
      Required Elements: mechanical_properties, permeability or hydraulic_conductivity, assigned_regions
      Optional Elements: comments, cap_pressure, rel_perm, sorption_isotherms 
  </material>
 
Mechanical_properties
---------------------

.. code-block:: xml

  <mechanical_properties>
      Required Elements: porosity, particle_density   (FILE OPTION NOT IMPLEMENTED) 
      Optional Elements: specific_storage, specific_yield, dispersion_tensor, tortuosity
  </mechanical_properties>

* `"mechanical_properties`" has six elements that can be either values or specified as files.  It has the following requirements.

    * `"porosity`" is defined in-line using attributes.  It is specified in oneof three ways: as a value between 0 and 1 using value="<value>", through a file using type="file" and filename="<filename>", or as a gslib file using type="gslib", parameter_file="<filename>", value="<value>" and (optionally) data_file="<filename>" (defaults to `"porosity_data`".  NOTE - FILE OPTION NOT IMPLEMENTED YET.

    * `"particle_density`" is defined in-line using attributes.  Either it is specified as a value greater than 0 using `"value`" or it specified through a file using `"filename`" and `"type`".  NOTE - FILE OPTION NOT IMPLEMENTED YET.

    * `"specific_storage`" is defined in-line using attributes.  Either it is specified as a value greater than 0 using `"value`" or it specified through a file using `"filename`" and `"type`".  NOTE - FILE OPTION NOT IMPLEMENTED YET.

    * `"specific_yield`" is defined in-line using attributes.  Either it is specified as a value using `"value`" or it specified through a file using `"filename`" and `"type`".  NOTE - FILE OPTION NOT IMPLEMENTED YET.

    * `"dispersion_tensor`" is defined in-line using attributes.  The attribute `"type`" is used to specify either the model to utilize of that a file is to be read.  The `"type`" options are: uniform_isotropic, burnett_frind, lichtner_kelkar_robinson, or file.  For `"uniform_isotropic`" values are specified using the attributes `"alpha_l`" and `"alpha_t`".  For `"burnett_frind`" values are specified using the attributes `"alpha_l`", `"alpha_th`", and `"alpha_tv`". For `"lichtner_kelkar_robinson`" values are specified using the attributes `"alpha_l`h", `"alpha_lv`", `"alpha_th`", and `"alpha_tv`".  For `"file`" the file name is specified using `"filename`".  NOTE - FILE OPTION NOT IMPLEMENTED YET.

    * `"tortuosity`" is defined in-line using attributes.  Either it is specified as a value using `"value`" or it specified through a file using `"filename`" and `"type`".  NOTE - FILE OPTION NOT IMPLEMENTED YET.


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

* `"assigned_regions`" is a comma separated list of region names for which this material is to be assigned.  Region names must be from the regions defined in the `"regions`" sections.  Region names can contain spaces.

.. code-block:: xml

    <assigned_regions>Region1, Region_2, Region 3</assigned_regions>

Permeability
------------

Permeability or hydraulic_conductivity must be specified but not both. If specified as constant values, permeability has the attributes `"x`", `"y`", and `"z`".  Permeability may also be extracted from the attributes of an Exodus II file, or generated as a gslib file.

.. code-block:: xml

  <permeability x="exponential" y="exponential" z="exponential" />
  or
  <permeability type="file" filename="file name" attribute="attribute name"/>
  or
  <permeability type="gslib" parameter_file="file name" value="exponential" data_file="file name"/>

Hydraulic_conductivity
----------------------

* `"hydraulic_conductivity`" is the hydraulic conductivity and has the attributes `"x`", `"y`", and `"z`". Permeability or hydraulic_conductivity must be specified but not both.

.. code-block:: xml

  <hydraulic_conductivity x="exponential" y="exponential" z="exponential" />
  or
  <hydraulic_conductivity type="gslib" parameter_file="file name" value="exponential" data_file="file name"/>

Cap_pressure
------------

*  `"cap_pressure`" is an optional element.  The available models are `"van_genuchten`", `"brooks_corey`", and `"none`".  The model name is specified in an attribute and parameters are specified in a subelement.  Model parameters are listed as attributes to the parameter element.

* `"van_genuchten`" parameters include `"alpha`", `"sr`", `"m`", and `"optional_krel_smoothing_interval`".  `"brooks_corey`" parameters include `"alpha`", `"sr`", `"m`", and `"optional_krel_smoothing_interval`".

.. code-block:: xml

  <cap_pressure model="van_genuchten | brooks_corey | none" >
      Required Elements: alpha, Sr, m (van_genuchten and brooks_corey only)
      Optional Elements: optional_krel_smoothing_interval (van_genuchten and brooks_corey only)
  </cap_pressure>

Rel_perm
--------

*  `"rel_perm`" is an optional element.  The available models are `"mualem`", `"burdine`", and `"none`".  The model name is specified in an attribute and parameters are specified in a subelement.  Model parameters are listed as attributes to the parameter element.

* `"mualem`" has no parameters.  `"burdine`" parameters include `"exp`".

.. code-block:: xml

  <rel_perm model="mualem | burdine | none )" >
      Required Elements: none 
      Optional Elements: exp (burdine only)
  </rel_perm>

Sorption_isotherms
------------------

*  `"sorption_isotherms`" is an optional element for providing Kd models and molecular diffusion values for individual solutes.  All solutes should be listed under each material.  Values of 0 indicate that the solute is not present/active in the current material.  The available Kd models are `"linear`", `"langmuir`", and `"freundlich`".  Different models and parameters are assigned per solute in sub-elements through attributes. The Kd and molecular diffusion parameters are specified in subelements.

.. code-block:: xml

    <sorption_isotherms>
	<solute name="string" />
            Required Elements: none
            Optional Elements: kd_model
    </sorption_isotherms>

.
    * `"kd_model`" takes the following form:

.. code-block:: xml
 
    <kd_model model="linear|langmuir|freundlich" kd="Value" b="Value (langmuir only)" n="Value (freundlich only)" />
  
    
Process Kernels
===============

.. code-block:: xml

  <process_kernels>
      Required Elements: flow, transport, chemistry
      Optional Elements: comments
  </process_kernels>

For each process kernel the element `"state`" indicates whether the solution is being calculated or not.  

Flow
----

* `"flow`" has the following attributes, 
      
      * `"state`" = "on | off"

      *  `"model`" = " richards | saturated | constant" 

Currently three scenarios are available for calculated the flow field.  `"richards`" is a single phase, variably saturated flow assuming constant gas pressure.  `"saturated`" is a single phase, fully saturated flow.  `"constant`" is equivalent to a flow model of single phase (saturated) with the time integration mode of transient with static flow in the version 1.2.1 input specification.  This flow model indicates that the flow field is static so no flow solver is called during time stepping. During initialization the flow field is set in one of two ways: (1) A constant Darcy velocity is specified in the initial condition; (2) Boundary conditions for the flow (e.g., pressure), along with the initial condition for the pressure field are used to solve for the Darcy velocity.

Note:  Unstructured options `"discretization_method`",  `"rel_perm_method`", and `"preconditioning_strategy`" have been moved to the `"unstr_flow_controls`" section under `"numerical_controls`"/

Transport
---------

* `"transport`" has the following attributes,
      
      * `"state`" = "on | off"

For `"transport`" the `"state`" must be specified.  

Note:  Unstructured options `"algorithm`" and `"sub_cycling`" have been moved to the `"unstr_transport_controls`" section under `"numerical_controls`"/

Chemistry
---------

* `"chemistry`" has the following attributes,
      
      * `"state`" = "on | off"
      
      * `"engine`" = "amanzi | pflotran | none"

      * `"process_model`" = "implicit operator split | none" 

For `"chemistry`" a combination of `"state`", `"engine`", and `"process_model`" must be specified.  If `"state`" is `"off`" then `"engine`" and `"process_model`" are set to `"none`".  Otherwise the `"engine`" and `"process_model`" model must be specified. 

Phases
======

Some general discussion of the `"Phases`" section goes here.

.. code-block:: xml

  <Phases>
      Required Elements: liquid_phase 
      Optional Elements: solid_phase
      Optional Elements: gas_phase [U]
  </Phases>

Liquid_phase
------------

* `"liquid_phase`" has the following elements

.. code-block:: xml

  <liquid_phase>
      Required Elements: viscosity, density
      Optional Elements: dissolved_components, eos [S]
  </liquid_phase>

Here is more info on the `"liquid_phase`" elements:

    * `"eos`" = "string" 

    * `"viscosity`" = "exponential"

    * `"density`" = "exponential"

    * `"dissolved_components`" has the required element

        * `"solutes`"

The subelement `"solutes`" can have an unbounded number of subelements `"solute`" which defines individual solutes present.  The `"solute`" element takes the following form:
  
    * `"solute`" = "string", containing the name of the solute

    * `"coefficient_of_diffusion`" = "exponential", this is an optional attribute

    * `"first_order_decay_constant`" = "exponential", this is an optional attribute

Solid_phase
-----------

* `"solid_phase`" has the following elements

.. code-block:: xml

  <solid_phase>
      Required Elements: minerals
      Optional Elements: NONE
  </solid_phase>

Here is more info on the `"solid_phase`" elements:

    * `"minerals`" has the element 

        * `"mineral`" which contains the name of the mineral

Initial Conditions
==================

Some general discussion of the `"initial_condition`" section goes here.

The `"initial_conditions`" section contains at least 1 and up to an unbounded number of `"initial_condition`" elements.  Each `"initial_condition`" element defines a single initial condition that is applied to one or more region.  The following is a description of the `"initial_condition`" element.

.. code-block:: xml

  <initial_condition>
      Required Elements: assigned_regions
      Optional Elements: liquid_phase (, comments, solid_phase - SKIPPED)
  </initial_condition>

Assigned_regions
----------------

* `"assigned_regions`" is a comma separated list of regions to apply the initial condition to.

Liquid_phase
------------

* `"liquid_phase`" has the following elements

.. code-block:: xml

  <liquid_phase>
      Required Elements: liquid_component
      Optional Elements: solute_component (, geochemistry  - SKIPPED)
  </liquid_phase>

*  Here is more info on the `"liquid_component`" elements:

    * `"uniform_pressure`" is defined in-line using attributes.  Uniform specifies that the initial condition is uniform in space.  Value specifies the value of the pressure.  
      
    * `"linear_pressure`" is defined in-line using attributes.  Linear specifies that the initial condition is linear in space.  Gradient specifies the gradient value in each direction in the form of a coordinate (grad_x, grad_y, grad_z).  Reference_coord specifies a reference location as a coordinate.  Value specifies the value of the pressure.
      
    * `"uniform_saturation`" is defined in-line using attributes.  See `"uniform_pressure`" for details.
      
    * `"linear_saturation`" is defined in-line using attributes. See `"linear_pressure`" for details.
      
    * `"velocity`" is defined in-line using attributes.  Specify the velocity is each direction using the appropriate attributes x, y, and z.

.. code-block:: xml

    <uniform_pressure name="some name" value="exponential" />
    <linear_pressure name="some name" value="exponential" reference_coord="coordinate" gradient="coordinate"/>
    <uniform_saturation name="some name" value="exponential" />
    <linear_saturation name="some name" value="exponential" reference_coord="coordinate" gradient="coordinate"/>
    <velocity name="some name" x="exponential" y="exponential" z="exponential"/>

*  Here is more info on the `"solute_component`" elements:

    * `"solute_component`" is defined in-line using attributes.  The attributes include "function", "value", and "name". Function specifies linear or constant temporal functional form during each time interval.  Value is the value of the `"solute_component`".  Name is the name of the solute component.

.. code-block:: xml

     <solute_component name="some name" value="exponential" function="uniform" />

..     <solute_component name="some name" (filename="filename" SKIPPED) value="exponential" function="uniform (|linear SKIPPED) " (reference_coord="coordinate" gradient="coordinate" - linear skipped) />

NOTE: Reading from a file is not yet implemented.  Also, the reference_coord and gradient attributes are only needed for the "linear" function type, which is also not yet implemented.

Geochemistry
------------

* `"geochemistry`" is an element with the following subelement: NOT IMPLEMENTED YET

   * `"constraint`" is an element with the following attributes: ONLY UNIFORM, for now

.. code-block:: xml

     <constraint name="some name" start="time" />

Solid_phase
-----------

* `"solid_phase`" has the following elements - Reminder this element has been SKIPPED

.. code-block:: xml

  <solid_phase>
      Required Elements: geochemistry - SKIPPED
      Optional Elements: mineral, geochemistry - BOTH SKIPPED 
  </solid_phase>

Here is more info on the `"solid_phase`" elements: - NOT IMPLEMENTED YET

    * `"mineral`" has the element - SKIPPED (EIB - I there's a typo in the schema here!)

        * `"mineral`" which contains the name of the mineral

    * `"geochemistry`" is an element with the following subelement: NOT IMPLEMENTED YET

        * `"constraint`" is an element with the following attributes: ONLY UNIFORM, for now

Boundary Conditions
===================

Some general discussion of the `"boundary_condition`" section goes here.

The `"boundary_conditions`" section contains at least 1 and up to an unbounded number of `"boundary_condition`" elements.  Each `"boundary_condition`" element defines a single initial condition that is applied to one or more region.  The following is a description of the `"boundary_condition`" element.

.. code-block:: xml

  <boundary_condition>
      Required Elements: assigned_regions, liquid_phase
      Optional Elements: comments - SKIPPED
  </boundary_condition>

Assigned_regions
----------------

* `"assigned_regions`" is a comma separated list of regions to apply the initial condition to.

Liquid_phase
------------

* `"liquid_phase`" has the following elements

.. code-block:: xml

  <liquid_phase>
      Required Elements: liquid_component
      Optional Elements: solute_component (, geochemistry - SKIPPED)
  </liquid_phase>

*  Here is more info on the `"liquid_component`" elements:

    * `"inward_mass_flux`" is defined in-line using attributes.  The attributes include "function", "start", and "value". Function specifies linear or constant temporal functional form during each time interval.  Start is a series of time values at which time intervals start.  Value is the value of the `"inward_mass_flux`" during the time interval. 

    * `"outward_mass_flux`" is defined in-line using attributes.  See `"inward_mass_flux`" for details.

    * `"inward_volumetric_flux`" is defined in-line using attributes.  See `"inward_mass_flux`" for details.

    * `"outward_volumetric_flux`" is defined in-line using attributes.  See `"inward_mass_flux`" for details.

    * `"uniform_pressure`" is defined in-line using attributes.  Uniform refers to uniform in spatial dimension.  See `"inward_mass_flux`" for details.

    * `"linear_pressure`" is defined in-line using attributes.  Linear refers to linear in spatial dimension. Gradient_value specifies the gradient value in each direction in the form of a coordinate (grad_x, grad_y, grad_z).  Reference_point specifies a reference location as a coordinate.  Reference_value specifies a reference value for the boundary condition. 

    * `"seepage_face`"is defined in-line using attributes.  The attributes include "function", "start", and "value". Function specifies linear or constant temporal functional form during each time interval.  Start is a series of time values at which time intervals start.  inward_mass_flux is the value of the inward_mass_flux during the time interval.
 
    * `"hydrostatic`" is an element with the attributes below.  By default the coordinate_system is set to "absolute".  Not specifying the attribute will result in the default value being used.  The attribute submodel is optional.  If not specified the submodel options will not be utilized.

    * `"linear_hydrostatic`" is defined in-line using attributes.  Linear refers to linear in spatial dimension. Gradient_value specifies the gradient value in each direction in the form of a coordinate (grad_x, grad_y, grad_z).  Reference_point specifies a reference location as a coordinate.  Reference_water_table_height specifies a reference value for the water table.  Optionally, the attribute "submodel" can be used to specify no flow above the water table height.

    * `"no_flow`" is defined in-line using attributes.  The attributes include "function" and "start". Function specifies linear or constant temporal functional form during each time interval.  Start is a series of time values at which time intervals start.  

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

Solute_component
----------------

*  Here is more info on the `"solute_component`" elements:

    * `"aqueous_conc`" is an element with the following attributes: ONLY CONSTANT, for now

.. code-block:: xml

     <aqueous_conc name="some name" value="exponential" function="linear | uniform | constant" start="time" />

*  Here is more info on the `"geochemistry`" elements:

    * `"constraint`" is an element with the following attributes: ONLY UNIFORM, for now

.. code-block:: xml

     <constraint name="some name" start="time" function="linear | uniform | constant"/>

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

* `"assigned_regions`" is a comma separated list of regions to apply the source to.

Liquid_phase
------------

* `"liquid_phase`" has the following elements

.. code-block:: xml

  <liquid_phase>
      Required Elements: liquid_component
      Optional Elements: solute_component (, geochemistry - SKIPPED)
  </liquid_phase>

*  Here is more info on the `"liquid_component`" elements:

    * `"volume_weighted`" is defined in-line using attributes.  The attributes include "function", "start", and "value". Function specifies linear or constant temporal functional form during each time interval.  Start is a series of time values at which time intervals start.  Value is the value of the `"volume_weighted`" during the time interval. 

    * `"perm_weighted`" is defined in-line using attributes.  See `"volume_weighted`" for details.

*  Here is more info on the `"solute_component`" elements:

    * `"uniform_conc`" is defined in-line using attributes.  The attributes include "name", "function", "start", and "value". Name is the name of a previously defined solute. Function specifies linear or constant temporal functional form during each time interval.  Start is a series of time values at which time intervals start.  Value is the value of the `"uniform_conc`" during the time interval. 

    * `"flow_weighted_conc`" is defined in-line using attributes.  See `"uniform_conc`" for details.

    * `"diffusion_dominated_release`" is defined in-line using attributes.  The attributes include "name", "start", "total_inventory", "mixing_length", and "effective_diffusion_coefficient". Name is the name of a previously defined solute. Start is a series of time values at which time intervals start.  Value is the value of the `"diffusion_dominated_release`" during the time interval. 

Output
======

Output data from Amanzi is currently organized into four specific elements: `"Vis`", `"Checkpoint`", `"Observations`", and `"Walkabout Data`".  Each of these is controlled in different ways, reflecting their intended use.

* `"Vis`" is intended to represent snapshots of the solution at defined instances during the simulation to be visualized.  The ''vis'' element defines the naming and frequencies of saving the visualization files.  The visualization files may include only a fraction of the state data, and may contain auxiliary "derived" information (see *elsewhere* for more discussion).

* `"Checkpoint`" is intended to represent all that is necessary to repeat or continue an Amanzi run.  The specific data contained in a Checkpoint Data dump is specific to the algorithm options and mesh framework selected.  Checkpoint is special in that no interpolation is performed prior to writing the data files; the raw binary state is necessary.  As a result, the user is allowed to only write Checkpoint at the discrete intervals of the simulation. The ''checkpoint'' element defines the naming and frequencies of saving the checkpoint files.

* `"Observations`" is intended to represent diagnostic values to be returned to the calling routine from Amanzi's simulation driver.  Observations are typically generated at arbitrary times, and frequently involve various point samplings and volumetric reductions that are interpolated in time to the desired instant.  Observations may involve derived quantities (see discussion below) or state fields.  The ''observations'' element may define one or more specific ''observation''.

* `"Walkabout Data`" is intended to be used as input to the particle tracking software Walkabout.

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

The ``vis`` element also includes an optional subelement ``write_regions``.  This was primarily implemented for debugging purposes but is also useful for visualizing fields only on specific regions.  The subelement accepts an arbitrary number of subelements named ``field``, with attibutes ``name`` (a string) and ``regions`` (a comma separated list of region names).  For each such subelement, a field will be created in the vis files using the name as a label.  The field will be initialized to 0, and then, for region list R1, R2, R3..., cells in R1 will be set to 1, cells in R2 will be set to 2, etc.  When regions in the list overlap, later ones in the list will take precedence.

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

The *base_filename* element contain the text component of the how the checkpoint files will be named.  The *base_filename* is appended with an index number to indicate the sequential order of the checkpoint files.  The *num_digits* elements indicates how many digits to use for the index. (*EIB NOTE* - verify if this is sequence index or iteration id)  Final the *cycle_macros* element indicates the previously defined cycle_macro to be used to determine the frequency at which to write the checkpoint files. Multiple cycle_macro may be specified in a comma seperated list. See the about NOTE about specifying a file location other than the current working directory.

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
                          solute_volumetric_flow_rate
     </liquid_phase>

The observation element identifies the field quantity to be observed.  Subelements identify the elements for a region, a model (functional) with which it will extract its source data, and a list of discrete times for its evaluation.  The observations are evaluated during the simulation and returned to the calling process through one of Amanzi arguments. The elements for each observation type are as follows:

.. code-block :: xml

   <observation_type>
     Required Elements: assigned_region, functional, time_macros or cycle_macros 
     Optional Elements: NONE
   </observation_type>

The only exceptions are aqueous_conc and solute_volumetric_flow_rate which both require a solute to be specified.  An additional subelement "solute" gives the name of the solute to calculate the aqueous concentration or volumetric flow rate for.  Be sure the name of given for the solute matches a defined solute elsewhere in the input file.  

NOTE: Previously individual observation elements had the subelement ''cycle_macro'' or ''time_macro''.  All output is moving away from only allowing a single macro to be specified to allowing multiple macros as a comma separated list.  To ease the transition for users both singular and plural are currently accepted.  However, the singular option will go away in the future.  Please update existing input files to use ''cycle_macros'' or ''time_macros''.


Example:

.. code-block :: xml

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

  <echo_translated_input format="some tag" file_name="some name"/>

* Write the input data after internal translation.  There are two specifiable attributes, `"format`" and `"file_name`".  If this parameter is omitted, no translated files are written.

  * `"format`" is a specific format tag, and can be `"v1`" (DEFAULT) or `"native`".  The actual format created for the `"native`" tag will depend on the value of the `"type`" specified under `"amanzi_input`" (see above).

  * `"file_name`" is the name of the translated output file.  If `"format`" = `"v1`", then `"file_name`" defaults to `"XXX_oldspec.xml`", where `"XXX.xml`" is the name of the original Amanzi input file.  If `"format`" = `"native`", then `"file_name`" defaults to `"translated_inpus.xml`".


Full Example
============

.. code-block:: xml

  <amanzi_input type="unstructured" version="2.1.0">
    <model_description name="example of full unstructured schema">
      <comments>comments here</comments>
      <model_id>XXX</model_id>
      <author>Erin Barker</author>
      <units>
        <length_unit>m</length_unit>
        <time_unit>s</time_unit>
        <mass_unit>kg</mass_unit>
        <conc_unit>molar</conc_unit>
      </units>
    </model_description>
    <echo_translated_input format="v1" file_name="my_translated_input.xml">
    <definitions>
      <macros>
        <time_macro name="time macro">
          <time>3.0e+10</time>
        </time_macro>
        <cycle_macro name="Every_20">
          <start>0</start>
          <timestep_interval>20</timestep_interval>
          <stop>-1</stop>
        </cycle_macro>
      </macros>
    </definitions>
    <process_kernels>
      <comments>Variably saturated flow</comments>
      <flow model="richards" state="on" discretization_method="fv-default" rel_perm_method="upwind-darcy_velocity"/>
      <transport algorithm="none" state="off" sub_cycling="off"/>
      <chemistry engine="none" process_model="none" state="off"/>
    </process_kernels>
    <phases>
      <liquid_phase name="water">
        <eos>false</eos>
        <viscosity>1.002E-03</viscosity>
        <density>998.2</density>
        <dissolved_components>
            <solutes>
                <solute coefficient_of_diffusion="1e-9" first_order_decay_constant="1.0">Tc-99</solute>
            </solutes>
        </dissolved_components>
      </liquid_phase>
      <solid_phase>
          <minerals>
              <mineral>Calcium</mineral>
          </minerals>
      </solid_phase>
    </phases>
    <execution_controls>
      <verbosity level="medium"/>
      <execution_control_defaults method="bdf1" mode="steady"/>
      <execution_control end="3.0e+10" init_dt="0.01" method="bdf1" mode="steady" reduction_factor="0.5" start="0.0"/>
    </execution_controls>
    <numerical_controls>
      <unstructured_controls>
        <unstr_linear_solver>
          <max_iterations>100</max_iterations>
          <tolerance>1.0e-17</tolerance>
          <method>gmres</method>
          <cfl>1</cfl>
          <preconditioner name="hypre_amg">
            <hypre_cycle_applications>5</hypre_cycle_applications>
            <hypre_smoother_sweeps>3</hypre_smoother_sweeps>
            <hypre_tolerance>0.0</hypre_tolerance>
            <hypre_strong_threshold>0.5</hypre_strong_threshold>
          </preconditioner>
        </unstr_linear_solver>
        <unstr_steady-state_controls>
          <initialize_with_darcy>true</initialize_with_darcy>
          <min_iterations>10</min_iterations>
          <max_iterations>15</max_iterations>
          <max_preconditioner_lag_iterations>5</max_preconditioner_lag_iterations>
          <nonlinear_tolerance>1.0e-5</nonlinear_tolerance>
          <limit_iterations>20</limit_iterations>
          <nonlinear_iteration_damping_factor>1</nonlinear_iteration_damping_factor>
          <nonlinear_iteration_divergence_factor>1000</nonlinear_iteration_divergence_factor>
          <max_divergent_iterations>3</max_divergent_iterations>
          <unstr_pseudo_time_integrator>
              <initialize_with_darcy>true</initialize_with_darcy>
              <clipping_saturation>0.9</clipping_saturation>
              <method>picard</method>
              <preconditioner>hypre_amg</preconditioner>
              <linear_solver>aztec00</linear_solver>
              <control_options>pressure</control_options>
              <convergence_tolerance>1.0e-8</convergence_tolerance>
              <max_iterations>100</max_iterations>
          </unstr_pseudo_time_integrator>
        </unstr_steady-state_controls>
      </unstructured_controls>
    </numerical_controls>
    <mesh framework="mstk">
      <comments>Two-dimensional box 499.872m x 73.152m</comments>
      <dimension>2</dimension>
      <generate>
        <number_of_cells nx="164" ny="120"/>
        <box high_coordinates="499.872, 73.152" low_coordinates="0.0, 0.0"/>
      </generate>
    </mesh>
    <regions>
      <comments/>
      <region name="Aquifer">
        <comments>One region comprising the entire domain</comments>
        <box high_coordinates="499.872, 73.152" low_coordinates="0.0, 0.0"/>
      </region>
      <region name="Left">
        <box high_coordinates="(0.0, 49.9872)" low_coordinates="(0.0, 0.0)"/>
      </region>
      <region name="Right">
        <box high_coordinates="(499.872, 73.152)" low_coordinates="(499.872, 0.0)"/>
      </region>
      <region name="Top">
        <box high_coordinates="(499.872, 73.152)" low_coordinates="(0.0, 73.152)"/>
      </region>
      <point coordinate="1.5240, 0.3048" name="Point5ft"/>
      <point coordinate="32.0040, 0.3048" name="Point105ft"/>
      <point coordinate="62.4840, 0.3048" name="Point205ft"/>
      <point coordinate="92.9640, 0.3048" name="Point305ft"/>
      <point coordinate="123.4440, 0.3048" name="Point405ft"/>
      <point coordinate="153.9240, 0.3048" name="Point505ft"/>
      <point coordinate="184.4040, 0.3048" name="Point605ft"/>
      <point coordinate="214.8840, 0.3048" name="Point705ft"/>
      <point coordinate="245.3640, 0.3048" name="Point805ft"/>
      <point coordinate="275.8440, 0.3048" name="Point905ft"/>
      <point coordinate="303.2760, 0.3048" name="Point1005ft"/>
      <point coordinate="336.8040, 0.3048" name="Point1105ft"/>
      <point coordinate="367.2840, 0.3048" name="Point1205ft"/>
      <point coordinate="397.7640, 0.3048" name="Point1305ft"/>
      <point coordinate="428.2440, 0.3048" name="Point1405ft"/>
      <point coordinate="458.7240, 0.3048" name="Point1505ft"/>
      <point coordinate="489.2040, 0.3048" name="Point1605ft"/>
      <point coordinate="498.3480, 0.3048" name="Point1635ft"/>
    </regions>
    <materials>
      <material name="Aquifer">
        <comments>Aquifer</comments>
        <mechanical_properties>
          <porosity value="0.43"/>
	  <particle_density value="2650.0"/>
        </mechanical_properties>
        <permeability x="1.1844e-12" y="1.1844e-12"/>
        <cap_pressure model="van_genuchten">
          <parameters alpha="1.46e-3" m="0.314" optional_krel_smoothing_interval="100.0" sr="0.052"/>
        </cap_pressure>
	<rel_perm model="mualem"/>
        <assigned_regions>Aquifer</assigned_regions>
        <sorption_isotherms>
            <solute name="Tc-99">
                <kd_model model="linear" kd="10.0"/>
            </solute>
        </sorption_isotherms>
      </material>
    </materials>
    <initial_conditions>
      <initial_condition name="Initial Condition">
        <comments>Aquifer</comments>
        <assigned_regions>Aquifer</assigned_regions>
        <liquid_phase name="water">
          <liquid_component name="water">
            <uniform_pressure value="101325.0"/>
          </liquid_component>
        </liquid_phase>
      </initial_condition>
    </initial_conditions>
    <boundary_conditions>
      <comments/>
      <boundary_condition name="LeftBC">
        <comments>Boundary condition at x=0</comments>
        <assigned_regions>Left</assigned_regions>
        <liquid_phase name="water">
          <liquid_component name="water">
            <hydrostatic function="constant" start="0.0" value="49.9872"/>
          </liquid_component>
        </liquid_phase>
      </boundary_condition>
      <boundary_condition name="TopBC">
        <comments>Boundary condition at y=73.152</comments>
        <assigned_regions>Top</assigned_regions>
        <liquid_phase name="water">
          <liquid_component name="water">
            <inward_mass_flux function="constant" start="0.0" value="1.1550e-4"/>
          </liquid_component>
        </liquid_phase>
      </boundary_condition>
    </boundary_conditions>
    <output>
       <vis>
        <base_filename>steady-flow</base_filename>
        <num_digits>5</num_digits>
        <time_macros>Steady State</time_macros>
      </vis>
    </output>
  </amanzi_input>

