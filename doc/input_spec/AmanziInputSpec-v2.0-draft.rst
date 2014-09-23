==============================================
Amanzi XML Input Specification (Version 2.0.x)
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

* `"Unstructured`": This instructs Amanzi to use data structures provided in the Trilinos software framework.  To the extent possible, the discretization algorithms implemented under this option are largely independent of the shape and connectivity of the underlying cells.  As a result, this option supports an arbitrarily complex computational mesh structure that enables users to work with numerical meshes that can be aligned with geometrically complex man-made or geostatigraphical features.  Under this option, the user typically provides a mesh file that was generated with an external software package.  The following mesh file formats are currently supported: `"Exodus II`".  Amanzi also provides a rudmentary capability to generate regular meshes within the unstructured framework internally.

An exmample root tag of an input file would look like the following.

.. code-block:: xml

  <amanzi_input version="2.0.0" type="unstructured"/>


Model Description
=================

This allows the users to provide a name and general description of model being devepoled.  This is also the section in which the units for the problem are stored. This entire section is optional but encouraged as documentation.

.. code-block:: xml

  <model_description name="Name of Model" >
      Required Elements: NONE
      Optional Elements: comment, author, created, modeified, model_id, description, purpose (units - NOT IMPLEMENTED YET)
  </model_description>

Units has the optional elements of length, time, mass, and concentration.  Each of those in turn have their own sturcture.  The structures are as follows.

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


Here is an overall example for the modle description element.

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

Definitions allows the user the define and name constants, times, and macros to be used in later sectons of the input file.  This is to streamline the look and readability of the input file.  The user should take care not to reuse names within this section or other sections.  This may have unindented consequences.

Named Times
-----------

Here the user can specify and name times to be used in other sections of the input file.   Note that if a name is repeated the last read value will be retained and all others will be overwritten.

A *time* requires the attributes `"name`" and `"value`".  If a unit is not specified with the value seconds is taken as the default.

Constants
---------

Here the user can define and name constants to be used in other sections of the input file.  Note that if a name is repeated the last read value will be retained and all others will be overwritten.

A *constant* has three attributes `"name`", `"type`", and `"value`".  The user can provide any name, but not it should not be repeated anywhere within the input to avoid confusion.  The available types include: `"none`", `"time`", `"constant`", `"numerical`", and `"area_mass_flux`".  Values assigned to constants of type `"time`" can include known units, otherwise seconds will be assumed as the default.

Macros
------

Three types of macros are currently available *time_macro*, *cycle_macro*, and *variable_macro*.

The *time_macro* requires an attribute `"name`".  The macro can then either take the form of one or more labeled time subelements or the subelements `"start`", `"timestep_interval`", and `"stop`" again containing labeled times.  A `"stop`" value of -1 will continue the cycle macro until the end of the simulation.  The labeled times can be time values assuming the default time unit of seconds or including a known time unit.

The *cycle_macro* requires an attribute `"name`" and the subelements `"start`", `"timestep_interval`", and `"stop`" with integer values.  A `"stop`" value of -1 will continue the cycle macro until the end of the simulation.

The *variable_macro* requires an attribute `"name`"  and one or more subelements `"variable`" containing strings.


An example *definitions* section would look as the following:

.. code-block:: xml

  <definitions>
    <constants>
      <constant name="zero" type="none" value="0.000"/>

      <constant name ="start"                   type="time" value="1956.0;y"/>
      <constant name ="B-18_release_end"        type="time" value ="1956.3288;y"/>
      <constant name="future_recharge"          type="area_mass_flux" value="1.48666E-6"/>

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


Execution Control
=================

Some general explaination of exection control goes here.

.. code-block:: xml
  
  <execution_controls>
      Required Elements: execution_control_defaults, execution_control
      Optional Elements: comments, verbosity
  </execution_controls>

Some explaination of each element goes here.

.. code-block:: xml
  
  <verbosity level="none | low | medium | high | extreme" />
 
Note, for debugging purposes use level="extreme". 

.. code-block:: xml

  <execution_control_defaults init_dt="labeled_time" max_dt="labeled_time" reduction_factor="exponential" increase_factor="exponential" mode="steady | transient" method=" bdf1 | picard" />

    * init_dt="labeled_time" 
      
    * max_dt="labeled_time" 
      
    * reduction_factor="exponential" 
      
    * increase_factor="exponential" 
      
    * mode="steady | transient" 
      
    * method=" bdf1 | picard" 

.. code-block:: xml

  <execution_control start="string" end="labeled_time" init_dt="labeled_time" max_dt="labeled_time" reduction_factor="exponential" increase_factor="exponential" mode="steady | transient" method=" bdf1 | picard" restart="string"/>

NOTE: start is REQUIRED
  
    * start="string", this attribute is required
      
    * end="labeled_time" 
      
    * init_dt="labeled_time" 
      
    * max_dt="labeled_time" 
      
    * reduction_factor="exponential" 
      
    * increase_factor="exponential" 
      
    * mode="steady | transient" 
      
    * method=" bdf1 | picard" 

    * restart="string"

Note, the value of the attribute ``restart`` is the name of the Amanzi checkpoint file previously created and to be used to initialize the current simulation.

Numerical Controls
==================

.. code-block:: xml

  <numerical_controls>
      Required Elements: NONE
      Optional Elements: comments, steady-state_controls, transient_controls, linear_solver, nonlinear_solver, chemistry_controls
  </numerical_controls>

Some discussion of the elements, what the minimum necessary for a simulation is goes here.  For now I have just listed the elements that are available.  

* `"comments`"="string" - SKIPPED 

    * Note: In many cases extra elements, such as comments, are not accommodated in the current input parsing. Therefore, for the most part `"comment`" elements are ignored.

* `"steady-state_controls`"  has the following elements

    * `"comments`"="string" - SKIPPED
 
    * `"min_iterations`"="integer" (min_iterations must be <= max_iterations)

    * `"max_iterations`"="integer"

    * `"max_preconditioner_lag_iterations`"="integer"

    * `"nonlinear_tolerance`"="exponential"

    * `"pseudo_time_integrator`"  has the following elements

        * `"method`"="string" (options: picard)

        * `"preconditioner`"="string" (options: trilinos_ml, hypre_amg, block_ilu) See below for subelements based on preconditioner name.

        * `"linear_solver`"="string" (options: aztec00)

        * `"control_options`"="string"

        * `"max_iterations`"="integer"

        * `"clipping_saturation`"="exponential"

        * `"convergence_tolerance`"="exponential"

        * `"initialize_with_darcy`"="boolean"

    * `"limit_iterations`"="integer"

    * `"nonlinear_iteration_damping_factor`"="exponential"

    * `"nonlinear_iteration_divergence_factor`"="exponential"

    * `"max_divergent_iterations`"="integer"

    * `"initialize_with_darcy`"="boolean"

* `"transient_controls`" has the elements `"comments`" and `"integration_method`". `"integration_method`" has the following elements

    * `"comments`"="string" - SKIPPED 
      
    * `"bdf1_integration_method`" has the following attributes

        * `"min_iterations`"="integer"

        * `"max_iterations`"="integer"

        * `"limit_iterations`"="integer"
 
        * `"nonlinear_tolerance`"="exponential"

        * `"max_preconditioner_lag_iterations`"="integer"

        * `"max_divergent_iterations`"="integer"

        * `"nonlinear_iteration_damping_factor`"="exponential"

        * `"nonlinear_iteration_divergence_factor`"="exponential"

        * `"restart_tolerance_factor`"="exponential"

        * `"restart_tolerance_relaxation_factor`"="exponential"

        * `"initialize_with_darcy`"="boolean"

    * `"preconditioner`" requires an attribute `"name`". (options: trilinos_ml, hypre_amg, block_ilu) See below for subelements based on preconditioner name.

* `"linear_solver`"  has the following elements

    * `"comments`"="string" - SKIPPED
 
    * `"method`"="string" (options: aztec00)

    * `"max_iterations`"="integer"

    * `"tolerance`"="exponential"

    * `"cfl`"="exponential"

    * `"preconditioner`" requires an attribute `"name`". (options: trilinos_ml, hypre_amg, block_ilu) See below for subelements based on preconditioner name.

* `"nonlinear_solver`"  has an attribute `"name`". (options: nka, newton, inexact newton)

* `"chemistry_controls`"  has the following elements

    * `"chem_tolerance`"="exponential" 
 
    * `"chem_max_newton_iterations`"="integer"

`"transient_controls`", `"linear_solver`", and `"pseudo_time_integrator`" accept a subelement for specifing the `"preconditioner`".  Current preconditioners available are Trilinos' ML, Hypre's AMG, and block ILU.  Below are the structures for each preconditioner.

* `"preconditioners`" with `"name = 'trilinos_ml'`" has the following optional elements

    * `"trilinos_smoother_type`"="string" (options: jacobi, gauss_seidel, ilu)

    * `"trilinos_threshold`"="exponential" 

    * `"trilinos_smoother_sweeps`"="integer"

    * `"trilinos_cycle_applications`"="integer"

* `"preconditioners`" with `"name = 'hypre_amg'`" has the following optional elements

    * `"hypre_cycle_applications`"="integer"

    * `"hypre_smoother_sweeps`"="integer"

    * `"hypre_tolerance`"="exponential" 

    * `"hypre_strong_threshold`"="exponential" 

* `"preconditioners`" with `"name = 'block_ilu'`" has the following optional elements

    * `"ilu_overlap`"="integer"

    * `"ilu_relax`"="exponential"

    * `"ilu_rel_threshold`"="exponential" 

    * `"ilu_abs_threshold`"="exponential" 

    * `"ilu_level_of_fill`"="integer" 


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

TODO: general description of what regions are

.. code-block:: xml

  <regions>
      Required Elements: NONE
      Optional Elements: comments, box, point, region
  </regions>

The regions block is required.  Within the region block no regions are required to be defined.  

The elements box and point allow for inline description of regions.  The region element uses a subelement to either define a box region or specify a region file.  

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

  <file name="filename" type=["color"|"labeled set"] format=["exodus ii"] entity=["cell"|"face"] label="integer"/>

Currently color functions and labeled sets can only be read from Exodus II files.  This will likely be the same file specified in the `"mesh`" element.  PLEASE NOTE the values listed within [] for attributes above are CASE SENSITIVE.  For many attributes within the Amanzi Input Schema the value is tested against a limited set of specific strings.  Therefore an user generated input file may generate errors due to a mismatch in cases.  Note that all specified names within this schema use lower case.

Geochemistry
============

Geochemistry allows users to define a reaction network and constraints to be associated with solutes defined under the `"dissolved_components`" section of the `"phases`" block. 

.. code-block:: xml

  <geochemistry>
      Required Elements: reaction_network, constraint
  </geochemistry>

PFLOTRAN Chemistry
------------------

For geochemisty simulated through PFLOTRAN, the user defines a reaction network and constraints.  These are defined within the same or separate text files through PFLOTRAN's input specification (see the CHEMISTRY and CONSTRAINT card definitions at https://bitbucket.org/pflotran/pflotran-dev/wiki/Documentation/QuickGuide).

`"reaction_network`" defines a file containing a PFLOTRAN CHEMISTRY block.

`"constraint`" defines a file containing a PFLOTRAN CONSTRAINT block.

.. code-block:: xml

  <geochemistry>
      <reaction_network file="calcite_flow_and_tran.in" format="simple"/>
      <constraint name="Initial" filename="calcite_flow_and_tran.in"/>
      <constraint name="Inlet" filename="calcite_flow_and_tran.in"/>
  </geochemistry>

Material
========

TODO - general description of the material section

Within the Materials block an unbounded number of `"material`" elements can be defined.  Each material has the following requirements.

.. code-block:: xml

  <material>
      Required Elements: mechanical_properties, permeability, hydraulic_conductivity, assigned_regions
      Optional Elements: comments, cap_pressure, rel_perm 
  </material>

`"mechanical_properties`" has two elements that can be either values or specified as files.  It has the following requirements.

.. code-block:: xml

  <mechanical_properties>
      Required Elements: porosity, particle_density   (FILE OPTION NOT IMPLEMENTED) 
  </mechanical_properties>

* `"porosity`" is defined inline using attributes.  Either it is specified as a value between 0 and 1 using `"value`" or it specified through a file using `"filename`" and `"type`". NOTE - FILE OPTION NOT IMPLEMENTED YET.

.. code-block:: xml

  <porosity value="decimal value"/>
  <porosity filename="file name" type="file"/>

* `"particle_density`" is defined inline using attributes.  Either it is specified as a value greater than 0 using `"value`" or it specified through a file using `"filename`" and `"type`".  See porosity for example.  NOTE - FILE OPTION NOT IMPLEMENTED YET.

* `"assigned_regions`" is a comma seperated list of region names for which this material is to be assigned.

* `"permeability`" is the permiability and has the attributes `"x`", `"y`", and `"z`".

.. code-block:: xml

  <permeability x="exponential" y="exponential" z="exponential" />

* `"hydraulic_conductivity`" is the hydraulic conductivity and has the attributes `"x`", `"y`", and `"z`".

.. code-block:: xml

  hydraulic_conductivity x="exponential" y="exponential" z="exponential" />

* `"cap_pressure`" is an optional element.  The available models are `"van_genuchten`", `"brooks_corey`", and `"none`".  The model name is specified in an attribute and parameters are specified in a subelement.  Model parameters are listed as attributes to the parameter element.

  * `"van_genuchten`" parameters include `"alpha`", `"sr`", `"m`", and `"optional_krel_smoothing_interval`".  `"brooks_corey`" parameters include `"alpha`", `"sr`", `"m`", and `"optional_krel_smoothing_interval`".

.. code-block:: xml

  <cap_pressure name="van_genuchten | brooks_corey | none )" >
      Required Elements: parameters
  </cap_pressure>

* `"rel_perm`" is an optional element.  The available models are `"mualem`", `"burdine`", and `"none`".  The model name is specified in an attribute and parameters are specified in a subelement.  Model parameters are listed as attributes to the parameter element.

  * `"mualem`" has no parameters.  `"burdine`" parameters include `"exp`".

.. code-block:: xml

  <rel_perm name="mualem | burdine | none )" >
      Required Elements: none 
      Optional Elements: exp (burdine only)
  </rel_perm>

* `"<sorption_isotherms>`" is an optional element for providing Kd models and molecular diffusion values for individual solutes.  All solutes should be listed under each material.  Values of 0 indicate that the solute is not present/active in the current material.  The available Kd models are `"linear`", `"langmuir`", and `"freundlich`".  Different models and parameters are assigned per solute in sub-elements through attributes. The Kd and molecular diffusion parameters are specified in subelements.

.. code-block:: xml

    <sorption_isotherms>
	<solute name="string" />
        model="linear | langmuir | langmuir" kd="exponential" b="exponential" n="exponential"/>
            Required Elements: none
            Optional Elements: kd_model, molecular_diffusion
    </sorption_isotherms>

The subelements kd_model and molecular_diffusion that the following forms:

.. code-block:: xml
 
    <kd_model model="linear|langmuir|freundlich" kd="Value" b="Value (langmuir only)" n="Value (freundlich only)" />
  
    
.. code-block:: xml
   
    <molecular_diffusion value="Value" />
    or
    <molecular_diffusion type="exodus ii" filename="file" />


Process Kernels
===============

.. code-block:: xml

  <process_kernels>
      Required Elements: flow, transport, chemistry
      Optional Elements: comments
  </process_kernels>

For each process kernel the element `"state`" indicates whether the solution is being calculated or not.  

* `"flow`" has two attributes, `"state`" and `"model`".
      
      * `"state`" = "on | off"

      *  `"model`" = " richards | saturated | constant" 

Currently three scenerios are avaiable for calculated the flow field.  `"richards`" is a single phase, variably saturated flow assuming constant gas pressure.  `"saturated`" is a single phase, fully saturated flow.  `"constant`" is equivalent to the a flow model of single phase (saturated) with the time integration mode of transient with static flow in the version 1.2.1 input specification.  This flow model indicates that the flow field is static so no flow solver is called during time stepping. During initialization the flow field is set in one of two ways: (1) A constant Darcy velocity is specified in the initial condition; (2) Boundary conditions for the flow (e.g., pressure), along with the initial condition for the pressure field are used to solve for the Darcy velocity.

* `"transport`" has two attributes, `"state`" and `"algorithm`".
      
      * `"state`" = "on | off"

      *  `"algorithm`" = " explicit first-order | explicit second-order | none " 

      * `"sub_cycling`" = "on | off"

For `"transport`" a combination of `"state`" and `"algorithm`" must be specified.  If `"state`" is `"off`" then `"algorithm`" is set to `"none`".  Otherwise the integration algorithm must be specified.  Whether sub-cycling is to be utilized within the transport algorithm is also specified here.

* `"chemistry`" has three attributes, `"state`", `"engine`", and `"process_model`".
      
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
      Optional Elements: solid_phase (comments - skipped)
  </Phases>

* `"liquid_phase`" has the following elements

.. code-block:: xml

  <liquid_phase>
      Required Elements: viscosity, density
      Optional Elements: dissolved_components, eos
  </liquid_phase>

Here is more info on the `"liquid_phase`" elements:

    * `"eos`"="string" 

    * `"viscosity`"="exponential"

    * `"density`"="exponential"

    * `"dissolved_components`" has the required element

        * `"solutes`"

The subelement `"solutes`" can have an unbounded number of subelements `"solute`" which defines individual solutes present.  The `"solute`" element takes the following form:
  
    * `"solute`"="string", containing the name of the solute

        * `"coefficient_of_diffusion`"="exponential", this is an optional attribute

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

* `"assigned_regions`" is a comma seperated list of regions to apply the initical condition to.

* `"liquid_phase`" has the following elements

.. code-block:: xml

  <liquid_phase>
      Required Elements: liquid_component
      Optional Elements: solute_component (, geochemistry  - SKIPPED)
  </liquid_phase>

Here is more info on the `"liquid_phase`" elements:

    * `"liquid_component`" is an element with the following subelement: 

        * `"uniform_pressure`" is an element with the following attributes: 

.. code-block:: xml

     <uniform_pressure name="some name" value="exponential" />

        * `"linear_pressure`" is an element with the following attributes: 

.. code-block:: xml

     linear_pressure name="some name" value="exponential" reference_coord="coordinate" gradient="coordinate"/>

        * `"velocity`" is an element with the following attributes: 

.. code-block:: xml

     <velocity name="some name" x="exponential" y="exponential" z="exponential"/>

    * `"solute_component`" is an element with the following attributes: 

.. code-block:: xml

     <solute_component name="some name" (filename="filename" SKIPPED) value="exponential" function="uniform (|linear SKIPPED) " (reference_coord="coordinate" gradient="coordinate" - linear skipped) />

NOTE: Reading from a file is not yet implemeneted.  Also, the reference_coord and gradient attributes are only needed for the "linear" function type, which is also not yet implemeneted.

    * `"geochemistry`" is an element with the following subelement: NOT IMPLEMENTED YET

        * `"constraint`" is an element with the following attributes: ONLY UNIFORM, for now

.. code-block:: xml

     <constraint name="some name" start="time" />

* `"solid_phase`" has the following elements - Remineder this element has been SKIPPED

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

* `"assigned_regions`" is a comma seperated list of regions to apply the initical condition to.

* `"liquid_phase`" has the following elements

.. code-block:: xml

  <liquid_phase>
      Required Elements: liquid_component
      Optional Elements: solute_component (, geochemistry - SKIPPED)
  </liquid_phase>

Here is more info on the `"liquid_phase`" elements:

    * `"liquid_component`" is an element with the following subelement: 

        * `"inward_mass_flux`" is an element with the following attributes: 

.. code-block:: xml

     <inward_mass_flux value="exponential" function="linear | uniform | constant" start="time" />

.
        * `"inward_volumetric_flux`" is an element with the following attributes: 

.. code-block:: xml

     <inward_volumetric_flux value="exponential" function="linear | uniform | constant" start="time" />

.
        * `"uniform_pressure`" is an element with the following attributes: 

.. code-block:: xml

     <uniform_pressure name="some name" value="exponential" function="uniform | constant" start="time" />

.
        * `"seepage_face`" is an element with the following attributes: 

.. code-block:: xml

     <seepage_face name="some name" inward_mass_flux="exponential" function="linear | uniform | constant" start="time" />

.
        * `"hydrostatic`" is an element with the following attributes: ONLY CONSTANT, for now

.. code-block:: xml

     <hydrostatic name="some name" value="exponential" function="uniform | constant" start="time" coordinate_system="absolute | relative to mesh top"/>

.
    * `"solute_component`" is an element with the following subelement: 

        * `"aqueous_conc`" is an element with the following attributes: ONLY CONTANT, for now

.. code-block:: xml

     <aqueous_conc name="some name" value="exponential" function="linear | uniform | constant" start="time" />

.
    * `"geochemistry`" is an element with the following subelement: NOT IMPLEMENTED YET

        * `"constraint`" is an element with the following attributes: ONLY UNIFORM, for now

.. code-block:: xml

     <constraint name="some name" start="time" function="linear | uniform | constant"/>

Output
======

Output data from Amanzi is currently organized into three specific elements: `"Vis`", `"Checkpoint`", `"Observations`", and `"Walkabout Data`".  Each of these is controlled in different ways, reflecting their intended use.

* `"Vis`" is intended to represent snapshots of the solution at defined instances during the simulation to be visualized.  The ''vis'' element defines the naming and frequencing of saving the visualization files.  The visualizatoin files may include only a fraction of the state data, and may contiain auxiliary "derived" information (see *elsewhere* for more discussion).

* `"Checkpoint`" is intended to represent all that is necesary to repeat or continue an Amanzi run.  The specific data contained in a Checkpoint Data dump is specific to the algorithm options and mesh framework selected.  Checkpoint is special in that no interpolation is perfomed prior to writing the data files; the raw binary state is necessary.  As a result, the user is allowed to only write Checkpoint at the discrete intervals of the simulation. The ''checkpoint'' element defines the naming and frequencing of saving the checkpoint files.

* `"Observations`" is intended to represent diagnostic values to be returned to the calling routine from Amanzi's simulation driver.  Observations are typically generated at arbitrary times, and frequently involve various point samplings and volumetric reductions that are interpolated in time to the desired instant.  Observations may involve derived quantities (see discussion below) or state fields.  The ''observations'' element may define one or more specific ''observation''.

* `"Walkabout Data`" is intended to be used as input to the particle tracking software Walkabout.

*EIB NOTE* - All three of the above are REQUIRED!!
For the obserservations I understand how to leave that empty.  But how do I execute without writing a checkpoint? If I'm running a dinky test am I really required to specify a checkpoint?  Will need to test this will validator.  Talk to Ellen about this.

Vis
---

The ''vis'' element defines the visualization filenaming scheme and how often to write out the files.  Thus, the ''vis'' element has the following requiements

.. code-block:: xml

  <vis>
      Required Elements: base_filename, num_digits 
      Optional Elements: time_macros, cycle_macros
  </vis>

The *base_filename* element contain the text component of the how the visualization files will be named.  The *base_filename* is appended with an index number to indicate the seqential order of the visualization files.  The *num_digits* elements indicates how many digits to use for the index. 

The presence of the ''vis'' element means that visualization files will be written out after cycle 0 and the final cycle of the simulation.  The optional elements *time_macros* or *cycle_macros* indicate additional points during the simulation at which visualization files are to be written out.  Both elements allow one or more of the appropriate type of macro to be listed.  These macros will be determine the appropriate times or cycles to write out visualization files.  See the `Definitions`_ section for defining individual macros.

(*EIB NOTE* - there should be a comment here about how the output is controlled, i.e. for each PK where do you go to turn on and off fields.  This will probably get filled in as the other sections fill out.)

Example:

.. code-block:: xml

  <vis>
     <base_filename>plot</base_filename>
     <num_digits>5</num_digits>
     <time_macros>Macro 1</time_macros>
  </vis>


Checkpoint
----------

The ''checkpoint'' element defines the filenaming scheme and frequency for writing out the checkpoint files.  As mentioned above, the user does not influence what is written to the checkpoint files.  Thus, the ''checkpoint'' element has the following requiements

.. code-block:: xml

  <checkpoint>
      Required Elements: base_filename, num_digits, cycle_macro
      Optional Elements: NONE
  </checkpoint>

The *base_filename* element contain the text component of the how the checkpoint files will be named.  The *base_filename* is appended with an index number to indicate the seqential order of the checkpoint files.  The *num_digits* elements indicates how many digits to use for the index. (*EIB NOTE* - verify if this is sequence index or interation id)  Final the *cycle_macro* element indicates the previously defined cycle_macro to be used to determine the frequency at which to write the checkpoint files.

Example:

.. code-block:: xml

  <checkpoint>
     <base_filename>chk</base_filename>
     <num_digits>5</num_digits>
     <cycle_macro>Every_100_steps</cycle_macro>
  </checkpoint>


Observations
------------

The Observations element holds all the observations that the user is requesting from Amanzi, as well as meta data, such as the name of the file that Amanzi will write observations to.  The observations are collected by their phase. Thus, the ''observations'' element has the following requirements

.. code-block:: xml

   <observations>
     Required Elements: filename, liquid_phase
     Optional Elements: NONE
   </observations>

The *filename* element contains the filename for the observation output, and may include the full path.  Currently, all observations are written to the same file.  

The *liquid_phase* element requires that the name of the phase be specified as an attribute and at least one observaton.  The observation element is named according to what is being observed.  The observations elements available are as follows:

.. code-block:: xml

     <liquid_phase name="Name of Phase (Required)">
       Required Elements: NONE 
       Optional Elements: integrated_mass, volumetric_water_content, gravimetric_water_content, aqueous_pressure, 
                          x_aqueous_volumetric_flux, y_aqueous_volumetric_flux, z_aqueous_volumetric_flux, material_id, 
                          hydraulic_head, aqueous_mass_flow_rate, aqueous_volumetric_flow_rate, aqueous_conc, drawdown
     </liquid_phase>

The observation element identifies the field quantity to be observed.  Subelements identify the elements for a region, a model (functional) with which it will extract its source data, and a list of discrete times for its evaluation.  The observations are evaluated during the simulation and returned to the calling process through one of Amanzi arguments. The elements for each observation type are as follows:

.. code-block :: xml

   <observation_type>
     Required Elements: assigned_region, functional, time_macro
     Optional Elements: NONE
   </observation_type>

The only exception is aqueous_conc requires an attribute Name="Solute Name".

Example:

.. code-block :: xml

    <observations>

      <filename>observation.out</filename>

      <liquid_phase name="water">
	<aqueous_pressure>
	  <assigned_regions>Obs_r1</assigned_regions>
	  <functional>point</functional>
	  <time_macro>Observation Times</time_macro>
	</aqueous_pressure>
	<aqueous_pressure>
	  <assigned_regions>Obs_r2</assigned_regions>
	  <functional>point</functional>
	  <time_macro>Observation Times</time_macro>
	</aqueous_pressure>
	<aqueous_pressure>
	  <assigned_regions>Obs_r2</assigned_regions>
	  <functional>point</functional>
	  <time_macro>Observation Times</time_macro>
	</aqueous_pressure>
      </liquid_phase>

    </observations>

Walkabout
----------

The ''walkabout'' element deines the filenaming scheme and frequency for writing out the walkabout files.  As mentioned above, the user does not influence what is written to the walkabout files only the writing frequency and naming scheme.  Thus, the ''walkabout'' element has the following requiements

.. code-block:: xml

  <walkabout>
      Required Elements: base_filename, num_digits, cycle_macro
      Optional Elements: NONE
  </walkabout>

The *base_filename* element contain the text component of the how the walkabout files will be named.  The *base_filename* is appended with an index number to indicate the seqential order of the walkabout files.  The *num_digits* elements indicates how many digits to use for the index.  Final the *cycle_macro* element indicates the previously defined cycle_macro to be used to determine the frequency at which to write the walkabout files.

Example:

.. code-block:: xml

  <walkabout>
     <base_filename>chk</base_filename>
     <num_digits>5</num_digits>
     <cycle_macro>Every_100_steps</cycle_macro>
  </walkabout>



