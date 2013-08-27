==============================================
Amanzi XML Input Specification (Version 2.0.x)
==============================================

.. contents:: **Table of Contents**

Overview
========

The Amanzi simulator evolves a system of conservation equations for
reacting flows in porous media, as detailed in the ASCEM report
entitled "Mathematical Formulation Requirements and Specifications for
the Process Models`" (hereafter referred to as the 'Model Requirements
Document (MRD)'). The purpose of the present document is to specify
the data required to execute Amanzi.  This specification should be
regarded as a companion to the MRD, and parameterizations of the
individual submodels are consistent between Amanzi, the MRD and this
document. Where applicable, the relevant sections of the MRD are
indicated.


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
  
  <verbosity level="low | medium | high" />
  
QUESTION: EIB - I don't understand what `"execution_control_defaults`" gets used for verses `"execution_control`"?  For now I am just skipping `"execution_control_defaults`".

.. code-block:: xml

  <execution_control_defaults init_dt="labeled_time" max_dt="labeled_time" reduction_factor="exponential" increase_factor="exponential" mode="stead | transient" method=" bdf1 | picard" />

NOTE: EIB - I don't understand how the method maps back to the old spec: bdf1 | picard? Is bdf1 the default? Does picard means, "Use Picard = true"?

    * init_dt="labeled_time" 
      
    * max_dt="labeled_time" 
      
    * reduction_factor="exponential" 
      
    * increase_factor="exponential" 
      
    * mode="stead | transient" 
      
    * method=" bdf1 | picard" 

.. code-block:: xml

  <execution_control start="string" end="labeled_time" init_dt="labeled_time" max_dt="labeled_time" reduction_factor="exponential" increase_factor="exponential" mode="stead | transient" method=" bdf1 | picard" />

NOTE: start is REQUIRED
  
    * start="string", this attribute is required
      
    * end="labeled_time" 
      
    * init_dt="labeled_time" 
      
    * max_dt="labeled_time" 
      
    * reduction_factor="exponential" 
      
    * increase_factor="exponential" 
      
    * mode="stead | transient" 
      
    * method=" bdf1 | picard" 

SKIPPED ATTRIBUTES: max_dt, reduction_factor, increase_factor

Numerical Controls
==================

.. code-block:: xml

  <numerical_controls>
      Required Elements: NONE????
      Optional Elements: steady-state_controls, transient_controls, comments, linear_solver (not specified)
  </numerical_controls>

NOTE: EIB - Currently `"linear_solver`" isn't listed in the schema with a min/max occurs.

Some discussion of the elements, what the minimum necessary for a simulation is goes here.  For now I have just listed the elements that are available.  

* `"comments`"="string" - SKIPPED 

    * Note: In many cases extra elements, such as comments, are not accommodated in the current input parsing. Therefore, for the most part `"comment`" elements are ignored.

* `"steady-state_controls`"  has the following elements

    * `"comments`"="string" - SKIPPED
 
    * `"min_iterations`"="integer"

    * `"max_iterations`"="integer"

    * `"max_preconditioner_lag_iterations`"="integer"

    * `"nonlinear_tolerance`"="exponential"

    * `"error_rel_tol`"="exponential"

    * `"error_abs_tol`"="exponential"

    * `"pseudo_time_integrator`"  has the following elements

        * `"method`"="string"

        * `"preconditioner`"="string"

        * `"linear_solver`"="string"

        * `"control_options`"="string"

        * `"divergent_max_iterations`"="integer"

        * `"clipping_saturation`"="exponential"

        * `"convergence_tolerance`"="exponential"

        * `"initialize_with_darcy`"="string"

* `"transient_controls`" has the elements `"comments`" and `"integration_method`". `"integration_method`" has the following elements

    * `"comments`"="string" - SKIPPED 
      
    * `"integration_method`" has the following elements

        * `"convergence_criteria`" has the following elements

            * `"error_rel_tol`"="exponential"

            * `"error_abs_tol`"="exponential"

        * `"nonlinear_solver_parameters`" has the following elements

            * `"min_iterations`"="integer"

            * `"max_iterations`"="integer"

            * `"limit_iterations`"="integer"
 
            * `"nonlinear_tolerance`"="exponential"

            * `"max_divergent_iterations`"="integer"

            * `"max_preconditioner_lag`"="integer"

* `"linear_solver`"  has the following elements

    * `"comments`"="string" - SKIPPED
 
    * `"method`"="string"

    * `"max_iterations`"="integer"

    * `"tolerance`"="exponential"

    * `"ml_cycle_applications`"="integer"

    * `"use_hypre_amg`"="string"

    * `"use_block_ilu`"="string"

    * `"hypre_amg_cycle_applications`"="integer"

    * `"hypre_amg_smoother_sweeps`"="integer"

    * `"hypre_amg_tolerance`"="exponential"

    * `"hypre_amg_threshold`"="exponential"

    * `"ml_smoother_type`"="string"

    * `"sub_cycling`"="string"

    * `"transport_sub_cycling`"="string"



Mesh
====

Amanzi supports both structured and unstructured numerical solution approaches.  This flexibility has a direct impact on the selection and design of the underlying numerical algorithms, the style of the software implementations, and, ultimately, the complexity of the user-interface.  "Mesh`" is used to select between the following options:

* `"Structured`": This instructs Amanzi to use BoxLib data structures and an associated paradigm to numerically represent the flow equations.  Data containers in the BoxLib software library, developed by CCSE at LBNL, are based on a hierarchical set of uniform Cartesian grid patches.  `"Structured`" requires that the simulation domain be a single coordinate-aligned rectangle, and that the "base mesh" consists of a logically rectangular set of uniform hexahedral cells.  This option supports a block-structured approach to dynamic mesh refinement, wherein successively refined subregions of the solution are constructed dynamically to track "interesting" features of the evolving solution.  The numerical solution approach implemented under the `"Structured`" framework is highly optimized to exploit regular data and access patterns on massively parallel computing architectures.

* `"Unstructured`": This instructs Amanzi to use data structures provided in the Trilinos software framework.  To the extent possible, the discretization algorithms implemented under this option are largely independent of the shape and connectivity of the underlying cells.  As a result, this option supports an arbitrarily complex computational mesh structure that enables users to work with numerical meshes that can be aligned with geometrically complex man-made or geostatigraphical features.  Under this option, the user typically provides a mesh file that was generated with an external software package.  The following mesh file formats are currently supported: `"Exodus 2`" (see example), `"MSTK`" (see example), `"MOAB`" (see example).  Amanzi also provides a rudmentary capability to generate unstructured meshes automatically.

.. code-block:: xml

   <mesh class=unstructured framework=["mstk"|"stk::mesh"|"moab"|"simple"]>

      <comments> May be included in the Mesh element </comments>

      <generate>
         <number_of_cells nx = "integer value"  ny = "integer value"  nz = "integer value"/>
         <box  low_coordinates = "x_low,y_low,z_low" high_coordinates = "x_high,y_high,z_high"/>
      </generate>

   </mesh>

testing.

.. code-block:: xml

  <mesh framework="mstk"> <!-- default is MSTK for unstructured -->
   <dimension>3</dimension>
   <generate>
     <number_of_cells nx = "64"  ny = "56"  nz = "107"/>
     <box  low_coordinates = "0.0,0.0,0.0" high_coordinates = "320.0,280.0,107.0"/>
   </generate>
  </mesh>


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
      Required Elements: region  ( OR file - NOT IMPLEMENTED YET)
      Optional Elements: comments
  </region>

A region is define as describe above.  A file is define as follows.

REMINDER - FILE OPTION NOT YET IMPLEMENTED

.. code-block:: xml

  <file name="file name" type="color | labeled set" format="exodus ii" entity="cell | face" label="integer"/>

Some discussion of reading a region file goes here. Talk about the color function or labeled set.  State we only read the ExodusII mesh format files.  State the region file must be specify cells or faces.  Explain what the label is for.

Geochemistry
============

Geochemistry allows users to define a reaction network and constraints to be associated with solutes defined under the dissolved_components section of the phases block. 

PFLOTRAN Chemistry
------------------

For geochemisty simulated through PFLOTRAN, the user defines a reaction network and constraints.  These are defined within the same or separate text files through PFLOTRAN's input specification (see the CHEMISTRY and CONSTRAINT card definitions at https://bitbucket.org/pflotran/pflotran-dev/wiki/Documentation/QuickGuide).

.. code-block::xml

  <geochemistry>
      Required elements: reaction_network, constraint
  </geochemistry>

`"reaction_network`" defines a file containing a PFLOTRAN CHEMISTRY block.

`"constraint`" defines a file containing a PFLOTRAN CONSTRAINT block.

Material
========

TODO - general description of the material section

Within the Materials block an unbounded number of `"material`" elements can be defined.  Each material has the following requirements.

.. code-block:: xml

  <material>
      Required Elements: mechanical_properties, permeability, assigned_regions
      Optional Elements: comments, cap_pressure (rel_perm - NOT YET IMPLEMENTED)
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

* `"cap_pressure`" is an optional element.  The available models are `"van_genuchten`", `"brooks_corey`", and `"none`".  The model name is specified in an attribute and parameters are specified in a subelement.  Model parameters are listed as attributes to the parameter element.

  * `"van_genuchten`" parameters include `"alpha`", `"sr`", and `"m`".  `"brooks_corey`" parameters include `"alpha`", `"sr`", and `"m`".

.. code-block:: xml

  <cap_pressure name="van_genuchten ( NOT IMPLEMENTED YET - | brooks_corey | none )" >
      Required Elements: parameters
  </cap_pressure>


REMINDER - REL_PERM IS NOT YET IMPLEMENTED

* `"rel_perm`" is an optional element.  The available models are `"mualem`", `"burdine`", and `"none`".  The model name is specified in an attribute and parameters are specified in a subelement.  Model parameters are listed as attributes to the parameter element.

  * `"mualem`" parameters include `"optional_krel_smoothing_interval`".  `"burdine`" parameters include `"optional_krel_smoothing_interval`", and `"exp`".

.. code-block:: xml

  <rel_perm name="mualem | burdine | none )" >
      Required Elements: parameters
  </rel_perm>



Process Kernels
===============

.. code-block:: xml

  <process_kernels>
      Required Elements: flow, transport, chemistry
      Optional Elements: comments
  </process_kernels>

* `"flow`" has two attributes, `"state`" and `"model`".
      
      * `"state`" = "on | off"

      *  `"model`" = " richards | saturated " 

* `"transport`" has two attributes, `"state`" and `"algorithm`".
      
      * SKIPPED FOR NOW

* `"chemistry`" has two attributes, `"state`" and `"process_model`".
      
      * SKIPPED FOR NOW

Phases
======

Some general discussion of the `"Phases`" section goes here.

.. code-block:: xml

  <Phases>
      Required Elements: liquid_phase
      Optional Elements: comments, solid_phase - SKIPPED
  </Phases>

* `"liquid_phase`" has the following elements

.. code-block:: xml

  <liquid_phase>
      Required Elements: eos, viscosity, density
      Optional Elements: dissolved_components - SKIPPED
  </liquid_phase>

Here is more info on the `"liquid_phase`" elements:

    * `"eos`"="string" - QUESTION - EIB: note sure what this translate to in the old spec.

    * `"viscosity`"="exponential"

    * `"density`"="exponential"

    * `"dissolved_components`" has the elements - SKIPPED

        * `"solutes`"

* `"solid_phase`" has the following elements - Remineder this element has been SKIPPED

.. code-block:: xml

  <solid_phase>
      Required Elements: minerals
      Optional Elements: NONE
  </solid_phase>

Here is more info on the `"solid_phase`" elements:

    * `"minerals`" has the element - SKIPPED

        * `"mineral`" which contains the name of the mineral

Initial Conditions
==================

Some general discussion of the `"initial_condition`" section goes here.

The `"initial_conditions`" section contains at least 1 and up to an unbounded number of `"initial_condition`" elements.  Each `"initial_condition`" element defines a single initial condition that is applied to one or more region.  The following is a description of the `"initial_condition`" element.

.. code-block:: xml

  <initial_condition>
      Required Elements: assigned_regions, liquid_phase
      Optional Elements: comments, solid_phase - SKIPPED
  </initial_condition>

* `"assigned_regions`" is a comma seperated list of regions to apply the initical condition to.

* `"liquid_phase`" has the following elements

.. code-block:: xml

  <liquid_phase>
      Required Elements: liquid_component
      Optional Elements: solute_component, geochemistry - BOTH SKIPPED
  </liquid_phase>

Here is more info on the `"liquid_phase`" elements:

    * `"liquid_component`" is an element with the following subelement: 

        * `"pressure`" is an element with the following attributes: ONLY UNIFORM, for now

.. code-block:: xml

     <pressure name="some name" value="exponential" function="linear | uniform" reference_coord="coordinate" gradient="coordinate"/>

.
    * `"solute_component`" is an element with the following attributes: NOT IMPLEMENTED YET

.. code-block:: xml

     <solute_component name="some name" filename="filename" value="exponential" function="linear | uniform" reference_coord="coordinate" gradient="coordinate"/>

.
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
      Optional Elements: solute_component, geochemistry - BOTH SKIPPED
  </liquid_phase>

Here is more info on the `"liquid_phase`" elements:

    * `"liquid_component`" is an element with the following subelement: 

        * `"inward_mass_flux`" is an element with the following attributes: ONLY CONSTANT, for now

.. code-block:: xml

     <inward_mass_flux value="exponential" function="linear | uniform | constant" start="time" />

.
        * `"inward_volumetric_flux`" is an element with the following attributes: ONLY CONSTANT, for now

.. code-block:: xml

     <inward_volumetric_flux value="exponential" function="linear | uniform | constant" start="time" />

.
        * `"uniform_pressure`" is an element with the following attributes: ONLY CONSTANT, for now

.. code-block:: xml

     <uniform_pressure name="some name" value="exponential" function="uniform | constant" start="time" />

.
        * `"hydrostatic`" is an element with the following attributes: ONLY CONSTANT, for now

.. code-block:: xml

     <hydrostatic name="some name" value="exponential" function="uniform | constant" start="time" />

.
    * `"solute_component`" is an element with the following subelement: NOT IMPLEMENTED YET

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

Output data from Amanzi is currently organized into three specific elements: `"Vis`", `"Checkpoint`", and `"Observations`".  
Each of these is controlled in different ways, reflecting their intended use.

* `"Vis`" is intended to represent snapshots of the solution at defined instances during the simulation to be visualized.  The ''vis'' element defines the naming and frequencing of saving the visualization files.  The visualizatoin files may include only a fraction of the state data, and may contiain auxiliary "derived" information (see *elsewhere* for more discussion).

* `"Checkpoint`" is intended to represent all that is necesary to repeat or continue an Amanzi run.  The specific data contained in a Checkpoint Data dump is specific to the algorithm options and mesh framework selected.  Checkpoint is special in that no interpolation is perfomed prior to writing the data files; the raw binary state is necessary.  As a result, the user is allowed to only write Checkpoint at the discrete intervals of the simulation. The ''checkpoint'' element defines the naming and frequencing of saving the checkpoint files.

* `"Observations`" is intended to represent diagnostic values to be returned to the calling routine from Amanzi's simulation driver.  Observations are typically generated at arbitrary times, and frequently involve various point samplings and volumetric reductions that are interpolated in time to the desired instant.  Observations may involve derived quantities (see discussion below) or state fields.  The ''observations'' element may define one or more specific ''observation''.

*EIB NOTE* - All three of the above are REQUIRED!!
For the obserservations I understand how to leave that empty.  But how do I execute without writing a checkpoint? If I'm running a dinky test am I really required to specify a checkpoint?  Will need to test this will validator.  Talk to Ellen about this.

Vis
---

The ''vis'' element defines the visualization filenaming scheme and how often to write out the files.  Thus, the ''vis'' element has the following requiements

.. code-block:: xml

  <vis>
      Required Elements: base_filename, num_digits, time_macro
      Optional Elements: NONE
  </vis>

The *base_filename* element contain the text component of the how the visualization files will be named.  The *base_filename* is appended with an index number to indicate the seqential order of the visualization files.  The *num_digits* elements indicates how many digits to use for the index. (*EIB NOTE* - verify if this is sequence index or interation id)  Final the *time_macro* element indicates the previously defined time_macro to be used to determin the frequency at which to write the visualization files.

(*EIB NOTE* - there should be a comment here about how the output is controlled, i.e. for each PK where do you go to turn on and off fields.  This will probably get filled in as the other sections fill out.)

Example:

.. code-block:: xml

  <vis>
     <base_filename>plot</base_filename>
     <num_digits>5</num_digits>
     <time_macro>Macro 1</time_macro>
  </vis>


Checkpoint
----------

The ''checkpoint'' element deines the filenaming scheme and frequency for writing out the checkpoint files.  As mentioned above, the user does not influence what is written to the checkpoint files.  Thus, the ''checkpoint'' element has the following requiements

.. code-block:: xml

  <checkpoint>
      Required Elements: base_filename, num_digits, time_macro
      Optional Elements: NONE
  </checkpoint>

The *base_filename* element contain the text component of the how the checkpoint files will be named.  The *base_filename* is appended with an index number to indicate the seqential order of the checkpoint files.  The *num_digits* elements indicates how many digits to use for the index. (*EIB NOTE* - verify if this is sequence index or interation id)  Final the *time_macro* element indicates the previously defined time_macro to be used to determin the frequency at which to write the checkpoint files.

Example:

.. code-block:: xml

  <checkpoint>
     <base_filename>chk</base_filename>
     <num_digits>5</num_digits>
     <time_macro>Every_100_timesteps</time_macro>
  </checkpoint>


Observations
------------

The Observations element holds all the observations that the user is requesting from Amanzi, as well as meta data, such as the name of the file that Amanzi will write observations to.  The observations are collected by their phase. Thus, the ''observations'' element has the following requirements

.. code-block:: xml

   <observations>
     Required Elements: filename, phase
     Optional Elements: NONE
   </observations>

The *filename* element contains the filename for the observation output, and may include the full path.  Currently, all observations are written to the same file.  

The *phase* element requires that the name of the phase be specified as an attribute and at least one observaton.

.. code-block:: xml

     <phase name="Name of Phase (Required)">
       Required Elements: observation (one observation element block for each observation)
       Optional Elements: NONE
     </phase>

In this release the only valid phase name is ''aqueous''.  The observation element requires a field quantity be given as an attribute, and elements for a region, a model (functional) with which it will extract its source data, and a list of discrete times for its evaluation.  The observations are evaluated during the simulation and returned to the calling process through one of Amanzi arguments. 

.. code-block :: xml

   <observation variable="Field Quantity (Required: see above for list of valid fields)">

     Required Elements: assigned_region, functional, time_macro 
     Optional Elements: NONE
     
   </observation>

Here the elements are ... 



Example:

.. code-block :: xml

   <observations>
     <filename>observation.out</filename> 
       <phase name="aqueous">
         <observation variable="H+ Aqueous concentration">
           <assigned_region>Well_1</assigned_region>
           <functional>point</functional>
           <time_macro>Every year</time_macro>
         </observation>
	 <observation variable="UO2++ Aqueous concentration">
	   <assigned_region>Well_3</assigned_region>
	   <functional>point</functional>
	   <time_macro>Every year</time_macro>
	 </observation>
       </phase>
     </observations>



