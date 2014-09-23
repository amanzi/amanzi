.. _Amanzi XML Schema :

============================================================
Input File XML Schema 
============================================================

Overview
++++++++

The present document describes how to specify the data required to execute Amanzi and perform a simulation.  This specification should be regarded as a companion to the mathematical requirements document entitled *Mathematical Formulation Requirements and Specifications for the Process Models* (see :ref:`ASCEM Overview <ASCEM Overview>`), and parameterizations of the individual submodels are consistent between Amanzi, the mathematical requirements document and this document.

The open-source, platform independent Akuna_ user environment can generate *Amanzi* models and generate corresponding valid, human-readable XML input files that can then be executed by *Amanzi*.  Example input files are available in the Amanzi source repository.

XML Schema 2.0
++++++++++++++

Amanzi solves a set of parameterized models for multiphase flow in porous media. An Amanzi simulation is specified by providing:

* values for a parameterized PDE-based transport model, including boundary and initial conditions, constitutive laws, and parameterized/phenomenological models for fluid and chemical sources and characterizations of the porous medium,
* parameters controlling the selection of key algorithmic options and output,
* a description of the (discrete) state of the computational system, including a list of the independent variables and instructions for obtaining or generating the discrete mesh, and a characterization of the (parallel) computing environment.

The primary input to *Amanzi* is through an XML file. The Amanzi input XML format is defined in terms of the XML schema that can found in the Amanzi source code repository.  Users can construct models and generate compliant XML input files using the Akuna_ tool suite.  Users can also choose to generate compliant file using a text editor or other method.

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

    <amanzi_input version="2.0.0" type="unstructured">
    </amanzi_input>

Model Description
-------------------

This section allows the user to provide information about the model being developed and how and when it was developed.  Default units for the model are also stored in this section.  This entire section is optional but encouraged for documentation purposes.

The opening tag ``model_description`` accepts an attribute ``name`` in which the user may give the current model a name.  The available elements within ``model_description`` include: ``comments``, ``author``, ``created``, ``modified``, ``model_id``, ``description``, ``purpose``, and ``units``.  Under the ``units`` element, the user may define the default units to used throughout the rest of the input file.  The options available for the units element are as follows:

.. code-block:: xml

    <units>
      <length_unit> m | cm </length_unit>
      <time_unit> y | d | h | s </time_unit>
      <mass_unit> kg </mass_unit>
      <conc_unit>molar</conc_unit>
    </units>


Here is an overall example for the model description element.

.. code-block:: xml

  <model_description name="BC Cribs">
    <comments>Some comment here</comments>
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
-----------

This section allows the user to provide useful definitions to be used throughout the other sections.  Definitions are grouped as elements constants, named times, and macros.

Constants can be one of four types: constant, time_constant, numerical_constant, and area_mass_flux_constant. Each of these types look as the follows:

.. code-block:: xml

  <constants>
    <constant name="Name of Constant" type="none | time | area_mass_flux" value="constant_value">
    <time_constant  name="Name of Time"  value="value,y|d|s">
    <numerical_constant name="Name of Constant" value="value_constant">
    <area_mass_flux_constant name="Name of Constant" value="value_of_flux">
  </constants>

Named_times allows the user to assign meaningful names to time values and define time values in a single location in the file.  Then the names are used throughout the file whenever needed by boundary conditions or execution controls, etc.  The named_times element contains an unbounded number of time ``time`` elements. The trailing character in the value attribute indicates the units of the time.

.. code-block:: xml

  <named_times>
    <time  name="Name of Time" value="time,y|d|s">
  </named_times>

The ``macro`` section defines time and cycle macros.  These specify a series of times or cycles for writing out visualization or checkpoint files.  Each ``time_macros`` requires a ``name`` attribute and one or more ``time`` elements.  An alternative option is to specify the start and stop times and interval time step, as shown in the ``cycle_macro``.

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

Using ``-1`` as the stop value will continue the interval until the simulation ends.

Here is an overall example for the ``definition`` element.

.. code-block:: xml

   <definitions>
	<constants>
		<constant name="zero" type="none" value="0.000"/>
		<constant name="start" type="time" value="1956.0;y"/>
		<constant name="future_recharge" type="area_mass_flux" value="1.48666E-6"/>
		<time_constant name="start_time" value="1956.0;y"/>
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

The ``execution_controls`` section defines the general execution of the Amanzi simulation.  Amanzi can execute in three modes: steady state, transient or initialize to a steady state and then continue it transient.  Default values for execution are defined in the ``execution_control_defaults`` element.  These values are used for any time period during the simulation for which the controls were not specified.  Individual time periods of the simulation are defined using ``execution_control`` elements.  For a steady state simulation, only one ``execution_control`` element will be defined.  However, for a transient simulation a series of controls may be defined during which different control values will be used.  For a valid ``execution_controls`` section the ``execution_control_defaults`` element and at least one ``execution_control`` element must appear.

The ``execution_control_defaults`` element has the following attributes.

.. code-block:: xml

  <execution_control_defaults init_dt="labeled_time" max_dt="labeled_time" reduction_factor="exponential" increase_factor="exponential" mode="steady | transient" method=" bdf1 | picard" />

The ``execution_control`` element has the following attributes.  

.. code-block:: xml

  <execution_control  restart="string" start="string" end="string" init_dt="labeled_time" max_dt="labeled_time" reduction_factor="exponential" increase_factor="exponential" mode="steady | transient" method=" bdf1 | picard" />

Each ``execution_control`` is required to define a ``start`` time.  The final control period must define and ``end`` time.  It is assumed that the start time of the next control period is the end time of the previous period.  Therefore, it is not required that each ``execution_control`` element have an ``end`` time defined.

The ``execution_control`` section also provides the elements ``comments`` and ``verbosity``.  Users may provide any text within the ``comment`` element to annotate this section.  ``verbosity`` takes the attribute level=`` extreme | high | medium | low | none``.  This triggers increasing levels of reporting from inside Amanzi.  For debugging purposes use the level extreme.

Restarting a simulation is available using the ``restart`` attribute.  The value given for the ``restart`` attribute is the name of the Amanzi checkpoint file to be read in and initialized from.

Here is an overall example for the ``execution_control`` element.

.. code-block:: xml

  <execution_controls>
    <execution_control_defaults init_dt= "3.168E-08"   max_dt="0.01"  reduction_factor="0.8"  increase_factor="1.25" mode="transient" method="bdf1"/>
    <execution_control  start="0.0;y"   end="1956.0,y"  init_dt= "0.01" max_dt="500.0" reduction_factor="0.8"  mode = "steady"   />
    <execution_control start="B-17_release_begin" />
    <execution_control start="B-17_release_end" />
    <execution_control start="B-18_release_begin" />
    <execution_control start="B-18_release_end"  end="3000.0,y" />
  </execution_controls>

Numerical Controls
------------------

This section allows the user to define control parameters associated with the underlying numerical implementation.  The list of available options is lengthy.  However, none are required for a valid input file.  The ``numerical_controls`` section is divided up into the subsections: ``steady-state_controls``, ``transient_controls``, ``linear_solver``, ``nonlinear_solver``, and ``chemistry_controls``.  The list of available options is as follows:

.. code-block:: xml

  <numerical_controls>

    <comments>Numerical controls comments here</comments>

    <steady-state_controls>
      <comments>Comment text here</comments>
      <min_iterations>Value</min_iterations>
      <max_iterations>Value</max_iterations>
      <max_preconditioner_lag_iterations>Value</max_preconditioner_lag_iterations>
      <nonlinear_tolerance>Value</nonlinear_tolerance>
      <pseudo_time_integrator>
        <initialize_with_darcy>Value</initialize_with_darcy>
        <clipping_saturation>Value</clipping_saturation>
        <method>picard</method>
        <preconditioner> trilinos_ml | hypre_amg | block_ilu </preconditioner>
             <!-- if trilinos_ml -->
          <trilinos_smoother_type> jacobi | gauss_seidel | ilu </trilinos_smoother_type>
          <trilinos_threshold> Value </trilinos_threshold>
          <trilinos_smoother_sweeps> Value </trilinos_smoother_sweeps>
          <trilinos_cycle_applications> Value </trilinos_cycle_applications>
             <!-- if hypre_amg -->
          <hypre_cycle_applications> Value </hypre_cycle_applications>
          <hypre_smoother_sweeps >Value </hypre_smoother_sweeps>
          <hypre_tolerance >Value </hypre_tolerance>
          <hypre_strong_threshold> Value </hypre_strong_threshold>
             <!-- if block_ilu -->
          <ilu_overlap> Value </ilu_overlap>
          <ilu_relax> Value </ilu_relax>
          <ilu_rel_threshold> Value </ilu_rel_threshold>
          <ilu_abs_threshold> Value </ilu_abs_threshold>
          <ilu_level_of_fill> Value </ilu_level_of_fill>
        </preconditioner>
        <linear_solver>aztec00</linear_solver>
        <control_options>Value</control_options>
        <convergence_tolerance>Value</convergence_tolerance>
        <max_iterations>Value</max_iterations>
      </pseudo_time_integrator>
      <limit_iterations>Value</limit_iterations>
      <nonlinear_iteration_damping_factor>Value</nonlinear_iteration_damping_factor>
      <nonlinear_iteration_divergence_factor>Value</nonlinear_iteration_divergence_factor>
      <restart_tolerance_factor>Value</restart_tolerance_factor>
      <restart_tolerance_relaxation_factor>Value</restart_tolerance_relaxation_factor>
      <max_divergent_iterations>Value</max_divergent_iterations>
    </steady-state_controls>

    <transient_controls>
      <comments>Comment text here</comments>
      <bdf1_integration_method min_iterations="Value" 
                               max_iterations="Value" 
                               limit_iterations="Value"
                               nonlinear_tolerance="Value"
                               nonlinear_iteration_damping_factor="Value"
                               max_preconditioner_lag_iterations="Value"
                               max_divergent_iterations="Value"
                               nonlinear_iteration_divergence_factor="Value"
                               restart_tolerance_factor="Value"
                               restart_tolerance_relaxation_factor="Value" />
      <preconditioner> trilinos_ml | hypre_amg | block_ilu </preconditioner>
           <!-- if trilinos_ml -->
        <trilinos_smoother_type> jacobi | gauss_seidel | ilu </trilinos_smoother_type>
        <trilinos_threshold> Value </trilinos_threshold>
        <trilinos_smoother_sweeps> Value </trilinos_smoother_sweeps>
        <trilinos_cycle_applications> Value </trilinos_cycle_applications>
           <!-- if hypre_amg -->
        <hypre_cycle_applications> Value </hypre_cycle_applications>
        <hypre_smoother_sweeps >Value </hypre_smoother_sweeps>
        <hypre_tolerance >Value </hypre_tolerance>
        <hypre_strong_threshold> Value </hypre_strong_threshold>
           <!-- if block_ilu -->
        <ilu_overlap> Value </ilu_overlap>
        <ilu_relax> Value </ilu_relax>
        <ilu_rel_threshold> Value </ilu_rel_threshold>
        <ilu_abs_threshold> Value </ilu_abs_threshold>
        <ilu_level_of_fill> Value </ilu_level_of_fill>
      </preconditioner>
    </transient_controls>

    <linear_solver>
      <comments>Comment text here</comments>
      <method> gmres </method>
      <max_iterations> Value </max_iterations>
      <tolerance> Value </tolerance>
      <preconditioner name="trilinos_ml | hypre_amg | block_ilu">
           <!-- if trilinos_ml -->
        <trilinos_smoother_type> jacobi | gauss_seidel | ilu </trilinos_smoother_type>
        <trilinos_threshold> Value </trilinos_threshold>
        <trilinos_smoother_sweeps> Value </trilinos_smoother_sweeps>
        <trilinos_cycle_applications> Value </trilinos_cycle_applications>
           <!-- if hypre_amg -->
        <hypre_cycle_applications> Value </hypre_cycle_applications>
        <hypre_smoother_sweeps >Value </hypre_smoother_sweeps>
        <hypre_tolerance >Value </hypre_tolerance>
        <hypre_strong_threshold> Value </hypre_strong_threshold>
           <!-- if block_ilu -->
        <ilu_overlap> Value </ilu_overlap>
        <ilu_relax> Value </ilu_relax>
        <ilu_rel_threshold> Value </ilu_rel_threshold>
        <ilu_abs_threshold> Value </ilu_abs_threshold>
        <ilu_level_of_fill> Value </ilu_level_of_fill>
      </preconditioner>
    </linear_solver>

    <nonlinear_solver name="nka | newton | inexact newton">

    <chemistry_controls>
      <chem_tolerance> Value </chem_tolerance>
      <chem_max_newton_iterations> Value </chem_max_newton_iterations>
    </chemistry_controls>

  </numerical_controls>

Here is an overall example for the ``numerical_controls`` element.

.. code-block:: xml

	<numerical_controls>

		<comments>Numerical controls comments here</comments>

		<steady-state_controls>
		        <comments>Note that this section contained data on timesteps, which was moved into the execution control section.</comments>
          		<min_iterations>10</min_iterations>
		      	<max_iterations>15</max_iterations>
          		<max_preconditioner_lag_iterations>30</max_preconditioner_lag_iterations>
          		<nonlinear_tolerance>1.0e-5</nonlinear_tolerance>
		</steady-state_controls>
		<transient_controls>
			<comments>Proposed comments section.</comments>
			<bdf1_integration_method min_iterations="10" max_iterations="15" max_preconditioner_lag_iterations="5" />
		</transient_controls>
		<linear_solver>
			<comments>Proposed comment section.</comments>
			<method>gmres</method>
			<max_iterations>20</max_iterations>
			<tolerance>1.0e-18</tolerance>
	                <preconditioner name = "hypre_amg">
	                     	<hypre_cycle_applications>10</hypre_cycle_applications>
	                	<hypre_smoother_sweeps>3</hypre_smoother_sweeps>
	                       	<hypre_tolerance>0.1</hypre_tolerance>
	                       	<hypre_strong_threshold>0.4</hypre_strong_threshold>
	                 </preconditioner>
 		</linear_solver>

	</numerical_controls>

Mesh
----

A mesh must be defined for the simulation to be conducted on.  The mesh can be structured or unstructured.  Structured meshes are always internally generated while unstructured meshes may be generated internally or imported from an existing `Exodus II <http://sourceforge.net/projects/exodusii/>`_ file. Generated meshes in both frameworks are always regular uniformly spaced meshes.

Mesh parameters are specified in the ``mesh`` section. If the mesh is unstructured the opening tag of the ``mesh`` section takes an attribute called called ``framework`` which can take the value of ``mstk``, ``stk::mesh``, ``moab`` or ``simple``. This specifies which mesh infrastructure library is to be used for managing the mesh queries under-the-hood. 

The ``mesh`` section takes a ``dimension`` element which indicates if the mesh is 2D or 3D. A 2D mesh can be given in 3D space with a third coordinate of 0. If a 2D mesh is specified this impacts other aspects of the input file.  It is up to the user to ensure consistency within the input file.  Other effected parts of the input file include region definitions and initial conditions which use coordinates, the material property permeability which must be specified using the correct subset of x, y, and z coordinates, and the initial condiction velocity which also requires the correct subset of x, y, and z coordinates.

This section also takes an element indicating if the mesh is to be internally generated (structured and unstructured) or read from an external file (unstructured only). If the mesh is to be generated internally, a ``generate`` element is specified with details about the number of cells in each direction and the low and high coordinates of the bounding box. If the mesh is to be read from a file, a ``read`` element is specified with the file name and file format. Currently only Exodus II files are supported.  Finally, as in other sections, a ``comments`` element is provide to include any comments or documentation the user wishes.

Here is an example specification for a structured ``mesh`` element.

.. code-block:: xml

  <mesh framework="mstk"> 
    <comments>3D block</comments>
    <dimension>3</dimension>
    <generate>
      <number_of_cells nx = "400"  ny = "200"  nz = "10"/>
      <box  low_coordinates = "0.0,0.0,0.0" high_coordinates = "200.0,200.0,1.0"/>
    </generate>
  </mesh>

The following is an example specification for a generated unstructured
mesh.

.. code-block:: xml

  <mesh framework="mstk"> 
    <comments>Pseudo 2D</comments>
    <dimension>3</dimension>
    <generate>
      <number_of_cells nx = "432"  ny = "1"  nz = "256"/>
      <box  low_coordinates = "0.0,0.0,0.0" high_coordinates = "216.0,1.0,107.52"/>
    </generate>
  </mesh>

Finally, an example of reading an unstructured mesh from a file is given below.

.. code-block:: xml

  <mesh framework="mstk"> 
    <comments>Read from Exodus II</comments>
    <dimension>3</dimension>
    <read>
      <file>dvz.exo</file>
      <format>exodus ii</format>
    </read>

  </mesh>

Regions
-------

Regions are geometrical constructs used in Amanzi to define subsets of the computational domain in order to specify the problem to be solved, and the output desired. Regions are commonly used to specify material properties, boundary conditions and obervation domains. Regions may represent zero-, one-, two- or three-dimensional subsets of physical space. For a three-dimensional problem, the simulation domain will be a three-dimensional region bounded by a set of two-dimensional regions. If the simulation domain is N-dimensional, the boundary conditions must be specified over a set of regions are (N-1)-dimensional.

Amanzi automatically defines the special region labeled "All", which is the entire simulation domain. Under the "Structured" option, Amanzi also automatically defines regions for the coordinate-aligned planes that bound the domain, using the following labels: "XLOBC", "XHIBC", "YLOBC", "YHIBC", "ZLOBC", "ZHIBC"

The ``regions`` block is required.  Within the region block no regions are required to be defined.  The optional elements include ``region``, ``box``, ``point``, and ``plane``.  As in other sections there is also an options ``comments`` element.

The elements ``box``, ``point``, and ``plane`` allow for inline description of regions.  The ``region`` element uses a subelement to either define a ``box`` or ``plane`` region or specify a region file.  Below are further descriptions of these elements.

Box
---

A box region region is defined by a low corner coordinates and high corner coordinates. Box regions can be degenerate in one or more directions.

.. code-block:: xml

  <box  name="box name" low_coordinates = "x_low,y_low,z_low"
  high_coordinates = "x_high,y_high,z_high"/>


Point
-----

A point region is defined by a point coordinates.

.. code-block:: xml

  <point name="point name" coordinate = "x,y,z" />

Plane
-----

A plane region is defined by a point on the plane and the normal direction of the plane

.. code-block:: xml

  <plane name="plane name" location="x,y,z" normal="dx,dy,dz" />

Labeled Set
-----------

A labeled set region is a predefined set of mesh entities defined in the Exodus II mesh file. This type of region is useful when applying boundary conditions on an irregular surface that has been tagged in the external mesh generator.  Please note that both the format and entity attribute values are case sensitive.

.. code-block:: xml

  <region name="region name">
      <region_file label="integer label" name="filename" type="labeled set" format="exodus ii" entity=["cell"|"face"] />
  </region>

Color function
--------------

A color function region defines a region based on a specified integer color in a structured color function file. The color values may be specified at the nodes or cells of the color function grid. A computational cell is assigned the color of the data grid cell containing its cell centroid or the data grid nearest its cell-centroid. Computational cell sets are then build from all cells with the specified color value. In order to avoid gaps and overlaps in specifying materials, it is strongly recommended that regions be defined using a single color function file.  At this time, Exodus II is the only file format available.   Please note that both the format and entity attribute values are case sensitive.


.. code-block:: xml

  <region name="region name">
      <region_file label="integer label" name="filename" type="color" format="exodus ii"  entity=["cell"|"face"]/>
  </region>

.. EIB:  The following are not exposed through the current XML Schema, only the OLD input spec.  I've commented out the text until a future date when they might be exposed.

.. Polygon
.. -------

.. A polygon region is used to define a bounded planar region and is specified by the number of points and a list of points

.. Logical
.. -------

.. Logical regions are compound regions formed from other primitive type regions using boolean operations. Supported operators are union, intersection, subtraction and complement.


.. Geochemistry
.. ------------

Material
--------

The ``material`` in this context is meant to represent the media through which fluid phases are transported. In the literature, this is also referred to as the "soil", "rock", "matrix", etc. Properties of the material must be specified over the entire simulation domain, and is carried out using the Region constructs defined above. For example, a single material may be defined over the "All" region (see above), or a set of materials can be defined over subsets of the domain via user-defined regions. If multiple regions are used for this purpose, they should be disjoint, but should collectively tile the entire domain. Each material requires a label and the following set of physical properties using the supported models described below.

A ``material`` element can contain the following:

.. code-block:: xml

  <material name="Name of material">
    <comments>Comment text here</comments>
    <mechanical_properties>
        <porosity value="Value"/>
        <particle_density value="Value"/>
        <specific_storage value="Value"/>
        <specific_yield value="Value"/>
        <dispersion_tensor alpha_l="Value" alpha_t="Value"/>
        <tortuosity value="Value"/>
    </mechanical_properties>
    <permeability x="Value" y="Value" z="Value"/>
    <hydraulic_conductivity x="Value" y="Value" z="Value"/>
    <cap_pressure model="van_genuchten | brooks_corey | none">
        <!-- for van_genuchten -->
 	<parameters m="Value" alpha="Value" sr="Value" optional_krel_smoothing_interval="Value"/>
        <!-- for brooks_corey -->
 	<parameters lambda="Value" alpha="Value" sr="Value" optional_krel_smoothing_interval="Value"/>
    </cap_pressure>
    <rel_perm model="mualem | burdine | none">
        <!-- burdine only -->
        <exp value="Value"/>
    </rel_perm>
    <sorption_isotherms>
        <!-- Three kd models available plus molecular diffusion -->
        <!-- Note: all solutes should be listed under all materials, value="0" indicates the solute isn't present/active -->
	<solute name="Name of Solute" >
            <kd_model model="linear" kd = "Value" />
            <molecular_diffusion value="Value" />
        </solute>
	<solute name="Name of Solute" >
            <kd_model model="langmuir" kd="Value" b="Value"/>
        </solute>
	<solute name="Name of Solute" >
            <kd_model model="freundlich" kd="Value" n="Value" />
        </solute>
    </sorption_isotherms>
    <assigned_regions>Comma seperated list of Regions</assigned_regions>
  </material>


While many material properties are available for the user to define, the minimum requirements for a valid material definition are specifying the ``assigned_regions`` and the ``porosity``.  However, if a capillary pressure model or relative permeability model is chosen (other than ``none``), the associated parameters must also be provided.

An example material would look like

.. code-block:: xml

  <material name="Facies_1">
    <comments>Material corresponds to region facies1</comments>
    <mechanical_properties>
      <porosity value="0.4082"/>
      <particle_density value="2720.0"/>
    </mechanical_properties>
    <permeability x="1.9976E-12" y="1.9976E-12" z="1.9976E-13"/>
    <cap_pressure model="van_genuchten">
      <parameters m="0.2294" alpha="1.9467E-04" sr="0.0"/>
    </cap_pressure>
    <rel_perm model="mualem"/>
    <assigned_regions>Between_Planes_1_and_2</assigned_regions>
  </material>

Process Kernels
---------------

Amanzi current employees three process kernels that need to be defined in the input file (1) flow, (2) transport, and (3) chemistry.  The ``process_kernels`` section allows the user to define which kernels are to be used during the section and select high level features of those kernels.  The ``process_kernels`` element is as follows:

.. code-block:: xml

  <process_kernels>
    <comments>Comment text here</comments>
    <flow state = "on | off" model = "richards | saturated | constant"/>
    <transport state = "on | off" algorithm = "explicit first-order | explicit second-order | implicit upwind | none" sub_cycling = "on | off"/>
    <chemistry state = "on | off" engine = "amanzi | pflotran | none" process_model="implicit operator split | none"/>
  </process_kernels>

Currently three scenerios are avaiable for calculated the flow field.  `"richards`" is a single phase, variably saturated flow assuming constant gas pressure.  `"saturated`" is a single phase, fully saturated flow.  `"constant`" is equivalent to the a flow model of single phase (saturated) with the time integration mode of transient with static flow in the version 1.2.1 input specification.  This flow model indicates that the flow field is static so no flow solver is called during time stepping. During initialization the flow field is set in one of two ways: (1) A constant Darcy velocity is specified in the initial condition; (2) Boundary conditions for the flow (e.g., pressure), along with the initial condition for the pressure field are used to solve for the Darcy velocity.

For `"transport`" a combination of `"state`" and `"algorithm`" must be specified.  If `"state`" is `"off`" then `"algorithm`" is set to `"none`".  Otherwise the integration algorithm must be specified.  Whether sub-cycling is to be utilized within the transport algorithm is also specified here.

For `"chemistry`" a combination of `"state`", `"engine`", and `"process_model`" must be specified.  If `"state`" is `"off`" then `"engine`" and `"process_model`" are set to `"none`".  Otherwise the `"engine`" and `"process_model`" model must be specified. 

An example ``process_kernels`` is as follows:

.. code-block:: xml

  <process_kernels>
    <comments>This is a proposed comment field for process_kernels</comments>
    <flow state = "on" model = "richards"/>
    <transport state = "on" algorithm = "explicit first-order" sub_cycling = "on"/>
    <chemistry state = "off" engine="none" process_model="none"/>
  </process_kernels>

Phases
------

The ``phases`` section is used to specify components of each of the phases that are mobile, and solutes that are contained within them. For each phase, the list identifies the set of all independent variables that are to be stored on each discrete mesh cell.

The terminology for flow in porous media can be somewhat ambiguous between the multiphase and groundwater communities, particularly in regards to "components", "solutes" and "chemicals". Since Amanzi is designed to handle a wide variety of problems, we must settle on a nomenclature for our use here. In the general problem, multiple "phases" may coexist in the domain (e.g. gaseous, aqueous/liquid, etc), and each is comprised of a number of "components" (section 2.2). In turn, each component may carry a number of "solutes" and some of these may participate in chemical reactions. As a result of reactions, a chemical source or sink term may appear for the solutes involved in the reaction, including solutes in other mobile phases or in the material matrix. Additionally, certain reactions such as precipitation may affect the flow properties of the material itself during the simulation, and some might affect the properties of the fluid (e.g. brines affect the liquid density). While Amanzi does not currently support chemical reactions and thermal processes, the specification here allows for the existence of the necessary data structures and input data framework. Note that if solute concentrations are significant, the system may be better modeled with that solute treated as a separate component. Clearly, these definitions are highly problem-dependent, so Amanzi provide a generalized interface to accommodate a variety of scenarios.

Currently in Amanzi, solutes are transported in the various phase components, and are treated in "complexes". Each complex is typically in chemical equilibrium with itself and does not undergo phase change. Under these conditions, knowledge of the local concentration of the "basis" or "primary" species (the terms are used here interchangeably) in a chemical complex is sufficient to determine the concentrations of all related secondary species in the phase. Each basis species has a total component concentration and a free ion concentration. The total component concentration for each basis species is a sum of the free ion concentrations in the phase components and its stoichiometric contribution to all secondary species. Amanzi splits the total component concentration into a set of totals for each of the transported phase components, and a total sorbed concentration. Given the free ion concentration of each basis species (and if there is more than one phase, a specification of the thermodynamic relationships that determine the partitioning between phase components (if mass transfer is allowed - not in current Amanzi), we can reconstruct the concentration of the primary and secondary species in each phase. As a result only the basis species are maintained in the state data structures for each phases component.

In addition to solutes in the transported phases, there may be various immobile chemical constituents within the porous media (material) matrix, such as "minerals" and "surface complexes". Bookkeeping for these constituents is managed in Amanzi data structures by generalizing the "solute" concept - a slot in the state is allocated for each of these immobile species, but their concentrations are not included in the transport/flow components of the numerical integration. To allow selective transport of the various solutes, Amanzi uses the concept of solute groups. The aqueous solute concentrations are typically treated together as a group, for example, and often represent the only chemical constituents that are mobile. Thus, the current Amanzi will assume that any other groups specified in an Aqueous phase are immobile.

This section specifies the phases present and specific properties about those phases.  The first grouping is by ``liquid_phase`` and ``solid_phase``.  The ``liquid_phase`` element is required to produce a valid input file.

The ``liquid_phase`` element requires an attribute *name*.  This is used by other sections to identify this phase.  Subelements are used to define the ``viscosity``, ``density``, and ``dissolved_components``. While ``viscosity`` and ``density`` are required elements, ``dissolved_components`` is optional.  ``dissolved_components`` contains a subelement ``solutes`` under which individual ``solute`` elements are used to specify any solutes present in the liquid phase.  The text of the ``solute`` contains the name of the solute while an attributes specifies the value *coefficient_of_diffusion*.  

The ``solid_phase`` element allows the user to define a ``minerals`` element under which a series of ``mineral`` elements can be listed to specify any minerals present in the solida phase.  The ``mineral`` elements contain the name of the mineral.

An example ``phases`` element looks like the following.

.. code-block:: xml

  <phases>
    <liquid_phase name = "water">
	<viscosity>1.002E-03</viscosity>
	<density>998.2</density>
	<dissolved_components> 
	    <solutes>
	       <solute coefficient_of_diffusion="1.0e-9">Tc-99</solute>
	    </solutes> 
	</dissolved_components>
    </liquid_phase>
  </phases>


Initial Conditions
------------------

The `"initial_conditions`" section contains at least 1 and up to an unbounded number of `"initial_condition`" elements.  Each `"initial_condition`" element defines a single initial condition that is applied to one or more region specified in the ``assigned_regions`` element.  The initial condition can be applied to a liquid phase or solid phase using the appropriate subelement.

To specify a liquid phase the ``liquid_phase`` element is used.  At least one ``liquid_component`` must be specified.  In addition an unbounded number of ``solute_component`` elements and a single ``geochemistry`` element can be specified.  Under the ``liquid_component`` and ``solute_component`` elements an initial condition can be defined.  Under the ``geochemistry`` element a geochemistry constraint is defined.

The initial conditions are defined using a specific elements.  The element name indicates the type of condition and the attributes define the necessary information.  Below is a table of the conditions available for the liquid phase and the attributes required to define them.

+-----------------------+------------------+---------------------------------+
| Initial Condition Type| Attributes       | Value Type                      |
+=======================+==================+=================================+
| uniform_pressure      | | name           | | string                        |
|                       | | value          | | double/time_constant/constant |
+-----------------------+------------------+---------------------------------+
| linear_pressure       | | name           | | string                        |
|                       | | value          | | double/time_constant/constant |
|                       | | reference_coord| | coordinate                    |
|                       | | gradient       | | coordinate                    |
+-----------------------+------------------+---------------------------------+
| velocity              | | name           | | string                        |
|                       | | x              | | double/constant               |
|                       | | y              | | double/constant               |
|                       | | (z)            | | double/constant               |
+-----------------------+------------------+---------------------------------+

.. | uniform_saturation    | | name           | | string                        |
.. |                       | | value          | | double/time_constant/constant |
.. +-----------------------+------------------+---------------------------------+
.. | linear_saturation     | | name           | | string                        |
.. |                       | | value          | | double/time_constant/constant |
.. |                       | | reference_coord| | coordinate                    |
.. |                       | | gradient       | | coordinate                    |
.. +-----------------------+------------------+---------------------------------+


For the solute_component the attributes available are *name*, *value*, *function*, *reference_coord*, and *gradient*.  The function options available are *uniform* and *linear*.  The attributes *reference_coord*, and *gradient* are only necessary for the *linear* function type.

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
		<solute_component name = "Tc-99" value = "0" function="uniform"/>
	    </liquid_phase>
          </initial_condition>
	</initial_conditions>

Boundary Conditions
-------------------

Boundary conditions are defined in a similar manor to the initial conditions.  Under the tag ``boundary_conditions`` and series of individual ``boundary_condition`` elements can be defined.  Within each ``boundary_condition`` element the ``assigned_regions`` and ``liquid_phase`` elements must appear.  The boundary condition can be applied to one or more region using a comma separated list of region names.  Under the ``liquid_phase`` element the ``liquid_component`` element must be define.  An unbounded number of ``solute_component`` elements and one ``geochemistry`` element may optionally be defined.

Under the ``liquid_component`` and ``solute_component`` elements a time series of boundary conditions is defined using the boundary condition elements available in the table below.  Each component element can only contain one type of boundary condition.  Both elements also accept a *name* attribute to indicate the phase associated with the boundary condition.

+-------------------------+--------------------+------------------------------------+
|Boundary Condition Type  | Attributes         | Value Type                         |
+=========================+====================+====================================+
|inward_mass_flux         | | name             | | string                           |
|inward_volumetric_flux   | | start            | | double/time_constant/constant    |
|outward_mass_flux        | | value            | | double                           |
|outward_volumetric_flux  | | function         | | 'linear','uniform','constant'    |
+-------------------------+--------------------+------------------------------------+
|uniform_pressure         | | name             | | string                           |
|                         | | start            | | double/time_constant/constant    |
|                         | | value            | | double                           |
|                         | | function         | | 'uniform','constant'             |
+-------------------------+--------------------+------------------------------------+
|hydrostatic              | | name             | | string                           |
|                         | | start            | | double/time_constant/constant    |
|                         | | value            | | double                           |
|                         | | function         | | 'uniform','constant'             |
|                         | | coordinate_system| | 'absolute','relative to mesh top'|
+-------------------------+--------------------+------------------------------------+ 
|seepage_face             | | name             | | string                           |
|                         | | start            | | double/time_constant/constant    |
|                         | | inward_mass_flux | | double/time_constant/constant    |
|                         | | function         | | 'linear','uniform','constant'    |
+-------------------------+--------------------+------------------------------------+
|no_flow                  | | name             | | string                           |
|                         | | start            | | double/time_constant/constant    |
|                         | | function         | | 'linear','uniform','constant'    |
+-------------------------+--------------------+------------------------------------+

For the solute component, the boundary condition available is ``aqueous_conc`` which has the attributes *name*, *value*, *function*, and *start*.  The function options available are *uniform*, *linear*, and *constant*.

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
	    <solute_component name = "solute">
		<aqueous_conc name = "Tc-99" start="0.0"     function= "constant"  value="zero"/>
		<aqueous_conc name = "Tc-99" start="1956.0,y"  function= "constant"  value="zero"/>
		<aqueous_conc name = "Tc-99" start="2012.0,y"  function= "constant"  value="zero"/>
		<aqueous_conc name = "Tc-99" start="3000.0,y"  function= "constant"  value="zero"/>
	    </solute_component>
	</liquid_phase>
    </boundary_condition>
  </boundary_conditions>


Sources
-------

Sources are defined in a similar manner to the boundary conditions.  Under the tag ``sources`` and series of individual ``source`` elements can be defined.  Within each ``source`` element the ``assigned_regions`` and ``liquid_phase`` elements must appear.  Sources can be applied to one or more region using a comma separated list of region names.  Under the ``liquid_phase`` element the ``liquid_component`` element must be define.  An unbounded number of ``solute_component`` elements and one ``geochemistry`` element may optionally be defined.

Under the ``liquid_component`` and ``solute_component`` elements a time series of boundary conditions is defined using the boundary condition elements available in the table below.  Each component element can only contain one type of source.  Both elements also accept a *name* attribute to indicate the phase associated with the source.

+-------------------------+--------------------+------------------------------------+
|Liquid Phase Source Type | Attributes         | Value Type                         |
+=========================+====================+====================================+
|volume_weighted          | | start            | | double/time_constant/constant    |
|perm_weighted            | | value            | | double                           |
|                         | | function         | | 'linear','uniform','constant'    |
+-------------------------+--------------------+------------------------------------+

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

Output data from Amanzi is currently organized into three specific elements: `"Vis`", `"Checkpoint`", `"Observations`", and `"Walkabout Data`".  Each of these is controlled in different ways, reflecting their intended use.

* `"Visualization`" is intended to represent snapshots of the solution at defined instances during the simulation to be visualized.  The ``vis`` element defines the naming and frequency of saving the visualization files.  The visualization files may include only a fraction of the state data, and may contain auxiliary "derived" information.

* `"Checkpoint`" is intended to represent all that is necessary to repeat or continue an Amanzi run.  The specific data contained in a checkpoint dump is specific to the algorithm options and mesh framework selected.  Checkpoint is special in that no interpolation is performed prior to writing the data files; the raw binary state is necessary.  As a result, the user is allowed to only write checkpoint at the discrete intervals of the simulation. The ``checkpoint`` element defines the naming and frequency of saving the checkpoint files.

* `"Observations`" is intended to represent diagnostic values to be returned to the calling routine from Amanzi's simulation driver.  Observations are typically generated at arbitrary times, and frequently involve various point samplings and volumetric reductions that are interpolated in time to the desired instant.  Observations may involve derived quantities or state fields.  The ``observations`` element may define one or more specific observation.

* `"Walkabout Data`" is intended to be used as input to the particle tracking software Walkabout.

Viz
---

The ``vis`` element defines the visualization file naming scheme and how often to write out the files.  The ``base_filename`` element contain the text component of the how the visualization files will be named.  The ``base_filename`` is appended with an index number to indicate the sequential order of the visualization files.  The ``num_digits`` elements indicates how many digits to use for the index.  Finally, the ``time_macros`` or ``cycle_macros`` element indicates previously defined time_macros or cycle_macros to be used to determine the frequency at which to write the visualization files.  One or more macro can be listed in a comma separated list.  Amanzi will converted the list of macros to a single list of times or cycles contained by all of the macros listed and output accordingly.

An example ``vis`` element looks like the following.

.. code-block:: xml

   <vis>
        <base_filename>plot</base_filename>
	<num_digits>5</num_digits>
	<time_macros>Macro 1</time_macros>
   </vis>

Checkpoint
----------

The ``checkpoint`` element defines the file naming scheme and frequency for writing out the checkpoint files.  The ``base_filename`` element contains the text component of the how the checkpoint files will be named.  The ``base_filename`` is appended with an index number to indicate the sequential order of the checkpoint files.  The ``num_digits`` elements indicates how many digits to use for the index.  Finally, the ``cycle_macro`` element indicates the previously defined cycle_macro to be used to determine the frequency at which to write the checkpoint files.

An example ``checkpoint`` element looks like the following.

.. code-block:: xml

    <checkpoint>
	<base_filename>chk</base_filename>
	<num_digits>5</num_digits>
	<cycle_macro>Every_1000_steps</cycle_macro>
    </checkpoint>

Observations
------------

The ``observations`` element defines the the file for writing observations to and specifies individual observations to be made.  At this time, all observations are written to a single file defined in the ``filename`` element.  Also, observations are only available for the liquid phases.  Therefore individual observations are defined in subelements under the ``liquid_phase`` tag.  The ``liquid_phase`` tag takes an attribute ``name`` to identify which phase the observations are associated with.

The element name of individual observations indicate the quantity being observed.  Below is a list of currently available observations.  Individual observations require the subelements ``assigned_regions``, ``functional``, and ``time_macro``.  ``aqueous_conc`` observations also take an attribute ``name`` which indicates the name of the solute being observed.

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

An example ``observations`` element looks like the following.

.. code-block:: xml

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

.. _Akuna : http://esd.lbl.gov/research/projects/ascem/thrusts/platform/
.. _Mathematical Formulation Requirements and Specifications for the Process Models: http://software.lanl.gov/ascem/trac/attachment/wiki/Documents/ASCEM-HPC-ProcessModels_2011-01-0a.pdf

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
