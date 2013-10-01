============================================================
Amanzi Input File XML Schema 2.0
============================================================

Overview
++++++++

The Amanzi simulator evolves a system of conservation equations for
reacting flows in porous media, as detailed in the ASCEM report
entitled "Mathematical Formulation Requirements and Specifications for
the Process Models`".  The purpose of the present document is to specify
the data required to execute Amanzi.  This specification should be
regarded as a companion to the MRD, and parameterizations of the
individual submodels are consistent between Amanzi, the MRD and this
document. 

Models can be setup and valid, human-readable XML input files can be generated using Akuna.  Example input files are available in the Amanzi source repository.   

XML Schema 2.0
++++++++++++++


Model Description
-------------------

This section allows the user to provide information about the model being developed and how and when it was developed.  Default units for the model are also stored in this section.  This entire section is optional but encourage for documentation.

.. code-block:: xml

  <model_description name="Name of Model" >
      Required Elements: NONE
      Optional Elements: comment, author, created, modeified, model_id, description, purpose, units
  </model_description>

Units has the optional elements of length, time, mass, and concentration.  Each of those in turn have their own sturcture.  The structures are as follows.

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

Named_times allows the user to assign meaningful names to time values and define time values in a single location in the file.  Then the names are used throoughout the file whenever needed by boundary conditions or execution controls, etc.  The named_times element contains an unbounded number of time ``time`` elements. The trailing character in the value attribute indicates the units of the time.

.. code-block:: xml

  <named_times>
    <time  name="Name of Time" value="time,y|d|s">
  </named_times>

The ``macro`` section defines time and cycle macros.  These specifiy a series of times or cycles for writing out visualization or checkpoint files.  Each ``time_macros`` requires a ``name`` attribute and one or more ``time`` elements.  An alternative option is to specify the start and stop times and interval timestep, as shown in the ``cycle_macro``.

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

The ``execution_control`` section defines the general execution of the Amanzi simulation.  Amanzi can execute in three modes: steady state, transient, and initialize to a steady state and then continue it transient.  The execution mode, numerical method, and initial time step information can be specificied in this section for multiple time series.  If a series of times of defined, default values can be defined using the ``execution_control_defaults`` element.  Any undefined elements in subsequent ``execution_control`` elements will be filled in from the ``execution_control_defaults`` element.

.. code-block:: xml

  <execution_control_defaults init_dt="labeled_time" max_dt="labeled_time" reduction_factor="exponential" increase_factor="exponential" mode="steady | transient" method=" bdf1 | picard" />

.. code-block:: xml

  <execution_control  restart="string" start="string" end="string" init_dt="labeled_time" max_dt="labeled_time" reduction_factor="exponential" increase_factor="exponential" mode="steady | transient" method=" bdf1 | picard" />

The ``execution_control`` section also provides the elements ``comments`` and ``verbosity``.  Users may provide any text within the ``comment`` element to annotate this section.  ``verbosity`` takes the attribute level=`` high | medium | low``.  This triggers increasing levels of reporting from inside Amanzi

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

Mesh
----

A mesh must be defined for the simulation to be conducted on.  The mesh can be structured or unstructured.  It can also be an existing mesh from an ExodusII file or be generated internally by Amanzi.  All of this is specificied in te ``mesh`` section.  The type of mesh is specified through the attribute ``framework``.  Current options are ``mstk``, ``moab``, and ``exodus ii``.  As in other section a ``comments`` element is provide to include any comments or documentation the user wishes.  A ``dimension`` element speicifies where the mesh is 2D or 3D.  A 2D mesh can be given in 3D space with a third coordinate of 0.  Finally, either a ``file`` or ``gnereate`` element must be present.  A filename can be specified in the ``file`` element.  An orthogonal mesh can be generated by specifing the number of elements in each direction and the coordinates of the low and high corners of a box region.

Here is an overall example for the ``mesh`` element.

.. code-block:: xml

  <mesh framework="mstk"> 
    <comments>Pseudo 2D</comments>
    <dimension>3</dimension>
    <generate>
      <number_of_cells nx = "432"  ny = "1"  nz = "256"/>
      <box  low_coordinates = "0.0,0.0,0.0" high_coordinates = "216.0,1.0,107.52"/>
    </generate>
  </mesh>


Regions
-------

Geochemistry
------------

Material
--------

Process Kernels
---------------

Phases
------

Initial Conditions
------------------

Boundary Conditions
-------------------

Sources
-------

Outputs
-------
