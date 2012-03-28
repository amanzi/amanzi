========================================
Amanzi Native XML Input Specification V2
========================================

.. contents:: **Table of Contents**



ParameterList XML
=================

The Amanzi input file is an ASCII text XML-formatted file that must be framed at the beginning and end by the following statements:


.. code-block:: xml

  <ParameterList name="Main">

  </ParameterList>

The value in the "name" can be anything ("Main" in this example).  A ParameterList consists of just two types of entries: Parameter and ParameterList.  ParameterLists are labeled with a `"name`" [string], while Parameters have a separate fields for `"name`" [string], `"type`" [string] and `"value`" [TYPE], where "TYPE" can be any of the following: double, float, short, int, bool, string, Array double, Array float, Array short, Array int, Array bool, Array string.  The value of the parameter is given in quotes (e.g. "2.7e3").  Array data is specified as a single comma-deliminated string bounded by {}'s (e.g. "{2.4, 2.1, 5.7}").

.. code-block:: xml

  <ParameterList name="Sub">
    <Parameter name="CFL" type="double" value="0.9"/>
    <Parameter name="ratio" type="Array int" value="{2, 2, 4}"/>
  </ParameterList>

In this example, the sublist "Sub" has a parameter named "CFL" that is a "double" and has the value of 0.9, and a Teuchos::Array<int>
parameter named "ratio" such that ratio[0] = 2. ratio[1]=2 and ratio[2]=4.


Syntax of the Specification
===========================

* Input specification for each ParameterList entry consists of two parts.  First, a bulleted list defines the usage syntax and available options.  This is followed by example snipets of XML code to demonstrate usage.

* In many cases, the input specifies data for a particular parameterized model, and Amanzi supports a number of parameterizations.  For example, initial data might be uniform (the value is required), or linear in y (the value and its gradient are required).  Where Amanzi supports a number of parameterized models for quantity Z, the available models will be listed by name, and then will be described in the subsequent section.  For example, the specification might begin with the following:


 * `"X`" [list] 

  * `"Y`" [string]

  * Z [list] Model for Z, choose exactly one of the following: (1) `"Z: z1`", or (2) `"Z: z2`" (see below) 

Here, an `"X`" is defined by a `"Y`" and a `"Z`".  The `"Y`" is a string parameter but the `"Z`" is given by a model (which will require its own set of parameters).
The optoins for `"Z`" will then be described:

 * `"Z: z1`" applies model z1.  Requires `"z1a`" [string]

 * `"Z: z2`" applies model z2.  Requires `"z2a`" [double] and `"z2b`" [int]

An example of using such a specification:

.. code-block:: xml

    <ParameterList name="X">
      <Parameter name="Y" type="string" value="hello"/>
      <ParameterList name="Z: z2">
        <Parameter name="z2a" type="double" value="0.7"/>
        <Parameter name="z2b" type="int" value="3"/>
      </ParameterList>   
    </ParameterList>   
 
Here, the user is defining X with Y="hello", and Z will be a z2 constructed with z2a=0.7 and z2b=3.

Conventions:

* Reserved keywords and labels are `"quoted and italicized`" -- these labels or values of parameters in user-generated input files must match (using XML matching rules) the specified or allowable values.  User-defined labels are indicated with ALL-CAPS, and are meant to represent a typical name given by a user - these can be names or numbers or whatever serves best the organization of the user input data.

* Where applicable, the relevant section of the MRD is referred to by section or chapter number in parentheses.



MPC
===

State
=====

Flow
====

Flow sublist includes only one sublist, either `"Darcy Problem`" or `"Richars Problem`".
The second sublist contains all objects of the first one.

Sublist `"Water retention models`"

Transport
=========

Mesh
====

Regions
=======================================

Regions are geometrical constructs used in Amanzi to define subsets of the computational domain in order to specify the problem
to be solved, and the output desired.  Regions may represents zero-, one-, two- or three-dimensional subsets of physical space.
for a three-dimensional problem, the simulation domain will be a three-dimensional region bounded by a set of two-dimensional 
regions.  If the simulation domain is N-dimensional, the boundary conditions must be specified over a set of regions are (N-1)-dimensional.

Amanzi automatically defines the special region labeled `"All`", which is the 
entire simulation domain. Currently, the unstructured framework does
not support the `"All`" region, but it is expected to do so in the
near future.

Under the `"Structured`" option, Amanzi also automatically defines regions for the coordinat-aligned planes that bound the domain,
using the following labels: `"XLOBC`", `"XHIBC`", `"YLOBC`", `"YHIBC`", `"ZLOBC`", `"ZHIBC`"

User-defined regions are constructed using the following syntax

 * [U][S] "Regions" [list] can accept a number of lists for named regions (REGION)

   * Shape [list] Geometric model primitive, choose exactly one of the following [see table below]: `"Region: Point`", `"Region: Box`", `"Region: Plane`", `"Region: Labeled Set`", `"Region: Layer`", `"Region: Surface`"

Amanzi supports parameterized forms for a number of analytic shapes, as well as more complex definitions based on triangulated surface files.  

+--------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
|  shape functional name         | parameters                              | type(s)                      | Comment                                                                |
+================================+=========================================+==============================+========================================================================+
| `"Region: Point"`  [SU]        | `"Coordinate`"                          | Array double                 | Location of point in space                                             |
+--------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
| `"Region: Box"` [SU]           | `"Low Coordinate`", `"High Coordinate`" | Array double, Array double   | Location of boundary points of box                                     |
+--------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
| `"Region: Plane"`  [SU]        | `"Direction`", `"Location`"             | string, double               | direction: `"X`", `"-X`", etc, and `"Location`" is coordinate value    |
+--------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
| `"Region: Labeled Set"`        | `"Label`", `"File`",                    | string, string,              | Set per label defined in mesh file (see below)                         |
|                                | `"Format`", `"Entity`"                  | string, string               |  (available for frameworks supporting the `"File`" keyword)            |
+--------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
| `"Region: Color Function"` [S] | `"File`", `"Value`"                     | string, int                  | Set defined by color in a tabulated function file (see below)          |
+--------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
| `"Region: Layer"`              | `"File#`", `"Label#`"                   | (#=1,2) string, string       | Region between two surfaces                                            |
+--------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
| `"Region: Surface"`            | `"File`" `"Label`"                      | string, string               | Labeled triangulated face set in file                                  |
+--------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+

Notes

* `"Region: Point`" defines a point in space. Using this definition, cell sets encompassing this point are retrieved inside Amanzi.

* `"Region: Box`" defines a region bounded by coordinate-aligned
  planes. Boxes are allowed to be of zero thickness in only one
  direction in which case they are equivalent to planes.

* Currently, `"Region: Plane`" is constrained to be coordinate-aligned.

* The `"Region: Labeled Set`" region defines a named set of mesh entities
  existing in an input mesh file. This is the same file that contains
  the computational mesh. The name of the entity set is given
  by `"Label`".  For example, a mesh file in the Exodus II
  format can be processed to tag cells, faces and/or nodes with
  specific labels, using a variety of external tools.  Regions based
  on such sets are assigned a user-defined label for Amanzi, which may
  or may not correspond to the original label in the exodus file.
  Note that the file used to express this labeled set may be in any
  Amanzi-supported mesh format (the mesh format is specified in the
  parameters for this option).  The `"entity`" parameter may be
  necessary to specify a unique set.  For example, an Exodus file
  requires `"Cell`", `"Face`" or `"Node`" as well as a label (which is
  an integer).  The resulting region will have the dimensionality 
  associated with the entities in the indicated set. 

  By definition, "Labeled Set" region is applicable only to the
  unstructured version of Amanzi. 

  Currently, Amanzi-U only supports mesh files in the Exodus II format.

* `"Region: Color Function`" defines a region based a specified
  integer color, `"Value`", in a structured color function file,
  `"File`". The format of the color function file is given below in
  the "Tabulated function file format" section. As
  shown in the file, the color values may be specified at the nodes or
  cells of the color function grid. A computational cell is assigned
  the 'color' of the data grid cell containing its cell centroid
  (cell-based colors) or the data grid nearest its cell-centroid
  (node-based colors). Computational cells sets are then built from
  all cells with the specified color `"Value`".

  In order to avoid, gaps and overlaps in specifying materials, it is
  strongly recommended that regions be defined using a single color
  function file. 

* Surface files contain labeled triangulated face sets.  The user is
  responsible for ensuring that the intersections with other surfaces
  in the problem, including the boundaries, are `"exact`" (*i.e.* that
  surface intersections are `"watertight`" where applicable), and that
  the surfaces are contained within the computational domain.  If
  nodes in the surface fall outside the domain, the elements they
  define are ignored.

  Examples of surface files are given in the `"Exodus II`" file 
  format here.

* Region names must NOT be repeated

Example:

.. code-block:: xml

  <ParameterList name="Regions">
    <ParameterList name="Top Section">
      <ParameterList name="Region: Box">
        <Parameter name="Low Coordinate" type="Array double" value="{2, 3, 5}"/>
        <Parameter name="High Coordinate" type="Array double" value="{4, 5, 8}"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="Middle Section">
      <ParameterList name="Region: Box">
        <Parameter name="Low Coordinate" type="Array double" value="{2, 3, 3}"/>
        <Parameter name="High Coordinate" type="Array double" value="{4, 5, 5}"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="Bottom Section">
      <ParameterList name="Region: Box">
        <Parameter name="Low Coordinate" type="Array double" value="{2, 3, 0}"/>
        <Parameter name="High Coordinate" type="Array double" value="{4, 5, 3}"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="Inflow Surface">
      <ParameterList name="Region: Labeled Set">
        <Parameter name="Label"  type="string" value="sideset_2"/>
	<Parameter name="File"   type="string" value="F_area_mesh.exo"/>
	<Parameter name="Format" type="string" value="Exodus II"/>
	<Parameter name="Entity" type="string" value="Face"/>
      </ParameterList>
    </ParamterList>
    <ParameterList name="Outflow plane">
      <ParameterList name="Region: Plane">
        <Parameter name="Location" type="Array double" value="{0.5, 0.5, 0.5}"/>
        <Parameter name="Direction" type="Array double" value="{0, 0, 1}"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="Sand">
      <ParameterList name="Region: Color Function">
        <Parameter name="File" type="string" value="F_area_col.txt"/>
        <Parameter name="Value" type="int" value="25"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>

In this example, "Top Section", "Middle Section" and "Bottom Section"
are three box-shaped volumetric regions. "Inflow Surface" is a
surface region defined in an Exodus II-formatted labeled set
file and "Outflow plane" is a planar region. "Sand" is a volumetric
region defined by the value 25 in color function file.


Time Functions
==============

Boundary condition functions utilize a parameterized model for time variations that is either piecewise constant or piecewise linear.  For example:

.. code-block:: xml

      <Parameter name="Times" type="Array double" value="{1, 2, 3}"/>
      <Parameter name="Time Values" type="Array double" value="{10, 20, 30}"/>
      <Parameter name="Time Functions" type="Array string" value="{Constant, Linear}"/>    


This defines four time intervals: (-inf,1), (1,2), (2,3), (3,+inf).  By assumption the function is constant over the first and last intervals.  The remaining 
two intervals are speicified by the `"Time Functions`" parameter.  Thus, the value here is 10 anytime prior to t=2. The value increases linearly from 10 to 
20 over the interval t=2 to t=3, and then is constant at 30 for t>3.



Output
======

Output data from Amanzi is currently organized into four specific groups: `"Observations`", `"Visualization Data`", `"Checkpoint Data`" `"Diagnostic Output`" and `"Log Data`".  
Each of these is controlled in different ways, reflecting their intended use.

* `"Checkpoint Data`" is intended to represent all that is necesary to repeat or continue an Amanzi run.  The specific data contained in a Checkpoint Data dump is specific to the algorithm optoins and mesh framework selected.  Checkpoint Data is special in that no interpolation is perfomed prior to writing the data files; the raw binary state is necessary.  As a result, the user is allowed to only write Checkpoint Data at the discrete intervals of the simulation.

* `"Visualization Data`" is intended to represent spatially complete snapshots of the solution at defined instances during the simulation.  Dependeing on the control parameters provided here, visualizatoin files may include only a fraction of the state data, and may contiain auxiliary "derived" information (see below for more discussion).

* `"Observation Data`" is intended to represent diagnostic values to be returned to the calling routine from Amanzi's simulation driver.  Observations are typically generated at arbitrary times, and frequently involve various point samplings and volumetric reductions that are interpolated in time to the desired instant.  Observations may involve derived quantities (see discussion below) or state fields.

* `"Diagnostic Output`" is intended to represent diagnostic values to be written to stdout during a simulation. The available diagnostics are for the most part analogous to what is available as observations under the Observation Data capability. 

* `"Log Data`" is intended to represent runtime diagnostics to indicate the status of the simulation in progress.  This data is typically written by the simulation code to the screen or some other stream or file pipe.  The volume of `"Log Data`" generated is a function of the `"Verbosity`" setting under `"Execution Control`".

"`Log Data`" is not explicitly controlled in this section, since it is easier to control in the context of specifying details of the algorithms.  The remaining data types are discussed in the section below.


Time and Cycle specification
----------------------------

The user must specify when the various types of output are desired.  For Observation Data, this can be in terms of physical time.  For Visualization Data or Checkpoint Data, this can only be in terms of cycle number.  We support the definition of useful macros to specify these quantities.  One must specify the quantity over which these operators must function.  For example, you may want the integral of Moisture Content at various times as Observation Data, or the molar concentration of a solute at periodic cycles as Visualization Data.  The quantities must be identified from the standardized set:

* Available field quantities

 * Volumetric water content [volume water / bulk volume]
 * Aqueous saturation [volume water / volume pore space]
 * Aqueous pressure [Pa]
 * XXX Aqueous concentration [moles of solute XXX / volume water in MKS] (name formed by string concatenation, given the definitions in `"Phase Definition`" section)
 * X-, Y-, Z- Aqueous volumetric fluxe [m/s]
 * MaterialID

Note that MaterialID will be treated as a double that is unique to each defined material.  Its value will be generated internal to Amanzi.  The log file will be appended with the (material name)->(integer) mapping used.  Also note that this list tacitly assumes the presence of Aqueous Water as one of the transported components.  Presently, it is an error if the `"Phase Definition`" above does not sufficiently define this component.


Time macros specify a rule to generate a list of time values.  They are defined in the parameter list `"Time Macros`":

* [SU] `"Time Macros`" [list] can accept multiple lists for user-named macros TMACRO

 * [S] TMACRO [list] can accept either `"Values`" or `"Start_Period_Stop`"

  * [SU] `"Values`" [Array double] values of time, or 

  * [SU] `"Start_Period_Stop`" [Array double] values of start time (ts), period (dt) and (optionally) end time (te) to generate times, t=ts + dt*i, for any integer i. If stop time is less than start time, the time intervals have no ending.


Cycle macros specify a rule to generate or list cycle values.  They are defined in the parameter list `"Cycle Macros`":

* [SU] `"Cycle Macros`" [list] can accept multiple lists for user-named macros CMACRO

 * [SU] CMACRO [list] can accept either `"Values`" or `"Start_Period_Stop`"

  * [SU] `"Values`" [Array int] values of cycle number, or 

  * [SU] `"Start_Period_Stop`" [Array int] values of start cycle (cs), period (dc) and (optionally) end cycle (ce) to generate cycle numbers, c=cs + dc*i, for any integer i. If stop cycle < 0, the cycle intervals will not end.



Observation Data
----------------

A user may request any number of specific observations from Amanzi.  Each labeled Observation Data quantity involves a field quantity, a model, a region from which it will extract its source data, and a list of discrete times 
for its evaluation.  The observations are evaluated during the simulation and returned to the calling process through one of Amanzi arguments.

* [SU] `"Observation Data`" [list] can accept multiple lists for named observations (OBSERVATION)

  * [SU] `"Observation Output Filename`" [string] user-defined name for the file that the observations are written to.

  * [SU] OBSERVATION [list] user-defined label, can accept values for `"Variables`", `"Functional`", `"Region`", `"Time Macro`", and `"Cycle Macro`".

    * [SU] `"Variables`" [Array string] a list of field quantities taken from the list of "Available field quantities" defined above

    * [SU] `"Functional`" [string] the label of a function to apply to each of the variables in the variable list (Function options detailed below)

    * [SU] `"Region`" [string] the label of a user-defined region

    * [SU] `"Time Macro`" [string] one of the labeled time macros (see below)

    * [SU] `"Cycle Macro`" [string] one of the labeled time macros (see below)


The following Observation Data functionals are currently supported.  All of them operate on the variables identified.

* [SU] `"Observation Data: Point`" returns the value of the field quantity at a point

* `"Observation Data: Mean`" returns the mean value of the field quantities over the region specified

* [SU] `"Observation Data: Integral`" returns the integral of the field quantity over the region specified

* `"Observation Data: Cummulative Integral`" returns the integral of the field quantity, accumulated over the intervals defined by the time macro

* `"Observation Data: Peak Value`" returns the peak value of the field quantity over the region


Example:

.. code-block:: xml

  <ParameterList name="Time Macros">
    <ParameterList name="Annual">
      <Parameter name="Start_Period_Stop" type="Array double" value="{0, 3.1536e7,-1}"/>
    </ParameterList>
  </ParameterList>

  <ParameterList name="Observation Data">
    <Parameter name="Observation Output Filename" type="string" value="obs_output.out"/>
    <ParameterList name="Integrated Mass">
      <Parameter name="Region" type="string" value="All"/>
      <Parameter name="Functional" type="string" value="Observation Data: Integral"/>
      <Parameter name="Variables" type="Array string" value="{Volumetric Water Content, Tc-99 Aqueous Concentration}"/>
      <Parameter name="Time Macro" type="string" value="Annual"/>
    </ParameterList>
  </ParameterList>

In this example, the user requests an annual report of the integrated volume of water and aqueous solute concentration over the entire domain.


Diagnostic Output
-----------------

A user may request any number of specific observations from Amanzi that are 
written to stdout at times or cycles that are specified by the user.  Each 
labeled Diagnostic Output quantity involves a field quantity, a model, a 
region from which it will extract its source data, and a list of discrete 
times or cycles for its evaluation.  The diagnostics are evaluated during 
the simulation and written to stdout while the simulation is running.

* `"Diagnostic Output`" [list] can accept multiple lists for named Diagnostics (DIAGNOSTIC)

  * DIAGNOSTIC [list] user-defined label, can accept values for `"Variables`", `"Functional`", `"Region`", `"Time Macro`", and `"Cycle Macro`".

    * `"Functional`" [string] the label of a function to apply to each of the variables in the variable list (Function options detailed below)

    * `"Region`" [string] the label of a user-defined region

    * `"Time Macro`" [string] one of the labeled time macros (see below)

    * `"Cycle Macro`" [string] one of the labeled cycle macros (see below)


The following Observation Data functionals are currently supported.  All of them operate on the variables identified.

* `"Diagnostic Output: Point`" returns the value of the field quantity at the specified point

* `"Diagnostic Output: Mean`" returns the mean value of the field quantity 

* `"Diagnostic Output: Integral`" returns the integral of the field quantity 

* `"Diagnostic Output: Cummulative Integral`" returns the integral of the field quantity , accumulated over the intervals defined by the time macro

* `"Diagnostic Output: Peak Value`" returns the peak value of the field quantity


Example:

.. code-block:: xml

  <ParameterList name="Cycle Macros">
    <ParameterList name="first 100">
      <Parameter name="Start_Period_Stop" type="Array int" value="{0, 1, 99}"/>
    </ParameterList>
  </ParameterList>

  <ParameterList name="Diagnostic Output">
    <ParameterList name="User specified name of this diagnostic output">
      <Parameter name="Region" type="string" value="Some user specified point region"/>
      <Parameter name="Functional" type="string" value="Diagnostic Output: Point"/>
      <Parameter name="Variables" type="Array string" value="{Volumetric Water Content, Tc-99 Aqueous Concentration}"/>
      <Parameter name="Cycle Macro" type="string" value="first 100"/>
    </ParameterList>
  </ParameterList>



In this example the simulation will make point observations of the water volume and
concentration of Tc-99 in every one of the first 100 cycles and write the result
of these to stdout. 


Checkpoint Data
---------------------------------

A user may request periodic dumps of Amanzi Checkpoint Data.  The user has no explicit control over the content of these files, but has the guarantee that the Amanzi run will be reproducible (with accuracies determined
by machine round errors and randomness due to execution in a parallel computing environment).  Therefore, output controls for Checkpoint Data are limited to file name generation and writing frequency, by numerical cycle number.

* [SU] `"Checkpoint Data`" [list] can accept a file name base [string] and cycle data [list] used to generate the file base name or directory base name that is used in writing Checkpoint Data. 

  * [SU] `"File Name Base`" [string]

  * [SU] `"Cycle Macro`" [string] can accept label of user-defined Cycle Macro (see above)


Example:

.. code-block:: xml

  <ParameterList name="Cycle Macros">
    <ParameterList name="Every-5">
      <Parameter name="Start_Period" type="Array int" value="{0, 5}"/>
    </ParameterList>
  </ParameterList>

  <ParameterList name="Checkpoint Data">
    <Parameter name="File Name Base" type="string" value="chk"/>
    <Parameter name="File Name Digits" type="int" value="5"/>
    <Parameter name="Cycle Macro" type="string" value="Every-5"/>
  </ParameterList>

In this example, Checkpoint Data files are written when the cycle number is evenly divisible by 5.


Visualization Data
---------------------------------

A user may request periodic writes of field data for the purposes of visualization.  The user will specify explicitly what is to be included in the file at each snapshot.  Visualization files can only be written 
at intervals corresponding to the numerical time step values; writes are controlled by timestep cycle number.

* [SU] `"Visualization Data`" [list] can accept a file name base [string] and cycle data [list] that is used to generate the file base name or directory base name that is used in writing visualization data.  It can also accept a set of lists to specify which field quantities to write

  * [SU] `"File Name Base`" [string]
  
  * [SU] `"Cycle Macro`" [string] can accept label of user-defined Cycle Macro (see above)

  * [SU] `"Variables`" [string] can accept a list of field quantities to include in the file


Example:

.. code-block:: xml

  <ParameterList name="Cycle Macros">
    <ParameterList name="Every-10">
      <Parameter name="Start_Period" type="Array int" value="{0, 10}"/>
    </ParameterList>
  </ParameterList>

  <ParameterList name="Visualization Data">
    <Parameter name="File Name Base" type="string" value="chk"/>
    <Parameter name="File Name Digits" type="int" value="5"/>
    <Parameter name="Cycle Macro" type="string" value="Every-10"}>
    <Parameter name="Variable Macro" type="string" value="{Aqueous Pressure, Moisture Content"}>
  </ParameterList>

In this example, the liquid pressure and moisture content are written when the cycle number is evenly divisble by 5.


Restart from Checkpoint Data File
---------------------------------

A user may request a restart from a Checkpoint Data file by including the sublist 
`"Restart from Checkpoint Data File`" in the Execution Control list. This mode of restarting
will overwrite all other initializations of data that are called out in the input file.
The purpose of restarting Amanzi in this fashion is mostly to continue a run that has been 
terminated because its allocation of time ran out.


* [S] `"Restart from Checkpoint Data File`" [list]

  * [S] `"Checkpoint Data File Name`" [string] file name of the specific Checkpoint Data file to restart from

Example:

.. code-block:: xml

  <ParameterList name="Restart from Checkpoint Data File">
     <Parameter name="Checkpoint Data File Name" type="string" value="chk00123.h5"/>
  </ParameterList>

In this example, Amanzi is restarted with all state data initialized from the Checkpoint 
Data file named chk00123.h5. All other initialization of field variables that might be called 
out in the input file is ignored.  Recall that the value of "time" is taken from the checkpoint, 
but may be overridden by the execution control parameters.

Output format of Observation Output File
========================================
ASCII format will be used.   The file is preceded by two header lines:

* `Observation Name, Region, Functional, Variable, Time, Value`
* `======================================`

the first line describes what information are being displayed in entries in subsequent lines.  Each subsequent line
consists of 6 entries separated by the delimiter ",":

* Entry 1: `"ParameterList name`" for a particular observation output. 

* Entry 2: `"Region`" in the above `"ParameterList`".

* Entry 3: `"Functional`" in the above `"ParameterList`".

* Entry 4: `"Variable`" in the above `"ParameterList`".

* Entry 5: Time at which the observation is requested.

* Entry 6: Value of the observation at time specified in Entry 5.

An example is given by the following:

`Integrated Mass, All, Observation Data: Integral, Water Mass Density, 1000, 1.00e3`

Tabulated Function File Format
==============================

The following ASCII input file format supports the definition of a tabulated function defined over a grid.  Several XML input Parameters refer to files in this format.  The file consists of the following records (lines).  Each record is on a single line, except for the DATAVAL record which may be split across multiple lines.

1. **DATATYPE**:  An integer value: 0 for integer data, 1 for real data.

  * An integer-valued file is used to define a 'color' function used in the definition of a region.

2. **GRIDTYPE**:  A string that specifies the type of grid used to define the function.  The format of the rest of the file is contingent upon this value.  The currently supported options are uniform rectilinear grids in 1, 2 and 3-D, which are indicated by the values `1DCoRectMesh`, `2DCoRectMesh` and `3DCoRectMesh`, respectively (names adopted from XDMF).

For the uniform rectilinear grids, the remaining records are as follows.  Several records take 1, 2 or 3 values depending on the space dimension of the grid.

3. **NXNYNZ**: 3 (or 2, 1) integer values (NX, NY, NZ) giving the number of zones in the x, y and z coordinate directions, respectively.

4. **CORNER1**: 3 (or 2, 1) floating point values (X1, Y1, Z1) giving the coordinate of the first corner of the domain.

5. **CORNER2**: 3 (or 2, 1) floating point values (X2, Y2, Z2) giving the coordinate of the second corner of the domain.  The grid points r_{i,j,k} = (x_i, y_j, z_j) are defined as:

      x_i = X1 + i*(X2-X1)/NX, 0 <= i <= NX

      y_j = Y1 + j*(Y2-Y1)/NY, 0 <= j <= NY

      z_k = Z1 + k*(Z2-Z1)/NZ, 0 <= k <= NZ

  The (i,j,k) grid cell is defined by the corner grid points r_{i-1,j-1,k-1} and r_{i,j,k}, for 1 <= i <= NX, 1 <= j <= NY, 1 <= k <= NZ.  Note that the corner points are any pair of opposite corner points; the ordering of grid points and cells starts at CORNER1 and ends at CORNER2.

6. **DATALOC**:  An integer value: 0 for cell-based data, 1 for point-based data.


7. **DATACOL**:  An integer (N) giving the number of "columns" in the data.  This is the number of values per grid cell/point.  N=1 for a scalar valued function; N>1 for a N-vector valued function.

  * [U] only a single column is currently supported.

8. **DATAVAL**: The values of the function on the cells/points of the grid.  The values should appear in Fortran array order were the values stored in the Fortran array A(N,NX,NY,NZ) (A(N,0:NX,0:NY,0:NZ) for point-based data).  That is, the column index varies fastest, x grid index next fastest, etc.
    
Example
-------

As an example, consider the following integer-valued function in 2-D:

::
 
                  +-----+-----+-----+ (2.0,3.0)
                  |     |     |     |
                  |  2  |  1  |  1  |
                  |     |     |     |
                  +-----+-----+-----+
                  |     |     |     |
                  |  5  |  1  |  2  |
                  |     |     |     |
        (0.0,0.0) +-----+-----+-----+


The corresponding input file would be:

.. code-block:: text

  0
  2DCoRectMesh
  3 2
  0.0 0.0
  2.0 3.0
  0
  1
  5 1 2 2 1 1


