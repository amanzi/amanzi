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



MPC (tbw)
=========

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



State (tbw)
===========


Flow
====

Flow sublist includes exactly one sublist, either `"Darcy Problem`" or `"Richards Problem`".
Structure of both sublists is quite similar. We make necessary comments on differences.

Water retention models
-----------------------

User defines water retention models in sublist `"Water retention models`". It contains as many sublists, 
e.g. `"Model 1`", `"Model 2`", etc, as there are different soils. 
These models are associated with non-overlapping regions. Each of the sublists `"Model N`" 
inludes a few mandatory parameters: a region name, model name, and parameters for the selected model.
The available models are `"van Genuchten`", `"Brooks Corey`", and `"fake`". 
The later is used to set up an analytical solution for convergence study. 
The available models for the relative permeability are `"Mualem`" (default) and `"Burdine`".
An example of the van Genuchten model specification is:

.. code-block:: xml

    <ParameterList name="Model 1">
       <Parameter name="Region" type="string" value="Top Half"/>
       <Parameter name="Water retention model" type="string" value="van Genuchten"/>
       <Parameter name="van Genuchten alpha" type="double" value="0.000194"/>
       <Parameter name="van Genuchten m" type="double" value="0.28571"/>
       <Parameter name="van Genuchten l" type="double" value="0.5"/>
       <Parameter name="residual saturation" type="double" value="0.103"/>
       <Parameter name="relative permeability model" type="string" value="Mualem"/>
    </ParameterList>

    <ParameterList name="Model 2">
       <Parameter name="Region" type="string" value="Bottom Half"/>
       <Parameter name="Water retention model" type="string" value="Brooks Corey"/>
       <Parameter name="Brooks Corey lambda" type="double" value="0.0014"/>
       <Parameter name="Brooks Corey alpha" type="double" value="0.000194"/>
       <Parameter name="Brooks Corey l" type="double" value="0.51"/>
       <Parameter name="residual saturation" type="double" value="0.103"/>
       <Parameter name="regularization interval" type="double" value="0.0"/>
       <Parameter name="relative permeability model" type="string" value="Burdine"/>
    </ParameterList>


Amanzi performs rudimentary checks of validity of the provided parameters.


Boundary conditions
-------------------

Boundary conditions are defined in sublist `"boundary conditions`". Four types of boundary 
conditions are supported:

* `"pressure`" [list] Dirichlet boundary condition, a pressure is prescribed on a surface region. 

* `"mass flux`" [list] Neumann boundary condition, an outward mass flux is prescribed on a surface region.
  This is the default boundary condtion. If no condition is specified on a mesh face, zero flux 
  boundary condition is used implicitly.

* `"static head`" [list] Dirichlet boundary condition, the hydrostatic pressure is prescribed on a surface region.

* `"seepage face`" [list] Seepage face boundary condition, a dynamic combination of the `"pressure`" and 
  `"mass flux`" boundary conditions on a region. 
  The atmospheric pressure is prescribed if internal pressure is higher. Otherwise, the outward mass flux is prescribed. 

The following example includes all four types of boundary conditions. The boundary of a square domain 
is split into six pieces. Constant finction is used for simplicity and can be replaced by any
of the other available functions:

.. code-block:: xml

     <ParameterList name="boundary conditions">
       <ParameterList name="pressure">
         <ParameterList name="BC 0">
           <Parameter name="regions" type="Array string" value="{West side Top, East side Top}"/>
           <ParameterList name="boundary pressure">
             <ParameterList name="function-constant">
               <Parameter name="value" type="double" value="101325.0"/>
             </ParameterList>
           </ParameterList>
         </ParameterList>
       </ParameterList>

       <ParameterList name="mass flux">
         <ParameterList name="BC 1">
           <Parameter name="regions" type="Array string" value="{North side, South side}"/>
           <ParameterList name="outward mass flux">
             <ParameterList name="function-constant">
               <Parameter name="value" type="double" value="0.0"/>
             </ParameterList>
           </ParameterList>
         </ParameterList>
       </ParameterList>

       <ParameterList name="static head">
         <ParameterList name="BC 2">
           <Parameter name="regions" type="Array string" value="{West side Bottom}"/>
           <ParameterList name="water table elevation">
             <ParameterList name="function-constant">
               <Parameter name="value" type="double" value="10.0"/>
             </ParameterList>
           </ParameterList>
         </ParameterList>
       </ParameterList>

       <ParameterList name="seepage face">
         <ParameterList name="BC 3">
           <Parameter name="regions" type="Array string" value="{East side Bottom}"/>
           <ParameterList name="outward mass flux">
             <ParameterList name="function-constant">
               <Parameter name="value" type="double" value="1.0"/>
             </ParameterList>
           </ParameterList>
         </ParameterList>
       </ParameterList>
     </ParameterList>


Sources and Sinks
-----------------

The external sources (e.g. wells) are supported only in sublist `"Darcy Problems`". The structure
of sublist `"source terms`" follows the specification of boundary conditions. 
Again, constant functions can be replaced by any of the available time-functions:

.. code-block:: xml

     <ParameterList name="source terms">
       <ParameterList name="SRC 0">
         <Parameter name="regions" type="Array string" value="{Well east}"/>
         <ParameterList name="sink">
           <ParameterList name="function-constant">
             <Parameter name="value" type="double" value="-0.1"/>
           </ParameterList>
         </ParameterList>
       </ParameterList>

       <ParameterList name="SRC 1">
         <Parameter name="regions" type="Array string" value="{Well west}"/>
         <ParameterList name="sink">
           <ParameterList name="function-constant">
             <Parameter name="value" type="double" value="-0.2"/>
           </ParameterList>
         </ParameterList>
       </ParameterList>
     </ParameterList>

The remaining `"Flow`" parameters are

* `"atmospheric pressure`" [double] defines the atmosperic pressure, [Pa].

* `"relative permeability`" [string] defines a method for calculating relative
  permeability. The available self-explanatory options `"upwind with gravity`",
  are `"upwind with Darcy flux`", `"arithmetic mean`" and `"cell centered`". 
  The first three calculate the relative permeability on mesh interfaces.

* `"discretization method`" [string] helps to test new discretization methods. 
  The available options are `"generic mfd`", `"optimized mfd`", `"monotone mfd`", and
  `"support operator`". The last option reproduces discretization method implemented in RC1. 
  The third option is recommended for orthogonal meshes and diagonal absolute permeability.
  The second option is still experimental (no papers were published) and produces 
  an optimal discretization.

* `"VerboseObject`" [list] defines default verbosity level for the process kernel.
  If it does not exists, it will be created on a fly and verbosity level will be set to `"high`".
  Here is an example:

.. code-block:: xml

    <ParameterList name="VerboseObject">
      <Parameter name="Verbosity Level" type="string" value="medium"/>
    </ParameterList>


Steady State Time Integratior
-----------------------------

The sublist `"steady state time integrator`" defines parameters controling linear and 
nonlinear solvers during steady state time integration. Here is an example:

.. code-block:: xml

    <ParameterList name="steady state time integrator">
      <Parameter name="time integration method" type="string" value="BDF1"/>
      <Parameter name="initialize with darcy" type="string" value="yes"/>
      <Parameter name="clipping saturation value" type="double" value="0.98"/>
      <Parameter name="preconditoner" type="string" value="Trilinos ML">
      <Parameter name="linear solver" type="string" value="AztecOO GMRES">
 
      <ParameterList name="nonlinear solver BDF1">
      ...
      </ParameterList>
    </ParameterList>

The parameters used here are

* `"time integration method`" [string] defines a time integration method.
  The available options are `"BDF1`", `"BDF2`", and `"Picard`".

* `"initialize with darcy`" [string] solves the fully saturated problem with the 
  boundary continious avaluated at time T=0. The solution defines a new pressure
  and saturation.

* `"clipping saturation value`" [double] is an experimental option. It is used 
  after re-initialization described in the previous bullet to cut-off small values 
  of pressure. By default, the pressure threshold is equal to the atmospheric pressure.
  The new pressure is calculated based of the defined saturation value. Default is 0.6.

* `"preconditioner`" [string] refferes to a preconditioner sublist of the list `"Precondtioners`".

* `"linear solver`" [string] refferes to a solver sublist of the list `"Solvers`".


Transient Time Integratior
-----------------------------

The sublist `"transient time integrator`" defines parameters controling linear and 
nonlinear solvers during transient time integration. Its parameters are similar to 
that in the sublist `"steady state time integrator`" except for parameters controling
pressure re-initialization.


Transport
=========

The boundary conditions sublist mimics specification of the boundary conditions in `"Flow`".
Its structure is slightly simple, which is unnecessary and be replaced in the future. 
For the advective transport, the boundary conditions must be specified on inflow parts of the
boundary. If no value is prescribed through the XML input, the zero inlux boundary condition
is used. Note that the boundary condition is set up separately for each component:

.. code-block:: xml

    <ParameterList name="Transport BCs">
      <Parameter name="number of BCs" type="int" value="2"/>
      <ParameterList name="BC 0">
        <Parameter name="Component 0" type="Array double" value="{1.0, 1.0}"/>
        <Parameter name="Regions" type="Array string" value="{Left side}"/>
        <Parameter name="Time Functions" type="Array string" value="{Constant}"/>
        <Parameter name="Times" type="Array double" value="{0.0, 0.1}"/>
      </ParameterList>  

      <ParameterList name="BC 1">
        <Parameter name="Component 1" type="Array double" value="{1.0, 1.0}"/>
        <Parameter name="Regions" type="Array string" value="{Bottom side}"/>
        <Parameter name="Time Functions" type="Array string" value="{Constant}"/>
        <Parameter name="Times" type="Array double" value="{0.0, 0.1}"/>
      </ParameterList>  
    </ParameterList>  

The remaining `"Transport`" parameters are:

* `"CFL`" [double] time step limiter, a number less than 1 with default of 1.
   
* `"spatial discretization order`" [int] the order of the spatial discretization, either
  1 or 2. The default is 1. 
  
* `"temporal discretization order`" [int] the order of temporar discretization, either
  1 or 2. The default is 1.

* `"VerboseObject`" [list] defines default verbosity level for the process kernel.
  If it does not exists, it will be created on a fly and verbosity level will be set to `"high`".
  See an example under `"Flow`".


The `"Transport`" parameters useful for developers are:

* `"enable internal tests`" [string] various internal tests will be executed during
  the run time. The default value is `no`.
   
* `"internal tests tolerance`" [double] tolerance for internal tests such as the 
  divergence-free condition. The defult value is 1e-6.


Linear and Nonlinear Solvers
============================

Version 2 of the native input spec introduces this list.
At the moment it constans sublists for various linear an nonlinear solvers such as AztecOO.
Here is and example:

.. code-block:: xml

     <ParameterList name="Solvers">
       <ParameterList name="GMRES via AztecOO">
         <Parameter name="error tolerance" type="double" value="1e-12"/>
         <Parameter name="iterative method" type="string" value="GMRES"/>
         <Parameter name="maximal number of iterations" type="int" value="400"/>
       </ParameterList>
     </ParameterList>

The name `"GMRES via AztecOO`" is selected by the user.
It can be used by a process kernel lists to define a solver.


Preconditioners
===============

Version 2 of the native input spec introduces this list. It contains sublists for various
preconditioners required by a simulation. At the moment, we support only Trilinos multilevel 
preconditioner. Here is an example:

.. code-block:: xml

     <ParameterList name="Preconditoners">
       <ParameterList name="Trilinos ML">
          <ParameterList name="ML Parameters">
            <Parameter name="ML output" type="int" value="0"/>
            <Parameter name="aggregation: damping factor" type="double" value="1.33333"/>
            ... 
         </ParameterList>
       </ParameterList>

       <ParameterList name="Trilinos ML 2">
       ...
       </ParameterList>

       <ParameterList name="External AMG">
       ...
       </ParameterList>
     </ParameterList>

Names `"Trilinos ML`", `"Trilinos ML 2`", and `"External AMG`" are selected by the user.
They can be used by a process kernel lists to define a preconditioner.


Mesh
====

Amanzi supports both structured and unstructured numerical solution approaches.  This flexibility has a direct impact on the selection and design of the underlying numerical algorithms, the style of the software implementations, and, ultimately, the complexity of the user-interface.  "Mesh`" is used to select between the following options:

* `"Structured`": This instructs Amanzi to use BoxLib data structures and an associated paradigm to numerically represent the flow equations.  Data containers in the BoxLib software library, developed by CCSE at LBNL, are based on a hierarchical set of uniform Cartesian grid patches.  `"Structured`" requires that the simulation domain be a single coordinate-aligned rectangle, and that the "base mesh" consists of a logically rectangular set of uniform hexahedral cells.  This option supports a block-structured approach to dynamic mesh refinement, wherein successively refined subregions of the solution are constructed dynamically to track "interesting" features of the evolving solution.  The numerical solution approach implemented under the `"Structured`" framework is highly optimized to exploit regular data and access patterns on massively parallel computing architectures.

* `"Unstructured`": This instructs Amanzi to use data structures provided in the Trilinos software framework.  To the extent possible, the discretization algorithms implemented under this option are largely independent of the shape and connectivity of the underlying cells.  As a result, this option supports an arbitrarily complex computational mesh structure that enables users to work with numerical meshes that can be aligned with geometrically complex man-made or geostatigraphical features.  Under this option, the user typically provides a mesh file that was generated with an external software package.  The following mesh file formats are currently supported: `"Exodus 2`" (see example), `"MSTK`" (see example), `"MOAB`" (see example).  Amanzi also provides a rudmentary capability to generate unstructured meshes automatically.

Usage:

* [SU] `"Mesh`" [list] accepts either (1) `"Structured`", or (2) `"Unstructured`" to indicate the meshing option that Amanzi will use

 * [S] `"Structured`" [list] accepts coordinates defining the extents of simulation domain, and number of cells in each direction.

  * [S] `"Domain Low Coordinate`" [Array double] Location of low corner of domain

  * [S] `"Domain High Coordinate`" [Array double] Location of high corner of domain

  * [S] `"Number Of Cells`" [Array int] the number of uniform cells in each coordinate direction

 * [U] `"Unstructured`" [list] accepts instructions to either (1) read or, (2) generate an unstructured mesh.

  * [U] `"Read Mesh File`" [list] accepts name, format of pre-generated mesh file

   * [U] `"File`" [string] name of pre-generated mesh file. Note that in the case of an Exodus II mesh file, the suffix of the serial mesh file must be .exo. When running in serial the code will read this file directly. When running in parallel, the code will instead read the partitioned files, that have been generated with a Nemesis tool. There is no need to change the file name in this case as the code will automatically load the proper files. 

   * [U] `"Format`" [string] format of pre-generated mesh file (`"MSTK`", `"MOAB`", or `"Exodus II`")

  * [U] `"Generate Mesh`" [list] accepts parameters of generated mesh (currently only `"Uniform`" supported)

   * [U] `"Uniform Structured`" [list] accepts coordinates defining the extents of simulation domain, and number of cells in each direction.

    * [U] `"Domain Low Coordinate`" [Array double] Location of low corner of domain

    * [U] `"Domain High Coordinate`" [Array double] Location of high corner of domain

    * [U] `"Number Of Cells`" [Array int] the number of uniform cells in each coordinate direction

   * [U] `"Expert`" [list] accepts parameters that control which particular mesh framework is to be used.

    * [U] `"Framework`" [string] one of "stk::mesh", "MSTK",
      "MOAB" or "Simple". 
    * [U] `"Verify Mesh`" [bool] true or false. 


Example of `"Structured`" mesh:

.. code-block:: xml

   <ParameterList name="Mesh">
     <ParameterList name="Structured"/>
       <Parameter name="Number of Cells" type="Array int" value="{100, 1, 100}"/>
       <Parameter name="Domain Low Corner" type="Array double" value="{0.0, 0.0, 0.0}" />
       <Parameter name="Domain High Corner" type="Array double" value="{103.2, 1.0, 103.2}" />
     </ParameterList>   
   </ParameterList>

Example of `"Unstructured`" mesh generated internally:

.. code-block:: xml

   <ParameterList name="Mesh">
     <ParameterList name="Unstructured"/>
       <ParameterList name="Generate Mesh"/>
         <ParameterList name="Uniform Structured"/>
           <Parameter name="Number of Cells" type="Array int" value="{100, 1, 100}"/>
           <Parameter name="Domain Low Corner" type="Array double" value="{0.0, 0.0, 0.0}" />
           <Parameter name="Domain High Corner" type="Array double" value="{103.2, 1.0, 103.2}" />
         </ParameterList>   
       </ParameterList>   
     </ParameterList>   
   </ParameterList>

Example of `"Unstructured`" mesh read from an external file:

.. code-block:: xml

    <ParameterList name="Mesh">
      <ParameterList name="Unstructured">
        <ParameterList name="Read Mesh File">
          <Parameter name="File" type="string" value="mesh_filename"/>
          <Parameter name="Format" type="string" value="Exodus II"/>
        </ParameterList>   
      </ParameterList>   
    </ParameterList>


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



Time and Cycle specification (must be reviewed)
================================================

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



Observation Data (must be reviewed)
===================================

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


Checkpoint Data
===============

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
==================

A user may request periodic writes of field data for the purposes of visualization.  The user will specify explicitly what is to be included in the file at each snapshot.  Visualization files can only be written 
at intervals corresponding to the numerical time step values; writes are controlled by timestep cycle number.

* `"Visualization Data`" [list] can accept a file name base [string] and cycle data [list] that is used to generate the file base name or directory base name that is used in writing visualization data.  It can also accept a set of lists to specify which field quantities to write

  * `"File Name Base`" [string]
  
  * `"cycle start period stop`" [list] this is a list of start period stop definitions for cycles, each of which must be a sublist. Currently there can only be one sublist.

   * CSPS [list] can accept the only the parameter `"start period stop`".
    
    *  `"start period stop`" [Array int] the first entry is the start cycle, the second is the cycle period, and the third is the stop cycle or -1 in which case there is no stop cycle. A visualization dump shall be written for such cycles that satisfy cycle = start + n*period, for n=0,1,2,... and cycle < stop if stop != -1.

  * `"time start period stop`" [list] this is a list of start period stop definitions, each of which must be a sublist

   * TSPS [list] can accept only the parameter `"start period stop`".

    * `"start period stop`" [Array double] the first entry is the start time, the second is the time period, and the third is the stop time or -1 in which case there is no stop time. A visualization dump shall be written at such times that satisfy time = start + n*period, for n=0,1,2,... and time < stop if stop != -1.0.

  * `"times`" an array of discrete times that at which a visualization dump shall be written.

Example:

.. code-block:: xml

  <ParameterList name="Visualization Data">
    <Parameter name="File Name Base" type="string" value="chk"/>
  
    <ParameterList name="cycle start period stop">
      <ParameterList name="some unique name">
        <Parameter name="start period stop" type="Array int" value="{0, 100, -1}"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="time start period stop">
      <ParameterList name="some unique name">
        <Parameter name="start period stop" type="Array double" value="{0.0, 10.0, -1.0}"/>
      </ParameterList>
    </ParameterList>
    <Parameter name="times" type="Array double" value="{100.0, 300.0, 450.0}"/>
  </ParameterList>

In this example, the liquid pressure and moisture content are written when the cycle number is evenly divisble by 5.


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


