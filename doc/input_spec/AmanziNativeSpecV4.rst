========================================
Amanzi Native XML Input Specification V4
========================================

.. contents:: **Table of Contents**


Overview
========
This is a continuously evolving specification format used by
the code developers. It is used for development of new capabilities.


ParameterList XML
=================

The Amanzi input file is an ASCII text XML-formatted file that must be framed 
at the beginning and end by the following statements:


.. code-block:: xml

  <ParameterList name="Main">
    ...
  </ParameterList>

The value in the "name" can be anything ("Main" in this example).  
A ParameterList consists of just two types of entries: Parameter and ParameterList.  
ParameterLists are labeled with a `"name`" [string], while Parameters have a separate 
fields for `"name`" [string], `"type`" [string] and `"value`" [TYPE], where "TYPE" can 
be any of the following: double, float, short, int, bool, string, Array(double), Array(float), 
Array(short), Array(int), Array(bool), Array(string).  
The value of the parameter is given in quotes (e.g. "2.7e3").  
Array data is specified as a single comma-deliminated string bounded by {}'s (e.g. "{2.4, 2.1, 5.7}").

.. code-block:: xml

  <ParameterList name="Main">
    <Parameter name="CFL" type="double" value="0.9"/>
    <Parameter name="ratio" type="Array(int)" value="{2, 1, 4}"/>
  </ParameterList>

In this example, the sublist "Main" has a parameter named "CFL" that is a "double" and has 
the value of 0.9, and a Array<int> parameter named "ratio" such that ratio[0] = 2, 
ratio[1]=1, and ratio[2]=4.


Syntax of the Specification
===========================

* Input specification for each ParameterList entry consists of two parts.  
  First, a bulleted list defines the usage syntax and available options.  
  This is followed by example snippets of XML code to demonstrate usage.

* In many cases, the input specifies data for a particular parameterized model, and Amanzi 
  supports a number of parameterizations.  
  For example, initial data might be uniform (the value is required), or linear in y (the value 
  and its gradient are required).  
  Where Amanzi supports a number of parameterized models for quantity Z, the available 
  models will be listed by name, and then will be described in the subsequent section.  
  For example, the specification might begin with the following:


  * `"X`" [list] 

  * `"Y`" [string]

  * Z [list] Model for Z, choose exactly one of the following: (1) `"Z: z1`", or (2) `"Z: z2`" (see below) 

Let an `"X`" is defined by a `"Y`" and a `"Z`".  
The `"Y`" is a string parameter but the `"Z`" is given by a model (which will require its own set of parameters).
The options for `"Z`" will then be described:

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

The MPC stands for MultiProcess Coordinators.
In the MPC sublist the user specifies which process kernels are on or off, which 
flow model is active, and the time integration mode that the MPC should run in.

To turn a particular process kernel on or off use these options:

 * `"disable Transport_PK`" [string], valid options are `"yes`" or `"no`".

 * `"disable Flow_PK`" [string], valid options are `"yes`" or `"no`".

 * `"Chemistry Model`" [string], valid options are `"On`", `"Off`", or `"Alquimia`".

To select a particular flow model, use this option:

 * `"Flow model`" [string], valid options are `"Darcy`", `"Steady State Saturated`" 
   (both will cause the instantiation of a Darcy_PK process kernel), `"Richards`", 
   `"Steady State Richards`" (both will cause the instantiation of a Richards_PK 
   process kernel.

The following parameters control MPC options related to particular process kernels:

 * `"transport subcycling`" [bool], default is `"false`".

 * `"max chemistry to transport time step ratio`" [double], default is 1.0.

 * `"time integration rescue reduction factor`" [double], default is 0.5.

Time Integration Mode
---------------------

The MPC list must have a sublist `"Time Integration Mode`" if flow is enabled.
This sublist must have exactly one of the following three sublists: `"Steady`",
`"Initialize To Steady`", `"Transient with Static Flow`", or `"Transient`". 
The first one is used to find a steady-state solution and terminate the simulation. 
The steady-state calculations starts at time `"Start`" and terminates at time `"End`". 
The corresponding values as well as the initial time step are given in seconds.

.. code-block:: xml

  <ParameterList name="MPC">
    <ParameterList name="Time Integration Mode">
      <ParameterList name="Steady">
        <Parameter name="Start" type="double" value="0.0"/>
        <Parameter name="End" type="double" value="5.0"/>
        <Parameter name="Initial Time Step" type="double" value="0.1"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>

The second time integration mode is used to find a steady-state solution and 
continue with a transient simulation. The transition from a steady-state phase
to a transient phase happens at time `"Switch`".

.. code-block:: xml

  <ParameterList name="MPC">
    <ParameterList name="Time Integration Mode">
      <ParameterList name="Initialize To Steady">
        <Parameter name="Start" type="double" value="0.0"/>
        <Parameter name="Switch" type="double" value="0.5"/>
        <Parameter name="End" type="double" value="5.0"/>
        <Parameter name="Steady Initial Time Step" type="double" value="0.1"/>
        <Parameter name="Transient Initial Time Step" type="double" value="0.1"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>

The third time integration mode is used for a transient simulation only.

.. code-block:: xml

  <ParameterList name="MPC">
    <ParameterList name="Time Integration Mode">
      <ParameterList name="Transient">
        <Parameter name="Start" type="double" value="0.0"/>
        <Parameter name="End" type="double" value="5.0"/>
        <Parameter name="Initial Time Step" type="double" value="0.1"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>

The fourth time integration mode is used for a transient simulation only.
The flow field is static so no flow solver is called during time stepping. 
During initialization the flow field is set in one of two ways: 
(1) A constant Darcy velocity is specified in the initial condition; 
(2) Boundary conditions for the flow (e.g., pressure), along with the 
initial condition for the pressure field are used to solve for the Darcy velocity. 
At present this mode only supports the "Single Phase" flow model.

.. code-block:: xml

  <ParameterList name="MPC">
    <ParameterList name="Time Integration Mode">
      <ParameterList name="Transient With Static Flow">
        <Parameter name="Start" type="double" value="0.0"/>
        <Parameter name="End" type="double" value="5.0"/>
        <Parameter name="Initial Time Step" type="double" value="0.1"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>


Restart from Checkpoint Data File
---------------------------------

A user may request a restart from a Checkpoint Data file by including the MPC sublist 
`"Restart from Checkpoint Data File`". This mode of restarting
will overwrite all other initialization of data that are called out in the input file.
The purpose of restarting Amanzi in this fashion is mostly to continue a run that has been 
terminated because its allocation of time ran out.


* `"Restart from Checkpoint Data File`" [list]

  * `"Checkpoint Data File Name`" [string] file name of the specific Checkpoint Data file to restart from

  * `"initialize from checkpoint data file and do not restart`" [bool] (optional) If this is set to false (default), then a restart is performed, if it is set to true, then all fields are initialized from the checkpoint data file.

.. code-block:: xml
  
  <ParameterList name="MPC">
    ...
    <ParameterList name="Restart from Checkpoint Data File">
      <Parameter name="Checkpoint Data File Name" type="string" value="chk00123.h5"/>
    </ParameterList>
    ...
  </ParameterList>


In this example, Amanzi is restarted with all state data initialized from the Checkpoint 
Data file named chk00123.h5. All other initialization of field variables that might be called 
out in the input file is ignored.  Recall that the value for the current time and current cycle
is read from the checkpoint. 


Example for a complete MPC list
-------------------------------

The following is an example of a complete MPC list:

.. code-block:: xml

  <ParameterList name="MPC">
    <ParameterList name="Time Integration Mode">
      <ParameterList name="Initialize To Steady">
        <Parameter name="Start" type="double" value="0.0"/>
        <Parameter name="Switch" type="double" value="0.5"/>
        <Parameter name="End" type="double" value="5.0"/>
        <Parameter name="Steady Initial Time Step" type="double" value="0.1"/>
        <Parameter name="Transient Initial Time Step" type="double" value="0.1"/>
      </ParameterList>
    </ParameterList>
    <Parameter name="disable Transport_PK" type="string" value="yes"/>
    <Parameter name="Chemistry Model" type="string" value="Off"/>
    <Parameter name="disable Flow_PK" type="string" value="no"/>
    <Parameter name="Flow model" type="string" value="Steady State Saturated"/>
    <ParameterList name="Restart from Checkpoint Data File">
      <Parameter name="Checkpoint Data File Name" type="string" value="steady-checkpoint.h5"/>
    </ParameterList>
    <ParameterList name="VerboseObject">
      <Parameter name="Verbosity Level" type="string" value="high"/>
    </ParameterList>
  </ParameterList>

In this example, we included `"VerboseObject`" sublist. Its parameter `"Verbosity Level`"
controls the number of output messages. Note that higher verbosity levels come with
additional (but usually small) computational overhead. 
Such a sublist can be added safely to various sublists of an XML file.


State
=====

Sublist `"State`" allows the user to initialize physical fields using a variety of 
tools. 

.. code-block:: xml

  <ParameterList name="State">
    <ParameterList name="initial conditions">
       ... initialization of fields
    </ParameterList>
  </ParameterList>


Initialization of constant scalars
----------------------------------

A constant scalar field is the global (with respect to the mesh) constant. 
At the moment, the set of such fields includes fluid density 
and fluid viscosity.
The initialization requires to provide a named sublist with a single
parameter `"value`".

.. code-block:: xml

  <ParameterList name="fluid_density">
    <Parameter name="value" type="double" value="998.0"/>
  </ParameterList>


Initialization of constant vectors
----------------------------------

A constant vector field is the global (with respect to the mesh) vector constant. 
At the moment, the set of such vector constants includes gravity.
The initialization requires to provide a named sublist with a single
parameter `"Array(double)`". In two dimensions, is looks like

.. code-block:: xml

  <ParameterList name="gravity">
    <Parameter name="value" type="Array(double)" value="{0.0, -9.81}"/>
  </ParameterList>


Initialization of scalar fields
-------------------------------

A variable scalar field is defined by a few functions (labeled for instance,
`"MESH BLOCK i`" with non-overlapping ranges. 
The required parameters for each function are `"region`", `"component`",
and the function itself.

* `"regions`" [Array(string)] list of mesh regions where the function
  should be applied.

* `"component`" [string] specifies a mesh object on which the discrete field 
  is defined.

.. code-block:: xml

  <ParameterList name="porosity"> 
    <ParameterList name="function">
      <ParameterList name="MESH BLOCK 1">
        <Parameter name="regions" type="Array(string)" value="DOMAIN 1"/>
        <Parameter name="component" type="string" value="cell"/>
        <ParameterList name="function">
          <ParameterList name="function-constant">
            <Parameter name="value" type="double" value="0.2"/>
          </ParameterList>
        </ParameterList>
      </ParameterList>
      <ParameterList name="MESH BLOCK 2">
        ...
      </ParameterList>
    </ParameterList>
  </ParameterList>

In this example, the discrete field `"porosity`" has constant value 0.2 in 
each mesh cell of region `"DOMAIN ``". The second mesh block will define
porosity in other mesh regions.


Initialization of vector and tensor fields
------------------------------------------

A variable tensor (or vector) field is defined similarly to 
a variable scalar field. 
The difference lies in the definition of the function which
is now a multi-value function.
The required parameters are `"Number of DoFs`" and `"Function type`". 

* `"dot with normal`" [bool] triggers special initialization of a
  vector field such as the darcy flux. This field is defined by
  projection of a vector field on face normals.

.. code-block:: xml

  <ParameterList name="darcy_flux">
    <Parameter name="dot with normal" type="bool" value="true"/>
    <ParameterList name="function">
      <ParameterList name="MESH BLOCK 1">
        <Parameter name="regions" type="Array(string)" value="{ALL DOMAIN}"/>
        <Parameter name="component" type="string" value="face"/>
        <ParameterList name="function">
          <Parameter name="Number of DoFs" type="int" value="2"/>
          <Parameter name="Function type" type="string" value="composite function"/>
          <ParameterList name="DoF 1 Function">
            <ParameterList name="function-constant">
              <Parameter name="value" type="double" value="0.002"/>
            </ParameterList>
          </ParameterList>
          <ParameterList name="DoF 2 Function">
            <ParameterList name="function-constant">
              <Parameter name="value" type="double" value="0.001"/>
            </ParameterList>
          </ParameterList>
        </ParameterList>
      </ParameterList>
    </ParameterList>
  </ParameterList>

In this example the constant vector (0.002, 0.001) is dotted with the face 
normal producing one number per mesh face.
changing value of `"dot with normal`" to false will produce a vector 


Mesh Partitioning
-----------------

Amanzi's state has a number of tools to verify completeness of initial data.
This is done using sublist `"mesh partitions`". 
Each sublist in there must have parameter `"region list`" specifying
regions that define unique partition of the mesh.

.. code-block:: xml

    <ParameterList name="mesh partitions">
      <ParameterList name="MATERIALS">
        <Parameter name="region list" type="Array(string)" value="{region1, region2, region3}"/>
      </ParameterList>
    </ParameterList>

In this example, we verify that three mesh regions cover the mesh without overlaps.
If so, all material fields, e.g. porosity, will be initialized properly.


Initialization from a file
--------------------------

Some data can be initialized from files. Additional sublist has to be added to
named sublist of the `"State`" list with the file name and the name of an attribute. 
For a serial run, the file extension must be `".exo`". 
For a parallel run, it must be `".par`".

.. code-block:: xml

  <ParameterList name="permeability">
    <ParameterList name="exodus file initialization">
      <Parameter name="file" type="string" value="mesh_with_data.exo"/>
      <Parameter name="attribute" type="string" value="perm"/>
    </ParameterList>
  </ParameterList>


Example
-------

The complete example of a state initialization is below. Note that
`"MATERIAL 1`" and `"MATERIAL 2`" must be valid labels of regions.

.. code-block:: xml

  <ParameterList name="state">
    <ParameterList name="initial conditions">
      <ParameterList name="fluid_density">
        <Parameter name="value" type="double" value="998.0"/>
      </ParameterList>

      <ParameterList name="fluid_viscosity">
        <Parameter name="value" type="double" value="0.001"/>
      </ParameterList>

      <ParameterList name="gravity">
        <Parameter name="value" type="Array(double)" value="{0.0, -9.81}"/>
      </ParameterList>

      <ParameterList name="porosity"> <!-- pressure is done similarly -->
        <ParameterList name="function">
          <ParameterList name="domain">
            <Parameter name="regions" type="Array(string)" value="Computational domain"/>
            <Parameter name="component" type="string" value="cell"/>
            <ParameterList name="function">
              <ParameterList name="function-constant">
                <Parameter name="value" type="double" value="0.2"/>
              </ParameterList>
            </ParameterList>
          </ParameterList>
        </ParameterList>
      </ParameterList>

      <ParameterList name="permeability">
        <ParameterList name="function">
          <ParameterList name="MESH BLOCK 1">
            <Parameter name="regions" type="Array(string)" value="MATERIAL 1"/>
            <Parameter name="component" type="string" value="cell"/>
            <ParameterList name="function">
              <Parameter name="Function type" type="string" value="composite function"/>
              <Parameter name="Number of DoFs" type="int" value="2"/>
              <ParameterList name="DoF 1 Function">
                <ParameterList name="function-constant">
                  <Parameter name="value" type="double" value="1e-12"/>
                </ParameterList>
              </ParameterList>
              <ParameterList name="DoF 2 Function">
                <ParameterList name="function-constant">
                  <Parameter name="value" type="double" value="1e-13"/>
                </ParameterList>
              </ParameterList>
            </ParameterList>
          </ParameterList>
          <ParameterList name="MESH BLOCK 2">
            <Parameter name="regions" type="Array(string)" value="MATERIAL 2"/>
            <Parameter name="component" type="string" value="cell"/>
            <ParameterList name="function">
              <Parameter name="Function type" type="string" value="composite function"/>
              <Parameter name="Number of DoFs" type="int" value="2"/>
              <ParameterList name="DoF 1 Function">
                <ParameterList name="function-constant">
                  <Parameter name="value" type="double" value="2e-13"/>
                </ParameterList>
              </ParameterList>
              <ParameterList name="DoF 2 Function">
                <ParameterList name="function-constant">
                  <Parameter name="value" type="double" value="2e-14"/>
                </ParameterList>
              </ParameterList>
            </ParameterList>
          </ParameterList>
        </ParameterList>
      </ParameterList>
    </ParameterList>
  </ParameterList>


Flow
====

Flow sublist includes exactly one sublist, either `"Darcy Problem`" or `"Richards Problem`".
Structure of both sublists is quite similar. We make necessary comments on their differences.

Water retention models
-----------------------

User defines water retention models in sublist `"Water retention models`". 
It contains as many sublists, 
e.g. `"Soil 1`", `"Soil 2`", etc, as there are different soils. 
This list is required for `"Richards Problem`" only.
 
The water retention models are associated with non-overlapping regions. Each of the sublists (e.g. `"Soil 1`") 
includes a few mandatory parameters: region name, model name, and parameters for the selected model.
The available models are `"van Genuchten`", `"Brooks Corey`", and `"fake`". 
The later is used only to set up a simple analytic solution for convergence study. 
The available models for the relative permeability are `"Mualem`" (default) and `"Burdine`".

Amanzi performs rudimentary checks of validity of the provided parameters. 
The relative permeability curves can be calculated and saved in an ASCI file 
if the list `"Output`" is provided. This list has two mandatory parameters:

* `"file`" [string] A user defined file name. It should be different for 
  each soil. 

* `"number of points`" [int] A number of data points. 
  Each file will contain a table with three columns: saturation, relative permeability, and
  capillary pressure. The data points are equidistributed between the residual saturation
  and 1.

.. code-block:: xml

  <ParameterList name="Water retention models">
    <ParameterList name="Soil 1">
      <Parameter name="region" type="string" value="Top Half"/>
      <Parameter name="water retention model" type="string" value="van Genuchten"/>
      <Parameter name="van Genuchten alpha" type="double" value="0.000194"/>
      <Parameter name="van Genuchten m" type="double" value="0.28571"/>
      <Parameter name="van Genuchten l" type="double" value="0.5"/>
      <Parameter name="residual saturation" type="double" value="0.103"/>
      <Parameter name="relative permeability model" type="string" value="Mualem"/>
      <ParameterList name="Output"/>
        <Parameter name="file" type="string" value="soil1.txt"/>
        <Parameter name="number of points" type="int" value="1000"/>
      </ParameterList>
    </ParameterList>

    <ParameterList name="Soil 2">
      <Parameter name="region" type="string" value="Bottom Half"/>
      <Parameter name="water retention model" type="string" value="Brooks Corey"/>
      <Parameter name="Brooks Corey lambda" type="double" value="0.0014"/>
      <Parameter name="Brooks Corey alpha" type="double" value="0.000194"/>
      <Parameter name="Brooks Corey l" type="double" value="0.51"/>
      <Parameter name="residual saturation" type="double" value="0.103"/>
      <Parameter name="regularization interval" type="double" value="0.0"/>
      <Parameter name="relative permeability model" type="string" value="Burdine"/>
    </ParameterList>
  </ParameterList>

In this example, we define two different water retention models in two soils.


Diffusion operators
-------------------

Operators sublist describes the PDE structure of the flow, specifies a discretization
scheme, and selects assembling schemas for matrices and preconditioners.

* `"discretization primary`" [string] specifies an advanced discretization method that
  has useful properties under some a priori conditions on the mesh and/or permeability tensor.
  The available options are `"mfd scaled`", `"optimized mfd scaled`",
  `"finite volume`" and `"support operator`". 
  The last option reproduces discretization method implemented in RC1. 
  The second option is recommended for general meshes.
  The third option is recommended for orthogonal meshes and diagonal absolute 
  permeability tensor. 

* `"discretization secondary`" [string] specifies the most robust discretization method
  that is used when the primary selection fails to satisfy the a priori conditions.

* `"schema`" [Array(string)] defines the operator stencil. It is a collection of 
  geometric objects.

* `"preconditioner schema`" [Array(string)] defines the preconditioner stencil.
  It is needed only when the default assembling procedure is not desirable. If skipped,
  the `"schema`" is used instead. 

* `"gravity`" [bool] specifies if flow is driven also by the gravity.

* `"nonstandard symbolic assembling`" [int] specifies a nonstandard treatment of schemas.
  It is used for experiments with preconditioners.
  Default is 0.

.. code-block:: xml

    <ParameterList name="operators">
      <ParameterList name="diffusion operator">
        <ParameterList name="matrix">
          <Parameter name="discretization primary" type="string" value="monotone mfd"/>
          <Parameter name="discretization secondary" type="string" value="optimized mfd scaled"/>
          <Parameter name="schema" type="Array(string)" value="{face, cell}"/>
          <Parameter name="preconditioner schema" type="Array(string)" value="{face}"/>
          <Parameter name="gravity" type="bool" value="true"/>
        </ParameterList>
        <ParameterList name="preconditioner">
          <Parameter name="discretization primary" type="string" value="monotone mfd"/>
          <Parameter name="discretization secondary" type="string" value="optimized mfd scaled"/>
          <Parameter name="schema" type="Array(string)" value="{face, cell}"/>
          <Parameter name="preconditioner schema" type="Array(string)" value="{face}"/>
          <Parameter name="gravity" type="bool" value="true"/>
        </ParameterList>
      </ParameterList>
    </ParameterList>

This example creates a p-lambda system, i.e. the pressure is
discretized in mesh cells and on mesh faces. 
The preconditioner is defined on faces only, i.e. cell-based unknowns
are elliminated explicitly and the preconditioner is applied to the
Schur complement.


Boundary conditions
-------------------

Boundary conditions are defined in sublist `"boundary conditions`". Four types of boundary 
conditions are supported:

* `"pressure`" [list] Dirichlet boundary condition, a pressure is prescribed on a surface region. 

* `"mass flux`" [list] Neumann boundary condition, an outward mass flux is prescribed on a surface region.
  This is the default boundary condition. If no condition is specified on a mesh face, zero flux 
  boundary condition is used. 

* `"static head`" [list] Dirichlet boundary condition, the hydrostatic pressure is prescribed on a surface region.

* `"seepage face`" [list] Seepage face boundary condition, a dynamic combination of the `"pressure`" and 
  `"mass flux`" boundary conditions on a region. 
  The atmospheric pressure is prescribed if internal pressure is higher. Otherwise, the outward mass flux is prescribed. 

  * `"reference pressure`" [double] defaults to the atmospheric pressure. 

.. code-block:: xml

     <ParameterList name="boundary conditions">
       <ParameterList name="pressure">
         <ParameterList name="BC 0">
           <Parameter name="regions" type="Array(string)" value="{West side Top, East side Top}"/>
           <ParameterList name="boundary pressure">
             <ParameterList name="function-constant">
               <Parameter name="value" type="double" value="101325.0"/>
             </ParameterList>
           </ParameterList>
         </ParameterList>
       </ParameterList>

       <ParameterList name="mass flux">
         <ParameterList name="BC 1">
           <Parameter name="regions" type="Array(string)" value="{North side, South side}"/>
           <Parameter name="rainfall" type="bool" value="false"/>
           <ParameterList name="outward mass flux">
             <ParameterList name="function-constant">
               <Parameter name="value" type="double" value="0.0"/>
             </ParameterList>
           </ParameterList>
         </ParameterList>
       </ParameterList>

       <ParameterList name="static head">
         <ParameterList name="BC 2">
           <Parameter name="regions" type="Array(string)" value="{West side Bottom}"/>
           <Parameter name="relative to top" type="bool" value="true"/>
           <ParameterList name="water table elevation">
             <ParameterList name="function-constant">
               <Parameter name="value" type="double" value="10.0"/>
             </ParameterList>
           </ParameterList>
         </ParameterList>
       </ParameterList>

       <ParameterList name="seepage face">
         <Parameter name="reference pressure" type="double" value="101325.0"/>
         <ParameterList name="BC 3">
           <Parameter name="regions" type="Array(string)" value="{East side Bottom}"/>
           <ParameterList name="outward mass flux">
             <ParameterList name="function-constant">
               <Parameter name="value" type="double" value="1.0"/>
             </ParameterList>
           </ParameterList>
         </ParameterList>
       </ParameterList>
     </ParameterList>

This example includes all four types of boundary conditions. The boundary of a square domain 
is split into six pieces. Constant function is used for simplicity and can be replaced by any
of the other available functions.

The above boundary conditions are the four major models supported by Amanzi. In addition to
that each model may support a few submodels. A submodel is defined by additional
parameters described below. Mix and match of parameters is allowed.

* `"rainfall`" [bool] indicates that the mass flux is defined with respect to the gravity 
  vector and the actual influx depends on boundary slope. Default value is `"false`".

* `"relative to top`" [bool] indicates that the static head is defined with respect
  to the top boundary (a curve in 3D) of the specified regions. Support of 2D is turned off.
  Default value is `"false`". 

* `"no flow above water table`" [bool] indicates that no-flow (Neumann) boundary condition 
  has to be used above the water table. This switch uses the pressure value at a face
  centroid. Default is `"false`".

* `"submodel`" [string] indicates different models for the seepage face boundary condition.
  It can take values `"PFloTran`", `"FACT`", and `"Amanzi`". The first option leads to a 
  discontinuous change of the boundary condition type from the infiltration to pressure. 
  The second option is described
  in the document on mathematical models. It employs a smooth transition from the infiltration 
  to mixed boundary condition. The third option combines the above two. Is uses a smooth transition
  from the infiltration to pressure boundary condition. 
  Default value is `"Amanzi`".

.. code-block:: xml

       <ParameterList name="seepage face">
         <ParameterList name="BC 3">
           <Parameter name="regions" type="Array(string)" value="{California}"/>
           <Parameter name="rainfall" type="bool" value="true"/>
           <Parameter name="submodel" type="string" value="pflotran"/>
           <ParameterList name="outward mass flux">
             <ParameterList name="function-constant">
               <Parameter name="value" type="double" value="1.0"/>
             </ParameterList>
           </ParameterList>
         </ParameterList>
       </ParameterList>


Sources and Sinks
-----------------

The external sources are typically pumping wells. The structure
of sublist `"source terms`" mimics that of boundary conditions. 
Again, constant functions can be replaced by any of the available time-functions.

* `"spatial distribution method`" [string] identifies a method for distributing
  source Q over the specified regions. The available options are `"volume`",
  `"none`", and `"permeability`". For option `"none`" the source term Q is measured
  in [kg/m^3/s]. For the other options, it is measured in [kg/s]. When the source function
  is defined over a few regions, Q will be distributed independently over each region.
  Default is `"none`".

.. code-block:: xml

     <ParameterList name="source terms">
       <ParameterList name="SRC 0">
         <Parameter name="regions" type="Array(string)" value="{Well east}"/>
         <Parameter name="spatial distribution method" type="string" value="volume"/>
         <ParameterList name="sink">
           <ParameterList name="function-constant">
             <Parameter name="value" type="double" value="-0.1"/>
           </ParameterList>
         </ParameterList>
       </ParameterList>

       <ParameterList name="SRC 1">
         <Parameter name="regions" type="Array(string)" value="{Well west}"/>
         <Parameter name="spatial distribution method" type="string" value="permeability"/>
         <ParameterList name="sink">
           <ParameterList name="function-constant">
             <Parameter name="value" type="double" value="-0.2"/>
           </ParameterList>
         </ParameterList>
       </ParameterList>
     </ParameterList>


Initial Guess Pseudo Time Integrator
-------------------------------------

The sublist `"initial guess pseudo time integrator`" defines parameters controlling linear and 
nonlinear solvers during calculation of an initial guess.
This sublist will disappear as nonlinear solvers become more mature.
Detailed description of parameters is in the next two subsections.

.. code-block:: xml

   <ParameterList name="initial guess pseudo time integrator">
     <Parameter name="time integration method" type="string" value="Picard"/>
     <Parameter name="error control options" type="Array(string)" value="{pressure}"/>
     <Parameter name="linear solver" type="string" value="GMRES with TrilinosML"/>

     <ParameterList name="initialization">
       <Parameter name="method" type="string" value="saturated solver"/>
       <Parameter name="linear solver" type="string" value="CG with HypreAMG"/>
       <Parameter name="clipping saturation value" type="double" value="0.9"/>
     </ParameterList>

     <ParameterList name="pressure-lambda constraints">
       <Parameter name="method" type="string" value="projection"/>
       <Parameter name="inflow krel correction" type="bool" value="false"/>
       <Parameter name="linear solver" type="string" value="CG with HypreAMG"/>
     </ParameterList>

     <ParameterList name="Picard">
       <ParameterList name="Picard parameters">
         <Parameter name="convergence tolerance" type="double" value="1e-08"/>
         <Parameter name="maximum number of iterations" type="int" value="400"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>



Steady State Time Integrator
----------------------------

The sublist `"steady state time integrator`" defines parameters controlling linear and 
nonlinear solvers during steady state time integration. 
We break this long sublist into smaller parts. 
The first part controls preliminary steps in the time integrator.

* `"error control options`" [Array(string)] lists various error control options. 
  A nonlinear solver is terminated when all listed options are passed. 
  The available options are `"pressure`", `"saturation`", and `"residual`". 
  All errors are relative, i.e. dimensionless. 
  The error in pressure is compared with capillary pressure plus atmospheric pressure. 
  The other two errors are compared with 1. 
  The option `"pressure`" is always active during steady-state time integration.
  The option  `"saturation`" is always active during transient time integration.

* `"preconditioner`" [string] specifies preconditioner for nonlinear solvers.

* `"initialization`" [list] defines parameters for calculating initial pressure guess.
  It can be used to obtain pressure field which is consistent with the boundary conditions.
  Default is empty list.

  * `"method`" [string] is a placeholder for different initialization methods. Now, the only
    available option is `"saturated solver`" which lead to solving a Darcy problem.

  * `"linear solver`" [string] refers to a solver sublist of the list `"Solvers`".

  * `"clipping saturation value`" [double] is an experimental option. It is used 
    after pressure initialization to cut-off small values of pressure.
    The new pressure is calculated based of the provided saturation value. Default is 0.6.

  * `"clipping pressure value`" [double] is an experimental option. It is used 
    after pressure initialization to cut-off small values of pressure below the provided
    value.

* `"enforce pressure-lambda constraints`" [list] each time the time integrator is 
  restarted, we need to re-enforce the pressure-lambda relationship for new boundary conditions. 
  Default is empty list.

  * `"method`" [string] is a placeholder for different algorithms. Now, the only 
    available option is `"projection`" which is default.

  * `"linear solver`" [string] refers to a solver sublist of the list `"Solvers`".

  * `"inflow krel correction`" [bool] estimates relative permeability on inflow 
    mesh faces. This estimate is more reliable than the upwinded relative permeability
    value, especially in steady-state calculations.

.. code-block:: xml

   <ParameterList name="steady state time integrator">
     <Parameter name="error control options" type="Array(string)" value="{pressure, saturation}"/>
     <Parameter name="linear solver" type="string" value="GMRES with HypreAMG"/>
     <Parameter name="preconditioner" type="string" value="Hypre AMG"/>

     <ParameterList name="initialization">
       <Parameter name="method" type="string" value="saturated solver"/>
       <Parameter name="linear solver" type="string" value="CG with HypreAMG"/>
       <Parameter name="clipping pressure value" type="double" value="50000.0"/>
     </ParameterList>

     <ParameterList name="pressure-lambda constraints">
       <Parameter name="method" type="string" value="projection"/>
       <Parameter name="inflow krel correction" type="bool" value="false"/>
       <Parameter name="linear solver" type="string" value="CG with HypreAMG"/>
     </ParameterList>
   </ParameterList>

A specific time integration method is invoked by parameter `"time integration method`".
The available options are `"BDF1`" and `"Picard`".
The time step change is controlled by parameter `"time step controller type`".
Available options are `"fixed`", `"standard`", `"smarter`", and `"adaptive`".
The later is under development and is based on a posteriori error estimates.

* `"max preconditioner lag iterations`" [int] specifies frequency of 
  preconditioner recalculation.

* `"extrapolate initial guess`" [bool] identifies forward time extrapolation
  of the initial guess. Default is `"true`".

* `"restart tolerance relaxation factor`" [double] changes the nonlinear
  tolerance. The time integrator is usually restarted when a boundary condition 
  changes drastically. It may be beneficial to loosen the nonlinear 
  tolerance on the first several time steps after the time integrator restart. 
  The default value is 1, while reasonable values maybe as large as 1000. 

* `"restart tolerance relaxation factor damping`" controls how fast the loosened 
  nonlinear tolerance will revert back to the one specified in `"nonlinear tolerance"`.
  If the nonlinear tolerance is `"tol`", the relaxation factor is `"factor`", and 
  the damping is `"d`", and the time step count is `"n`" then the actual nonlinear 
  tolerance is `"tol * max(1.0, factor * d ** n)`".
  The default value is 1, while reasonable values are between 0 and 1.

* `"time step increase factor`" [double] defines geometric grow rate for the
  initial time step. This factor is applied when nonlinear solver converged
  in less than `"min iterations`" iterations. Default is 1.0.

* `"time step reduction factor`" [double] defines abrupt time step reduction
  when nonlinear solver failed or did not converge in  `"max iterations`" iterations.

* `"max time step`" [double] is the maximum allowed time step.

* `"min time step`" [double] is the minimum allowed time step.

.. code-block:: xml

   <ParameterList name="steady state time integrator">
     <Parameter name="max preconditioner lag iterations" type="int" value="5"/>
     <Parameter name="extrapolate initial guess" type="bool" value="true"/>
     <Parameter name="restart tolerance relaxation factor" type="double" value="1000.0"/>
     <Parameter name="restart tolerance relaxation factor damping" type="double" value="0.9"/>

     <Parameter name="time integration method" type="string" value="BDF1"/>
     <ParameterList name="BDF1">
       <Parameter name="timestep controller type" type="string" value="standard"/>
       <ParameterList name="timestep controller standard parameters">
         <Parameter name="min iterations" type="int" value="10"/>
         <Parameter name="max iterations" type="int" value="15"/>
         <Parameter name="time step increase factor" type="double" value="1.2"/>
         <Parameter name="time step reduction factor" type="double" value="0.5"/>
         <Parameter name="max time step" type="double" value="1e+9"/>
         <Parameter name="min time step" type="double" value="0.0"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>

In this example, the time step is increased by factor 1.2 when the nonlinear
solver converges in 10 or less iterations. 
The time step is not changed when the number of nonlinear iterations is
between 11 and 15.
The time step will be cut twice if the number of nonlinear iterations exceeds 15.

Amanzi supports a few nonlinear solvers described in details in a separate section.
Here, we recall parameters used in the NKA solver.

* `"solver type`" [string] defines nonlinear solver used on each time step for
  a nonlinear algebraic system :math:`F(x) = 0`. 
  The available options `"nka`" and `"Newton`".

* `"nka parameters`" [list] internal parameters for the nonlinear solver NKA.

  * `"nonlinear tolerance`" [double] is the convergence tolerance.

  * `"limit iterations`" [int] is the maximum allowed number of iterations.

  * `"diverged tolerance`" [double] is the maximum allowed error norm.

  * `"diverged l2 tolerance`" [double] is the maximum allowed relative L2 error norm.
    At the moment it is to prevent overflow only in the first NKA increment.

  * `"max du growth factor`" [double] limits the maximum change of the norm of
    the increment `du` during one nonlinear iteration step. 

  * `"max divergent iterations`" [int] limits the number of times the error
    can jump up during sequence of nonlinear iterations.

  * `"max nka vectors`" [int] is the size of the Krylov space.

  * `"modify correction`" [bool] allows to change (e.g. clip or damp) 
    the NKA or Newton correction. This is the experimental option with default `"false`".


.. code-block:: xml

     <ParameterList name="BDF1">
       <Parameter name="solver type" type="string" value="nka"/>
       <ParameterList name="nka parameters">
         <Parameter name="nonlinear tolerance" type="double" value="1e-5"/>
         <Parameter name="limit iterations" type="int" value="30"/>
         <Parameter name="diverged tolerance" type="double" value="1e+10"/>
         <Parameter name="diverged l2 tolerance" type="double" value="1e+5"/>
         <Parameter name="max du growth factor" type="double" value="1e+5"/>
         <Parameter name="max divergent iterations" type="int" value="3"/>
         <Parameter name="max nka vectors" type="int" value="10"/>
         <Parameter name="modify correction" type="bool" value="false"/>
         <ParameterList name="VerboseObject">
         <Parameter name="Verbosity Level" type="string" value="high"/>
         </ParameterList>
       </ParameterList>
     </ParameterList>

The remaining parameters in the time integrator sublist include 
those needed for unit tests, and future code development. 

.. code-block:: xml

   <ParameterList name="steady state time integrator">
     <ParameterList name="obsolete parameters">
       <Parameter name="start time" type="double" value="0.0"/>
       <Parameter name="end time" type="double" value="100.0"/>
       <Parameter name="maximum number of iterations" type="int" value="400"/>
       <Parameter name="error abs tol" type="double" value="1"/>
       <Parameter name="error rel tol" type="double" value="0"/>
     </ParameterList>
   </ParameterList>


Transient Time Integrator
-------------------------

The sublist `"transient time integrator`" defines parameters controlling linear and 
nonlinear solvers during transient time integration. Its parameters are similar to 
that in the sublist `"steady state time integrator`".
Here is a short example:
Note that the transient time integrator can be restarted multiple times, 
preferably every time a simulation goes through a stress test (e.g. external sourced 
are turning on and off abruptly).
If a non-empty `"initialization`" list is specified, it will be executed only once.

.. code-block:: xml

   <ParameterList name="transient time integrator">
     <Parameter name="time integration method" type="string" value="BDF1"/>
     <Parameter name="error control options" type="Array(string)" value="{pressure, saturation}"/>
     <Parameter name="linear solver" type="string" value="GMRES with HypreAMG"/>
     <Parameter name="time stepping strategy" type="string" value="adaptive"/>

     <ParameterList name="pressure-lambda constraints">
       <Parameter name="method" type="string" value="projection"/>
       <Parameter name="linear solver" type="string" value="CG with HypreAMG"/>
     </ParameterList>

     <ParameterList name="BDF1">
       ...
     </ParameterList>
   </ParameterList>



Other Parameters
----------------

The remaining `"Flow`" parameters are

* `"atmospheric pressure`" [double] defines the atmospheric pressure, [Pa].

* `"relative permeability`" [string] defines a method for calculating relative
  permeability. The available self-explanatory options `"upwind with gravity`",
  are `"upwind with Darcy flux`", `"arithmetic mean`" and `"cell centered`". 
  The first three calculate the relative permeability on mesh interfaces.

* `"upwind update`" [string] defines frequency of recalculating Darcy flux inside
  nonlinear solver. The available options are `"every time step`" and `"every nonlinear iteration`".
  The first option freezes the Darcy flux for the whole time step. The second option
  updates it on each iteration of a nonlinear solver. The second option is recommended
  for the New ton solver. It may impact significantly upwinding of the relative permeability 
  and convergence rate of this solver.

* `"clipping parameters`"[list] defines how corrections in nonlinear solver modified (clipped)

.. code-block:: xml

   <ParameterList name="clipping parameters">
      <Parameter name="max sat change" type="double" value="0.25"/>
      <Parameter name="damping factor" type="double" value="0.5"/>
   </ParameterList>	

* `"plot time history`" [bool] produces an ASCII file with time history when exists.

* `"VerboseObject`" [list] defines default verbosity level for the process kernel.
  If it does not exists, it will be created on a fly and verbosity level will be set to `"high`".

.. code-block:: xml

  <ParameterList name="Flow">
    <ParameterList name="Richards Problem">
      <Parameter name="atmospheric pressure" type="double" value="101325.0"/>
      <Parameter name="relative permeability" type="string" value="upwind with Darcy flux"/>
      <Parameter name="upwind update" type="string" value="every timestep"/>

      <ParameterList name="VerboseObject">
        <Parameter name="Verbosity Level" type="string" value="medium"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>



Transport
=========

The main parameters control temporal stability, spatial 
and temporal accuracy, and verbosity:

* `"CFL`" [double] time step limiter, a number less than 1 with default of 1.
   
* `"spatial discretization order`" [int] the order of the spatial discretization, either
  1 or 2. The default is 1. 
  
* `"temporal discretization order`" [int] the order of temporal discretization, either
  1 or 2. The default is 1.

* `"VerboseObject`" [list] defines default verbosity level for the process kernel.
  If it does not exists, it will be created on a fly and verbosity level will be set to `"high`".
  See an example under `"Flow`".

.. code-block:: xml

   <ParameterList name="Transport">
     <Parameter name="CFL" type="double" value="1.0"/>
     <Parameter name="spatial discretization order" type="int" value="1"/>
     <Parameter name="temporal discretization order" type="int" value="1"/>
     <Parameter name="advection limiter" type="string" value="Tensorial"/>

     <ParameterList name="VerboseObject">
       <Parameter name="Verbosity Level" type="string" value="high"/>
     </ParameterList>
   </ParameterList>  


Dispersivity models
-------------------
Two dispersivity models have been implemented: `"isotropic`" and `"Bear`". 
The anisotropic model `"Lichtner`" is pending for a more detailed 
description in the Process Models document.

.. code-block:: xml

   <ParameterList name="Dispersivity">
     <Parameter name="solver" type="string" value="Dispersive Solver"/>

     <ParameterList name="Brown Sugar">
       <Parameter name="regions" type="Array(string)" value="{top region, bottom region}"/>
       <Parameter name="model" type="string" value="Bear"/>
       <Parameter name="alphaL" type="double" value="1e-2"/>
       <Parameter name="alphaT" type="double" value="1e-5"/>
       <Parameter name="D" type="double" value="1e-8"/>
       <Parameter name="tortuosity" type="double" value="1e-4"/>       
     </ParameterList>  
     
     <ParameterList name="Grey Soil">
       <Parameter name="regions" type="Array(string)" value="{middle region}"/>
       <Parameter name="model" type="string" value="Bear"/>
       <Parameter name="alphaL" type="double" value="1e-2"/>
       <Parameter name="alphaT" type="double" value="1e-5"/>
       <Parameter name="D" type="double" value="1e-8"/>
       <Parameter name="tortuosity" type="double" value="1e-4"/>
     </ParameterList>  
   </ParameterList>  


Boundary Conditions
-------------------

For the advective transport, the boundary conditions must be specified on inflow parts of the
boundary. If no value is prescribed through the XML input, the zero influx boundary condition
is used. Note that the boundary condition is set up separately for each component.
The structure of boundary conditions is aligned with that used for Flow and
allows us to define spatially variable boundary conditions. 

.. code-block:: xml

   <ParameterList name="boundary conditions">
     <ParameterList name="concentration">
       <ParameterList name="H+"> 
         <ParameterList name="EAST CRIB">   <!-- user defined name -->
           <Parameter name="regions" type="Array(string)" value="{Top, Bottom}"/>
             <ParameterList name="boundary concentration">
               <ParameterList name="function-constant">  <!-- any time function -->
                 <Parameter name="value" type="double" value="0.0"/>
               </ParameterList>
             </ParameterList>
           </ParameterList>
         </ParameterList>
         <ParameterList name="WEST CRIB">   <!-- user defined name -->
           ...
         </ParameterList>
       </ParameterList>

       <ParameterList name="CO2"> <!-- Next component --> 
         ...
       </ParameterList>
     </ParameterList>
   </ParameterList>


Geochemical boundary conditions are concentration-type boundary conditions
but require special treatment. 

.. code-block:: xml

   <ParameterList name="boundary conditions">
     <ParameterList name="geochemical conditions">
       <ParameterList name="EAST CRIB">   <!-- user defined name -->
         <Parameter name="regions" type="Array(string)" value="{Crib1}"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>


Sources and Sinks
-----------------

The external sources are typically located at pumping wells. The structure
of sublist `"source terms`" includes only sublists named after components. 
Again, constant functions can be replaced by any available time-function.
Note that the source values are set up separately for each component.

* `"spatial distribution method`" [string] identifies a method for distributing
  source Q over the specified regions. The available options are `"volume`",
  `"none`", and `"permeability`". For option `"none`" the source term Q is measured
  in [mol/m^3/s]. For the other options, it is measured in [mol/s]. When the source function
  is defined over a few regions, Q will be distributed independently over each region.
  Default is `"none`".

.. code-block:: xml

     <ParameterList name="source terms">
       <ParameterList name="concentration">
         <ParameterList name="H+"> 
           <ParameterList name="source for east well">   <!-- user defined name -->
	     <Parameter name="regions" type="Array(string)" value="{Well east}"/>
             <Parameter name="spatial distribution method" type="string" value="volume"/>
             <ParameterList name="sink">   <!-- keyword, do not change -->
             <ParameterList name="function-constant">
               <Parameter name="value" type="double" value="-0.01"/>
             </ParameterList>
           </ParameterList>
           <ParameterList name="source for west well">
              ...
           </ParameterList>
         </ParameterList>
     
         <ParameterList name="Sugar syrop">   <!-- next component name -->
           <ParameterList name="source for west well">   <!-- user defined name -->
             <Parameter name="regions" type="Array(string)" value="{Well west}"/>
             <Parameter name="spatial distribution method" type="string" value="permeability"/>
             <ParameterList name="sink">  
               <ParameterList name="function-constant">
               <Parameter name="value" type="double" value="-0.02"/>
             </ParameterList>
           </ParameterList>
         </ParameterList>
       </ParameterList>
     </ParameterList>
    

Other parameters
-----------------

The other parameters useful for developers are:

* `"enable internal tests`" [string] various internal tests will be executed during
  the run time. The default value is `"no`".
   
* `"internal tests tolerance`" [double] tolerance for internal tests such as the 
  divergence-free condition. The default value is 1e-6.


Chemistry
=========

Geochemical engines
-------------------

This chemistry list specifies the default and the third-party geochemical engines. 
In the case of the third-party engine most details are provided in the trimmed 
PFloTran file `"1d-tritium-trim.in`".

The Alquimia chemistry process kernel only requires the `"Engine`" and `"Engine Input File`"
entries, but will also accept and respect the value given for `"Max Time Step (s)`". 
The rest are only used by the native chemistry kernel.

.. code-block:: xml

  <ParameterList name="Chemistry">
    <ParameterList name="Thermodynamic Database">
      <Parameter name="Format" type="string" value="simple"/>
      <Parameter name="File" type="string" value="tritium.bgd"/>
    </ParameterList>
    <Parameter name="Engine" type="string" value="PFloTran"/>
    <Parameter name="Engine Input File" type="string" value="1d-tritium-trim.in"/>
    <Parameter name="Verbosity" type="Array(string)" value="{verbose}"/>
    <Parameter name="Activity Model" type="string" value="unit"/>
    <Parameter name="Tolerance" type="double" value="1.5e-12"/>
    <Parameter name="Maximum Newton Iterations" type="int" value="25"/>
    <Parameter name="Max Time Step (s)" type="double" value="1.5778463e+07"/>
    <Parameter name="Number of component concentrations" type="int" value="1"/>
  </ParameterList>


Initial conditions
------------------

This sublist completes initialization of state variable, see list `"State`" for 
more detail. This section is only required for the native chemistry kernel, the
Alquimia chemistry kernel reads initial conditions from the `"State`" list.

.. code-block:: xml

    <ParameterList name="initial conditions">
      <ParameterList name="free_ion_species">
        <ParameterList name="function">
          <ParameterList name="ENTIRE DOMAIN">
            <Parameter name="region" type="string" value="Entire Domain"/>
            <Parameter name="component" type="string" value="cell"/>
            <ParameterList name="function">
              <Parameter name="Number of DoFs" type="int" value="1"/>
              <Parameter name="Function type" type="string" value="composite function"/>
              <ParameterList name="DoF 1 Function">
                <ParameterList name="function-constant">
                  <Parameter name="value" type="double" value="1.0e-09"/>
                </ParameterList>
              </ParameterList>
            </ParameterList>
          </ParameterList>
        </ParameterList>
      </ParameterList>
    </ParameterList>


Operators
=========

Operators are discrete forms of linearized PDEs operators.
They form a layer between physical process kernels and solvers
and include diffusion, advection, and source operators.
At the moment, a PK decides which collection of operators to be used to build 
a preconditioner.

* `"discretization primary`" [string] identifies a primary discretization method.
  Advanced discretization methods may have limitations due to mesh geometry and/or
  problem coefficients. In such a case the second discretization method is needed.

* `"discretization secondary`" [string] identifies a fallback discretization method.

* `"schema`" [Array(string)] defines the operator stencil. It is a collection of 
  geometric objects.

* `"preconditioner schema`" [Array(string)] defines the preconditioner stencil.
  It is needed only when the default assembling procedure is not desirable.

.. code-block:: xml

  <ParameterList name="PK operator">
    <Parameter name="preconditioner" type="string" value="Hypre AMG"/>

    <ParameterList name="diffusion operator">
      <Parameter name="discretization primary" type="string" value="mfd monotone"/>
      <Parameter name="discretization secondary" type="string" value="mfd optimized scaled"/>
      <Parameter name="schema" type="Array(string)" value="{cell,face}"/>
      <Parameter name="preconditioner schema" type="Array(string)" value="{cell,face}"/>
    </ParameterList>

    <ParameterList name="advection operator">
      <Parameter name="discretization primary" type="string" value="upwind"/>
      <Parameter name="reconstruction order" type="int" value="0"/>
    </ParameterList>
  </ParameterList>

In this example, the diffusion matrix is defined on mesh faces and cells.
The corresponding preconditioner is defined on the same objects. 
This discretization scheme corresponds to a p-lambda system.


Functions
=========

To set up non-trivial boundary conditions and/or initial fields, `Amanzi`
supports a few mathematical functions. 
New function types can added easily.
Each function is defined by a list:

.. code-block:: xml

  <ParameterList name="NAME">
    function-specification
  </ParameterList>

The parameter list name string `"NAME`" is arbitrary and meaningful only to the
parent parameter list.
This list is given as input to the Amanzi::FunctionFactory::Create
method which instantiates a new Amanzi::Function object.
The function-specification is one of the following parameter lists.


Constant function
-----------------
Constant function is defined as `f(x) = a`, for all `x`. 
The specification of this function needs only one parameter.
For example, when `a = 1`, we have:

.. code-block:: xml

  <ParameterList name="function-constant">
    <Parameter name="value" type="double" value="1.0"/>
  </ParameterList>
  

Tabular function
----------------
Given values :math:`x_i, y_i, i=0, ... n-1`, a tabular function :math:`f(x)` is 
defined piecewise: 

.. math::
  \begin{matrix}
  f(x) &=& x_0, & x \le x_0,\\
  f(x) &=& f(x_{i-1}) + (x - x_{i-1}) \frac{f(x_i) - f(x_{i-1})}{x_i - x_{i-1}},
  & x \in (x_{i-1}, x_i],\\
  f(x) &=& x_{n-1}, & x > x_{n-1}.
  \end{matrix}

This function is continuous and linear between two consecutive points.
This behavior can be changed using parameter `"forms`".
This parameter is optional.
If specified it must be an array of length equal to one less than the length 
of `x values`.  
Each value is either `"linear`" to indicate linear interpolation on that 
interval, or `"constant`" to use the left endpoint value for that interval.
The example defines function that is zero on interval :math:`(-\infty,\,0]`,
linear on interval :math:`(0,\,1]`, constant (`f(x)=1`) on interval :math:`(1,\,2]`,
and constant (`f(x)=2`) on interval :math:`(2,\,\infty]`.
The parameter `"x coordinate`" defines whether the `"x values`" refers to time `"t`",
x-coordinate `"x`", y-coordinate `"y`", or z-coordinate `"z`".
The default value for `"x coordinate`" is `"t`".

.. code-block:: xml

  <ParameterList name="function-tabular">
    <Parameter name="x values" type="Array(double)" value="{0.0, 1.0, 2.0}"/>
    <Parameter name="x coordinate" type="string" value="t"/>
    <Parameter name="y values" type="Array(double)" value="{0.0, 1.0, 2.0}"/>
    <Parameter name="forms" type="Array(string)" value="{linear, constant}"/>
  </ParameterList>
  

Bilinear function
-----------------
The bilinear function provides an extension of the linear form of the tabular function 
to a function with 2 variables `f(x,y)`.
A 2x2 matrix of values for `f(x,y)` and arrays of associated values for `x`
and `y` are read in from datasets in an HDF5 file. The dataset headers are indicated
by parameters `"row header`", `"column header`", and `"value header`" for `x`, `y`, 
and `f(x,y)`, respectively. The `x` and `y` arrays in the HDF5 file are expected to be
strictly increasing.
The parameters `"row coordinate`" and `"column coordinate`" define the model 
coordinate for `x` and `y` in the function, respectively, where
`"t`" refers to time, `"x`" to the x-coordinate, `"y`" to the y-coordinate, 
and `"z`" to the z-coordinate. 

The following code block defines a bilinear interpolation function for pressures
that vary in time and the x dimension.

.. code-block:: xml

  <ParameterList name="function-bilinear">
    <Parameter name="file" type="string" value="pressure_face.h5" />
    <Parameter name="row header" type="string" value="/time" />
    <Parameter name="row coordinate" type="string" value="time" />
    <Parameter name="column header" type="string" value="/x" />
    <Parameter name="column coordinate" type="string" value="x" />
    <Parameter name="value header" type="string" value="/pressures" />
  </ParameterList>
  

Smooth step function
--------------------
A smooth :math:`C^2` function `f(x)` on interval :math:`[x_0,\,x_1]` is 
defined such that `f(x) = y_0` for `x < x0`, `f(x) = y_1` for `x > x_1`, 
and monotonically increasing for :math:`x \in [x_0, x_1]`.
Here is an example:

.. code-block:: xml

  <ParameterList name="function-smooth-step">
    <Parameter name="x0" type="double" value="0.0"/>
    <Parameter name="y0" type="double" value="0.0"/>
    <Parameter name="x1" type="double" value="1.0"/>
    <Parameter name="y1" type="double" value="2.0"/>
  </ParameterList>


Polynomial function
-------------------
A generic polynomial function is given by the following expression:

.. math::
  f(x) = \sum_{j=0}^n c_j (x - x_0)^{p_j}

where :math:`c_j` are coefficients of monomials,
:math:`p_j` are integer exponents, and :math:`x_0` is the reference point.
Here i san example of a quartic polynomial:

.. code-block:: xml

  <ParameterList name="function-polynomial">
    <Parameter name="coefficients" type="Array(double)" value="{1.0, 1.0}"/>
    <Parameter name="exponents" type="Array(int)" value="{2, 4}"/>
    <Parameter name="reference point" type="double" value="0.0"/>
  </ParameterList>
  

Multi-variable linear function
------------------------------
A multi-variable linear function is formally defined by
 
.. math::
  f(x) = y_0 + \sum_{j=0}^{n-1} g_j (x_j - x_{0,j}) 

with the constant term "math:`y_0` and gradient :math:`g_0,\, g_1\,..., g_{n-1}`.
If the reference point :math:`x_0` is specified, it must have the same
number of values as the gradient.  Otherwise, it defaults to zero.
Note that one of the parameters in a multi-valued linear function can be time.
Here is an example:

.. code-block:: xml

  <ParameterList name="function-linear">
    <Parameter name="y0" type="double" value="1.0"/>
    <Parameter name="gradient" type="Array(double)" value="{1.0, 2.0, 3.0}"/>
    <Parameter name="x0" type="Array(double)" value="{2.0, 3.0, 1.0}"/>
  </ParameterList>
  

Separable function
------------------
A separable function is defined as the product of other functions such as

.. math::
  f(x_0, x_1,...,x_{n-1}) = f_1(x_0)\, f_2(x_1,...,x_{n-1})

where :math:`f_1` is defined by the `"function1`" sublist, and 
:math:`f_2` by the `"function2`" sublist:

.. code-block:: xml

  <ParameterList name="function-separable">
    <ParameterList name="function1">
      function-specification
    </ParameterList>
    <ParameterList name="function2">
      function-specification
    </ParameterList>
  </ParameterList>



Linear Solvers
==============

This list contains sublists for various linear solvers such as PCG, GMRES, and NKA.

* `"preconditioner`" [string] is name in the list of preconditioners. If it is missing, 
  the identity preconditioner is employed.

* `"iterative method`" [string] defines a Krylov-based method. The available options
  include `"pcg`" and `"gmres`".

* `"xxx parameters`" [sublist] provides parameters for the iterative method specified 
  by variable `"iterative method`".
 
.. code-block:: xml

     <ParameterList name="Solvers">
       <ParameterList name="GMRES with HYPRE AMG">
         <Parameter name="preconditioner" type="string" value="Hypre AMG"/>
         <Parameter name="iterative method" type="string" value="gmres"/>

         <ParameterList name="gmres parameters">
           ...
         </ParameterList>
       </ParameterList>

       <ParameterList name="PCG with HYPRE AMG">
         <Parameter name="preconditioner" type="string" value="Hypre AMG"/>
         <Parameter name="iterative method" type="string" value="pcg"/>
         <ParameterList name="pcg parameters">
           ...
         </ParameterList>
       </ParameterList>
     </ParameterList>

The names `"GMRES with HYPRE AMG`" and similar are chosen by the user.


GMRES
-----

Internal parameters for GMRES include

* `"error tolerance`" [double] is used in the convergence test. The default value is 1e-6.

* `"maximum number of iterations`" [int] is used in the convergence test. The default is 100.

* `"convergence criteria`" [Array(string)] specifies multiple convergence criteria. The list
  may include `"relative residual`", `"relative rhs`" (default), `"absolute residual`", and
  `"make one iteration`". The latter enforces the solver to perform at least one iteration
  which may be critical for extremely small time steps.

* `"size of Krylov space`" [int] defines the maximum size of the Krylov space. The default value is 10.

* `"overflow tolerance`" [double] defines the maximum allowed jump in residual. The default
  value is 3.0e+50.

.. code-block:: xml

    <ParameterList name="gmres parameters">
      <Parameter name="error tolerance" type="double" value="1e-12"/>
      <Parameter name="maximum number of iterations" type="int" value="400"/>
      <Parameter name="convergence criteria" type="Array(string)" value="{relative residual}"/>
      <Parameter name="size of Krylov space" type="int" value="10"/>
      <Parameter name="overflow tolerance" type="double" value="3.0e+50"/>

      <ParameterList name="VerboseObject">
        <Parameter name="Verbosity Level" type="string" value="high"/>
      </ParameterList>
    </ParameterList>


PCG
---

Internal parameters for PCG include

* `"error tolerance`" [double] is used in the convergence test. The default value is 1e-6.

* `"maximum number of iterations`" [int] is used in the convergence test. The default is 100.

* `"convergence criteria`" [Array(string)] specifies multiple convergence criteria. The list
  may include `"relative residual`", `"relative rhs`" (default), `"absolute residual`", and
  `"make one iteration`". The latter enforces the solver to perform at least one iteration
  which may be critical for extremely small time steps.

* `"overflow tolerance`" [double] defines the maximum allowed jump in residual. The default
  value is 3.0e+50.

.. code-block:: xml

    <ParameterList name="pcg parameters">
      <Parameter name="error tolerance" type="double" value="1e-12"/>
      <Parameter name="maximum number of iterations" type="int" value="400"/>
      <Parameter name="convergence criteria" type="Array(string)" value="{relative residual,make one iteration}"/>
      <Parameter name="overflow tolerance" type="double" value="3.0e+50"/>

      <ParameterList name="VerboseObject">
        <Parameter name="Verbosity Level" type="string" value="high"/>
      </ParameterList>
    </ParameterList>


NKA
---

This is a variation of the GMRES solver. Internal parameters for NKA include

* `"error tolerance`" [double] is used in the convergence test. The default value is 1e-6.

* `"maximum number of iterations`" [int] is used in the convergence test. The default is 100.

* `"convergence criteria`" [Array(string)] specifies multiple convergence criteria. The list
  may include `"relative residual`", `"relative rhs`" (default), `"absolute residual`", and
  `"make one iteration`". The latter enforces the solver to perform at least one iteration
  which may be critical for extremely small time steps.

* `"overflow tolerance`" [double] defines the maximum allowed jump in residual. The default
  value is 3.0e+50.

* `"max nka vectors`" [int] defines the maximum number of consecutive vectors used for 
  a local space.  The default value is 10.

* `"nka vector tolerance`" [int] defines the minimum allowed orthogonality between vectors in 
  the local space. If a new vector does not satisfy this requirement, the space is modified. 
  The default value is 0.05.

.. code-block:: xml

    <ParameterList name="nka parameters">
      <Parameter name="error tolerance" type="double" value="1e-12"/>
      <Parameter name="maximum number of iterations" type="int" value="400"/>
      <Parameter name="convergence criteria" type="Array(string)" value="{relative residual}"/>
      <Parameter name="overflow tolerance" type="double" value="3.0e+50"/>
      <Parameter name="max nka vectors" type="int" value="10"/>
      <Parameter name="nka vector tolerance" type="double" value="0.05"/>

      <ParameterList name="VerboseObject">
        <Parameter name="Verbosity Level" type="string" value="high"/>
      </ParameterList>
    </ParameterList>


Preconditioners
===============

This sublist contains entries for various
preconditioners required by a simulation. At the moment, we support Trilinos multilevel
preconditioner, Hypre BoomerAMG preconditioner, ILU preconditioner, Euclid ILU
preconditioner, and identity preconditioner. 

* `"type`" [string] defines preconditioner name.

* `"xxx parameters`" [sublist] provides parameters for the preconditioner specified 
  by variable `"type`".
 
.. code-block:: xml

     <ParameterList name="Preconditoners">
       <ParameterList name="TRILINOS ML">
          <Parameter name="type" type="string" value="ml"/>
          <ParameterList name="ml parameters">
            ... 
         </ParameterList>
       </ParameterList>

       <ParameterList name="HYPRE AMG">
          <Parameter name="type" type="string" value="boomer amg"/>
          <ParameterList name="boomer amg parameters">
            ...
          </ParameterList>
       </ParameterList>

       <ParameterList name="BLOCK ILU">
          <Parameter name="type" type="string" value="block ilu"/>
          <ParameterList name="block ilu parameters">
            ...
          </ParameterList>
       </ParameterList>
     </ParameterList>

Names `"TRILINOS ML`", `"HYPRE AMG`", and `"BLOCK ILU`" are choosen by the user.


Hypre AMG
---------

Internal parameters for Boomer AMG include

* `"tolerance`" [double] if is not zero, the preconditioner is dynamic 
  and approximate the inverse matrix with the prescribed tolerance (in
  the energy norm ???).

* `"relaxation type`" [int] defines the smoother to be used. Default is 6 
  which specifies a symmetric hybrid Gauss-Seidel / Jacobi hybrid method.

* `"smoother sweeps`" [int] defines the number of smoothing loops. Default is 3.

* `"cycle applications`" [int] defines the number of V-cycles. Default is 5.

* `"max multigrid levels`" [int] defined the maximum number of multigrid levels.

* `"number of function`" [int] the value > 1 tells Boomer AMG to use the "systems 
  of PDEs" code.  Note that, to use this approach, unknowns must be ordered with 
  DoF fastest varying (i.e. not the native Epetra_MultiVector order).  By default, it
  uses the "unknown" approach in which each equation is coarsened and
  interpolated independently.  Comments below are taken from Allison
  Baker's email to the PETSc mailing list, 25 Apr 2007, as these features
  of BoomerAMG are not documented very well.  Here we ignore her option
  2, as she warns it is inefficient and likely not useful.
  http://lists.mcs.anl.gov/pipermail/petsc-users/2007-April/001487.html

  * `"nodal strength of connection norm`" [int] tells AMG to coarsen such
    that each variable has the same coarse grid - sometimes this is more
    "physical" for a particular problem. The value chosen here for nodal
    determines how strength of connection is determined between the
    coupled system.  I suggest setting nodal = 1, which uses a Frobenius
    norm.  This does NOT tell AMG to use nodal relaxation.
    Default is 0.

* `"verbosity`" [int] prints a summary of run time settings and timing 
  information to stdout.  Default is 0.

.. code-block:: xml

   <ParameterList name="boomer amg parameters">
     <Parameter name="tolerance" type="double" value="0.0"/>
     <Parameter name="smoother sweeps" type="int" value="3"/>
     <Parameter name="cycle applications" type="int" value="5"/>
     <Parameter name="coarsen type" type="int" value="0"/>
     <Parameter name="strong threshold" type="double" value="0.5"/>
     <Parameter name="relaxation type" type="int" value="3"/>
     <Parameter name="verbosity" type="int" value="0"/>
   </ParameterList>


Euclid ILU
----------

The Euclid Parallel ILU algorithm was presented at SC99 and published in expanded 
form in the SIAM Journal on Scientific Computing. 
Scalability means that the factorization (setup) and application (triangular solve) timings remain
nearly constant when the global problem size is scaled in proportion to the number of processors.
As with all ILU preconditioning methods, the number of iterations is expected to increase with
global problem size.
Internal parameters for this preconditioner include

* `"ILU(k) fill level`" [int] is the factorization level. Default is 1.

* `"ILUT drop tolerance`" defines a drop tolerance relative to the largest 
  absolute value of any entry in the row being factored.

* `"rescale row`" [bool] if true, values are scaled prior to factorization 
  so that largest value in any row is +1 or -1. Note that this can destroy 
  matrix symmetry. 

* `"verbosity`" [int] prints a summary of runtime settings and timing 
  information to stdout.  Default is 0.

.. code-block:: xml

   <ParameterList name="euclid parameters">
     <Parameter name="ILU(k) fill level" type="int" value="6"/>
     <Parameter name="ILUT drop tolerance" type="double" value="0.01"/>
     <Parameter name="rescale rows" type="bool" value="true"/>
     <Parameter name="verbosity" type="int" value="0"/>
   </ParameterList>


Trilinos ML
-----------

Internal parameters for Trilinos ML include

.. code-block:: xml

   <ParameterList name="ml parameters">
     <Parameter name="ML output" type="int" value="0"/>
     <Parameter name="aggregation: damping factor" type="double" value="1.33"/>
     <Parameter name="aggregation: nodes per aggregate" type="int" value="3"/>
     <Parameter name="aggregation: threshold" type="double" value="0.0"/>
     <Parameter name="aggregation: type" type="string" value="Uncoupled"/>
     <Parameter name="coarse: type" type="string" value="Amesos-KLU"/>
     <Parameter name="coarse: max size" type="int" value="128"/>
     <Parameter name="coarse: damping factor" type="double" value="1.0"/>
     <Parameter name="cycle applications" type="int" value="2"/>
     <Parameter name="eigen-analysis: iterations" type="int" value="10"/>
     <Parameter name="eigen-analysis: type" type="string" value="cg"/>
     <Parameter name="max levels" type="int" value="40"/>
     <Parameter name="prec type" type="string" value="MGW"/>
     <Parameter name="smoother: damping factor" type="double" value="1.0"/>
     <Parameter name="smoother: pre or post" type="string" value="both"/>
     <Parameter name="smoother: sweeps" type="int" value="2"/>
     <Parameter name="smoother: type" type="string" value="Gauss-Seidel"/>
   </ParameterList>


Block ILU
---------

The internal parameters for block ILU are as follows:

.. code-block:: xml

   <ParameterList name="block ilu parameters">
     <Parameter name="fact: relax value" type="double" value="1.0"/>
     <Parameter name="fact: absolute threshold" type="double" value="0.0"/>
     <Parameter name="fact: relative threshold" type="double" value="1.0"/>
     <Parameter name="fact: level-of-fill" type="int" value="0"/>
     <Parameter name="overlap" type="int" value="0"/>
     <Parameter name="schwarz: combine mode" type="string" value="Add"/>
   </ParameterList>


Indentity
---------

The identity preconditioner is instantiated if either no preconditioner is
specified or the specified preconditioner list does not exists.


Mesh
====

Amanzi supports both structured and unstructured numerical solution approaches.
This flexibility has a direct impact on the selection and design of the underlying 
numerical algorithms, the style of the software implementations, and, ultimately, 
the complexity of the user-interface.  
`"Mesh`" is used to select between the following options:

* `"Structured`": This instructs Amanzi to use BoxLib data structures and an associated paradigm to numerically represent the flow equations.  Data containers in the BoxLib software library, developed by CCSE at LBNL, are based on a hierarchical set of uniform Cartesian grid patches.  `"Structured`" requires that the simulation domain be a single coordinate-aligned rectangle, and that the "base mesh" consists of a logically rectangular set of uniform hexahedral cells.  This option supports a block-structured approach to dynamic mesh refinement, wherein successively refined subregions of the solution are constructed dynamically to track "interesting" features of the evolving solution.  The numerical solution approach implemented under the `"Structured`" framework is highly optimized to exploit regular data and access patterns on massively parallel computing architectures.

* `"Unstructured`": This instructs Amanzi to use data structures provided in the Trilinos software framework.  To the extent possible, the discretization algorithms implemented under this option are largely independent of the shape and connectivity of the underlying cells.  As a result, this option supports an arbitrarily complex computational mesh structure that enables users to work with numerical meshes that can be aligned with geometrically complex man-made or geostatigraphical features.  Under this option, the user typically provides a mesh file that was generated with an external software package.  The following mesh file formats are currently supported: `"Exodus 2`" (see example), `"MSTK`" (see example), `"MOAB`" (see example).  Amanzi also provides a rudimentary capability to generate unstructured meshes automatically.

Usage:

* [SU] `"Mesh`" [list] accepts either (1) `"Structured`", or (2) `"Unstructured`" to indicate the meshing option that Amanzi will use

 * [S] `"Structured`" [list] accepts coordinates defining the extents of simulation domain, and number of cells in each direction.

  * [S] `"Domain Low Coordinate`" [Array(double)] Location of low corner of domain

  * [S] `"Domain High Coordinate`" [Array(double)] Location of high corner of domain

  * [S] `"Number Of Cells`" [Array(int)] the number of uniform cells in each coordinate direction

 * [U] `"Unstructured`" [list] accepts instructions to either (1) read or, (2) generate an unstructured mesh.

  * [U] `"Read Mesh File`" [list] accepts name, format of pre-generated mesh file

   * [U] `"File`" [string] name of pre-generated mesh file. Note that in the case of an Exodus II mesh file, the suffix of the serial mesh file must be .exo. When running in serial the code will read this file directly. When running in parallel, the code will instead read the partitioned files, that have been generated with a Nemesis tool. There is no need to change the file name in this case as the code will automatically load the proper files. 

   * [U] `"Format`" [string] format of pre-generated mesh file (`"MSTK`", `"MOAB`", or `"Exodus II`")

  * [U] `"Generate Mesh`" [list] accepts parameters of generated mesh (currently only `"Uniform`" supported)

   * [U] `"Uniform Structured`" [list] accepts coordinates defining the extents of simulation domain, and number of cells in each direction.

    * [U] `"Domain Low Coordinate`" [Array(double)] Location of low corner of domain

    * [U] `"Domain High Coordinate`" [Array(double)] Location of high corner of domain

    * [U] `"Number Of Cells`" [Array(int)] the number of uniform cells in each coordinate direction

   * [U] `"Expert`" [list] accepts parameters that control which particular mesh framework is to be used.

    * [U] `"Framework`" [string] one of "stk::mesh", "MSTK",
      "MOAB" or "Simple". 
    * [U] `"Verify Mesh`" [bool] true or false. 


Example of `"Structured`" mesh:

.. code-block:: xml

   <ParameterList name="Mesh">
     <ParameterList name="Structured"/>
       <Parameter name="Number of Cells" type="Array(int)" value="{100, 1, 100}"/>
       <Parameter name="Domain Low Corner" type="Array(double)" value="{0.0, 0.0, 0.0}" />
       <Parameter name="Domain High Corner" type="Array(double)" value="{103.2, 1.0, 103.2}" />
     </ParameterList>   
   </ParameterList>

Example of `"Unstructured`" mesh generated internally:

.. code-block:: xml

   <ParameterList name="Mesh">
     <ParameterList name="Unstructured"/>
       <ParameterList name="Generate Mesh"/>
         <ParameterList name="Uniform Structured"/>
           <Parameter name="Number of Cells" type="Array(int)" value="{100, 1, 100}"/>
           <Parameter name="Domain Low Corner" type="Array(double)" value="{0.0, 0.0, 0.0}" />
           <Parameter name="Domain High Corner" type="Array(double)" value="{103.2, 1.0, 103.2}" />
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
=======

Regions are geometrical constructs used in Amanzi to define subsets of the computational domain in order to specify the problem
to be solved, and the output desired.  Regions may represents zero-, one-, two- or three-dimensional subsets of physical space.
for a three-dimensional problem, the simulation domain will be a three-dimensional region bounded by a set of two-dimensional 
regions.  If the simulation domain is N-dimensional, the boundary conditions must be specified over a set of regions are (N-1)-dimensional.

Amanzi automatically defines the special region labeled `"All`", which is the 
entire simulation domain. Currently, the unstructured framework does
not support the `"All`" region, but it is expected to do so in the
near future.

Under the `"Structured`" option, Amanzi also automatically defines regions for the coordinate-aligned planes that bound the domain,
using the following labels: `"XLOBC`", `"XHIBC`", `"YLOBC`", `"YHIBC`", `"ZLOBC`", `"ZHIBC`"

User-defined regions are constructed using the following syntax

 * [U][S] "Regions" [list] can accept a number of lists for named regions (REGION)

   * Shape [list] Geometric model primitive, choose exactly one of the following [see table below]: `"Region: Point`", `"Region: Box`", `"Region: Plane`", `"Region: Labeled Set`", `"Region: Layer`", `"Region: Surface`"

Amanzi supports parameterized forms for a number of analytic shapes, as well as more complex definitions based on triangulated surface files.  

+--------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
|  shape functional name         | parameters                              | type(s)                      | Comment                                                                |
+================================+=========================================+==============================+========================================================================+
| `"Region: Point"`  [SU]        | `"Coordinate`"                          | Array(double)                | Location of point in space                                             |
+--------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
| `"Region: Box"` [SU]           | `"Low Coordinate`", `"High Coordinate`" | Array(double), Array(double) | Location of boundary points of box                                     |
+--------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
| `"Region: Plane"`  [SU]        | `"Direction`", `"Location`"             | string, double               | direction: `"X`", `"-X`", etc, and `"Location`" is coordinate value    |
+--------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
| `"Region: Polygon"`  [U]       | `"Number of points`", `"Points`"        | int, Array double            | Number of polygon points and point coordinates in linear array         |
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

* `"Region: Polygon`" defines a polygonal region on which mesh faces and
  nodes can be queried. NOTE that one cannot ask for cells in a polygonal
  region.In 2D, the "polygonal" region is a line and is specified by 2 points
  In 3D, the "polygonal" region is specified by an arbitrary number of points.
  In both cases the point coordinates are given as a linear array. The polygon
  can be non-convex.

  The polygonal region can be queried for a normal. In 2D, the normal is
  defined as [Vy,-Vx] where [Vx,Vy] is the vector from point 1 to point 2.
  In 3D, the normal of the polygon is defined by the order in which points 
  are specified.

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

.. code-block:: xml

  <ParameterList name="Regions">
    <ParameterList name="Top Section">
      <ParameterList name="Region: Box">
        <Parameter name="Low Coordinate" type="Array(double)" value="{2, 3, 5}"/>
        <Parameter name="High Coordinate" type="Array(double)" value="{4, 5, 8}"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="Middle Section">
      <ParameterList name="Region: Box">
        <Parameter name="Low Coordinate" type="Array(double)" value="{2, 3, 3}"/>
        <Parameter name="High Coordinate" type="Array(double)" value="{4, 5, 5}"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="Bottom Section">
      <ParameterList name="Region: Box">
        <Parameter name="Low Coordinate" type="Array(double)" value="{2, 3, 0}"/>
        <Parameter name="High Coordinate" type="Array(double)" value="{4, 5, 3}"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="Inflow Surface">
      <ParameterList name="Region: Labeled Set">
        <Parameter name="Label"  type="string" value="sideset_2"/>
	<Parameter name="File"   type="string" value="F_area_mesh.exo"/>
	<Parameter name="Format" type="string" value="Exodus II"/>
	<Parameter name="Entity" type="string" value="Face"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="Outflow plane">
      <ParameterList name="Region: Plane">
        <Parameter name="Location" type="Array(double)" value="{0.5, 0.5, 0.5}"/>
        <Parameter name="Direction" type="Array(double)" value="{0, 0, 1}"/>
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


Output
======

VerboseObject
-------------

Output of all components of Amanzi is controlled by a standard verbose 
object sublist. If this list is not specified, the default verbosity
value is used.

.. code-block:: xml

    <ParameterList name="VerboseObject">
      <Parameter name="Verbosity Level" type="string" value="high"/>
    </ParameterList>



Time Functions
--------------

Boundary condition functions utilize a parameterized model for time variations that is either piecewise constant or piecewise linear.  For example:

.. code-block:: xml

      <Parameter name="Times" type="Array(double)" value="{1, 2, 3}"/>
      <Parameter name="Time Values" type="Array(double)" value="{10, 20, 30}"/>
      <Parameter name="Time Functions" type="Array(string)" value="{Constant, Linear}"/>    


This defines four time intervals: (-inf,1), (1,2), (2,3), (3,+inf).  By assumption the function is constant over the first and last intervals.  The remaining 
two intervals are specified by the `"Time Functions`" parameter.  Thus, the value here is 10 anytime prior to t=2. The value increases linearly from 10 to 
20 over the interval t=2 to t=3, and then is constant at 30 for t>3.


Observation Data
----------------

A user may request any number of specific observations from Amanzi.  Each labeled Observation Data quantity involves a field quantity, a model, a region from which it will extract its source data, and a list of discrete times 
for its evaluation.  The observations are evaluated during the simulation and returned to the calling process through one of Amanzi arguments.

* `"Observation Data`" [list] can accept multiple lists for named observations (OBSERVATION)

  * `"Observation Output Filename`" [string] user-defined name for the file that the observations are written to.

  * OBSERVATION [list] user-defined label, can accept values for `"variables`", `"functional`", `"region`", `"times`", and TSPS (see below).

    * `"variables`" [Array(string)] a list of field quantities taken from the list of 
      available field quantities:

      * Volumetric water content [volume water / bulk volume]
      * Aqueous saturation [volume water / volume pore space]
      * Aqueous pressure [Pa]
      * Hydraulic Head [m] 
      * Drawdown [m] 
      * XXX Aqueous concentration [moles of solute XXX / volume water in MKS] (name formed by string concatenation, given the definitions in `"Phase Definition`" section)
      * X-, Y-, Z- Aqueous volumetric flux [m/s]
      * MaterialID

    Observation "Drawdown" is calculated with respect to the value registered at the first time
    it was requested.

    * `"functional`" [string] the label of a function to apply to each of the variables in the variable list (Function options detailed below)

    * `"region`" [string] the label of a user-defined region

    * `"cycles start period stop`" [Array(int)] the first entry is the start cycle, the second is the cycle period, and the third is the stop cycle or -1 in which case there is no stop cycle. A visualization dump shall be written at such cycles that satisfy cycle = start + n*period, for n=0,1,2,... and cycle < stop if stop != -1.0.

    * `"cycles start period stop n`" [Array(int)] if multiple cycles start period stop parameters are needed, then use these parameters with n=0,1,2,..., and not the single `"cycles start period stop`" parameter.

    * `"cycles`" [Array(int)] an array of discrete cycles that at which a visualization dump shall be written. 

    * `"times start period stop`" [Array(double)] the first entry is the start time, the second is the time period, and the third is the stop time or -1 in which case there is no stop time. A visualization dump shall be written at such times that satisfy time = start + n*period, for n=0,1,2,... and time < stop if stop != -1.0.

    * `"times start period stop n`" [Array(double) if multiple start period stop parameters are needed, then use this these parameters with n=0,1,2,..., and not the single  `"times start period stop`" parameter.

    * `"times`" [Array(double)] an array of discrete times that at which a visualization dump shall be written.

    * `"delimiter`" [string] the string used to delimit columns in the observation file output, default is `",`"

The following Observation Data functionals are currently supported.  All of them operate on the variables identified.

* `"Observation Data: Point`" returns the value of the field quantity at a point

* `"Observation Data: Integral`" returns the integral of the field quantity over the region specified

.. code-block:: xml

  <ParameterList name="Observation Data">
    <Parameter name="Observation Output Filename" type="string" value="obs_output.out"/>
    <ParameterList name="some observation name">
      <Parameter name="Region" type="string" value="some point region name"/>
      <Parameter name="Functional" type="string" value="Observation Data: Point"/>
      <Parameter name="Variable" type="string" value="Volumetric water content"/>
      <Parameter name="times" type="Array(double)" value="{100000.0, 200000.0}"/>

      <Parameter name="cycles" type="Array(int)" value="{100000, 200000, 400000, 500000}"/>
      <Parameter name="cycles start period stop" type="Array(int)" value="{0, 100, -1}" />

      <Parameter name="times start period stop 0" type="Array(double)" value="{0.0, 10.0, 100.0}"/>
      <Parameter name="times start period stop 1" type="Array(double)" value="{100.0, 25.0, -1.0}"/>
      <Parameter name="times" type="Array(double)" value="{101.0, 303.0, 422.0}"/>

    </ParameterList>
  </ParameterList>


Checkpoint Data
---------------

A user may request periodic dumps of Amanzi Checkpoint Data.  The user has no explicit control over the content of these files, but has the guarantee that the Amanzi run will be reproducible (with accuracies determined
by machine round errors and randomness due to execution in a parallel computing environment).  Therefore, output controls for Checkpoint Data are limited to file name generation and writing frequency, by numerical cycle number.

* `"Checkpoint Data`" [list] can accept a file name base [string] and cycle data [list] 
  used to generate the file base name or directory base name that is used in writing Checkpoint Data. 

  * `"file name base`" [string] ("checkpoint")
  
  * `"file name digits`" [int] (5)

  * `"cycles start period stop`" [Array(int)] the first entry is the start cycle, the second is the cycle period, and the third is the stop cycle or -1 in which case there is no stop cycle. A visualization dump shall be written at such cycles that satisfy cycle = start + n*period, for n=0,1,2,... and cycle < stop if stop != -1.0.

  * `"cycles start period stop n`" [Array(int)] if multiple cycles start period stop parameters 
    are needed, then use these parameters with n=0,1,2,..., and not the single `"cycles start period stop`" parameter.

  * `"cycles`" [Array(int)] an array of discrete cycles that at which a visualization dump shall be written. 

  * `"times start period stop`" [Array(double)] the first entry is the start time, 
    the second is the time period, and the third is the stop time or -1 in which 
    case there is no stop time. A visualization dump shall be written at such times 
    that satisfy time = start + n*period, for n=0,1,2,... and time < stop if stop != -1.0.

  * `"times start period stop n`" [Array(double) if multiple start period stop parameters 
    are needed, then use this these parameters with n=0,1,2,..., and not the 
    single  `"times start period stop`" parameter.

  * `"times`" [Array(double)] an array of discrete times that at which a visualization dump shall be written.

.. code-block:: xml

  <ParameterList name="Checkpoint Data">
    <Parameter name="file name base" type="string" value="chkpoint"/>
    <Parameter name="file name digits" type="int" value="5"/>

    <Parameter name="cycles start period stop" type="Array(int)" value="{0, 100, -1}" />
    <Parameter name="cycles" type="Array(int)" value="{999, 1001}" />

    <Parameter name="times start period stop 0" type="Array(double)" value="{0.0, 10.0, 100.0}"/>
    <Parameter name="times start period stop 1" type="Array(double)" value="{100.0, 25.0, -1.0}"/>
    <Parameter name="times" type="Array(double)" value="{101.0, 303.0, 422.0}"/>
  </ParameterList>

In this example, Checkpoint Data files are written when the cycle number is 
a multiple of 100.


Walkabout Data
--------------

A user may request periodic dumps of Walkabout Data. Output controls for Walkabout Data are limited to file name generation and writing frequency, by numerical cycle number or time.

* `"Walkabout Data`" [list] can accept a file name base [string] and cycle data [list] 
  used to generate the file base name or directory base name that is used in writing Checkpoint Data. 

  * `"file name base`" [string] ("walkabout")
  
  * `"file name digits`" [int] (5)

  * `"cycles start period stop`" [Array(int)] the first entry is the start cycle, the second is the cycle period, and the third is the stop cycle or -1 in which case there is no stop cycle. A visualization dump shall be written at such cycles that satisfy cycle = start + n*period, for n=0,1,2,... and cycle < stop if stop != -1.0.

  * `"cycles start period stop n`" [Array(int)] if multiple cycles start period stop parameters are needed, then use these parameters with n=0,1,2,..., and not the single `"cycles start period stop`" parameter.

  * `"cycles`" [Array(int)] an array of discrete cycles that at which a visualization dump shall be written. 

  * `"times start period stop`" [Array(double)] the first entry is the start time, the second is the time period, and the third is the stop time or -1 in which case there is no stop time. A visualization dump shall be written at such times that satisfy time = start + n*period, for n=0,1,2,... and time < stop if stop != -1.0.

  * `"times start period stop n`" [Array(double) if multiple start period stop parameters are needed, then use this these parameters with n=0,1,2,..., and not the single  `"times start period stop`" parameter.

  * `"times`" [Array(double)] an array of discrete times that at which a visualization dump shall be written.

.. code-block:: xml

  <ParameterList name="Walkabout Data">
    <Parameter name="file name base" type="string" value="walkabout"/>
    <Parameter name="file name digits" type="int" value="5"/>

    <Parameter name="cycles start period stop" type="Array(int)" value="{0, 100, -1}" />
    <Parameter name="cycles" type="Array(int)" value="{999, 1001}" />

    <Parameter name="times start period stop 0" type="Array(double)" value="{0.0, 10.0, 100.0}"/>
    <Parameter name="times start period stop 1" type="Array(double)" value="{100.0, 25.0, -1.0}"/>
    <Parameter name="times" type="Array(double)" value="{101.0, 303.0, 422.0}"/>
  </ParameterList>

In this example, walkabout data files are written when the cycle number is 
a multiple of 100.



Visualization Data
------------------

A user may request periodic writes of field data for the purposes of visualization.  The user will specify explicitly what is to be included in the file at each snapshot.  Visualization files can only be written 
at intervals corresponding to the numerical time step values or intervals corresponding to the cycle number; writes are controlled by time step cycle number.

* `"Visualization Data`" [list] can accept a file name base [string] and cycle data [list] that is used to generate the file base name or directory base name that is used in writing visualization data.  It can also accept a set of lists to specify which field quantities to write

  * `"file name base`" [string] ("amanzi_vis")
  
  * `"cycles start period stop`" [Array(int)] the first entry is the start cycle, the second is the cycle period, and the third is the stop cycle or -1 in which case there is no stop cycle. A visualization dump shall be written at such cycles that satisfy cycle = start + n*period, for n=0,1,2,... and cycle < stop if stop != -1.0.

  * `"cycles start period stop n`" [Array(int)] if multiple cycles start period stop parameters are needed, then use these parameters with n=0,1,2,..., and not the single `"cycles start period stop`" parameter.

  * `"cycles`" [Array(int)] an array of discrete cycles that at which a visualization dump shall be written. 

  * `"times start period stop`" [Array(double)] the first entry is the start time, the second is the time period, and the third is the stop time or -1 in which case there is no stop time. A visualization dump shall be written at such times that satisfy time = start + n*period, for n=0,1,2,... and time < stop if stop != -1.0.

  * `"times start period stop n`" [Array(double) if multiple start period stop parameters are needed, then use this these parameters with n=0,1,2,..., and not the single  `"times start period stop`" parameter.

  * `"times`" [Array(double)] an array of discrete times that at which a visualization dump shall be written.

  * `"dynamic mesh`" [bool] (false) write mesh data for every visualization dump, this facilitates visualizing deforming meshes.

  * `"regions`" [Array(string)] (empty array) write an array into the visualization file that can be used to identify a region or regions. The first entry in the regions array is marked with the value 1.0 in the array, the second with the value 2.0, and so forth. The code ignores entries in the regions array that are not valid regions that contain cells.

  * `"write partition`" [bool] (false) if this parameter is true, then write an array into the visualization file that contains the rank number of the processor that owns a mesh cell. 

.. code-block:: xml

  <ParameterList name="Visualization Data">
    <Parameter name="file name base" type="string" value="chk"/>
  
    <Parameter name="cycles start period stop" type="Array(int)" value="{0, 100, -1}" />
    <Parameter name="cycles" type="Array(int)" value="{999, 1001}" />

    <Parameter name="times start period stop 0" type="Array(double)" value="{0.0, 10.0, 100.0}"/>
    <Parameter name="times start period stop 1" type="Array(double)" value="{100.0, 25.0, -1.0}"/>
    <Parameter name="times" type="Array(double)" value="{101.0, 303.0, 422.0}"/>

    <Parameter name="dynamic mesh" type="bool" value="false"/>
  </ParameterList>


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


