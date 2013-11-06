========================================
Amanzi Native XML Input Specification V3
========================================

.. contents:: **Table of Contents**



ParameterList XML
=================

The Amanzi input file is an ASCII text XML-formatted file that must be framed 
at the beginning and end by the following statements:


.. code-block:: xml

  <ParameterList name="Main">

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

  <ParameterList name="Sub">
    <Parameter name="CFL" type="double" value="0.9"/>
    <Parameter name="ratio" type="Array(int)" value="{2, 2, 4}"/>
  </ParameterList>

In this example, the sublist "Sub" has a parameter named "CFL" that is a "double" and has 
the value of 0.9, and a Teuchos::Array<int> parameter named "ratio" such that ratio[0] = 2, 
ratio[1]=2, and ratio[2]=4.


Syntax of the Specification
===========================

* Input specification for each ParameterList entry consists of two parts.  
  First, a bulleted list defines the usage syntax and available options.  
  This is followed by example snipets of XML code to demonstrate usage.

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

Here, an `"X`" is defined by a `"Y`" and a `"Z`".  
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

In the MPC sublist the user specifies which process kernels are on or off, which 
flow model is active, and the time integration mode that the MPC should run in.

To turn a particular process kernel on or off use these options:

 * `"disable Transport_PK`" [string], valid options are `"yes`" or `"no`".

 * `"disable Flow_PK`" [string], valid options are `"yes`" or `"no`".

 * `"Chemistry Model`" [string], valid options are `"On`" or `"Off`".

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

The MPC list must have a sublist named `"Time Integration Mode`" if flow is enabled.
This list must have exactly one of the following three sublists

.. code-block:: xml

      <ParameterList name="Steady">
        <Parameter name="Start" type="double" value="0.0"/>
        <Parameter name="End" type="double" value="5.0"/>
        <Parameter name="Initial Time Step" type="double" value="0.1"/>
      </ParameterList>

or

.. code-block:: xml

      <ParameterList name="Initialize To Steady">
        <Parameter name="Start" type="double" value="0.0"/>
        <Parameter name="Switch" type="double" value="0.5"/>
        <Parameter name="End" type="double" value="5.0"/>
        <Parameter name="Steady Initial Time Step" type="double" value="0.1"/>
        <Parameter name="Transient Initial Time Step" type="double" value="0.1"/>
      </ParameterList>

or

.. code-block:: xml

      <ParameterList name="Transient">
        <Parameter name="Start" type="double" value="0.0"/>
        <Parameter name="End" type="double" value="5.0"/>
        <Parameter name="Initial Time Step" type="double" value="0.1"/>
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

Example

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
        <Parameter name="Start" type="double" value="0.00000000000000000e+00"/>
        <Parameter name="Switch" type="double" value="5.00000000000000000e-01"/>
        <Parameter name="End" type="double" value="5.00000000000000000e+00"/>
        <Parameter name="Steady Initial Time Step" type="double" value="1.00000000000000006e-01"/>
        <Parameter name="Transient Initial Time Step" type="double" value="1.00000000000000006e-01"/>
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



State
=====

State allows the user to initialize physical fields using a variety of 
tools. 

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
`"Mesh Block i`" with non-overlapping ranges. 
The required parameters for each function are `"region`", `"component`",
and the function itself.

.. code-block:: xml

  <ParameterList name="porosity"> 
    <ParameterList name="function">
      <ParameterList name="Mesh Block 1">
        <Parameter name="region" type="string" value="Computational domain"/>
        <Parameter name="component" type="string" value="cell"/>
        <ParameterList name="function">
          <ParameterList name="function-constant">
            <Parameter name="value" type="double" value="0.2"/>
          </ParameterList>
        </ParameterList>
      </ParameterList>
      <ParameterList name="Mesh Block 2">
        ...
      </ParameterList>
    </ParameterList>
  </ParameterList>


Initialization of tensor fields
-------------------------------

A variable tensor (or vector) field is defined similarly to 
a variable scalar field. 
The difference lies in the definition of the function which
is now a multi-values function.
The required parameters are `"Number of DoFs`" and `"Function type`". 

.. code-block:: xml

  <ParameterList name="function">
    <Parameter name="Number of DoFs" type="int" value="2"/>
    <Parameter name="Function type" type="string" value="composite function"/>
    <ParameterList name="DoF 1 Function">
      <ParameterList name="function-constant">
        <Parameter name="value" type="double" value="1.9976e-12"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="DoF 2 Function">
      <ParameterList name="function-constant">
        <Parameter name="value" type="double" value="1.9976e-13"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>

Example
-------

The complete example of a state initialization is below.

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
            <Parameter name="region" type="string" value="Computational domain"/>
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
          <ParameterList name="Mesh Block 1">
            <Parameter name="region" type="string" value="Material 1 Region"/>
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
          <ParameterList name="Mesh Block 2">
            <Parameter name="region" type="string" value="Material 2 Region"/>
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


Initialization from a file
--------------------------

Some data can be initialized from files. Additional sublist has to be added to
named sublist of the `"state`" list with the file name and names of attributes. 
The provided data will override results of other initialization tools. 
Here is an example:

.. code-block:: xml

  <ParameterList name="exodus file initialization">
    <Parameter name="file" type="string" value="mesh_with_data.exo"/>
    <Parameter name="attribute" type="string" value="perm"/>
  </ParameterList>


Flow
====

Flow sublist includes exactly one sublist, either `"Darcy Problem`" or `"Richards Problem`".
Structure of both sublists is quite similar. We make necessary comments on differences.

Water retention models
-----------------------

User defines water retention models in sublist `"Water retention models`". 
It contains as many sublists, 
e.g. `"Soil 1`", `"Soil 2`", etc, as there are different soils. 
These models are associated with non-overlapping regions. Each of the sublists `"Model N`" 
includes a few mandatory parameters: a region name, model name, and parameters for the selected model.
The available models are `"van Genuchten`", `"Brooks Corey`", and `"fake`". 
The later is used to set up an analytic solution for convergence study. 
The available models for the relative permeability are `"Mualem`" (default) and `"Burdine`".
An example of the van Genuchten model specification is:

.. code-block:: xml

    <ParameterList name="Soil 1">
       <Parameter name="region" type="string" value="Top Half"/>
       <Parameter name="water retention model" type="string" value="van Genuchten"/>
       <Parameter name="van Genuchten alpha" type="double" value="0.000194"/>
       <Parameter name="van Genuchten m" type="double" value="0.28571"/>
       <Parameter name="van Genuchten l" type="double" value="0.5"/>
       <Parameter name="residual saturation" type="double" value="0.103"/>
       <Parameter name="relative permeability model" type="string" value="Mualem"/>
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


Amanzi performs rudimentary checks of validity of the provided parameters. 
The relative permeability curves can be calculated and saved in the file krel_pc.txt
and krel_sat.txt using the following optional commands (that go to `"Water Retention Models`" list):

.. code-block:: xml

    <Parameter name="plot krel-pc curves" type="Array(double)" value="{0.0, 0.1, 3000.0}"/>
    <Parameter name="plot krel-sat curves" type="Array(double)" value="{0.0001, 0.01, 1.0}"/>

The triple of doubles means the starting capillary pressure (resp., saturation), the period, and 
the final capillary pressure (resp., saturation).
Each line in the output file will contain the capillary pressure (resp., saturation) and relative 
permeability values for all water retention models in the order they appear in the input spec.


Boundary conditions
-------------------

Boundary conditions are defined in sublist `"boundary conditions`". Four types of boundary 
conditions are supported:

* `"pressure`" [list] Dirichlet boundary condition, a pressure is prescribed on a surface region. 

* `"mass flux`" [list] Neumann boundary condition, an outward mass flux is prescribed on a surface region.
  This is the default boundary condition. If no condition is specified on a mesh face, zero flux 
  boundary condition is used implicitly. 

* `"static head`" [list] Dirichlet boundary condition, the hydrostatic pressure is prescribed on a surface region.

* `"seepage face`" [list] Seepage face boundary condition, a dynamic combination of the `"pressure`" and 
  `"mass flux`" boundary conditions on a region. 
  The atmospheric pressure is prescribed if internal pressure is higher. Otherwise, the outward mass flux is prescribed. 

The following example includes all four types of boundary conditions. The boundary of a square domain 
is split into six pieces. Constant function is used for simplicity and can be replaced by any
of the other available functions:

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

The above boundary conditions are the four major models supported by Amanzi. In addition to
that each model may support a few submodels. A submodel is defined by additional
parameters described below. Mix and match of parameters is allowed.

* `"rainfall`" [bool] indicates that the mass flux is defined with respect to the gravity 
  vector and the actual influx depends on boundary slope. Default value is `"false`".

* `"relative to top`" [bool] indicates that the static head is defined with respect
  to the top boundary (a curve in 3D) of the specified regions. Support of 2D is turned off.
  Default value is `"false`". 

* `"submodel`" [string] indicates different models for the seepage face boundary condition.
  It can take values `"PFloTran`", `"FACT`", and `"Amanzi`". The first option leads to a 
  discontinuous change of the boundary condition type from the infiltration to pressure. 
  The second option is described
  in the document on mathematical models. It employs a smooth transition from the infiltration 
  to mixed boundary condition. The third option combines the above two. Is uses a smooth transition
  from the infiltration to pressure boundary condition. 
  Default value is `"Amanzi`".

Here is an example:

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
of sublist `"source terms`" follows the specification of boundary conditions. 
Again, constant functions can be replaced by any of the available time-functions:

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

* `"spatial distribution method`" [string] identifies a method for distributing
  source Q over the specified regions. The available options are `"volume`",
  `"none`", and `"permeability`". For option `"none`" the source term Q is measured
  in [kg/m^3/s]. For the other options, it is measured in [kg/s]. When the source function
  is defined over a few regions, Q will be distributed independently over each region.
  Default is `"none`".


Initial Guess Pseudo Time Integrator
-------------------------------------

The sublist `"initial guess pseudo time integrator`" defines parameters controlling linear and 
nonlinear solvers during calculation of the initial guess time integration. Here is an example:

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
       <Parameter name="linear solver" type="string" value="CG with HypreAMG"/>
     </ParameterList>

     <ParameterList name="Picard">
       <ParameterList name="Picard parameters">
         <Parameter name="convergence tolerance" type="double" value="1e-08"/>
         <Parameter name="maximum number of iterations" type="int" value="400"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>

Detailed description of parameters is in the next two subsection.


Steady State Time Integrator
----------------------------

The sublist `"steady state time integrator`" defines parameters controlling linear and 
nonlinear solvers during steady state time integration. Here is an example:

.. code-block:: xml

   <ParameterList name="steady state time integrator">
     <Parameter name="time integration method" type="string" value="BDF1"/>
     <Parameter name="error control options" type="Array(string)" value="{pressure, saturation}"/>
     <Parameter name="linear solver" type="string" value="GMRES with HypreAMG"/>

     <ParameterList name="initialization">
       <Parameter name="method" type="string" value="saturated solver"/>
       <Parameter name="linear solver" type="string" value="CG with HypreAMG"/>
       <Parameter name="clipping pressure value" type="double" value="50000.0"/>
     </ParameterList>

     <ParameterList name="pressure-lambda constraints">
       <Parameter name="method" type="string" value="projection"/>
       <Parameter name="linear solver" type="string" value="CG with HypreAMG"/>
     </ParameterList>

     <ParameterList name="BDF1">
       <Parameter name="max iterations" type="int" value="15"/>
       <Parameter name="min iterations" type="int" value="10"/>
       <Parameter name="limit iterations" type="int" value="20"/>
       <Parameter name="nonlinear tolerance" type="double" value="1e-05"/>
       <Parameter name="time step reduction factor" type="double" value="0.8"/>
       <Parameter name="time step increase factor" type="double" value="1.25"/>
       <Parameter name="max time step" type="double" value="6e+10"/>
       <Parameter name="max preconditioner lag iterations" type="int" value="20"/>
       <Parameter name="error abs tol" type="double" value="1.0"/>
       <Parameter name="error rel tol" type="double" value="0.0"/>
       <Parameter name="time step increase factor" type="double" value="1.2"/>
       <Parameter name="max divergent iterations" type="int" value="3"/>
       <Parameter name="nonlinear iteration damping factor" type="double" value="1.0"/>
       <Parameter name="nonlinear iteration initial guess extrapolation order" type="int" value="1"/>
       <Parameter name="restart tolerance relaxation factor" type="double" value="1.0"/>
       <Parameter name="restart tolerance relaxation factor damping" type="double" value="1.0"/>
       <Parameter name="nonlinear iteration divergence factor" type="double" value="1e+03"/>

       <Parameter name="initial time step" type="double" value="1e-07"/>
       <Parameter name="maximum time step" type="double" value="1e+10"/>
       <Parameter name="maximum number of iterations" type="int" value="400"/>
       <Parameter name="convergence tolerance" type="double" value="1e-12"/>
       <Parameter name="maximal number of iterations" type="int" value="200"/>
       <Parameter name="start time" type="double" value="0.0"/>
       <Parameter name="end time" type="double" value="100.0"/>
     </ParameterList>
   </ParameterList>

The parameters used here are

* `"time integration method`" [string] defines a time integration method.
  The available options are `"BDF1`", `"BDF2`", `"Picard`", and `"backward Euler`".

* `"error control options`" [Array(string)] lists various error control options. 
  A nonlinear solver is terminated when all listed options are passed. 
  The available options are `"pressure`", `"saturation`", and `"residual`". 
  All errors are relative, i.e. dimensionless. 
  The error in pressure is compared with capillary pressure plus atmospheric pressure. 
  The other two error are compared with 1. 
  The option `"pressure`" is always active during steady-state time integration.
  The option  `"saturation`" is always active during transient time integration.

* `"time stepping strategy`" [string] allows one to define an adaptive time step increment 
  through an error estimator. The only available option is `"adaptive`". It is supported
  for the Darcy flow only. 
  The error estimator can be controlled via two parameters in the list `"time integration method`" 
  called `"absolute error tolerance`" and `"relative error tolerance`". The default values
  for these parameters are 0.001. 

* `"BDF1`" [list] list specified in `"time integration method`".
  It includes the following parameters.

  * `"time step increase factor`" [double] defines geometric grow rate for the
    initial time step. If adaptive time stepping strategy is specified, this
    parameter is ignored. Default is 1.0.

  * Other parameters will be described later.

* `"initialization`" [list] defines parameters for calculating initial pressure guess.
  It can be used to obtain pressure field which is consistent with the boundary conditions.
  Default is empty list.

  * `"method`" [string] refers to a constraint enforcement method. The only 
    available option is `"projection`" which is default.

  * `"linear solver`" [string] refers to a solver sublist of the list `"Solvers`".

  * `"clipping saturation value`" [double] is an experimental option. It is used 
    after pressure initialization to cut-off small values of pressure. By default, the 
    pressure threshold is equal to the atmospheric pressure.
    The new pressure is calculated based of the provided saturation value. Default is 0.6.

  * `"clipping pressure value`" [double] is an experimental option. It is used 
    after pressure initialization to cut-off small values of pressure below the provided
    value.

* `"enforce pressure-lambda constraints`" [list] each time the time integrator is 
  restarted, we may re-enforce the pressure-lambda relationship for new boundary conditions. 
  Default is empty list.

  * `"method`" [string] refers to a constraint enforcement method. The only 
    available option is `"projection`" which is default.

  * `"linear solver`" [string] refers to a solver sublist of the list `"Solvers`".

* `"BFD1`" [list] the named list used to control the nonlinear solver.


Transient Time Integrator
-------------------------

The sublist `"transient time integrator`" defines parameters controlling linear and 
nonlinear solvers during transient time integration. Its parameters are similar to 
that in the sublist `"steady state time integrator`" except for parameters controlling
pressure re-initialization. Here is an example:

.. code-block:: xml

   <ParameterList name="transient time integrator">
     <Parameter name="time integration method" type="string" value="BDF1"/>
     <Parameter name="error control options" type="Array(string)" value="{pressure, saturation}"/>
     <Parameter name="linear solver" type="string" value="GMRES with HypreAMG"/>
     <Parameter name="time stepping strategy" type="string" value="adaptive"/>

     <ParameterList name="initialization">
       <Parameter name="method" type="string" value="projection"/>
       <Parameter name="linear solver" type="string" value="CG with HypreAMG"/>
     </ParameterList>

     <ParameterList name="pressure-lambda constraints">
       <Parameter name="method" type="string" value="projection"/>
       <Parameter name="linear solver" type="string" value="CG with HypreAMG"/>
     </ParameterList>

     <ParameterList name="BDF1">
       ...
     </ParameterList>
   </ParameterList>

The parameters were defined above. A non-empty `"initialization`" list 
may be useful for a transient saturated simulation.


Other Parameters
-----------------------------

The remaining `"Flow`" parameters are

* `"atmospheric pressure`" [double] defines the atmospheric pressure, [Pa].

* `"relative permeability`" [string] defines a method for calculating relative
  permeability. The available self-explanatory options `"upwind with gravity`",
  are `"upwind with Darcy flux`", `"arithmetic mean`" and `"cell centered`". 
  The first three calculate the relative permeability on mesh interfaces.

* `"discretization method`" [string] helps to test new discretization methods. 
  The available options are `"mfd scaled`", `"optimized mfd scaled`",
  `"two-point flux approximation`", `"two point flux approximation`" and
  `"support operator`". The last option reproduces discretization method implemented in RC1. 
  The third option is recommended for orthogonal meshes and diagonal absolute permeability.
  The second option is still experimental (no papers were published) and produces 
  an optimal discretization.

* `"plot time history`" [bool] produces an ASCII file with time history when exists.

* `"VerboseObject`" [list] defines default verbosity level for the process kernel.
  If it does not exists, it will be created on a fly and verbosity level will be set to `"high`".
  Here is an example:

.. code-block:: xml

    <ParameterList name="VerboseObject">
      <Parameter name="Verbosity Level" type="string" value="medium"/>
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

Here is an example:

.. code-block:: xml

   <ParameterList name="Transport">
     <Parameter name="CFL" type="double" value="1.0"/>
     <Parameter name="spatial discretization order" type="int" value="1"/>
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

Two discretization methods that preserve the maximum principles are 
`"two point flux approximation`" and `"nonliner finite volume`". 
The first one may show significant numerical dispersion on unstructured meshes, 
the second-one is more accurate but also is a few times more expensive.

.. code-block:: xml

   <ParameterList name="Dispersivity">
     <Parameter name="numerical method" type="string" value="two point flux approximation"/>
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

Parameter `"preconditioner`" will be replaced with more appropriate `"linear solver`".
 

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
         <ParameterList name="source for east well">   <!-- user defined name -->
           <Parameter name="regions" type="Array(string)" value="{Top, Bottom}"/>
             <ParameterList name="boundary concentration">
               <ParameterList name="function-constant">  <!-- any time function -->
                 <Parameter name="value" type="double" value="0.0"/>
               </ParameterList>
             </ParameterList>
           </ParameterList>
         </ParameterList>
         <ParameterList name="source for west well">   <!-- user defined name -->
           ...
         </ParameterList>
       </ParameterList>

       <ParameterList name="Sugar syrop"> <!-- Next component --> 
         ...
       </ParameterList>
     </ParameterList>

     <ParameterList name="outward flux">  <!-- Future boundary conditions -->
     </ParameterList>
   </ParameterList>

Sources and Sinks
-----------------

The external sources are typically located at pumping wells. The structure
of sublist `"source terms`" includes only sublists named after components. 
Again, constant functions can be replaced by any available time-function:
Note that the source values are set up separately for each component:

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
           <ParameterList name="source for Well west">   <!-- user defined name -->
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
    

* `"spatial distribution method`" [string] identifies a method for distributing
  source Q over the specified regions. The available options are `"volume`",
  `"none`", and `"permeability`". For option `"none`" the source term Q is measured
  in [mol/m^3/s]. For the other options, it is measured in [mol/s]. When the source function
  is defined over a few regions, Q will be distributed independently over each region.
  Default is `"none`".


Other parameters
-----------------

The `"Transport`" parameters useful for developers are:

* `"enable internal tests`" [string] various internal tests will be executed during
  the run time. The default value is `no`.
   
* `"internal tests tolerance`" [double] tolerance for internal tests such as the 
  divergence-free condition. The default value is 1e-6.


Linear Solvers
==============

This list contains sublists for various linear solvers such as PCG and GMRES.
Here is and example:

.. code-block:: xml

     <ParameterList name="Solvers">
       <ParameterList name="GMRES with HypreAMG">
         <Parameter name="iterative method" type="string" value="gmres"/>
         <Parameter name="error tolerance" type="double" value="1e-12"/>
         <Parameter name="maximum number of iterations" type="int" value="400"/>
         <Parameter name="convergence criteria" type="Array(string)" value="{relative residual}"/>
         <Parameter name="preconditioner" type="string" value="Hypre AMG"/>
         <Parameter name="size of Krylov space" type="int" value="10"/>

         <ParameterList name="VerboseObject">
           <Parameter name="Verbosity Level" type="string" value="high"/>
         </ParameterList>
       </ParameterList>
     </ParameterList>

The name `"GMRES with Hypre AMG`" is selected by the user.
It can be used by a process kernel lists to define a solver.
The verbose object is discussed below.

* `"iterative method`" [string] defines a type of Krylov-based method. The parameters
  include `"pcg'" and `"gmres`".

* `"error tolerance`" [double] is used in the convergence test. The default value is 1e-6.

* `"maximum number of iterations`" [int] is used in the convergence test. The default is 100.

* `"convergence criteria`" [Array(string)] specifies multiple convergence criteria. The list
  may include `"relative residual`", `"relative rhs`" (default), and `"absolute residual`".

* `"preconditioner`" [string] is name in the list of preconditioners. If missing, the identity
  preconditioner will be employed. Support of the identity preconditioner is the work in progress.

* `"size of Krylov space`" [int] is used in GMRES iterative method. The default value is 10.


Nonlinear Solvers
=================
This list contains the name of a nonlinear solver. Curently, there are
two available options: nka and newton. Using of the Newton method
assumes that two-point flux discretization will be implemented.

.. code-block:: xml

  <ParameterList name="Nonlinear solvers">
    <Parameter name="solver" type="string" value="newton"/>
  </ParameterList>

Preconditioners
===============

Version 2 of the native input spec introduces this list. It contains sublists for various
preconditioners required by a simulation. At the moment, we support Trilinos multilevel 
preconditioner and Hypre BoomerAMG preconditioner. Here is an example:

.. code-block:: xml

     <ParameterList name="Preconditoners">
       <ParameterList name="Trilinos ML">
          <Parameter name="discretization method" type="string" value="optimized mfd scaled"/>
          <Parameter name="type" type="string" value="trilinos ml"/>
          <ParameterList name="ML Parameters">
            ... 
         </ParameterList>
       </ParameterList>

       <ParameterList name="Hypre AMG">
          <Parameter name="discretization method" type="string" value="optimized mfd scaled"/>
          <Parameter name="type" type="string" value="boomer amg"/>
          <ParameterList name="BoomerAMG Parameters">
            ...
          </ParameterList>
       </ParameterList>

       <ParameterList name="Block ILU">
          <Parameter name="discretization method" type="string" value="optimized mfd scaled"/>
          <Parameter name="type" type="string" value="block ilu"/>
          <ParameterList name="Block ILU Parameters">
            ...
          </ParameterList>
       </ParameterList>
     </ParameterList>

Names `"Trilinos ML`" and `"Hypre AMG`" are selected by the user.
They can be used by a process kernel lists to define a preconditioner.

Hypre AMG
---------

Internal parameters of Boomer AMG includes

.. code-block:: xml

   <ParameterList name="BoomerAMG Parameters">
     <Parameter name="tolerance" type="double" value="0.0"/>
     <Parameter name="smoother sweeps" type="int" value="3"/>
     <Parameter name="cycle applications" type="int" value="5"/>
     <Parameter name="strong threshold" type="double" value="0.5"/>
   </ParameterList>

Trilinos ML
-----------

Internal parameters of Trilinos ML includes

.. code-block:: xml

   <ParameterList name="ML Parameters">
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

The internal parameters of the block ILU are as follows:

.. code-block:: xml

   <ParameterList name="Block ILU Parameters">
     <Parameter name="fact: relax value" type="double" value="1.00000000000000000e+00"/>
     <Parameter name="fact: absolute threshold" type="double" value="0.00000000000000000e+00"/>
     <Parameter name="fact: relative threshold" type="double" value="1.00000000000000000e+00"/>
     <Parameter name="fact: level-of-fill" type="int" value="0"/>
     <Parameter name="overlap" type="int" value="0"/>
     <Parameter name="schwarz: combine mode" type="string" value="Add"/>
   </ParameterList>


Indentity
---------

The identity preconditioner is instantiated if either no preconditioner is
pecified or the specified preconditoner list does not exists.


Mesh
====

Amanzi supports both structured and unstructured numerical solution approaches.
This flexibility has a direct impact on the selection and design of the underlying 
numerical algorithms, the style of the software implementations, and, ultimately, 
the complexity of the user-interface.  
`"Mesh`" is used to select between the following options:

* `"Structured`": This instructs Amanzi to use BoxLib data structures and an associated paradigm to numerically represent the flow equations.  Data containers in the BoxLib software library, developed by CCSE at LBNL, are based on a hierarchical set of uniform Cartesian grid patches.  `"Structured`" requires that the simulation domain be a single coordinate-aligned rectangle, and that the "base mesh" consists of a logically rectangular set of uniform hexahedral cells.  This option supports a block-structured approach to dynamic mesh refinement, wherein successively refined subregions of the solution are constructed dynamically to track "interesting" features of the evolving solution.  The numerical solution approach implemented under the `"Structured`" framework is highly optimized to exploit regular data and access patterns on massively parallel computing architectures.

* `"Unstructured`": This instructs Amanzi to use data structures provided in the Trilinos software framework.  To the extent possible, the discretization algorithms implemented under this option are largely independent of the shape and connectivity of the underlying cells.  As a result, this option supports an arbitrarily complex computational mesh structure that enables users to work with numerical meshes that can be aligned with geometrically complex man-made or geostatigraphical features.  Under this option, the user typically provides a mesh file that was generated with an external software package.  The following mesh file formats are currently supported: `"Exodus 2`" (see example), `"MSTK`" (see example), `"MOAB`" (see example).  Amanzi also provides a rudmentary capability to generate unstructured meshes automatically.

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

Example:

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

  * OBSERVATION [list] user-defined label, can accept values for `"Variables`", `"Functional`", `"Region`", `"times`", and TSPS (see below).

    * `"Variables`" [Array(string)] a list of field quantities taken from the list of 
      available field quantities:

      * Volumetric water content [volume water / bulk volume]
      * Aqueous saturation [volume water / volume pore space]
      * Aqueous pressure [Pa]
      * Hydraulic Head [m] 
      * XXX Aqueous concentration [moles of solute XXX / volume water in MKS] (name formed by string concatenation, given the definitions in `"Phase Definition`" section)
      * X-, Y-, Z- Aqueous volumetric fluxe [m/s]
      * MaterialID

    * `"Functional`" [string] the label of a function to apply to each of the variables in the variable list (Function options detailed below)

    * `"Region`" [string] the label of a user-defined region

    * `"cycles start period stop`" [Array(int)] the first entry is the start cycle, the second is the cycle period, and the third is the stop cycle or -1 in which case there is no stop cycle. A visualization dump shall be written at such cycles that satisfy cycle = start + n*period, for n=0,1,2,... and cycle < stop if stop != -1.0.

    * `"cycles start period stop n`" [Array(int)] if multiple cycles start period stop parameters are needed, then use these parameters with n=0,1,2,..., and not the single `"cycles start period stop`" parameter.

    * `"cycles`" [Array(int)] an array of discrete cycles that at which a visualization dump shall be written. 

    * `"times start period stop`" [Array(double)] the first entry is the start time, the second is the time period, and the third is the stop time or -1 in which case there is no stop time. A visualization dump shall be written at such times that satisfy time = start + n*period, for n=0,1,2,... and time < stop if stop != -1.0.

    * `"times start period stop n`" [Array(double) if multiple start period stop parameters are needed, then use this these parameters with n=0,1,2,..., and not the single  `"times start period stop`" parameter.

    * `"times`" [Array(double)] an array of discrete times that at which a visualization dump shall be written.


The following Observation Data functionals are currently supported.  All of them operate on the variables identified.

* `"Observation Data: Point`" returns the value of the field quantity at a point

* `"Observation Data: Integral`" returns the integral of the field quantity over the region specified


Example:

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

  * `"cycles start period stop n`" [Array(int)] if multiple cycles start period stop parameters are needed, then use these parameters with n=0,1,2,..., and not the single `"cycles start period stop`" parameter.

  * `"cycles`" [Array(int)] an array of discrete cycles that at which a visualization dump shall be written. 

  * `"times start period stop`" [Array(double)] the first entry is the start time, the second is the time period, and the third is the stop time or -1 in which case there is no stop time. A visualization dump shall be written at such times that satisfy time = start + n*period, for n=0,1,2,... and time < stop if stop != -1.0.

  * `"times start period stop n`" [Array(double) if multiple start period stop parameters are needed, then use this these parameters with n=0,1,2,..., and not the single  `"times start period stop`" parameter.

  * `"times`" [Array(double)] an array of discrete times that at which a visualization dump shall be written.

  * `"walkabout`" [bool] (false) include postprocessed output for walkabout in the checkpoint files.

Example:

.. code-block:: xml

  <ParameterList name="Checkpoint Data">
    <Parameter name="file name base" type="string" value="chkpoint"/>
    <Parameter name="file name digits" type="int" value="5"/>

    <Parameter name="cycles start period stop" type="Array(int)" value="{0, 100, -1}" />
    <Parameter name="cycles" type="Array(int)" value="{999, 1001}" />

    <Parameter name="times start period stop 0" type="Array(double)" value="{0.0, 10.0, 100.0}"/>
    <Parameter name="times start period stop 1" type="Array(double)" value="{100.0, 25.0, -1.0}"/>
    <Parameter name="times" type="Array(double)" value="{101.0, 303.0, 422.0}"/>

    <Parameter name="walkabout" type="bool" value="false"/>
  </ParameterList>

In this example, Checkpoint Data files are written when the cycle number is 
a multiple of 100.

Additional data are written to this file when parameter `"walkabout`"
is set to true.  


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

Example:

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


