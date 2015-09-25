==========================================
Amanzi-U Native XML Input Specification V7
==========================================

.. contents:: **Table of Contents**


Overview
========

This is a continuously evolving specification format used by the code developers. 
It is main purpose is to develop and test new capabilities without disruption of end-users.
Parameters labeled by [WIP] (Work-In-Progress) are under development.
Parameters labeled by [O] (Obsolete) are old capabilities and will be removed soon.


Changes V6 -> V7
================

* Observations use lower-case names.
* Added dual porosity model to flow and transport.


ParameterList XML
=================

The Amanzi input file is an ASCII text XML-formatted file that must be framed 
at the beginning and end by the following statements:

.. code-block:: xml

  <ParameterList name="Main">
    various lists and sublists
  </ParameterList>

The value of *name* can be anything (*Main* in this example).  
A ParameterList consists of just two types of entries: Parameter and ParameterList.  
ParameterLists are labeled with *name* [string], while Parameters have a separate 
fields called *name* [string], *type* [string] and *value* [TYPE], where TYPE can 
be any of the following: double, int, bool, string, Array(double), Array(int), 
and Array(string).  
The value of the parameter is given in quotes (e.g. value="2.7e3").  
Array data is specified as a single comma-delimited string bounded by {}'s (e.g. value="{2.4, 2.1, 5.7}").

.. code-block:: xml

  <ParameterList name="Main">
    <Parameter name="cfl" type="double" value="0.9"/>
    <Parameter name="ratio" type="Array(int)" value="{2, 1, 4}"/>
  </ParameterList>

In this example, the list *Main* has parameter *cfl* that is the double with 
value 0.9, and parameter *ratio* that is the integer array such that ratio[0] = 2, 
ratio[1]=1, and ratio[2]=4.


Syntax of the specification
---------------------------

Input specification for each ParameterList entry consists of two parts.  
First, a bulleted list defines the usage syntax and available options.  
This is followed by example snippets of XML code to demonstrate usage.

In many cases, the input specifies data for a particular parameterized model, and Amanzi 
supports a number of parameterizations.  
For example, initial data might be uniform (the value is required), or linear in y (the value 
and its gradient are required).  
Where Amanzi supports a number of parameterized models for parameter *model*, the available 
models will be listed by name, and then will be described in the subsequent section.  
In the manufactured example below, the specification looks as follows:

* SOIL [list] accepts parameters that describes properties of this soil.

  * `"region`" [string] defines a subdomain of the computational domain.

  * `"model`" [list] specifies a model for the soil. Available options are `"van Genuchten`" 
    and `"Brooks-Corey`".

Here SOIL is defined by a *region* and a *model*.  
The *region* is a string parameter but the *model* is given by a sublist with its own set of parameters.
The parameter for *model* can be described in the same section or in a separate section
of this document. For instance, the local description may look like:

* `"model`" [list] specifies a model for the soil. Available options are `"van Genuchten`"
  and `"Brooks-Corey`".
  The option `"van Genuchten`" requires `"m`" [double].
  The option `"Brooks-Corey`" requires `"lambda`" [double] and `"alpha`" [double].

Each part of the spec is illustrated by an example followed by optional comments:

.. code-block:: xml

   <ParameterList name="water retention models">
     <ParameterList name="SOIL">
       <Parameter name="region" type="string" value="TOP_DOMAIN"/>
       <ParameterList name="Brooks-Corey">
         <Parameter name="lambda" type="double" value="0.7"/>
         <Parameter name="alpha" type="double" value="1e-3"/>
       </ParameterList>   
     </ParameterList>   
   </ParameterList>   
 
This defines soil properties in region TOP_DOMAIN using the
Brocks-Corey model with parameters *lambda=0.7* and *alpha=1e-3*.

Additional conventions:

* Reserved keywords and labels are *italicized* in discussions and `"quoted and italicized`" in the spec.
  These are usually labels or values of parameters 
  in the input file and must match (using XML matching rules) the specified or allowable values.

* User-defined labels are marked with ALL_CAPS in this documents.
  In practice, no rules are imposed on these names.

* For developers: we are gradually migrating to low-case naming convention for most parameters
  with exceptions of proper name (e.g. CO2, Linux), personal names (e.g. van Genuchten), and high-level
  lists (e.g. Mesh). Abbreviations are capitalized only in the list names.

* Lists with too many parameters are described using multiple sections and multiple examples.
  For most examples we show name of the parent sublist.


Verbose output
--------------

Output of all components of Amanzi is controlled by a standard verbose 
object sublist. This sublist can be inserted in almost any significant
component of this spec to produce a verbose output, see the embedded examples.
If this list is not specified, the default verbosity value is used.
The name *Verbosity Level* is reserved by Trilinos.

.. code-block:: xml

   <ParameterList name="VerboseObject">
     <Parameter name="Verbosity Level" type="string" value="high"/>
   </ParameterList>


Units
-----

Amanzi's internal default units are SI units.


Cycle driver
============

New multi-processor cycle driver which provides more flexibility
to handle multiphysics process kernels (PKs) and multiple time periods.

* `"component names`" [Array(string)] provides the list of species names.
  It is required for reactive transport.

* `"number of liquid components`" [int] is the number of liquid components. 
   
* `"time periods`" [list] contains the list of time periods involved in the simulation.
  the number of periods is not limited.

  * `"TP #`" [list] defines a particular time period. The numbering
    should be sequential starting with 0.

    * `"PK tree`" [list] describes a hierarchical structure of the process kernels
      that reflect their weak and strong coupling.

      * `"PKNAME`"  [list] name of PK which is used in the
        simulation. Name can be arbitrary but the sublist with the same name
        should exist in the list of PKs (see below).

      * `"PK type`" [string] specifies the type of PK supported by Amanzi. At the moment
        available options are (`"darcy`", `"richards`", `"transport`", `"reactive
        transport`", `"flow reactive transport`", and `"chemistry`").
 
      * `"start period time`" [double] is the start time of the current time period.

      * `"end period time`" [double] is the end time of the current time period.

      * `"maximum cycle number`" [int] is the maximum allowed number of cycles in 
        the current time period. Special value -1 means unlimited number of cycles.

      * `"initial time step`" is the initial time step for the current time period.

.. code-block:: xml

  <ParameterList>  <!-- parent list -->
    <ParameterList name="Cycle Driver">
      <Parameter name="component names" type="Array(string)" value="{H+, Na+, NO3-, Zn++}"/>
      <Parameter name="number of liquid components" type="int" value="4"/>
      <ParameterList name="time periods">
        <ParameterList name="TP 0">
          <ParameterList name="PK Tree">
            <ParameterList name="FLOW and REACTIVE TRANSPORT">
              <Parameter name="PK type" type="string" value="flow reactive transport"/>
              <ParameterList name="REACTIVE TRANSPORT">
                <Parameter name="PK type" type="string" value="reactive transport"/>
                <ParameterList name="TRANSPORT">
                  <Parameter name="PK type" type="string" value="transport"/>
                </ParameterList>
                <ParameterList name="CHEMISTRY">
                  <Parameter name="PK type" type="string" value="chemistry"/>
                </ParameterList>
              </ParameterList>
              <ParameterList name="FLOW">
                <Parameter name="PK type" type="string" value="darcy"/>
              </ParameterList>
            </ParameterList>
          </ParameterList>
          <Parameter name="start period time" type="double" value="0.0"/>
          <Parameter name="end period time" type="double" value="1.5778463e+09"/>
          <Parameter name="maximum cycle number" type="int" value="-1"/>
          <Parameter name="initial time step" type="double" value="1.57680e+05"/>
        </ParameterList>

        <ParameterList name="TP 1">
        ... 
        </ParameterList>
      </ParameterList>
    </ParameterList>
  </ParameterList>

In this simulation, we use the PK called *flow reactive transport*. It is
defined internally as sequential application of two PKs, *flow* and *reactive transport*.
The latter is defined as sequential application of two PKs, *transport* and *chemistry*.
Process kernel *reactive transport* can susbcycle with respect to *flow*.
Process kernel *chemistry* can susbcycle with respect to *transport*.


Time period control
-------------------

A set of times that simulation hits exactly can be used to avoid problems with
sudden change of boundary conditions or source/sink terms.
This list must *NOT* include start times for time periods *TP #*.

* `"start times`" [Array(double)] is the list of particular times that we want to hit exactly.

* `"initial time step`" [Array(double)] is the first time step after we hit a special
  time specified above.

* `"maximum time step`" [Array(double)] allows the user to limit the time step between
  two particular times.

.. code-block:: xml

   <ParameterList name="Cycle Driver">  <!-- parent list -->
     <ParameterList name="time period control">
       <Parameter name="start times" type="Array(double)" value="{3.16e+10, 6.32e+10}"/>
       <Parameter name="initial time step" type="Array(double)" value="{100.0, 100.0}"/>
       <Parameter name="maximum time step" type="Array(double)" value="{3.16e+8, 4.0e+17}"/>
     </ParameterList>
   </ParameterList>

Between approximately 1000 [y] and 2000 [y], we limit the maximum time step to 10 [y]. 


Restart from checkpoint data file
---------------------------------

A user may request a restart from a checkpoint data file by including sublist 
*restart*. In this scenario, the code will overwrite data initialized using the input XML file.
The purpose of restart is to continue the simulation that has been terminated before for some reasons,
e.g. because its allocation of time ran out.

The value for the current time and current cycle is read from the checkpoint file.
currently, the checkpoint file does not save the list of observations.

* `"restart`" [list]

  * `"file name`" [string] provides name of the existing checkpoint data file to restart from.

.. code-block:: xml
  
  <ParameterList name="Cycle Driver">  <!-- parent list -->
    <ParameterList name="restart">
      <Parameter name="file name" type="string" value="CHECK00123.h5"/>
    </ParameterList>
  </ParameterList>

In this example, Amanzi is restarted with all state data initialized from file
CHECK00123.h5. 


State
=====

List *State* allows the user to initialize various fields and field evaluators 
using a variety of tools. 
A field evaluator is a node in the Phalanx-like (acyclic) dependency tree. 
The corresponding sublist of *State* is named *field evaluators*.
The initialization sublist of *State* is named *initial conditions*.

* `"initialization filename`" [string] (optional) provides name of the existing checkpoint data 
  file. The initialization sequence is as follows. First, we try to initialize a
  field using the provided check-point file. Second, regardless of the outcome of the
  previous step, we try to initialize the field using the sublist `"initial conditions`".
  By design, the second step allows us to overwrite only part for the field.

.. code-block:: xml

  <ParameterList>  <!-- parent list -->
    <ParameterList name="State">
      <Parameter name="initialization filename" type="string" value="CHECK00123.h5"/>
      <ParameterList name="field evaluators">
         ... list of field evaluators
      </ParameterList>
      <ParameterList name="initial conditions">
         ... initialization of fields
      </ParameterList>
    </ParameterList>
  </ParameterList>


Field evaluators
----------------

There are three different types of field evaluators.

Independent field evaluator
...........................

An independent field evaluator has no dependencies and is specified by a function.
Typically, it is evaluated once per simulation.
The evaluator has the following fields.

* `"field evaluator type`" [string] The value of this parameter is used by the factory
  of evaluators. The available option are `"independent variable`", `"primary variable`",
  and `"secondary variable`".

* `"function`" [list] defines a piecewise function for calculating the independent variable.
  In may contain multiple sublists `"DOMAIN`" with identical structure.
  
  * `"DOMAIN`" [list]

    * `"region`" [string] specifies domain of the function, a single region.

    * `"regions`" [Array(string)] is the alternative to option `"region`", domain on 
      the function consists of many regions.

    * `"component`" [string] specifies geometric object associated with the mesh function.
      Available options are `"cell`", `"face`", and `"node`".

    * `"function`" [list] defines an analytic function for calculation. Its structure
      is described in the separate section below.

.. code-block:: xml

  <ParameterList name="field_evaluators">  <!-- parent list -->
    <ParameterList name="SATURATION_LIQUID">
      <Parameter name="field evaluator type" type="string" value="independent variable"/>
      <ParameterList name="function">
        <ParameterList name="DOMAIN">
          <Parameter name="region" type="string" value="COMPUTATIONAL DOMAIN"/>
          <Parameter name="component" type="string" value="cell"/>
          <ParameterList name="function">
            <ParameterList name="function-constant">
              <Parameter name="value" type="double" value="0.8"/>
            </ParameterList>
          </ParameterList>
        </ParameterList>
      </ParameterList>
      <ParameterList name="VerboseObject">
        <Parameter name="Verbosity Level" type="string" value="extreme"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>

The independent variable *SATURATION_LIQUID* is defined as a cell-based variable with
constant value 0.8. 
Note that the user-defined name for this field cannot have spaces.


Primary field evaluator
.......................

The primary field evaluator has no dependencies solved for by a PK.
Examples of independent field evaluators are primary variable of PDEs, such as
pressure and temperature.
Typically this evaluator is used internally to inform the dependency tree about 
a new state of the primary variable.


Secondary field evaluators
..........................

Secondary fields are derived either from primary fields or other secondary fields.
There are two types of secondary fields evaluators.
The first type is used to evaluate a single field.
The second type is used to evaluate efficiently (in one call of an evaluator) multiple fields.

Typically, secondary fields are created by high-level PKs during the setup phase and
inserted automatically in the list of evaluators.
The related XML syntax can provide various parameters needed for evaluation as explained in two
examples below.
The developer can create a secondary field evaluator using common parameters as well
as custom parameters (see the examples).

* `"evaluator dependencies`" [Array(string)] provides a list of fields on which this evaluator
  depends.

* `"check derivatives`" [bool] allows the develop to check derivatives with finite differences.
  This is the expensive option involving finite difference approximations and is recommended for
  code debugging only. Default is *false*.

* `"finite difference epsilon`" [double] defines the finite difference epsilon.
  Default is 1e-10.

.. code-block:: xml

  <ParameterList name="field_evaluators">  <!-- parent list -->
    <ParameterList name="molar_density_liquid">
      <Parameter name="field evaluator type" type="string" value="eos"/>
      <Parameter name="eos basis" type="string" value="both"/>
      <Parameter name="molar density key" type="string" value="molar_density_liquid"/>
      <Parameter name="mass density key" type="string" value="mass_density_liquid"/>
      <ParameterList name="EOS parameters">
        <Parameter name="eos type" type="string" value="liquid water"/>
      </ParameterList>
      <ParameterList name="VerboseObject">
        <Parameter name="Verbosity Level" type="string" value="extreme"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>

In this example the molar density of liquid is evaluated using the *eos* evaluator.
The secondary field name is *molar_density_liquid*.
It is evaluated simultaneously with the secondary field *mass_density_liquid*.
The internal *eos* evaluator knows that these fields depend on fields *temperature* 
and *pressure*; hence, this information is not provided in the input list.
The *eos* evaluator requires one-parameter list to select the proper model for evaluation.

.. code-block:: xml

  <ParameterList name="field_evaluators">  <!-- parent list -->
    <ParameterList name="internal_energy_rock">
      <Parameter name="field evaluator type" type="string" value="iem"/>
      <Parameter name="internal energy key" type="string" value="internal_energy_rock"/>
      <ParameterList name="IEM parameters">
        <Parameter name="iem type" type="string" value="linear"/>
        <Parameter name="heat capacity [J/kg-K]" type="double" value="620.0"/>
      </ParameterList>
      <ParameterList name="VerboseObject">
        <Parameter name="Verbosity Level" type="string" value="extreme"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>

In this example, the secondary field *internal_energy_rock* is evaluated using one of the 
internal *iem* evaluators. 
A particular evaluator is selected dynamically using parameter *iem type*.


Initial conditions
------------------

Constant scalar field
.....................

A constant field is the global (with respect to the mesh) constant. 
At the moment, the set of such fields includes *fluid_density*
and *fluid_viscosity*.
The initialization requires to provide a named sublist with a single
parameter *value*.

.. code-block:: xml

   <ParameterList name="initial conditions">  <!-- parent list -->
     <ParameterList name="fluid_density">
       <Parameter name="value" type="double" value="998.0"/>
     </ParameterList>
   </ParameterList>


Constant vector field
.....................

A constant vector field is the global (with respect to the mesh) vector constant. 
At the moment, the set of such vector constants includes *gravity*.
The initialization requires to provide a named sublist with a single
parameter *value*. In two dimensions, is looks like

.. code-block:: xml

   <ParameterList name="initial conditions">  <!-- parent list -->
     <ParameterList name="gravity">
       <Parameter name="value" type="Array(double)" value="{0.0, -9.81}"/>
     </ParameterList>
   </ParameterList>


A scalar field
..............

A variable scalar field is defined by a few functions (labeled *MESH BLOCK #* in our
example) with non-overlapping domains. 
The required parameters for each function are *region*, *component*,
and the function itself.

* `"regions`" [Array(string)] is list of mesh regions where the function
  should be applied, the domain of the function.

* `"component`" [string] specifies a mesh object on which the discrete field 
  is defined.

Optional parameters are *write checkpoint* and *write vis*.
These parameters define whether the field has to be written into
checkpoints of visualization files. Default values are *true*.

.. code-block:: xml

   <ParameterList name="initial conditions">  <!-- parent list -->
     <ParameterList name="pressure"> 
       <Parameter name="write checkpoint" type="bool" value ="false">   
       <Parameter name="write vis" type="bool" value ="true">
       <ParameterList name="function">
         <ParameterList name="MESH BLOCK 1">
           <Parameter name="regions" type="Array(string)" value="DOMAIN 1"/>
           <Parameter name="component" type="string" value="cell"/>
           <ParameterList name="function">
             <ParameterList name="function-constant">
               <Parameter name="value" type="double" value="90000.0"/>
             </ParameterList>
           </ParameterList>
         </ParameterList>
         <ParameterList name="MESH BLOCK 2">
           ... 
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>

In this example, the field *pressure* has constant value 90000 [Pa] in 
each mesh cell of region *DOMAIN 1*. The second mesh block will define
the pressure in the second mesh region and so on.


A vector or tensor field
........................

A variable tensor (or vector) field is defined similarly to a variable scalar field. 
The difference lies in the definition of the function which is now a multi-value function.

* `"Number of DoFs`" [int] is the number of components in the vector or tensor.

* `"Function type`" [string] defines the function type. The only available option 
  is `"composite function`".

* `"dot with normal`" [bool] triggers the special initialization of a
  vector field such as the `"darcy_flux`". This field is defined by
  projection of the velocity (a vector field) on face normals.
  Changing value to *false* will produce the vector field.

Optional parameters are *write checkpoint*,  and *write vis*.
These parameters define whether the field has to be written into
checkpoints of vis files. Default values are *true*.

.. code-block:: xml

   <ParameterList name="initial conditions">  <!-- parent list -->
     <ParameterList name="darcy_flux">
       <Parameter name="write checkpoint" type="bool" value="true"/>
       <Parameter name="write vis" type="bool" value="false"/>
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
   </ParameterList>

In this example the constant Darcy velocity (0.002, 0.001) [m/s] is dotted with the face 
normal producing one number per mesh face.


Geochemical constraint
......................

We can define geochemical contraint as follows:

.. code-block:: xml

   <ParameterList name="initial conditions">  <!-- parent list -->
     <ParameterList name="geochemical conditions">
       <ParameterList name="initial">
         <Parameter name="regions" type="Array(string)" value="{Entire Domain}"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>

Mesh partitioning
-----------------

Amanzi's state has a number of tools to verify completeness of initial data.
This is done using list *mesh partitions*. 
Each sublist in there must have parameter *region list* specifying
regions that define unique partition of the mesh.

.. code-block:: xml

   <ParameterList name="State">  <!-- parent list -->
     <ParameterList name="mesh partitions">
       <ParameterList name="MATERIALS">
         <Parameter name="region list" type="Array(string)" value="{region1, region2, region3}"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>

In this example, we verify that three mesh regions cover completely the mesh without overlaps.
If so, all material fields, such as *porosity*, will be initialized properly.


Initialization from a file
--------------------------

Some fields can be initialized from files. 
For each field, an additional sublist has to be added to the
named sublist of *State* list with the file name and the name of an attribute. 
For a serial run, the file extension must be *.exo*. 
For a parallel run, it must be *.par*.

.. code-block:: xml

   <ParameterList name="initial conditions">  <!-- parent list -->
     <ParameterList name="permeability">
       <ParameterList name="exodus file initialization">
         <Parameter name="file" type="string" value="mesh_with_data.exo"/>
         <Parameter name="attribute" type="string" value="perm"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>


Example
-------

The complete example of a state initialization is below. Note that
*MATERIAL 1* and *MATERIAL 2* must be labels of the existing regions.

.. code-block:: xml

  <ParameterList name="state">
    <ParameterList name="field evaluators">
      <ParameterList name="porosity">
        <ParameterList name="function">
          <ParameterList name="ALL">
            <Parameter name="regions" type="Array(string)" value="{COMPUTATIONAL DOMAIN}"/>
            <Parameter name="component" type="string" value="cell"/>
            <ParameterList name="function">
              <ParameterList name="function-constant">
                <Parameter name="value" type="double" value="0.408"/>
              </ParameterList>
            </ParameterList>
          </ParameterList>
        </ParameterList>
      </ParameterList>
    </ParameterList>

    <ParameterList name="initial conditions">
      <ParameterList name="fluid_density">
        <Parameter name="value" type="double" value="998.0"/>
      </ParameterList>

      <ParameterList name="gravity">
        <Parameter name="value" type="Array(double)" value="{0.0, -9.81}"/>
      </ParameterList>

      <ParameterList name="pressure">
        <ParameterList name="function">
          <ParameterList name="domain">
            <Parameter name="regions" type="Array(string)" value="{COMPUTATIONAL DOMAIN}"/>
            <Parameter name="component" type="string" value="cell"/>
            <ParameterList name="function">
              <ParameterList name="function-constant">
                <Parameter name="value" type="double" value="90000.0"/>
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


Process kernels (PKs)
=====================

The process kernels list describes all PKs used in a simulation.
The name of the PKs in this list must match *PKNAMEs* in *Cycle Driver* list.

.. code-block:: xml

  <ParameterList>  <!-- parent list -->
    <ParameterList name="PKs">
      <ParameterList name="FLOW and TRANSPORT">
        <Parameter name="PK type" type="string" value="flow reactive transport"/>      
        <Parameter name="PKs order" type="Array(string)" value="{FLOW, TRANSPORT}"/> 
        <Parameter name="master PK index" type="int" value="0"/>
      </ParameterList>
      <ParameterList name="FLOW">
        ...
      </ParameterList>
      <ParameterList name="TRANSPORT">
        ...
      </ParameterList>
    </ParameterList>
  </ParameterList>


Flow PK
-------

The conceptual PDE model for the fully saturated flow is

.. math::
  \phi (s_s + s_y) \frac{\partial p_l}{\partial t} 
  =
  \boldsymbol{\nabla} \cdot (\rho_l \boldsymbol{q}_l) + Q,
  \quad
  \boldsymbol{q}_l 
  = -\frac{\boldsymbol{K}}{\mu} 
  (\boldsymbol{\nabla} p - \rho_l \boldsymbol{g}),

where 
:math:`\phi` is porosity,
:math:`s_s` and :math:`s_y` are specific storage and specific yield, respectively,
:math:`\rho_l` is fluid density,
:math:`Q` is source or sink term,
:math:`\boldsymbol{q}_l` is the Darcy velocity,
and :math:`\boldsymbol{g}` is gravity.

The conceptual PDE model for the partially saturated flow is

.. math::
  \frac{\partial \theta}{\partial t} 
  =
  \boldsymbol{\nabla} \cdot (\eta_l \boldsymbol{q}_l) + Q,
  \quad
  \boldsymbol{q}_l 
  = -\frac{\boldsymbol{K} k_r}{\mu} 
  (\boldsymbol{\nabla} p - \rho_l \boldsymbol{g})

where 
:math:`\theta` is total water content,
:math:`\eta_l` is molar density of liquid,
:math:`\rho_l` is fluid density,
:math:`Q` is source or sink term,
:math:`\boldsymbol{q}_l` is the Darcy velocity,
:math:`k_r` is relative permeability,
and :math:`\boldsymbol{g}` is gravity.
We define 

.. math::
  \theta = \phi \eta_l s_l

where :math:`s_l` is liquid saturation,
and :math:`\phi` is porosity.

Based on these two models, the flow sublist includes exactly one sublist, either 
*Darcy problem* or *Richards problem*.
Structure of both sublists is quite similar.

.. code-block:: xml

   <ParameterList name="Flow">  <!-- parent list -->
     <ParameterList name="Richards problem">
       ...
     </ParameterList>
   </ParameterList>


Physical models and assumptions
...............................

This list is used to summarize physical models and assumptions, such as
coupling with other PKs.
This list is often generated or extended by a high-level MPC PK.

In the code development, this list plays a two-fold role. 
First, it provides necessary information for coupling different PKs such 
as flags for adding a vapor diffusion to Richards' equations.
Second, developers may use it instead of a factory of evaluators such as
creation of primary and secondary evaluators for rock porosity models.
Combination of both approaches may lead to a more efficient code.

* `"vapor diffusion`" [bool] is set up automatically by a high-level PK,
  e.g. by EnergyFlow PK. The default value is `"false`".

* `"multiscale model`" [string] specifies a multiscale model.
  Available options are `"single porosity`" (default) and `"dual porosity`".

* `"water content model`" [string] changes the evaluator for water
  content. Available options are `"generic`" and `"constant density`" (default).

* `"viscosity model`" [string] changes the evaluator for liquid viscosity.
  Available options are `"generic`" and `"constant viscosity`" (default).

* `"porosity model`" [string] specifies an isothermal porosity model.
  Available options are `"compressible: storativity coefficient`",
  `"compressible: pressure function`", and `"constant porosity`" (default).

.. code-block:: xml

   <ParameterList name="Richards problem">  <!-- parent list -->
     <ParameterList name="physical models and assumptions">
       <Parameter name="vapor diffusion" type="bool" value="false"/>
       <Parameter name="water content model" type="string" value="constant density"/>
       <Parameter name="viscosity model" type="string" value="constant viscosity"/>
       <Parameter name="porosity model" type="string" value="compressible: pressure function"/>
       <Parameter name="multiscale model" type="string" value="single porosity"/>
     </ParameterList>
   </ParameterList>


Water retention models
......................

User defines water retention models in sublist *water retention models*. 
It contains as many sublists, e.g. *SOIL_1*, *SOIL_2*, etc, as there are different soils. 
This list is required for *Richards problem* only.
 
The water retention models are associated with non-overlapping regions. Each of the sublists (e.g. *Soil 1*) 
includes a few mandatory parameters: region name, model name, and parameters for the selected model.

* `"water retention model`" [string] specifies a model for the soil.
  The available models are `"van Genuchten`", `"Brooks Corey`", and `"fake`". 
  The later is used only to set up a simple analytic solution for convergence study. 

  * The model `"van Genuchten`" requires `"van Genuchten alpha`" [double],
    `"van Genuchten m`" [double], `"van Genuchten l`" [double], `"residual saturation`" [double],
    and `"relative permeability model`" [string].

  * The model `"Brooks-Corey`" requires `"Brooks Corey lambda`" [double], `"Brooks Corey alpha`" [double],
    `"Brooks Corey l`" [double], `"residual saturation`" [double],
    and `"relative permeability model`" [string].

* `"relative permeability model`" [string] The available options are `"Mualem`" (default) 
  and `"Burdine`".

* `"regularization interval`" [double] removes the kink in the water retention curve at the
  saturation point using a cubic spline. The parameter specifies the regularization region [Pa].
  Default value is 0.

Amanzi performs rudimentary checks of validity of the provided parameters. 
The relative permeability curves can be calculated and saved in an ASCI file 
if the list *output* is provided. This list has two mandatory parameters:

* `"file`" [string] is the user defined file name. It should be different for 
  each soil. 

* `"number of points`" [int] is the number of data points. 
  Each file will contain a table with three columns: saturation, relative permeability, and
  capillary pressure. The data points are equidistributed between the residual saturation
  and 1.

.. code-block:: xml

   <ParameterList name="Richards problem">  <!-- parent list -->
     <ParameterList name="water retention models">
       <ParameterList name="SOIL_1">
         <Parameter name="region" type="string" value="TOP HALF"/>
         <Parameter name="water retention model" type="string" value="van Genuchten"/>
         <Parameter name="van Genuchten alpha" type="double" value="0.000194"/>
         <Parameter name="van Genuchten m" type="double" value="0.28571"/>
         <Parameter name="van Genuchten l" type="double" value="0.5"/>
         <Parameter name="residual saturation" type="double" value="0.103"/>
         <Parameter name="regularization interval" type="double" value="100.0"/>
         <Parameter name="relative permeability model" type="string" value="Mualem"/>
         <ParameterList name="output">
           <Parameter name="file" type="string" value="soil1.txt"/>
           <Parameter name="number of points" type="int" value="1000"/>
         </ParameterList>
       </ParameterList>

       <ParameterList name="SOIL_2">
         <Parameter name="region" type="string" value="BOTTOM HALF"/>
         <Parameter name="water retention model" type="string" value="Brooks Corey"/>
         <Parameter name="Brooks Corey lambda" type="double" value="0.0014"/>
         <Parameter name="Brooks Corey alpha" type="double" value="0.000194"/>
         <Parameter name="Brooks Corey l" type="double" value="0.51"/>
         <Parameter name="residual saturation" type="double" value="0.103"/>
         <Parameter name="regularization interval" type="double" value="0.0"/>
         <Parameter name="relative permeability model" type="string" value="Burdine"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>

In this example, we define two different water retention models in two soils.


Porosity models
...............

User defines porosity models in sublist *porosity models*. 
It contains as many sublists, e.g. *SOIL_1*, *SOIL_2*, etc, as there are different soils. 

The porosity models are associated with non-overlapping regions. Each of the sublists (e.g. *Soil 1*) 
includes a few mandatory parameters: *region name*, *model name*, and parameters for the selected model.

* `"porosity model`" [string] specifies a model for the soil.
  The available models are `"compressible`" and `"constant`". 

  * The model `"compressible`" requires `"underformed soil porosity"`" [double],
    `"reference pressure`" [double], and `"pore compressibility`" [string] [Pa^-1].
    Default value for `"reference pressure`" is 101325.0 [Pa].

  * The model `"constant`" requires `"value`" [double].

.. code-block:: xml

   <ParameterList name="Richards problem">  <!-- parent list -->
     <ParameterList name="porosity models">
       <ParameterList name="SOIL_1">
         <Parameter name="region" type="string" value="TOP HALF"/>
         <Parameter name="porosity model" type="string" value="constant"/>
         <Parameter name="value" type="double" value="0.2"/>
       </ParameterList>

       <ParameterList name="SOIL_2">
         <Parameter name="region" type="string" value="BOTTOM HALF"/>
         <Parameter name="porosity model" type="string" value="compressible"/>
         <Parameter name="underformed soil porosity" type="double" value="0.2"/>
         <Parameter name="reference pressure" type="double" value="101325.0"/>
         <Parameter name="pore compressibility" type="double" value="1e-8"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>

In this example, we define two different porosity models in two soils.


Multiscale continuum models
...........................

The list *multiscale models* is the placeholder for various multiscale models.
The list is extension of the list *water retention models*. 
Its ordered by soil regions and includes parameters for the multiscale,
capillary pressure, and relative permebility models.
This list is optional. 

* `"multiscale model`" [string] is the model name. Available option is `"dual porosity`".

* `"mass transfer coefficient`" [double] is the mass transfer coefficient.

* `"tolerance`" [double] defines tolerance for iterative methods used to solve
  secondary equations. Default is 1e-8.

* `"water retention model`" [string] specifies a model for the soil.
  The available models are `"van Genuchten`" and `"Brooks Corey`". 
  Parameters for each model are described above.

* `"relative permeability model`" [string] The available options are `"Mualem`" (default) 
  and `"Burdine`".

.. code-block:: xml

   <ParameterList name="Richards problem">  <!-- parent list -->
     <ParameterList name="multiscale models"> 
       <ParameterList name="SOIL_1">
         <Parameter name="region" type="string" value="TOP HALF"/>
         <Parameter name="multiscale model" type="string" value="dual porosity"/> 
         <Paramater name="mass transfer coefficient" type="double" value="4.0e-5"/>
         <Paramater name="tolerance" type="double" value="1e-8"/>

         <Parameter name="water retention model" type="string" value="van Genuchten"/>
         <Parameter name="van Genuchten alpha" type="double" value="0.000194"/>
         <Parameter name="van Genuchten m" type="double" value="0.28571"/>
         <Parameter name="van Genuchten l" type="double" value="0.5"/>
         <Parameter name="residual saturation" type="double" value="0.103"/>
         <Parameter name="relative permeability model" type="string" value="Mualem"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>


Upwind 
......

This section discusses interface treatment of cell-centered fields such as 
relative permeability, density and viscosity.

* `"upwind`" [list] collects information required for treatment of
  relative permeability, density and viscosity on mesh faces.

  * `"relative permeability`" [string] defines a method for calculating the *upwinded* 
    relative permeability. The available options are: `"upwind: gravity`", 
    `"upwind: darcy velocity`" (default), `"upwind: amanzi``", 
    `"other: harmonic average`", and `"other: arithmetic average`".

  * `"upwind update`" [string] defines frequency of recalculating Darcy flux inside
    nonlinear solver. The available options are `"every time step`" and `"every nonlinear iteration`".
    The first option freezes the Darcy flux for the whole time step. The second option
    updates it on each iteration of a nonlinear solver. The second option is recommended
    for the New ton solver. It may impact significantly upwinding of the relative permeability 
    and convergence rate of this solver.

  * `"upwind method`" [string] specifies a method for treating nonlinear diffusion coefficient.
    Available options are `"standard`", `"divk`" (default), and `"second-order`" (experimental). 

  * `"upwind NAME parameters`" [list] defines parameters for upwind method `"NAME`".

    * `"tolerance`" [double] specifies relative tolerance for almost zero local flux. In such
      a case the flow is assumed to be parallel to a mesh face. Default value is 1e-12.

    * [WIP] `"reconstruction method`" [string] defines a reconstruction method for the second-order upwind.

    * [WIP] `"limiting method`" [string] defines limiting method for the second-order upwind.

.. code-block:: xml

   <ParameterList name="Richards problem">  <!-- parent list -->
     <ParameterList name="upwind">
       <Parameter name="relative permeability" type="string" value="upwind with Darcy flux"/>
       <Parameter name="upwind update" type="string" value="every timestep"/>

       <Parameter name="upwind method" type="string" value="standard"/>
       <ParameterList name="upwind standard parameters">
          <Parameter name="tolerance" type="double" value="1e-12"/>
       </ParameterList>
     </ParameterList>  
   </ParameterList>  


Diffusion operators
...................

List *operators* describes the PDE structure of the flow, specifies a discretization
scheme, and selects assembling schemas for matrices and preconditioners.

* `"operators`" [list] 

  * `"diffusion operator`" [list] defines parameters for generating and assembling diffusion matrix.

    * `"matrix`" [list] defines parameters for generating and assembling diffusion matrix. See section
      describing operators. 
      When `"Richards problem`" is selected, Flow PK sets up proper value for parameter `"upwind method`" of 
      this sublist.

    * `"preconditioner`" [list] defines parameters for generating and assembling diffusion 
      matrix that is used to create preconditioner. 
      This sublist is ignored inside sublist `"Darcy problem`".
      Since update of preconditioner can be lagged, we need two objects called `"matrix`" and `"preconditioner`".
      When `"Richards problem`" is selected, Flow PK sets up proper value for parameter `"upwind method`" of 
      this sublist.

.. code-block:: xml

  <ParameterList name="Richards problem">  <!-- parent list -->
    <ParameterList name="operators">
      <ParameterList name="diffusion operator">
        <ParameterList name="matrix">
          <Parameter name="discretization primary" type="string" value="mfd: optimized for monotonicity"/>
          <Parameter name="discretization secondary" type="string" value="mfd: optimized for sparsity"/>
          <Parameter name="schema" type="Array(string)" value="{face, cell}"/>
          <Parameter name="preconditioner schema" type="Array(string)" value="{face}"/>
          <Parameter name="gravity" type="bool" value="true"/>
          <Parameter name="gravity term discretization" type="string" value="hydraulic head"/>
        </ParameterList>
        <ParameterList name="preconditioner">
          <Parameter name="discretization primary" type="string" value="mfd: optimized for monotonicity"/>
          <Parameter name="discretization secondary" type="string" value="mfd: optimized for sparsity"/>
          <Parameter name="schema" type="Array(string)" value="{face, cell}"/>
          <Parameter name="preconditioner schema" type="Array(string)" value="{face}"/>
          <Parameter name="newton correction" type="string" value="approximate jacobian"/>
        </ParameterList>
      </ParameterList>
    </ParameterList>
  </ParameterList>

This example creates a p-lambda system, i.e. the pressure is
discretized in mesh cells and on mesh faces. 
The preconditioner is defined on faces only, i.e. cell-based unknowns
are eliminated explicitly and the preconditioner is applied to the
Schur complement.


Boundary conditions
...................

Boundary conditions are defined in sublist *boundary conditions*. 
Four types of boundary conditions are supported.
Each type has a similar structure: a list of identical elements that contain
information about a part of the boundary where it is prescribed, a function
to calculate it, and optional parameters to modify it slightly.
This modification is referred to as a submodel and requires additional parameters as described below. 

* `"pressure`" [list] is the Dirichlet boundary condition where the pressure is prescribed on a part of the 
  boundary surface. No submodels is available.

* `"mass flux`" [list] is the Neumann boundary condition where an outward mass flux is prescribed on a 
  part of the boundary surface.
  This is the default boundary condition. If no condition is specified on a mesh face, the zero flux 
  boundary condition is used. 

  * `"rainfall`" [bool] indicates the submodel where the mass flux is defined with respect to the gravity 
    vector and the actual flux depends on the boundary slope. Default is `"false`".

* `"static head`" [list] is the Dirichlet boundary condition where the hydrostatic pressure is prescribed 
  on a part of the boundary.

  * `"relative to top`" [bool] indicates the submodel where the static head is defined with respect
    to the top boundary (a curve in 3D) of the specified regions. Support of 2D is turned off.
    Default value is `"false`". 

  * `"no flow above water table`" [bool] indicates the submodel where the no-flow boundary condition 
    has to be used above the water table. This switch uses the pressure value at a face
    centroid. Default is `"false`".

* `"seepage face`" [list] is the seepage face boundary condition, a dynamic combination of the `"pressure`" and 
  `"mass flux`" boundary conditions over the specified region. 
  The atmospheric pressure is prescribed if internal pressure is higher it. 
  Otherwise, the outward mass flux is prescribed. 

  * `"reference pressure`" [double] defaults to the atmospheric pressure. 

  * `"rainfall`" [bool] indicates the submodel where the mass flux is defined with respect to the gravity 
    vector and the actual influx depends on the boundary slope. Default is `"false`".

  * `"submodel`" [string] indicates different models for the seepage face boundary condition.
    It can take values `"PFloTran`" and `"FACT`". The first option leads to a 
    discontinuous change of the boundary condition type from the infiltration to pressure. 
    The second option is described in the document on mathematical models. 
    It employs a smooth transition from the infiltration 
    to mixed boundary condition. Default is `"PFloTran`".

.. code-block:: xml

   <ParameterList name="Richards problem">  <!-- parent list -->
     <ParameterList name="boundary conditions">
       <ParameterList name="pressure">
         <ParameterList name="BC 0">
           <Parameter name="regions" type="Array(string)" value="{WEST_SIDE}"/>
           <ParameterList name="boundary pressure">
             <ParameterList name="function-constant">
               <Parameter name="value" type="double" value="101325.0"/>
             </ParameterList>
           </ParameterList>
         </ParameterList>
       </ParameterList>

       <ParameterList name="mass flux">
         <ParameterList name="BC 1">
           <Parameter name="regions" type="Array(string)" value="{NORTH_SIDE, SOUTH_SIDE}"/>
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
           <Parameter name="regions" type="Array(string)" value="{EAST_SIDE}"/>
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
           <Parameter name="regions" type="Array(string)" value="{EAST_SIDE_BOTTOM}"/>
           <Parameter name="rainfall" type="bool" value="true"/>
           <Parameter name="submodel" type="string" value="PFloTran"/>
           <ParameterList name="outward mass flux">
             <ParameterList name="function-constant">
               <Parameter name="value" type="double" value="1.0"/>
             </ParameterList>
           </ParameterList>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>

This example includes all four types of boundary conditions. The boundary of a square domain 
is split into six pieces. Constant function is used for simplicity and can be replaced by any
of the other available functions.


Sources and sinks
.................

The sources and sinks are typically associated with pumping wells. The structure
of list *source terms* mimics that of list *boundary conditions*. 
Again, constant functions can be replaced by any of the available functions.

* `"regions`" [Array(string)] is the list of regions where the source is defined.

* `"spatial distribution method`" [string] is the method for distributing
  source Q over the specified regions. The available options are `"volume`",
  `"none`", and `"permeability`". For option `"none`", the source term Q is measured
  in [kg/m^3/s]. For the other options, it is measured in [kg/s]. When the source function
  is defined over a few regions, Q is distributed independently over each region.
  Default is `"none`".

* `"submodel`" [string] refines definition of the source. Available options are `"rate`"
  and `"integrated source`". The first option defines the source in a natural way as the rate 
  of change `q`. The second option defines the indefinite integral `Q` of the rate 
  of change, i.e. the source term is calculated as `q = dQ/dt`. Default is `"rate`".

.. code-block:: xml

   <ParameterList name="Richards problem">  <!-- parent list -->
     <ParameterList name="source terms">
       <ParameterList name="SRC 0">
         <Parameter name="regions" type="Array(string)" value="{WELL_EAST}"/>
         <Parameter name="spatial distribution method" type="string" value="volume"/>
         <Parameter name="submodel" type="string" value="rate"/>
         <ParameterList name="sink">
           <ParameterList name="function-constant">
             <Parameter name="value" type="double" value="-0.1"/>
           </ParameterList>
         </ParameterList>
       </ParameterList>

       <ParameterList name="SRC 1">
         <Parameter name="regions" type="Array(string)" value="{WELL_WEST}"/>
         <Parameter name="spatial distribution method" type="string" value="permeability"/>
         <ParameterList name="sink">
           <ParameterList name="function-constant">
             <Parameter name="value" type="double" value="-0.2"/>
           </ParameterList>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>


Time integrator
...............

The list *time integrator* defines a generic time integrator used
by the cycle driver. 
This driver assumes that each PK has only one time integrator.
The list *time integrator* defines parameters controlling linear and 
nonlinear solvers during a time integration period.
We break this long sublist into smaller parts. 


Initialization and constraints
``````````````````````````````

* `"error control options`" [Array(string)] lists various error control options. 
  A nonlinear solver is terminated when all listed options are passed. 
  The available options are `"pressure`", `"saturation`", and `"residual`". 
  All errors are relative, i.e. dimensionless. 
  The error in pressure is compared with capillary pressure plus atmospheric pressure. 
  The other two errors are compared with 1. 
  The option `"pressure`" is always active during steady-state time integration.
  The option  `"saturation`" is always active during transient time integration.

* `"linear solver`" [string] refers to a generic linear solver from list `"Solvers`".
  It is used in all cases except for `"initialization`" and `"enforce pressure-lambda constraints`".

* `"preconditioner`" [string] specifies preconditioner for linear and nonlinear solvers.

* `"preconditioner enhancement`" [string] specifies a linear solver that binds 
  the above preconditioner to improve spectral properties. Default is `"none`".

* `"initialization`" [list] defines parameters for calculating initial pressure guess.
  It can be used to obtain pressure field which is consistent with the boundary conditions.
  Default is empty list.

  * `"method`" [string] specifies an optional initialization methods. The available 
    options are `"picard`" and `"saturated solver`". The latter option leads to solving 
    a Darcy problem. The former option uses sublist `"picard parameters`".
    *Picard works better if a bounded initial pressure guess is provided.* 

  * `"picard parameters`" [list] defines control parameters for the Picard solver.

    * `"convergence tolerance`" [double] specifies nonlinear convergence tolerance. 
      Default is 1e-8.
    * `"maximum number of iterations`" [int] limits the number of iterations. Default is 400. 

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

   <ParameterList name="Richards problem">  <!-- parent list -->
     <ParameterList name="time integrator">
       <Parameter name="error control options" type="Array(string)" value="{pressure, saturation}"/>
       <Parameter name="linear solver" type="string" value="GMRES_with_AMG"/>
       <Parameter name="linear solver as preconditioner" type="string" value="GMRES_with_AMG"/>
       <Parameter name="preconditioner" type="string" value="HYPRE_AMG"/>

       <ParameterList name="initialization">  <!-- first method -->
         <Parameter name="method" type="string" value="saturated solver"/>
         <Parameter name="linear solver" type="string" value="PCG_with_AMG"/>
         <Parameter name="clipping pressure value" type="double" value="50000.0"/>
       </ParameterList>

       <ParameterList name="initialization">  <!-- alternative method -->
         <Parameter name="method" type="string" value="picard"/>
         <Parameter name="linear solver" type="string" value="PCG_with_AMG"/>
         <ParameterList name="picard parameters">
           <Parameter name="convergence tolerance" type="double" value="1e-8"/> 
           <Parameter name="maximum number of iterations" type="int" value="20"/> 
         </ParameterList>
       </ParameterList>

       <ParameterList name="pressure-lambda constraints">
         <Parameter name="method" type="string" value="projection"/>
         <Parameter name="inflow krel correction" type="bool" value="false"/>
         <Parameter name="linear solver" type="string" value="PCG_with_AMG"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>


Time step controller and nonlinear solver
`````````````````````````````````````````

The time step is controlled by parameter *time step controller type*
and the related list of options.
Nonlinear solver is controlled by parameter *solver type*  and related list of options.
Amanzi supports a few nonlinear solvers described in details in a separate section.

* `"time step controller type`" [list]
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
    The default value is 1, while a reasonable value may be as large as 1000. 

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

  * `"solver type`" [string] defines nonlinear solver used on each time step for
    a nonlinear algebraic system :math:`F(x) = 0`. 
    The available options `"aa`", `"nka`" and `"Newton`".

  * `"nka parameters`" [list] internal parameters for the nonlinear
    solver NKA.

  * `"aa parameters`" [list] internal parameters for the nonlinear
    solver AA (Anderson acceleration).

.. code-block:: xml

   <ParameterList name="Richards problem">  <!-- parent list -->
     <ParameterList name="time integrator">
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

         <Parameter name="solver type" type="string" value="nka"/>
         <ParameterList name="nka parameters">
           <Parameter name="nonlinear tolerance" type="double" value="1e-5"/>
           <Parameter name="limit iterations" type="int" value="30"/>
           <Parameter name="diverged tolerance" type="double" value="1e+10"/>
           <Parameter name="diverged l2 tolerance" type="double" value="1e+10"/>
           <Parameter name="diverged pc tolerance" type="double" value="1e+10"/>
           <Parameter name="max du growth factor" type="double" value="1e+5"/>
           <Parameter name="max divergent iterations" type="int" value="3"/>
           <Parameter name="max nka vectors" type="int" value="10"/>
           <Parameter name="modify correction" type="bool" value="false"/>
           <ParameterList name="VerboseObject">
             <Parameter name="Verbosity Level" type="string" value="high"/>
           </ParameterList>
         </ParameterList>

         <!-- alternative solver 
         <Parameter name="solver type" type="string" value="aa"/>
         <ParameterList name="aa parameters">
           <Parameter name="nonlinear tolerance" type="double" value="1e-5"/>
           <Parameter name="limit iterations" type="int" value="30"/>
           <Parameter name="diverged tolerance" type="double" value="1e+10"/>
           <Parameter name="diverged l2 tolerance" type="double" value="1e+10"/>
           <Parameter name="diverged pc tolerance" type="double" value="1e+10"/>
           <Parameter name="max du growth factor" type="double" value="1e+5"/>
           <Parameter name="max divergent iterations" type="int" value="3"/>
           <Parameter name="max aa vectors" type="int" value="10"/>
           <Parameter name="modify correction" type="bool" value="false"/>
           <Parameter name="relaxation parameter" type="double" value="0.75"/>
         </ParameterList-->

       </ParameterList>
     </ParameterList>
   </ParameterList>

In this example, the time step is increased by factor 1.2 when the nonlinear
solver converges in 10 or less iterations. 
The time step is not changed when the number of nonlinear iterations is
between 11 and 15.
The time step will be cut twice if the number of nonlinear iterations exceeds 15.


Developer parameters
````````````````````

The remaining parameters in the time integrator sublist include 
those needed for unit tests, and future code development. 

.. code-block:: xml

   <ParameterList name="time integrator">
     <ParameterList name="obsolete parameters">
       <Parameter name="start time" type="double" value="0.0"/>
       <Parameter name="end time" type="double" value="100.0"/>
       <Parameter name="maximum number of iterations" type="int" value="400"/>
       <Parameter name="error abs tol" type="double" value="1"/>
       <Parameter name="error rel tol" type="double" value="0"/>
     </ParameterList>
   </ParameterList>


Other parameters
................

The remaining *Flow* parameters are

* `"atmospheric pressure`" [double] defines the atmospheric pressure, [Pa].

* `"absolute permeability coordinate system`" [string] defines coordinate system
  for calculating absolute permeability. The available options are `"cartesian`"
  and `"layer`".

* `"clipping parameters`"[list] defines how solution increment calculated by a nonlinear 
  solver is modified e.g., clipped.

* `"plot time history`" [bool] produces an ASCII file with the time history. Default is `"false`".

.. code-block:: xml

   <ParameterList name="Richards problem">  <!-- parent list -->
     <ParameterList name="clipping parameters">
        <Parameter name="maximum saturation change" type="double" value="0.25"/>
        <Parameter name="pressure damping factor" type="double" value="0.5"/>
     </ParameterList>	

     <Parameter name="plot time history" type="bool" value="false"/>
   </ParameterList>	


Explanation of verbose output
.............................

When verbosity is set to *high*, this PK reports information about 
current status of the simulation.
Here after keyword *global* refers to the whole simulation including
all time periods, keyword *local* refers to the current time period.
The incomplete list is

 * [global] cycle number, time before the step, and time step dt (in years)
 * [local] step number, time T, and dT inside the time integrator (in seconds)
 * [local] frequency of preconditioner updates
 * [local] number of performed nonlinear steps and value of the nonlinear residual
 * [local] total number of successful time steps (TS), failed time steps (FS),
   preconditioner updates (PC/1) and preconditioner applies (PC/2),
   linear solves insides preconditioner (LS)
 * [local] amount of liquid (water) in the reservoir and amount of water entering
   and living domain through its boundary (based on darcy flux).
 * [global] current simulation time (in years)

.. code-block:: xml

  CycleDriver      |   Cycle 40: time(y) = 0.953452, dt(y) = 0.238395
  TI::BDF1         |    step 40 T = 3.00887e+07 [sec]  dT = 7.52316e+06
  TI::BDF1         |    preconditioner lag is 20 out of 20
  TI::BDF1         |    success: 4 nonlinear itrs error=7.87642e-08
  TI::BDF1         |    TS:40 FS:0 NS:64 PC:42 64 LS:0 dt:1.0000e+03 7.5232e+06
  FlowPK::Richards |    reservoir water mass=1.36211e+06 [kg], total influx=897.175 [kg]
  CycleDriver      |   New time(y) = 1.19185


Transport PK
------------

The conceptual PDE model for the fully saturated flow is

.. math::
  \frac{\partial (\phi s_l C_l)}{\partial t} 
  =
  - \boldsymbol{\nabla} \cdot (\boldsymbol{q}_l C_l) 
  + \boldsymbol{\nabla} \cdot (\phi s_l (\boldsymbol{D}_l + \tau \boldsymbol{M}_l) \boldsymbol{\nabla} C_l) + Q,

where 
:math:`\phi` is porosity,
:math:`s_l` is liquid saturation, 
:math:`Q` is source or sink term,
:math:`\boldsymbol{q}_l` is the Darcy velocity,
:math:`\boldsymbol{D}_l` is dispersion tensor,
:math:`\boldsymbol{M}_l` is diffusion coefficient,
and :math:`\tau` is tortuosity.
For an isotropic medium with no preferred axis of symmetry the dispersion 
tensor has the following form:

.. math::
  \boldsymbol{D}_l 
  = \alpha_T \|\boldsymbol{v}\| \boldsymbol{I} 
  + \left(\alpha_L-\alpha_T \right) 
    \frac{\boldsymbol{v} \boldsymbol{v}}{\|\boldsymbol{v}\|},

where
:math:`\alpha_L` is longitudinal dispersivity,
:math:`\alpha_T` is  transverse dispersivity,
and :math:`\boldsymbol{v}` is average pore velocity.


Physical models and assumptions
...............................

This list is used to summarize physical models and assumptions, such as
coupling with other PKs.
This list is often generated or extended by a high-level MPC PK.

* `"gas diffusion`" [bool] indicates that air-water partitioning coefficients
  are used to distribute components between liquid and as phases. Default is *false*.

* `"permeability field is required`" [bool] indicates if some transport features
  require absolute permeability. Default is *false*.

* `"multiscale model`" [string] specifies a multiscale model.
  Available options are `"single porosity`" (default) and `"dual porosity`".

.. code-block:: xml

   <ParameterList name="Transport">  <!-- parent list -->
     <ParameterList name="physical models and assumptions">
       <Parameter name="gas diffusion" type="bool" value="false"/>
       <Parameter name="permeability field is required" type="bool" value="false"/>
       <Parameter name="multiscale model" type="string" value="single porosity"/>
     </ParameterList>
   </ParameterList>


Global parameters
.................

This list is used to summarize physical models and assumptions, such as
The transport component of Amanzi performs advection of aqueous and gaseous
components and their dispersion and diffusion. 
The main parameters control temporal stability, spatial 
and temporal accuracy, and verbosity:

* `"cfl`" [double] Time step limiter, a number less than 1. Default value is 1.
   
* `"spatial discretization order`" [int] defines accuracy of spatial discretization.
  It allows values 1 or 2. Default value is 1. 
  
* `"temporal discretization order`" [int] defines accuracy of temporal discretization.
  It allows values 1 or 2. Default value is 1.

* `"reconstruction`" [list] collects reconstruction parameters. The available options are
  describe in the separate section below.

* `"solver`" [string] Specifies the dispersion/diffusion solver.

* `"number of aqueous components`" [int] The total number of aqueous components. 
  Default value is the total number of components.

* [WIP] `"number of gaseous components`" [int] The total number of gaseous components. 
  Default value is 0.
   
.. code-block:: xml

  <ParameterList>  <!-- parent list -->
    <ParameterList name="Transport">
      <Parameter name="cfl" type="double" value="1.0"/>
      <Parameter name="spatial discretization order" type="int" value="1"/>
      <Parameter name="temporal discretization order" type="int" value="1"/>
      <Parameter name="solver" type="string" value="PCG_SOLVER"/>

      <ParameterList name="reconstruction">
        <Parameter name="method" type="string" value="cell-based"/>
        <Parameter name="polynomial order" type="int" value="1"/>
        <Parameter name="limiter" type="string" value="tensorial"/>
        <Parameter name="limiter extension for transport" type="bool" value="true"/>
      </ParameterList>

      <ParameterList name="VerboseObject">
        <Parameter name="Verbosity Level" type="string" value="high"/>
      </ParameterList>
    </ParameterList>  
  </ParameterList>  


Material properties
...................

The material properties include dispersivity model and diffusion parameters 
for aqueous and gaseous phases.
The dispersivity is defined as a soil property. 
The diffusivity is defined independently for each solute.

* SOIL [list] Defines material properties.
  
  * `"region`" [Array(string)] Defines geometric regions for material SOIL.
  * `"model`" [string] Defines dispersivity model, choose exactly one of the following: `"scalar`", `"Bear`",
    `"Burnett-Frind`", or `"Lichtner-Kelkar-Robinson`".
  * `"parameters for MODEL`" [list] where `"MODEL`" is the model name.
    For model `"scalar`", the following options must be specified:

      * `"alpha`" [double] defines dispersion in all directions. 

    For model `"Bear`", the following options must be specified:

      * `"alphaL`" [double] defines dispersion in the direction of Darcy velocity.
      * `"alphaT`" [double] defines dispersion in the orthogonal direction.
    
    For model `"Burnett-Frind`", the following options must be specified:

      * `"alphaL`" [double] defines the longitudinal dispersion in the direction of Darcy velocity.
      * `"alphaTH`" [double] Defines the transverse dispersion in the horizonla direction orthogonal directions.
      * `"alphaTV`" [double] Defines dispersion in the orthogonal directions.
        When `"alphaTH`" equals to `"alphaTV`", we obtain dispersion in the direction of the Darcy velocity.
        This and the above parameters must be defined for `"Burnett-Frind`" and `"Lichtner-Kelkar-Robinson`" models.

    For model `"Lichtner-Kelker-Robinson`", the following options must be specified:

      * `"alphaLH`" [double] defines the longitudinal dispersion in the horizontal direction.
      * `"alphaLV`" [double] Defines the longitudinal dispersion in the vertical direction.
        When `"alphaLH`" equals to `"alphaLV`", we obtain dispersion in the direction of the Darcy velocity.
        This and the above parameters must be defined for `"Burnett-Frind`" and `"Lichtner-Kelker-Robinson`" models.
      * `"alphaTH`" [double] Defines the transverse dispersion in the horizontal direction orthogonal directions.
      * `"alphaTV`" [double] Defines dispersion in the orthogonal directions.
        When `"alphaTH`" equals to `"alphaTV`", we obtain dispersion in the direction of the Darcy velocity.
        This and the above parameters must be defined for `"Burnett-Frind`" and `"Lichtner-Kelker-Robinson`" models.

  * `"aqueous tortuosity`" [double] Defines tortuosity for calculating diffusivity of liquid solutes.
  * `"gaseous tortuosity`" [double] Defines tortuosity for calculating diffusivity of gas solutes.
 
Three examples are below:

.. code-block:: xml

  <ParameterList name="Transport">  <!-- parent list -->
    <ParameterList name="material properties">
      <ParameterList name="WHITE SOIL">
        <Parameter name="regions" type="Array(string)" value="{TOP_REGION, BOTTOM_REGION}"/>
        <Parameter name="model" type="string" value="Bear"/>
        <ParameterList name="parameters for Bear">
          <Parameter name="alphaL" type="double" value="1e-2"/>
          <Parameter name="alphaT" type="double" value="1e-5"/>
        <ParameterList>
        <Parameter name="aqueous tortuosity" type="double" value="1.0"/>       
        <Parameter name="gaseous tortuosity" type="double" value="1.0"/>       
      </ParameterList>  
     
      <ParameterList name="GREY SOIL">
        <Parameter name="regions" type="Array(string)" value="{MIDDLE_REGION}"/>
        <Parameter name="model" type="string" value="Burnett-Frind"/>
        <ParameterList name="parameters for Burnett-Frind">
          <Parameter name="alphaL" type="double" value="1e-2"/>
          <Parameter name="alphaTH" type="double" value="1e-3"/>
          <Parameter name="alphaTV" type="double" value="2e-3"/>
        <ParameterList>
        <Parameter name="aqueous tortuosity" type="double" value="0.5"/>
        <Parameter name="gaseous tortuosity" type="double" value="1.0"/>       
      </ParameterList>  
    </ParameterList>  
  </ParameterList>  


* `"molecular diffusion`" [list] defines names of solutes in aqueous and gaseous phases and related
  diffusivity values.

.. code-block:: xml

  <ParameterList name="Transport">  <!-- parent list -->
    <ParameterList name="molecular diffusion">
      <Parameter name="aqueous names" type=Array(string)" value="{CO2(l),Tc-99}"/>
      <Parameter name="aqueous values" type=Array(double)" value="{1e-8,1e-9}"/>

      <Parameter name="gaseous names" type=Array(string)" value="{CO2(g)}"/>
      <Parameter name="gaseous values" type=Array(double)" value="{1e-8}"/>
      <Parameter name="air-water partitioning coefficient" type=Array(double)" value="{0.03}"/> 
    </ParameterList>  
  </ParameterList>  


Dispersion operator
...................

List *operators* describes the PDE structure of the flow, specifies a discretization
scheme, and selects assembling schemas for matrices and preconditioners.

* `"operators`" [list] 

  * `"diffusion operator`" [list] 

    * `"matrix`" [list] defines parameters for generating and assembling dispersion matrix.
      See section describing operators. 

.. code-block:: xml

  <ParameterList name="Transport">  <!-- parent list -->
    <ParameterList name="operators">
      <ParameterList name="diffusion operator">
        <ParameterList name="matrix">
          <Parameter name="discretization primary" type="string" value="mfd: optimized for monotonicity"/>
          <Parameter name="discretization secondary" type="string" value="mfd: two-point flux approximation"/>
          <Parameter name="schema" type="Array(string)" value="{face, cell}"/>
          <Parameter name="preconditioner schema" type="Array(string)" value="{face}"/>
        </ParameterList>
      </ParameterList>
    </ParameterList>
  </ParameterList>

This example creates a p-lambda system, i.e. the concentation is
discretized in mesh cells and on mesh faces. The later unknowns are auxiliary unknwons.


Multiscale continuum models
...........................

The list of multiscale models is the placeholder for various multiscale models.
Its ordered by materials and includes parameters for the assigned multiscale model
This list is optional.

* `"multiscale model`" [string] is the model name. Available option is "dual porosity".

* `"regions`" [Array(string)] is the list of regions where this model should be applied.

* `"solute transfer coefficient`" [list] defines diffusive solute transport due to
  convetration gradients.

.. code-block:: xml

  <ParameterList name="Transport">  <!-- parent list -->
    <ParameterList name="multiscale models">
      <ParameterList name="WHITE SOIL">
        <Parameter name="multiscale model" type="string" value="dual porosity"/>
        <Parameter name="regions" type="Array(string)" value="{TOP_REGION, BOTTOM_REGION}"/>
        <Paramater name="solute transfer coefficient" type="double" value="4.0e-5"/>
      </ParameterList>  

      <ParameterList name="GREY SOIL">
         ...  
      </ParameterList>  
    </ParameterList>  
  </ParameterList>  


Boundary conditions
...................

For the advective transport, the boundary conditions must be specified on inflow parts of the
boundary. If no value is prescribed through the XML input, the zero influx boundary condition
is used. Note that the boundary condition is set up separately for each component.
The structure of boundary conditions is aligned with that used for Flow and
allows us to define spatially variable boundary conditions. 

* `"boundary conditions`" [list]

  * `"concentration`" [list] This is a reserved keyword.
   
    * "COMP" [list] contains a few sublists (e.g. BC_1, BC_2) for boundary conditions.
      The name *COMP* must be the name in the list of solutes.
 
      * "BC_1" [list] defines boundary conditions using arrays of boundary regions and attached
        functions.
   
      * `"regions`" [Array(string)] defines a list of boundary regions where a boundary condition
        must be applied.
      * `"boundary concentration`" [list] defines a function for calculating boundary conditions.
        The function specification is described in subsection Functions.

The example below sets constant boundary condition 1e-5 for the duration of transient simulation.

.. code-block:: xml

  <ParameterList name="Transport">  <!-- parent list -->
    <ParameterList name="boundary conditions">
      <ParameterList name="concentration">
        <ParameterList name="H+"> 
          <ParameterList name="EAST CRIB">   <!-- user defined name -->
            <Parameter name="regions" type="Array(string)" value="{TOP, LEFT}"/>
            <ParameterList name="boundary concentration">
              <ParameterList name="function-constant">  <!-- any time function -->
                <Parameter name="value" type="double" value="1e-5"/>
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
  </ParameterList>


Geochemical boundary conditions are concentration-type boundary conditions
but require special treatment. 
Note that the number of *forms* below is one less than the number of times
and geochemical conditions.

.. code-block:: xml

  <ParameterList name="Transport">  <!-- parent list -->
    <ParameterList name="boundary conditions">
      <ParameterList name="geochemical conditions">
        <ParameterList name="H+"> 
          <ParameterList name="EAST CRIB">   <!-- user defined name -->
            <Parameter name="times" type="Array(double)" value="{0.0, 100.0}"/>
            <Parameter name="geochemical conditions" type="Array(string)" value="{cond1, cond2}"/>
            <Parameter name="time functions" type="Array(string)" value="{constant}"/>
            <Parameter name="regions" type="Array(string)" value="{CRIB1}"/>
          </ParameterList>
        </ParameterList>
      </ParameterList>
    </ParameterList>
  </ParameterList>


Sources and sinks
.................

The external sources are typically located at pumping wells. The structure
of list *source terms* includes only sublists named after components. 
Again, constant functions can be replaced by any available time-function.
Note that the source values are set up separately for each component.

* `"concentration`" [list] This is a reserved keyword.

 * "COMP" [list] contains a few sublists (e.g. SRC_1, SRC_2) for multiple sources and sinks.

  * "SRC_1" [list] defines a source using arrays of domain regions, a function, and 
    a distribution method.
   
   * `"regions`" [Array(string)] defines a list of domain regions where a source term
     must be applied.

   * `"sink`" [list] is a function for calculating a source term.
     The function specification is described in subsection Functions.

    * `"spatial distribution method`" [string] identifies a method for distributing
      source Q over the specified regions. The available options are `"volume`",
      `"none`", and `"permeability`". For option `"none`" the source term Q is measured
      in [mol/m^3/s]. For the other options, it is measured in [mol/s]. When the source function
      is defined over a few regions, Q will be distributed independently over each region.
      Default value is `"none`".

    * `"submodel`" [string] refines definition of source. Available options are `"rate`"
      and `"integrand`". The first option defines rate of change `q`, the second one 
      defines integrand `Q` of a rate `Q = dq/dt`. Default is `"rate`".

This example defines one well and one sink.

.. code-block:: xml

   <ParameterList name="Transport">  <!-- parent list -->
     <ParameterList name="source terms">
       <ParameterList name="concentration">
         <ParameterList name="H+"> 
           <ParameterList name="SOURCE: EAST WELL">   <!-- user defined name -->
	     <Parameter name="regions" type="Array(string)" value="{EAST_WELL}"/>
             <Parameter name="spatial distribution method" type="string" value="volume"/>
             <Parameter name="submodel" type="string" value="rate"/>
             <ParameterList name="sink">   <!-- keyword, do not change -->
               <ParameterList name="function-constant">
                 <Parameter name="value" type="double" value="-0.01"/>
               </ParameterList>
             </ParameterList>
           </ParameterList>
           <ParameterList name="source for west well">
              ...
           </ParameterList>
         </ParameterList>
     
         <ParameterList name="CO2(g)">   <!-- next component, a gas -->
           <ParameterList name="SOURCE: WEST WELL">   <!-- user defined name -->
             <Parameter name="regions" type="Array(string)" value="{WEST_WELL}"/>
             <Parameter name="spatial distribution method" type="string" value="permeability"/>
             <ParameterList name="sink">  
               <ParameterList name="function-constant">
                 <Parameter name="value" type="double" value="0.02"/>
               </ParameterList>
             </ParameterList>
           </ParameterList>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
    

Developer parameters
....................

The remaining parameters that can be used by a developer include

* `"enable internal tests`" [string] various internal tests will be executed during
  the run time. The default value is `"no`".
   
* `"internal tests tolerance`" [double] tolerance for internal tests such as the 
  divergence-free condition. The default value is 1e-6.

* `"runtime diagnostics: solute names`" [Array(string)] defines solutes that will be 
  tracked closely each time step if verbosity `"high`". Default value is the first 
  solute in the global list of `"aqueous names`" and the first gas in the global list 
  of `"gaseous names`".

* `"runtime diagnostics: regions`" [Array(string)] defines a boundary region for 
  tracking solutes. Default value is a seepage face boundary, see Flow PK.


Explanation of verbose output
.............................

When verbosity is set to *high*, this PK reports information about 
current status of the simulation.
Here after keyword *global* refers to the whole simulation including
all time periods, keyword *local* refers to the current time period.
The incomplete list is

 * [global] cycle number, time before step, and time step dt (in years)
 * [local] cell id and position with the smallest time step
 * [local] convergence of a linear solver for dispersion, PCG here
 * [local] number of subcycles, stable time step, and global time step (in seconds)
 * [local] species's name, concentration extrema, total amount of it in the 
   reservoir, and amount escaped through the outflow boundary
 * [global] current simulation time (in years)

.. code-block:: xml

  CycleDriver      |   Cycle 10: time(y) = 0.803511, dt(y) = 0.089279
  TransportPK      |    cell 0 has smallest dt, (-270, -270)
  TransportPK      |    dispersion solver (pcg) ||r||=8.33085e-39 itrs=2
  TransportPK      |    1 sub-cycles, dt_stable=2.81743e+06 [sec]  dt_MPC=2.81743e+06 [sec]
  TransportPK      |    Tc-99: min/max=7.111e-21 0.001461 [mol/m^3], total/out=2.2957 1.4211e-14 [mol]
  CycleDriver      |   New time(y) = 0.89279


Chemistry PK
------------

The chemistry header includes three parameters:

* `"chemistry model`" [string] defines chemical model. The available options are `"Alquimia`"
  and `"Amanzi`" (default).

.. code-block:: xml

  <ParameterList name="Chemistry">
    <Parameter name="chemistry model" type="string" value="Amanzi"/>
  </ParameterList>


Geochemical engines
...................

Here we specify either the default or the third-party geochemical engine. 


Alquimia
````````

The Alquimia chemistry process kernel only requires the *Engine* and *Engine Input File*
entries, but will also accept and respect the value given for *max time step (s)*. 
Most details are provided in the trimmed PFloTran file *1d-tritium-trim.in*.

* `"Minerals`" [Array(string)] is the list of mineral names.

* `"max time step (s)`" [double] is the maximum time step that chemistry will allow the MPC to take.

* `"min time step (s)`" [double] is the minimum time step that chemistry will allow the MPC to take.

* `"initial time step (s)`" [double] is the initial time step that chemistry will ask the MPC to take.

* `"time step control method`" [string] specifies time step control method for chemistry subcycling. 
  Choose either "fixed" (default) or "simple".  For option "fixed", time step is fixed.
  For option "simple", the time step is adjusted in response to stiffness of system of equations 
  based on a simple scheme. This option require the following parameters: `"time step cut threshold`",
  `"time step cut factor`", `"time step increase threshold`", and `"time step increase factor`".

* `"time step cut threshold`" [int] is the number of Newton iterations that if exceeded
  will trigger a time step cut. Default is 8.

* `"time step cut factor`" [double] is the factor by which the time step is cut. Default is 2.0

* `"time step increase threshold`" [int] is the number of consecutive successful time steps that
  will trigger a time step increase. Default is 4.

* `"time step increase factor`" [double] is the factor by which the time step is increased. Default is 1.2

.. code-block:: xml

  <ParameterList>  <!-- parent list -->
    <ParameterList name="Chemistry">
      <Parameter name="Engine" type="string" value="PFloTran"/>
      <Parameter name="Engine Input File" type="string" value="1d-tritium-trim.in"/>
      <Parameter name="Verbosity" type="Array(string)" value="{verbose}"/>
      <Parameter name="Minerals" type="Array(string)" value="{quartz, kaolinite, goethite, opal}"/>
      <Parameter name="min time step (s)" type="double" value="1.5778463e-07"/>
      <Parameter name="max time step (s)" type="double" value="1.5778463e+07"/>
      <Parameter name="initial time step (s)" type="double" value="1.0e-02"/>
      <Parameter name="time step control method" type="string" value="simple"/>
      <Parameter name="time step cut threshold" type="int" value="8"/>
      <Parameter name="time step cut factor" type="double" value="2.0"/>
      <Parameter name="time step increase threshold" type="int" value="4"/>
      <Parameter name="time step increase factor" type="double" value="1.2"/>
    </ParameterList>
  </ParameterList>


Amanzi
``````

The Amanzi chemistry process kernel uses the following parameters.

* `"activity model`" [string] is the type of model used for activity corrections. 
  Valid options are `"unit`" and `"debye-huckel`".

* `"tolerance`" [double] defines tolerance in Newton solves inside the chemistry library.

* `"maximum Newton iterations`" [int] is the maximum number of iteration the chemistry 
  library can take.

* `"auxiliary data`" [Array(string)] defines additional chemistry related data that the user 
  can request be saved to vis files. Currently `"pH`" is the only variable supported.

.. code-block:: xml

  <ParameterList>  <!-- parent list -->
    <ParameterList name="Chemistry">
      <ParameterList name="Thermodynamic Database">
        <Parameter name="Format" type="string" value="simple"/>
        <Parameter name="File" type="string" value="tritium.bgd"/>
      </ParameterList>
      <Parameter name="Verbosity" type="Array(string)" value="{verbose}"/>
      <Parameter name="activity model" type="string" value="unit"/>
      <Parameter name="tolerance" type="double" value="1.5e-12"/>
      <Parameter name="maximum Newton iterations" type="int" value="25"/>
      <Parameter name="max time step (s)" type="double" value="1.5e+07"/>
      <Parameter name="min time step (s)" type="double" value="1.0e-6"/>
      <Parameter name="number of component concentrations" type="int" value="1"/>
      <Parameter name="auxiliary data" type="Array(string)" value="{pH}"/>
    </ParameterList>
  </ParameterList>


Initial conditions
..................

This sublist completes initialization of state variable, see list `"State`" for 
more detail. This section is only required for the native chemistry kernel, the
Alquimia chemistry kernel reads initial conditions from the `"State`" list.
The following cell-based fields can be initialized here:

* `"mineral_volume_fractions`" (Alquimia only)
* `"mineral_specific_surface_area`" (Alqumia only)
* `"ion_exchange_sites`"
* `"ion_exchange_ref_cation_conc`"
* `"isotherm_kd`"
* `"isotherm_langmuir_b`"
* `"surface_complexation`"
* `"free_ion_species`"

.. code-block:: xml

  <ParameterList name="Chemistry">  <!-- parent list -->
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
  </ParameterList>


Format of chemistry database (.bgd) file
........................................

A section header starts with token `"<`". 
A comment line starts with token `"#`". 
Data fields are separated by semicolumns.


Primary Species
```````````````

Each line in this section has four data fields: 
name of a primary component, ion size parameter, charge, and atomic mass [u].

.. code-block:: txt

   <Primary Species
   Al+++    ;   9.0 ;   3.0 ;  26.9815
   Ca++     ;   6.0 ;   2.0 ;  40.078
   Cl-      ;   3.0 ;  -1.0 ;  35.4527
   CO2(aq)  ;   3.0 ;   0.0 ;  44.01
   Cs137    ;   2.5 ;   1.0 ; 132.9054
   F-       ;   3.5 ;  -1.0 ;  18.9984
   Fe++     ;   6.0 ;   2.0 ;  55.847
   H+       ;   9.0 ;   1.0 ;   1.0079
   HCO3-    ;   4.0 ;  -1.0 ;  61.0171
   HPO4--   ;   4.0 ;  -2.0 ;  95.9793
   K+       ;   3.0 ;   1.0 ;  39.0983
   Mg++     ;   8.0 ;   2.0 ;  24.30
   Na+      ;   4.0 ;   1.0 ;  22.9898
   N2(aq)   ;   3.0 ;   0.0 ;  28.0135
   NO3-     ;   3.0 ;  -1.0 ;  62.0049
   O2(aq)   ;   3.0 ;   0.0 ;  31.9988
   Pb_210   ;   1.0 ;   0.0 ; 210.00
   Pu_238   ;   1.0 ;   0.0 ; 238.00
   Ra_226   ;   1.0 ;   0.0 ; 226.00
   SiO2(aq) ;   3.0 ;   0.0 ;  60.0843
   SO4--    ;   4.0 ;  -2.0 ;  96.0636
   Sr90     ;   5.0 ;   2.0 ;  87.6200
   Tc_99    ;   1.0 ;   0.0 ;  99.00
   Th_230   ;   1.0 ;   0.0 ; 230.00
   Tritium  ;   9.0 ;   0.0 ;   1.01
   U_234    ;   1.0 ;   0.0 ; 234.00
   UO2++    ;   4.5 ;   2.0 ; 270.028
   Zn++     ;   6.0 ;   2.0 ;  65.39


Isotherms
`````````

Each line in this section has three fields: primary species name, 
adsorption isotherm model, and parameters. The adsorption model is one of: linear, langmuir, or freundlich.
The parameters is a space delimited list of numbers. The number of  parameters and 
their meaning depends on the model; although the first one is always *kd*.

.. code-block:: txt

   <Isotherm
   B      ; langmuir   ;      30.0 0.1
   C      ; freundlich ;       1.5 1.25
   Pu_238 ; linear     ;  461168.4
   U_234  ; linear     ;  329406.0
   Th_230 ; linear     ; 1482327.0
   Ra_226 ; linear     ;   41175.75
   Pb_210 ; linear     ; 3294060.0
   Tc_99  ; linear     ;     988.218


General Kinetics
````````````````

Each line in this section has five data fields.
The first field is the reaction string that has format 
"30 A(aq) + 2 B(aq) <-> C(aq) + .3 D(aq) +- 4 E(aq)"
where number (stoichiometires) is followed by species name. 
The second and fourth fields contain information about reactanct and products.
The fourth and fifth columns contain rate constants.

.. code-block:: txt

   <Primary Species
   Tritium  ;   9.00 ;   0.00 ;   1.01

   <General Kinetics
   1.00 Tritium <->  ;   1.00 Tritium ;  1.78577E-09 ; ; 


Aqueous Equilibrium Complexes
`````````````````````````````

Each line in this section has five 
fields for secondary species: name = coeff reactant, log Keq, size parameter, charge, and 
gram molecular weight.

.. code-block:: txt

   <Aqueous Equilibrium Complexes
   AlHPO4+    =  1.0 Al+++ 1.0 HPO4-- ;  -7.4    ;  4.0 ;   1.0 ; 122.961
   CaCl+      =  1.0 Ca++  1.0 Cl-    ;   0.6956 ;  4.0 ;   1.0 ;  75.5307
   CaCl2(aq)  =  1.0 Ca++  2.0 Cl-    ;   0.6436 ;  3.0 ;   0.0 ; 110.9834
   CaCO3(aq)  =  1.0 Ca+2  1.0 CO3-2  ;  -3.151  ;  3.5 ;   0.0 ;  17.0073
   CaHCO3+    =  1.0 Ca++  1.0 HCO3-  ;  -1.0467 ;  4.0 ;   1.0 ; 101.0951
   CaHPO4(aq) =  1.0 Ca++  1.0 HPO4-- ;  -2.7400 ;  3.0 ;   0.0 ; 136.0573
   CO3--      = -1.0 H+    1.0 HCO3-  ;  10.3288 ;  4.5 ;  -2.0 ;  60.0092
   FeCl+      =  1.0 Fe++  1.0 Cl-    ;   0.1605 ;  4.0 ;   1.0 ;  91.2997 
   FeCl2(aq)  =  1.0 Fe++  2.0 Cl-    ;   2.4541 ;  3.0 ;   0.0 ; 126.752 
   FeCl4--    =  1.0 Fe++  4.0 Cl-    ;  -1.9    ;  4.0 ;  -2.0 ; 197.658 
   FeHCO3+    =  1.0 Fe++  1.0 HCO3-  ;  -2.72   ;  4.0 ;   1.0 ; 116.864 
   FeHPO4(aq) =  1.0 Fe++  1.0 HPO4-- ;  -3.6    ;  3.0 ;   0.0 ; 151.826 
   FeF+       =  1.0 Fe++  1.0 F-     ;  -1.36   ;  4.0 ;   1.0 ;  74.8454 
   FeSO4(aq)  =  1.0 Fe++  1.0 SO4--  ;  -2.2    ;  3.0 ;   0.0 ; 151.911 
   H2PO4-     =  1.0 H+    1.0 HPO4-- ;  -7.2054 ;  4.0 ;  -1.0 ;  96.9872
   H3PO4(aq)  =  2.0 H+    1.0 HPO4-- ;  -9.3751 ;  3.0 ;   0.0 ;  97.9952
   H2SO4(aq)  =  2.0 H+    1.0 SO4--  ;   1.0209 ;  3.0 ;   0.0 ;  98.0795 
   HCl(aq)    =  1.0 H+    1.0 Cl-    ;  -0.67   ;  3.0 ;   0.0 ;  36.4606 
   HNO3(aq)   =  1.0 H+    1.0 NO3-   ;   1.3025 ;  3.0 ;   0.0 ;  63.0129 
   HSO4-      =  1.0 H+    1.0 SO4--  ;  -1.9791 ;  4.0 ;  -1.0 ;  97.0715 
   KCl(aq)    =  1.0 K+    1.0 Cl-    ;   1.4946 ;  3.0 ;   0.0 ;  74.551
   KHPO4-     =  1.0 K+    1.0 HPO4-- ;  -0.78   ;  4.0 ;  -1.0 ; 135.078
   KSO4-      =  1.0 K+    1.0 SO4--  ;  -0.8796 ;  4.0 ;  -1.0 ; 135.162
   NaHCO3(aq) =  1.0 Na+   1.0 HCO3-  ;  -0.1541 ;  3.0 ;   0.0 ;  84.0069 
   NaCl(aq)   =  1.0 Na+   1.0 Cl-    ;   0.777  ;  3.0 ;   0.0 ;  58.4425 
   NaF(aq)    =  1.0 Na+   1.0 F-     ;   0.9976 ;  3.0 ;   0.0 ;  41.9882  
   NaHPO4-    =  1.0 Na+   1.0 HPO4-- ;  -0.92   ;  4.0 ;  -1.0 ; 118.969
   NaNO3(aq)  =  1.0 Na+   1.0 NO3-   ;   1.044  ;  3.0 ;   0.0 ;  84.9947
   NaSO4-     =  1.0 Na+   1.0 SO4--  ;  -0.82   ;  4.0 ;  -1.0 ; 119.053
   MgCO3(aq)  =  1.0 Mg+2  1.0 CO3-2  ;  -2.928  ;  3.5 ;   0.0 ;  17.0073
   OH-        =  1.0 H2O  -1.0 H+     ;  13.9951 ;  3.5 ;  -1.0 ;  17.00730
   P2O7----   = -1.0 H2O   2.0 HPO4-- ;   3.7463 ;  4.0 ;  -4.0 ; 173.9433
   PO4---     = -1.0 H+    1.0 HPO4-- ;  12.3218 ;  4.0 ;  -3.0 ;  94.9714
   UO2Cl+     =  1.0 Cl-   1.0 UO2++  ;  -0.1572 ;  4.0 ;   1.0 ; 305.48 
   UO2Cl2(aq) =  2.0 Cl-   1.0 UO2++  ;   1.1253 ;  3.0 ;   0.0 ; 340.933 
   UO2F+      =  1.0 F-    1.0 UO2++  ;  -5.0502 ;  4.0 ;   1.0 ; 289.026 
   UO2F2(aq)  =  2.0 F-    1.0 UO2++  ;  -8.5403 ;  3.0 ;   0.0 ; 308.024 
   UO2F3-     =  3.0 F-    1.0 UO2++  ; -10.7806 ;  4.0 ;  -1.0 ; 327.023 
   UO2F4--    =  4.0 F-    1.0 UO2++  ; -11.5407 ;  4.0 ;  -2.0 ; 346.021 
   UO2HPO4(aq)= 1.0 HPO4-- 1.0 UO2++  ;  -8.4398 ;  3.0 ;   0.0 ; 366.007 
   UO2NO3+    =  1.0 NO3-  1.0 UO2++  ;  -0.2805 ;  4.0 ;   1.0 ; 332.033 
   UO2SO4(aq) =  1.0 SO4-- 1.0 UO2++  ;  -3.0703 ;  3.0 ;   0.0 ; 366.091 
   UO2(SO4)2-- = 2.0 SO4-- 1.0 UO2++  ;  -3.9806 ;  4.0 ;  -2.0 ; 462.155 

   Al2(OH)2++++	 = -2.0 H+    2.0 Al+++   2.0 H2O     ;   7.6902 ;  5.5	;  4.0 ;  87.9778
   Al3(OH)4(5+)	 = -4.0 H+    3.0 Al+++   4.0 H2O     ;  13.8803 ;  6.0 ;  5.0 ; 148.9740
   Al(OH)2+      =  2.0 H2O   1.0 Al+++  -2.0 H+      ;  10.5945 ;  4.0 ;  1.0 ;  60.9962 
   Al(OH)3(aq)   =  3.0 H2O   1.0 Al+++  -3.0 H+      ;  16.1577 ;  3.0 ;  0.0 ;  78.0034 
   Al(OH)4-      =  4.0 H2O   1.0 Al+++  -4.0 H+      ;  22.8833 ;  4.0 ; -1.0 ;  95.0107 
   AlH2PO4++     =  1.0 Al+++ 1.0 H+      1.0 HPO4--  ;  -3.1    ;  4.5 ;  2.0 ; 123.969
   AlO2-         =  2.0 H2O   1.0 Al+++  -4.0 H+      ;  22.8833 ;  4.0 ; -1.0 ;  58.9803
   AlOH++        =  1.0 H2O   1.0 Al+++  -1.0 H+      ;   4.9571 ;  4.5 ;  2.0 ;  43.9889 
   CaCO3(aq)     = -1.0 H+    1.0 Ca++    1.0 HCO3-   ;   7.0017 ;  3.0 ;  0.0 ; 100.0872
   CaH2PO4+      =  1.0 Ca++  1.0 H+      1.0 HPO4--  ;  -1.4000 ;  4.0 ;  1.0 ; 137.0652
   CaP2O7--      = -1.0 H2O   1.0 Ca++    2.0 HPO4--  ;  -3.0537 ;  4.0 ; -2.0 ; 214.0213
   CaPO4-        = -1.0 H+    1.0 Ca++    1.0 HPO4--  ;   5.8618 ;  4.0 ; -1.0 ; 135.0494
   CaOH+         = -1.0 H+    1.0 Ca++    1.0 H2O     ;  12.8500 ;  4.0 ;  1.0 ;  57.0853
   CO2(aq)       = -1.0 H2O   1.0 H+      1.0 HCO3-   ;  -6.3447 ;  3.0 ;  0.0 ;  44.0098
   H2P2O7--      = -1.0 H2O   2.0 H+      2.0 HPO4--  ; -12.0709 ;  4.0 ; -2.0 ; 175.9592
   H2S(aq)       =  2.0 H+    1.0 SO4--  -2.0 O2(aq)  ; 131.329  ;  3.0 ;  0.0 ;  34.0819 
   H3P2O7-       = -1.0 H2O   2.0 HPO4--  3.0 H+      ; -14.4165 ;  4.0 ; -1.0 ; 176.9671
   H4P2O7(aq)    = -1.0 H2O   2.0 HPO4--  4.0 H+      ; -15.9263 ;  3.0 ;  0.0 ; 177.9751
   HAlO2(aq)     =  2.0 H2O   1.0 Al+++  -3.0 H+      ;  16.4329 ;  3.0 ;  0.0 ;  59.9883
   HCO3-         =  1.0 H2O  -1.0 H+      1.0 CO2(aq) ;  6.34470 ;  4.0 ; -1.0 ;  61.01710
   HO2-          =  1.0 H2O  -1.0 H+      0.5 O2(aq)  ;  28.302  ;  4.0 ; -1.0 ;  33.0067 
   HP2O7---      = -1.0 H2O   1.0 H+      2.0 HPO4--  ;  -5.4498 ;  4.0 ; -3.0 ; 174.9513
   HS-           =  1.0 H+    1.0 SO4--  -2.0 O2(aq)  ; 138.317  ;  3.5 ; -1.0 ;  33.0739 
   Fe2(OH)2++++  =  1.0 H2O   2.0 Fe++    0.5 O2(aq)  ; -14.0299 ;  5.5 ;  4.0 ; 145.709 
   FeCO3(aq)     =  1.0 Fe++ -1.0 H+      1.0 HCO3-   ;   5.5988 ;  3.0 ;  0.0 ; 115.856 
   FeH2PO4+      =  1.0 Fe++  1.0 H+      1.0 HPO4--  ;  -2.7    ;  4.0 ;  1.0 ; 152.834 
   Fe(OH)2(aq)   =  2.0 H2O   1.0 Fe++   -2.0 H+      ;  20.6    ;  3.0 ;  0.0 ;  89.8617 
   Fe(OH)3-      =  3.0 H2O   1.0 Fe++   -3.0 H+      ;  31.0    ;  4.0 ; -1.0 ; 106.869 
   Fe(OH)4--     =  4.0 H2O   1.0 Fe++   -4.0 H+      ;  46.0    ;  4.0 ; -2.0 ; 123.876 
   FeOH+         =  1.0 H2O   1.0 Fe++   -1.0 H+      ;   9.5    ;  4.0 ;  1.0 ;  72.8543 
   FeOH++        =  0.5 H2O   1.0 Fe++    0.25 O2(aq) ;  -6.3    ;  4.5 ;  2.0 ;  72.8543 
   FePO4-        =  1.0 Fe++ -1.0 H+      1.0 HPO4--  ;   4.3918 ;  4.0 ; -1.0 ; 150.818 
   KHSO4(aq)     =  1.0 K+    1.0 H+      1.0 SO4--   ;  -0.8136 ;  3.0 ;  0.0 ; 136.17
   KOH(aq)       =  1.0 H2O   1.0 K+     -1.0 H+      ;  14.46   ;  3.0 ;  0.0 ;  56.1056
   KP2O7---      = -1.0 H2O   1.0 K+      2.0 HPO4--  ;   1.4286 ;  4.0 ; -3.0 ; 213.042
   MgOH+         =  1.0 H2O  -1.0 H+      1.0 Mg++    ; 11.78510 ;  4.0 ;  1.0 ;  41.3123
   NaCO3-        = -1.0 H+    1.0 HCO3-   1.0 Na+     ;   9.8144 ;  4.0 ; -1.0 ;  82.9990
   NaOH(aq)      =  1.0 H2O   1.0 Na+    -1.0 H+      ;  14.7948 ;  3.0 ;  0.0 ;  39.9971
   NH3(aq)       =  1.5 H2O   0.5 N2(aq) -0.75 O2(aq) ;  58.2305 ;  3.0 ;  0.0 ;  17.0306 
   UO2CO3(aq)    = -1.0 H+    1.0 HCO3-   1.0 UO2++   ;   0.6634 ;  3.0 ;  0.0 ; 330.037 
   UO2(CO3)2--   = -2.0 H+    2.0 HCO3-   1.0 UO2++   ;   3.7467 ;  4.0 ; -2.0 ; 390.046 
   UO2(CO3)3---- = -3.0 H+    3.0 HCO3-   1.0 UO2++   ;   9.4302 ;  4.0 ; -4.0 ; 450.055 
   UO2H2PO4+     =  1.0 H+    1.0 HPO4--  1.0 UO2++   ; -11.6719 ;  4.0 ;  1.0 ; 367.015 
   UO2H3PO4++    =  2.0 H+    1.0 HPO4--  1.0 UO2++   ; -11.3119 ;  4.5 ;  2.0 ; 368.023 
   UO2OH+        =  1.0 H2O  -1.0 H+      1.0 UO2++   ;   5.2073 ;  4.0 ;  1.0 ; 287.035 
   UO2PO4-       = -1.0 H+    1.0 HPO4--  1.0 UO2++   ;  -2.0798 ;  4.0 ; -1.0 ; 364.999 
   UO2(OH)2(aq)  =  2.0 H2O  -2.0 H+      1.0 UO2++   ;  10.3146 ;  3.0 ;  0.0 ; 304.042 
   UO2(OH)3-     =  3.0 H2O  -3.0 H+      1.0 UO2++   ;  19.2218 ;  4.0 ; -1.0 ; 321.05 
   UO2(OH)4--    =  4.0 H2O  -4.0 H+      1.0 UO2++   ;  33.0291 ;  4.0 ; -2.0 ; 338.057 
   (UO2)2OH+++   =  1.0 H2O  -1.0 H+      2.0 UO2++   ;   2.7072 ;  5.0 ;  3.0 ; 557.063 
   (UO2)2(OH)2++ =  2.0 H2O  -2.0 H+      2.0 UO2++   ;   5.6346 ;  4.5 ;  2.0 ; 574.07 
   (UO2)3(OH)4++ =  4.0 H2O  -4.0 H+      3.0 UO2++   ;  11.929  ;  4.5 ;  2.0 ; 878.112 
   (UO2)3(OH)5+  =  5.0 H2O  -5.0 H+      3.0 UO2++   ;  15.5862 ;  4.0 ;  1.0 ; 895.12 
   (UO2)3(OH)7-  =  7.0 H2O  -7.0 H+      3.0 UO2++   ;  31.0508 ;  4.0 ; -1.0 ; 929.135 
   (UO2)4(OH)7+  =  7.0 H2O  -7.0 H+      4.0 UO2++   ;  21.9508 ;  4.0 ;  1.0 ; 1199.16 
   UO2(H2PO4)(H3PO4)+ = 3.0 H+ 2.0 HPO4-- 1.0 UO2++   ; -22.7537 ;  4.0 ;  1.0 ; 465.01 
   UO2(H2PO4)2(aq) =    2.0 H+ 2.0 HPO4-- 1.0 UO2++   ; -21.7437 ;  3.0 ;  0.0 ; 464.002 
   Zn(OH)2(aq)   =  2.0 H2O  -2.0 H+      1.0 Zn++    ;  17.3282 ;  3.0 ;  0.0 ;  99.4047
   Zn(OH)3-      =  3.0 H2O  -3.0 H+      1.0 Zn++    ;  28.8369 ;  4.0 ; -1.0 ; 116.41200
   Zn(OH)4--     =  4.0 H2O  -4.0 H+      1.0 Zn++    ;  41.6052 ;  4.0 ; -2.0 ; 133.41940
   ZnOH+         =  1.0 H2O  -1.0 H+      1.0 Zn++    ;   8.9600 ;  4.0 ;  1.0 ;  82.39730

   Ca2UO2(CO3)3(aq) =  2.0 Ca++ -3.0 H+     3.0 HCO3-   1.0 UO2++    ;   0.2864 ; 4.0 ;  0.0 ; 530.215 
   CaUO2(CO3)3--    =  1.0 Ca++ -3.0 H+     3.0 HCO3-   1.0 UO2++    ;   3.8064 ; 4.0 ; -2.0 ; 530.215 
   CH4(aq)          =  1.0 H2O   1.0 H+     1.0 HCO3-  -2.0 O2(aq)   ; 144.141  ; 3.0 ;  0.0 ;   0.0
   NaAlO2(aq)       =  2.0 H2O   1.0 Na+    1.0 Al+++  -4.0 H+       ;  23.6266 ; 3.0 ;  0.0 ;  81.9701
   NaHP2O7--        = -1.0 H2O   1.0 Na+    1.0 H+      2.0 HPO4--   ;  -6.8498 ; 4.0 ; -2.0 ; 197.941
   NaHSiO3(aq)      =  1.0 H2O   1.0 Na+   -1.0 H+      1.0 SiO2(aq) ;   8.304  ; 3.0 ;  0.0 ; 100.081
   Fe+++            = -0.5 H2O   1.0 Fe++   1.0 H+      0.25 O2(aq)  ;  -8.49   ; 9.0 ;  3.0 ;  55.847 
   Fe3(OH)4(5+)     =  2.5 H2O   3.0 Fe++  -1.0 H+      0.75 O2(aq)  ; -19.1699 ; 6.0 ;  5.0 ; 235.57 
   Fe(OH)2+         =  1.5 H2O   1.0 Fe++  -1.0 H+      0.25 O2(aq)  ;  -2.82   ; 4.0 ;  1.0 ;  89.8617 
   Fe(OH)3(aq)      =  2.5 H2O   1.0 Fe++  -2.0 H+      0.25 O2(aq)  ;   3.51   ; 3.0 ;  0.0 ; 106.869 
   Fe(OH)4-         =  3.5 H2O   1.0 Fe++  -3.0 H+      0.25 O2(aq)  ;  13.11   ; 4.0 ; -1.0 ; 123.876 
   FeCO3+           = -0.5 H2O   1.0 Fe++   1.0 HCO3-   0.25 O2(aq)  ;  -7.8812 ; 4.0 ;  1.0 ; 115.856 
   MgHCO3+          =  1.0 H2O  -1.0 H+     1.0 CO2(aq) 1.0 Mg++     ;   5.309  ; 4.0 ;  1.0 ;  85.3221
   N3-              =  0.5 H2O  -1.0 H+     1.5 N2(aq) -0.25 O2(aq)  ;  77.7234 ; 4.0 ; -1.0 ;  42.0202 
   NH4+             =  1.5 H2O   1.0 H+     0.5 N2(aq) -0.75 O2(aq)  ;  48.9895 ; 2.5 ;  1.0 ;  18.0385 
   U+++             = -0.5 H2O   1.0 H+     1.0 UO2++  -0.75 O2(aq)  ;  64.8028 ; 5.0 ;  3.0 ; 238.029 
   U++++            = -1.0 H2O   2.0 H+     1.0 UO2++  -0.5 O2(aq)   ;  33.949  ; 5.5 ;  4.0 ; 238.029 
   UO2+             =  0.5 H2O  -1.0 H+     1.0 UO2++  -0.25 O2(aq)  ;  20.0169 ; 4.0 ;  1.0 ; 270.028 
   UO2OSi(OH)3+     =  2.0 H2O  -1.0 H+     1.0 SiO2(aq) 1.0 UO2++   ;   2.4810 ; 9.0 ;  1.0 ; 365.135
   (UO2)2CO3(OH)3-  =  3.0 H2O  -4.0 H+     1.0 HCO3-   2.0 UO2++    ;  11.2229 ; 4.0 ; -1.0 ; 651.087 

   Fe(SO4)2- = -0.5 H2O  1.0 Fe++  1.0 H+      2.0 SO4--   0.25 O2(aq) ; -11.7037 ; 4.0 ; -1.0 ; 247.974 
   FeCl++    = -0.5 H2O  1.0 Fe++  1.0 H+      1.0 Cl-     0.25 O2(aq) ;  -7.6792 ; 4.5 ;  2.0 ;  91.2997 
   FeCl2+    = -0.5 H2O  1.0 Fe++  1.0 H+      2.0 Cl-     0.25 O2(aq) ; -10.62   ; 4.0 ;  1.0 ; 126.752 
   FeCl4-    = -0.5 H2O  1.0 Fe++  1.0 H+      4.0 Cl-     0.25 O2(aq) ;  -7.7    ; 4.0 ; -1.0 ; 197.658 
   FeF++     = -0.5 H2O  1.0 Fe++  1.0 H+      1.0 F-      0.25 O2(aq) ; -12.6265 ; 4.5 ;  2.0 ;  74.8454 
   FeF2+     = -0.5 H2O  1.0 Fe++  1.0 H+      2.0 F-      0.25 O2(aq) ; -16.8398 ; 4.0 ;  1.0 ;  93.8438 
   FeH2PO4++ = -0.5 H2O  1.0 Fe++  2.0 H+      1.0 HPO4--  0.25 O2(aq) ; -12.66   ; 4.5 ;  2.0 ; 152.834 
   FeHPO4+   = -0.5 H2O  1.0 Fe++  1.0 H+      1.0 HPO4--  0.25 O2(aq) ; -18.67   ; 4.0 ;  1.0 ; 151.826 
   FeNO3++   = -0.5 H2O  1.0 Fe++  1.0 H+      1.0 NO3-    0.25 O2(aq) ;  -9.49   ; 4.5 ;  2.0 ; 117.852 
   FeSO4+    = -0.5 H2O  1.0 Fe++  1.0 H+      1.0 SO4--   0.25 O2(aq) ; -10.4176 ; 4.0 ;  1.0 ; 151.911 
   MgUO2(CO3)3-- = 3.0 H2O -6.0 H+ 3.0 CO2(aq) 1.0 Mg++    1.00 UO2++  ;  23.9105 ; 3.0 ; -2.0 ; 500.0
   NH4SO4-   =  1.5 H2O  1.0 H+    0.5 N2(aq)  1.0 SO4--  -0.75 O2(aq) ;  57.2905 ; 4.0 ; -1.0 ; 114.102


Minerals
````````

Each line in this section has five fields for secondary species:
Name = coeff reactant, log Keq, gram molecular weight [g/mole], molar volume [cm^3/mole],
and specific surface area [cm^2 mineral / cm^3 bulk].

.. code-block:: txt

   <Minerals
   Chalcedony = 1.0 SiO2(aq) ; -3.7281 ; 60.0843 ; 22.688 ; 1.0
   Opal       = 1.0 SiO2(aq) ; -3.005  ; 60.084  ; 29.0   ; 1.0
   Quartz     = 1.0 SiO2(aq) ; -3.9993 ; 60.0843 ; 22.688 ; 1.0

   Halite = 1.0 Na+  1.0 Cl- ; 1.58550 ; 58.4425 ; 27.0150 ; 1.0

   Gypsum    =  2.0 H2O   1.0 SO4-2   1.0 Ca++  ; -4.581  ; 172.1722 ; 74.21216 ; 1.0
   Calcite   = -1.0 H+    1.0 HCO3-   1.0 Ca++  ;  1.8487 ; 100.087  ; 36.934   ; 1.0
   Gibbsite  =  3.0 H2O   1.0 Al+++  -3.0 H+    ;  7.756  ;  78.0036 ; 31.956   ; 1.0
   Schoepite =  3.0 H2O  -2.0 H+      1.0 UO2++ ;  4.8333 ; 322.058  ; 66.08    ; 1.0

   Basaluminite = 15.0 H2O -10.0 H+     4.0 Al+++    1.0 SO4--    ;  22.2511 ;  464.140 ; 218.934 ; 1.0
   Ferrihydrite = 2.5 H2O   1.0 Fe++   -2.0 H+       0.25 O2(aq)  ;  -3.594  ;  106.869 ;  23.99  ; 1.0
   Jurbanite    = 6.0 H2O   1.0 Al+++  -1.0 H+       1.0 SO4--    ;  -3.23   ;  230.064 ; 126.0   ; 1.0
   Kaolinite    = 5.0 H2O  -6.0 H+      2.0 Al+++    2.0 SiO2(aq) ;   7.570  ;  258.160 ;  99.520 ; 1.0
   Soddyite     = 4.0 H2O  -4.0 H+      1.0 SiO2(aq) 2.0 UO2++    ;   0.392  ;  668.169 ; 131.27  ; 1.0
   (UO2)3(PO4)2.4H2O = 4.0 H2O -2.0 H+  2.0 HPO4--   3.0 UO2++    ; -27.0349 ; 1072.09  ; 500.0   ; 1.0

   K-Feldspar = 2.0 H2O  1.0 K+   1.0 Al+++  -4.0 H+  3.0 SiO2(aq) ;  -0.2753 ; 278.332 ; 108.87   ; 1.0
   Polyhalite = 2.0 H2O  1.0 Mg++ 2.0 Ca++    2.0 K+  4.0 SO4-2    ; -13.7440 ; 218.1   ; 100.9722 ; 1.0


Mineral Kinetics
````````````````

Each line in this section has four fields.
The first field contains mineral name that is assumed to have the same stoichiometry 
as the mineral definition.
The second field is the rate name.

.. code-block:: txt

   <Mineral Kinetics
   Basaluminite ; TST ; log10_rate_constant    -8.0 moles/cm^2/sec
   Calcite      ; TST ; log10_rate_constant   -10.0 moles/m^2/sec
   Chalcedony   ; TST ; log10_rate_constant   -14.0 moles/m^2/sec
   Ferrihydrite ; TST ; log10_rate_constant   -14.0 moles/m^2/sec
   Gibbsite     ; TST ; log10_rate_constant   -14.0 moles/m^2/sec
   Halite       ; TST ; log10_rate_constant   -40.0 moles/cm^2/sec
   Jurbanite    ; TST ; log10_rate_constant   -14.0 moles/m^2/sec
   K-Feldspar   ; TST ; log10_rate_constant   -16.699 moles/m^2/sec
   Kaolinite    ; TST ; log10_rate_constant   -16.699 moles/m^2/sec
   Opal         ; TST ; log10_rate_constant   -12.135 moles/cm^2/sec
   Quartz       ; TST ; log10_rate_constant   -18.0 moles/m^2/sec
   Schoepite    ; TST ; log10_rate_constant   -10.0 moles/m^2/sec
   Soddyite     ; TST ; log10_rate_constant   -10.0 moles/m^2/sec
   (UO2)3(PO4)2.4H2O ; TST ; log10_rate_constant  -10.0 moles/m^2/sec


Ion Exchange Sites
``````````````````

Each line in this section has three fields: 
exchanger name, exchanger change, and exchanger location. 
The location is the mineral where the exchanger is located, i.e. kaolinite.

.. code-block:: txt

  <Ion Exchange Sites
   X- ; -1.0 ; Halite


Ion Exchange Complexes
``````````````````````

Each line in this section has two fields.
The first field has format "name = coeffient and primary name followed by coefficient 
and exchanger name. the second field is Keq.
The following assumptions are made:

   - The coefficient of the ion exchange complex is one.
   - Each complexation reaction is written between a single
     primary species and a single exchange site.

.. code-block:: txt

   <Ion Exchange Complexes
   Al+++X = 1.0 Al+++ 3.0 X- ;  1.71133
   Ca++X  = 1.0 Ca++  2.0 X- ;  0.29531
   Ca0.5X = 0.5 Ca++  1.0 X- ; -0.99
   H+X    = 1.0 H+    1.0 X- ;  0.0251189
   Mg++X  = 1.0 Mg++  2.0 X- ;  0.1666
   Na+X   = 1.0 Na+   1.0 X- ;  1.0
   NaX    = 1.0 Na+   1.0 X- ;  0.0


Surface Complex Sites
`````````````````````

Each line in this section has two fields: species name and surface density.

.. code-block:: txt

   <Surface Complex Sites
   >AlOH   ; 6.3600E-03
   >FeOH   ; 6.3600E-03
   >FeOH_w ; 7.6355E+04
   >FeOH_s ; 1.9080E+03
   >SiOH   ; 6.3600E-03
   >davis_OH ; 1.56199E-01


Surface Complexes
`````````````````

Each line in this section has three fields
for secondary species. The first field has format "name = coefficient primary_name coefficient exchanger site".
The second field is Keq. The third field is charge.

.. code-block:: txt

   <Surface Complexes
   >FeOH2+_w  = 1.0 >FeOH_w   1.0 H+    ; -7.18 ;  1.0
   >FeOH2+_s  = 1.0 >FeOH_s   1.0 H+    ; -7.18 ;  1.0
   >FeO-_w    = 1.0 >FeOH_w  -1.0 H+    ;  8.82 ; -1.0
   >FeO-_s    = 1.0 >FeOH_s  -1.0 H+    ;  8.82 ; -1.0
   >FeOHUO2++ = 1.0 >FeOH     1.0 UO2++ ; -6.63 ;  2.0
   >SiO-      =-1.0 H+        1.0 >SiOH ;  0.0
   >SiOH2+    = 1.0 H+        1.0 >SiOH ;  0.0

   >AlOUO2+    = 1.0 >AlOH   -1.0 H+   1.0 UO2++ ; -3.13 ; 1.0
   >FeOHZn+_w  = 1.0 >FeOH_w -1.0 H+   1.0 Zn++  ;  2.32 ; 1.0
   >FeOHZn+_s  = 1.0 >FeOH_s -1.0 H+   1.0 Zn++  ; -0.66 ; 1.0
   >SiOUO3H3++ = 1.0 >SiOH    1.0 H2O  1.0 UO2++ ;  5.18 ; 2.0
   >UO2++      = 1.0 UO2++   -1.0 Ca++ 1.0 >Ca++ ; -5.12 ; 0.0
   (>davis_O)UO2+ = 1.0 >davis_OH -1.0 H+ 1.0 UO2++; -0.444 ; 1.0

   >SiOUO3H2+    = 1.0 >SiOH  1.0 H2O  -1.0 H+  1.0 UO2++ ;  5.18 ;  1.0
   >SiOUO3H      = 1.0 >SiOH  1.0 H2O  -2.0 H+  1.0 UO2++ ;  5.18 ;  0.0
   >SiOUO3-      = 1.0 >SiOH  1.0 H2O  -3.0 H+  1.0 UO2++ ; 12.35 ; -1.0
   >SiOUO2(OH)2- = 1.0 >SiOH  2.0 H2O  -3.0 H+  1.0 UO2++ ; 12.35 ; -1.0
   >FeOHUO3      = 1.0 >FeOH  1.0 H2O  -2.0 H+  1.0 UO2++ ;  3.05 ;  0.0


Radiactive Decay
````````````````

Each line in this section has two fields.
The first field has format "parent name --> stoichiometric coefficient and species name.
The second fields is half-life time with units.
The stoichiometric coefficient of the parent should always be one.
The units is one of the following: years, days, hours, minutes, or seconds.

.. code-block:: txt

   <Radioactive Decay
   Cs137  -->  1.0 Cs137  ; half_life 30.2 years
   Pb_210 -->             ; half_life 22.2 years
   Pu_238 -->  1.0 U_234  ; half_life 87.7 years
   Ra_226 -->  1.0 Pb_210 ; half_life 1.6e3 years
   Th_230 -->  1.0 Ra_226 ; half_life 7.54e4 years
   U_234  -->  1.0 Th_230 ; half_life 2.45e5 years
   Tc_99  -->             ; half_life 2.111e5 years
   Sr90   -->  1.0 Sr90   ; half_life 28.8 years


Energy PK
---------

The conceptual PDE model for the energy equation is 

.. math::
  \frac{\partial \varepsilon}{\partial t} 
  =
  \boldsymbol{\nabla} \cdot (\kappa \nabla T) -
  \boldsymbol{\nabla} \cdot (\eta_l H_l \boldsymbol{q}_l) + Q

where 
:math:`\varepsilon` is the internal energy,
:math:`\eta_l` is molar density of liquid,
:math:`Q` is source or sink term,
:math:`\boldsymbol{q}_l` is the Darcy velocity,
:math:`\kappa` is thermal conductivity,
and :math:`H_l` is molar enthalpy of liquid.
We define 

.. math::
   \varepsilon = \phi (\eta_l s_l U_l + \eta_g s_g U_g) + 
   (1 - \phi) \rho_r c_r T

where
:math:`s_l` is liquid saturation,
:math:`s_g` is gas saturation (water vapor),
:math:`\eta_l` is molar density of liquid,
:math:`\eta_g` is molar density of gas,
:math:`U_l` is molar internal energy of liquid,
:math:`U_g` is molar internal energy of gas (water vapor),
:math:`\phi` is porosity,
:math:`\rho_r` is rock density,
:math:`c_r` is specific heat of rock,
and :math:`T` is temperature.

Energy sublist includes exactly one sublist, either `"Single-phase problem`" or `"Two-phase problem`".
Structure of both sublists is quite similar. We make necessary comments on their differences.


Physical models and assumptions
...............................

This list is used to summarize physical models and assumptions, such as
coupling with other PKs.
This list is often generated on a fly by a high-level MPC PK.

* `"vapor diffusion`" [bool] is set up automatically by a high-level PK,
  e.g. by EnergyFlow PK. The default value is `"false`".

* `"water content model`" [string] changes the evaluator for water
  content. Available options are `"generic`" and `"constant density`" (default).

.. code-block:: xml

   <ParameterList name="Energy">  <!-- parent list -->
     <ParameterList name="physical models and assumptions">
       <Parameter name="vapor diffusion" type="bool" value="false"/>
       <Parameter name="water content model" type="string" value="constant density"/>
     </ParameterList>
   </ParameterList>


Internal energy
...............

Internal energy list has a few parameters that allows us to run this PK
in a variety of regimes, e.g. with or without gas phase.

* `"energy key`" [string] specifies name for the internal energy field.
  The default value is `"energy`".

* `"evaluator type`" [string] changes the evaluator for internal energy.
  Available options are `"generic`" and `"constant liquid density`" (default).

* `"vapor diffusion`" [bool] specifies presence of a gas phase.
  The default value is `"true`".

.. code-block:: xml

   <ParameterList name="Energy">  <!-- parent list -->
     <ParameterList name="energy evaluator">
       <Parameter name="energy key" type="string" value="energy"/>
       <Parameter name="evaluator type" type="string" value="constant liquid density"/>
       <Parameter name="vapor diffusion" type="bool" value="true"/>
       <ParameterList name="VerboseObject">
         <Parameter name="Verbosity Level" type="string" value="high"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>


Molar enthalpy
..............

.. code-block:: xml

   <ParameterList name="Energy">  <!-- parent list -->
     <ParameterList name="enthalpy evaluator">
       <Parameter name="enthalpy key" type="string" value="enthalpy_liquid"/>
       <Parameter name="internal energy key" type="string" value="internal_energy_liquid"/>

       <Parameter name="include work term" type="bool" value="true"/>
       <Parameter name="pressure key" type="string" value="pressure"/>
       <Parameter name="molar density key" type="string" value="molar_density_liquid"/>
     </ParameterList>
   </ParameterList>


Thermal conductivity
....................

Evaluator for thermal conductivity allows us to select a proper model. 
The variety of available models allows to run the energy PK by itself or in
coupling with flow PK. 
The structure of the thermal conductivity list resembles that of a field
evaluator list in state. 
The two-phase model accepts the following parameters.

* `"thermal conductivity parameters`" [list] defines a model and its parameters.

* `"thermal conductivity type`" [string] is the name of a conductivity model in the
  list of registered models. Available two-phase models are `"two-phase Peters-Lidard`",
  and `"two-phase wet/dry`". Available one-phase model is `"one-phase polynomial`".

* `"thermal conductivity of rock`" [double] defines constant conductivity of rock.

* `"thermal conductivity of gas`" [double] defines constant conductivity of gas.

* `"thermal conductivity of liquid`" [double] defines constant conductivity of fluid.
  Default value is 0.6065 [W/m/K].

* `"unsaturated alpha`" [double] is used to define the Kersten number to interpolate
  between saturated and dry conductivities.

* `"epsilon`" [double] is needed for the case of zero saturation. Default is `"1.0e-10`".

.. code-block:: xml

   <ParameterList name="Energy">  <!-- parent list -->
     <ParameterList name="thermal conductivity evaluator">
       <ParameterList name="thermal conductivity parameters">
         <Parameter name="thermal conductivity type" type="string" value="two-phase Peters-Lidard"/>
         <Parameter name="thermal conductivity of rock" type="double" value="0.2"/>
         <Parameter name="thermal conductivity of gas" type="double" value="0.02"/>
         <Parameter name="thermal conductivity of liquid" type="double" value="0.6065"/>

         <Parameter name="unsaturated alpha" type="double" value="1.0"/>
         <Parameter name="epsilon" type="double" value="1.e-10"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>

The single-phase model accepts some of the parameters defined above (see the example) 
and a few additional parameters.

* `"reference temperature`" [double] defines temperature at which reference conductivity
  of liquid is calculated. Default value is 298.15 [K].

* `"polynomial expansion`" [Array(double)] collect coefficients in the quadratic representation of the 
  thermal conductivity of liquid with respect to the dimensionless parameter T/Tref.

.. code-block:: xml

   <ParameterList name="Energy">  <!-- parent list -->
     <ParameterList name="thermal conductivity evaluator">
       <ParameterList name="thermal conductivity parameters">
         <Parameter name="thermal conductivity type" type="string" value="one-phase polynomial"/>
         <Parameter name="thermal conductivity of rock" type="double" value="0.2"/>
         <Parameter name="reference temperature" type="double" value="298.15"/>
         <Parameter name="polynomial expansion" type="Array(double)" value="{-1.48445, 4.12292, -1.63866}"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>


Operators
.........

This section contains sublist for diffusion and advection operators.
It also has one global parameters.

* `"operators`" [list] 
  
  * `"include enthalpy in preconditioner`" [bool] allows us to study impact (usually positive) 
    of including enthalpy term in the preconditioner. Default value is *true*.


Diffusion operator
``````````````````

Operators sublist describes the PDE structure of the flow, specifies a discretization
scheme, and selects assembling schemas for matrices and preconditioners.

* `"diffusion operator`" [list] defines parameters for generating and assembling diffusion matrix.

  * `"matrix`" [list] defines parameters for generating and assembling diffusion matrix. See section
    describing operators. 
    When `"Richards problem`" is selected, Flow PK sets up proper value for parameter `"upwind method`" of 
    this sublist.

  * `"preconditioner`" [list] defines parameters for generating and assembling diffusion 
    matrix that is used to create preconditioner. 
    This sublist is ignored inside sublist `"Darcy problem`".
    Since update of preconditioner can be lagged, we need two objects called `"matrix`" and `"preconditioner`".
    When `"Richards problem`" is selected, Flow PK sets up proper value for parameter `"upwind method`" of 
    this sublist.

.. code-block:: xml

     <ParameterList name="operators">
       <Parameter name="include enthalpy in preconditioner" type="boll" value="true"/>
       <ParameterList name="diffusion operator">
         <ParameterList name="matrix">
           <Parameter name="discretization primary" type="string" value="mdf: optimized for monotonicity"/>
           <Parameter name="discretization secondary" type="string" value="mfd: optimized for sparsity"/>
           <Parameter name="schema" type="Array(string)" value="{face, cell}"/>
           <Parameter name="preconditioner schema" type="Array(string)" value="{face}"/>
           <Parameter name="gravity" type="bool" value="false"/>
           <Parameter name="upwind method" type="string" value="standard: cell"/> 
         </ParameterList>
         <ParameterList name="preconditioner">
           <Parameter name="discretization primary" type="string" value="mfd: optimized for monotonicity"/>
           <Parameter name="discretization secondary" type="string" value="mfd: optimized for sparsity"/>
           <Parameter name="schema" type="Array(string)" value="{face, cell}"/>
           <Parameter name="preconditioner schema" type="Array(string)" value="{face}"/>
           <Parameter name="gravity" type="bool" value="true"/>
           <Parameter name="newton correction" type="string" value="approximate jacobian"/>
           <Parameter name="upwind method" type="string" value="standard: cell"/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>

This example uses cell-centered discretization for 


Advection operator
``````````````````

This section to be written.

.. code-block:: xml

   <ParameterList name="operators">  <!-- parent list -->
     <ParameterList name="advection operator">
       <Parameter name="discretization primary" type="string" value="upwind"/>
     <Parameter name="reconstruction order" type="int" value="0"/>
   </ParameterList>


Coupled process kernels
=======================

Coupling of process kernels requires additional parameters for PK 
described above.


Flow and Energy PK
------------------

The conceptual PDE model of the coupled flow and energy equations is

.. math::
  \begin{array}{l}
  \frac{\partial \theta}{\partial t} 
  =
  - \boldsymbol{\nabla} \cdot (\eta_l \boldsymbol{q}_l)
  - \boldsymbol{\nabla} \cdot (\phi s_g \tau_g D_g \boldsymbol{\nabla} X_g) + Q_1,
  \quad
  \boldsymbol{q}_l 
  = -\frac{\boldsymbol{K} k_r}{\mu} 
  (\boldsymbol{\nabla} p - \rho_l \boldsymbol{g}) \\
  %
  \frac{\partial \varepsilon}{\partial t} 
  =
  \boldsymbol{\nabla} \cdot (\kappa \nabla T) -
  \boldsymbol{\nabla} \cdot (\eta_l H_l \boldsymbol{q}_l) + Q_2
  \end{array}

In the first equation,
:math:`\theta` is total water content,
:math:`\eta_l` is molar density of liquid,
:math:`\rho_l` is fluid density,
:math:`Q_1` is source or sink term,
:math:`\boldsymbol{q}_l` is the Darcy velocity,
:math:`k_r` is relative permeability,
:math:`\boldsymbol{g}` is gravity,
:math:`\phi` is porosity,
:math:`s_g` is gas saturation (water vapor),
:math:`\tau_g` is tortuosity of gas,
:math:`D_g` is diffusion coefficient,
and :math:`X_g` is molar fraction of water in the gas phase.
We define 

.. math::
   \theta = \phi (s_g \eta_g X_g + s_l \eta_l)

where
:math:`s_l` is liquid saturation,
and :math:`\eta_g` is molar density of gas.

In the second equation,
:math:`\varepsilon` is the internal energy,
:math:`Q_2` is source or sink term,
:math:`\kappa` is thermal conductivity,
:math:`H_l` is molar enthalphy of liquid,
and :math:`T` is temperature.
We define 

.. math::
   \varepsilon = \phi (\eta_l s_l U_l + \eta_g s_g U_g) + 
   (1 - \phi) \rho_r c_r T

where
:math:`U_l` is molar internal energy of liquid,
:math:`U_g` is molar internal energy of gas (water vapor),
:math:`\rho_r` is rock density,
and :math:`c_r` is specific heat of rock.


Diffusion operator
..................

.. code-block:: xml

   <ParameterList name="Flow">  <!-- parent lists -->
   <ParameterList name="operator"> 
   <ParameterList name="diffusion operator">
     <ParameterList name="vapor matrix">
       <Parameter name="discretization primary" type="string" value="mfd: optimized for sparsity"/>
       <Parameter name="discretization secondary" type="string" value="mfd: optimized for sparsity"/>
       <Parameter name="schema" type="Array(string)" value="{face, cell}"/>
       <Parameter name="nonlinear coefficient" type="string" value="standard: cell"/>
       <Parameter name="exclude primary terms" type="bool" value="false"/>
       <Parameter name="scaled constraint equation" type="bool" value="false"/>
       <Parameter name="gravity" type="bool" value="false"/>
       <Parameter name="newton correction" type="string" value="none"/>
     </ParameterList>
   </ParameterList>
   </ParameterList>
   </ParameterList>


Generic capabilities
====================

Collection of generic tools used by PKs.

Operators
---------

Operators are discrete forms of linearized PDEs operators.
They form a layer between physical process kernels and solvers
and include diffusion, advection, and source operators.
A PK decides which collection of operators must be used to build a preconditioner.

Operators use a few generic tools that are generic in nature and can be used 
independently by PKs. 
The list includes reconstruction and limiting algorithms. 


Diffusion operator
..................

* `"OPERATOR_NAME`" [list] a PK specific name for the diffusion operator.

  * `"discretization primary`" [string] specifies an advanced discretization method that
    has useful properties under some a priori conditions on the mesh and/or permeability tensor.
    The available options are `"mfd: optimized for sparsity`", `"mfd: optimized for monotonicity`",
    `"mfd: default`", `"mfd: support operator`", `"mfd: two-point flux approximation`",
    and `"fv: default`". 
    The first option is recommended for general meshes.
    The second option is recommended for orthogonal meshes and diagonal absolute 
    permeability tensor. 

  * `"discretization secondary`" [string] specifies the most robust discretization method
    that is used when the primary selection fails to satisfy all a priori conditions.

  * `"nonlinear coefficient`" [string] specifies a method for treating nonlinear diffusion
    coefficient, if any. Available options are `"upwind: face`", `"divk: cell-face`" (default),
    `"divk: face`", `"standard: cell`", `"divk: cell-face-twin`" and `"divk: cell-grad-face-twin`".
    Symmetry preserving methods are the divk-family of methods and the classical cell-centered
    method (`"standard: cell`"). The first part of the name indicates the base scheme.
    The second part (after the semi-column) indicates required components of the composite vector
    that must be provided by a physical PK.

  * `"schema`" [Array(string)] defines the operator stencil. It is a collection of 
    geometric objects.

  * `"preconditioner schema`" [Array(string)] defines the preconditioner stencil.
    It is needed only when the default assembling procedure is not desirable. If skipped,
    the `"schema`" is used instead. 

  * `"gravity`" [bool] specifies if flow is driven also by the gravity.

  * `"gravity term discretization`" [string] selects a model for discretizing the 
    gravity term. Available options are `"hydraulic head`" [default] and `"finite volume`". 
    The first option starts with equation for the shifted solution, i.e. the hydraulic head,
    and derives gravity discretization by the reserve shifting.
    The second option is based on the divergence formula.

  * `"newton correction`" [string] specifies a model for non-physical terms 
    that must be added to the matrix. These terms represent Jacobian and are needed 
    for the preconditioner. Available options are `"true jacobian`" and `"approximate jacobian`".

  * `"scaled constraint equation`" [bool] rescales flux continuity equations on mesh faces.
    These equations are divided by the nonlinear coefficient. This option allows us to 
    treat the case of zero nonlinear coefficient. At moment this feature does not work 
    with non-zero gravity term. Default is *false*.

  * `"constraint equation scaling cutoff"`" [double] specifies the cutoff value for
    applying rescaling strategy described above.  

  * `"consistent faces`" [list] may contain a `"preconditioner`" and
    `"linear operator`" list (see sections Preconditioners_ and LinearSolvers_
    respectively).  If these lists are provided, and the `"discretization
    primary`" is of type `"mfd: *`", then the diffusion method
    UpdateConsistentFaces() can be used.  This method, given a set of cell
    values, determines the faces constraints that satisfy the constraint
    equation in MFD by assembling and inverting the face-only system.  This is
    not currently used by any Amanzi PKs.

.. code-block:: xml

    <ParameterList name="OPERATOR_NAME">
      <Parameter name="discretization primary" type="string" value="mfd: optimized for monotonicity"/>
      <Parameter name="discretization secondary" type="string" value="mfd: two-point flux approximation"/>
      <Parameter name="schema" type="Array(string)" value="{face, cell}"/>
      <Parameter name="preconditioner schema" type="Array(string)" value="{face}"/>
      <Parameter name="gravity" type="bool" value="true"/>
      <Parameter name="gravity term discretization" type="string" value="hydraulic head"/>
      <Parameter name="upwind method" type="string" value="standard: cell"/>
      <Parameter name="newton correction" type="string" value="true jacobian"/>

      <ParameterList name="consistent faces">
        <ParameterList name="linear solver">
          ...
        </ParameterList>
        <ParameterList name="preconditioner">
          ...
        </ParameterList>
      </ParameterList>
    </ParameterList>

This example creates a p-lambda system, i.e. the pressure is
discretized in mesh cells and on mesh faces. 
The preconditioner is defined on faces only, i.e. cell-based unknowns
are eliminated explicitly and the preconditioner is applied to the
Schur complement.


Advection operator
..................

This section is under construction.

* `"OPERATOR_NAME`" [list] a PK specific name for the advection operator.

  * [WIP] `"discretization primary`" defines a discretization method. The only available option is `"upwind`".

  * [WIP] `"reconstruction order`" defines accuracy of this discrete operator.

.. code-block:: xml

  <ParameterList name="OPERATOR_NAME">
    <Parameter name="discretization primary" type="string" value="upwind"/>
    <Parameter name="reconstruction order" type="int" value="0"/>
  </ParameterList>


Reconstruction and limiters
...........................

A reconstruction of discrete fields is used to increase accuracy of discrete models.
The reconstruction can be either unconstrained or limited. 
Amanzi supports a variety of state-of-the-art reconstruction and limiting algorithms 
and their extensions for various PKs.

* `"reconstruction`" [list] describes parameters used by reconstruction algorithms.

 * [WIP] `"method`" [string] specifies a reconstruction method. Available option is
   `"cell-based`" (default).

 * [WIP] `"polynomial order`" [int] defines the polynomial order of a reconstructed function. 
   Default is 1.

 * `"limiter`" [string] specifies limiting method. Available options are 
   `"Barth-Jespersen`" (default), `"tensorial`", and `"Kuzmin`". 

 * `"limiter extension for transport`" [bool] adds additional corrections to 
   limiters required by the transport PK. Default value is *false*.

.. code-block:: xml

  <ParameterList name="reconstruction">
    <Parameter name="method" type="string" value="cell-based"/>
    <Parameter name="order" type="int" value="1"/>
    <Parameter name="limiter" type="string" value="tensorial"/>
    <Parameter name="limiter extension for transport" type="bool" value="false"/>
  </ParameterList>


Functions
---------

To set up non-trivial boundary conditions and/or initial fields, Amanzi
supports a few mathematical functions. 
New function types can added easily.
Each function is defined by a list:

.. code-block:: xml

  <ParameterList name="NAME">
    function-specification
  </ParameterList>

The parameter list name string NAME is arbitrary and meaningful only to the
parent parameter list.
This list is given as input to the Amanzi::FunctionFactory::Create
method which instantiates a new Amanzi::Function object.
The function-specification is one of the following parameter lists.


Constant function
.................

Constant function is defined as `f(x) = a`, for all `x`. 
The specification of this function needs only one parameter.
For example, when `a = 1`, we have:

.. code-block:: xml

  <ParameterList name="function-constant">
    <Parameter name="value" type="double" value="1.0"/>
  </ParameterList>
  

Tabular function
................

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
This behavior can be changed using parameter *forms*.
This parameter is optional.
If specified it must be an array of length equal to one less than the length 
of *x values*.  
Each value in *forms* is either *linear* to indicate linear interpolation on that 
interval, *constant* to use the left endpoint value for that interval, or *FUNCTION*
to indicate an arbitrary user function, usually a math function. 
The default value for parameter *x coordinate* is *t*.

.. code-block:: xml

  <ParameterList name="function-tabular">
    <Parameter name="x values" type="Array(double)" value="{0.0, 1.0, 2.0, 3.0}"/>
    <Parameter name="x coordinate" type="string" value="t"/>
    <Parameter name="y values" type="Array(double)" value="{0.0, 1.0, 2.0, 2.0}"/>
    <Parameter name="forms" type="Array(string)" value="{linear, constant, USER_FUNC}"/>

    <ParameterList name="USER_FUNC">
      <ParameterList name="function-standard-math">
        <Parameter name="operator" type="string" value="sqrt"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>
  
The example defines function that is zero on interval :math:`(-\infty,\,0]`,
linear on interval :math:`(0,\,1]`, constant (`f(x)=1`) on interval :math:`(1,\,2]`, 
square root of `t` on interval :math:`(2,\,3]`,
and constant (`f(x)=2`) on interval :math:`(3,\,\infty]`.
The parameter *x coordinate* defines whether the *x values* refers to time *t*,
x-coordinate *x*, y-coordinate *y*, or z-coordinate *z*.


Bilinear function
.................

The bilinear function provides an extension of the linear form of the tabular function 
to a function with 2 variables `f(x,y)`.
A 2x2 matrix of values for `f(x,y)` and arrays of associated values for `x`
and `y` are read in from datasets in an HDF5 file. The dataset headers are indicated
by parameters *row header*, *column header*, and *value header* for `x`, `y`, 
and `f(x,y)`, respectively. The `x` and `y` arrays in the HDF5 file are expected to be
strictly increasing.
The parameters *row coordinate* and *column coordinate* define the model 
coordinate for `x` and `y` in the function, respectively, where
`t` refers to time, `x` to the x-coordinate, `y` to the y-coordinate, 
and `z` to the z-coordinate. 

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
....................

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
...................

A generic polynomial function is given by the following expression:

.. math::
  f(x) = \sum_{j=0}^n c_j (x - x_0)^{p_j}

where :math:`c_j` are coefficients of monomials,
:math:`p_j` are integer exponents, and :math:`x_0` is the reference point.
Here is an example of a quartic polynomial:

.. code-block:: xml

  <ParameterList name="function-polynomial">
    <Parameter name="coefficients" type="Array(double)" value="{1.0, 1.0}"/>
    <Parameter name="exponents" type="Array(int)" value="{2, 4}"/>
    <Parameter name="reference point" type="double" value="0.0"/>
  </ParameterList>
  

Multi-variable linear function
..............................

A multi-variable linear function is formally defined by
 
.. math::
  f(x) = y_0 + \sum_{j=0}^{n-1} g_j (x_j - x_{0,j}) 

with the constant term :math:`y_0` and gradient :math:`g_0,\, g_1\,..., g_{n-1}`.
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
..................

A separable function is defined as the product of other functions such as

.. math::
  f(x_0, x_1,...,x_{n-1}) = f_1(x_0)\, f_2(x_1,...,x_{n-1})

where :math:`f_1` is defined by the list *function1*, and 
:math:`f_2` by the list *function2*:

.. code-block:: xml

  <ParameterList name="function-separable">
    <ParameterList name="function1">
      function-specification
    </ParameterList>
    <ParameterList name="function2">
      function-specification
    </ParameterList>
  </ParameterList>


Standard math functions
.......................

Amanzi supports a set of standard functions `f(x) = f(x[0])`. 
In Amanzi, the first index of vector `x` corresponds to time.
These functions allow to set up non-trivial time-dependent boundary conditions 
which increases a set of analytic solutions that can be used in convergence 
analysis tests.

* `"operator`" [string] specifies the name of a standard mathematical function.
  Available options are `"cos`", `"sin`", `"tan`", `"acos`", `"asin`", `"atan`", 
  `"cosh`", `"sinh`", `"tanh`", `"exp`", `"log`", `"log10`", `"sqrt`", `"ceil`",
  `"fabs`", `"floor`", `"mod`", and `"pow`".

* `"amplitude`" [double] specifies a multiplication factor `a` in formula `a f(x)`. 
  The multiplication factor is ignored by function `mod`. Default value is 1.

* `"parameter`" [double] specifies additional parameter `p` for math functions 
  with two arguments. These functions are `"a pow(x[0], p)`" and `"a mod(x[0], p)`".
  Defualt value is 0.

* `"shift`" [double] specifies shift of the function argument. Default is 0.

.. code-block:: xml

  <ParameterList name="function-standard-math">
    <Parameter name="operator" type="string" value="sqrt"/>
    <Parameter name="amplitude" type="double" value="1e-7"/>
    <Parameter name="parameter" type="double" value="0.5"/>
    <Parameter name="shift" type="double" value="0.1"/>
  </ParameterList>

This example defines function `1e-7 sqrt(t-0.1)`.


Additive function
.................

To increase calculus of standard math functions, we support a few basic operations
with them. The first one is the sum of two functions, `f(t) = f1(t) + f2(t)`.
This function requires two sublists *function1* and *function2*.

.. code-block:: xml

  <ParameterList name="function-additive">
    <ParameterList name="function1">
      <ParameterList name="function-standard-math">
        <Parameter name="operator" type="string" value="sqrt"/>
        <Parameter name="parameter" type="double" value="0.5"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="function2">
      <ParameterList name="function-standard-math">
        <Parameter name="operator" type="string" value="sin"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>

This example defines function `srqt(t) + sin(t)`.


Composition function
....................

To increase calculus of standard math functions, we support a few basic operations
with them. The second one is the composition of two functions, `f(t) = f1(f2(t))`.
This function requires two sublists *function1* and *function2*.

.. code-block:: xml

  <ParameterList name="function-composition">
    <ParameterList name="function1">
      <ParameterList name="function-standard-math">
        <Parameter name="operator" type="string" value="sqrt"/>
        <Parameter name="parameter" type="double" value="0.5"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="function2">
      <ParameterList name="function-linear">
        <Parameter name="y0" type="double" value="1.0"/>
        <Parameter name="gradient" type="Array(double)" value="{1.0, 2.0, 1.0}"/>
        <Parameter name="x0" type="Array(double)" value="{3.0, 2.0, 1.0}"/>
    </ParameterList>
  </ParameterList>

In two dimensions, this example defines function `srqt(1 + (t-3) + 2(x-2) + (y-1))`.
In three dimension, we have to add one additional argument to parameters *gradient* and *x0*.


Multiplicative function
.......................

To increase calculus of standard math functions, we support a few basic operations
with them. The third one is the multiplication of two functions, `f(t) = f1(t) * f2(t)`.
This function requires two sublists *function1* and *function2*.

.. code-block:: xml

  <ParameterList name="function-multiplicative">
    <ParameterList name="function1">
      <ParameterList name="function-standard-math">
        <Parameter name="operator" type="string" value="sqrt"/>
        <Parameter name="parameter" type="double" value="0.5"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="function2">
      <ParameterList name="function-standard-math">
        <Parameter name="operator" type="string" value="sin"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>

This example defines function `srqt(t) * sin(t)`.


Time functions
--------------

Boundary condition functions utilize a parameterized model for time variations that is either piecewise constant or piecewise linear.  For example:

.. code-block:: xml

   <Parameter name="times" type="Array(double)" value="{1, 2, 3}"/>
   <Parameter name="time values" type="Array(double)" value="{10, 20, 30}"/>
   <Parameter name="time functions" type="Array(string)" value="{Constant, Linear}"/>    

This defines four time intervals: (-inf,1), (1,2), (2,3), (3,+inf).  
By assumption the function is constant over the first and last intervals.
The remaining two intervals are specified by the parameter *time functions*. 
Thus, the value here is 10 anytime prior to `t=2`. The value increases linearly from 10 to 
20 over the interval `t=2` to `t=3`, and then is constant at 30 for `t>3`.


Solvers
=======

This section describes generic solvers and preconditioners that can be used
by various PKs.

.. _LinearSolvers:

Linear solvers
--------------

This list contains sublists for various linear solvers such as PCG, GMRES, and NKA.

* `"preconditioner`" [string] is name in the list of preconditioners. If it is missing, 
  the identity preconditioner is employed.

* `"iterative method`" [string] defines a Krylov-based method. The available options
  include `"pcg`" and `"gmres`".

* `"xxx parameters`" [list] provides parameters for the iterative method specified 
  by variable `"iterative method`".
 
.. code-block:: xml

   <ParameterList>  <!-- parent list -->
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
   </ParameterList>

The names *GMRES with HYPRE AMG* and similar are chosen by the user.


Generalized minimal residuals (GMRES)
.....................................

Not all scientists know that idea of GMRES method was formulated first in 1968.
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

   <ParameterList name="GMRES with HYPRE AMG">  <!-- parent list -->
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
   </ParameterList>


Preconditioner conjugate gradient (PCG)
.......................................

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

  <ParameterList name="PCG with HYPRE AMG">  <!-- parent list -->
    <ParameterList name="pcg parameters">
      <Parameter name="error tolerance" type="double" value="1e-12"/>
      <Parameter name="maximum number of iterations" type="int" value="400"/>
      <Parameter name="convergence criteria" type="Array(string)" value="{relative residual,make one iteration}"/>
      <Parameter name="overflow tolerance" type="double" value="3.0e+50"/>

      <ParameterList name="VerboseObject">
        <Parameter name="Verbosity Level" type="string" value="high"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>


Newton-Krylov acceleration (NKA)
................................

This is a variation of the GMRES solver for nonlinear problems. 
Internal parameters for NKA include

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

  <ParameterList name="NKA">  <!-- parent list -->
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
  </ParameterList>


Nonlinear solvers
-----------------

Amanzi supports a few nonlinear solvers. 
Typically, a process kernel uses a factory to select a nonlinear solver.
This factory uses parameter *solver type* to find parameters for 
the selected solver.


Newton-Krylov acceleration (NKA)
................................

* `"nonlinear tolerance`" [double] defines the required error tolerance. 
  The error is calculated by a PK. Default is 1e-6. 

* `"monitor`" [string] specifies control of the nonlinear residual. The available 
  options are `"monitor update`" (default), `"monitor residual`", and 
  `"monitor preconditioned residual`".

* `"limit iterations`" [int] defines the maximum allowed number of iterations.
  Default is 20.

* `"diverged tolerance`" [double] defines the error level indicating divergence 
  of the solver. The error is calculated by a PK. Default is 1e+10.

* `"diverged l2 tolerance`" [double] defines another way to identify divergence
  of the solver. If the relative L2 norm of the solution increment is above this
  value, the solver is terminated. Default is 1e+10.

* `"diverged pc tolerance`" [double] defines another way to identify divergence
  of the solver. If the relative maximum norm of the solution increment (with respect
  to the initial increment) is above this value, the solver is terminated.
  Default is 1e+10.

* `"diverged residual tolerance`" [double] defines another way to identify divergence
  of the solver. If the relative L2 norm of the residual (with respect
  to the initial residual) is above this value, the solver is terminated.
  Default is 1e+10.

* `"max du growth factor`" [double] allows the solver to identify divergence 
  pattern on earlier iterations. If the maximum norm of the solution increment
  changes drastically on two consecutive iterations, the solver is terminated.
  Default is 1e+5.

* `"max error growth factor`" [double] defines another way to identify divergence 
  pattern on earlier iterations. If the PK-specific error changes drastically on 
  two consecutive iterations, the solver is terminated. Default is 1e+5.

* `"max divergent iterations`" [int] defines another way to identify divergence
  pattern on earlier iterations. If the maximum norm of the solution increment grows 
  on too many consecutive iterations, the solver is terminated. Default is 3.

* `"modify correction`" [bool] allows a PK to modify the solution increment.
  One example is a physics-based clipping of extreme solution values. Default is *false*.

* `"lag iterations`" [int] delays the NKA acceleration, but updates the Krylov space.
  Default is 0.

* `"max nka vectors`" [int] defines the maximum number of consecutive vectors used for 
  a local space. Default is 10.

* `"nka vector tolerance`" [int] defines the minimum allowed orthogonality between vectors in 
  the local space. If a new vector does not satisfy this requirement, the space is modified. 
  Default is 0.05.

.. code-block:: xml

   <Parameter name="solver type" type="string" value="nka"/>
   <ParameterList name="nka parameters">
     <Parameter name="nonlinear tolerance" type="double" value="1.0e-06"/>
     <Parameter name="monitor" type="string" value="monitor update"/>
     <Parameter name="limit iterations" type="int" value="20"/>
     <Parameter name="diverged tolerance" type="double" value="1.0e+10"/>
     <Parameter name="diverged l2 tolerance" type="double" value="1.0e+10"/>
     <Parameter name="diverged pc tolerance" type="double" value="1.0e+10"/>
     <Parameter name="diverged residual tolerance" type="double" value="1.0e+10"/>
     <Parameter name="max du growth factor" type="double" value="1.0e+03"/>
     <Parameter name="max error growth factor" type="double" value="1.0e+05"/>
     <Parameter name="max divergent iterations" type="int" value="3"/>
     <Parameter name="max nka vectors" type="int" value="10"/>
     <Parameter name="nka vector tolerance" type="double" value="0.05"/>
     <Parameter name="modify correction" type="bool" value="false"/>
     <Parameter name="lag iterations" type="int" value="0"/>

     <ParameterList name="VerboseObject">
       <Parameter name="Verbosity Level" type="string" value="high"/>
     </ParameterList>
   </ParameterList>


Anderson acceleration (AA)
..........................

This is a variation of the GMRES solver for nonlinear problems. 
Internal parameters for AA include

* `"nonlinear tolerance`" [double] defines the required error tolerance. 
  The error is calculated by a PK. Default is 1e-6. 

* `"limit iterations`" [int] defines the maximum allowed number of iterations.
  Default is 20.

* `"diverged tolerance`" [double] defines the error level indicating divergence 
  of the solver. The error is calculated by a PK. Default is 1e+10.

* `"diverged l2 tolerance`" [double] defines another way to identify divergence
  of the solver. If the relative L2 norm of the solution increment is above this
  value, the solver is terminated. Default is 1e+10.

* `"diverged pc tolerance`" [double] defines another way to identify divergence
  of the solver. If the relative maximum norm of the solution increment (with respect
  to the initial increment) is above this value, the solver is terminated.
  Default is 1e+10.

* `"max du growth factor`" [double] allows the solver to identify divergence 
  pattern on earlier iterations. If the maximum norm of the solution increment
  changes drastically on two consecutive iterations, the solver is terminated.
  Default is 1e+5.

* `"max divergent iterations`" [int] defines another way to identify divergence
  pattern on earlier iterations. If the maximum norm of the solution increment grows 
  on too many consecutive iterations, the solver is terminated. Default is 3.

* `"max aa vectors`" [int] defines the maximum number of consecutive vectors used for 
  a local space. Default is 10.

* `"modify correction`" [bool] allows a PK to modify the solution increment.
  One example is a physics-based clipping of extreme solution values. Default is *false*.

* `"relaxation parameter`" [double] damping factor for increment. Default is 1.

.. code-block:: xml

  <ParameterList name="AA">  <!-- parent list -->
    <ParameterList name="aa parameters">
      <Parameter name="nonlinear tolerance" type="double" value="1e-5"/>
      <Parameter name="limit iterations" type="int" value="30"/>
      <Parameter name="diverged tolerance" type="double" value="1e+10"/>
      <Parameter name="diverged l2 tolerance" type="double" value="1e+10"/>
      <Parameter name="diverged pc tolerance" type="double" value="1e+10"/>
      <Parameter name="max du growth factor" type="double" value="1e+5"/>
      <Parameter name="max divergent iterations" type="int" value="3"/>
      <Parameter name="max aa vectors" type="int" value="10"/>
      <Parameter name="modify correction" type="bool" value="false"/>
      <Parameter name="relaxation parameter" type="double" value="0.75"/>
    </ParameterList>
  </ParameterList>


Newton
......

The classical Newton method works well for cases where Jacobian is available and
corresponds to a stable (e.g. upwind) discretization.

* `"nonlinear tolerance`" [double] defines the required error tolerance. 
  The error is calculated by a PK. Default is 1e-6. 

* `"monitor`" [string] specifies control of the nonlinear residual. The available 
  options are `"monitor update`" (default) and `"monitor residual`".

* `"limit iterations`" [int] defines the maximum allowed number of iterations.
  Default is 50.

* `"diverged tolerance`" [double] defines the error level indicating divergence 
  of the solver. The error is calculated by a PK. Default is 1e+10.

* `"max du growth factor`" [double] allows the solver to identify divergence 
  pattern on earlier iterations. If the maximum norm of the solution increment
  changes drastically on two consecutive iterations, the solver is terminated.
  Default is 1e+5.

* `"max error growth factor`" [double] defines another way to identify divergence 
  pattern on earlier iterations. If the PK-specific error changes drastically on 
  two consecutive iterations, the solver is terminated. Default is 1e+5.

* `"max divergent iterations`" [int] defines another way to identify divergence
  pattern on earlier iterations. If the maximum norm of the solution increment grows 
  on too many consecutive iterations, the solver is terminated. Default is 3.

* `"modify correction`" [bool] allows a PK to modify the solution increment.
  One example is a physics-based clipping of extreme solution values. Default is *true*.

* `"stagnation iteration check`" determines the number of iterations before the
  stagnation check is turned on. The stagnation happens when the current L2-error
  exceeds the initial L2-error. Default is 8.

.. code-block:: xml

   <Parameter name="solver type" type="string" value="Newton"/>
   <ParameterList name="Newton parameters">
     <Parameter name="nonlinear tolerance" type="double" value="1.0e-05"/>
     <Parameter name="diverged tolerance" type="double" value="1.0e+10"/>
     <Parameter name="max du growth factor" type="double" value="1.0e+03"/>
     <Parameter name="max divergent iterations" type="int" value="3"/>
     <Parameter name="limit iterations" type="int" value="20"/>
     <Parameter name="modify correction" type="bool" value="true"/>
   </ParameterList>


Jacobian-free Newton-Krylov (JFNK)
..................................

JFNK is the example of an inexact Newton solver. 
It requires three sublists for a nonlinear solver (NKA, Newton, etc), 
a preconditioner, and a linear operator that uses this preconditioner.
We describe parameters of the second sublist only.

* `"typical solution value`" [double] Default is 1.

* `"nonlinear solver`" [list] specifies the base nonlinear solvers.

* `"linear operator`" [list] specifies the linear solver for inverting 
  the approximate Jacobian.

* `"finite difference epsilon`" [double] defines the base finite difference epsilon.
  Default is 1e-8.

* `"method for epsilon`" [string] defines a method for calculating finite difference epsilon.
  Available option is `"Knoll-Keyes`", `"Knoll-Keyes L2`", `"Brown-Saad`".

.. code-block:: xml

   <Parameter name="solver type" type="string" value="JFNK"/>
     <ParameterList name="JFNK parameters">
       <Parameter name="typical solution value" type="double" value="1.0"/>

       <ParameterList name="JF matrix parameters">
         <Parameter name="finite difference epsilon" type="double" value="1.0e-8"/>
         <Parameter name="method for epsilon" type="string" value="Knoll-Keyes L2"/>
       </ParameterList>

       <ParameterList name="nonlinear solver">
         <Parameter name="solver type" type="string" value="Newton"/>
         <ParameterList name="Newton parameters">
           ...
         </ParameterList>
       </ParameterList>

       <ParameterList name="linear operator">
         <Parameter name="iterative method" type="string" value="gmres"/>
         <ParameterList name="gmres parameters">
           ...
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>


Preconditioners
---------------

This sublist contains entries for various
preconditioners required by a simulation. At the moment, we support Trilinos multilevel
preconditioner, Hypre BoomerAMG preconditioner, ILU preconditioner, Euclid ILU
preconditioner, and identity preconditioner. 

* `"type`" [string] defines preconditioner name.

* `"xxx parameters`" [list] provides parameters for the preconditioner specified 
  by variable `"type`".
 
.. code-block:: xml

   <ParameterList>  <!-- parent list -->
     <ParameterList name="Preconditioners">
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
   </ParameterList>

Names *TRILINOS ML*, *HYPRE AMG*, and *BLOCK ILU* are chosen by the user.


Hypre's algebraic multigrid (AMG)
.................................

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

   <ParameterList name="HYPRE AMG">  <!-- parent list -->
     <ParameterList name="boomer amg parameters">
       <Parameter name="tolerance" type="double" value="0.0"/>
       <Parameter name="smoother sweeps" type="int" value="3"/>
       <Parameter name="cycle applications" type="int" value="5"/>
       <Parameter name="coarsen type" type="int" value="0"/>
       <Parameter name="strong threshold" type="double" value="0.5"/>
       <Parameter name="relaxation type" type="int" value="3"/>
       <Parameter name="verbosity" type="int" value="0"/>
     </ParameterList>
   </ParameterList>


Euclid ILU
..........

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

   <ParameterList name="MY EUCLID">  <!-- parent list -->
     <ParameterList name="euclid parameters">
       <Parameter name="ILU(k) fill level" type="int" value="6"/>
       <Parameter name="ILUT drop tolerance" type="double" value="0.01"/>
       <Parameter name="rescale rows" type="bool" value="true"/>
       <Parameter name="verbosity" type="int" value="0"/>
     </ParameterList>
   </ParameterList>


Trilinos ML
...........

Internal parameters for Trilinos ML include

.. code-block:: xml

   <ParameterList name="MY ML">  <!-- parent list -->
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
   </ParameterList>


Block ILU
.........

The internal parameters for block ILU are as follows:

.. code-block:: xml

   <ParameterList name="MY ILU">  <!-- parent list -->
     <ParameterList name="block ilu parameters">
       <Parameter name="fact: relax value" type="double" value="1.0"/>
       <Parameter name="fact: absolute threshold" type="double" value="0.0"/>
       <Parameter name="fact: relative threshold" type="double" value="1.0"/>
       <Parameter name="fact: level-of-fill" type="int" value="0"/>
       <Parameter name="overlap" type="int" value="0"/>
       <Parameter name="schwarz: combine mode" type="string" value="Add"/>
     </ParameterList>
   </ParameterList>


Identity
.........

The identity preconditioner is instantiated if either no preconditioner is
specified or the specified preconditioner list does not exists.


Mesh
====

Amanzi supports both structured and unstructured numerical solution approaches.
This flexibility has a direct impact on the selection and design of the underlying 
numerical algorithms, the style of the software implementations, and, ultimately, 
the complexity of the user-interface.  
This specification format uses and describes the unstructured mesh only.

* `"Mesh`" [list] accepts `"Unstructured`" to indicate the meshing option that Amanzi will use.
  This instructs Amanzi to use data structures provided in the Trilinos or MSTK software frameworks.
  To the extent possible, the discretization algorithms implemented under this option 
  are largely independent of the shape and connectivity of the underlying cells.
  As a result, this option supports an arbitrarily complex computational mesh structure
  that enables users to work with numerical meshes that can be aligned with geometrically
  complex man-made or geostatigraphical features.
  Under this option, the user typically provides a mesh file that was generated with 
  an external software package.
  The following mesh file formats are currently supported: `"Exodus II`" (see example),
  `"MSTK`" (see example), and `"MOAB`" (obsolete).
  Amanzi also provides a rudimentary capability to generate unstructured meshes automatically.

  * `"Unstructured`" [list] accepts instructions to either (1) read or, (2) generate an unstructured mesh.

    * `"Read Mesh File`" [list] accepts name, format of pre-generated mesh file

      * `"File`" [string] name of pre-generated mesh file. Note that in the case of an
        Exodus II mesh file, the suffix of the serial mesh file must be .exo and 
        the suffix of the parallel mesh file must be .par.
        When running in serial the code will read this the indicated file directly.
        When running in parallel and the suffix is .par, the code will instead read
        the partitioned files, that have been generated with a Nemesis tool and
        named as filename.par.N.r where N is the number of processors and r is the rank.
        When running in parallel and the suffix is .exo, the code will partition automatically
        the serial file.
     
      * `"Format`" [string] format of pre-generated mesh file (`"MSTK`", `"MOAB`", or `"Exodus II`")

    * `"Generate Mesh`" [list] accepts parameters of generated mesh (currently only `"Uniform`" supported)

      * `"Uniform Structured`" [list] accepts coordinates defining the extents of simulation domain, and number of cells in each direction.

        * `"Domain Low Coordinate`" [Array(double)] Location of low corner of domain
        * `"Domain High Coordinate`" [Array(double)] Location of high corner of domain
        * `"Number Of Cells`" [Array(int)] the number of uniform cells in each coordinate direction

      * `"Expert`" [list] accepts parameters that control which particular mesh framework is to be used.

        * `"Framework`" [string] one of `"stk::mesh`", `"MSTK`", `"MOAB`" or `"Simple`". 
        * `"Verify Mesh`" [bool] true or false. 

Example of *Unstructured* mesh generated internally:

.. code-block:: xml

   <ParameterList>  <!-- parent list -->
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
   </ParameterList>

Example of *Unstructured* mesh read from an external file:

.. code-block:: xml

   <ParameterList name="Mesh">  <!-- parent list -->
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

Amanzi automatically defines the special region labeled *All*, which is the 
entire simulation domain. Currently, the unstructured framework does
not support the *All* region, but it is expected to do so in the
near future.

User-defined regions are constructed using the following syntax

 * "Regions" [list] can accept a number of lists for named regions (REGION)

   * Shape [list] Geometric model primitive, choose exactly one of the following [see table below]: `"Region: Point`", `"Region: Box`", `"Region: Plane`", `"Region: Labeled Set`", `"Region: Layer`", `"Region: Surface`"

Amanzi supports parameterized forms for a number of analytic shapes, as well as more complex definitions based on triangulated surface files.  

+------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
|  shape functional name       | parameters                              | type(s)                      | Comment                                                                |
+==============================+=========================================+==============================+========================================================================+
| `"Region: Point"`            | `"Coordinate`"                          | Array(double)                | Location of point in space                                             |
+------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
| `"Region: Box"`              | `"Low Coordinate`", `"High Coordinate`" | Array(double), Array(double) | Location of boundary points of box                                     |
+------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
| `"Region: Plane"`            | `"Direction`", `"Location`"             | string, double               | direction: `"X`", `"-X`", etc, and `"Location`" is coordinate value    |
+------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
| `"Region: Polygonal Surface"`| `"Number of points`", `"Points`"        | int, Array(double)           | Number of polygon points and point coordinates in linear array. This   |
|                              |                                         |                              | provides a set of faces with a normal for computing flux               |    
+------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
| `"Region: Logical"`          | `"Operation`", `"RegionList`"           | string, Array(string)        | Operation can be Union, Intersection, Subtraction, Complement          |
+------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
| `"Region: Labeled Set"`      | `"Label`", `"File`",                    | string, string,              | Set per label defined in mesh file (see below)                         |
|                              | `"Format`", `"Entity`"                  | string, string               |  (available for frameworks supporting the `"File`" keyword)            |
+------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
| `"Region: Color Function"`   | `"File`", `"Value`"                     | string, int                  | Set defined by color in a tabulated function file (see below)          |
+------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+


Point
-----

List *Region: Point* defines a point in space. 
Using this definition, cell sets encompassing this point are retrieved inside Amanzi.

.. code-block:: xml

   <ParameterList name="Dnwind150"> <!-- parent list -->
     <ParameterList name="Region: Point">
       <Parameter name="Coordinate" type="Array(double)" value="{-150.0, 0.0, 0.0}"/>
     </ParameterList>
   </ParameterList>


Box
---

List *Region: Box* defines a region bounded by coordinate-aligned
planes. Boxes are allowed to be of zero thickness in only one
direction in which case they are equivalent to planes.

.. code-block:: xml

   <ParameterList name="Well">  <!-- parent list -->
     <ParameterList name="Region: Box">
       <Parameter name="Low Coordinate" type="Array(double)" value="{-5.0,-5.0, -5.0}"/>
       <Parameter name="High Coordinate" type="Array(double)" value="{5.0, 5.0,  5.0}"/>
     </ParameterList>
   </ParameterList>


Plane
-----

List *Region: Plane* defines a plane using a point lying on the plane, called *Location*,
and normal to the plane, called *Direction*.

.. code-block:: xml

   <ParameterList name="Top Section"> <!-- parent list -->
     <ParameterList name="Region: Plane">
       <Parameter name="Location" type="Array(double)" value="{2, 3, 5}"/>
       <Parameter name="Direction" type="Array(double)" value="{1, 1, 0}"/>
       <ParameterList name="Expert Parameters">
         <Parameter name="Tolerance" type="double" value="1.0e-05"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>


Labeled Set
-----------

The list *Region: Labeled Set* defines a named set of mesh entities
existing in an input mesh file. This is the same file that contains
the computational mesh. The name of the entity set is given
by *Label*.  For example, a mesh file in the Exodus II
format can be processed to tag cells, faces and/or nodes with
specific labels, using a variety of external tools.  Regions based
on such sets are assigned a user-defined label for Amanzi, which may
or may not correspond to the original label in the exodus file.
Note that the file used to express this labeled set may be in any
Amanzi-supported mesh format (the mesh format is specified in the
parameters for this option).  The *entity* parameter may be
necessary to specify a unique set.  For example, an Exodus file
requires *Cell*, *Face* or *Node* as well as a label (which is
an integer).  The resulting region will have the dimensionality 
associated with the entities in the indicated set. 

Currently, Amanzi-U only supports mesh files in the Exodus II format.

.. code-block:: xml

   <ParameterList name="Aquifer">
     <ParameterList name="Region: Labeled Set">
       <Parameter name="Entity" type="string" value="Cell"/>
       <Parameter name="File" type="string" value="porflow4_4.exo"/>
       <Parameter name="Format" type="string" value="Exodus II"/>
       <Parameter name="Label" type="string" value="1"/>
     </ParameterList>
   </ParameterList>


Color Function
--------------

The list *Region: Color Function* defines a region based a specified
integer color, *Value*, in a structured color function file,
*File*. The format of the color function file is given below in
the "Tabulated function file format" section. As
shown in the file, the color values may be specified at the nodes or
cells of the color function grid. A computational cell is assigned
the 'color' of the data grid cell containing its cell centroid
(cell-based colors) or the data grid nearest its cell-centroid
(node-based colors). Computational cells sets are then built from
all cells with the specified color *Value*.

In order to avoid, gaps and overlaps in specifying materials, it is
strongly recommended that regions be defined using a single color
function file. 

.. code-block:: xml

   <ParameterList name="SOIL1">
     <ParameterList name="Region: Color Function">
       <Parameter name="File" type="string" value="geology_resamp_2D.tf3"/>
       <Parameter name="Value" type="int" value="1"/>
     </ParameterList>
   </ParameterList>


Polygon
-------

The list *Region: Polygon* defines a polygonal region on which mesh faces and
nodes can be queried. NOTE that one cannot ask for cells in a polygonal surface
region. In 2D, the "polygonal surface" region is a line and is specified by 2 points.
In 3D, the "polygonal surface" region is specified by an arbitrary number of points.
In both cases the point coordinates are given as a linear array. The polygon
can be non-convex.

The polygonal surface region can be queried for a normal. In 2D, the normal is
defined as [Vy,-Vx] where [Vx,Vy] is the vector from point 1 to point 2.
In 3D, the normal of the polygon is defined by the order in which points 
are specified.

.. code-block:: xml

   <ParameterList name="XY pentagon">
     <ParameterList name="Region: Polygon">
       <Parameter name="Number of points" type="int" value="5"/>
       <Parameter name="Points" type="Array(double)" value="{-0.5, -0.5, -0.5, 
                                                              0.5, -0.5, -0.5,
                                                              0.8, 0.0, 0.0,
                                                              0.5,  0.5, 0.5,
                                                             -0.5, 0.5, 0.5}"/>
       <ParameterList name="Expert Parameters">
         <Parameter name="Tolerance" type="double" value="1.0e-3"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>


Logical
-------

The list *Region: Logical* defines logical operations on regions allow for more
advanced region definitions. At this time the Logical Region allows
for logical operations on a list of regions.  In the case of Union
the result is obvious, it is the union of all regions.  Similarly
for Intersection. In the case of Subtraction, subtraction is
performed from the first region in the list.  The Complement is a
special case in that it is the only case that operates on single
region, and returns the complement to it within the domain *Entire Domain*.
Currently, multi-region booleans are not supported in the same expression.

.. code-block:: xml

  <ParameterList name="Lower Layers">
    <ParameterList name="Region: Logical">
      <Parameter name="Operation" type="string" value="Union"/>
      <Parameter name="RegionList" type="Array(string)" value="{Middle1, Middle2, Bottom}"/>
    </ParameterList>
  </ParameterList>


Notes and example
-----------------

* Surface files contain labeled triangulated face sets.  The user is
  responsible for ensuring that the intersections with other surfaces
  in the problem, including the boundaries, are *exact* (*i.e.* that
  surface intersections are *watertight* where applicable), and that
  the surfaces are contained within the computational domain.  If
  nodes in the surface fall outside the domain, the elements they
  define are ignored.

  Examples of surface files are given in the *Exodus II* file 
  format here.

* Region names must NOT be repeated.

.. code-block:: xml

   <ParameterList>  <!-- parent list -->
     <ParameterList name="Regions">
       <ParameterList name="TOP SECTION">
         <ParameterList name="Region: Box">
           <Parameter name="Low Coordinate" type="Array(double)" value="{2, 3, 5}"/>
           <Parameter name="High Coordinate" type="Array(double)" value="{4, 5, 8}"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="MIDDLE SECTION">
         <ParameterList name="Region: Box">
           <Parameter name="Low Coordinate" type="Array(double)" value="{2, 3, 3}"/>
           <Parameter name="High Coordinate" type="Array(double)" value="{4, 5, 5}"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="BOTTOM SECTION">
         <ParameterList name="Region: Box">
           <Parameter name="Low Coordinate" type="Array(double)" value="{2, 3, 0}"/>
           <Parameter name="High Coordinate" type="Array(double)" value="{4, 5, 3}"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="INFLOW SURFACE">
         <ParameterList name="Region: Labeled Set">
           <Parameter name="Label"  type="string" value="sideset_2"/>
           <Parameter name="File"   type="string" value="F_area_mesh.exo"/>
           <Parameter name="Format" type="string" value="Exodus II"/>
           <Parameter name="Entity" type="string" value="Face"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="OUTFLOW PLANE">
         <ParameterList name="Region: Plane">
           <Parameter name="Location" type="Array(double)" value="{0.5, 0.5, 0.5}"/>
           <Parameter name="Direction" type="Array(double)" value="{0, 0, 1}"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="BLOODY SAND">
         <ParameterList name="Region: Color Function">
           <Parameter name="File" type="string" value="F_area_col.txt"/>
           <Parameter name="Value" type="int" value="25"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="FLUX PLANE">
         <ParameterList name="Region: Polygon">
           <Parameter name="Number of points" type="int" value="5"/>
           <Parameter name="Points" type="Array(double)" value="{-0.5, -0.5, -0.5, 
                                                                  0.5, -0.5, -0.5,
                                                                  0.8, 0.0, 0.0,
                                                                  0.5,  0.5, 0.5,
                                                                 -0.5, 0.5, 0.5}"/>
          </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>

In this example, *TOP SESCTION*, *MIDDLE SECTION* and *BOTTOM SECTION*
are three box-shaped volumetric regions. *INFLOW SURFACE* is a
surface region defined in an Exodus II-formatted labeled set
file and *OUTFLOW PLANE* is a planar region. *BLOODY SAND* is a volumetric
region defined by the value 25 in color function file.


Output data
===========

Amanzi uses a few ways to communicate simulation data to the user that includes
a short file with observations and full-scale visualization files.

Observation file
----------------

A user may request any number of specific observations from Amanzi.  Each labeled Observation Data quantity involves a field quantity, a model, a region from which it will extract its source data, and a list of discrete times 
for its evaluation.  The observations are evaluated during the simulation and returned to the calling process through one of Amanzi arguments.

* `"Observation Data`" [list] can accept multiple lists for named observations.

  * `"observation output filename`" [string] user-defined name for the file that the observations are written to.
    The file name can contain relative or absolute path to an *existing* directory only. 

  * `"precision`" [int] defined the number of significant digits. Default is 16.

  * OBSERVATION [list] user-defined label, can accept values for `"variables`", `"functional`",
    `"region`", `"times`", and TSPS (see below).

    * `"variables`" [Array(string)] a list of field quantities taken from the list of 
      available field quantities:

      * volumetric water content [-] (volume water / bulk volume)
      * aqueous saturation [-] (volume water / volume pore space)
      * aqueous pressure [Pa]
      * hydraulic head [m] 
      * drawdown [m] 
      * SOLUTE Aqueous concentration [mol/m^3]
      * SOLUTE gaseous concentration [mol/m^3]
      * x-, y-, z- aqueous volumetric flux [m/s]
      * material id [-]
      * aqueous mass flow rate [kg/s] (when funtional="integral")
      * aqueous volumetric flow rate [m^3/s] (when functional="integral")
      * SOLUTE volumetric flow rate [mol/s] (when functional="integral")

    Observation *drawdown* is calculated with respect to the value registered at the first time
    it was requested.

    * `"functional`" [string] the label of a function to apply to each of the variables
      in the variable list (Function options detailed below)

    * `"region`" [string] the label of a user-defined region

    * `"cycles start period stop`" [Array(int)] the first entry is the start cycle,
      the second is the cycle period, and the third is the stop cycle or -1 in which case 
      there is no stop cycle. A visualization dump shall be written at such cycles that 
      satisfy cycle = start + n*period, for n=0,1,2,... and cycle < stop if stop != -1.0.

    * `"cycles start period stop n`" [Array(int)] if multiple cycles start-period-stop
      parameters are needed, then use these parameters with n=0,1,2,..., and not the single 
      `"cycles start period stop`" parameter.

    * `"cycles`" [Array(int)] an array of discrete cycles that at which a visualization dump shall be written. 

    * `"times start period stop`" [Array(double)] the first entry is the start time,
      the second is the time period, and the third is the stop time or -1 in which case 
      there is no stop time. A visualization dump shall be written at such times that 
      satisfy time = start + n*period, for n=0,1,2,... and time < stop if stop != -1.0.

    * `"times start period stop n`" [Array(double) if multiple start-period-stop parameters
      are needed, then use this these parameters with n=0,1,2,..., and not the 
      single  `"times start period stop`" parameter.

    * `"times`" [Array(double)] an array of discrete times that at which a visualization
      dump shall be written.

    * `"delimiter`" [string] the string used to delimit columns in the observation file
      output, default is `",`".

The following observation functionals are currently supported.
All of them operate on the variables identified.

* `"Observation Data: Point`" returns the value of the field quantity at a point.

* `"Observation Data: Integral`" returns the integral of the field quantity over the region specified.

* `"Observation Data: Extensive Integral`" returns the integral of an extensive variable
  over the region specified.  Note that this should be used over the above Integral when 
  the variable to be integrated is an extensive quantity, i.e. water content or flux.

* `"Observation Data: Minimum`" and `"Observation Data: Maximum`" returns the minimum 
  (respectively maximum) of the field quantity over the region specified.

.. code-block:: xml

   <ParameterList>  <!-- parent list -->
     <ParameterList name="Observation Data">
       <Parameter name="observation output filename" type="string" value="obs_output.out"/>
       <Parameter name="precision" type="int" value="10"/>
       <ParameterList name="SOME OBSERVATION NAME">
         <Parameter name="region" type="string" value="some point region name"/>
         <Parameter name="functional" type="string" value="Observation Data: Point"/>
         <Parameter name="variable" type="string" value="volumetric water content"/>
         <Parameter name="times" type="Array(double)" value="{100000.0, 200000.0}"/>

         <Parameter name="cycles" type="Array(int)" value="{100000, 200000, 400000, 500000}"/>
         <Parameter name="cycles start period stop" type="Array(int)" value="{0, 100, -1}" />

         <Parameter name="times start period stop 0" type="Array(double)" value="{0.0, 10.0, 100.0}"/>
         <Parameter name="times start period stop 1" type="Array(double)" value="{100.0, 25.0, -1.0}"/>
         <Parameter name="times" type="Array(double)" value="{101.0, 303.0, 422.0}"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>


Checkpoint file
---------------

A user may request periodic dumps of Amanzi Checkpoint Data.  
The user has no explicit control over the content of these files, but has the guarantee that 
the Amanzi run will be reproducible (with accuracy determined
by machine round errors and randomness due to execution in a parallel computing environment).
Therefore, output controls for Checkpoint Data are limited to file name generation and writing 
frequency, by numerical cycle number.

* `"Checkpoint Data`" [list] can accept a file name base [string] and cycle data [list] 
  used to generate the file base name or directory base name that is used in writing Checkpoint Data. 

  * `"file name base`" [string] ("checkpoint")
  
  * `"file name digits`" [int] Default value is 5.

  * `"cycles start period stop`" [Array(int)] the first entry is the start cycle, the second 
    is the cycle period, and the third is the stop cycle or -1 in which case there is no stop cycle.
    A visualization dump shall be written at such cycles that satisfy cycle = start + n*period, 
    for n=0,1,2,... and cycle < stop if stop != -1.0.

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

   <ParameterList>  <!-- parent list -->
     <ParameterList name="Checkpoint Data">
       <Parameter name="file name base" type="string" value="chkpoint"/>
       <Parameter name="file name digits" type="int" value="5"/>

       <Parameter name="cycles start period stop" type="Array(int)" value="{0, 100, -1}" />
       <Parameter name="cycles" type="Array(int)" value="{999, 1001}" />

       <Parameter name="times start period stop 0" type="Array(double)" value="{0.0, 10.0, 100.0}"/>
       <Parameter name="times start period stop 1" type="Array(double)" value="{100.0, 25.0, -1.0}"/>
       <Parameter name="times" type="Array(double)" value="{101.0, 303.0, 422.0}"/>
     </ParameterList>
   </ParameterList>

In this example, Checkpoint Data files are written when the cycle number is 
a multiple of 100.


Visualization file
------------------

A user may request periodic writes of field data for the purposes of visualization.  
The user will specify explicitly what is to be included in the file at each snapshot.
Visualization files can only be written at intervals corresponding to the numerical 
time step values or intervals corresponding to the cycle number; writes are controlled by time step cycle number.

* `"Visualization Data`" [list] can accept a file name base [string] and cycle data [list] 
  that is used to generate the file base name or directory base name that is used in writing visualization data.
  It can also accept a set of lists to specify which field quantities to write.
  The file name can contain relative or absolute path to an *existing* directory only. 

  * `"file name base`" [string] ("amanzi_vis")
  
  * `"cycles start period stop`" [Array(int)] the first entry is the start cycle, 
    the second is the cycle period, and the third is the stop cycle or -1 in which case 
    there is no stop cycle. A visualization dump shall be written at such cycles that 
    satisfy cycle = start + n*period, for n=0,1,2,... and cycle < stop if stop != -1.0.

  * `"cycles start period stop n`" [Array(int)] if multiple cycles start-period-stop parameters 
    are needed, then use these parameters with n=0,1,2,..., and not the single 
    `"cycles start period stop`" parameter.

  * `"cycles`" [Array(int)] an array of discrete cycles that at which a visualization dump shall be written. 

  * `"times start period stop`" [Array(double)] the first entry is the start time, 
    the second is the time period, and the third is the stop time or -1 in which case 
    there is no stop time. A visualization dump shall be written at such times that 
    satisfy time = start + n*period, for n=0,1,2,... and time < stop if stop != -1.0.

  * `"times start period stop n`" [Array(double) if multiple start-period-stop parameters 
    are needed, then use this these parameters with n=0,1,2,..., and not the single
    `"times start period stop`" parameter.

  * `"times`" [Array(double)] an array of discrete times that at which a visualization dump shall be written.

  * `"dynamic mesh`" [bool] (false) write mesh data for every visualization dump, 
    this facilitates visualizing deforming meshes.

  * `"write regions`" [list] contains auxiliary fields with region ids to write into the visualization file.

    * `"REGION_NAME`" [Array(string)] the user-defined field name and the list of assigned regions. 
      The first entry in the regions array is marked with the value 1.0 in the array, 
      the second with the value 2.0, and so forth. 
      The code ignores entries in the regions array that are not valid regions that contain cells.

  * `"write partitions`" [bool] (false) if this parameter is true, then write an array into 
    the visualization file that contains the rank number of the processor that owns a mesh cell. 

.. code-block:: xml

   <ParameterList>  <!-- parent list -->
     <ParameterList name="Visualization Data">
       <Parameter name="file name base" type="string" value="chk"/>
  
       <Parameter name="cycles start period stop" type="Array(int)" value="{0, 100, -1}" />
       <Parameter name="cycles" type="Array(int)" value="{999, 1001}" />

       <Parameter name="times start period stop 0" type="Array(double)" value="{0.0, 10.0, 100.0}"/>
       <Parameter name="times start period stop 1" type="Array(double)" value="{100.0, 25.0, -1.0}"/>
       <Parameter name="times" type="Array(double)" value="{101.0, 303.0, 422.0}"/>

       <Parameter name="dynamic mesh" type="bool" value="false"/>

       <ParameterList name="write regions">
         <Parameter name="regions" type="Array(string)" value="{Obs_r1, Obs_r1, Obs_r3}"/>
         <Parameter name="wells" type="Array(string)" value="{Obs_r1}"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>


Walkabout file
--------------

A user may request periodic dumps of Walkabout Data. Output controls for Walkabout Data are limited to file name generation and writing frequency, by numerical cycle number or time.

* `"Walkabout Data`" [list] can accept a file name base [string] and cycle data [list] 
  used to generate the file base name or directory base name that is used in writing Checkpoint Data. 

  * `"file name base`" [string] The file name can contain relative or absolute path to an *existing* 
    directory only.  Default is `"walkabout`".
  
  * `"file name digits`" [int] specify the number of digits that should be appended to the file 
    name for the cycle number. Default is 5.

  * `"cycles start period stop`" [Array(int)] the first entry is the start cycle, 
    the second is the cycle period, and the third is the stop cycle or -1 in which case 
    there is no stop cycle. A visualization dump shall be written at such cycles that 
    satisfy cycle = start + n*period, for n=0,1,2,... and cycle < stop if stop != -1.0.

  * `"cycles start period stop n`" [Array(int)] if multiple cycles start-period-stop parameters 
    are needed, then use these parameters with n=0,1,2,..., and not the single
    `"cycles start period stop`" parameter.

  * `"cycles`" [Array(int)] an array of discrete cycles that at which a visualization dump shall be written. 

  * `"times start period stop`" [Array(double)] the first entry is the start time, 
    the second is the time period, and the third is the stop time or -1 in which case 
    there is no stop time. A visualization dump shall be written at such times that 
    satisfy time = start + n*period, for n=0,1,2,... and time < stop if stop != -1.0.

  * `"times start period stop n`" [Array(double) if multiple start-period-stop parameters 
    are needed, then use this these parameters with n=0,1,2,..., and not the single
    `"times start period stop`" parameter.

  * `"times`" [Array(double)] an array of discrete times that at which a visualization dump shall be written.

.. code-block:: xml

   <ParameterList>  <!-- parent list -->
     <ParameterList name="Walkabout Data">
       <Parameter name="file name base" type="string" value="walkabout"/>
       <Parameter name="file name digits" type="int" value="5"/>

       <Parameter name="cycles start period stop" type="Array(int)" value="{0, 100, -1}" />
       <Parameter name="cycles" type="Array(int)" value="{999, 1001}" />

       <Parameter name="times start period stop 0" type="Array(double)" value="{0.0, 10.0, 100.0}"/>
       <Parameter name="times start period stop 1" type="Array(double)" value="{100.0, 25.0, -1.0}"/>
       <Parameter name="times" type="Array(double)" value="{101.0, 303.0, 422.0}"/>
     </ParameterList>
   </ParameterList>

In this example, walkabout data files are written when the cycle number is 
a multiple of 100.


Input data
==========

This section describes format and purpose of various input files.
In addition it explain how to verify the input information (e.g. regions) using special sublists.


Input analysis
--------------

This list contains data collected by the input parser of a higher-level spec. 

* `"used boundary condition regions`" [Array(string)] provides list of boundary regions 
  for analysis. The simulator will print number of faces and total area of these regions
  if verbosity level is equal to or above *high*.

* `"used source and sink regions`" [Array(string)] provides list of source and sink regions
  for analysis. The simulator will print number of cells and the total volume of these regions
  if verbosity level is equal to or above *high*.

* `"used observation regions`" [Array(string)] provides list of observation regions
  for analysis. The simulator will print number of faces(or cells) and the total area 
  (or volume) of these regions if verbosity level is equal to or above *high*.

.. code-block:: xml

   <ParameterList>  <!-- parent list -->
     <ParameterList name="Analysis">
       <Parameter name="used boundary condition regions" type="Array(string)" value="{region1,region2}"/>
       <Parameter name="used source and sink regions" type="Array(string)" value="{region3,region4}"/>
       <Parameter name="used observation regions" type="Array(string)" value="{region5}"/>
       <ParameterList name="VerboseObject">
         <Parameter name="Verbosity Level" type="string" value="high"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>
  

Tabulated function file format
------------------------------

The following ASCII input file format supports the definition of a tabulated function defined over a grid.
Several XML input Parameters refer to files in this format.
The file consists of the following records (lines).
Each record is on a single line, except for the DATAVAL record which may be split across multiple lines.

1. **DATATYPE**:  An integer value: 0 for integer data, 1 for real data.

  * An integer-valued file is used to define a 'color' function used in the definition of a region.

2. **GRIDTYPE**:  A string that specifies the type of grid used to define the function.
The format of the rest of the file is contingent upon this value.
The currently supported options are uniform rectilinear grids in 1, 2 and 3-D, 
which are indicated by the values `1DCoRectMesh`, `2DCoRectMesh` and `3DCoRectMesh`, 
respectively (names adopted from XDMF).

For the uniform rectilinear grids, the remaining records are as follows.  
Several records take 1, 2 or 3 values depending on the space dimension of the grid.

3. **NXNYNZ**: 3 (or 2, 1) integer values (NX, NY, NZ) giving the number of zones in 
the x, y and z coordinate directions, respectively.

4. **CORNER1**: 3 (or 2, 1) floating point values (X1, Y1, Z1) giving the coordinate 
of the first corner of the domain.

5. **CORNER2**: 3 (or 2, 1) floating point values (X2, Y2, Z2) giving the coordinate 
of the second corner of the domain.  The grid points r_{i,j,k} = (x_i, y_j, z_j) are defined as:

      x_i = X1 + i*(X2-X1)/NX, 0 <= i <= NX

      y_j = Y1 + j*(Y2-Y1)/NY, 0 <= j <= NY

      z_k = Z1 + k*(Z2-Z1)/NZ, 0 <= k <= NZ

The (i,j,k) grid cell is defined by the corner grid points r_{i-1,j-1,k-1} and 
r_{i,j,k}, for 1 <= i <= NX, 1 <= j <= NY, 1 <= k <= NZ. 
Note that the corner points are any pair of opposite corner points; the ordering of grid 
points and cells starts at CORNER1 and ends at CORNER2.

6. **DATALOC**:  An integer value: 0 for cell-based data, 1 for point-based data.


7. **DATACOL**:  An integer (N) giving the number of "columns" in the data.  
This is the number of values per grid cell/point.  
N=1 for a scalar valued function; N>1 for a N-vector valued function.

  * only a single column is currently supported.

8. **DATAVAL**: The values of the function on the cells/points of the grid.  
The values should appear in Fortran array order were the values stored in 
the Fortran array A(N,NX,NY,NZ) (A(N,0:NX,0:NY,0:NZ) for point-based data).  
That is, the column index varies fastest, x grid index next fastest, etc.
    

Example
.......

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










