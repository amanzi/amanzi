==========================================
Amanzi-U Native XML Input Specification V8
==========================================

.. contents:: **Table of Contents**


Overview
========

This is a continuously evolving specification format used by the code developers. 
Its main purpose is to develop and test new capabilities without disruption of end-users.


ParameterList XML
=================

The Amanzi input file is an ASCII text XML-formatted file that must be framed 
at the beginning and end by the following statements:

.. code-block:: xml

  <ParameterList name="transport">
    various parameters and sublists
  </ParameterList>

The value of *name* can be anything (*transport* in this example).  
A ParameterList consists of just two types of entries: Parameter and ParameterList.  
ParameterLists are labeled with *name* [string], while Parameters have a separate 
fields called *name* [string], *type* [string] and *value* [TYPE], where TYPE can 
be any of the following: double, int, bool, string, Array(double), Array(int), 
and Array(string).  
The value of the parameter is given in quotes (e.g. value="2.7e3").  
Array data is specified as a single comma-delimited string bounded by {}'s (e.g. value="{2.4, 2.1, 5.7}").

.. code-block:: xml

  <ParameterList name="transport">
    <Parameter name="cfl" type="double" value="0.9"/>
    <Parameter name="ratio" type="Array(int)" value="{2, 1, 4}"/>
  </ParameterList>

In this example, the list *transport* has parameter *cfl* that is the double with 
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

* User-defined labels are marked with ALL_CAPS in this document.
  In practice, no rules are imposed on these names.

* Lists with too many parameters are described using multiple sections and multiple examples.
  For most examples we show name of the parent sublist.


Naming convention rule
----------------------

It is hard to overestimate importance of a reasonable naming convention rule for efficient
code development and its daily usage in reasearch.

* Camel-case names should *not* be used as names for fixed keywords (parameters and parameter lists).  
  The following is a short list of allowed exceptions. 
 
  * The names created by the user are not fixed/reserved keywords and are exempt from the above
    rule. In this documents, we always prefix user-defined names with 
    the underscore symbol.

  * Proper names such as an individual person, place, or organization, including their derivatives
    *should* be spelled using capital letters. Examples: *van Genuchten m*, *Brocks-Corey lambda*, 
    *Jacobian matrix*, and *Newton correction*.

  * Names of chemical species (inside fixed keywords) should be capitalized. Examples: *CO2*, *H+*.

  * A few well-established abbreviations. Their complete list is here: *PK*, *MPC*, *BDF1*, *EOS*,
    *IEM*, *PFloTran*, *pH*, *TP*. Note that names of linear and nonlinear solvers and preconditioners are 
    not included in this list. Thus, we have to use *pcg*, *gmres*, *nka*, *amg*, *ml*, and *ilu*.

  * Units such as energy [J] and temperature [K].

  * The Hilbert spaces *L2* and *H1*. Note that *L2* and *l2* are different spaces and should be used
    appropriately.

  * Trilinos parameters. There are a few camel-case parameters that
    go directly to Trilinos functions and therefore outside of our control, e.g. *ML output*.


Verbose output
--------------

Output of all components of Amanzi is controlled by a standard verbose 
object list. This list can be inserted in almost any significant
component of this spec to produce a verbose output, see the embedded examples.
If this list is not specified, the default verbosity value is used.

* `"verbosity level`" [string] Available options are *none*, *low*, *medium*, *high*, and *extreme*.
  Option *extreme is used by the developers only. For communication between users and developers, 
  the recommended option is *high*. 

* `"hide line prefix`" [bool] defines prefix for output messages. Default value is *true*.

* `"name`" [string] is the name of the prefix.

* `"write on rank`" [int] is processor rank on which the output is performed. Default is 0.

.. code-block:: xml

   <ParameterList name="verbose object">
     <Parameter name="verbosity level" type="string" value="medium"/>
     <Parameter name="name" type="string" value="my header"/>
     <Parameter name="hide line prefix" type="bool" value="false"/>
     <Parameter name="write on rank" type="int" value="0"/>
   </ParameterList>


Residual debugger
-----------------

Some components (currently just nonlinear solver, this may change)
leverage a *residual debugger* object for writing, to file, residuals,
corrections, and internal iterates of a solve process for solver
debugging/work.  Control of when these iterates are written is
controlled by a few parameters.  This should be written sparingly --
each attempt at a timestep and each cycle is its own file, and writes
its own mesh file, so this should be considered i/o and runtime
expensive.

  * `"cycles start period stop`" [Array(int)] the first entry is the start cycle, 
    the second is the cycle period, and the third is the stop cycle or -1 in which case 
    there is no stop cycle. All iterations shall be written at such cycles that 
    satisfy formula cycle = start + n*period, for n=0,1,2,... and cycle < stop if stop != -1.0.

  * `"cycles start period stop n`" [Array(int)] if multiple cycles start-period-stop parameters 
    are needed, then use these parameters with n=0,1,2,..., and not the single 
    `"cycles start period stop`" parameter.

  * `"cycles`" [Array(int)] an array of single cycles at which all iterations shall be written. 

Note: *cycle* here means the current time integration step and *not* the global cycle.

.. code-block:: xml

  <ParameterList name="BDF1">  <!-- parent list -->
  <ParameterList name="residual debugger">
    <Parameter name="cycles start period stop" type="Array(int)" value="{0,100,-1}"/>
    <Parameter name="cycles" type="Array(int)" value="{999,1001}"/>
  </ParameterList>
  </ParameterList>
   

Units
-----

Amanzi's internal default units are SI units except for the concentration.

* `"concentration`" [string] defines units for concentration. Available options
  are `"molar`" (default) which is `"mol/L`" and `"SI`" which is `"mol/m^3`". 

.. code-block:: xml

  <ParameterList>  <!-- parent list -->
  <ParameterList name="units">
    <Parameter name="length" type="string" value="m"/>
    <Parameter name="time" type="string" value="s"/>
    <Parameter name="mass" type="string" value="kg"/>
    <Parameter name="temperature" type="string" value="K"/>
    <Parameter name="concentration" type="string" value="molar"/>
  </ParameterList>
  </ParameterList>


Cycle driver
============

The new multi-processor cycle driver provides more flexibility
to handle multiphysics process kernels (PKs) and multiple time periods.

* `"component names`" [Array(string)] provides the list of species names.
  It is required for reactive transport.

* `"component molar masses`" [Array(string)] provides the list of 
  molar masses of species. It is required for proper conversion to and from 
  dimensionless units. Default is 1. 

* `"number of liquid components`" [int] is the number of liquid components. 
   
* `"time periods`" [list] contains the list of time periods involved in the simulation.
  The number of time periods is not limited.

  * `"TP #`" [list] defines a particular time period. The numbering
    should be sequential starting with 0.

    * `"PK tree`" [list] describes a hierarchical structure of the process kernels
      that reflect their weak and strong coupling.

      * `"PKNAME`"  [list] name of PK which is used in the
        simulation. Name can be arbitrary but the sublist with the same name
        should exist in the list of PKs (see below).

      * `"PK type`" [string] specifies the type of PK supported by Amanzi. At the moment
        available options are (`"darcy`", `"richards`", `"transport`", `"one-phase energy`", 
        `"two-phase energy`", `"reactive transport`", `"flow reactive transport`", 
        `"thermal richards`", `"chemistry`", `"transport implicit`", `"transport matrix fracture`",
        `"transport matrix fracture implicit`", `"flow`", and `"darcy matrix fracture`").
 
      * `"start period time`" [double] is the start time of the current time period.

      * `"end period time`" [double] is the end time of the current time period.

      * `"maximum cycle number`" [int] is the maximum allowed number of cycles in 
        the current time period. Special value -1 means unlimited number of cycles.

      * `"initial time step`" is the initial time step for the current time period.

.. code-block:: xml

  <ParameterList>  <!-- parent list -->
  <ParameterList name="cycle driver">
    <Parameter name="component names" type="Array(string)" value="{H+, Na+, NO3-, Zn++}"/>
    <Parameter name="component molar masses" type="Array(double)" value="{1.0e-3, 23.0e-3, 62.0e-3, 65.4e-3}"/>
    <Parameter name="number of liquid components" type="int" value="4"/>
    <ParameterList name="time periods">
      <ParameterList name="TP 0">
        <ParameterList name="PK tree">
          <ParameterList name="_FLOW and REACTIVE TRANSPORT">
            <Parameter name="PK type" type="string" value="flow reactive transport"/>
            <ParameterList name="_REACTIVE TRANSPORT">
              <Parameter name="PK type" type="string" value="reactive transport"/>
              <ParameterList name="_TRANSPORT">
                <Parameter name="PK type" type="string" value="transport"/>
              </ParameterList>
              <ParameterList name="_CHEMISTRY">
                <Parameter name="PK type" type="string" value="chemistry"/>
              </ParameterList>
            </ParameterList>
            <ParameterList name="_FLOW">
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

In this simulation, we use the PK labeled as *flow reactive transport*. It is
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

* `"initial time step`" [Array(double)] is the size of the first time step after we hit a special
  time specified above.

* `"maximum time step`" [Array(double)] allows the user to limit the time step between
  two particular times.

.. code-block:: xml

  <ParameterList name="cycle driver">  <!-- parent list -->
  <ParameterList name="time period control">
    <Parameter name="start times" type="Array(double)" value="{3.16e+10, 6.32e+10}"/>
    <Parameter name="initial time step" type="Array(double)" value="{100.0, 100.0}"/>
    <Parameter name="maximum time step" type="Array(double)" value="{3.2e+8, 4e+17}"/>
   </ParameterList>
   </ParameterList>

Between approximately 1000 and 2000 years, we limit the maximum time step to 10 years. 


Restart from checkpoint data file
---------------------------------

A user may request to restart a simulation from a checkpoint data file by creating list 
*restart*. In this scenario, the code will overwrite data initialized using the input XML file.
The purpose of restart is to continue the simulation that has been terminated before for some reasons,
e.g. because its allocation of time ran out.
The value for the current time and current cycle is read from the checkpoint file.

* `"restart`" [list]

  * `"file name`" [string] provides name of the existing checkpoint data file to restart from.

.. code-block:: xml
  
  <ParameterList name="cycle driver">  <!-- parent list -->
  <ParameterList name="restart">
    <Parameter name="file name" type="string" value="_CHECK00123.h5"/>
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
  field using the provided checkpoint file. Second, regardless of the outcome of the
  previous step, we try to initialize the field using the sublist `"initial conditions`".
  By design, the second step allows us to overwrite only part for the
  field. There are several options available to initialize field using
  the sublist `"initial conditions`": `"restart file`" - read field
  from existing HDF55 file, `"exodus file initialization`" - read field
  from existing Exodus file, `"cells from file`" - read cell
  components from HDF5 file, `"constant`" - set field values to constant, `"initialize
  from 1D column`" - initialize 1D column from file and `"function`" -
  field is initialized by function.


.. code-block:: xml

  <ParameterList>  <!-- parent list -->
  <ParameterList name="state">
    <Parameter name="initialization filename" type="string" value="_CHECK00123.h5"/>
    <ParameterList name="field evaluators">
       ... list of field evaluators
    </ParameterList>
    <ParameterList name="initial conditions">
      ... initialization of fields
    </ParameterList>
  </ParameterList>
  </ParameterList>


Primary and derived fields
--------------------------

* Primary fields [default units]

  * pressure [Pa]
  * total component concentration [mol/L] or [mol/m^3]
  * temperature [K]

* Secondary fields

  * saturation [-]
  * hydraulic_head [m]
  * darcy_flux (more precisely, volumetric flow rate) [m^3/s] 
  * permeability [m^2]
  * porosity [-]
  * specific_storage [m^-1]
  * specific_yield [-]
  * transport_porosity [-] 


Field evaluators
----------------

There are three different types of field evaluators.

Independent variable field evaluator
....................................

An independent ivariable field evaluator has no dependencies and is specified by a function.
Typically, it is evaluated once per simulation.
The evaluator has the following control parameters.

* `"field evaluator type`" [string] The value of this parameter is used by the factory
  of evaluators. The available option are `"independent variable`", `"primary variable`",
  and `"secondary variable`".

* `"constant in time`" [bool] specifies time-dependence nature of the field.

* `"function`" [list] defines a piecewise continuous function for calculating the independent variable.
  In may contain multiple sublists `"_DOMAIN`" with identical structure.
  
  * `"_DOMAIN`" [list]

    * `"region`" [string] specifies domain of the function, a single region.

    * `"regions`" [Array(string)] is the alternative to option `"region`", domain on 
      the function consists of many regions.

    * `"component`" [string] specifies geometric object associated with the mesh function.
      Available options are `"cell`", `"face`", and `"node`".

    * `"function`" [list] defines an analytic function for calculation. Its structure
      is described in the Functions_ section below.

    * `"initialize faces from cells`" [bool] instructs state to initialize face-component
      and boundary face-component (if any) of a composite vector from a cell-component 
      using simple averaging. Default is false.

.. code-block:: xml

  <ParameterList name="field_evaluators">  <!-- parent list -->
  <ParameterList name="saturation_liquid">
    <Parameter name="field evaluator type" type="string" value="independent variable"/>
    <Parameter name="constant in time" type="bool" value="true"/>
    <ParameterList name="function">
      <ParameterList name="_DOMAIN">
        <Parameter name="region" type="string" value="_ALL DOMAIN"/>
        <Parameter name="component" type="string" value="cell"/>
        <ParameterList name="function">
          <ParameterList name="function-constant">
            <Parameter name="value" type="double" value="0.8"/>
          </ParameterList>
        </ParameterList>
      </ParameterList>
    </ParameterList>
    <ParameterList name="verbose object">
      <Parameter name="verbosity level" type="string" value="extreme"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>

In this example field *saturation_liquid* is defined as a cell-based variable with constant value 0.8. 
Note that the user-defined name for this field cannot have spaces.


Independent variable field evaluator from file
..............................................

An independent variable field evaluator from file has no dependencies and is specified by 
data at specific time moments.

* `"filename`" [string] defines name of a data file.
  
* `"domain name`" [string] specifies mesh. Default is `"domain`".

* `"variable name`" [string] defines variable name in the data file.

* `"component name`" [string] defines component name in a composite vector.

* `"mesh entity`" [string] specifies geometric object associated with the mesh function.
  Available options are `"cell`", `"face`", and `"node`".

* `"number of dofs`" [string] defines the number of degrees of freedom. Default is 1.

* `"time function`" [list] defines a time function to interpolate data. This is the 
  optional parameter.

.. code-block:: xml

  <ParameterList name="field_evaluators">  <!-- parent list -->
  <ParameterList name="porosity">
    <Parameter name="field evaluator type" type="string" value="independent variable from file"/>
    <Parameter name="filename" type="string" value="_DATA_FILE.h5"/>
    <Parameter name="domain name" type="string" value="domain"/>
    <Parameter name="variable name" type="string" value="porosity"/>
    <Parameter name="component name" type="string" value="cell"/>
    <Parameter name="mesh entity" type="string" value="cell"/>
    <Parameter name="number of dofs" type="int" value="1"/>

    <ParameterList name="time function">  
      <Parameter name="times" type="Array(double)" value="{1.0, 2.0, 3.0}"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>

The field *porosity* is defined as a cell-based variable and
interpolated between three time intervals.


Constant variable field evaluator
.................................

Constant variable field evaluator as a simplified version of independent field evaluator from
file which allows one to define constant in time field. Initialization of the field 
has to be done in the initial conditions sublist of state.

.. code-block:: xml

  <ParameterList name="initial conditions">  <!-- parent list -->
  <ParameterList name="porosity"> 
    <ParameterList name="function">
      <ParameterList name="_ANY NAME">
        <Parameter name="regions" type="Array(string)" value="_ALL DOMAIN"/>
        <Parameter name="component" type="string" value="cell"/>
        <ParameterList name="function">
          <ParameterList name="function-constant">
            <Parameter name="value" type="double" value="90000.0"/>
          </ParameterList>
        </ParameterList>
      </ParameterList>
    </ParameterList>
  </ParameterList>
  </ParameterList>

  <ParameterList name="field_evaluators">  <!-- parent list -->
  <ParameterList name="porosity">
    <Parameter name="field evaluator type" type="string" value="constant variable"/>
  </ParameterList>
  </ParameterList>


Primary variable field evaluator
................................

The primary variable field evaluator has no dependencies solved for by a PK.
Examples of independent field evaluators are primary variable of PDEs, such as
pressure and temperature.
Typically this evaluator is used internally to inform the dependency tree about 
a new state of the primary variable.


Secondary variable field evaluators
...................................

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
    <ParameterList name="verbose object">
      <Parameter name="verbosity level" type="string" value="extreme"/>
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
      <ParameterList name="SOIL1">
        <Parameter name="regions" type="Array(string)" value="{TopRegion}"/>
        <ParameterList name="IEM parameters">
          <Parameter name="iem type" type="string" value="linear"/>
          <Parameter name="heat capacity [J/kg-K]" type="double" value="620.0"/>
        </ParameterList>
      </ParameterList>
      <ParameterList name="SOIL2">
        ...
      </ParameterList>
    </ParameterList>
    <ParameterList name="verbose object">
      <Parameter name="verbosity level" type="string" value="extreme"/>
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
parameter *value*. In three dimensions, it looks like

.. code-block:: xml

  <ParameterList name="initial conditions">  <!-- parent list -->
  <ParameterList name="gravity">
    <Parameter name="value" type="Array(double)" value="{0.0, 0.0, -9.81}"/>
  </ParameterList>
  </ParameterList>


A scalar field
..............

A variable scalar field is defined by a few functions (labeled with _MESH BLOCK #
in our example) with non-overlapping domains. 
The required parameters for each function are *region*, *component*,
and *function*.

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
      <ParameterList name="_MESH BLOCK 1">
        <Parameter name="regions" type="Array(string)" value="_DOMAIN 1"/>
        <Parameter name="component" type="string" value="cell"/>
        <ParameterList name="function">
          <ParameterList name="function-constant">
            <Parameter name="value" type="double" value="90000.0"/>
          </ParameterList>
        </ParameterList>
      </ParameterList>
      <ParameterList name="_MESH BLOCK 2">
        <Parameter name="regions" type="Array(string)" value="_DOMAIN 2, _DOMAIN 3"/>
        ... 
      </ParameterList>
    </ParameterList>
  </ParameterList>
  </ParameterList>

In this example, the field *pressure* has constant value 90000 [Pa] in 
each mesh cell of the first region. The second mesh block will define
the pressure in the second mesh region and so on.
Note that the names of functions may coinsider with the names of regions.


A vector or tensor field
........................

A variable tensor (or vector) field is defined similarly to a variable scalar field. 
The difference lies in the definition of the function which is now a multi-valued function.

* `"number of dofs`" [int] is the number of components in the vector or tensor.

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
      <ParameterList name="_MESH BLOCK 1">
        <Parameter name="regions" type="Array(string)" value="{_ALL DOMAIN}"/>
        <Parameter name="component" type="string" value="face"/>
        <ParameterList name="function">
          <Parameter name="number of dofs" type="int" value="2"/>
          <Parameter name="function type" type="string" value="composite function"/>
          <ParameterList name="dof 1 function">
            <ParameterList name="function-constant">
              <Parameter name="value" type="double" value="0.002"/>
            </ParameterList>
          </ParameterList>
          <ParameterList name="dof 2 function">
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
      <Parameter name="regions" type="Array(string)" value="{_ENTIRE DOMAIN}"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>


Initialization from Exodus II file
-------------------------------------

Some fields can be initialized from Exodus II files. 
For each field, an additional sublist has to be added to the
named sublist of *State* list with the file name and the name of attributes. 
For a serial run, the file extension must be *.exo*. 
For a parallel run, it must be *.par*.

* `"attributes`" [Array(string)] defines names of attributes. The number of names
  must be equal to the number of components in the field. The names can be repeated.
  Scalar fields (e.g. porosity) require one name, tensorial fields (e.g. permeability)
  require two or three names.

.. code-block:: xml

  <ParameterList name="initial conditions">  <!-- parent list -->
  <ParameterList name="permeability">
    <ParameterList name="exodus file initialization">
      <Parameter name="file" type="string" value="_MESH_AND_DATA.exo"/>
      <Parameter name="attributes" type="Array(string)" value="{permx, permx, permz}"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>


Initialization from HDF5 file
-------------------------------

Some field can be initialized from HDF5 file. The field has to written
to HDF5 file as 2D array (number_elements, number_of_components) and
has to name as field_name.entity.component, e.g
transport_porosity.cell.0. Parameter `"cell from file`" initializes
only cell part of the field.

.. code-block:: xml

  <ParameterList name="initial conditions">  <!-- parent list -->
  <ParameterList name="transport_porosity">
    <Parameter name="restart file" type="string" value="_TEST1.h5"/>
    </ParameterList>
  <ParameterList name="porosity">
    <Parameter name="cells from file" type="string" value="_TEST3.h5"/>
  </ParameterList>
  </ParameterList>


Initialization of 1D column
-----------------------------

It is possible to initialize only 1D column portion of a particular field.

.. code-block:: xml

  <ParameterList name="initial conditions">  <!-- parent list -->
  <ParameterList name="temperature">
    <ParameterList name="initialize from 1D column">
      <Parameter name="file" type="string" value="_COLUMN_DATA.h5"/>
      <Parameter name="z header" type="string" value="/z"/>
      <Parameter name="f header" type="string" value="/temperature"/>
      <Parameter name="coordinate orientation" type="string" value="depth"/>
      <Parameter name="surface sideset" type="string" value="surface"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>


Mesh partitioning
-----------------

Amanzi's state has a number of tools to verify completeness of initial data.
This is done using list *mesh partitions*. 
Each sublist must have parameter *region list* specifying
regions that define unique partition of the mesh.

.. code-block:: xml

  <ParameterList name="state">  <!-- parent list -->
  <ParameterList name="mesh partitions">
    <ParameterList name="_MATERIALS">
      <Parameter name="region list" type="Array(string)" value="{_SAND1,_CLAY,_SAND2}"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>

In this example, we verify that three mesh regions representing sand and clay cover completely
the mesh without overlaps.
If so, all material fields, such as *porosity* and *permeability*, will be initialized properly.


Data IO control
---------------

Two parameters below allow us to control fields that will go into a visuzalization file.
First, we remove all fields matching the patterns specified by *blacklist*.
Second, we add all fields matching the patterns specified by *whitelist*.
Both parameters are optional.

* `"blacklist`" [Array(string)] list of fields that should *not* be written to the visualization file.
  Standard regular expressuion rules can be used, e.g. `"(secondary_)(.*)`" skips all fields 
  those names start with `"secondary_`".

* `"whitelist`" [Array(string)] list of fields that should *be* written to the visualization file.
  Standard regular expressuion rules can be used, e.g. `"(primary_)(.*)`" adds all fields 
  those names start with `"primary_`".


Example
-------

The complete example of a state initialization is below. Note that
_MATERIAL1 and _MATERIAL2 must be labels of the existing regions that cover
the computational domain.
The fields *porosity* and *pressure* are constant over the whole domain. 
The field *permeability* is the piecewise constant diagonal tensor.

.. code-block:: xml

  <ParameterList name="state">
  <ParameterList name="field evaluators">
    <ParameterList name="porosity">
      <ParameterList name="function">
        <ParameterList name="_ANY NAME ">
          <Parameter name="regions" type="Array(string)" value="{_ALL DOMAIN}"/>
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
        <ParameterList name="_ANY NAME">
          <Parameter name="regions" type="Array(string)" value="{_ALL DOMAIN}"/>
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
        <ParameterList name="_ANY NAME">
          <Parameter name="regions" type="Array(string)" value="_MATERIAL1"/>
          <Parameter name="component" type="string" value="cell"/>
          <ParameterList name="function">
            <Parameter name="function type" type="string" value="composite function"/>
            <Parameter name="number of dofs" type="int" value="2"/>
            <ParameterList name="dof 1 function">
              <ParameterList name="function-constant">
                <Parameter name="value" type="double" value="1e-12"/>
              </ParameterList>
            </ParameterList>
            <ParameterList name="dof 2 function">
              <ParameterList name="function-constant">
                <Parameter name="value" type="double" value="1e-13"/>
              </ParameterList>
            </ParameterList>
          </ParameterList>
        </ParameterList>
        <ParameterList name="_ANY_NAME">
          <Parameter name="regions" type="Array(string)" value="_MATERIAL2"/>
          <Parameter name="component" type="string" value="cell"/>
          <ParameterList name="function">
            <Parameter name="function type" type="string" value="composite function"/>
            <Parameter name="number of dofs" type="int" value="2"/>
            <ParameterList name="dof 1 function">
              <ParameterList name="function-constant">
                <Parameter name="value" type="double" value="2e-13"/>
              </ParameterList>
            </ParameterList>
            <ParameterList name="dof 2 function">
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
The name of the PKs in this list must match *PKNAMEs* in *cycle driver* list.

.. code-block:: xml

  <ParameterList>  <!-- parent list -->
  <ParameterList name="PKs">
    <ParameterList name="_FLOW and TRANSPORT">
      <Parameter name="PK type" type="string" value="flow reactive transport"/>      
      <Parameter name="PKs order" type="Array(string)" value="{_FLOW, _TRANSPORT}"/> 
      <Parameter name="master PK index" type="int" value="0"/>
    </ParameterList>
    <ParameterList name="_FLOW">
      ...
      ... flow parameters, lists, and sublists
      ...
    </ParameterList>
    <ParameterList name="_TRANSPORT">
      ...
      ... transport parameters, lists, and sublists
      ...
    </ParameterList>
  </ParameterList>
  </ParameterList>


Flow PK
-------

Mathematical models
...................

A few PDE models can be instantiated using the parameters described below.


Fully saturated flow
````````````````````

The conceptual PDE model for the fully saturated flow is

.. math::
  \phi \left(\frac{S_s}{g} + \frac{S_y}{Lg}\right)\frac{\partial p_l}{\partial t} 
  =
  -\boldsymbol{\nabla} \cdot (\rho_l \boldsymbol{q}_l) + Q,
  \quad
  \boldsymbol{q}_l 
  = -\frac{\boldsymbol{K}}{\mu} 
  (\boldsymbol{\nabla} p - \rho_l \boldsymbol{g}),

where 
:math:`\phi` is porosity [-],
:math:`s_s` and :math:`s_y` are specific storage [m] and specific yield [-], respectively,
:math:`L` is characteristic length [m],
:math:`\rho_l` is fluid density [:math:`kg / m^3`],
:math:`Q` is source or sink term [:math:`kg / m^3 / s`],
:math:`\boldsymbol{q}_l` is the Darcy velocity [:math:`m/s`],
and :math:`\boldsymbol{g}` is gravity [:math:`m/s^2`].


Partially saturated flow
````````````````````````

The conceptual PDE model for the partially saturated flow is

.. math::
  \frac{\partial \theta}{\partial t} 
  =
  -\boldsymbol{\nabla} \cdot (\eta_l \boldsymbol{q}_l) + Q,
  \qquad
  \boldsymbol{q}_l 
  = -\frac{\boldsymbol{K} k_r}{\mu} 
  (\boldsymbol{\nabla} p - \rho_l \boldsymbol{g})

where 
:math:`\theta` is total water content [:math:`mol/m^3`],
:math:`\eta_l` is molar density of liquid [:math:`mol/m^3`],
:math:`\rho_l` is fluid density [:math:`kg/m^3`],
:math:`Q` is source or sink term [:math:`mol/m^3/s`],
:math:`\boldsymbol{q}_l` is the Darcy velocity [:math:`m/s`],
:math:`k_r` is relative permeability [-],
and :math:`\boldsymbol{g}` is gravity [:math:`m/s^2`].
We define 

.. math::
  \theta = \phi \eta_l s_l

where :math:`s_l` is liquid saturation [-],
and :math:`\phi` is porosity [-].


Partially saturated flow with water vapor
`````````````````````````````````````````

The conceptual PDE model for the partially saturated flow with water vapor 
includes liquid phase (liquid water) and gas phase (water vapor):

.. math::
  \frac{\partial \theta}{\partial t} 
  =
  - \boldsymbol{\nabla} \cdot (\eta_l \boldsymbol{q}_l)
  - \boldsymbol{\nabla} \cdot (\boldsymbol{K}_g \boldsymbol{\nabla} \big(\frac{p_v}{p_g}\big)) + Q,
  \quad
  \boldsymbol{q}_l 
  = -\frac{\boldsymbol{K} k_r}{\mu} 
  (\boldsymbol{\nabla} p - \rho_l \boldsymbol{g})

where 
:math:`\theta` is total water content [:math:`mol/m^3`],
:math:`\eta_l` is molar density of liquid (water) [:math:`mol/m^3`],
:math:`\rho_l` is fluid density [:math:`kg/m^3`],
:math:`Q` is source or sink term [:math:`mol/m^3/s`],
:math:`\boldsymbol{q}_l` is the Darcy velocity [:math:`m/s`],
:math:`k_r` is relative permeability [-],
:math:`\boldsymbol{g}` is gravity [:math:`m/s^2`],
:math:`p_v` is the vapor pressure [Pa],
:math:`p_g` is the gas pressure [Pa],
and :math:`\boldsymbol{K}_g` is the effective diffusion coefficient of the water vapor.
We define 

.. math::
  \theta = \phi \eta_l s_l + \phi \eta_g (1 - s_l) X_g

where :math:`s_l` is liquid saturation [-],
:math:`\phi` is porosity [-],
:math:`\eta_g` is molar density of water vapor [:math:`mol/m^3`],
and :math:`X_g` is molar fraction of water vapor.
The effective diffusion coefficient of the water vapor is given by

.. math::
  \boldsymbol{K}_g = \phi s_g \tau_g \eta_g \boldsymbol{D}_g

where :math:`s_g` is gas saturation [-],
:math:`\tau_g` is the tortuosity of the gas phase [-],
:math:`\eta_g` is the molar density of gas [:math:`kg/m^3`],
and :math:`\boldsymbol{D}_g` is the diffusion coefficient of the gas phase [:math:`m^2/s`],
The gas pressure :math:`p_g` is set to the atmosperic pressure and the vapor pressure
model assumes thermal equlibrium of liquid and gas phases:

.. math::
  p_v = P_{sat}(T) \exp\left(\frac{P_{cgl}}{\eta_l R T}\right)

where
:math:`R` is the ideal gas constant [:math:`kg m^2/K/mol/s^2`],
:math:`P_{cgl}` is the liquid-gas capillary pressure [Pa],
:math:`P_{sat}` is the saturated vapor pressure [Pa],
and :math:`T` is the temperature [K].
The diffusion coefficient is based of TOUGHT2 model

.. math::
   D_g = D_0 \frac{P_{ref}}{p} \left(\frac{T}{273.15}\right)^a

where
:math:`D_0 = 2.14 \cdot 10^{-5}`,
:math:`P_{ref}` is atmospheric pressure,
and :math:`a = 1.8`. 
finally we need a model for the gas tortuosity. We use the Millington and Quirk model:

.. math::
   \tau_g = \phi^\beta s_g^\gamma

where
:math:`\beta = 1/3` and 
:math:`\gamma = 7/3`.


Isothermal partially saturated flow with dual porosity model
````````````````````````````````````````````````````````````

The conceptual model for the partially saturated flow with dual porosity model
assumes that water flow is restricted to the fractures and the water in the matrix does not move.
The rock matrix represents immobile pockets that can exchange, retain and store water
but do not permit convective flow.
This leads to dual-porosity type flow and transport models that partition the liquid
phase into mobile and immobile regions.
The Richards equation in the mobile region is augmented by the water exchange
term :math:`\Sigma_w`:
 
.. math::
  \frac{\partial \theta_{lf}}{\partial t} 
  = -\boldsymbol{\nabla} \cdot (\eta_l \boldsymbol{q}_l) 
    -\frac{K_m\,k_{rm}\,\eta_l}{\mu\,L_m}\, \nabla p_m + Q_f,
  \qquad
  \boldsymbol{q}_l 
  = -\frac{\boldsymbol{K}_f\, k_{rf}}{\mu} 
  (\boldsymbol{\nabla} p_f - \rho_l \boldsymbol{g})

where 
:math:`p_f` is fracture pressure [Pa],
:math:`p_m` is matrix pressure [Pa],
:math:`L_m` is the characteristic matrix depth defined typically as the ratio of a matrix block [m],
and :math:`Q_f` is source or sink term [:math:`kg \cdot m^{-3} \cdot s^{-1}`].
The equation for water balance in the matrix is

.. math::
  \frac{\partial \theta_{lm}}{\partial t} 
  = Q_m
    +\nabla\cdot \left(\frac{K_m\, k_{rm}\,\eta_l}{\mu}\, \nabla p_{m}\right),

where 
:math:`Q_m` is source or sink term [:math:`kg / m^3 / s`].
The volumetric volumetric water contents are defined as

.. math::
  \theta_f = \phi_f\, \eta_l\, s_{lf},\quad
  \theta_m = \phi_m\, \eta_l\, s_{lm},

where saturations :math:`s_{lf}` and :math:`s_{lm}` may use different capillary 
pressure - saturation models.
In the simplified model, the rate of water exchange between the fracture and matrix regions 
is proportional to the difference in hydraulic heads:

.. math::
  \frac{K_m\,k_{rm}\,\eta_l}{\mu\,L_m}\, \nabla p_m 
  \approx
  \alpha_w (h_f - h_m),

where :math:`\alpha_w` is the mass transfer coefficient.
Since hydraulic heads are needed for both regions, this equation requires
retention curves for both regions and therefore is nonlinear.
 

Model specifications and assumptions
....................................

This list is used to summarize physical models and assumptions, such as
coupling with other PKs.
This list is often generated or extended by a high-level MPC PK.
In the code development, this list plays a two-fold role. 
First, it provides necessary information for coupling different PKs such 
as flags for adding a vapor diffusion to Richards' equation.
Second, the developers may use it instead of a factory of evaluators such as
creation of primary and secondary evaluators for rock porosity models.
Combination of both approaches may lead to a more efficient code.

* `"vapor diffusion`" [bool] is set up automatically by a high-level PK,
  e.g. by EnergyFlow PK. The default value is `"false`".

* `"flow in fractures`" [bool] indicates that Darcy flow is calculated in fractures. 
  This option is ignored is mesh dimentionaly equals to manifold dimensionality.

* `"multiscale model`" [string] specifies a multiscale model.
  Available options are `"single porosity`" (default) and `"dual continuum discontinum matrix`".

* `"water content model`" [string] changes the evaluator for water
  content. Available options are `"generic`" and `"constant density`" (default).

* `"viscosity model`" [string] changes the evaluator for liquid viscosity.
  Available options are `"generic`" and `"constant viscosity`" (default).

* `"porosity model`" [string] specifies an isothermal porosity model.
  Available options are `"compressible: storativity coefficient`",
  `"compressible: pressure function`", and `"constant porosity`" (default).

* `"coupled matrix fracture flow`" [string] specifies PK's role in the strong 
  coupling of two flow PKs. The value is either `"matrix`" or `"fracture`".

.. code-block:: xml

  <ParameterList name="flow">  <!-- parent list -->
  <ParameterList name="physical models and assumptions">
    <Parameter name="vapor diffusion" type="bool" value="false"/>
    <Parameter name="water content model" type="string" value="constant density"/>
    <Parameter name="viscosity model" type="string" value="constant viscosity"/>
    <Parameter name="porosity model" type="string" value="compressible: pressure function"/>
    <Parameter name="multiscale model" type="string" value="single porosity"/>
    <Parameter name="coupled matrix fracture flow" type="string" value="matrix"/>
  </ParameterList>
  </ParameterList>


Global parameters
.................

* `"domain name`" [string] specifies mesh name that defined domain of this PK.
  Default is `"domain`".


Water retention models
......................

User defines water retention models in sublist *water retention models*. 
It contains as many sublists, e.g. *SOIL_1*, *SOIL_2*, etc, as there are different soils. 
This list is required for the Richards problem only.
 
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

  <ParameterList name="flow">  <!-- parent list -->
  <ParameterList name="water retention models">
    <ParameterList name="_SOIL_1">
      <Parameter name="regions" type="Array(string)" value="{_TOP HALF}"/>
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

    <ParameterList name="_SOIL_2">
      <Parameter name="regions" type="Array(string)" value="{_BOTTOM HALF}"/>
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
It contains as many sublists, e.g. _SOIL1 and _SOIL2, as there are different soils. 
The porosity models are associated with non-overlapping regions. Each of the sublists (e.g. _SOIL1) 
includes a few mandatory parameters: *regions names*, *model name*, and parameters for the selected model.

* `"porosity model`" [string] specifies a model for the soil.
  The available models are `"compressible`" and `"constant`". 

  * The model `"compressible`" requires `"underformed soil porosity"`" [double],
    `"reference pressure`" [double], and `"pore compressibility`" [string] [Pa^-1].
    Default value for `"reference pressure`" is 101325.0 [Pa].

  * The model `"constant`" requires `"value`" [double].

.. code-block:: xml

  <ParameterList name="flow">  <!-- parent list -->
  <ParameterList name="porosity models">
    <ParameterList name="_SOIL1">
      <Parameter name="regions" type="Array(string)" value="{_TOP HALF}"/>
      <Parameter name="porosity model" type="string" value="constant"/>
      <Parameter name="value" type="double" value="0.2"/>
    </ParameterList>

    <ParameterList name="_SOIL2">
      <Parameter name="regions" type="Array(string)" value="{_BOTTOM HALF}"/>
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

The list *multiscale models* is the place for various multiscale models.
The list is extension of the list *water retention models*. 
Its ordered by soil regions and includes parameters for the multiscale,
capillary pressure, and relative permebility models.
This list is optional. 

* `"multiscale model`" [string] is the model name. Available options are `"dual porosity`"
  and `"generalized dual porosity`".

* `"xxx parameters`" [sublist] provides parameters for the model specified by variable `"multiscale model`".

* `"water retention model`" [string] specifies a model for the soil.
  The available models are `"van Genuchten`" and `"Brooks Corey`". 
  Parameters for each model are described above.

* `"relative permeability model`" [string] The available options are `"Mualem`" (default) 
  and `"Burdine`".


Dual porosity model
```````````````````

* `"mass transfer coefficient`" [double] is the mass transfer coefficient.

* `"tolerance`" [double] defines tolerance for iterative methods used to solve
  secondary equations. Default is 1e-8.


Generalized dual porosity model
```````````````````````````````

* `"number of matrix nodes`" [int] defines number of matrix layers.
* `"matrix depth`" [double] is the characteristic length for matrix continuum.
* `"matrix volume fraction`" [double] defines relative volume of matrix continuum.

.. code-block:: xml

  <ParameterList name="flow">  <!-- parent list -->
  <ParameterList name="multiscale models"> 
    <ParameterList name="_SOIL1">
      <Parameter name="regions" type="Array(string)" value="{_TOP HALF}"/>
      <Parameter name="multiscale model" type="string" value="dual porosity"/> 
      <ParameterList name="dual porosity parameters">
        <Paramater name="mass transfer coefficient" type="double" value="4.0e-5"/>
        <Paramater name="tolerance" type="double" value="1e-8"/>
      </ParameterList>

      <Parameter name="water retention model" type="string" value="van Genuchten"/>
      <Parameter name="van Genuchten alpha" type="double" value="0.000194"/>
      <Parameter name="van Genuchten m" type="double" value="0.28571"/>
      <Parameter name="van Genuchten l" type="double" value="0.5"/>
      <Parameter name="residual saturation" type="double" value="0.103"/>
      <Parameter name="relative permeability model" type="string" value="Mualem"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>


Absolute permeability
.....................

* `"coordinate system`" [string] defines coordinate system
  for calculating absolute permeability. The available options are `"cartesian`"
  and `"layer`".

* `"off-diagonal components`" [int] defines additional (typically off-diagonal) 
  components of the absolute permeability. Deafult is 0.

.. code-block:: xml

  <ParameterList name="flow">  <!-- parent list -->
  <ParameterList name="absolute permeability">
    <Parameter name="coordinate system" type="string" value="cartesian"/>
    <Parameter name="off-diagonal components" type="int" value="0"/>
  </ParameterList>
  </ParameterList>


Relative permeability
.....................

This section discusses interface treatment of cell-centered fields such as 
relative permeability, density and viscosity.

* `"relative permeability`" [list] collects information required for treatment of
  relative permeability, density and viscosity on mesh faces.

  * `"upwind method`" [string] defines a method for calculating the *upwinded* 
    relative permeability. The available options are: `"upwind: gravity`", 
    `"upwind: darcy velocity`" (default), `"upwind: second-order`", `"upwind: amanzi`" (experimental), 
    `"upwind: amanzi new`" (experiemental), `"other: harmonic average`", and `"other: arithmetic average`".

  * `"upwind frequency`" [string] defines frequency of recalculating Darcy flux inside
    nonlinear solver. The available options are `"every timestep`" and `"every nonlinear iteration`".
    The first option freezes the Darcy flux for the whole time step. The second option
    updates it on each iteration of a nonlinear solver. The second option is recommended
    for the Newton solver. It may impact significantly upwinding of the relative permeability 
    and convergence rate of this solver.

  * `"upwind parameters`" [list] defines parameters for upwind method specified by `"relative permeability`".

    * `"tolerance`" [double] specifies relative tolerance for almost zero local flux. In such
      a case the flow is assumed to be parallel to a mesh face. Default value is 1e-12.

    * `"method`" [string] specifies a reconstruction method. Available option is
      `"cell-based`" (default).

    * `"polynomial order`" [int] defines the polynomial order of a reconstructed function. Default is 1.

    * `"limiter`" [string] specifies limiting method for a high-order reconstruction. 
      Available options are `"Barth-Jespersen`" (default), `"Michalak-Gooch`", `"tensorial`",
      and `"Kuzmin`". 

.. code-block:: xml

  <ParameterList name="flow">  <!-- parent list -->
  <ParameterList name="relative permeability">
    <Parameter name="upwind method" type="string" value="upwind: darcy velocity"/>
    <Parameter name="upwind frequency" type="string" value="every timestep"/>
    <ParameterList name="upwind parameters">
       <Parameter name="tolerance" type="double" value="1e-12"/>
       <Parameter name="method" type="string" value="cell-based"/>
       <Parameter name="polynomial order" type="int" value="1"/>
       <Parameter name="limiter" type="string" value="Barth-Jespersen"/>
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
      When the Richards problem is set up, Flow PK sets up proper value for parameter `"upwind method`" of 
      this sublist.

    * `"preconditioner`" [list] defines parameters for generating and assembling diffusion 
      matrix that is used to create preconditioner. 
      This sublist is ignored for the saturated problem.
      Since update of preconditioner can be lagged, we need two objects called `"matrix`" and `"preconditioner`".
      When the Richards problem is set up, Flow PK sets up proper value for parameter `"upwind method`" of 
      this sublist.

.. code-block:: xml

  <ParameterList name="flow">  <!-- parent list -->
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
        <Parameter name="preconditioner schema" type="Array(string)" value="{face, cell}"/>
        <Parameter name="Newton correction" type="string" value="approximate Jacobian"/>
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

  * `"relative to bottom`" [bool] indicates the submodel where the static head is defined with respect
    to the bottom boundary (a curve in 3D) of the specified regions. Support of 2D is turned off.
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
    to mixed boundary condition. The recommended value is `"PFloTran`".

  * `"seepage flux threshold`" [double] sets up the threshold for switching from the pressure 
    to influx boundary condition in submodel `"PFloTran`". The pressure condition remains 
    for a small influx value until it exceeds the certain fraction of the `"mass flux`" specified 
    by this parameter. The admissible range is from 0 to 0.1. Default value is 0. 

Each boundary condition accepts three parameters: `"regions`", 
`"use area fractions`", and `"spatial distribution method`". Parameter `"regions`"
specifies the list of regions where the boundary condition is defined. 
The boolen parameter `"use area fractions`" instructs the code to use all available volume fractions. 
Default value is *false*, it corresponds to :math:`f=1` in the formulas below.
Parameter `"spatial distribution method`" defines the method for distributing
data (e.g. the total mass flux) over the specified regions. The available options 
are `"volume`", `"permeability`", `"domain coupling`", `"subgrid`", `"simple well`", or `"none`". 
For instance, for a given boundary function :math:`g(x)`, these options correspond to 
different boundary conditions for the Darcy velocity in the original PDE:

.. math::
  {\boldsymbol q} \cdot {\boldsymbol n} = g(x)\, f\, \frac{1}{|B|},\quad\mbox{and}\quad
  {\boldsymbol q} \cdot {\boldsymbol n} = g(x)\, f,

where :math:`f` is the folume fraction function, and :math:`|B|` is the area of the
specified regions calculated using the folume fraction function.

.. code-block:: xml

  <ParameterList name="flow">  <!-- parent list -->
  <ParameterList name="boundary conditions">
    <ParameterList name="pressure">
      <ParameterList name="_BC 0">
        <Parameter name="regions" type="Array(string)" value="{_WEST_SIDE}"/>
        <Parameter name="spatial distribution method" type="string" value="none"/>
        <ParameterList name="boundary pressure">
          <ParameterList name="function-constant">
            <Parameter name="value" type="double" value="101325.0"/>
          </ParameterList>
        </ParameterList>
      </ParameterList>
    </ParameterList>

    <ParameterList name="mass flux">
      <ParameterList name="_BC 1">
        <Parameter name="regions" type="Array(string)" value="{_NORTH_SIDE, _SOUTH_SIDE}"/>
        <Parameter name="spatial distribution method" type="string" value="volume"/>
        <Parameter name="rainfall" type="bool" value="false"/>
        <ParameterList name="outward mass flux">
          <ParameterList name="function-constant">
            <Parameter name="value" type="double" value="0.0"/>
          </ParameterList>
        </ParameterList>
      </ParameterList>
    </ParameterList>

    <ParameterList name="static head">
      <ParameterList name="_BC 2">
        <Parameter name="regions" type="Array(string)" value="{_EAST_SIDE_TOP}"/>
        <Parameter name="spatial distribution method" type="string" value="none"/>
        <Parameter name="relative to top" type="bool" value="true"/>
        <Parameter name="relative to bottom" type="bool" value="true"/>
        <ParameterList name="static head">
          <ParameterList name="function-static-head">
            <Parameter name="p0" type="double" value="101325.0"/>
            <Parameter name="density" type="double" value="1000.0"/>
            <Parameter name="gravity" type="double" value="9.8"/>
            <Parameter name="space dimension" type="int" value="3"/>
            <ParameterList name="water table elevation">
              <ParameterList name="function-constant">
                <Parameter name="value" type="double" value="10.0"/>
              </ParameterList>
            </ParameterList>
          </ParameterList>
        </ParameterList>
      </ParameterList>
    </ParameterList>

    <ParameterList name="seepage face">
      <ParameterList name="_BC 3">
        <Parameter name="regions" type="Array(string)" value="{_EAST_SIDE_BOTTOM}"/>
        <Parameter name="spatial distribution method" type="string" value="none"/>
        <Parameter name="rainfall" type="bool" value="true"/>
        <Parameter name="submodel" type="string" value="PFloTran"/>
        <Parameter name="reference pressure" type="double" value="101325.0"/>
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

The sources and sinks are typically associated with wells. 
Negative source means a producing well. 
Positive source means an injecting well. 
The structure of list *source terms* mimics that of list *boundary conditions*. 
Again, constant functions can be replaced by any of the available functions.

* `"regions`" [Array(string)] is the list of regions where the source is defined.

* `"spatial distribution method`" [string] is the method for distributing
  source Q over the specified regions. The available options are `"volume`",
  `"none`", `"permeability`" and `"simple well`".
  For option `"none`", the source term function Q is measured in [kg/m^3/s]. 
  For the other options, it is measured in [kg/s]. 
  When the source function is defined over a few regions, Q is distributed over their union.
  Option `"volume fraction`" can be used when the region geometric
  model support volume fractions. Option `"simple well`" implements the Peaceman model. 
  The well flux is defined as `q_w = WI (p - p_w)` [kg/s], where `WI` is the well index 
  and `p_w` is the well pressure. The pressure in a well is assumed to be hydrostatic.

* `"use volume fractions`" instructs the code to use all available volume fractions. 
  Note that the region geometric model supports volume fractions only for a few regions.

* `"submodel`" [string] refines definition of the source. Available options are `"rate`",
  `"integrated source`" and `"bhp"` (bottom hole pressure). The first option defines the source 
  in a natural way as the rate of change `q`. The second option defines the indefinite
  integral `Q` of the rate of change, i.e. the source term is calculated as `q = dQ/dt`. 
  For most distributions methods, two submodles are available: `"rate`" and `"integrated source`".
  For distribution method `"simple well`", two submodels are available: `"rate`" and
  `"bhp`". Submodel `"bhp`" requires `"depth"`, `"well radius`" and 
  `"bhp`" function. Submodel `"rate`" requires only rate function.
  Default is `"rate`". 

.. code-block:: xml

  <ParameterList name="flow">  <!-- parent list -->
  <ParameterList name="source terms">
    <ParameterList name="_SRC 0">
      <Parameter name="regions" type="Array(string)" value="{_WELL_EAST}"/>
      <Parameter name="spatial distribution method" type="string" value="volume"/>
      <Parameter name="submodel" type="string" value="rate"/>
      <ParameterList name="well">
        <ParameterList name="function-constant">
          <Parameter name="value" type="double" value="-0.1"/>
        </ParameterList>
      </ParameterList>
    </ParameterList>

    <ParameterList name="_SRC 1">
      <Parameter name="regions" type="Array(string)" value="{_WELL_WEST}"/>
      <Parameter name="spatial distribution method" type="string" value="permeability"/>
      <ParameterList name="well">
        <ParameterList name="function-constant">
          <Parameter name="value" type="double" value="-0.2"/>
        </ParameterList>
      </ParameterList>
    </ParameterList>

    <ParameterList name="_SRC 2">
      <Parameter name="regions" type="Array(string)" value="{_WELL_NORTH}"/>
        <Parameter name="spatial distribution method" type="string" value="simple well"/>    
        <ParameterList name="well">
          <Parameter name="submodel" type="string" value="bhp"/>
          <Parameter name="depth" type="double" value="-2.5"/>
          <Parameter name="well radius" type="double" value="0.1"/>
          <ParameterList name="bhp">
            <ParameterList name="function-constant">
              <Parameter name="value" type="double" value="10.0"/>
            </ParameterList>
          </ParameterList>
        </ParameterList>
      </ParameterList>
    </ParameterList>

    <ParameterList name="_SRC 3">
      <Parameter name="regions" type="Array(string)" value="{_WELL_SOUTH}"/>
        <Parameter name="spatial distribution method" type="string" value="simple well"/>
        <ParameterList name="well">
          <Parameter name="submodel" type="string" value="rate"/>
          <ParameterList name="rate">
            <ParameterList name="function-constant">
              <Parameter name="value" type="double" value="100.0"/>
            </ParameterList>
          </ParameterList>
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

* `"linear solver`" [string] refers to a generic linear solver from list `"solvers`".
  It is used in all cases except for `"initialization`" and `"enforce pressure-lambda constraints`".
  Currently, it is used by the Darcy PK only.

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

  * `"active wells`" [bool] specifies if wells are active or turned off. Default is *false*.

  * `"picard parameters`" [list] defines control parameters for the Picard solver.

    * `"convergence tolerance`" [double] specifies nonlinear convergence tolerance. 
      Default is 1e-8.
    * `"maximum number of iterations`" [int] limits the number of iterations. Default is 400. 

  * `"linear solver`" [string] refers to a solver sublist of the list `"solvers`". 

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

  * `"linear solver`" [string] refers to a solver sublist of the list `"solvers`".

  * `"inflow krel correction`" [bool] estimates relative permeability on inflow 
    mesh faces. This estimate is more reliable than the upwinded relative permeability
    value, especially in steady-state calculations.

.. code-block:: xml

  <ParameterList name="flow">  <!-- parent list -->
  <ParameterList name="time integrator">
    <Parameter name="error control options" type="Array(string)" value="{pressure, saturation}"/>
    <Parameter name="linear solver" type="string" value="_GMRES_WITH_AMG"/>
     <Parameter name="preconditioner" type="string" value="_HYPRE_AMG"/>
     <Parameter name="preconditioner enhancement" type="string" value="none"/>

     <ParameterList name="initialization">  <!-- first method -->
       <Parameter name="method" type="string" value="saturated solver"/>
       <Parameter name="linear solver" type="string" value="_PCG_WITH_AMG"/>
       <Parameter name="clipping pressure value" type="double" value="50000.0"/>
     </ParameterList>

     <ParameterList name="initialization">  <!-- alternative method -->
       <Parameter name="method" type="string" value="picard"/>
       <Parameter name="linear solver" type="string" value="_PCG_WITH_AMG"/>
       <ParameterList name="picard parameters">
         <Parameter name="convergence tolerance" type="double" value="1e-8"/> 
         <Parameter name="maximum number of iterations" type="int" value="20"/> 
       </ParameterList>
     </ParameterList>

     <ParameterList name="pressure-lambda constraints">
       <Parameter name="method" type="string" value="projection"/>
       <Parameter name="inflow krel correction" type="bool" value="false"/>
       <Parameter name="linear solver" type="string" value="_PCG_WITH_AMG"/>
     </ParameterList>
   </ParameterList>
   </ParameterList>


Time step controller and nonlinear solver
`````````````````````````````````````````

The time step is controlled by parameter *time step controller type*
and the related list of options, see section TimeStepController_ for the list
of supported parameter.
Nonlinear solver is controlled by parameter *solver type*  and related list of options.
Amanzi supports a few nonlinear solvers described in details in a separate section.

.. code-block:: xml

  <ParameterList name="flow">  <!-- parent list -->
  <ParameterList name="time integrator">
    <Parameter name="time integration method" type="string" value="BDF1"/>
    <ParameterList name="BDF1">
      <Parameter name="max preconditioner lag iterations" type="int" value="5"/>
      <Parameter name="extrapolate initial guess" type="bool" value="true"/>
      <Parameter name="restart tolerance relaxation factor" type="double" value="1000.0"/>
      <Parameter name="restart tolerance relaxation factor damping" type="double" value="0.9"/>

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
        <ParameterList name="verbose object">
          <Parameter name="verbosity level" type="string" value="high"/>
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


Other parameters
................

The remaining *flow* parameters are

* `"clipping parameters`" [list] defines how solution increment calculated by a nonlinear 
  solver is modified e.g., clipped.

  * `"maximum saturation change`" [double] Default is 0.25.

  * `"pressure damping factor`" [double] Default is 0.5.

* `"plot time history`" [bool] produces an ASCII file with the time history. Default is `"false`".

* `"algebraic water content balance`" [bool] uses algebraic correction to enforce consistency of 
  water content and Darcy fluxes. It leads to a monotone transport. Default is *false*.

.. code-block:: xml

  <ParameterList name="flow">  <!-- parent list -->
  <ParameterList name="clipping parameters">
     <Parameter name="maximum saturation change" type="double" value="0.25"/>
     <Parameter name="pressure damping factor" type="double" value="0.5"/>
  </ParameterList>	

  <Parameter name="plot time history" type="bool" value="false"/>
  <Parameter name="algebraic water content balance" type="bool" value="false"/>
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

Mathematical models
...................

A few PDE models can be instantiated using the parameters described below.


Single-phase transport
``````````````````````

The conceptual PDE model for the transport in partially saturated media is

.. math::
  \frac{\partial (\phi s_l C_l)}{\partial t} 
  =
  - \boldsymbol{\nabla} \cdot (\boldsymbol{q}_l C_l) 
  + \boldsymbol{\nabla} \cdot (\phi_e s_l\, (\boldsymbol{D}_l + \tau \boldsymbol{M}_l) \boldsymbol{\nabla} C_l) + Q,

where 
:math:`\phi` is total porosity [-],
:math:`\phi_e` is effective transport porosity [-],
:math:`s_l` is liquid saturation [-], 
:math:`Q` is source or sink term,
:math:`\boldsymbol{q}_l` is the Darcy velocity [m/s],
:math:`\boldsymbol{D}_l` is dispersion tensor,
:math:`\boldsymbol{M}_l` is diffusion coefficient,
and :math:`\tau` is tortuosity [-].
For an isotropic medium with no preferred axis of symmetry the dispersion 
tensor has the following form:

.. math::
  \boldsymbol{D}_l 
  = \alpha_t \|\boldsymbol{v}\| \boldsymbol{I} 
  + \left(\alpha_l-\alpha_t \right) 
    \frac{\boldsymbol{v} \boldsymbol{v}}{\|\boldsymbol{v}\|}, \qquad
  \boldsymbol{v} = \frac{\boldsymbol{q}}{\phi_e}

where
:math:`\alpha_l` is longitudinal dispersivity [m],
:math:`\alpha_t` is  transverse dispersivity [m],
and :math:`\boldsymbol{v}` is average pore velocity [m/s].
Amanzi supports two additional models for dispersivity with 3 and 4 parameters.


Single-phase transport with dual porosity model
```````````````````````````````````````````````

The dual porosity formulation of the solute transport consists of two equations
for the fracture and matrix regions. 
In the fracture region, we have \citep{simunek-vangenuchten_2008}

.. math::
  \frac{\partial (\phi_f\, s_{lf}\, C_{lf})}{\partial t} 
  =
  - \boldsymbol{\nabla} \cdot (\boldsymbol{q}_l C_{lf}) 
  + \boldsymbol{\nabla} \cdot (\phi_f\, s_{lf}\, (\boldsymbol{D}_l + \tau_f M) \boldsymbol{\nabla} C_{lf}) 
  - \frac{\phi_m\,\tau_m}{L_m}\, M \nabla C_m - \Sigma_w C^* + Q_f,

where 
:math:`\phi_f` is fracture porosity [-],
:math:`\phi_m` is matrix porosity [-],
:math:`s_{lf}` is liquid saturation in fracture [-], 
:math:`\boldsymbol{q}_l` is the Darcy velocity [m/s],
:math:`\boldsymbol{D}_l` is dispersion tensor,
:math:`\tau_f` is fracture tortuosity [-],
:math:`\tau_m` is matrix tortuosity [-],
:math:`M` is molecular diffusion coefficient [:math:`m^2/s`], and
:math:`L_m` is the characteristic matrix depth defined typically as the ratio of a matrix block [m],
:math:`\Sigma_w` is transfer rate due to flow from the matrix to the fracture, 
:math:`C^*` is equal to :math:`C_{lf}` if :math:`\Sigma_w > 0` and :math:`C_{lm}` is :math:`\Sigma_w < 0`,
and :math:`Q_f` is source or sink term.
In the matrix region, we have

.. math::
  \frac{\partial (\phi_m\, s_{lm}\, C_{lm})}{\partial t}
  = \nabla\cdot (\phi_m\, \tau_m\, M_m \nabla C_{lm}) + \Sigma_w C^* + Q_m,

where 
:math:`\phi_m` is matrix porosity [-],
:math:`s_{lm}` is liquid saturation in matrix [-], 
:math:`Q_m` is source or sink term.
The simplified one-node dual porosity model uses a finite difference approximation of the 
solute gradient:

.. math::
  \nabla C_{lm} \approx WR \, \frac{C_{lf} - C_{lm}}{L_m},

where 
:math:`WR` is the Warren-Root coefficient that estimates the poro-space geometry, [-]


Model specifications and assumptions
....................................

This list is used to summarize physical models and assumptions, such as
coupling with other PKs.
This list is often generated or extended by a high-level MPC PK.

* `"gas diffusion`" [bool] indicates that air-water partitioning coefficients
  are used to distribute components between liquid and as phases. Default is *false*.

* `"permeability field is required`" [bool] indicates if some transport features
  require absolute permeability. Default is *false*.

* `"multiscale model`" [string] specifies a multiscale model.
  Available options are `"single porosity`" (default) and `"dual porosity`".

* `"effective transport porosity`" [bool] If *true*, effective transport porosity
  will be used by dispersive-diffusive fluxes instead of total porosity. 
  Default is *false*.

.. code-block:: xml

  <ParameterList name="_TRANSPORT">  <!-- parent list -->
  <ParameterList name="physical models and assumptions">
    <Parameter name="gas diffusion" type="bool" value="false"/>
    <Parameter name="permeability field is required" type="bool" value="false"/>
    <Parameter name="multiscale model" type="string" value="single porosity"/>
    <Parameter name="effective transport porosity" type="bool" value="false"/>
  </ParameterList>
  </ParameterList>


Global parameters
.................

This list is used to summarize physical models and assumptions, such as
The transport component of Amanzi performs advection of aqueous and gaseous
components and their dispersion and diffusion. 
The main parameters control temporal stability, spatial 
and temporal accuracy, and verbosity:


* `"domain name`" [string] specifies mesh name that defined domain of this PK.
  Default is `"domain`".

* `"cfl`" [double] Time step limiter, a number less than 1. Default value is 1.
   
* `"spatial discretization order`" [int] defines accuracy of spatial discretization.
  It permits values 1 or 2. Default value is 1. 
  
* `"temporal discretization order`" [int] defines accuracy of temporal discretization.
  It permits values 1 or 2 and values 3 or 4 when expert parameter 
  `"generic RK implementation`" is set to true. Note that RK3 is not monotone.
  Default value is 1.

* `"reconstruction`" [list] collects reconstruction parameters. The available options are
  describe in the separate section below.

* `"solver`" [string] Specifies the dispersion/diffusion solver.

* `"preconditioner`" [string] specifies preconditioner for dispersion solver.

* `"number of aqueous components`" [int] The total number of aqueous components. 
  Default value is the total number of components.

* `"number of gaseous components`" [int] The total number of gaseous components. 
  Default value is 0.
   
.. code-block:: xml

  <ParameterList>  <!-- parent list -->
  <ParameterList name="_TRANSPORT">
    <Parameter name="domain name" type="string" value="domain"/>
    <Parameter name="cfl" type="double" value="1.0"/>
    <Parameter name="spatial discretization order" type="int" value="1"/>
    <Parameter name="temporal discretization order" type="int" value="1"/>
    <Parameter name="solver" type="string" value="_PCG_SOLVER"/>

    <ParameterList name="reconstruction">
      <Parameter name="method" type="string" value="cell-based"/>
      <Parameter name="polynomial order" type="int" value="1"/>
      <Parameter name="limiter" type="string" value="tensorial"/>
      <Parameter name="limiter extension for transport" type="bool" value="true"/>
    </ParameterList>

    <ParameterList name="verbose object">
      <Parameter name="verbosity level" type="string" value="high"/>
    </ParameterList>
  </ParameterList>  
  </ParameterList>  


Material properties
...................

The material properties include dispersivity model and diffusion parameters 
for aqueous and gaseous phases.
The dispersivity is defined as a soil property. 
The diffusivity is defined independently for each solute.

* _SOIL [list] Defines material properties.
  
  * `"region`" [Array(string)] Defines geometric regions for material SOIL.
  * `"model`" [string] Defines dispersivity model, choose exactly one of the following: `"scalar`", `"Bear`",
    `"Burnett-Frind`", or `"Lichtner-Kelkar-Robinson`".
  * `"parameters for MODEL`" [list] where `"MODEL`" is the model name.
    For model `"scalar`", *only* one of the following options must be specified:

      * `"alpha`" [double] defines dispersivity in all directions, [m].
      * `"dispersion coefficient`" [double] defines dispersion coefficient [m^2/s].

    For model `"Bear`", the following options must be specified:

      * `"alpha_l`" [double] defines dispersion in the direction of Darcy velocity, [m].
      * `"alpha_t`" [double] defines dispersion in the orthogonal direction, [m].
    
    For model `"Burnett-Frind`", the following options must be specified:

      * `"alphaL`" [double] defines the longitudinal dispersion in the direction of Darcy velocity, [m].
      * `"alpha_th`" [double] Defines the transverse dispersion in the horizonla direction orthogonal directions, [m].
      * `"alpha_tv`" [double] Defines dispersion in the orthogonal directions, [m].
        When `"alpha_th`" equals to `"alpha_tv`", we obtain dispersion in the direction of the Darcy velocity.
        This and the above parameters must be defined for `"Burnett-Frind`" and `"Lichtner-Kelkar-Robinson`" models.

    For model `"Lichtner-Kelker-Robinson`", the following options must be specified:

      * `"alpha_lh`" [double] defines the longitudinal dispersion in the horizontal direction, [m].
      * `"alpha_lv`" [double] Defines the longitudinal dispersion in the vertical direction, [m].
        When `"alpha_lh`" equals to `"alpha_lv`", we obtain dispersion in the direction of the Darcy velocity.
        This and the above parameters must be defined for `"Burnett-Frind`" and `"Lichtner-Kelker-Robinson`" models.
      * `"alpha_th`" [double] Defines the transverse dispersion in the horizontal direction orthogonal directions, [m].
      * `"alpha_tv" [double] Defines dispersion in the orthogonal directions.
        When `"alpha_th`" equals to `"alpha_tv`", we obtain dispersion in the direction of the Darcy velocity.
        This and the above parameters must be defined for `"Burnett-Frind`" and `"Lichtner-Kelker-Robinson`" models.

  * `"aqueous tortuosity`" [double] Defines tortuosity for calculating diffusivity of liquid solutes, [-].
  * `"gaseous tortuosity`" [double] Defines tortuosity for calculating diffusivity of gas solutes, [-].
 
Three examples are below:

.. code-block:: xml

  <ParameterList name="_TRANSPORT">  <!-- parent list -->
  <ParameterList name="material properties">
    <ParameterList name="_WHITE SOIL">
      <Parameter name="regions" type="Array(string)" value="{_TOP_REGION, _BOTTOM_REGION}"/>
      <Parameter name="model" type="string" value="Bear"/>
      <ParameterList name="parameters for Bear">
        <Parameter name="alpha_l" type="double" value="1e-2"/>
        <Parameter name="alpha_t" type="double" value="1e-5"/>
      <ParameterList>
      <Parameter name="aqueous tortuosity" type="double" value="1.0"/>
      <Parameter name="gaseous tortuosity" type="double" value="1.0"/>       
    </ParameterList>  
     
    <ParameterList name="_GREY SOIL">
      <Parameter name="regions" type="Array(string)" value="{_MIDDLE_REGION}"/>
      <Parameter name="model" type="string" value="Burnett-Frind"/>
      <ParameterList name="parameters for Burnett-Frind">
        <Parameter name="alpha_l" type="double" value="1e-2"/>
        <Parameter name="alpha_th" type="double" value="1e-3"/>
        <Parameter name="alpha_tv" type="double" value="2e-3"/>
      <ParameterList>
      <Parameter name="aqueous tortuosity" type="double" value="0.5"/>
      <Parameter name="gaseous tortuosity" type="double" value="1.0"/>       
    </ParameterList>  
  </ParameterList>  
  </ParameterList>  


* `"molecular diffusion`" [list] defines names of solutes in aqueous and gaseous phases and related
  diffusivity values.

.. code-block:: xml

  <ParameterList name="_TRANSPORT">  <!-- parent list -->
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

  <ParameterList name="_TRANSPORT">  <!-- parent list -->
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
discretized in mesh cells and on mesh faces. The later unknowns are auxiliary unknowns.


Multiscale continuum models
...........................

The list of multiscale models is the place for various subscale models that coul 
be mixed and matched.
Its ordered by materials and includes parameters for the assigned multiscale model
This list is optional.

* `"multiscale model`" [string] is the model name. Available option is `"dual porosity`"
  and `"generalized dual porosity`".

* `"regions`" [Array(string)] is the list of regions where this model should be applied.

* `"xxx parameters`" [sublist] provides parameters for the model specified by variable `"multiscale model`".


Dual porosity model
```````````````````

* `"Warren Root parameter`" [list] scales diffusive solute transport due to
  concentration gradient.
* `"tortousity`" [double] defines tortuosity to correct diffusivity of a liquid solute.


Generalized dual porosity model
```````````````````````````````

* `"number of matrix nodes`" [int] defines number of matrix layers.
* `"matrix depth`" [double] is the characteristic length for matrix continuum.
* `"tortousity`" [double] defines tortuosity to correct diffusivity of a liquid solute.
* `"matrix volume fraction`" [double] defines relative volume of matrix continuum.

.. code-block:: xml

  <ParameterList name="_TRANSPORT">  <!-- parent list -->
  <ParameterList name="multiscale models">
    <ParameterList name="_WHITE SOIL">
      <Parameter name="multiscale model" type="string" value="dual porosity"/>
      <Parameter name="regions" type="Array(string)" value="{_TOP_REGION, _BOTTOM_REGION}"/>
      <ParameterList name="dual porosity parameters">
        <Paramater name="Warren Root parameter" type="double" value="4.0e-5"/>
        <Paramater name="matrix tortuosity" type="double" value="0.95"/>
        <Paramater name="matrix volume fraction" type="double" value="0.9999"/>
      </ParameterList>  
    </ParameterList>  

    <ParameterList name="_GREY SOIL">
      <Parameter name="multiscale model" type="string" value="generalized dual porosity"/>
      <Parameter name="regions" type="Array(string)" value="{_MIDDLE_REGION}"/>
      <ParameterList name="generalized dual porosity parameters">
        <Paramater name="number of matrix nodes" type="int" value="2"/>
        <Paramater name="matrix depth" type="double" value="0.01"/>
        <Paramater name="matrix tortuosity" type="double" value="1.0"/>
      </ParameterList>  
    </ParameterList>  
  </ParameterList>  
  </ParameterList>  


Boundary conditions
...................

For the advective transport, the boundary conditions must be specified on inflow parts of the
boundary. If no value is prescribed through the XML input, the zero influx boundary condition
is used. Note that the boundary condition is set up separately for each component.
The structure of boundary conditions is aligned with that used for flow and
allows us to define spatially variable boundary conditions. 

* `"boundary conditions`" [list]

  * `"concentration`" [list] This is a reserved keyword.
   
    * "_SOLUTE" [list] contains a few sublists (e.g. _EAST_CRIB) for boundary conditions.
      The name _SOLUTE must be the name in the list of solutes.
 
      * "_BC1" [list] defines boundary conditions using arrays of boundary regions and attached
        functions.
   
      * `"regions`" [Array(string)] defines a list of boundary regions where a boundary condition
        must be applied.
      * `"spatial distribution method`" [string] defines the method for distributing
        data  over the specified regions. The available options are `"area`" or `"none`". 
      * `"boundary concentration`" [list] defines a function for calculating boundary conditions.
        The function specification is described in subsection Functions.

The example below sets constant boundary condition 1e-5 for the duration of transient simulation.

.. code-block:: xml

  <ParameterList name="_TRANSPORT">  <!-- parent list -->
  <ParameterList name="boundary conditions">
    <ParameterList name="concentration">
      <ParameterList name="NO3-"> 
        <ParameterList name="_EAST_CRIB">   <!-- user defined name -->
          <Parameter name="regions" type="Array(string)" value="{_TOP, _LEFT}"/>
          <Parameter name="spatial distribution method" type="string" value="none"/>
          <ParameterList name="boundary concentration">
            <ParameterList name="function-constant">  <!-- any time function -->
              <Parameter name="value" type="double" value="1e-5"/>
            </ParameterList>
          </ParameterList>
        </ParameterList>
        <ParameterList name="_WEST CRIB">   <!-- user defined name -->
          ...
        </ParameterList>
      </ParameterList>

      <ParameterList name="CO2"> <!-- Next component --> 
        ...
      </ParameterList>
    </ParameterList>
  </ParameterList>
  </ParameterList>

Geochemical boundary condition is the Dirichlet boundary condition
which requires calculation of a geochemical balance.
Note that the number of *time functions* below is one less than the number of times
and geochemical conditions.

.. code-block:: xml

  <ParameterList name="boundary conditions">  <!-- parent list -->
  <ParameterList name="geochemical">
    <ParameterList name="_EAST_CRIB">
      <Parameter name="solutes" type="Array(string)" value={H+,HCO3-,Ca++}"/>
      <Parameter name="times" type="Array(double)" value="{0.0, 100.0}"/>
      <Parameter name="geochemical conditions" type="Array(string)" value="{cond1, cond2}"/>
      <Parameter name="time functions" type="Array(string)" value="{constant}"/>
      <Parameter name="regions" type="Array(string)" value="{EAST_CRIB}"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>


Sources and sinks
.................

The sources and sinks are typically located at wells. 
Stability condition requires to distinguish between injecting and producing wells.
A source function used for an injecting well specifies concentration of solute.
A source function used for a producing well specifies volumetric flow rate [m^3/s] of water. 

The structure of list *source terms* includes only sublists named after components. 
Again, constant functions can be replaced by any available time-function.
Note that the source values are set up separately for each component.

* `"concentration`" [list] This is a reserved keyword.

 * "_SOLUTE" [list] contains a few sublists (e.g. _SRC1 and _SRC2) for various sources and sinks.
   The name _SOLUTE must exist in the list of solutes.

  * "_SRC1" [list] defines a source using arrays of domain regions, a function, and 
    a distribution method.
   
   * `"regions`" [Array(string)] defines a list of domain regions where a source term
     must be applied.

   * `"spatial distribution method`" [string] identifies a method for distributing
     source Q over the specified regions. The available options are `"volume`",
     `"none`", and `"permeability`". For option `"none`" the source term Q is measured
     in [mol/L/s] (if units for concetration is mol/L) or [mol/m^3/s] (otherwise). 
     For the other options, it is measured in [mol/s]. When the source function
     is defined over a few regions, Q will be distributed over their union.

   * `"submodel`" [string] refines definition of source. Available options are `"rate`"
     and `"integrand`". The first option defines rate of change `q`, the second one 
     defines integrand `Q` of a rate `Q = dq/dt`. Default is `"rate`".

   * `"sink`" [list] is a function for calculating a source term.
     The function specification is described in subsection Functions.


This example defines one well and one sink.

.. code-block:: xml

  <ParameterList name="source terms"> <!-- parent list -->
  <ParameterList name="concentration">
    <ParameterList name="H+"> 
      <ParameterList name="_SOURCE: EAST WELL">   <!-- user defined name -->
        <Parameter name="regions" type="Array(string)" value="{_EAST_WELL}"/>
        <Parameter name="spatial distribution method" type="string" value="volume"/>
        <Parameter name="submodel" type="string" value="rate"/>
        <ParameterList name="injector">   <!-- reserved keyword -->
          <ParameterList name="function-constant">
            <Parameter name="value" type="double" value="0.01"/>
          </ParameterList>
        </ParameterList>
      </ParameterList>
      <ParameterList name="_SOURCE: WEST_WELL">
         ...
      </ParameterList>
    </ParameterList>
     
    <ParameterList name="CO2(g)">   <!-- next component, a gas -->
      <ParameterList name="_SOURCE: WEST WELL">   <!-- user defined name -->
        <Parameter name="regions" type="Array(string)" value="{_WEST_WELL}"/>
        <Parameter name="spatial distribution method" type="string" value="permeability"/>
        <ParameterList name="injector">  
          <ParameterList name="function-constant">
            <Parameter name="value" type="double" value="0.02"/>
          </ParameterList>
        </ParameterList>
      </ParameterList>
    </ParameterList>
  </ParameterList>
  </ParameterList>
    

Developer parameters
....................

The remaining parameters that can be used by a developer include

* `"enable internal tests`" [bool] turns on various internal tests during
  run time. Default value is `"false`".
   
* `"generic RK implementation`" [bool] leads to generic implementation of 
  all Runge-Kutta methods. Default value is `"false`".
   
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
 * [local] species's name, concentration extrema, and total amount of it in the reservoir
 * [global] current simulation time (in years)

.. code-block:: xml

  CycleDriver      |   Cycle 10: time(y) = 0.803511, dt(y) = 0.089279
  TransportPK      |    cell 0 has smallest dt, (-270, -270)
  TransportPK      |    dispersion solver (pcg) ||r||=8.33085e-39 itrs=2
  TransportPK      |    1 sub-cycles, dt_stable=2.81743e+06 [sec]  dt_MPC=2.81743e+06 [sec]
  TransportPK      |    Tc-99: min=8.08339e-06 mol/L max=0.0952948 mol/L, total=9.07795 mol
  CycleDriver      |   New time(y) = 0.89279


Chemistry PK
------------

The chemistry header includes three parameters:

* `"chemistry model`" [string] defines chemical model. The available options are `"Alquimia`"
  and `"Amanzi`" (default).

.. code-block:: xml

  <ParameterList name="_CHEMISTRY">
    <Parameter name="chemistry model" type="string" value="Amanzi"/>
  </ParameterList>


Geochemical engines
...................

Here we specify either the default or the third-party geochemical engine. 

Common parameters
`````````````````

The following parameters are common for all supported engines.

* `"time step control method`" [string] specifies time step control method for chemistry subcycling. 
  Choose either "fixed" (default) or "simple".  For option "fixed", time step is fixed.
  For option "simple", the time step is adjusted in response to stiffness of system of equations 
  based on a simple scheme. This option require the following parameters: `"time step cut threshold`",
  `"time step cut factor`", `"time step increase threshold`", and `"time step increase factor`".

* `"time step cut threshold`" [int] is the number of Newton iterations that if exceeded
  will trigger a time step cut. Default is 8.

* `"max time step (s)`" [double] is the maximum time step that chemistry will allow the MPC to take.

* `"initial time step (s)`" [double] is the initial time step that chemistry will ask the MPC to take.

* `"time step cut factor`" [double] is the factor by which the time step is cut. Default is 2.0

* `"time step increase threshold`" [int] is the number of consecutive successful time steps that
  will trigger a time step increase. Default is 4.

* `"time step increase factor`" [double] is the factor by which the time step is increased. Default is 1.2

* `"free ion initial guess`" [double] provides an estimate of the free ion concentration for solutes.
  It used to help convergence of the initial solution of the chemistry. If this parameter is absent, 
  a fraction (10%) of the total component concentration is used.

* `"initial conditions time`" [double] specifies time for applying initial conditions. This parameter
  is useful for simulation restart. Default value is the state time when chemistry PK is instantiated. 

Alquimia
````````

The Alquimia chemistry process kernel only requires the *Engine* and *Engine Input File*
entries, but will also accept and respect the value given for *max time step (s)*. 
Most details are provided in the trimmed PFloTran file *1d-tritium-trim.in*.

* `"minerals`" [Array(string)] is the list of mineral names.

* `"sorption sites`" [Array(string)] is the list of sorption sites.

* `"auxiliary data`" [Array(string)] defines additional chemistry related data that the user 
  can request be saved to vis files. 

* `"min time step (s)`" [double] is the minimum time step that chemistry will allow the MPC to take.

.. code-block:: xml

  <ParameterList>  <!-- parent list -->
  <ParameterList name="_CHEMISTRY">
    <Parameter name="engine" type="string" value="PFloTran"/>
    <Parameter name="engine input file" type="string" value="_TRITIUM.in"/>
    <Parameter name="minerals" type="Array(string)" value="{quartz, kaolinite, goethite, opal}"/>
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

* `"thermodynamic database`" [list] 

  * `"file`" [string] is the name of the chemistry database file, relative to the execution directory.

  * `"format`" [string] is the format of the database file. Actual database format is not XML and 
    is the same as described for the 2010 demo with additions for the new chemical processes. 
    Valid values: "simple".

* `"minerals`" [Array(string)] is the list of mineral names.

* `"sorption sites`" [Array(string)] is the list of sorption sites.

* `"activity model`" [string] is the type of model used for activity corrections. 
  Valid options are `"unit`", `"debye-huckel`", and `"pitzer-hwm`",

* `"tolerance`" [double] defines tolerance in Newton solves inside the chemistry library.

* `"maximum Newton iterations`" [int] is the maximum number of iteration the chemistry 
  library can take.

* `"auxiliary data`" [Array(string)] defines additional chemistry related data that the user 
  can request be saved to vis files. Currently `"pH`" is the only variable supported.

.. code-block:: xml

  <ParameterList>  <!-- parent list -->
  <ParameterList name="_CHEMISTRY">
    <ParameterList name="thermodynamic database">
      <Parameter name="file" type="string" value="_TRITIUM.bgd"/>
      <Parameter name="format" type="string" value="simple"/>
    </ParameterList>
    <Parameter name="activity model" type="string" value="unit"/>
    <Parameter name="tolerance" type="double" value="1.5e-12"/>
    <Parameter name="maximum Newton iterations" type="int" value="25"/>
    <Parameter name="max time step (s)" type="double" value="1.5e+07"/>
    <Parameter name="auxiliary data" type="Array(string)" value="{pH}"/>
    <Parameter name="number of component concentrations" type="int" value="1"/>
    <Parameter name="time step control method" type="string" value="simple"/>
  </ParameterList>
  </ParameterList>


Format of chemistry database (.bgd) file
........................................

A section header starts with token `"<`". 
A comment line starts with token `"#`". 
Data fields are separated by semicolumns.


Primary species
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


General kinetics
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


Aqueous equilibrium complexes
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
Name = coeff reactant, log Keq, gram molecular weight [g/mol], molar volume [cm^3/mol],
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


Mineral kinetics
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


Ion exchange sites
``````````````````

Each line in this section has three fields: 
exchanger name, exchanger change, and exchanger location. 
The location is the mineral where the exchanger is located, i.e. kaolinite.

.. code-block:: txt

  <Ion Exchange Sites
   X- ; -1.0 ; Halite


Ion exchange complexes
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


Surface complex sites
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


Surface complexes
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


Radiactive decay
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
:math:`\varepsilon` is the energy density [:math:`J/m^3`],
:math:`\eta_l` is molar density of liquid [:math:`mol/m^3`],
:math:`Q` is heat source term,
:math:`\boldsymbol{q}_l` is the Darcy velocity [m/s],
:math:`\kappa` is thermal conductivity,
and :math:`H_l` is molar enthalpy of liquid [J/mol].
We define 

.. math::
   \varepsilon = \phi (\eta_l s_l U_l + \eta_g s_g U_g) + 
   (1 - \phi) \rho_r c_r T

where
:math:`s_l` is liquid saturation [-],
:math:`s_g` is gas saturation (water vapor),
:math:`\eta_l` is molar density of liquid [:math:`mol/m^3`],
:math:`\eta_g` is molar density of gas,
:math:`U_l` is molar internal energy of liquid [J/mol],
:math:`U_g` is molar internal energy of gas (water vapor) [J/mol],
:math:`\phi` is porosity [-],
:math:`\rho_r` is rock density [:math:`kg/m^3`],
:math:`c_r` is specific heat of rock [J/kg/K],
and :math:`T` is temperature [K].


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

  <ParameterList>  <!-- parent list -->
  <ParameterList name="_ENERGY"> 
    <ParameterList name="physical models and assumptions">
      <Parameter name="vapor diffusion" type="bool" value="false"/>
      <Parameter name="water content model" type="string" value="constant density"/>
    </ParameterList>
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

  <ParameterList name="_ENERGY">  <!-- parent list -->
  <ParameterList name="energy evaluator">
    <Parameter name="energy key" type="string" value="energy"/>
    <Parameter name="evaluator type" type="string" value="constant liquid density"/>
    <Parameter name="vapor diffusion" type="bool" value="true"/>
    <ParameterList name="verbose object">
      <Parameter name="verbosity level" type="string" value="high"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>


Molar enthalpy
..............

.. code-block:: xml

  <ParameterList name="_ENERGY">  <!-- parent list -->
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

  <ParameterList name="_ENERGY">  <!-- parent list -->
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

  <ParameterList name="_ENERGY">  <!-- parent list -->
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
It also has one global parameter.

* `"operators`" [list] 

  * `"include enthalpy in preconditioner`" [bool] allows us to study impact (usually positive) 
    of including enthalpy term in the preconditioner. Default value is *true*.
  
  * `"diagonal shift`" [double] allows for a constant shift to be applied to
    the diagonal of the assembled operator, which can b useful for dealing
    with singular or near-singular matrices.  Default is *0.0*.


Diffusion operator
``````````````````

Operators sublist describes the PDE structure of the flow, specifies a discretization
scheme, and selects assembling schemas for matrices and preconditioners.

* `"diffusion operator`" [list] defines parameters for generating and assembling diffusion matrix.

  * `"matrix`" [list] defines parameters for generating and assembling diffusion matrix. See section
    describing operators. 

  * `"preconditioner`" [list] defines parameters for generating and assembling diffusion 
    matrix that is used to create preconditioner. 
    Since update of preconditioner can be lagged, we need two objects called `"matrix`" and `"preconditioner`".

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
        <Parameter name="Newton correction" type="string" value="approximate Jacobian"/>
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
    <Parameter name="method" type="string" value="upwind"/>
  </ParameterList>
  </ParameterList>


Sources and sinks
.................

The sources and sinks for injecting and removing energy from the system. 
Negative source removes energy. 
Positive source inject energy.
The structure of list *source terms* mimics that of list *boundary conditions*. 
Again, constant functions can be replaced by any of the available functions.

* `"regions`" [Array(string)] is the list of regions where the source is defined.

* `"spatial distribution method`" [string] is the method for distributing
  source Q over the specified regions. The available options are `"volume`" and
  `"none`".
  For option `"none`", the source term function Q is measured in [J/m^3/s]. 
  For option `"volume`", it is measured in [J/s]. 
  When the source function is defined over a few regions, Q is distributed over their union.
  Option `"volume fraction`" can be used when the region geometric
  model support volume fractions. 

* `"use volume fractions`" instructs the code to use all available volume fractions. 
  Note that the region geometric model supports volume fractions only for a few regions.

* `"submodel`" [string] refines definition of the source. Available options are `"rate`",
  `"integrated source`". The first option defines the source 
  in a natural way as the rate of change `q`. The second option defines the indefinite
  integral `Q` of the rate of change, i.e. the source term is calculated as `q = dQ/dt`. 
  Default is `"rate`". 

.. code-block:: xml

  <ParameterList name="energy">  <!-- parent list -->

    <ParameterList name="source terms">
      <ParameterList name="_SRC 0">
        <Parameter name="regions" type="Array(string)" value="{_WELL_EAST}"/>
        <Parameter name="spatial distribution method" type="string" value="volume"/>
        <Parameter name="submodel" type="string" value="rate"/>
        <ParameterList name="source">
          <ParameterList name="function-constant">
            <Parameter name="value" type="double" value="-0.1"/>
          </ParameterList>
        </ParameterList>
      </ParameterList>
    </ParameterList>
    
   </ParameterList>


  

Navier Stokes PK
----------------

The conceptual PDE model for the incompressible Navier Stokes equations are

.. math::
  \frac{\partial (\rho \boldsymbol{u})}{\partial t} 
  + \boldsymbol{\nabla} \cdot (\rho \boldsymbol{u} \otimes \boldsymbol{u})
  =
  - \boldsymbol{\nabla} p 
  + \boldsymbol{\nabla} \cdot \boldsymbol{\sigma} 
  + \rho \boldsymbol{g}

where 
:math:`\rho` is the fluid density [kg/m^3],
:math:`p` is the pressure [Pa],
:math:`\boldsymbol{\sigma}` is the deviatoric stress tensor,
:math:`\boldsymbol{g}` is the gravity vector [:math:`m/s^2`], 
and :math:`u \otimes v = u \times v^T`.
The Stokes stress contitutive law for incompressible viscous fluid is

.. math::
  \boldsymbol{\sigma} = 
  \mu \left(\boldsymbol{\nabla} \boldsymbol{u} + 
            \boldsymbol{\nabla} \boldsymbol{u}^{T}\right),

where 
:math:`\mu` is the dynamic viscosity [:math:`Pa \cdot s`]. It can depend on density and pressure. 


Physical models and assumptions
...............................

This list is used to summarize physical models and assumptions, such as
coupling with other PKs.
This list is often generated on a fly by a high-level MPC PK.

* `"gravity`" [bool] is set up automatically by a high-level PK.
  The default value is `"false`".

.. code-block:: xml

  <ParameterList name="_NAVIER_STOKES">  <!-- parent list -->
  <ParameterList name="physical models and assumptions">
    <Parameter name="gravity" type="bool" value="false"/>
  </ParameterList>
  </ParameterList>


Operators
.........

This section contains sublist for diffusion and advection operators.
It also has one global parameter.

* `"operators`" [list] 


Elasticity operator
```````````````````

.. code-block:: xml

  <ParameterList name="operators">  <!-- parent list -->
  <ParameterList name="elasticity operator">
    <Parameter name="method" type="string" value="BernardiRaugel"/>
    <ParameterList name="schema">
      <Parameter name="base" type="string" value="cell"/>
      <Parameter name="location" type="Array(string)" value="{node, face}"/>
      <Parameter name="type" type="Array(string)" value="{vector, normal component}"/>
      <Parameter name="number" type="Array(int)" value="{2, 1}"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>


Convection operator
```````````````````

This section to be written.

.. code-block:: xml

  <ParameterList name="operators">  <!-- parent list -->
  <ParameterList name="convection operator">
    <Parameter name="flux formula" type="string" value="Rusanov"/>
  </ParameterList>
  </ParameterList>


Divergence operator
```````````````````

This section to be written.

.. code-block:: xml

  <ParameterList name="operators">  <!-- parent list -->
  <ParameterList name="divergence operator">
    <Parameter name="method" type="string" value="BernardiRaugel"/>
    <ParameterList name="schema domain">
      <Parameter name="base" type="string" value="cell"/>
      <Parameter name="location" type="Array(string)" value="{node, face}"/>
      <Parameter name="type" type="Array(string)" value="{scalar, normal component}"/>
      <Parameter name="number" type="Array(int)" value="{2, 1}"/>
    </ParameterList>
    <ParameterList name="schema range">
      <Parameter name="base" type="string" value="cell"/>
      <Parameter name="location" type="Array(string)" value="{cell}"/>
      <Parameter name="type" type="Array(string)" value="{scalar}"/>
      <Parameter name="number" type="Array(int)" value="{1}"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>


Multiphase PK
-------------

Mathematical models
...................

The conceptual PDE model for the isothermal multiphase flow inlcude
transport equations for components and nonlinear algebraic constraints for the phase presence. 
At the moment we consider two phases (liquid and gas), multiple components, and one 
constraint.
Each transport equation has the following form:

.. math::
  \frac{\partial \Theta}{\partial t}
  + \nabla \cdot \boldsymbol{\Psi} = Q,

where 
:math:`\Theta` is the storage and 
:math:`\boldsymbol{\Psi}` is the flux.
The storage term sums up component amount across two phases, :math:`\alpha=l` for liquid
phase and :math:`\alpha=g` for gas phase:

.. math::
  \Theta = \phi \sum_\alpha \eta_\alpha x_\alpha s_\alpha

where
:math:`\phi` is porosity [-],
:math:`\eta` is molar density [mol/m^3],
:math:`x` is molar fraction of component [-], and
:math:`s` is phase saturation [-].

The flux includes advective and diffusion terms:

.. math::
  \boldsymbol{\Psi} 
  = -\sum_\alpha \eta_\alpha \left(\boldsymbol{q}_\alpha + D_\alpha \nabla x \right)

where
:math:`\boldsymbol{q}` is Darcy phase velocity,
:math:`D` is molecular duffusion coefficient.

The nonlinear algebraic constraint may have different forms. One of the available forms is

.. math::
  min (s_g, 1 - x_l - x_g) = 0.

It implies that if gas compounent is present then we must have :math:`x_l + x_g = 1`.

The PK provides three choices of primary variables. 
The first one includes pressure liquid, mole gas fraction, and saturation liquid.
The second one includes pressure liquid, molar gas density, and saturation liquid.
The third one is used for verification purposes and is based on the model in Jaffre's paper. 
This model describes two-phase two-component system with water and hydrogen. 


Shallow water PK
----------------

The mathematical model describing two-dimensional shallow water flow is

.. math::
  \begin{align*}
  & h_t + (hu)_x + (hv)_y = 0, \\
  & (hu)_t + (hu^2 + \frac{1}{2} gh^2)_x + (huv)_y = -ghB_x \\
  & (hv)_t + (huv)_x + (hv^2 + \frac{1}{2} gh^2)_y = -ghB_y
  \end{align*}  

Here
:math:`h` [m] is water depth, 
:math:`g` [m/s^2] is gravity acceleration,
:math:`u` [m/s] is depth averaged velocity in x direction,
:math:`v` [m/s] is depth averaged velocity in y direction,
:math:`B` [m] is bottom elevation (bathymetry),
:math:`H = h + B` [m] is water surface elevation.


Global parameters
.................

Global parameters are placed in the sublist `"shallow water`". 
The list of global parameters include:

* `"domain name`" [string] specifies mesh name that defined domain of this PK.
  Default is `"domain`".

* `"limiter cfl`" [double] is a safety factor (less than 1) applied to the limiter.
  Default value is 1.


Reconstruction and limiters
...........................

The control of the second-order numerical scheme is done via `"reconstruction`"
sublist, described in Reconstruction_. Here is the example:


.. code-block:: xml

  <ParameterList name="shallow water">  <!-- parent list -->
  <ParameterList name="reconstruction">
    <Parameter name="method" type="string" value="cell-based"/>
    <Parameter name="polynomial order" type="int" value="1"/>
    <Parameter name="limiter" type="string" value="Barth-Jespersen"/>
    <Parameter name="limiter stencil" type="string" value="cell to closest cells"/>
    <Parameter name="limiter location" type="string" value="node"/>
    <Parameter name="limiter points" type="int" value="0"/>
    <Parameter name="limiter cfl" type="double" value="0.1"/>
  </ParameterList>
  </ParameterList>


Coupled process kernels
=======================

Coupling of process kernels requires additional parameters for PK 
described above.


Reactive transport PK
---------------------

Reactive transport can be setup using a steady-state flow.
The two PKs are executed consequitively. 
The input spec requires new keyword *reactive transport*.

.. code-block:: xml

  <ParameterList name="PK tree">  <!-- parent list -->
  <ParameterList name="_REACTIVE TRANSPORT">
    <Parameter name="PK type" type="string" value="reactive transport"/>
    <ParameterList name="_TRANSPORT">
      <Parameter name="PK type" type="string" value="transport"/>
    </ParameterList>
    <ParameterList name="_CHEMISTRY">
      <Parameter name="PK type" type="string" value="chemistry amanzi"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>


Flow and reactive transport PK
------------------------------

Amanzi uses operator splitting approach for coupled physical kernels.
The coupling of PKs is described as a tree where flow and reactive 
transport are executed consequitively.
The input spec requires new keyword *flow reactive transport*.

.. code-block:: xml

  <ParameterList name="PK tree">  <!-- parent list -->
  <ParameterList name="_FLOW and REACTIVE TRANSPORT">
    <Parameter name="PK type" type="string" value="flow reactive transport"/>
    <ParameterList name="_FLOW">
      <Parameter name="PK type" type="string" value="darcy"/>
    </ParameterList>
    <ParameterList name="_REACTIVE TRANSPORT">
      <Parameter name="PK type" type="string" value="reactive transport"/>
      <ParameterList name="_TRANSPORT">
      <Parameter name="PK type" type="string" value="transport"/>
      </ParameterList>
      <ParameterList name="_CHEMISTRY">
        <Parameter name="PK type" type="string" value="chemistry amanzi"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>
  </ParameterList>

This example describe four PKs identified by keywords *darcy*, *reactive transport*,
*transport*, and *chemistry amanzi*. 
The flow is fully saturated. 
The transport of reactive chemicals is based on the native chemistry package *chemistry amanzi*.

Details of PKs are organized as a plain list of ParameterLists.
Note that *reactive transport* is MPC-PK and hence its description is short.

.. code-block:: xml

  <ParameterList name="PKs">
  <ParameterList name="_FLOW and REACTIVE TRANSPORT">
    <Parameter name="PK type" type="string" value="flow reactive transport"/>
    <Parameter name="PKs order" type="Array(string)" value="{_FLOW, _REACTIVE TRANSPORT}"/>
    <Parameter name="master PK index" type="int" value="0"/>
  </ParameterList>

  <ParameterList name="_REACTIVE TRANSPORT">
    <Parameter name="PK type" type="string" value="reactive transport"/>
    <Parameter name="PKs order" type="Array(string)" value="{_CHEMISTRY, _TRANSPORT}"/>
  </ParameterList>

  <ParameterList name="_FLOW">
    ...
  </ParameterList>

  <ParameterList name="_TRANSPORT">
    ...
  </ParameterList>

  <ParameterList name="_CHEMISTRY">
    ...
  </ParameterList>
  </ParameterList>


Thermal flow PK
---------------

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
:math:`\theta` is total water content (we use non-conventional definition) [:math:`mol/m^3`],
:math:`\eta_l` is molar density of liquid [:math:`mol/m^3`],
:math:`\rho_l` is fluid density [:math:`kg/m^3`],
:math:`Q_1` is source or sink term,
:math:`\boldsymbol{q}_l` is the Darcy velocity [m/s],
:math:`k_r` is relative permeability [-],
:math:`\boldsymbol{g}` is gravity [:math:`m/s^2`],
:math:`\phi` is porosity [-],
:math:`s_g` is gas saturation (water vapor) [-],
:math:`\tau_g` is tortuosity of gas [-],
:math:`D_g` is diffusion coefficient,
and :math:`X_g` is molar fraction of water in the gas phase [-].
We define 

.. math::
   \theta = \phi (s_g \eta_g X_g + s_l \eta_l)

where
:math:`s_l` is liquid saturation [-],
and :math:`\eta_g` is molar density of gas.

In the second equation,
:math:`\varepsilon` is the energy density [:math:`J/mol^3`],
:math:`Q_2` is source or sink term,
:math:`\kappa` is thermal conductivity [W/m/K],
:math:`H_l` is molar enthalphy of liquid [J/mol],
and :math:`T` is temperature [K].
We define 

.. math::
   \varepsilon = \phi (\eta_l s_l U_l + \eta_g s_g U_g) + 
   (1 - \phi) \rho_r c_r T

where
:math:`U_l` is molar internal energy of liquid [J/mol],
:math:`U_g` is molar internal energy of gas (water vapor) [J/mol],
:math:`\rho_r` is rock density [kg/m^3],
and :math:`c_r` is specific heat of rock [J/kg/K].



Diffusion operator
..................

.. code-block:: xml

  <ParameterList name="_NAVIER_STOKES">  <!-- parent lists -->
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
       <Parameter name="Newton correction" type="string" value="none"/>
     </ParameterList>
   </ParameterList>
   </ParameterList>



Coupled matrix-fracture Darcy flow PK
-------------------------------------

Mathematical models
...................

Let subscripts :math:`m` and :math:`f` correspond to matrix and fracture, respectively.
The conceptual PDE model of the stationary coupled matrix-fracture flow is

.. math::
  \begin{array}{l}
  \phi_m \frac{S_{s,m}}{g} \frac{\partial p_m}{\partial t}
  - \boldsymbol{\nabla} \cdot (\rho_l \boldsymbol{q}_m) = Q_m,
  \quad
  \boldsymbol{q}_m = -\frac{\boldsymbol{K}_m}{\mu} 
  (\boldsymbol{\nabla} p_m - \rho_l \boldsymbol{g}) \\
  %
  \phi_f \frac{S_{s,f}}{g} \frac{\partial p_f}{\partial t}
  -\boldsymbol{\nabla} \cdot (\rho_l \boldsymbol{q}_f) = 
    -\rho_l [[ \tilde{\boldsymbol{q}}_m \cdot \boldsymbol{n} ]],
  \quad
  \boldsymbol{q}_f = -\frac{\boldsymbol{K}_f}{\mu} 
  (\boldsymbol{\nabla} p_f - \rho_l \boldsymbol{g}) \\
  %
  \tilde{\boldsymbol{q}}_m \cdot \boldsymbol{n} = \frac{k}{g} (p_f - p_m)
  \end{array}

subject to convential boundary conditions for both matrix and fracture domains expect for 
the matrix-fracture boundary where the boundary condition is

.. math::
  \boldsymbol{q}_m \cdot \boldsymbol{n} = \tilde{\boldsymbol{q}}_m \cdot \boldsymbol{n}

Here
:math:`\rho_l` is fluid density [kg/m^3],
:math:`\phi` is porosity [-],
:math:`S_s` is ispecific storage [m],
:math:`p` is aqueous pressure [Pa],
:math:`\boldsymbol{K}` is absolute permeability [m^2] for matrix domain and [m^3] for fracture domain,
:math:`Q_m` is source or sink term,
:math:`\boldsymbol{q}` is the Darcy velocity [m/s] for matrix domain and [m^2/s] for fracture domain,
:math:`k` is effective normal premeability [s^-1],
and
:math:`\boldsymbol{g}` is gravity [:math:`m/s^2`].


Main parameters and sublists 
............................

* `"PKs order`" [array(string)] defines user names for two flow PKs. The matrix PK should be
  defined *first*.

* `"time integrator`" [list] defines a generic time integrator used by the cycle driver. 

.. code-block:: xml

  <ParameterList name="PKs">  <!-- parent list -->
  <ParameterList name="_COUPLED DARCY FLOW">
    <Parameter name="PKs order" type="Array(string)" value="{_FLOW MATRIX, _FLOW FRACTURE}"/>
    <Parameter name="master PK index" type="int" value="0"/>
    <ParameterList name="time integrator">
      ...
    </ParameterList>
  </ParameterList>
  </ParameterList>


Coupled matrix-fracture tansport PK
-----------------------------------

Mathematical models
...................

Let subscripts :math:`m` and :math:`f` correspond to matrix and fracture, respectively.
The conceptual PDE model of the coupled matrix-fracture advective transport is

.. math::
  \begin{array}{l}
  \frac{\partial(\phi_m C_m)}{\partial t} = 
    -\boldsymbol{\nabla} \cdot (\boldsymbol{q}_m C_m) = Q_m,\\
  %
  \frac{\partial(\varepsilon_f\phi_f C_f)}{\partial t} = 
    -\boldsymbol{\nabla} \cdot (\boldsymbol{q}_f C_f)
    +[[ \tilde{C} (\boldsymbol{q}_m \cdot \boldsymbol{n}) ]] + Q_f,
  \end{array}

subject to the  Dirichlet boundary conditions on the inflow part of the computational domain.
Here
:math:`\phi` is porosity [-],
:math:`Q` is source or sink term,
:math:`\boldsymbol{q}` is the Darcy velocity [m/s] for matrix domain and [m^2/s] for fracture domain,
:math:`\varepsilon_f` is fracture cross-sectional length [m],
and
:math:`\tilde{C}` is concentration in fracture for outflow and in matrix for inflow.


Main parameters and sublists 
............................

* `"PKs order`" [array(string)] defines user names for two transport PKs.

.. code-block:: xml

  <ParameterList name="PKs">  <!-- parent list -->
  <ParameterList name="_COUPLED DARCY FLOW">
    <Parameter name="PKs order" type="Array(string)" value="{_TRANSPORT MATRIX, _TRANSPORT FRACTURE}"/>
    <Parameter name="master PK index" type="int" value="0"/>
  </ParameterList>
  </ParameterList>


Generic capabilities
====================

Collection of generic tools used by PKs.

Operators
---------

Operators are discrete forms of linearized PDEs operators.
They form a layer between physical process kernels and solvers
and include accumulation, diffusion, advection, elasticity, reaction, 
and source operators.
The residual associated with an operator :math:`L_h` helps to 
understand the employed sign convention:

.. math::
  r = f - L_h u.

A PK decides how to bundle operators in a collection of operators.
For example, an advection-diffusion problem may benefit from using
a single operator that combines two operators representing diffusion and advection process.
Collection of operators must be used for implicit solvers and for building preconditioners.
In such a case, the collections acts as a single operator.

Operators use a few tools that are generic in nature and can be used independently by PKs. 
The list includes reconstruction and limiting algorithms. 


Schema
......

The operators use notion of schema to describe operator's abstract structure.
Old operators use a simple schema which is simply the list of geometric objects where
scalar degrees of freedom are defined.
New operators use a list to define location, type, and number of degrees of freedom.
In addition, the base of local stencil is either *face* or *cell*.
A rectangular operator needs two schemas do describe its domain (called `"schema domain`") 
and its range (called `"schema range`").
A square operator may use either two identical schema lists or a single list called `"schema`".

.. code-block:: xml

  <ParameterList name="pks operator name">  <!-- parent list-->
  <ParameterList name="schema domain">
    <Parameter name="base" type="string" value="cell"/>
    <Parameter name="location" type="Array(string)" value="{node, face}"/>
    <Parameter name="type" type="Array(string)" value="{scalar, normal component}"/>
    <Parameter name="number" type="Array(int)" value="{2, 1}"/>
  </ParameterList>
  <ParameterList name="schema domain">
    <Parameter name="base" type="string" value="cell"/>
    <Parameter name="location" type="Array(string)" value="{node, face}"/>
    <Parameter name="type" type="Array(string)" value="{scalar, normal component}"/>
    <Parameter name="number" type="Array(int)" value="{2, 1}"/>
  </ParameterList>
  </ParameterList>

This example describes a square operator with two degrees of freedom per mesh node and one
degree of freedom per mesh face. 
The face-based degree of freedom is the normal component of a vector field. 
Such set of degrees of freedom is used in the Bernardi-Raugel element for discretizing 
Stokes equations.
Parameter `"base`" indicates that local matrices are associated with mesh cells. 


Diffusion operator
..................

Diffusion is the most frequently used operator. It employs the old schema.

* `"pks operator name`" [list] a PK specific name for the diffusion operator.

  * `"discretization primary`" [string] specifies an advanced discretization method that
    has useful properties under some a priori conditions on the mesh and/or permeability tensor.
    The available options are `"mfd: optimized for sparsity`", `"mfd: optimized for monotonicity`",
    `"mfd: default`", `"mfd: support operator`", `"mfd: two-point flux approximation`",
    `"fv: default`", and `"nlfv: default`".
    The first option is recommended for general meshes.
    The second option is recommended for orthogonal meshes and diagonal absolute 
    permeability tensor. 

  * `"discretization secondary`" [string] specifies the most robust discretization method
    that is used when the primary selection fails to satisfy all a priori conditions.
    Default value is equal to that for the primary discretization.

  * `"diffusion tensor`" [string] specifies additional properties of the diffusion tensor.
    It allows us to solve problems with non-symmetric but positive definite tensors. 
    Available options are *symmetric* (default) and *nonsymmetric*.

  * `"nonlinear coefficient`" [string] specifies a method for treating nonlinear diffusion
    coefficient, if any. Available options are `"none`", `"upwind: face`", `"divk: cell-face`" (default),
    `"divk: face`", `"standard: cell`", and `"divk: cell-face-twin`".
    Symmetry preserving methods are the divk-family of methods and the classical cell-centered
    method (`"standard: cell`"). The first part of the name indicates the base scheme.
    The second part (after the semi-column) indicates required components of the composite vector
    that must be provided by a physical PK.
    Default is `"none`".

  * `"schema`" [Array(string)] defines the operator stencil. It is a collection of 
    geometric objects. It equals to `"{cell}`" for finite volume schemes. 
    It is typically `"{face, cell}`" for mimetic discretizations.

  * `"preconditioner schema`" [Array(string)] defines the preconditioner stencil.
    It is needed only when the default assembling procedure is not desirable. 
    If skipped, the `"schema`" is used instead. 

  * `"gravity`" [bool] specifies if flow is driven also by the gravity.

  * `"gravity term discretization`" [string] selects a model for discretizing the 
    gravity term. Available options are `"hydraulic head`" [default] and `"finite volume`". 
    The first option starts with equation for the shifted solution, i.e. the hydraulic head,
    and derives gravity discretization by the reserve shifting.
    The second option is based on the divergence formula.

  * `"gravity magnitude`" [double] defined magnitude of the gravity vector.

  * `"Newton correction`" [string] specifies a model for correction (non-physical) terms 
    that must be added to the preconditioner. These terms approximate some Jacobian terms.
    Available options are `"true Jacobian`" and `"approximate Jacobian`".
    The FV scheme accepts only the first options. The othre schemes accept only the second option.

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

  * `"fracture`" [Array(string)] provides list of regions that defines a fracture network.
    This parameter is used only by the coupled flow PK.

.. code-block:: xml

  <ParameterList name="pks operator name">
    <Parameter name="discretization primary" type="string" value="mfd: optimized for monotonicity"/>
    <Parameter name="discretization secondary" type="string" value="mfd: two-point flux approximation"/>
    <Parameter name="schema" type="Array(string)" value="{face, cell}"/>
    <Parameter name="preconditioner schema" type="Array(string)" value="{face}"/>
    <Parameter name="gravity" type="bool" value="true"/>
    <Parameter name="gravity term discretization" type="string" value="hydraulic head"/>
    <Parameter name="gravity magnitude" type="double" value="9.81"/>
    <Parameter name="nonlinear coefficient" type="string" value="upwind: face"/>
    <Parameter name="Newton correction" type="string" value="true Jacobian"/>

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

A high-order advection operator may have different domain and range and therefore requires two schemas.
The structure of the new schema is described in the previous section.
A high-order advection operator has two terms in a weak formulation, corresponding to 
volume and surface integrals. These two terms are discretixed using two operators with
matrix of types *advection* and *flux*, respectively.


* `"pks operator name`" [list] a PK specific name for the advection operator.

  * `"method`" [string] defines a discretization method. The available option is `"dg modal`".

  * `"method order`" [int] defines method order. For example, the classical low-order finite 
    volume scheme is equivalent to DG of order 0.

  * `"matrix type`" [string] defines matrix type. The supported options are `"advection`"
    and `"flux`".

  * `"dg basis`" [string] defines bases for DG schemes. The available options are 
    `"regularized`" (recommended), `"normalized`", `"orthonormalized`", and `"natural`" 
    (not recommended).

  * `"gradient operator on test function`" [bool] defines place of the gradient operator.
    For integration by parts schemes, the gradient is transfered to a test function.
    This option is needed for discretizing volumetric integrals.

  * `"jump operator on test function`" [bool] defines place of the jump operator.
    For integration by parts schemes, the jump operator is applied to a test function.
    This option is needed for discretizing surface fluxes.

  * `"flux formula`" [string] defines type of the flux. The available options 
    are `"Rusanov`" (default), `"upwind`", `"downwind`", and `"NavierStokes`".

  * `"schema domain`" [list] defines a discretization schema for the operator domain.

  * `"schema range`" [list] defines a discretization schema for the operator range. 

.. code-block:: xml

  <ParameterList name="pks operator name">
    <Parameter name="method" type="string" value="dg modal"/>
    <Parameter name="method order" type="int" value="2"/>
    <Parameter name="flux formula" type="string" value="Rusanov"/>
    <Parameter name="matrix type" type="string" value="flux"/>
    <Parameter name="jump operator on test function" type="bool" value="true"/>

    <ParameterList name="schema domain">
      <Parameter name="base" type="string" value="cell"/>
      <Parameter name="location" type="Array(string)" value="{node, face}"/>
      <Parameter name="type" type="Array(string)" value="{scalar, normal component}"/>
      <Parameter name="number" type="Array(int)" value="{2, 1}"/>
    </ParameterList>
    <ParameterList name="schema range">
      <Parameter name="base" type="string" value="cell"/>
      <Parameter name="location" type="Array(string)" value="{cell}"/>
      <Parameter name="type" type="Array(string)" value="{scalar}"/>
      <Parameter name="number" type="Array(int)" value="{1}"/>
    </ParameterList>
  </ParameterList>

In this example, we construct an operator for volumetric integrals in a weak formulation
of advection problem.

The only low-order advection operator in Amanzi is the upwind operator. 
It employes the old schema.

.. code-block:: xml

  <ParameterList name="pks operator name">
    <Parameter name="base" type="string" value="face"/>
    <Parameter name="schema" type="Array(string)" value="{cell}"/>
    <Parameter name="method order" type="int" value="0"/>
    <Parameter name="matrix type" type="string" value="advection"/>
  </ParameterList>


Reaction operator
.................

A reaction operator may represent either reaction of identity operator.
It is symmetric so far and requires one schema.
The structure of the schema is described in the previous section.

* `"pks operator name`" [list] a PK specific name for the advection operator.

  * `"method`" [string] defines a discretization method. The only supported
    option is `"dg nodal`".

  * `"schema`" [list] defines a discretization schema for the operator domain.

.. code-block:: xml

  <ParameterList name="pks operator name">
    <Parameter name="method" type="string" value="dg modal"/>
    <Parameter name="method order" type="int" value="1"/>
    <Parameter name="matrix type" type="string" value="mass"/>
    <ParameterList name="schema">
      <Parameter name="base" type="string" value="cell"/>
      <Parameter name="location" type="Array(string)" value="{cell}"/>
      <Parameter name="type" type="Array(string)" value="{scalar}"/>
      <Parameter name="number" type="Array(int)" value="{3}"/>
    </ParameterList>
  </ParameterList>


Elasticity operator
...................

Elasticity operator is used for describing soil deformation or fluid flow (Stokes 
and Navier-Stokes).

* `"method`" [string] defines a discretization method. The available
  options are `"BernardiRaugel`".

* `"schema`" [list] defines a discretization schema.

  * `"location`" [Array(string)] defines geometric location of degrees of freedom.

  * `"type`" [Array(string)] defines type of degrees of freedom. The available options 
    are `"scalar`" and `"normal component`".

  * `"number`" [Array(int)] indicates how many time this degree of freedom is repeated.

.. code-block:: xml

  <ParameterList name="elasticity operator">
    <Parameter name="method" type="string" value="BernardiRaugel"/>
    <ParameterList name="schema">
      <Parameter name="base" type="string" value="cell"/>
      <Parameter name="location" type="Array(string)" value="{node, face}"/>
      <Parameter name="type" type="Array(string)" value="{scalar, normal component}"/>
      <Parameter name="number" type="Array(int)" value="{2, 1}"/>
    </ParameterList>
  </ParameterList>


Abstract operator
.................
An abstract operator is designed for testing new discretization methods. 
It uses the factory of discretization methods and a few control parameters
required by this factory and/or particular method in it.

* `"method`" [string] defines a discretization method. The available
  options are `"diffusion`", `"diffusion generalized`", `"BernardiRaugel`",
  `"CrouzeixRaviart`", `"CrouzeixRaviart serendipity`", `"Lagrange`", 
  `"Lagrange serendipity`", and `"dg modal`".

* `"method order`" [int] defines disretization order. It is used by 
  high-order discretization methods such as the discontinuous Galerkin.

* `"matrix type`" [string] defines type of local matrix. Available options are
  `"mass`", `"mass inverse`", `"stiffness`", `"divergence`", and `"advection`".

.. code-block:: xml

  <ParameterList name="_ABSTRACT OPERATOR">
    <Parameter name="method" type="string" value="dg modal"/>
    <Parameter name="method order" type="int" value="2"/>
    <Parameter name="dg basis" type="string" value="regularized"/>
    <Parameter name="matrix type" type="string" value="flux"/>

    <ParameterList name="schema domain">
      ...
    </ParameterList>
    <ParameterList name="schema range">
      ...
    </ParameterList>
  </ParameterList>


Diffusion is the most frequently used operator.


.. _Reconstruction:

Reconstruction and limiters
...........................

A reconstruction of discrete fields is used to increase accuracy of discrete models.
The reconstruction can be either unconstrained or limited. 
Amanzi supports a variety of state-of-the-art reconstruction and limiting algorithms 
and their extensions for various PKs.

* `"reconstruction`" [list] describes parameters used by reconstruction algorithms.

 * `"method`" [string] specifies a reconstruction method. Available option is
   `"cell-based`" (default).

 * `"polynomial order`" [int] defines the polynomial order of the reconstructed function. 
   Default is 1.

 * `"limiter`" [string] specifies limiting method. Available options are 
   `"Barth-Jespersen`" (default), `"Michalak-Gooch`", `"tensorial`", and `"Kuzmin`". 

 * `"limiter stencil`" [string] specifies stencil for calculating local bounds. Available 
   options are `"face to cells`", `"cell to closets cells`", `"cell to all cells`",
   and `"node to cells`".
   For a square mesh, the above options define stencils of size 2, 5, 9, and 4,
   respectively.
   Option `"face to cells`" is default for `"Barth-Jespersen`", `"Michalak-Gooch`", 
   and `"tensorial`".  Option `"node to cells`" is default for `"Kuzmin`".

 * `"limiter points`" [int] specifies the number of integration points (Gauss points in 2D) 
   on face where limiting occurs. Default is 1. Limited to 2D.

 * `"limiter location`" [string] defines geometry entity where the *limiter points*
   are located. Available options are `"node`", `"face`", and `"cell`".
   Option `"node`" is default for `"node to cells`" stencil.
   Option `"face`" is default for other stencils.

 * `"limiter cfl`" [double] is a safety factor (less than 1) applied to the limiter.
   Default value is 1.

 * `"use external bounds`" [bool] specifies if bounds for limiters are provided by 
   the hosting application. Default is `"false`".`

 * `"limiter extension for transport`" [bool] adds additional corrections to 
   limiters required by the transport PK. Default value is *false*.

.. code-block:: xml

  <ParameterList name="reconstruction">
    <Parameter name="method" type="string" value="cell-based"/>
    <Parameter name="polynomial order" type="int" value="1"/>
    <Parameter name="limiter" type="string" value="tensorial"/>
    <Parameter name="limiter extension for transport" type="bool" value="false"/>
    <Parameter name="limiter stencil" type="string" value="face to cells"/>
    <Parameter name="limiter points" type="int" value="0"/>
  </ParameterList>


Time integrator
---------------

There exists only one time integrator, *BDF1*.
The standard parameters for this time integrator include sublists of control
options for the timestep controller, preconditioner, and nonlinear solver.
Each PK may have additional parameters.

* `"max preconditioner lag iterations`" [int] specifies frequency of 
  preconditioner recalculation.

* `"freeze preconditioner`" [bool] enforces preconditioner to be updated only
  once per non-linear solver. When set to *true*, the above parameter is ignored.
  Default value is *false*.

* `"extrapolate initial guess`" [bool] identifies forward time extrapolation
  of the initial guess. Default is `"true`".

* `"nonlinear iteration initial guess extrapolation order`" [int] defines
  extrapolation algorithm. Defualt is 1. Zero value implies no extrapolation.

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

* `"timestep controller type`" [string] defines one of a few available controllers.
  This parameter typically requires a sublist with controller's parameters, e.g.
  `"timestep controller standard parameters`".

* `"solver type`" [string] defines nonlinear solver used on each time step for
  a nonlinear algebraic system :math:`F(x) = 0`. 
  The available options `"aa`", `"nka`" and `"Newton`".
  This parameter typically requires a sublist with solver's control parameters, e.g. 
  `"nka parameters`".

* `"residual debugger`" [list] a residual debugger specification.

  * `"file name base`" [string] specifies file name with the debug data. Default is `"amanzi_dbg`".
    
.. code-block:: xml

  <ParameterList name="time integrator">
    <Parameter name="time integration method" type="string" value="BDF1"/>
    <ParameterList name="BDF1">
      <Parameter name="max preconditioner lag iterations" type="int" value="5"/>
      <Parameter name="freeze preconditioner" type="bool" value="false"/>
      <Parameter name="extrapolate initial guess" type="bool" value="true"/>
      <Parameter name="nonlinear iteration initial guess extrapolation order" type="int" value="1"/>
      <Parameter name="restart tolerance relaxation factor" type="double" value="1.0"/>
      <Parameter name="restart tolerance relaxation factor damping" type="double" value="1.0"/>

      <Parameter name="timestep controller type" type="string" value="standard"/>
      <ParameterList name="timestep controller standard parameters">
        ...
      </ParameterList>

      <Parameter name="solver type" type="string" value="nka"/>
      <ParameterList name="nka parameters">
        ... 
      </ParameterList>
    </ParameterList>
  </ParameterList>


.. _TimeStepController:

Time step controller
--------------------

The time step is controlled by parameter *time step controller type*
and the related list of options.
Nonlinear solver is controlled by parameter *solver type*  and related list of options.
Amanzi supports a few nonlinear solvers described in details in a separate section.

The time step controller *standard* is a simple timestep control mechanism
which sets the next timestep based upon the previous timestep and how many
nonlinear iterations the previous timestep took to converge.
The next time step is given by the following rule:

* if :math:`N_k > N^{max}` then :math:`\Delta t_{k+1} = f_{reduction} \Delta t_{k}`
* if :math:`N_k < N^{min}` then :math:`\Delta t_{k+1} = f_{increase} \Delta t_{k}`
* otherwise :math:`\Delta t_{k+1} = \Delta t_{k}`

where :math:`\Delta t_{k}` is the previous timestep and :math:`N_k` is the number of nonlinear 
iterations required to solve step :math:`k`.

The time step controller *smart* is based on *standard*, but also tries to be a bit 
smarter to avoid repeated increase/decrease loops where the step size decreases, 
converges in few iterations, increases, but then fails again.  It also tries to grow 
the time step geometrically to more quickly recover from tricky nonlinearities.

The time step controller *from file* loads a timestep history from a file, then
advances the step size with those values.  This is mostly used for testing
purposes, where we need to force the same timestep history as previous runs to
do regression testing.  Otherwise roundoff errors can eventually alter
number of iterations enough to alter the timestep history, resulting in
solutions which are enough different to cause doubt over their correctness.

* `"time step controller type`" [list]
  Available options are `"fixed`", `"standard`", `"smarter`",  `"adaptive`", and `"from file`"
  The later is under development and is based on a posteriori error estimates.

  * `"time step increase factor`" [double] defines geometric grow rate for the
    initial time step. This factor is applied when nonlinear solver converged
    in less than `"min iterations`" iterations. 
    This value can be modified geometrically by the smart controller in the 
    case of repeated successful steps. Default is 1.

  * `"min iterations`" [int]  triggers increase of the time step if the previous step 
    took less than this.

  * `"time step reduction factor`" [double] defines abrupt time step reduction
    when nonlinear solver failed or did not converge in `"max iterations`" iterations.

  * `"max iterations`" [int] itriggers decrease of the time step if the previous step 
    took more than this.

  * `"max time step`" [double] is the maximum allowed time step.

  * `"min time step`" [double] is the minimum allowed time step.

.. code-block:: xml

  <ParameterList name="BDF1"> <!-- parent list -->
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

In this example, the time step is increased by factor 1.2 when the nonlinear
solver converges in 10 or less iterations. 
The time step is not changed when the number of nonlinear iterations is
between 11 and 15.
The time step will be cut twice if the number of nonlinear iterations exceeds 15.

Parameters for other controllers are:

* `"file name`" [string] is the path to hdf5 file containing timestep information. 
  The parameter is used by only one controller.

* `"timestep header`" [string] is the name of the dataset containing the history 
  of timestep sizes. The parameter is used by only one controller.

* `"max time step increase factor`" [double] specifies the maximum value for 
  parameter `"time step increase factor`" in the smart controller. Default is 10.

* `"growth wait after fail`" [int] defined the number of skipped timesteps before
  attempting to grow the timestep after a failed timestep. 
  This parameter is used by the smart controller only.

* `"count before increasing increase factor`" [int] defines the number of successive 
  increasions before multiplying parameter `"time step increase factor`". 
  This parameter is used by the smart controller only.
 

.. _Functions:

Functions
---------

To set up non-trivial boundary conditions and/or initial fields, Amanzi
supports a few mathematical functions. 
New function types can added easily.
Each function is defined by a list:

.. code-block:: xml

  <ParameterList name="function name">
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
    <Parameter name="forms" type="Array(string)" value="{linear, constant, _USER_FUNC}"/>

    <ParameterList name="_USER_FUNC">
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

It is possible to populate *x values* and *y values* parameters in a tabular 
function from an HDF5 file. The parameter *forms* is optional.

* `"file`" [string] is the path to hdf5 file containing tabulr function information. 

* `"x header`" [string] name of the dataset containing *x values*.

* `"y header`" [string] name of the dataset containing *y values*.

.. code-block:: xml

  <ParameterList name="function-tabular">
    <Parameter name="file" type="string" value="surface.h5"/>
    <Parameter name="x header" type="string" value="times"/>
    <Parameter name="x coordinate" type="string" value="t"/>
    <Parameter name="y header" type="string" value="recharge"/>
    <Parameter name="forms" type="Array(string)" value="{linear, constant, constant}"/>
  </ParameterList>
  

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
    <Parameter name="file" type="string" value="pressure.h5"/>
    <Parameter name="row header" type="string" value="/time"/>
    <Parameter name="row coordinate" type="string" value="time"/>
    <Parameter name="column header" type="string" value="/x"/>
    <Parameter name="column coordinate" type="string" value="x"/>
    <Parameter name="value header" type="string" value="/pressures"/>
  </ParameterList>


Bilinear-and-time function
..........................

The bilinear-and-time function does trilinear interpolation between
two spatial dimensions and one temporal dimension, but does so with
lazy loading of an HDF5 file for the temporal dimension.  This allows,
for instance, efficient interpolation of precipitation data that is
both temporally varying and provided as a raster.  By lazily loading
the coefficients, long time series can be safely interpolated.

.. code-block:: xml

  <ParameterList name="function-bilinear-and-time">
    <Parameter name="file" type="string" value="precipitation_rain.h5"/>
    <Parameter name="time header" type="string" value="/time"/>
    <Parameter name="x header" type="string" value="/x"/>
    <Parameter name="y header" type="string" value="/y"/>
    <Parameter name="value header" type="string" value="/values"/>
  </ParameterList>

Note this expects and HDF5 file laid out as:

.. code-block:: text

   precipitation_rain.h5
   | time (NTIMES 1D array)
   | x (NX 1D array)
   | y (NY 1D array)
   | values (group)
   |  | 0 (NX x NY 2D array)
   |  | 1 (NX x NY 2D array)
   |  | ...
   |  | NTIMES (NX x NY 2D array)


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


Distance function
..................

A distance function calculates distance from reference point :math:`x_0`
using by the following expression:

.. math::
  f(x) = \sum_{j=0}^{n} m_j (x_j - x_{0,j})^2

Note that the first parameter in :math:`x` can be time.
Here is an example of a distance function using isotropic metric:

.. code-block:: xml

  <ParameterList name="function-distance">
    <Parameter name="x0" type="Array(double)" value="{1.0, 3.0, 0.0}"/>
    <Parameter name="metric" type="Array(double)" value="{1.0, 1.0, 1.0}"/>
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
  

Polynomial function
...................

A generic polynomial function of one argument is given by the following expression:

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
  

Multi-variable monomial function
................................

A multi-variable monomial function is given by the following expression:

.. math::
  f(x) = c \prod_{j=0}^{n} (x_j - x_{0,j})^{p_j}

with the constant factor :math:`c`, the reference point :math:`x_0`, and
integer exponents :math:`p_j`. 
Note that the first parameter in :math:`x` can be time.
Here is an example of monomial of degree 6 in three variables:

.. code-block:: xml

  <ParameterList name="function-monomial">
    <Parameter name="c" type="double" value="1.0"/>
    <Parameter name="x0" type="Array(double)" value="{1.0, 3.0, 0.0}"/>
    <Parameter name="exponents" type="Array(int)" value="{2, 3, 1}"/>
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
  Default value is 0.

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


Static head function
....................

This function requires a few parameters as well as one of the standard math functions.

.. code-block:: xml

  <ParameterList name="function-static-head">
    <Parameter name="p0" type="double" value="101325.0"/>
    <Parameter name="density" type="double" value="1000.0"/>
    <Parameter name="gravity" type="double" value="9.8"/>
    <Parameter name="space dimension" type="int" value="3"/>
    <ParameterList name="water table elevation">
      <ParameterList name="function-constant">
        <Parameter name="value" type="double" value="1.0"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>


Expression function
...................

The expression function is defined by a string expression.
The function has min(N, D + 1) arguments t, x, y, and z. The argument t is required. 
D is the space dimension, and N is the user specified number of arguments which could 
be less than D + 1.

Example of a quadratic function in 2D:

.. code-block:: xml

  <ParameterList name="function-exprtk">
    <Parameter name="number of arguments" type="int" value="3"/>
    <Parameter name="formula" type="string" value="t + x + y^2"/>
  </ParameterList>


Time functions
--------------

Boundary condition functions utilize a parameterized model for time variations that is either piecewise constant or piecewise linear.  For example:

.. code-block:: xml

  <Parameter name="times" type="Array(double)" value="{1, 2, 3}"/>
  <Parameter name="time values" type="Array(double)" value="{10, 20, 30}"/>
  <Parameter name="time functions" type="Array(string)" value="{constant, linear}"/>    

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

Iterative solvers
-----------------

This list contains sublists for various linear algebra solvers such as PCG, GMRES, and NKA.
Note that only PK can provide a preconditioner for a linear solver; hence, we cannot
specify it here.

* `"iterative method`" [string] defines a Krylov-based method. The available options
  include `"pcg`" and `"gmres`".

* `"direct method`" [string] defines a direct method. The available option is `"amesos`".

* `"xxx parameters`" [list] provides parameters for the iterative method specified 
  by variable `"iterative method`".
 
.. code-block:: xml

  <ParameterList>  <!-- parent list -->
  <ParameterList name="solvers">
    <ParameterList name="_GMRES WITH HYPRE AMG">
      <Parameter name="iterative method" type="string" value="gmres"/>

      <ParameterList name="gmres parameters">
        ...
      </ParameterList>
    </ParameterList>

    <ParameterList name="_PCG with HYPRE AMG">
      <Parameter name="iterative method" type="string" value="pcg"/>
      <ParameterList name="pcg parameters">
        ...
      </ParameterList>
    </ParameterList>
  </ParameterList>
  </ParameterList>



Generalized minimal residuals (GMRES)
.....................................

Not all scientists know that idea of GMRES method was formulated first in 1968.
We support two implementations (Amanzi and Tilinos). 
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

* `"preconditioning strategy`" [string] defines either `"left`" or `"right`" preconditioner.
  Default is `"left`".

* `"controller training start`" [int] defines the iteration number when the stagnation controller 
  starts to collect data of the convergence history. Default is 0.

* `"controller training end`" [int] defines the iteration number when the stagnation controller
  stops to collect data of the convergence history. The cotroller becomes active on the next
  iteration. Default is 3.

* `"maximum size of deflation space`" [int] defines the size of deflation space. It should be 
  smaller than the size of the Krylov space. Default is 0. This is experimental feature.

* `"release Krylov vectors`" [bool] cleans the stack of Krylov vectors. Default is *false*.
  This may be useful to reduce internal work memory of a persistent solver.

.. code-block:: xml

  <ParameterList name="_GMRES with HYPRE AMG">  <!-- parent list -->
  <ParameterList name="gmres parameters">
    <Parameter name="error tolerance" type="double" value="1e-12"/>
    <Parameter name="maximum number of iterations" type="int" value="400"/>
    <Parameter name="convergence criteria" type="Array(string)" value="{relative residual}"/>
    <Parameter name="size of Krylov space" type="int" value="10"/>
    <Parameter name="overflow tolerance" type="double" value="3.0e+50"/>
    <Parameter name="maximum size of deflation space" type="int" value="0"/>
    <Parameter name="preconditioning strategy`" type="string" value="left"/>
    <Parameter name="release Krylov vectors" type="bool" value="false"/>

    <ParameterList name="verbose object">
      <Parameter name="verbosity level" type="string" value="high"/>
    </ParameterList>
  </ParameterList>

  <!-- Alternative implementation
  <ParameterList name="belos gmres parameters">
    <Parameter name="error tolerance" type="double" value="1e-12"/>
    <Parameter name="maximum number of iterations" type="int" value="400"/>
    <Parameter name="convergence criteria" type="Array(string)" value="{relative residual}"/>
    <Parameter name="size of Krylov space" type="int" value="10"/>
    <Parameter name="overflow tolerance" type="double" value="3.0e+50"/>
  </ParameterList-->
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

  <ParameterList name="_PCG with HYPRE AMG">  <!-- parent list -->
  <ParameterList name="pcg parameters">
    <Parameter name="error tolerance" type="double" value="1e-12"/>
    <Parameter name="maximum number of iterations" type="int" value="400"/>
    <Parameter name="convergence criteria" type="Array(string)" value="{relative residual,make one iteration}"/>
    <Parameter name="overflow tolerance" type="double" value="3.0e+50"/>

    <ParameterList name="verbose object">
      <Parameter name="verbosity level" type="string" value="high"/>
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

  <ParameterList name="_NKA">  <!-- parent list -->
  <ParameterList name="nka parameters">
    <Parameter name="error tolerance" type="double" value="1e-12"/>
    <Parameter name="maximum number of iterations" type="int" value="400"/>
    <Parameter name="convergence criteria" type="Array(string)" value="{relative residual}"/>
    <Parameter name="overflow tolerance" type="double" value="3.0e+50"/>
    <Parameter name="max nka vectors" type="int" value="10"/>
    <Parameter name="nka vector tolerance" type="double" value="0.05"/>

    <ParameterList name="verbose object">
      <Parameter name="verbosity level" type="string" value="high"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>


Direct solvers from Amesos library 
..................................

Amesos library of Trilinos package provides interfaces to a few direct solvers.
List `"amesos parameters`" contains parameters that understood by this library.
These parameters may violate the camel-case convention employed by this spec.
Additional parameters are:

* `"solver name`" [string] declares name of one of the supported direct solvers. 
  Available options are `"klu`", `"superludist`", `"basker`", etc, see Amesos and 
  Amesos2 manuals for details. The default value is serial solver `"klu`".

* `"amesos version`" [int] specifies version of Amesos. Available options are 1 and 2.
  The default value is 1.

.. code-block:: xml

  <ParameterList name="_AMESOS KLU">  <!-- parent list -->
  <Parameter name="direct method" type="string" value="amesos"/>
  <ParameterList name="amesos parameters">
    <Parameter name="solver name" type="string" value="klu"/>
    <Parameter name="amesos version" type="int" value="1"/>
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
  options are `"monitor update`" (default), `"monitor residual`", 
  `"monitor preconditioned residual`", `"monitor l2 residual`", and 
  `"monitor preconditioned l2 residual`".

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

  <ParameterList name="BDF1">  <!-- typical parent list -->
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

    <ParameterList name="verbose object">
      <Parameter name="verbosity level" type="string" value="high"/>
    </ParameterList>
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

  <ParameterList name="BDF1">  <!-- typical parent list -->
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

  <ParameterList name="BDF1">  <!-- typical parent list -->
  <Parameter name="solver type" type="string" value="Newton"/>
  <ParameterList name="Newton parameters">
    <Parameter name="nonlinear tolerance" type="double" value="1.0e-05"/>
    <Parameter name="diverged tolerance" type="double" value="1.0e+10"/>
    <Parameter name="max du growth factor" type="double" value="1.0e+03"/>
    <Parameter name="max divergent iterations" type="int" value="3"/>
    <Parameter name="limit iterations" type="int" value="20"/>
    <Parameter name="modify correction" type="bool" value="true"/>
  </ParameterList>
  </ParameterList>


Inexact Newton
..............

The inexact Newton methods work for cases where the discrete Jacobian is either 
*not* available, or not stable, or computationally expensive. The discrete
Jacobian is replaced by a stable approximation of the continuum Jacobian.
This solver has the same list of parameters as the Newton solver. 

The difference between these solvers is in the preconditioner parameters.
Here is the list of selected parameters for the Newton-Picard solver.

.. code-block:: xml

  <ParameterList name="operators">
    <ParameterList name="diffusion operator">
      <ParameterList name="preconditioner">
        <Parameter name="discretization primary" type="string" value="mfd: optimized for monotonicity"/>
        <Parameter name="discretization secondary" type="string" value="mfd: optimized for sparsity"/>
        <Parameter name="schema" type="Array(string)" value="{face, cell}"/>
        <Parameter name="preconditioner schema" type="Array(string)" value="{face, cell}"/>
        <Parameter name="Newton correction" type="string" value="approximate Jacobian"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>
       
  <Parameter name="solver type" type="string" value="Newton"/>
  <ParameterList name="Newton parameters">
    <Parameter name="nonlinear tolerance" type="double" value="1.0e-05"/>
    <Parameter name="diverged tolerance" type="double" value="1.0e+10"/>
    <Parameter name="max du growth factor" type="double" value="1.0e+03"/>
    <Parameter name="max divergent iterations" type="int" value="3"/>
    <Parameter name="limit iterations" type="int" value="20"/>
    <Parameter name="modify correction" type="bool" value="false"/>
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
      <Parameter name="nonlinear tolerance" type="double" value="1.0e-05"/>
      <Parameter name="diverged tolerance" type="double" value="1.0e+10"/>
      <Parameter name="limit iterations" type="int" value="20"/>
      <Parameter name="max divergent iterations" type="int" value="3"/>
    </ParameterList>

    <ParameterList name="linear operator">
      <Parameter name="iterative method" type="string" value="gmres"/>
      <ParameterList name="gmres parameters">
        ...
      </ParameterList>
    </ParameterList>
  </ParameterList>
  </ParameterList>


Nonlinear Object-Oriented Solution (NOX)
........................................

The interface to Trilinos NOX solver is as follows:

.. code-block:: xml

  <Parameter name="solver type" type="string" value="nox"/>
    <ParameterList name="nox parameters">
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
preconditioner, Hypre BoomerAMG preconditioner, ILU preconditioner, Hypre's Euclid ILU
preconditioner, and identity preconditioner. 

* `"preconditioning method`" [string] defines preconditioner algorithm.

* `"xxx parameters`" [list] provides parameters for the preconditioner specified 
  by parameter `"preconditioning method`".
 
.. code-block:: xml

  <ParameterList>  <!-- parent list -->
  <ParameterList name="preconditioners">
    <ParameterList name="_TRILINOS ML">
      <Parameter name="preconditioning method" type="string" value="ml"/>
      <ParameterList name="ml parameters">
        ... 
      </ParameterList>
    </ParameterList>

    <ParameterList name="_HYPRE AMG">
      <Parameter name="preconditioning method" type="string" value="boomer amg"/>
      <ParameterList name="boomer amg parameters">
        ...
      </ParameterList>
    </ParameterList>

    <ParameterList name="_BLOCK ILU">
      <Parameter name="preconditioning method" type="string" value="block ilu"/>
      <ParameterList name="block ilu parameters">
        ...
      </ParameterList>
    </ParameterList>

    <ParameterList name="_DIAGONAL">
      <Parameter name="preconditioning method" type="string" value="diagonal"/>
    </ParameterList>
  </ParameterList>


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

* `"max multigrid levels`" [int] defines the maximum number of multigrid levels.

* `"strong threshold`" [double] sets AMG strength threshold. The default is 0.25. 
  For 2D Laplace operators, 0.25 is a good value, for 3D Laplace operators, 0.5 or 
  0.6 is a better value. For elasticity problems, a large strength threshold,
  such as 0.9, is often better.

* `"coarsen type`" [int] defines which parallel coarsening algorithm is used. 
  The following options for coarsen type:

  * 0 - (default) CLJP-coarsening, a parallel coarsening algorithm using independent sets.
  * 3 - classical Ruge-Stueben coarsening on each processor, followed by a third pass, which adds coarse points on the boundaries.
  * 8 - PMIS-coarsening, a parallel coarsening algorithm using independent sets, generating lower complexities than CLJP, might also lead to slower convergence.
  * 21 - CGC coarsening by M. Griebel, B. Metsch and A. Schweitzer.
  * 22 - CGC-E coarsening by M. Griebel, B. Metsch and A.Schweitzer.

* `"max coarse size`" [int] sets maximum size of coarsest grid. The default is 9.

* `"use block indices`" [bool] If true, uses the `"systems of PDEs`" code with blocks 
  given by the SuperMap, or one per DoF per entity type. Default is *false*.
  Note Hypre's BoomerAMG cannot be given both `"use block indices`" and 
  `"number of functions`" (see later) options as these are two ways of specifying 
  the same thing.

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

  <ParameterList name="_HYPRE AMG">  <!-- parent list -->
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


Hypre's Euclid ILU
..................

The Euclid Parallel ILU algorithm was presented at SC99 and published in expanded 
form in the SIAM Journal on Scientific Computing. 
Scalability means that the factorization (setup) and application (triangular solve) timings remain
nearly constant when the global problem size is scaled in proportion to the number of processors.
As with all ILU preconditioning methods, the number of iterations is expected to increase with
global problem size.
Internal parameters for this preconditioner include

* `"ilu(k) fill level`" [int] is the factorization level. Default is 1.

* `"ilut drop tolerance`" defines a drop tolerance relative to the largest 
  absolute value of any entry in the row being factored.

* `"rescale row`" [bool] if true, values are scaled prior to factorization 
  so that largest value in any row is +1 or -1. Note that this can destroy 
  matrix symmetry. 

* `"verbosity`" [int] prints a summary of runtime settings and timing 
  information to stdout.  Default is 0.

.. code-block:: xml

  <ParameterList name="_MY EUCLID">  <!-- parent list -->
  <ParameterList name="euclid parameters">
    <Parameter name="ilu(k) fill level" type="int" value="6"/>
    <Parameter name="ilut drop tolerance" type="double" value="0.01"/>
    <Parameter name="rescale rows" type="bool" value="true"/>
    <Parameter name="verbosity" type="int" value="0"/>
  </ParameterList>
  </ParameterList>


Hypre's Block ILU
..................

The internal parameters for block ILU are as follows (see 
https://docs.trilinos.org/dev/packages/ifpack/doc/html/index.html for more detail):

* `"fact: level-of-fill`" [int] sets the level of fill used by the level-based ilu(k) strategy.
  This is based on powers of the graph, so the value 0 means no-fill. Default is 0.

* `"fact: absolute threshold`" [double] defines the value to add to each diagonal element 
  (times the sign of the actual diagonal element). Default is 0.

* `"fact: relative threshold`" [double] multiplies the diagonal by this value before checking 
  the threshold. Default is 1.

* `"fact: relax value`" [double] if nonzero, dropped values are added to the diagonal (times 
  this factor). Default is 0.

* `"overlap`" [int] defines overlap of the additive Schwarz. Default is 0.

* `"schwarz: combine mode`" [string] defines how values corresponding to overlapping nodes are 
  handled. Users should set this value to `"Add`" if interested in a symmetric preconditioner. 
  Otherwise, the default value of `"Zero`" usually results in better convergence.
  If the level of overlap is set to zero, the rows of the user matrix that are stored on a 
  given processor are treated as a self-contained local matrix and all column entries that 
  reach to off-processor entries are ignored. Setting the level of overlap to one tells Ifpack 
  to increase the size of the local matrix by adding rows that are reached to by rows owned by 
  this processor. Increasing levels of overlap are defined recursively in the same way. 
  For sufficiently large levels of overlap, the entire matrix would be part of each processor's 
  local ILU factorization process. 
  Default is `"Add`".

.. code-block:: xml

  <ParameterList name="_MY ILU">  <!-- parent list -->
  <ParameterList name="block ilu parameters">
    <Parameter name="fact: relax value" type="double" value="1.0"/>
    <Parameter name="fact: absolute threshold" type="double" value="0.0"/>
    <Parameter name="fact: relative threshold" type="double" value="1.0"/>
    <Parameter name="fact: level-of-fill" type="int" value="10"/>
    <Parameter name="overlap" type="int" value="0"/>
    <Parameter name="schwarz: combine mode" type="string" value="Add"/>
  </ParameterList>
  </ParameterList>


Trilinos ML
...........

Internal parameters for Trilinos ML include

.. code-block:: xml

  <ParameterList name="_MY ML">  <!-- parent list -->
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


Identity
........

The identity preconditioner is instantiated if either no preconditioner is
specified or the specified preconditioner list does not exists.


Diagonal
........

The diagonal preconditioner is instantiated by the name and needs no additional
parameters.


Mesh
====

Amanzi supports both structured and unstructured numerical solution approaches.
This flexibility has a direct impact on the selection and design of the underlying 
numerical algorithms, the style of the software implementations, and, ultimately, 
the complexity of the user-interface.  
This specification format uses and describes the unstructured mesh only.

* `"mesh`" [list] accepts `"unstructured`" to indicate the meshing option that Amanzi will use.
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

  * `"unstructured`" [list] accepts instructions to either (1) read or, (2) generate an unstructured mesh.

    * `"read mesh file`" [list] accepts name, format of pre-generated mesh file

      * `"file`" [string] name of pre-generated mesh file. Note that in the case of an
        Exodus II mesh file, the suffix of the serial mesh file must be .exo and 
        the suffix of the parallel mesh file must be .par.
        When running in serial the code will read this the indicated file directly.
        When running in parallel and the suffix is .par, the code will instead read
        the partitioned files, that have been generated with a Nemesis tool and
        named as filename.par.N.r where N is the number of processors and r is the rank.
        When running in parallel and the suffix is .exo, the code will partition automatically
        the serial file.
     
      * `"format`" [string] format of pre-generated mesh file (`"MSTK`", `"MOAB`", or `"Exodus II`")

    * `"generate mesh`" [list] accepts parameters of generated mesh

      * `"domain low coordinate`" [Array(double)] Location of low corner of domain
      * `"domain high coordinate`" [Array(double)] Location of high corner of domain
      * `"number of cells`" [Array(int)] the number of uniform cells in each coordinate direction

    * `"expert`" [list] accepts parameters that control which particular mesh framework is to be used.

      * `"framework`" [string] one of `"stk::mesh`", `"MSTK`", `"MOAB`" or `"Simple`". 
      * `"verify mesh`" [bool] true or false. 

      * `"partitioner`" [string] defines the partitioning algorithm for parallel unstructured meshes.
        The available options are `"metis"` (default), `"zoltan_graph"` and `"zoltan_rcb"`. `"metis"`
        and `"zoltan_graph"` perform a graph partitioning of the mesh with no regard to the geometry 
        of the mesh. `"zoltan_rcb"` partitions meshes using Recursive Coordinate Bisection which 
        can lead to better partitioning in meshes that are thin in a particular direction. 
        Additionally, the use of `"zoltan_rcb"` with the MSTK framework triggers an option to 
        detect columns of elements in a mesh and adjust the partitioning such that no column is 
        split over multiple partitions. If no partitioner is specified, the default one is used.

      * `"request edges`" [bool] builds support for mesh edges. Only in 3D.

      * `"contiguous global ids`" [bool] enforces contigoud global ids. Default is *true*.

    * `"submesh`" [list] parameters for extracted meshes

      * `"extraction method`" [string] specifies the extraction method. The only available option
        is `"manifold mesh`". If this parameter is missing, the parent mesh framework is used 
        for submesh extraction..

      * `"regions`" [Array(string)] defines a list of regions for submesh. Parameter 
        `"extraction method`" requires a single name in this list.

Example of *Unstructured* mesh generated internally:

.. code-block:: xml

  <ParameterList>  <!-- parent list -->
  <ParameterList name="mesh">
    <ParameterList name="unstructured"/>
      <ParameterList name="generate mesh">
        <Parameter name="number of cells" type="Array(int)" value="{100, 1, 100}"/>
        <Parameter name="domain low cooordinate" type="Array(double)" value="{0.0, 0.0, 0.0}"/>
        <Parameter name="domain high coordinate" type="Array(double)" value="{103.2, 1.0, 103.2}"/>
      </ParameterList>   

      <ParameterList name="expert">
        <Parameter name="framework" type="string" value="MSTK"/>
        <Parameter name="partitioner" type="string" value="metis"/>
      </ParameterList>
    </ParameterList>   
  </ParameterList>
  </ParameterList>

Example of *Unstructured* mesh read from an external file:

.. code-block:: xml

  <ParameterList name="mesh">  <!-- parent list -->
  <ParameterList name="unstructured">
    <ParameterList name="read mesh file">
      <Parameter name="file" type="string" value="mesh_filename"/>
      <Parameter name="format" type="string" value="Exodus II"/>
    </ParameterList>   
  </ParameterList>   
  </ParameterList>


Geometric model
===============

Domain
------

It is not always possible to extract space dimension from provided data.
Therefore, we require the user to provide simple list *domain* with only 
one parameter *spatial dimension*.

* `"spatial dimension`" [int] defined space dimenstion. The available values are 2 or 3.

.. code-block:: xml

  <ParameterList name="domain">
    <Parameter name="spatial dimension" type="int" value="2"/>
  </ParameterList>


Regions
-------

Regions are geometrical constructs used in Amanzi to define subsets of the computational domain in order to specify the problem
to be solved, and the output desired. Regions may represents zero-, one-, two- or three-dimensional subsets of physical space.
for a three-dimensional problem, the simulation domain will be a three-dimensional region bounded by a set of two-dimensional 
regions.  If the simulation domain is N-dimensional, the boundary conditions must be specified over a set of regions are (N-1)-dimensional.

Amanzi automatically defines the special region labeled *All*, which is the 
entire simulation domain. Currently, the unstructured framework does
not support the *All* region, but it is expected to do so in the
near future.

User-defined regions are constructed using the following syntax

 * `"regions`" [list] can accept a number of lists for named regions (REGION)

   * Shape [list] Geometric model primitive, choose exactly one of the following: 
     `"region: point`", `"region: box`", `"region: plane`", `"region: halfspace`",
     `"region: cylinder`", `"region: labeled set`", `"region: layer`", `"region: surface`",
     `"region: boundary`", `"region: box volume fractions`", `"region: line segment"`, 
     or `"region: all"`.

Amanzi supports parameterized forms for a number of analytic shapes, as well as more complex 
definitions based on triangulated surface files.  


Point
.....

List *region: point* defines a point in space. 
Using this definition, cell sets encompassing this point are retrieved inside Amanzi.

* `"coordinate`" [Array(double)] Location of point in space.

.. code-block:: xml

  <ParameterList name="_DOWNWIND"> <!-- parent list -->
    <ParameterList name="region: point">
      <Parameter name="coordinate" type="Array(double)" value="{-150.0, 0.0, 0.0}"/>
    </ParameterList>
  </ParameterList>


Box
...

List *region: box* defines a region bounded by coordinate-aligned
planes. Boxes are allowed to be of zero thickness in only one
direction in which case they are equivalent to planes.

* `"low coordinate`" [Array(double)] Location of the boundary point with the lowest coordinates.

* `"high coordinate`" [Array(double)] Location of the boundary points with the highest coordinates.

.. code-block:: xml

  <ParameterList name="_WELL">  <!-- parent list -->
    <ParameterList name="region: box">
      <Parameter name="low coordinate" type="Array(double)" value="{-5.0,-5.0, -5.0}"/>
      <Parameter name="high coordinate" type="Array(double)" value="{5.0, 5.0,  5.0}"/>
    </ParameterList>
  </ParameterList>


Plane
.....

List *region: plane* defines a plane using a point lying on the plane and normal to the plane.

* `"normal`" [Array(double)] Normal to the plane.

* `"point`" [Array(double)] Point in space.

.. code-block:: xml

  <ParameterList name="_TOP_SECTION"> <!-- parent list -->
    <ParameterList name="region: plane">
      <Parameter name="point" type="Array(double)" value="{2, 3, 5}"/>
      <Parameter name="normal" type="Array(double)" value="{1, 1, 0}"/>
      <ParameterList name="expert parameters">
        <Parameter name="tolerance" type="double" value="1.0e-05"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>


Halfspace
.........

List *region: halfspace* defines a half-space using a point lying on the plane 
and outward normal to the plane. Single expert parameter provides tolerance
for geometric predicates.

* `"normal`" [Array(double)] Normal to the plane.

* `"point`" [Array(double)] Point in space.

.. code-block:: xml

  <ParameterList name="_TOP_SECTION"> <!-- parent list -->
    <ParameterList name="region: halfspace">
      <Parameter name="point" type="Array(double)" value="{2, 3, 5}"/>
      <Parameter name="normal" type="Array(double)" value="{1, 1, 0}"/>
      <ParameterList name="expert parameters">
        <Parameter name="tolerance" type="double" value="1.0e-05"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>


Labeled Set
...........

The list *region: labeled set* defines a named set of mesh entities
existing in an input mesh file. This is the same file that contains
the computational mesh. The name of the entity set is given
by *Label*.  For example, a mesh file in the Exodus II
format can be processed to tag cells, faces and/or nodes with
specific labels, using a variety of external tools. Regions based
on such sets are assigned a user-defined label for Amanzi, which may
or may not correspond to the original label in the exodus file.
Note that the file used to express this labeled set may be in any
Amanzi-supported mesh format (the mesh format is specified in the
parameters for this option).  The *entity* parameter may be
necessary to specify a unique set.  For example, an Exodus file
requires *Cell*, *Face* or *Node* as well as a label (which is
an integer).  The resulting region will have the dimensionality 
associated with the entities in the indicated set. 

* `"label`" [string] Set per label defined in the mesh file.

* `"file`" [string] File name.

* `"format`" [string] Currently, we only support mesh files in the Exodus II format.

* `"entity`" [string] Type of the mesh object (cell, face, etc).

.. code-block:: xml

  <ParameterList name="_AQUIFER">
    <ParameterList name="region: labeled set">
      <Parameter name="entity" type="string" value="cell"/>
      <Parameter name="file" type="string" value="porflow4_4.exo"/>
      <Parameter name="format" type="string" value="Exodus II"/>
      <Parameter name="label" type="string" value="1"/>
    </ParameterList>
  </ParameterList>


Color function
..............

The list *region: color function* defines a region based a specified
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

* `"file`" [string] File name.

* `"value`" [int] Color that defines the set in a tabulated function file.

.. code-block:: xml

  <ParameterList name="SOIL_TOP">
    <ParameterList name="region: color function">
      <Parameter name="file" type="string" value="_GEOLOGY_RESAMP_2D.tf3"/>
      <Parameter name="value" type="int" value="1"/>
    </ParameterList>
  </ParameterList>


Polygon
.......

The list *region: polygon* defines a polygonal region on which mesh faces and
nodes can be queried. NOTE that one cannot ask for cells in a polygonal surface
region. In 2D, the polygonal region is a line and is specified by 2 points.
In 3D, the polygonal region is specified by an arbitrary number of points.
In both cases the point coordinates are given as a linear array. The polygon
can be non-convex.

This provides a set of faces with a normal for computing flux.

The polygonal surface region can be queried for a normal. In 2D, the normal is
defined as [Vy,-Vx] where [Vx,Vy] is the vector from point 1 to point 2.
In 3D, the normal of the polygon is defined by the order in which points 
are specified.

* `"number of points`" [int] Number of polygon points.

* `"points`" [Array(double)] Point coordinates in a linear array. 

.. code-block:: xml

  <ParameterList name="_XY_PENTAGON">
    <ParameterList name="region: polygon">
      <Parameter name="number of points" type="int" value="5"/>
      <Parameter name="points" type="Array(double)" value="{-0.5, -0.5, -0.5, 
                                                             0.5, -0.5, -0.5,
                                                             0.8, 0.0, 0.0,
                                                             0.5,  0.5, 0.5,
                                                            -0.5, 0.5, 0.5}"/>
      <ParameterList name="expert parameters">
        <Parameter name="tolerance" type="double" value="1.0e-3"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>


Cylinder
........

List *region: cylinder* defines a oblique cylinder using the axis of symmetry, 
a point lying on this axis, and radius.

* `"axis`" [Array(double)] Axis of symmetry.

* `"point`" [Array(double)] Point on the axis of symmetry.

* `"radius`" [double] Radius of the cylinder.

.. code-block:: xml

  <ParameterList name="_TOP_SECTION"> <!-- parent list -->
    <ParameterList name="region: halfspace">
      <Parameter name="point" type="Array(double)" value="{2, 3, 5}"/>
      <Parameter name="axis" type="Array(double)" value="{1, 1, 0}"/>
      <Parameter name="radius" type="double" value="1.4"/>
      <ParameterList name="expert parameters">
        <Parameter name="tolerance" type="double" value="1.0e-05"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>


Logical
.......

The list *region: logical* defines logical operations on regions allow for more
advanced region definitions. At this time the Logical Region allows
for logical operations on a list of regions.  In the case of Union
the result is obvious, it is the union of all regions.  Similarly
for Intersection. In the case of Subtraction, subtraction is
performed from the first region in the list.  The Complement is a
special case in that it is the only case that operates on single
region, and returns the complement to it within the domain ENTIRE_DOMAIN.
Currently, multi-region booleans are not supported in the same expression.

* `"operation`" [string] defines operation on the list of regions.
  Availale options are *union*, *intersect*, *subtract*, *complement*

* `"regions`" [Array(string)] specifies the list of involved regions.

.. code-block:: xml

  <ParameterList name="_LOWER_LAYERs">
    <ParameterList name="region: logical">
      <Parameter name="operation" type="string" value="union"/>
      <Parameter name="regions" type="Array(string)" value="{_LAYER1, _LAYER2, _LAYER3}"/>
    </ParameterList>
  </ParameterList>


Boundary
........

List *region: boundary* defines a set of all boundary faces. 
Using this definition, faces located on the domain boundary are extracted.

* `"entity`" [string] Type of the mesh object.  

.. code-block:: xml

  <ParameterList name="_DOMAIN_BOUNDARY"> <!-- parent list -->
    <ParameterList name="region: boundary">
      <Parameter name="entity" type="string" value="face"/>
    </ParameterList>
  </ParameterList>


Enumerated set
..............

List *region: enumerated set* defines a set of mesh entities via the list 
of input global ids..

* `"entity`" ``[string]`` Type of the mesh object.  Valid are *cell*, *face*, *edge*, *node*
* `"entity gids`" ``[Array(int)]`` List of the global IDs of the entities.
  
.. code-block:: xml

  <ParameterList name="_WELL"> <!-- parent list -->
    <ParameterList name="region: enumerated set">
      <Parameter name="entity" type="string" value="face"/>
      <Parameter name="entity gids" type="Array(int)" value="{1, 12, 23, 34}"/>
    </ParameterList>
  </ParameterList>


Box volume fractions
....................

List *region: box volume fraction* defines a region bounded by a box *not* 
aligned with coordinate axes. 
Boxes are allowed to be of zero thickness in only one direction in which case 
they are equivalent to rectangles on a plane or segments on a line.

* `"corner coordinate`" [Array(double)] Location of one box corner.

* `"opposite corner coordinate`" [Array(double)] Location of the opposite box corner.

* `"normals`" [Array(double)] Normals to sides in a linear array. Default is columns of
  the identity matrix. The normals may be scaled arbitrarily but must be orthogonal to
  one another and form the right coordinate frame.

.. code-block:: xml

  <ParameterList name="_BASIN">  <!-- parent list -->
    <ParameterList name="region: box volume fractions">
      <Parameter name="corner coordinate" type="Array(double)" value="{-1.0,-1.0, 1.0}"/>
      <Parameter name="opposite corner coordinate" type="Array(double)" value="{1.0, 1.0, 1.0}"/>
      <Parameter name="normals" type="Array(double)" value="{1.0, 0.0, 0.0
                                                             0.0, 2.0, 0.0,
                                                             0.0, 0.0, 3.0}"/>
    </ParameterList>
  </ParameterList>

This example defines a degenerate box, a square on a surface *z=1*.

Line Segment
............

List *region: line segment* desribes a region defined by a line
segment. This region is a set of cells which intersect with a line
segment.  The line segment is allowed to intersect with one or more cells. Zero length
line segments are allowed. The line segment is defined by its ends
points.

* `"end coordinate"` [Array(double)] Location of one end of a line
  segment.

* `"opposite end coordinate`" [Array(double)] Location of the opposite
  end of a line segment.

.. code-block:: xml

  <ParameterList name="_WELL"> <!-- parent list -->
    <ParameterList name="region: line segment">
      <Parameter name="end coordinate" type="Array(double)" value="{49754.44, 53937.77, 0.0}"/>
      <Parameter name="opposite end coordinate" type="Array(double)" value="{49754.44, 53937.77, 100.0}"/>
    </ParameterList>
  </ParameterList>     


All
...

List *region: all* desribes a region which matches all entities in the
mesh.  No parameters are required.

.. code-block:: xml

  <ParameterList name="_ENTIRE MESH"> <!-- parent list -->
    <ParameterList name="region: all">
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
  <ParameterList name="regions">
    <ParameterList name="_TOP SECTION">
      <ParameterList name="region: box">
        <Parameter name="low coordinate" type="Array(double)" value="{2, 3, 5}"/>
        <Parameter name="high coordinate" type="Array(double)" value="{4, 5, 8}"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="_MIDDLE SECTION">
      <ParameterList name="region: box">
        <Parameter name="low coordinate" type="Array(double)" value="{2, 3, 3}"/>
        <Parameter name="high coordinate" type="Array(double)" value="{4, 5, 5}"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="_BOTTOM SECTION">
      <ParameterList name="region: box">
        <Parameter name="low coordinate" type="Array(double)" value="{2, 3, 0}"/>
        <Parameter name="high coordinate" type="Array(double)" value="{4, 5, 3}"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="_INFLOW SURFACE">
      <ParameterList name="region: labeled set">
        <Parameter name="label"  type="string" value="_SIDESET2"/>
        <Parameter name="file"   type="string" value="_MESH.exo"/>
        <Parameter name="format" type="string" value="Exodus II"/>
        <Parameter name="entity" type="string" value="face"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="_OUTFLOW PLANE">
      <ParameterList name="region: plane">
        <Parameter name="point" type="Array(double)" value="{0.5, 0.5, 0.5}"/>
        <Parameter name="normal" type="Array(double)" value="{0, 0, 1}"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="_BLOODY SAND">
      <ParameterList name="region: color function">
        <Parameter name="file" type="string" value="_FAREA_COLOR.txt"/>
        <Parameter name="value" type="int" value="25"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="_FLUX PLANE">
      <ParameterList name="region: polygon">
        <Parameter name="number of points" type="int" value="5"/>
        <Parameter name="points" type="Array(double)" value="{-0.5, -0.5, -0.5, 
                                                               0.5, -0.5, -0.5,
                                                               0.8, 0.0, 0.0,
                                                               0.5,  0.5, 0.5,
                                                              -0.5, 0.5, 0.5}"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="_ENTIRE MESH">
      <ParameterList name="region: all"\>
    </ParameterList>
  </ParameterList>
  </ParameterList>

In this example, _TOP SECTION, _MIDDLE SECTION and _BOTTOM SECTION
are three box-shaped volumetric regions. _INFLOW SURFACE is a
surface region defined in an Exodus II-formatted labeled set
file and _OUTFLOW PLANE is a planar region. _BLOODY SAND is a volumetric
region defined by the value 25 in color function file.


Output data
===========

Amanzi uses a few ways to communicate simulation data to the user that includes
a short file with observations and full-scale visualization files.

Observation file
----------------

A user may request any number of specific observations from Amanzi.  Each labeled observation data quantity involves a field quantity, a model, a region from which it will extract its source data, and a list of discrete times 
for its evaluation.  The observations are evaluated during the simulation and returned to the calling process through one of Amanzi arguments.

* `"observation data`" [list] can accept multiple lists for named observations.

  * `"observation output filename`" [string] user-defined name for the file that the observations are written to.
    The file name can contain relative or absolute path to an *existing* directory only. 

  * `"time unit`" [string] defines time unit for output data.
    Available options are `"s`", `"h`", `"d`", and `"y`". Default is `"s`".

  * `"mass unit`" [string] defines mass unit for output data. 
    Available options are `"g`", `"lb`", and `"lb`". Default is `"kg`".

  * `"length unit`" [string] defines length unit for output data.
     Available options are `"cm`", `"in`", `"ft`", `"yd`" , `"m`", and `"km`". Default is `"m`".

  * `"concentration unit`" [string] defines concentration unit for output data.
     Available options are `"molar`", and `"SI`". Default is `"molar`".

  * `"precision`" [int] defines the number of significant digits. Default is 16.

  * OBSERVATION [list] user-defined label, can accept values for `"variables`", `"functional`",
    `"region`", `"times`", and TSPS (see below).

    * `"domain name`" [string] name of the domain. Typically, it is either `"domain`" for
      the matrix/subsurface or `"fracture`" for the fracture network.

    * `"variables`" [Array(string)] a list of field quantities taken from the list of 
      available field quantities:

      * volumetric water content [-] (volume water / bulk volume)
      * aqueous saturation [-] (volume water / volume pore space)
      * aqueous pressure [Pa]
      * hydraulic head [m] 
      * permeability-weighted hydraulic head [m] 
      * drawdown [m] 
      * permeability-weighted drawdown [m] 
      * volumetric water content [-]
      * gravimetric water content [-]
      * water table [m]
      * SOLUTE aqueous concentration [mol/m^3]
      * SOLUTE gaseous concentration [mol/m^3]
      * SOLUTE sorbed concentration [mol/kg] 
      * SOLUTE free ion concentration
      * x-, y-, z- aqueous volumetric flux [m/s]
      * material id [-]
      * aqueous mass flow rate [kg/s] (when funtional="integral")
      * aqueous volumetric flow rate [m^3/s] (when functional="integral")
      * fractures aqueous volumetric flow rate [m^3/s] (when functional="integral")
      * SOLUTE volumetric flow rate [mol/s] (when functional="integral")
      * SOLUTE breakthrough curve [mol] (when functional="integral")
      * pH [-] 

    Observations *drawdown* and *permeability-weighted* are calculated with respect to the value 
    registered at the first time it was requested.

    The following observations are point-type obervations: "water table", "drawdown".

    The following observations are integrated continuously in time but saved only at specified
    times: "SOLUTE breakthrough curve". 

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

    * `"interpolation`" [string] the string which defines
      the interpolation method to compute observation. Works ONLY with
      Line Segment region at the moment. Available options `"linear`"
      and `"constant`". Default is `"linear`" 

    * `"weighting`"  [string] the string defined the weighting
      function applied to compute observation. Works ONLY with
      Line Segment region at the moment. Available options `"flux norm`"
      and `"none`". Default is `"none`". `"flux norm`" is the absolute
      value of the Darcy flux in a cell.

The following observation functionals are currently supported.
All of them operate on the variables identified.

* `"observation data: point`" returns the value of the field quantity at a point.

* `"observation data: integral`" returns the integral of the field quantity over the region specified.

* `"observation data: extensive integral`" returns the integral of an extensive variable
  over the region specified.  Note that this should be used over the above Integral when 
  the variable to be integrated is an extensive quantity, i.e. water content or flux.

* `"observation data: minimum`" and `"observation data: maximum`" returns the minimum 
  (respectively maximum) of the field quantity over the region specified.

.. code-block:: xml

  <ParameterList>  <!-- parent list -->
  <ParameterList name="observation data">
    <Parameter name="observation output filename" type="string" value="_OUTPUT.out"/>
    <Parameter name="precision" type="int" value="10"/>
    <ParameterList name="_ANY OBSERVATION NAME">
      <Parameter name="region" type="string" value="_REGION"/>
      <Parameter name="functional" type="string" value="observation data: point"/>
      <Parameter name="variable" type="string" value="volumetric water content"/>

      <Parameter name="cycles" type="Array(int)" value="{100000, 200000, 400000, 500000}"/>
      <Parameter name="cycles start period stop" type="Array(int)" value="{0, 100, -1}"/>

      <Parameter name="times" type="Array(double)" value="{101.0, 303.0, 422.0}"/>
      <Parameter name="times start period stop 0" type="Array(double)" value="{0.0, 10.0, 100.0}"/>
      <Parameter name="times start period stop 1" type="Array(double)" value="{100.0, 25.0, -1.0}"/>
    </ParameterList>

    <ParameterList name="_ANY OBSERVATION NAME B">  <!-- another observation -->
      ...
    </ParameterList>
  </ParameterList>
  </ParameterList>

In this example, we collect `"volumetric water content`" on four selected cycles,
every 100 cycles, three selected times, every 10 seconds from 0 to 100, and every 25 seconds after that.


Checkpoint file
---------------

A user may request periodic dumps of Amanzi checkpoint data.  
The user has no explicit control over the content of these files, but has the guarantee that 
the Amanzi run will be reproducible (with accuracy determined
by machine round errors and randomness due to execution in a parallel computing environment).
Therefore, output controls for *checkpoint data* are limited to file name generation and writing 
frequency, by numerical cycle number.

* `"checkpoint data`" [list] can accept a file name base [string] and cycle data [list] 
  used to generate the file base name or directory base name that is used in writing checkpoint data. 

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
  <ParameterList name="checkpoint data">
    <Parameter name="file name base" type="string" value="_CHKPOINT"/>
    <Parameter name="file name digits" type="int" value="5"/>

    <Parameter name="cycles start period stop" type="Array(int)" value="{0, 100, -1}"/>
    <Parameter name="cycles" type="Array(int)" value="{999, 1001}"/>

    <Parameter name="times start period stop 0" type="Array(double)" value="{0.0, 10.0, 100.0}"/>
    <Parameter name="times start period stop 1" type="Array(double)" value="{100.0, 25.0, -1.0}"/>
    <Parameter name="times" type="Array(double)" value="{101.0, 303.0, 422.0}"/>
  </ParameterList>
  </ParameterList>

In this example, checkpoint data files are written when the cycle number is 
a multiple of 100.


Visualization file
------------------

A user may request periodic writes of field data for the purposes of visualization.  
The user will specify explicitly what is to be included in the file at each snapshot.
Visualization files can only be written at intervals corresponding to the numerical 
time step values or intervals corresponding to the cycle number; writes are controlled by time step cycle number.

* `"visualization data`" [list] can accept a file name base [string] and cycle data [list] 
  that is used to generate the file base name or directory base name that is used in writing visualization data.
  It can also accept a set of lists to specify which field quantities to write.
  The file name can contain relative or absolute path to an *existing* directory only. 

  * `"file name base`" [string] ("amanzi_vis")

  * `"file format`" [string] ("XDMF") Amanzi supports two types of
    visualization files.  XDMF is the default and preferred method, but does
    not correctly handle general polyhedra.  Serial, 3D general polyhedral
    support is supported by the "SILO" option.  This will eventually be
    extended to parallel, 2/3D support, but this is not yet implemented.

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
  <ParameterList name="visualization data">
    <Parameter name="file name base" type="string" value="_PLOT"/>
  
    <Parameter name="cycles start period stop" type="Array(int)" value="{0, 100, -1}"/>
    <Parameter name="cycles" type="Array(int)" value="{999, 1001}"/>

    <Parameter name="times start period stop 0" type="Array(double)" value="{0.0, 10.0, 100.0}"/>
    <Parameter name="times start period stop 1" type="Array(double)" value="{100.0, 25.0, -1.0}"/>
    <Parameter name="times" type="Array(double)" value="{101.0, 303.0, 422.0}"/>

    <Parameter name="dynamic mesh" type="bool" value="false"/>

    <ParameterList name="write regions">
      <Parameter name="regions" type="Array(string)" value="{_OBS1, _OBS2, _OBS3}"/>
   A  <Parameter name="wells" type="Array(string)" value="{_OBS4}"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>


Walkabout file
--------------

A user may request periodic dumps of Walkabout data. Output controls for Walkabout data are 
limited to file name generation and writing frequency, by numerical cycle number or time.

* `"walkabout data`" [list] can accept a file name base [string] and cycle data [list] 
  used to generate the file base name or directory base name that is used in writing Walkabout data. 

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

  * `"write regions`" [list] contains three lists of equal size with region names,
    material names, and material ids to write into the output file.

    * `"region names`" [Array(string)] specifies names of regions.
    * `"material names`" [Array(int)] specifies names of materials. 
    * `"material ids`" [Array(int)] specifies material ids. 

.. code-block:: xml

  <ParameterList>  <!-- parent list -->
  <ParameterList name="walkabout data">
    <Parameter name="file name base" type="string" value="_WALKABOUT"/>
    <Parameter name="file name digits" type="int" value="5"/>

    <Parameter name="cycles start period stop" type="Array(int)" value="{0, 100, -1}"/>
    <Parameter name="cycles" type="Array(int)" value="{999, 1001}"/>

    <Parameter name="times start period stop 0" type="Array(double)" value="{0.0, 10.0, 100.0}"/>
    <Parameter name="times start period stop 1" type="Array(double)" value="{100.0, 25.0, -1.0}"/>
    <Parameter name="times" type="Array(double)" value="{101.0, 303.0, 422.0}"/>

    <ParameterList name="write regions">
      <Parameter name="region names" type="Array(string)" value="{_REGION1, _REGION2}"/>
      <Parameter name="material names" type="Array(string)" value="{_MAT1, _MAT2}"/>
      <Parameter name="material ids" type="Array(int)" value="{1000, 2000}"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>

In this example, walkabout data files are written when the cycle number is 
a multiple of 100.


Mesh info
---------

A user may request to dump mesh information. Mesh information includes coordinates of cell centroids
written is the order consistent with all output fields.


* `"filename`"[string] - name of the HDF5 file where coordinates of the centroids are dumped.

.. code-block:: xml
                  
  <ParameterList>  <!-- parent list -->                
  <ParameterList name="mesh info">
    <Parameter name="filename" type="string" value="centroids"/>
  </ParameterList>

  <ParameterList name="mesh info fracture">
    <Parameter name="filename" type="string" value="centroids_fracture"/>
  </ParameterList>
  </ParameterList>


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
  <ParameterList name="analysis">
    <Parameter name="used boundary condition regions" type="Array(string)" value="{_REG1,_REG2}"/>
    <Parameter name="used source and sink regions" type="Array(string)" value="{_REG3,_REG4}"/>
    <Parameter name="used observation regions" type="Array(string)" value="{_REG5}"/>
    <ParameterList name="verbose object">
      <Parameter name="verbosity level" type="string" value="high"/>
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










