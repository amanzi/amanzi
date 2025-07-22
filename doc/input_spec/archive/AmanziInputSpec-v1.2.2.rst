========================================
Amanzi XML Input Specification
========================================

.. contents:: **Table of Contents**


Overview
========

The Amanzi simulator evolves a system of conservation equations for reacting flow in porous media, as detailed in the ASCEM report entitled "Mathematical Formulation Requirements and Specifications for the Process Models`" (hereafter referred to as the 'Model Requirements Document (MRD)'). The purpose of the present document is to specify the data required to execute Amanzi.  This specification should be regarded as a companion to the MRD, and parameterizations of the individual submodels are consistent between Amanzi, the MRD and this document. Where applicable, the relevant sections of the MRD are indicated.


Preliminary Concepts
--------------------

Amanzi solves a set of parameterized models for multiphase flow in porous media.  An Amanzi simulation is specified by providing:

* values for a parameterized PDE-based transport model, including boundary and initial conditions, constitutive laws, and parameterized/phenomenological models for fluid and chemical sources and characterizations of the porous medium,

* parameters controlling the selection of key algorithmic options and output, 

* a description of the (discrete) state of the computational system, including a list of the independent variables and instructions for obtaining or generating the discrete mesh, and a characterization of the (parallel) computing environment.

The primary software interface to Amanzi is a compiled C++ function, and much of the input data required is communicated through a single `Teuchos::ParameterList <http://trilinos.sandia.gov/packages/docs/r7.0/packages/teuchos/doc/html/index.html>`_.
A ParameterList consist of a simple hierarchy of parameters and lists of parameters, and is constructed directly from a similarly structured XML file.  The Amanzi input specification is defined in terms of the XML file format
used to construct a `Teuchos::ParameterList <http://trilinos.sandia.gov/packages/docs/r7.0/packages/teuchos/doc/html/index.html>`_.

In practice, Amanzi is called by a "simulation coordinator" which manages the simulation instructions and orchestrates the flow of data.  A basic simulation coordinator is
provided with the Amanzi source code distribution.  This simple stand-alone coordinator can be used to drive a simple sequence of Amanzi runs, or can serve as a template for user-generated extensions supporting more intricate workflows.  


Model Characterization
~~~~~~~~~~~~~~~~~~~~~~

For each fluid phase in the model system, Amanzi requires the specification of the fluid composition.  If a component exists in multiple phases, a relationship is required to compute its phase distribution as a function of the state of the system.
Equation 2.11 of the MRD governs the conservation and transport of each component, where the volumetric flow rate has been specified via Darcy's law (equation 2.10).  For Darcy flow, the properties of the porous medium must be specified over the entire simulation domain.  All transported phases
require an appropriate set of boundary conditions along the edge of the simulation domain.  Source/sink terms and initial data are provided for each phase component and any solutes they contain.  Additionally, the extent of physical domain
is specified, along with its discrete representation (i.e., the mesh).

Categories of Output Data
~~~~~~~~~~~~~~~~~~~~~~~~~

Output data from Amanzi is currently organized into four specific groups:

* `"Observation Data`": values generated during a run that characterize a particular feature of interest in the solution.  Such data would be used, for example, when quantifying the "response" of Amanzi to a set of targeted parameter variations.  `"Observation data`" is evaluated according to instructions given via the input parameter list, and is assembled into a special data structure that is returned to the calling routine through Amanzi's simulation driver function arguments.

* `"Visualization Data`": field data intended for post-processing, and may include various component subsets of the simulation state data.  Future extensions may provide for including point-wise derived fields and material properties.

* `"Checkpoint Data`": the complete set of information that is required to duplicate a given execution of Amanzi.  Typically `"Checkpoint Data`" is created at periodic intervals during a long run in order to facilitate a repeatable simulation restart capability and archiving procedures. Unlike the first two data groups, there are no user-settable parameters to control the contents of a `"Checkpoint Data`" file - the implementation will determine the quantity, format and precision of this data in order to guarantee reproducibility.


* `"Walkabout Data`": post-processed field data specific to be used as input to the particle tracking software Walkabout.

* `"Log Data`": running commentary on the performance and status of Amanzi's execution.  Typically such data is written to a C++ stream which may be directed to a pipe or file.  The amount and detail of log data is determined by a range of verbosity controls.

Generally, `"Visualization Data`" and `"Checkpoint Data`" consists of high-dimensional field data representing snapshots of the evolving discrete variables.  These are large datasets, relative to the other types, and are most often written to disk in a file format that allows a direct repesentation of the underlying discrete mesh and parallel data distribution.


ParameterList XML
-----------------

The Amanzi input file is an ASCII text XML-formatted file that must be framed at the beginning and end by the following statements:


.. code-block:: xml

  <ParameterList name="Main">

  </ParameterList>

The value in the "name" can be anything ("Main" in this example).  A ParameterList consists of just two types of entries: Parameter and ParameterList.  ParameterLists are labeled with a `"name`" [string], while Parameters have a separate fields for `"name`" [string], `"type`" [string] and `"value`" [TYPE], where "TYPE" can be any of the following: double, int, bool, string, Array(double), Array(int), Array(bool), Array(string).  The value of the parameter is given in quotes (e.g. "2.7e3").  Array data is specified as a single comma-deliminated string bounded by {}'s (e.g. "{2.4, 2.1, 5.7}").

.. code-block:: xml

  <ParameterList name="Sub">
    <Parameter name="CFL" type="double" value="0.9"/>
    <Parameter name="ratio" type="Array(int)" value="{2, 2, 4}"/>
  </ParameterList>

In this example, the sublist "Sub" has a parameter named "CFL" that is a "double" and has the value of 0.9, and a Teuchos::Array<int>
parameter named "ratio" such that ratio[0] = 2. ratio[1]=2 and ratio[2]=4.


Syntax of the Specification
~~~~~~~~~~~~~~~~~~~~~~~~~~~

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



Version
=======

Each input set contains at the top level a string variable `"Amanzi Input Format Version`".  As of the most recent update of this specification, the current version of the Amanzi input is `"1.2.1`".  If the version is unspecified, it is assumed to be earlier than `"0.9.0`".  Release notes documenting the evolving input specification version can be found *here*.

* [SU] "Amanzi Input Format Version" [string] Three part version string

Example:

.. code-block:: xml

  <ParameterList name="Main">
    <Parameter name="Amanzi Input Format Version" type="string" value="1.2.1"/>
  </ParameterList>


PETSc
=====

Amanzi uses PETSc (TODO: insert link here to Amanzi TPL list) for access to specialized solvers for flow and chemistry evolution.  PETSc implements its own flexible approach to runtime selection of solver options; the list of available options changes considerably between library versions.  Input of PETSc-specific  options at runtime is managed via an optional string parameter, which names an auxiliary text file.  By default, no file is included.

* [S] "Petsc Options File" [string] Relative path naming a text file with PETSc parameters included

Example:

.. code-block:: xml

  <ParameterList name="Main">
    <Parameter name="Petsc Options File" type="string" value=".petsc"/>
  </ParameterList>

See the PETSc documentation for details of the format of this file, and the available parameters.  At this time, PETSc is used by the flow solvers supporting structured-grid AMR Amanzi runs, and the pFloTran chemistry process kernal (see usage notes below).

General Description
===================

The `"General Description`" parameter list can be used to provide a brief description of the problem specified in the file.  ANY number of string entries can be provided
with ANY label that may be useful for the user own purposes.  This section is not parsed by Amanzi and thus optional.

* [S] LABEL [string] A descriptive string

Example:

.. code-block:: xml

   <ParameterList name="General Description">
     <Parameter name="Model ID" type="string" value="Transient Richards"/>
     <Parameter name="Model name" type="string" value="BC Cribs PE Template"/>
     <Parameter name="Description" type="string" value="Unsat flow and transport"/>
     <Parameter name="Purpose" type="string" value="Provide input req. for Phase II Demo"/>
     <Parameter name="Creation date" type="string" value="09.25.11 01:28"/>
     <Parameter name="Last modified" type="string" value="09.25.11 01:28"/>
  </ParameterList>
  
Unstructured Amanzi ignores this list.


Execution Control
=================

Amanzi supports both single-phase saturated and variably saturated groundwater flow and solute transport on structured and unstructured grids.  As part of the execution control, the user must specify the process models to be employed for each simulation.  There are currently three process models or modes that need to be defined in the input file (1) flow, (2) transport, and (3) chemistry (chemistry is currently a placeholder).  Additionally, the user must indicate whether a time-accurate or steady solution is requested.

Usage:

* [SU] `"Execution Control`"

 * [SU] `"Flow Model`" [string]: flow process model

  * [SU] `"Off`" [string]: No flow model

  * [SU] `"Richards`" [string]: Single phase, variably saturated flow (assume constant gas pressure)

  * [SU] `"Single Phase`" [string]: Single phase, fully saturated flow

 * [SU] `"Transport Model`" [string]: Transport of phases.  Accepts `"Off`" or `"On`" [string]

 * [SU] `"Chemistry Model`" [string]: Chemical interface and engine for reaction of constituents.

  * [SU] `"Off`" [string]: No chemistry model

  * [SU] `"Amanzi`" [string]: Original Amanzi geochemistry engine, which supports only primary species concentrations in initial and boundary conditions, and source terms. 

  * [SU] `"Alquimia`" [string]: Alquimia interface to a geochemistry engine, supporting geochemical constraints in initial and boundary conditions, and in source terms.

 * [SU] `"Time Integration Mode`" [list]: accepts one of three integration modes:

  * [SU] `"Steady`" [list] - Amanzi is run in steady mode.

   * [SU] `"Start`" [double] (S: Optional) Initial value for psuedo time (used as a continuation parameter) to generate a steady solution.

   * [SU] `"End`" [double]: Time that defines a steady solution.  (stopping criteria may be generalized in future releases).

   * [SU] `"Initial Time Step`" [double]: The initial time step for the steady calculation.

  * [SU] `"Transient`" [list] - A time-accurate evolution is desired

   * [SU] `"Start`" [double] (S: Optional) Start time for integration (if a steady mode exists then this time must equal the steady end time)

   * [SU] `"End`" [double]: End of integration period
   
   * [SU] `"Initial Time Step`" [double] (Optional) The intitial time step for the transient calculation. (see S Note below)


   * [SU] `"Maximum Cycle Number`" [int]: (Optional) The maximum allowed cycle number.

  * [SU] `"Transient with Static Flow`" [list] - The flow field is static so no flow solver is called during time stepping. During initialization the flow 
    field is set in one of two ways: (1) A constant Darcy velocity is specified in the initial condition; 
    (2) Boundary conditions for the flow    (e.g., pressure), along with the initial condition for the pressure field are used to solve for the 
    Darcy velocity.  At present this mode only supports the "Single Phase" flow model.

   * [SU] `"Start`" [double] (S: Optional) Start time for integration (if a steady mode exists then this time must equal the steady end time)

   * [SU] `"End`" [double]: End of integration period
   
   * [SU] `"Initial Time Step`" [double]  (Optional) The intitial time step for the transient calculation. (see S Note below)

   * [SU] `"Maximum Cycle Number`" [int]: (Optional) The maximum allowed cycle number.

  * [SU] `"Initialize To Steady`" [list] - Amanzi is run in steady mode with `"Chemistry Model`" = `"Transport Model`" = `"Off`" until a steady solution is obtained.  Any solutes defined below are ignored.  When the solution is steady, the transport and chemistry models are set to user input and the transient integration mode is employed.  Integration continues forward in time.  Method for detection of a steady solution is specified.

   * [SU] `"Start`" [double]: Initial value for time to generate a steady solution

   * [SU] `"Switch`" [double]: Time when Chemistry Model and Transport Model are set to user specified input and Amanzi switches to time-accurate solution approach.

   * [SU] `"End`" [double]: The end of the time-integration period
    
   * [SU] `"Steady Initial Time Step`" [double]: (Optional) The intitial time step for the steady state initialization calculation.

   * [SU] `"Transient Initial Time Step`" [double]: (Optional)  The intitial time step for the transient calculation after "Switch" time.  (see S Note below)

S Note: If unspecified, Amanzi will compute this value based on numerical stability limitations, scaled by the parameter `"Initial Time Step Multiplier`"

 * [SU] `"Time Period Control`" (Optional)

  * [SU] `"Start Times`" [Array(double)]: List of times at which the current time-integrator will be reinitialized.
  * [SU] `"Initial Time Step`"[Array(double)]: The initial time step for each time period. If unspecified, Amanzi will compute this value based on numerical stability limitations, scaled by the parameter `"Initial Time Step Multiplier`"
  * [SU] `"Maximum Time Step`"[Array(double)]: (Optional) The maximum time step for each time period. 
  * [SU] `"Default Initial Time Step`" [double]: (Optional) set the default initial time step, this is used for time integrator restarts that are required by boundary conditions and sources, but are not specified in this list under Start Times, the default value is 1.0. 

 * [SU] `"Verbosity`" [string]: (default: `"Medium`") Choose one of `"None"`, `"Low"`, `"Medium"`, `"High`", or `"Extreme`".

  * [SU] `"None`": No output is written to run log

  * [SU] `"Low`": Minimal logging output, includes information about time stepsizes attempted, and notification of I/O operations

  * [SU] `"Medium`": Includes summary-level activity of each process kernel

  * [SU] `"High`": Includes numerical performance statistics of each process kernal, and miscellaneous status of primary variables

  * [SU] `"Extreme`": Includes detailed iteration-level convergence properties of process kernal sovlers

 * [SU] `"Numerical Control Parameters`" [list]: Any common control parameters are in the `"Common Controls`" list while detailed control parameters associated with the underlying numerical implementation are seperated into specific lists.
 
 For both unstructured and structured, the following common list of parameters is valid:

  * [SU] `"Common Controls"` [list]: Control parameters associtated with both algorithms.  This section will be filled as the input parameters for both algorithms are brought closer in line with each other.

  If the unstructured option is active, the following list of parameters is valid:

  * [U] `"Unstructured Algorithm"` [list]: Control parameters associtated with the unstructured algorithm.

   * [U] `"Flow Process Kernel`" [list]: Control parameters for the flow methods

     * [U] `"Discretization Method`" [string]: Specifies the spatial discretization 
       method. The Available options are: `"mfd scaled`", `"optimized mfd scaled`"
       (default), `"two point flux approximation`", and `"support operator`".
       The second option is recommended for orthogonal meshes and diagonal absolute permeability.

     * [U] `"Relative Permeability`" [string]: Defines a method for calculating the *upwinded*
       relative permeability. The available options are: `"upwind with gravity`", 
       `"upwind with Darcy flux`" (default), `"cell centered"`, and `"upwind amanzi`"
       (experimental).  The first three calculate the relative permeability on mesh interfaces.

     * [U] `"atmospheric pressure`" [double]: Defines the atmospheric pressure, [Pa].   

     * [U] `"Use Picard`" [bool]: Use the Picard solver to find a good initial guess for the steady state solver. (default: `"false`")

   * [U] `"Transport Process Kernel`" [list]: Control parameters for the transport methods

     * [U] `"Transport Integration Algorithm`" [string]: Accepts `"Explicit First-Order`" or `"Explicit Second-Order`" (default: `"Explicit First-Order`")

     * [U] `"CFL`" [double]: Time step limiter, a number less than 1 with default of 1.

     * [U] `"transport subcycling`" [bool]: Accepts `"true`" or `"false`" which corresponds to transport subcycling on or off, respectively. (default: `"true`") Note that setting this parameter to false does not always preclude transport from subcycling. Since the estimate for the transport time step is based on velocities from the previous time step, the actual time step that transport can take after the current flow step might be different from its initial estimate.

   * [U] `"Chemistry Process Kernel`" [list]: Control parameters for the reactive transport methods

     * [U] `"max chemistry to transport timestep ratio`" [double] when both chemistry and transport process kernels are on, the chemistry time step will be limited such that the ratio of (chemistry time step)/(transport time step) < this parameter. By default this parameter equals 1.0. If this parameter is set for example to 10.0, then we limit the chemistry time step to 10 times what the current transport time step is, such that for each chemistry sub-cycle, there will be at most 10 transport sub cycles. (default: `"1.0`", suggested range: 0.2 ... 10.0)

   * [U] `"Steady-State Implicit Time Integration`" [list] Parameters for BDF1 time integration to reach steady-state

     * [U] `"steady max iterations"` [int] If during the steady state calculation, the number of iterations of the nonlinear solver exceeds this number, the subsequent time step is reduced by the factor specified in `"steady time step reduction factor"`. (default: `"15`", suggested range: 10 ... 20)

     * [U] `"steady min iterations"` [int] If during the steady state calculation, the number of iterations of the nonlinear solver exceeds this number, the subsequent time step is increased by the factor specified in `"steady time step increase factor"`. (default: `"10`", suggested range: 5 ... 15)

     * [U] `"steady limit iterations"` [int] If during the steady state calculation, the number of iterations of the nonlinear solver exceeds this number, the current time step is cut in half and the current time step is repeated. (default: `"20`", suggested range: 20 ... 50)

     * [U] `"steady nonlinear tolerance"` [double] The tolerance for the nonlinear solver during the steady state computation. (default: `"1.0e-5`", suggested range: 1.0e-8 ... 1.0e-6)

     * [U] `"steady nonlinear iteration damping factor"` [double] Damp the nonlinear iteration (fixed point iteration) by this factor, the default is 1.0 (no damping). (default: `"1.0`", suggested range: 0.1 ... 1.0)
 
     * [U] `"steady time step reduction factor"` [double] When time step reduction is necessary during the steady calculation, use this factor. (default: `"0.8`", suggested range: 0.5 ... 0.9)

     * [U] `"steady time step increase factor"` [double] When time step increase is possible during the steady calculation, use this factor. (default: `"1.2`", suggested range: 1.1 ... 2.0)

     * [U] `"steady max time step"` [double] During the steady state solve, the time step is limited to the value specified here. (default: `"1.0e+10`", suggested range: 1.0e+8 ... 1.0e+10)

     * [U] `"steady max preconditioner lag iterations"` [int] During the steady state solve, the preconditioner is lagged this amount of iterations during the nonlinear solve. For example, if a value of 4 is specified here, the preconditioner is updated at the beginning of each nonlinear solve and then in each fourth iteration during each nonlinear solve. To force a preconditioner in each iteration of each nonlinear solve, set this parameter to one (very expensive, but also very robust), and to disable updates of the preconditioner, except at the beginning of each nonlinear solve, set this parameter to a value larger than `"steady limit iterations"`. (default: `"5`", suggested range: 0 ... 10)

     * [U] `"steady max divergent iterations`" [int] the BDF1 time integrator will tolerate one less than that many subsequent divergent nonlinear iterations. if there are `"steady max divergent iterations`" then the time iterator will give up on this time step and will cause the current time step to be cut by 50% and the current time step to be repeated. (default: `"3`", suggested range: 3 ... 8)

     * [U] `"steady nonlinear iteration divergence factor`" [double] If during the nonlinear solve, the inf norm of the nonlinear update is larger by this factor than the inf norm of the update in the prior iteration, we abort the nonlinear solve to protect against a runaway divergent iteration that causes numerical overflow. As a result the current time step will repeated with a smaller delta T. (default: `"1000.0`", suggested range: 100.0 ... 10000.0)

     * [U] `"steady restart tolerance relaxation factor`" [double] when the time integrator is started, it may be beneficial to set this parameter to something > 1.0 to loosen the nonlinear tolerance on the first several time steps. The parameter `"steady restart tolerance relaxation factor damping`" controls how fast the this loosened nonlinear tolerance will revert back to the one specified in  `"steady nonlinear tolerance"`: If the nonlinear tolerance is ntol, the initial timestep factor is ntol_factor, and the damping is ntol_damping, then the actual nonlinear tolerance is ntol*ntol_factor, and after every time step, ntol_factor = max(1.0,ntol_factor*ntol_damping), such that a few iterations after a time integrator start, the actual tolerance equals ntol, again. The default for this paramameter is 1.0, while reasonable values are > 1.0, maybe as large as 1000.0. The default for the damping factor is 1.0, while reasonable values are between 0 and 1. (default: `"1.0`", suggested range: 1.0 ... 1000.0)

     * [U] `"steady restart tolerance relaxation factor damping`" [double] see `"steady nonlinear iteration initial timestep factor`" for a detailed explanation of this parameter. (default: `"1.0`", suggested range: 0.001 ... 1.0)

     * [U] `"steady preconditioner`" [string] select the preconditioner to be used in the nonlinear solver for the steady state problem, choose one of `"Trilinos ML`", `"Hypre AMG`", or `"Block ILU`". (default: `"Hypre AMG`")

     * [U] `"steady initialize with darcy`" [bool] Initialize the flow field using a Darcy solve. (default `"true`")  

     * [U] `"steady nonlinear iteration initial guess extrapolation order`" [int] defines how the initial guess (predictor) for a new time step is calculated. If set to zero, the previous solution is used as the initial guess. (default: 1)  

   * [U] `"Transient Implicit Time Integration`" [list] Parameters for BDF1 transient time integration 

     * [U] `"transient max iterations"` [int] If during the transient calculation, the number of iterations of the nonlinear solver exceeds this number, the subsequent time step is reduced by the factor specified in `"transient time step reduction factor"`. (default: `"15`", suggested range: 10 ... 20)

     * [U] `"transient min iterations"` [int] If during the transient calculation, the number of iterations of the nonlinear solver exceeds this number, the subsequent time step is increased by the factor specified in `"transient time step increase factor"`. (default: `"10`", suggested range: 5 ... 15)

     * [U] `"transient limit iterations"` [int] If during the transient calculation, the number of iterations of the nonlinear solver exceeds this number, the current time step is cut in half and the current time step is repeated. (default: `"20`", suggested range: 20 ... 50)

     * [U] `"transient nonlinear tolerance"` [double] The tolerance for the nonlinear solver during the transient computation. (default: `"1.0e-5`", suggested range: 1.0e-6 ... 1.0e-5)

     * [U] `"transient nonlinear iteration damping factor"` [double] Damp the nonlinear iteration (fixed point iteration) by this factor, the default is 1.0 (no damping). (default: `"1.0`", suggested range: 0.1 ... 1.0)

     * [U] `"transient time step reduction factor"` [double] When time step reduction is necessary during the transient calculation, use this factor. (default: `"0.8`", suggested range: 0.5 ... 0.9)

     * [U] `"transient time step increase factor"` [double] When time step increase is possible during the transient calculation, use this factor. (default: `"1.2`", suggested range: 1.1 ... 2.0) Note that this paramter also works in the case where the flow mode `"Single Phase`" was selected. In that case, the default is `"1.0`".

     * [U] `"transient max time step"` [double] During the transient solve, the time step is limited to the value specified here. (default: `"1.0e+8`", suggested range: 1.0e+8 ... 10e+10)

     * [U] `"transient max preconditioner lag iterations"` [int] During the transient solve, the preconditioner is lagged this amount of iterations during the nonlinear solve. For example, if a value of 4 is specified here, the preconditioner is updated at the beginning of each nonlinear solve and then in each fourth iteration during each nonlinear solve. To force a preconditioner in each iteration of each nonlinear solve, set this parameter to one (very expensive, but also very robust), and to disable updates of the preconditioner, except at the beginning of each nonlinear solve, set this parameter to a value larger than `"transient limit iterations"`. (default: `"5`", suggested range: 0 ... 10)

     * [U] `"transient max divergent iterations`" [int] the BDF1 time integrator will tolerate one less than that many subsequent divergent nonlinear iterations. if there are `"transient max divergent iterations`" then the time iterator will give up on this time step and will cause the current time step to be cut by 50% and the current time step to be repeated. (default: `"3`", suggested range: 3 ... 8)

     * [U] `"transient nonlinear iteration divergence factor`" [double] If during the nonlinear solve, the inf norm of the nonlinear update is larger by this factor than the inf norm of the update in the prior iteration, we abort the nonlinear solve to protect against a runaway divergent iteration that causes numerical overflow. As a result the current time step will repeated with a smaller delta T. (default: `"1000.0`", suggested range: 100.0 ... 10000.0)

     * [U] `"transient restart tolerance relaxation factor`" [double] when the time integrator is restarted, at a time when a boundary condition drastically changes, it may be beneficial to set this parameter to something > 1.0 to loosen the nonlinear tolerance on the first several time steps after the time integrator restart. The parameter `"transient restart tolerance relaxation factor damping`" controls how fast the this loosened nonlinear tolerance will revert back to the one specified in  `"transient nonlinear tolerance"`: If the nonlinear tolerance is ntol, the initial timestep factor is ntol_factor, and the damping is ntol_damping, then the actual nonlinear tolerance is ntol*ntol_factor, and after every time step, ntol_factor = max(1.0,ntol_factor*ntol_damping), such that a few iterations after a time integrator restart, the actual tolerance equals ntol, again. The default for this paramameter is 1.0, while reasonable values are > 1.0, maybe as large as 1000.0. The default for the damping factor is 1.0, while reasonable values are between 0 and 1. (default: `"1.0`", suggested range: 1.0 ... 1000.0)

     * [U] `"transient restart tolerance relaxation factor damping`" [double] see `"transient nonlinear iteration initial timestep factor`" for a detailed explanation of this parameter. (default: `"1.0`", suggested range: 0.001 ... 1.0)

     * [U] `"transient preconditioner`" [string] select the preconditioner to be used in the nonlinear solver for the steady state problem, choose one of `"Trilinos ML`", `"Hypre AMG`", or `"Block ILU`". (default: `"Hypre AMG`")

     * [U] `"transient initialize with darcy`" [bool] Initialize the flow field using a Darcy solve. (default `"false`") 

     * [U] `"transient nonlinear iteration initial guess extrapolation order`" [int] defines how the initial guess (predictor) for a new time step is calculated. If set to zero, the previous solution is used as the initial guess. (default: 1)  


   * [U] `"Steady-State Pseudo-Time Implicit Solver`" [list] Parameters for Damped Picard iteration to reach steady-state

     * [U] `"pseudo time integrator initialize with darcy`" [bool] Initialize the pseudo time integrator (Picard) with a Darcy solution. (default: `"true`")

     * [U] `"pseudo time integrator clipping saturation value`" [double] (default: 0.9, suggested range: 0.7 ... 0.95)

     * [U] `"pseudo time integrator time integration method`" [double] select the pseudo time integration method (currrently only Picard is supported). (default: `"Picard`")

     * [U] `"pseudo time integrator preconditioner`" [string] select the preconditioner to be used in the pseudo time integration method, choose one of `"Trilinos ML`", `"Hypre AMG`", or `"Block ILU`". (default: `"Hypre AMG`")

     * [U] `"pseudo time integrator linear solver`" [string] select the linear solver to be used in the pseudo time integration method. (default: `"AztecOO`")

     * [U] `"pseudo time integrator error control options`" [Array(string)] (default: `"pressure`")

     * [U] `"pseudo time integrator picard convergence tolerance`" [double] Picard convergence tolerance. (default: `"1.0e-8`", suggested range: 1.0e-10 ... 1.0e-4)

     * [U] `"pseudo time integrator picard maximum number of iterations`" [int] Picard maximum number of iterations. (default: `"400`", suggested range: 50 ... 1000)

   * [U] `"Linear Solver`" [list] Parameters for the linear solver used in single-phase steady-state solves, and in the damped Picard iteration to reach steady-state.

     * [U] `"linear solver tolerance`" [double] Set the tolerance for the AztecOO linear solver that may be used in a saturated steady state computation. (default: `"1.0e-16`", suggested range: 1.0e-20 ... 1.0e-14)

     * [U] `"linear solver maximum iterations`" [int] Set the maximum number of iterations for the linear solver that may be used in a saturated steady state computation. (default: `"100`", suggested range: 50 ... 1000)
 
     * [U] `"linear solver preconditioner`" [string] select the preconditioner to be used in the nonlinear solver for linear problems, choose one of `"Trilinos ML`", `"Hypre AMG`", or `"Block ILU`". (default: `"Hypre AMG`")

     * [U] `"linear solver iterative method`" [string] select the iterative method to be used in linear solvers, choose one of `"pcg`", or `"gmres`". (default: `"gmres`")

   * [U] `"MPC`" [list] Parameters for the multiprocess coordinator

     * [U] `"time integration rescue reduction factor`" [double] when the time integrator threatens to fail, for example, due to exceeding the number of limit iterations, or by threatening to diverge, the multiprocess coordinator will repeat the current time step with a time step that is reduced by this factor (default: `"0.5`").

   * [U] `"Nonlinear Solver`" [list] Parameters for the nonlinear solver used in time-integration.

     * [U] `"Nonlinear Solver Type`" [string] select the nonlinear solver type from `"NKA`", `"Newton`", and `"inexact Newton`".

     * [U] `"modify correction`" [bool] allows a process kernel to modify correction to a solution.(default: `"false`")

     * [U] `"update upwind frequency`" [string] define frequency of the updates for upwind direction: `"every nonlinear iteration`", `"every timestep`".(default: `"every timestep`")

   * [U] `"Preconditioners`" [list] Parameters to control the linear solver algorithms used in the preconditioner.

     * [U] `"Trilinos ML`" Parameter used by Trilinos multi-level solver, ML

       * [U] `"ML smoother type`" [string] The smoother to be used by ML, valid paramters are `"Jacobi`" (default), `"Gauss-Seidel`", and `"ILU`". (default: `"Jacobi`")

       * [U] `"ML aggregation threshold`" [double] This parameter influences the coarsening strategy of ML. The default is 0.0, which is a good choice for regular meshes. For meshes that have high aspect ratio cells, it is worth trying to set this parameter to something positive, but small, for example 0.0001. (default: `"0.0`", suggested range: 0.0 ... 0.1)

       * [U] `"ML smoother sweeps`" [int] The smoother will be called this many times before and after the coarse grid correction in the multilevel algorithm. (default: `"3`", suggested range: 1 ... 5) 

       * [U] `"ML cycle applications`" [int] This is the number of V-cycles that are performed in each preconditioner invocation. (default: `"2`", suggested range: 1 ... 5).

     * [U] `"Hypre AMG`" Parameters used by Hypre Algebraic Multigrid solver, BoomerAMG

       * [U] `"Hypre AMG tolerance`" [double] set a tolerance stopping criterion for the Hypre BoomerAMG preconditioner. If this is greater zero, then the preconditioner will run as many V-cycles as necessary to reach this prescribed accuracy, up to the maximum number of cycles that can also be specified as a parameter (see `"Hypre AMG cycle applications`"). (default: `"0.0`", suggested range: 0.0 ... 0.1)

       * [U] `"Hypre AMG cycle applications`" [int] the maximum number of V-cycles that are performed per preconditioner invocation. Note that if  `"Hypre AMG tolerance`" is zero, then this is the exact number of V-cycles that are performed per preconditioner invocation. (default: `"5`", suggested range: 1 ... 5)

       * [U] `"Hypre AMG smoother sweeps`" [int] the number of both pre and post smoothing sweeps. (default: `"3`", suggested range: 1 ... 5)

       * [U] `"Hypre AMG strong threshold`" [double] set this to 0.25 for a 2D problem, and to 0.5 for a 3D problem. (default: `"0.5`", suggested range: 0.2 ... 0.8) 

     * [U] `"Block ILU`" Parameters used by Trilinos Block ILU

       * [U] `"Block ILU overlap`" [int] specify the domain decomposition overlap that will be used in constructing the additive Schwarz block ILU preconditioner. (default: `"0`", suggested range: 0 ... 3)

       * [U] `"Block ILU relax value`" [double] corresponds to the Trilinos Ifpack ILU parameter `"fact: relax value`". (default: `"1.0`", suggested range: )

       * [U] `"Block ILU relative threshold`" [double] corresponds to the Trilinos Ifpack ILU parameter `"fact: relative threshold`". (default: `"1.0`", see Ifpack manual) 

       * [U] `"Block ILU absolute threshold`" [double] corresponds to the Trilinos Ifpack ILU parameter `"fact: absolute threshold`". (default: `"0.0`", see Ifpack manual)

       * [U] `"Block ILU level of fill`" [int] corresponds to the Trilinos Ifpack ILU parameter `"fact: level-of-fill`". (default: `"0`", suggested range: 0 ... 2)




  If the structured option is active, the following list of parameters is valid (Note: all lists here accept an optional sublist `"Expert Settings`".  Parameters listed in the expert area are not checked for validity/relevance during input reading stage, but are simply passed to the underlying implementation.)

  * [S] `"Structured Algorithm`" [list] (Optional) Additional controls for details of the structured-grid algorithm. 

   * [S] `"Expert Settings`" [int] Options passed to Amanzi that are not specifically checked for validity/relevance

     * [S] `"steady_limit_iterations"` [int] Maximum number of Newton iterations to attempt when solving for a single time step evolution of Richards equation. (default: "20", suggested range: 5 ... 200)

     * [S] `"steady_time_step_reduction_factor"` [double] When time step reduction is necessary during the steady calculation, use this factor. (default: `"0.8`", suggested range: 0.5 ... 0.9)

     * [S] `"steady_min_iterations"` [int] Maximum iteration count of successful Newton solve leading to time step increase of "steady_time_increase_factor". (default: "10", suggested range: 5 ... 100)

     * [S] `"steady_time_step_increase_factor"` [double] When time step increase is possible during the steady calculation, use this factor. (default: `"1.2`", suggested range: 1.1 ... 2.0)

     * [S] `"transient_time_step_reduction_factor"` [double] When time step reduction is necessary during the transient calculation, use this factor. (default: `"0.8`", suggested range: 0.5 ... 0.9)

     * [S] `"transient_time_step_increase_factor"` [double] When time step increase is possible during the steady calculation, use this factor. (default: `"1.2`", suggested range: 1.1 ... 2.0)

     * [S] `"maximum_time_step_size`" [double]: The maximum time step size allowed.

     * [S] `"initial_time_step_multiplier`" [double] (Optional) If internally computed time step used, it will be scaled by this factor (default value: 1)

     * [S] `"do_richard_init_to_steady`" [int]  If 1, triggers a psuedo-transient time-evolution of the initial data, prior to entering the `"Execution Mode`" phases descussed above.  (default: `"0`", suggested range: 0 ... 1)

     * [S] `"richard_init_to_steady_verbose`" [int]  Verbosity level of psuedo-transient time-evolution of the initial data, prior to entering the `"Execution Mode`" phases descussed above.  (default: `"0`", suggested range: 0 ... 4)

     * [S] `"steady_max_pseudo_time`" [double]  Stopping time for the psuedo-transient time-evolution of the initial data, prior to entering the `"Execution Mode`" phases descussed above.  (default: `"1.e10`", suggested range: 0 ... 1.e12)

     * [S] `"steady_time_step_reduction_factor`" [double]  Scale factor to reduce time step size for retry if Newton iterations fail.  (default: `"0.8`", suggested range: 0.1 ... 0.99)

     * [S] `"steady_time_step_increase_factor`" [double]  Scale factor to increase next step after successful solve with less than `"steady_min_iterations`" newton iterations.  (default: `"1.25`", suggested range: 1.1 ... 10)

     * [S] `"steady_min_iterations_2`" [int]  Iteration count of successful Newton solve leading to time step increase of `"steady_time_increase_factor_2`".  (default: `"0`", suggested range: 5 ... 100)

     * [S] `"steady_time_step_increase_factor_2`" [double]  Scale factor to increase next step after successful solve if iteration count of successful Newton solve is less than `"steady_min_iterations_2`".  (default: `"10`", suggested range: 1.1 ... 10)

     * [S] `"steady_max_consecutive_failures_1`" [int]  Number of failed time step attempts before reducing time step size by factor of `"steady_time_step_retry_factor_1`" (default: `"3`", suggested range: 5 ... 10)

     * [S] `"steady_time_step_retry_factor_1`" [double]  Scale factor to decrease time step after `"steady_max_consecutive_failures_1`" failed time steps.  (default: `"0.5`", suggested range: 0.1 ... 0.5)

     * [S] `"steady_max_consecutive_failures_2`" [int]  Number of failed time step attempts before reducing time step size by factor of `"steady_time_step_retry_factor_2`" (default: `"4`", suggested range: 5 ... 10)

     * [S] `"steady_time_step_retry_factor_2`" [double]  Scale factor to decrease time step after `"steady_max_consecutive_failures_2`" failed time steps.  (default: `"0.01`", suggested range: 0.01 ... 0.1)

     * [S] `"steady_time_step_retry_factor_f`" [double]  Scale factor to decrease time step after `"steady_max_consecutive_failures_2`" + 1 failed time steps.  (default: `"0.001`", suggested range: 0.001 ... 0.01)

     * [S] `"steady_max_num_consecutive_success`" [int]  Number of consecutive successful time step attempts, after which the time step will be increased by factor of `"steady_extra_time_step_increase_factor`" (default: `"15`", suggested range: 5 ... 100)

     * [S] `"steady_extra_time_step_increase_factor`" [double]  Scale factor to increase time step after `"steady_max_num_consecutive_success`" successful time steps.  (default: `"10`", suggested range: 5 ... 100)

     * [S] `"steady_abort_on_psuedo_timestep_failure`" [int]  If > 0, abort the run when the solver fails to successfully complete a time step.  (default: `"0`", suggested values: 0, 1)

     * [S] `"steady_use_PETSc_snes`" [bool]  If true, use a backward Euler discretization of Richards equation, and use the PETSC SNES software to drive the solution of the system.  (default: `"True`")

     * [S] `"steady_limit_function_evals`" [int]  If > 0, the maximum number of function evaluations during a single PETSC SNES time step solve.  Aborts if more are attempted.  (default: `"-1`", suggested values: -1, 1 ... 1.e10)

     * [S] `"richard_solver_verbose`" [int]  Verbosity of Richard solve. (default: `"1`", suggested values: 0 ... 3)

     * [S] `"richard_max_ls_iterations`" [int]  Maximum number of line search attempts before declaring Newton solver failure. (default: `"10`", suggested values: 5 ... 15)

     * [S] `"richard_ls_reduction_factor`" [double]  Factor to scale line search parameter for subsequent line search attempt. (default: `"0.1`", suggested values: .01 ... 0.9)

     * [S] `"richard_min_ls_factor`" [double]  Smallest allowable line search factor before declaring Newton solver failure. (default: `"1.e-8`", suggested values: .001 ... 1.e-10)

     * [S] `"richard_ls_acceptance_factor`" [double]  Maximum factor by which residual from previous Newton iterate is reduced by scaled Newton update (default: `"1.4`", suggested values: .9 ... 200)

     * [S] `"richard_monitor_line_search`" [int]  If > 0, print progress of line search. (default: `"0`", suggested values: 0, 1)

     * [S] `"richard_monitor_linear_solve`" [int]  If > 0, print progress of linear solve for Newton systems. (default: `"0`", suggested values: 0, 1)

     * [S] `"richard_use_fd_jac`" [bool]  If True, use finite-difference approximation for Jacobian in Newton system. (default: `"True`") - ANALYTIC JACOBIAN NOT CURRENTLY SUPPORTED

     * [S] `"richard_perturbation_scale_for_J`" [double]  Perturbation on scaled pressure values used to compute finite-difference Jacobian. (default: `"1.e-8`", suggested values: 1.e-12 ... 1.e-6)

     * [S] `"richard_use_dense_Jacobian`" [bool]  If True, use dense storage methods for Newton system. (default: `"False`")

     * [S] `"richard_upwind_krel`" [bool]  If True, use upwind saturation values to evaluate the relative permeability at a cell face.  (default: `"True`")

     * [S] `"richard_pressure_maxorder`" [int]  Polynomial order used to construct pressure gradients at coarse-fine interfaces. (default: `"3`", suggested values: 1 ... 4)

     * [S] `"richard_scale_solution_before_solve`" [bool] If True, scale pressure variable in SNES prior to solve. (default: `"True`")

     * [S] `"richard_semi_analytic_J`" [bool] If True, form numerical Jacobian by finite-differencing divergence of Darcy flux but using analytic form of time derivative.  (default: `"True`")

     * [S] `"steady_do_grid_sequence`" [bool] If True and richard_init_to_steady, psuedo-evolve coarsest only level solution, then interpolate solution to next finer level and repeat.  (default: `"True`")

     * [S] `"steady_grid_sequence_new_level_dt_factor`" [Array(double)] Factor by which to scale final psuedo time step from previous (coarser) steady solve in order to compute initial psuedo time step for next steady solve.  If more than one value given, each will be used in successive solves.

     * [S] `"max_n_subcycle_transport`" [int] Maximum number of level-0 subcycled transport time steps for each flow step.  Transport will be limited by an advective CFL stability constriant, so this will contribute to limiting the over step size taken. (default: `"10`", suggested values: 1 ... 20)

     * [S] `"cfl`" [double]: Fraction of stability-limited maximum time step allowed by the advective transport scheme. (default: "1", suggested values: .01 ... 1)

   * [S] `"Adaptive Mesh Refinement Control`" [list] (Optional) Additional details related to the adaptive mesh refinement algorithm. 

     * [S] `"Number Of AMR Levels`" [int] Maximum number of adaptive levels, including the base grid (default=1)

     * [S] `"Refinement Ratio`" [Array(int)] Grid spacing ratio between adjacent refinement levels.  One value required for each coarse level. Only values of 2 or 4 are supported.

     * [S] `"Do AMR Subcycling`" [bool] For integration of transport and chemistry, AMR subcycling time-steps each level with the same ratio of dx/dt, the levels are integrated and synchronized recursively.  If "`False"`, the time step is identical across levels.

     * [S] `"Regrid Interval`" [Array(int)] Number of base (coarse) grid time steps between regrid operations (one value > 0 required for each coarse level)

     * [S] `"Blocking Factor`" [Array(int)] Number by which each grid per level is evenly divisable in each dimension (typically used to guarantee multigrid hierachy depth).  A single value implies that the same is to be used for all levels, otherwise one value is required for each fine level.

     * [S] `"Number Error Buffer Cells`" [Array(int)] Number of coarse cells automatically tagged to surround user-tagged cells prior to generation of fine grids.  Used to guarantee buffer between refinement levels.

     * [S] `"Maximum Grid Size`" [Array(int)] Size of largest dimension of any mesh generated at each level.  A single value implies that the same value is to be used for all levels.

     * [S] `"Refinement Indicators`" [list] A list of user-labeled refinement indicators, REFINE.  Criteria will be applied in the order listed.

      * [S] REFINE [list] A user-defined label for a single refinement criteria indicator function.  For all criteria except `"Inside Region`" definition must indicate a `"Field Name`" selected from the list of known available field quantities (defined under "Time and Cycle specification").  Additionally, one must specify `"Regions`" (a list of user-named regions over which this criteria is to apply) and one of the following parameters refinement rules:

       * [S] `"Value Greater`" [double] The threshold value.  For each coarse cell, if the value of the given field is larger than this value, tag the cell for refinement

       * [S] `"Value Less`" [double] The threshold value.  For each coarse cell, if the value of the given field is smaller than this value, tag the cell for refinement

       * [S] `"Adjacent Difference Greater`" [double] The threshold value.  For each coarse cell, if the maximum difference of the values for the given field between adjacent neighbors is larger than this value, tag the cell for refinement.

       * [S] `"Inside Region`" [bool] Set this TRUE if all coarse cells in the identified list of regions should be tagged for refinement.

       The following optional parameters are available in order to fine-tune the refinement strategy:

       * [S] `"Maximum Refinement Level`" [int] If set, this identifies the highest level of refinement that will be triggered by this indicator

       * [S] `"Start Time`" [double] If set, this identifies the time after which this criteria will be applied

       * [S] `"End Time`" [double] If set, this identifies the time before which this criteria will be applied

  * [S] `"Diffusion Discretization Control`" [list] (Optional) Additional details related to the parabolic diffusion solver. Details to be added.

  * [S] `"Pressure Discretization Control`" [list] (Optional) Algorithmic options for pressure solve. Details to be added.

  * [S] `"Iterative Linear Solver Control`" [list] Detailed controls for linear solvers. Details to be added.

   * [S] `"Conjugate Gradient Algorithm`" [list] (Optional) Algorithmic options for CG Solver. Details to be added.

   * [S] `"Multigrid Algorithm`" [list] (Optional) Algorithmic options for Multigrid Solver. Details to be added.



Example:

.. code-block:: xml

  <ParameterList name="Execution Control">

    <Parameter name="Flow Model" type="string" value="Richards"/>
    <Parameter name="Transport Model" type="string" value="On"/>
    <Parameter name="Chemistry Model" type="string" value="Off"/>

    <ParameterList name="Time Integration Mode">
      <ParameterList name="Transient">
         <Parameter name="Start" type="double" value="0"/>
         <Parameter name="End" type="double" value="1.5768e9"/>
      </ParameterList>
    </ParameterList>

    <ParameterList name="Time Period Control">
      <Parameter name="Period Start Times" type="Array(double)" value="{6.1726667E10, 6.1731787E10, 6.1737054E10, 9.4672798E10}"/>
      <Parameter name="Initial Time Step" type="Array(double)" value="{60.0, 60.0, 60.0, 800.0}"/>
    </ParameterList>

    <Parameter name="Verbosity" type="string" type="High"/>

    <ParameterList name="Numerical Control Parameters">
      <ParameterList name="Adaptive Mesh Refinement Control">
        <Parameter name="Number Of AMR Levels" type="int" value="3"/>
        <Parameter name="Refinement Ratio" type="Array(int)" value="{4, 4}"/>
        <Parameter name="Regrid Interval" type="Array(int)" value="{2}"/>
        <Parameter name="Blocking Factor" type="Array(int)" value="{8, 8, 8}"/>
        <Parameter name="Maximum Grid Size" type="Array(int)" value="{16, 16, 16}"/>
        <Parameter name="Numbers Error Buffer Cells" type="Array(int)" value="{2, 1}"/>

        <Parameter name="Refinement Indicators" type="Array(string)" value="{Pc ref, Region ref}"/>
        <ParameterList name="Pc ref">
          <Parameter name="Maximum Refinement Level" type="int" value="1"/>
          <Parameter name="Field Name" type="string" value="Capillary Pressure"/>
          <Parameter name="Regions" type="Array(string)" value="{CCugr}"/>
          <Parameter name="Value Greater" type="double" value="1.e6"/>
        </ParameterList>
        <ParameterList name="Region ref">
          <Parameter name="Regions" type="Array(string)" value="{Hgr, CCugr}"/>
          <Parameter name="Inside Region" type="bool" value="TRUE"/>
        </ParameterList>
      </ParameterList>

      <ParameterList name="Basic Algorithm Settings">
        <ParameterList name="Expert Settings">
          <Parameter name="visc_abs_tol" type="double" value="1.e-14"/>
        </ParameterList>
      </ParameterList>
    </ParameterList>

  </ParameterList>


This example specifies that a time-dependent evolution of Richards equation is desired, evolving over the physical time interval, 0 to 1.5768e9 seconds.  Here a 3-level AMR hierarchy is desired, where refinement up to level 1 is based on a threshold value of capillary pressure in the CCugr region.  Additionally, fine grid up to the maximum allowed (3) is generated over the Hgr and CCugr regions.  The user has also set the expert setting for a parameter called "visc_abs_tol".


Domain
======

[SU] The `"Domain`" parameter list contains the spatial dimension.

Example:

.. code-block:: xml

  <ParameterList name="Domain">
    <Parameter name="Spatial Dimension" type="int" value="2"/>
  </ParameterList>

For unstructured Amanzi, this parameter can equal 2 or 3.

Mesh
====

Amanzi supports both structured and unstructured numerical solution approaches.  This flexibility has a direct impact on the selection and design of the underlying numerical algorithms, the style of the software implementations, and, ultimately, the complexity of the user-interface.  "Mesh`" is used to select between the following options:

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
     <ParameterList name="Structured">
       <Parameter name="Number of Cells" type="Array(int)" value="{100, 1, 100}"/>
       <Parameter name="Domain Low Coordinate" type="Array(double)" value="{0.0, 0.0, 0.0}" />
       <Parameter name="Domain High Coordinate" type="Array(double)" value="{103.2, 1.0, 103.2}" />
     </ParameterList>   
   </ParameterList>

Example of `"Unstructured`" mesh generated internally:

.. code-block:: xml

   <ParameterList name="Mesh">
     <ParameterList name="Unstructured">
       <ParameterList name="Generate Mesh">
         <ParameterList name="Uniform Structured">
           <Parameter name="Number of Cells" type="Array(int)" value="{100, 1, 100}"/>
           <Parameter name="Domain Low Coordinate" type="Array(double)" value="{0.0, 0.0, 0.0}" />
           <Parameter name="Domain High Coordinate" type="Array(double)" value="{103.2, 1.0, 103.2}" />
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

Under the `"Structured`" option, Amanzi also automatically defines regions for the coordinate-aligned planes that bound the domain,
using the following labels: `"XLOBC`", `"XHIBC`", `"YLOBC`", `"YHIBC`", `"ZLOBC`", `"ZHIBC`"

User-defined regions are constructed using the following syntax

 * [SU] "Regions" [list] can accept a number of lists for named regions (REGION)

   * Shape [list] Geometric model primitive, choose exactly one of the
     following [see table below]: `"Region: Point`", `"Region: Box`",
     `"Region: Plane`", `"Region: Layer`", `"Region: Polygon`", `"Region: Circle`", `"Region: Rotated Polygon`", `"Region: Swept Polygon`", `"Region: Logical`"

Amanzi supports parameterized forms for a number of analytic shapes, as well as more complex definitions based on triangulated surface files.  

+-----------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
|  shape functional name            | parameters                              | type(s)                      | Comment                                                                |
+===================================+=========================================+==============================+========================================================================+
| `"Region: Point"`  [SU]           | `"Coordinate`"                          | Array(double)                | Location of point in space                                             |
+-----------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
| `"Region: Box"` [SU]              | `"Low Coordinate`", `"High Coordinate`" | Array(double), Array(double) | Location of boundary points of box                                     |
+-----------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
| `"Region: Plane"`  [SU]           | `"Direction`", `"Location`"             | string, double               | direction: `"X`", `"-X`", etc, and `"Location`" is coordinate value    |
+-----------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
| `"Region: Polygonal Surface"` [SU]| `"Number of points`", `"Points`"        | int, Array(double)           | Number of polygon points and point coordinates in linear array. This   |
|                                   |                                         |                              | provides a set of faces with a normal for computing flux               |    
+-----------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
| `"Region: Polygon"`  [S]          |                                         | Array(double), Array(double) | V1=(x1,x2,...) and V2=(y1,y2,...) coordinates of an ordered sequence of|
|       (2D-only)                   | `"VerticesV1`", `"VerticesV2`"          |                              | points (x1,y1), (x2,y2),... defining a polygon.                        |
|                                   |                                         |                              | The first and last points are connected.                               |    
+-----------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
| `"Region: Ellipse"`  [S]          | `"Center`", `"Radius`"                  | Array(double), Array(double) | Coordinate (x,y,z), of center, and radii (rx, ry) of an ellipse        |
|       (2D-only)                   |                                         |                              | in the plane.                                                          |    
+-----------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
| `"Region: Rotated Polygon"`  [S]  |                                         | Array(double), Array(double),| V1=(x1,x2,...) and V2=(y1,y2,...) coordinates of an ordered sequence of|
|  (3D-only)                        | `"VerticesV1`", `"VerticesV2`",         | string, string, Array(double)| points (x1,y1), (x2,y2),... defining a polygon in the specified plane  |
|                                   | `"Plane`", `"Axis`", `"Reference Point`"|                              | (`"XY`", `"YZ`", `"XZ`"), rotated about the                            |
|                                   |                                         |                              | given axis (`"X`", `"Y`", `"Z`") that is located at the given          |
|                                   |                                         |                              | reference point, (x,y,z).                                              |
+-----------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
| `"Region: Swept Polygon"`  [S]    | `"VerticesV1`", `"VerticesV2`",         | Array(double), Array(double),| V1=(x1,x2,...) and V2=(y1,y2,...) coordinates of an ordered sequence of|
|  (3D-only)                        | `"Plane`", `"Extent`"                   | string, Array(double)        | points (x1,y1), (x2,y2),... defining a polygon in the specified plane  |
|                                   |                                         |                              | (`"XY`", `"YZ`", `"XZ`"), and swept over the extent, (min, max) in the |
|                                   |                                         |                              | direction normal to the plane                                          |
+-----------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
| `"Region: Logical"` [SU]          | `"Operation`", `"RegionList`"           | string, Array(string)        | Operation can be Union, Intersection, Subtraction, Complement          |
+-----------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
| `"Region: Labeled Set"` [SU]       | `"Label`", `"File`",                    | string, string,              | Set per label defined in mesh file (see below)                         |
|                                   | `"Format`", `"Entity`"                  | string, string               |  (available for frameworks supporting the `"File`" keyword)            |
+-----------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+
| `"Region: Color Function"` [SU]   | `"File`", `"Value`"                     | string, int                  | Set defined by color in a tabulated function file (see below)          |
+-----------------------------------+-----------------------------------------+------------------------------+------------------------------------------------------------------------+

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

* `"Region: Polygonal Surface`" defines a polygonal region on which mesh faces and
  nodes can be queried. NOTE that one cannot ask for cells in a polygonal surface
  region. In 2D, the "polygonal surface" region is a line and is specified by 2 points.
  In 3D, the "polygonal surface" region is specified by an arbitrary number of points.
  In both cases the point coordinates are given as a linear array. The polygon
  can be non-convex.

  The polygonal surface region can be queried for a normal. In 2D, the normal is
  defined as [Vy,-Vx] where [Vx,Vy] is the vector from point 1 to point 2.
  In 3D, the normal of the polygon is defined by the order in which points 
  are specified.

* `"Region: Logical`" Logical operations on regions allow for more
  advanced region definitions. At this time the Logical Region allows
  for logical operations on a list of regions.  In the case of Union
  the result is obvious, it is the union of all regions.  Similarly
  for Intersection. In the case of Subtraction, subtraction is
  performed from the first region in the list.  The Complement is a
  special case in that it is the only case that operates on single
  region, and returns the complement to it within the domain 'Entire
  Domain'.  Currently, multi-region booleans are not supported in the same expression.

.. code-block:: xml

  <ParameterList name="Lower Layers">
    <ParameterList name="Region: Logical">
      <Parameter name="Operation" type="string" value="Union"/>
      <Parameter name="RegionList" type="Array(string)" value="{Middle1, Middle2, Bottom}"/>
    </ParameterList>
  </ParameterList>

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
    </ParamterList>
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
    <ParameterList name="Flux plane">
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

In this example, "Top Section", "Middle Section" and "Bottom Section"
are three box-shaped volumetric regions. "Inflow Surface" is a
surface region defined in an Exodus II-formatted labeled set
file and "Outflow plane" is a planar region. "Sand" is a volumetric
region defined by the value 25 in color function file.



Material Properties
===================

The "material" in this context is meant to represent the media through with  fluid phases are transported.  In the literature, this is also referred to as the "soil", "rock", "matrix", etc.
Properties of the material must be specified over the entire simulation domain, and is carried out using the Region constructs defined above. For example, a single material 
may be defined over the `"All`" region (see above), or a set of materials can be defined over subsets of the domain via user-defined regions.
If multiple regions are used for this purpose, they should be disjoint, but should collectively tile the entire domain.  Each material requires (Section 2.6) a label and 
the following set of physical properties using the supported models described below.

* [SU] "Material Properties" [list] can accept multiple lists for named material types (MATERIAL)

 * [SU] MATERIAL [list] can accept lists to specify models, and `"Assigned Regions`" to specify where this model applies

  The flow related material properties *Intrinsic Permeability* or *Hydraulic Conductivity* must be specified, but not both:  

  * [SU] Intrinsic Permeability [list] Parameterized model for intrinsic permeability.  Choose exactly one of the following: `"Intrinsic Permeability: Uniform`", `"Intrinsic Permeability: Anisotropic Uniform`" (see below)

  * Hydraulic Conductivity [list] Parameterized model for intrinsic permeability.  Choose exactly one of the following: `"Hydraulic Conductivity: Uniform`", `"Hydraulic Conductivity: Anisotropic Uniform`" (see below)

  Additional ''Material Properties'' related to flow are:

  * [SU] Porosity [list] Parameterized model for porosity.  Choose exactly one of the following: `"Porosity: Uniform`" (see below)

  * [SU] Capillary Pressure [list] Parameterized mass density model.  Choose exactly one of the following: `"van Genuchten`" or `"Brooks Corey`"

  * [SU] Particle Density [list] Choose exatly one of the following: `"Particle Density: Uniform`". 

  * [SU] Specific Storage [list] Parameterized model for Specific Storage [L^-1]. Choose exactly one of the following: `"Specific Storage: Uniform`".

  * [SU] Specific Yield [list] Parameterized model for Specific Yield [-]. Choose exactly one of the following: `"Specific Yield: Uniform`".

  Material properties related to transport (dispersion and diffusion):

  * [SU] Dispersion Tensor [list] Parameterized model for Dispersion Tensor. Choose exactly one of the following: `"Dispersion Tensor:  Uniform Isotropic`".

..  * [SU] Molecular Diffusion [list] Parameterized model for
    a single molecular diffusion coefficient for all primary species [L^2 / time = m^2 / s]. Choose exactly one of the following: `"Molecular Diffusion: Uniform`".

  * [SU] Tortuosity [list] Parameterized model for the Tortuosity [-]. Choose exactly one of the following: `"Tortuosity: Uniform`".

  Material properties related to geochemistry are constant over the assigned regions and constant in time:

  * [SU] Mineralogy [list] List of minerals in the system.  Any mineral present in the phase definition but 
    not listed here must still be allocated in memory with default to zero parameter values.

    * [SU] `"MINERAL NAME`" [list] Name of a mineral from the phase definitions "Minerals" list.

      * [SU] `"Volume Fraction`" (double) [-] Uniform over the assigned region, constant in time (default: 0.0).
      * [SU] `"Specific Surface Area`" (double) [m^2 / m^3 bulk] Uniform over the assigned region, constant in time (default: 0.0).

  * [SU] `"Surface Complexation Sites`" [list]

    * [SU] `"SITE NAME`" [list] Name of a site from the phase definitions "Sorption Sites" list.

      * [SU] `"Site Density`" (double) [mol / m^3 bulk] Uniform over the assigned region, constant in time (default: 0.0)

  * [SU] `"Cation Exchange Capacity`" (double) [equivalent / m^3 bulk] Uniform over the assigned region, constant in time (default: 0.0)

  * [SU] `"Sorption Isotherms`" [list]

    * [SU] `"Solute Name`" [list] The name of one of the solutes from the phase definitions "Component Solutes" list.
 
      * [SU] `"Kd`" (double) [Kg H2O / m^3 bulk] molality-based distribution coefficient for this solute on this material (default: 0.0). If Kd is available in the more conventional units of mL/g or L/Kg, one needs to multiply that value by the water density [Kg water/L] and bulk density [Kg/m3 bulk]. Note: `"Kd`" is also used to enter the distribution coefficient in the Langmuir and Freundlich models. In the empirical Freundlich model, units will depend on the choice of n, i.e. [ mol^n / (m^3 bulk * Kg H2O)^n) ] . In the Langmuir model, units will be in [L H2O / mol]
      * [SU] `"Langmuir b`" (double) (Optional) [mol/m^3 bulk] Langmuir isotherm "b" coefficient, (default: 0.0)
      * [SU] `"Freundlich n`" (double) (Optional) [-] Freundlich isotherm "n" coefficient, (default to 1.0).

  Assigned regions are typically specified last:

  * [SU] `"Assigned Regions`" (Array(string)) a set of labels corresponding to volumetric regions defined above.  If any regions specified here are not three-dimensional, an error is thrown. (NOTE: [S] if layers in this list overlap spatially, this list implies the precedence ordering, right to left)

The following models can be specified for porosity (only `"Porosity: Uniform`" is supported at the moment):

* [SU] `"Porosity: Uniform`" [list] requires 
 
 * [SU] `"Value`" [double] to specify the constant value of porosity.

The following models can be specified for the intrinsic permeability of the material:

* [SU] `"Intrinsic Permeability: Uniform`" [list] requires 
 
 * [SU] `"Value`" [double] to specify the constant value of the intrinsic permeability

* [U] `"Intrinsic Permeability: File`" [list] requires 
 
 * [U] `"File`" [string] provides the name of the file containing the permeability field

 * [U] `"Format`" [string] specifies the format of the file (`"exodus`" is the only supported format at this time)

 * [U] `"Attribute`" [string] to specify the attribute name used to identify the permeability values

* [SU] `"Intrinsic Permeability: Anisotropic Uniform`" [list] requires
 
 * [SU] `"x`" [double] to specify the constant value of the intrinsic permeability in the x-direction; and

 * [SU] `"y`" [double] to specify the constant value of the intrinsic permeability in the y-direction; and

 * [SU] `"z`" [double] to specify the constant value of the intrinsic permeability in the z (vertical) direction.

 where the directions refer to the global cartesian coordinates.


The following models can be specified for the Hydraulic Conductivity of the material:

* [SU] `"Hydraulic Conductivity: Uniform`" [list] requires 
 
 * [SU] `"Value`" [double] to specify the constant value of the intrinsic permeability

* [SU] `"Hydraulic Conductivity: Anisotropic Uniform`" [list] requires
 
 * [SU] `"x`" [double] to specify the constant value of the intrinsic permeability in the x-direction; and

 * [SU] `"y`" [double] to specify the constant value of the intrinsic permeability in the y-direction; and

 * [SU] `"z`" [double] to specify the constant value of the intrinsic permeability in the z (vertical) direction.

where the directions refer to the global cartesian coordinates.  Note that internally Amanzi works with a pressure formulation and uses intrinsic permeability.  If Hydraulic Conductivity is specified, constant density and viscosity values are used to convert it to intrinsic permeability
(see Equation 3.25).  Hence, either Intrinsic Permeability or Hydraulic Conductivity must be specified, but not both.

The following models are currently supported for capillary pressure (Section 3.3.2):

* `"Capillary Pressure: None`" [list] requires no parameters, pc = 0

* [SU] `"Capillary Pressure: van Genuchten`" [list] requires 

 * [SU] `"alpha`" [double] to specify alpha in Equation 3.7.

 * [SU] `"Sr`" [double] to specify the residual saturation, s^r_l, in Equation 3.5.

 * [SU] `"m`" [double] to specify m in Equation 3.7.

 * [SU] `"ell`" [double] ''l'' in Equation 3.11 (default = 0.5)

 * [SU] `"Relative Permeability`" [string] (either (0) [U] `"Burdine`", or (2) [SU] `"Mualem`") determines n
   in Equation 3.10, and the form of relative permeability (either Equation 3.12, or Equation 3.11, respectively).

 * [SU] `"krel smoothing interval`" [double] If this parameter is positive, a cubic hermite interpolant in used in place of the van Genuchten relative permeability function when the capillary pressure is in the interval [0.0, krel smoothing interval]. The default for this parameter is 0.0, such that there is no relative premeability smoothing.  

 * [S] `"WRM Plot File`" [string] (Optional) name of ASCII text file to write 3-column, space-delimited data for water saturation, capillary pressure and relative permeability.  If not give, no file is created.  Also, if file exists, it will be overwritten.

 * [S] `"WRM Plot File Number Of Points`" [int] (Optional, defaults to 1000) Number of evaluation points used to create the WRM Plot File, spaced uniformly in saturation between Sr+epsilon and 1.

* [SU] `"Capillary Pressure: Brooks Corey`" [list] requires

 * [SU] `"lambda`" [double] to specify lambda in Equation 3.9

 * [SU] `"alpha`" [double]  to specify alpha in Equation 3.9 

 * [SU] `"ell`" [double] to specify ''l'' in Equation 3.12 (default is 2.0)

 * [SU] `"Sr`" [double] to specify residual saturation, s^r_l, in Equation 3.5

 * [SU] `"Relative Permeability`" [string] (either (0) `"Burdine`", or (2) `"Mualem`") chooses the form of the
   relative permeability (either Equation 3.15, or Equation 3.14, respectively)

 * [SU] `"krel smoothing interval`" [double] (default value gives no relative permeability smoothing).

 * [S] `"WRM Plot File`" [string] (Optional) name of ASCII text file to write 3-column, space-delimited data for water saturation, capillary pressure and relative permeability.  If not give, no file is created.  Also, if file exists, it will be overwritten.

 * [S] `"WRM Plot File Number Of Points`" [int] (Optional, defaults to 1000) Number of evaluation points used to create the WRM Plot File, spaced uniformly in saturation between Sr+epsilon and 1.


The following models can be specified for particle density (only `"Particle Density: Uniform`" is supported at the moment):

* [SU] `"Particle Density: Uniform`" [list] requires 
 
 * [SU] `"Value`" [double] to specify the constant value of rock density.


The following models are currently supported for Specific Yield.

* [SU] `"Specific Yield: Uniform`" [list] requires

 * [SU] `"Value`" [double] to specify specific yield.


The following models are currently supported for Specific Storage.

* [SU] `"Specific Storage: Uniform`" [list] requires

 * [SU] `"Value`" [double] to specify specific storage.

The following models are currently supported for the dispersion tensor
in transport

* [SU] `"Dispersion Tensor: Uniform Isotropic`" (see Equation 4.9) [list] requires

 * [SU] `"alphaL`" [m]  the longitudinal dispersion  (default 0)
 * [SU] `"alphaT`" [m]  the transverse dispersion    (default 0)

The following models are currently supported for the Molecular
Diffusion coefficient.

* [SU] `"Molecular Diffusion: Uniform`" [list] requires

 * [SU] `"Value`" [double] to specify diffusion coefficient [m^2/s] (see Equation 4.15).

The following models are currently supported for Tortuosity.

* [SU] `"Tortuosity: Uniform`" [list] requires

 * [SU] `"Value`" [double] to specify Tortuosity [-] (see Equation 4.18).


Example:

.. code-block:: xml

  <ParameterList name="Material Properties">
    <ParameterList name="Backfill">
      <ParameterList name="Particle Density: Uniform">
        <Parameter name="Value" type="double" value="2.8e3"/>
      </ParameterList>
      <ParameterList name="Intrinsic Permeability: Anisotropic Uniform">
        <Parameter name="Horizontal" type="double" value="2.05e-8"/>
        <Parameter name="Vertical" type="double" value="2.05e-9"/>
      </ParameterList>
      <ParameterList name="Porosity: Uniform">
        <Parameter name="Value" type="double" value="0.38"/>
      </ParameterList>
      <ParameterList name="Capillary Pressure: van Genuchten">
        <Parameter name="alpha" type="double" value="2.14e-4"/> 
        <Parameter name="sr" type="double" value="0"/> 
        <Parameter name="m" type="double" value=".601"/> 
        <Parameter name="Relative Permeability" type="string" value="Mualem"/>
      </ParameterList>
      <Parameter name="Assigned regions" type="Array(string)" value="{Top Region, Bottom Region}"/>
    </ParameterList>

In this example, the material `"Backfill`" (which fills `"Bottom Region`" and `"Top Region`") has a
van Genuchten model for capillary pressure and a Mualem closure for relative permeability.  It also has an
anisotropic permeability which is uniform throughout the domain.


Fluid Phases
====================================

The "Fluid Phases" parameter list is used to specify fluid phases. A phase is defined as a homogeneous mixture of its chemical constituents. In the current version of Amanzi the aqueous phase serves as a reference phase in terms of which the composition all other fluid phases are derived through chemical equilibrium relations in the form of mass action equations. For the aqueous phase, the `"Fluid Phases`" parameter list identifies a set of independent variables through a flow mode (pressure equation) and a list of primary species (also referred to as basis species or components) that fully determine the chemical composition of each fluid phase in the system.  In the current version of the Amanzi the flow mode corresponds to a single liquid phase in a variably saturated porous medium, commonly referred to as Richards equation. The flow equation and primary species reactive transport equations are sequentially coupled.

Primary species, basis species or chemical components
-----------------------------------------------------

The primary species must be chosen from chemical constituents in the aqueous reference phase, but their choice is otherwise arbitrary except that they must form a linearly independent set of species, i.e. no linear combination of the primary species can exist which forms a valid chemical reaction. The concentrations of the remaining chemical constituents in the various fluid phases, referred to as secondary species, are obtained from the primary species concentrations through appropriate mass action relations under conditions of chemical equilibrium for given temperature and pressure conditions.

Each primary species has associated with it a total component concentration and a free ion concentration. The total concentration for each primary species is a sum of its free ion concentration in the aqueous phase and its stoichiometric contribution to all secondary species, which may also include other fluid phases for which it is in equilibrium. Amanzi splits the total primary species concentrations into a set of total concentrations for each fluid phase, and a total sorbed concentration. Mineral concentrations are not included in the total primary species concentrations.

In a general problem, multiple fluid phases may coexist in a mesh cell (e.g. aqueous/liquid, gaseous, etc.), with each phase comprised of a number of chemical constituents. The chemical constituents making up a fluid phase are typically divided into the solvent, the dominant species in the phase such as H2O in an aqueous phase, and the remaining "solute" species. All of these species may participate in various chemical reactions either as homogeneous reactions within a particular phase, or heterogeneous reactions involving more than one phase, for example, aqueous, solid and gas phases. Mineral reactions are treated as kinetically controlled with a reaction rate term appearing in the primary species transport equations. For each mineral an additional mass transfer equation is solved to obtain its spatial distribution throughout the computational domain. Sorbed species involving ion exchange and surface complexation reactions are treated as local equilibrium reactions with the sorbed concentration obtained through a mass action relation.

During initialization, Amanzi performs a distribution of species calculation that partitions the primary species concentrations among the secondary species within each fluid phase and equilibrates the aqueous solution with any specified minerals or gases. Various options may be used to constrain the speciation calculation, such as specifying charge balance, pH, total or free ion primary species concentration, total aqueous plus sorbed concentration, equilibrium with minerals and gases, and other options. 

In addition, certain reactions such as mineral precipitation and dissolution may affect the flow properties of the porous medium itself during the simulation through changes in porosity, permeability and tortuosity. Fluid properties (e.g. fluid density) may be affected through changes in species concentrations, temperature and pressure. While Amanzi does not currently support the effect of chemical reactions on material or fluid properties - the specification here, however, allows for the existence of the necessary input data framework and data structures to include such processes. Clearly, these specifications are highly problem dependent, so Amanzi attempts to provide a generalized interface to accommodate a variety of scenarios.

Given the free ion concentration of each primary species (and if there is more than one phase, a specification of the thermodynamic relationships that determine the partitioning between fluid phases, one can reconstruct the concentration of the primary and secondary species in each fluid phase. As a result only the primary species are maintained in the state data structures for each fluid phase. In addition, mineral concentrations and corresponding specific surface areas must also be stored in a state data structure.

Specification of Amanzi's numerical state is organized fundamentally around the list of fluid and solid phases that are present. Each fluid phase requires a specification of its physical properties (Section 4.6), and a list of its primary species. For each phase, Amanzi requires a label, and a list of chemical constituents. For each species, a group membership is specified. Note that Amanzi will eventually support the use of a master chemistry database, where a list of chemical species including aqueous, gaseous, surface complexes and mineral species together with their reaction stoichiometry, equilibrium constants over a range of temperatures and pressures, charge and other properties are defined. In that case, inclusion of a particular species in the Amanzi input file is conditioned on its presence in the appropriate section of the master thermodynamic database.

Sources and Initial and Boundary Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Fluid phases and the chemical constituents contained in them, require boundary conditions over the surface bounding the computational domain (Sections 3.3, 3.6, 3.10 and 4.3). 
Generally, boundary conditions are determined by specifying the phase pressure (Dirichlet condition), Darcy velocity (Neumann condition), or the phase saturation (Dirichlet condition) at the boundary. 
The fluid composition at a boundary may be specified either through Dirichlet or Neumann conditions. For simplicity, any boundary conditions not explicitly set in the input are defaulted to outflow with a zero gradient applied to each primary species. Volumetric source terms, used to model infiltration (Section 3.7) and a wide variety of production and loss processes, are defined for each phase, if applicable, and include the concentration or flux of any species that are carried into the domain with that phase. However, sources and sinks are not currently supported in Amanzi.

In order to support the rather general specification requirements (involving combinations of different fluid phases), it is necessary to first define the composition of the "state" of the system being simulated by identifying all fluid phases and chemical constituents that will be present in the system. We do this hierarchically, first by fluid phase then by chemical constituent:


Generalized Specification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Phases
=======================================

The `"Phases`" parameter list is used to specify components of each of the phases that are mobile, and solutes that are contained within them.  For each
phase, the list identifies the set of all independent variables that are to be stored on each discrete mesh cell.

Phases, components and solutes
------------------------------

The terminology for flow in porous media can be somewhat ambiguous between the multiphase and groundwater communities, particurly in regards to "components", "solutes" and "chemicals".  Since Amanzi is designed to handle a 
wide variety of problems, we must settle on a nomenclature for our use here.  In the general problem, multiple "phases" may coexist in the domain (e.g. gaseous, aqueous/liquid, etc), and each is
comprised of a number of "components" (section 2.2).  In turn, each component may carry a number of "solutes" and some of these may participate
in chemical reactions.  As a result of reactions, a chemical source or sink term may appear for the solutes involved in the reaction, including solutes in other mobile phases or in the material matrix.  
Additionally, certain reactions such as precipitation may affect the flow properties of the material itself during the simulation, and 
some might affect the properties of the fluid (e.g. brines affect the liquid density). While Amanzi does not currently support chemical reactions and thermal processes, the specification here allows for the existence of
the necessary data structures and input data framework.  Note that if solute concentrations are significant, the system may be better modeled with that solute treated as a separate component.  Clearly, these definitions
are highly problem-dependent, so Amanzi provide a generalized interface to accommodate a variety of scenarios.

Currently in Amanzi, solutes are transported in the various phase components, and are treated in "complexes".  Each complex is typically in chemical equilibrium with itself and does not undergo phase change.
Under these conditions, knowledge of the local concentration of the "basis" or "primary" species (the terms are used here interchangeably) in a chemical complex is sufficient to determine the concentrations of all related secondary species
in the phase. Each basis species has a total component concentration and a free ion concentration. The total component concentration for each basis species is a sum of the
free ion concentrations in the phase components and its stoichiometric contribution to all secondary species. Amanzi splits the total component concentration into a set of totals for each of the transported phase components,
and a total sorbed concentration. Given the free ion concentration of each basis species (and if there is more than one phase, a specification of the thermodynamic relationships that determine the partitioning 
between phase components (if mass transfer is allowed - not in current Amanzi), we can reconstruct the concentration of the primary and secondary species in each phase. As a result only the basis species are maintained in the state
data structures for each phases component.

In addition to solutes in the transported phases, there may be various immobile chemical constituents within the
porous media (material) matrix, such as "minerals" and "surface complexes". Bookkeeping for these constituents is managed in Amanzi
data structures by generalizing the "solute" concept - a slot in the state is allocated for each of these immobile species, but their concentrations are
not included in the transport/flow components of the numerical integration.  To allow selective transport of the various solutes, Amanzi
uses the concept of solute groups.   The aqueous solute concentrations are typically treated together as a group, for example, and often represent the only 
chemical constituents that are mobile.  Thus, the current Amanzi will assume that any other groups specified in an Aqueous phase are immobile.

Specification of Amanzi's state is organized fundamentally around the list of phases that are present.  Each phase requires a 
a specification of its physical properties (Section 4.6), and a list of its components.  For each component,
Amanzi requires a label, and a list of solutes.  For each solute, a group membership is specified.
Note that Amanzi will eventually support the use of a master chemistry database, where the solute complexes and their chemical activity are defined.  In that case, inclusion of a particular solute in the
Amanzi input file will be conditioned on its presence in the appropriate section of the master list.

Sources and Initial and Boundary Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Mobile phase components, and solutes contained in them, require boundary conditions along the entire surface bounding the computational domain (Sections 3.3, 3.6, 3.10 and 4.3).  Generally, boundary conditions are
specified in porous media systems by giving either the phase pressure or Darcy velocity on the boundary, and/or the component saturations.  Since mobile solutes are carried with the resulting flow,
inflowing boundary conditions for solutes are typically specified using Dirichlet conditions that define the effective solute concentration in the incoming flow.
For simplicity here, any boundary conditions not explicitly set in the input are defaulted to outflow with a zero gradient applied to any transport solutes. 
Volumetric source terms, used to model infiltration (Section 3.7) and a wide variety of production and loss processes, are defined for each phase component, if applicable, and include the distribution of any solutes that are carried into the domain with the phase component.  However, sources are not currently supported in Amanzi.

In order to support the rather general specification requirements (involving combinations of phase pressures and component saturations), we must first define the composition of the "state" of the simulations by identifying all phases, components and solutes that will be present in the system.  We do this hierarchically, first by phase then by component:
/
* [SU] `"Phase Definitions`" [list] can accept lists of named phases (currently PHASE can be `"Aqueous`" or `"Solid`").

 * [SU] `"Aqueous`" phase [list] can accept the following lists: `"Phase Properties`", `"Phase Components`"

  * [SU] `"Phase Properties`" can accept models for viscosity and density

   * [SU] Density [list] Parameterized model for phase mass density.  Choose exactly one of the following: `"Density: Uniform`", `"Density: File`" (see below)

   * [SU] Viscosity [list] Parameterized model for phase viscosity.  Choose exactly one of the following: `"Viscosity: Uniform`", `"Viscosity: File`" (see below)

  * [SU] `"Phase Components`" can accept COMP [list] named after a user-defined phase component (e.g., Water).

   * [SU] COMP [list] can accept a list of solutes carried by the component.

    * [SU] `"Component Solutes`" [Array(string)] List of primary or basis species for the aqueous solutes in the system. The order of this list must be the same as the order in the chemistry database file.

    * [SU] SOLUTE [list] can accept a solute name from the above list. 

      * [SU] `"Molecular Diffusivity`" [double] (Optional) specify the molecular diffusivity coefficient [m^2/s] for the current solute
        
      * [SU] `"First Order Decay Constant`" [double] (Optional) specify the decay constant for the current solute

 * [SU] `"Solid`" phase [list] can accept the following parameters: `"Minerals`", and `"Sorption Sites`"

  * [SU] `"Minerals`" [string array] list of minerals present in the system. The order of this list must be the same as the order in the chemistry database file.
  * [SU] `"Sorption Sites`" [string array] list of surface sites present in the system. The order of this list must be the same as the order in the chemistry database file.

Next, we specify the initial conditions.  Note that support is provided for specifying initial data on the phases and/or components simultaneously (the capillary pressure relationships are used to convert between the various options).  Thus, boundary conditions on the phases and components are specified together.  The solutes are specified afterward, organized first by phase then component.  If a solute exists in more than one phase/component, a thermodynamic relationship is required to partition the distribution - Amanzi does not currently support such a situation.

* [SU] `"Initial Conditions`" [list] accepts labels, IC, of named initial condition specifications 

 * [U] `"Init from Checkpoint File`" [string] (Optional) specify the checkpoint file that all fields are to be initialized from. If this parameter is present, all initial conditions are ignored.

 * [SU] IC [list] label for an initial condition, accepts initial condition function names, and parameters to specify assigned regions and solute initial conditions

  * [SU] Function [list] Parameterized model to specify initial profiles.  Choose exactly one of the following: `"IC: Uniform Pressure`", `"IC: Linear Pressure`", `"IC: Uniform Saturation`", `"IC: Linear Saturation`", `"IC: Uniform Velocity`" (see below)

  * [SU] `"Assigned Regions`" [Array(string)] list of regions to which this condition is assigned.  Note [S] when multiple regions specified overlap, this list implies a precedence, ordered right to left.

  * [SU] `"Solute IC`" can accept PHASE (labels of phases defined above)

   * [SU] PHASE [list] can accept COMPONENT (labels of components defined above), or keyword `"Alquimia`" to support geochemical conditions.

    * [SU] COMPONENT [list] can accept SOLUTE (label of solute defined above)

     * [SU] Component IC [list] Parameterized model for initial component conditions; only `"IC: Uniform Concentration`" is supported (see below) in units of molarity (moles/volume).

Next, we specify boundary conditions.  Again, support is provided for specifying boundary conditions on the aqueous phase (flow), and on the solutes 
(total component concentration of primary species in reactive transport).

* [SU] `"Boundary Conditions`" [list] accepts labels, BC, of named boundary condition specifications 

 * [SU] BC [list] label for a boundary condition, accepts boundary condition function names, and parameters to specify assigned regions and solute boundary conditions

  * [SU see below] Function [list] Parameterized model to specify boundary conditions.  Choose exactly one of the following: `"BC: Uniform Pressure`", `"BC: Linear Pressure`", `"BC: Hydrostatic`", `"BC: Linear Hydrostatic`", `"BC: Flux`", `"BC: Inflow`", `"BC: Impermeable`", `"BC: Zero Flow`" (see below)

  * [SU] `"Assigned Regions`" [Array(string)] list of regions to which this condition is assigned

  * [SU] `"Solute BC`" can accept PHASE (labels of phases defined above)

   * [SU] PHASE [list] can accept COMPONENT (labels of components defined above)

    * [SU] COMPONENT [list] can accept SOLUTE (label of solute defined above), or keyword `"Alquimia`" to support geochemical conditions

     * [SU] BC function [list] Parameterized model to specify the contcentration profile, only `"BC: Uniform Concentration`" is supported (see below) in units of molarity (moles/volume).

Finally, we specify sources.  Support is provided for specifying sources on the aqueous phase (flow), and for the solutes (total component concentration of primary species in reactive transport).

* [SU] `"Sources"` [list] accepts labels, SOURCE, of named source specifications

 * [SU] SOURCE [list] label for a source term, accepts source function names, and parameters to specify assigned regions and solute source conditions.

  * [SU] Function [list] Parameterized model to specify source. Choose exactly one of the following: `"Source: Uniform`", `"Source: Volume Weighted`", `"Source: Permeability Weighted`" (see below).
  
  * [SU] `"Assigned Regions`" [Array(string)] list of regions to which this condition is assigned

  * `"Solute SOURCE`" can accept PHASE (labels of phases defined above)

   * PHASE [list] can accept COMPONENT (labels of components defined above)

    * COMPONENT [list] can accept SOLUTE (label of solute defined above), or keyword `"Alquimia`" to support geochemical conditions

     * Source function [list] Parameterized model to specify the concentration profile, `"Source: Uniform Concentration`" and `"Source: Flow Weighted Concentration`" are supported (see below) in units of molarity (moles/volume).

The following initial condition parameterizations are supported:

* [U] `"IC: Uniform Saturation`" requires `"Value`" [double] OR `"Geochemical Condition`" [string] if Alquimia is providing initial conditions.

* [U] `"IC: Linear Saturation`" requires `"Reference Point`" (Array(double)), `"Reference Value`" [double], and  `"Gradient Value`" (Array(double)) 

* [SU] `"IC: Uniform Pressure`" requires `"Value`" [double]

* [SU] `"IC: Linear Pressure`" requires `"Reference Point`" (Array(double)), `"Reference Value`" [double], and  `"Gradient Value`" (Array(double)) 

* [SU] `"IC: Uniform Velocity`" requires `"Velocity Vector`" (Array(double)).

* [SU] `"IC: Uniform Concentration`" [list] requires `"Value`" and `"Free Ion Guess`", OR `"Geochemical Condition`"

  * [SU] `"Value`" [double]  total component concentration in units of molarity (moles/volume)
  * [SU] `"Free Ion Guess`" [double]  estimate of the free ion concentration for this solute; used to help convergence of the initial solution of the chemistry.
  * [SU] `"Geochemical Condition`" [String]: name of a geochemical condition defined in Alquimia's chemistry engine input file or in the Chemistry block.

The following boundary condition parameterizations are supported:

* [SU] `"BC: Flux`" requires a TIME FUNCTION (see below), where VALUE_NAME is one of the following: 

    * [SU]  `"Inward Volumetric Flux`" [Array(double)], 

    * [SU] `"Inward Mass Flux`" [Array(double)], 

    * [SU]  `"Outward Volumetric Flux`" [Array(double)], or

    * [SU] `"Outward Mass Flux`" [Array(double)]. 

  Here volumetriuc flux is interpreted as meters cubed per meters squared per second, and mass
  flux is interpreted as kilogramms per meter squared per
  second. Inward or outward refers to the flux being in the direction
  of the inward or outward normal to each face of the boundary region,
  respectively. With `"Inward`" fluxes there is an additional
  parameter to describe how sloping topography is handled:

    * [U] "rainfall" [bool] indicates that the mass flux is defined with respect
      to the gravity vector and the actual influx depends on boundary
      slope (default value is "false").

* [SU] `"BC: Uniform Pressure`" requires a TIME FUNCTION (see below), where VALUE_NAME is `"Values`" [Array(double)].

* [SU] `"BC: Linear Pressure`" [list] requires `"Reference Value`" [double] `"Reference Point`" [Array(double)] `"Gradient Value`" [Array(double)]

* [U] `"BC: Seepage`" [list] requires a TIME FUNCTION (see below), where VALUE_NAME is one of `"Inward Mass Flux`" [Array(double)] or `"Inward Volumetric Flux`" [Array(double)].  Here volumetriuc flux is interpreted as meters cubed per meters squared per second, and mass flux is interpreted as kilogramms per meter squared per second. Inward refers to the flux being in the direction of the inward normal to each face of the boundary region, respectively.
    * [U] "rainfall" [bool] indicates that the mass flux is defined with respect
      to the gravity vector and the actual influx depends on boundary
      slope (default value is "false").

* [SU] `"BC: Hydrostatic`" [list] requires a TIME FUNCTION (see below), where VALUE_NAME is `"Water Table Height`" [Array(double)].  Optional parameters supported by the unstructured framework include: `"Coordinate System`" [string] (available options: `"Relative`", `"Absolute`" [DEFAULT]), and `"Submodel`" [string] (available options: `"No Flow Above Water Table`", `"None`" [DEFAULT])

* [U] `"BC: Linear Hydrostatic`" [list] requires `"Reference Water Table Height`" [double] `"Reference Point`" [Array(double)] `"Gradient Value`" [Array(double)]

* `"BC: Impermeable`"  requires no parameters

* [SU] `"BC: Zero Flow`"  [list] takes no additional parameters

* [S] `"BC: Zero Gradient`" [list] takes no additional parameters

* [SU] `"BC: Uniform Concentration`" [list] requires a TIME FUNCTION (see below), where VALUE_NAME is `"Geochemical Condition`" [Array(string)] (if using Alquimia, names of geochemical conditions defined in Alquimia's chemistry engine input file), or `"Values`" [Array(double)] otherwise.  

The following source parameterizations are supported.

* [SU] `"Source: Uniform`" [kg/s/m^3] requires a TIME FUNCTION (see below), where VALUE_NAME is `"Values`" [Array(double)].

* [SU] `"Source: Volume Weighted`" [kg/s] requires a TIME FUNCTION (see below), where VALUE_NAME is `"Values`" [Array(double)].

* [SU] `"Source: Permeability Weighted`" [kg/s] requires a TIME FUNCTION (see below), where VALUE_NAME is `"Values`" [Array(double)].

* [SU] `"Source: Uniform Concentration`" [mol/s/m^3] uses a volume weighting to distribute the source uniformally over the specified region(s).  Requires a TIME FUNCTION (see below), where VALUE_NAME is `"Geochemical Condition`" [Array(string)] (if using Alquimia, names of geochemical conditions defined in Alquimia's chemistry engine input file), or `"Values`" [Array(double)] otherwise.  

* [SU] `"Source: Flow Weighted Concentration`" aligns the spatial distribution of the concentration with the distribution selected for the flow. Requires a TIME FUNCTION (see below), where VALUE_NAME is `"Geochemical Condition`" [Array(string)] (if using Alquimia, names of geochemical conditions defined in Alquimia's chemistry engine input file), or `"Values`" [Array(double)] otherwise.



TIME FUNCTIONS
~~~~~~~~~~~~~~

Boundary condition functions utilize a parameterized model for time variations that is either piecewise constant or piecewise linear, where data values are specified via the parameter VALUE_NAME (see above).  For example, in the case that VALUE_NAME=`"Values`" [Array(double)]:

.. code-block:: xml

      <Parameter name="Times" type="Array(double)" value="{1, 2, 3}"/>
      <Parameter name="Values" type="Array(double)" value="{10, 20, 30}"/>
      <Parameter name="Time Functions" type="Array(string)" value="{Constant, Linear}"/>    

This defines two time intervals, [1,2) and [2,3), three <time,time value> pairs, <1,10>, <2,20>, and <3,30>, and
two time functions, `"Constant`" and `"Linear`" that correspond to the intervals [1,2) and [2,3) respectively.

The type  of the function (linear or constant) on the intervals is specified by the `"Time Functions`" parameter. 
For a particular interval it can be either `"Linear`" or `"Constant`". 

- If the function type is `"Constant`", the function has the value that is specified as the time value at the left 
  end of the interval. In this example, we have f(t) = 10, for t in the interval [1,2).

- If the function type is `"Linear`", the function is the linear interpolant that is defined by the <time,time value>
  pairs at the two endpoints of the interval. In this example, these two  <time,time value> pairs are <2,20> and <3,30> and 
  we have f(x) = (30-20)*(x-2)/(3-2) + 20, for t in the interval [2,3). Note that this is indeed the linear function 
  for which f(2) = 20 and f(3) = 30.

- By convention, the function is constant for t in (-inf,1) and for t in [3,+inf) such that 

  - f(x) = 10, for x <1, specified by the first <time,time value> pair, and
  - f(x) = 30, for x >= 3, specified by the last <time,time value> pair. 

NOTE: If a single constant value is to be applied over all time, an abbreviated form is available:

.. code-block:: xml

      <Parameter name="Values" type="Array(double)" value="{10}"/>


Example Phase Definition
~~~~~~~~~~~~~~~~~~~~~~~~
Due to its length, an XML example of the `"Phases`" parameter list appears in the example appended to this specification.

Chemistry
=========

The chemistry list is needed if the Chemistry model is set to `"Alquimia`" or `"Amanzi`". The parameters within this list vary according to which chemistry process kernel is being used.

* [SU] `"Chemistry`" [list] (*`"Alquimia`" chemistry process kernel*)

  * [SU] `"Engine`" [string] The name of the backend chemistry engine (e.g. `"PFloTran`").
  * [SU] `"Engine Input File`" [string] The file specifying the input for the backend chemistry engine.
  * [SU] `"Max Time Step (s)`" [double] The maximum time step that chemistry will allow the MPC to take.

  * [SU] `"Geochemical Conditions`" [list] (*optional*, allows definition of geochemical conditions within XML.)

    * [SU] CONDITION NAME [list] The geochemical condition, defined in terms of aqueous and mineral constraints.

      * [SU] AQUEOUS CONSTRAINT NAME [list] Entry for an aqueous constraint involving a species and/or a mineral.

        * [SU] TYPE [double] The type of aqueous constraint (total_aqueous, total_sorb, free, mineral, gas, pH, charge) and its associated value.

        * [SU] SPECIES [string] The name of any associated mineral species for a mineral constraint.

      * [SU] MINERAL SPECIES [list] Entry for a mineral constraint mentioned in an aqueous constraint.

        * [SU] `"Volume Fraction`" [double] Volume fraction for the mineral.

        * [SU] `"Specific Surface Area`" [double] Specific surface area for the mineral.

* [SU] `"Chemistry`" [list] (*Native `"Amanzi`" chemistry process kernel*)

  * [SU] `"Thermodynamic Database`" [list]
   
    * [SU] `"File`" [string] Name of the chemistry database file, relative to the execution directory.
    * [SU] `"Format`" [string] The format of the database file. Actual database format is not XML and is the same as described for the 2010 demo with additions for the new chemical processes. Valid values: `"simple`".

  * [SU] `"Verbosity`" [string array] Control how much information the chemistry library will log. Valid values: "silent", "terse", "verbose", "warnings", "errors".
  * [SU] `"Activity Model`" [string] The type of model used for activity corrections. Valid values: "unit", "debye-huckel".
  * [SU] `"Tolerance`" [float] Tolerance in newton solves inside the chemistry library.
  * [SU] `"Maximum Newton Iterations`" [int] Maximum number of iteration the chemistry library can take.
  * [SU] `"Output File Name`" [string] (Optional) A file name that the chemistry library should use to write simulation information and debugging info. An empty string (default) indicates that the chemistry library should not write to a file.
  * [SU] `"Use Standard Out`" [bool] A flag indicating whether the chemistry library can write simulation information and debugging info to standard out. Default is true, so the amanzi u/s drivers will need to set this appropriately on mpi/openmp processes.
  * [SU] `"Auxiliary Data`" [string array] Additional chemistry related data that the user can request be saved to vis files. Currently "pH" is the only variable supported.
  * [SU] `"Max Time Step (s)`" [double] The maximum time step that chemistry will allow the MPC to take.


Output
======

Output data from Amanzi is currently organized into four specific groups: `"Observation Data`", `"Visualization Data`", `"Checkpoint Data`", `"Walkabout Data`", and `"Log Data`".  Each of these is controlled in different ways, reflecting their intended use.

* `"Checkpoint Data`" is intended to represent all that is necesary to repeat or continue an Amanzi run.  The specific data contained in a Checkpoint Data dump is specific to the algorithm optoins and mesh framework selected.  Checkpoint Data is special in that no interpolation is perfomed prior to writing the data files; the raw binary state is necessary.  As a result, the user is allowed to only write Checkpoint Data at the discrete intervals of the simulation.

* `"Visualization Data`" is intended to represent spatially complete snapshots of the solution at defined instances during the simulation.  Dependeing on the control parameters provided here, visualizatoin files may include only a fraction of the state data, and may contiain auxiliary "derived" information (see below for more discussion).

* `"Observation Data`" is intended to represent diagnostic values to be returned to the calling routine from Amanzi's simulation driver.  Observations are typically generated at arbitrary times, and frequently involve various point samplings and volumetric reductions that are interpolated in time to the desired instant.  Observations may involve derived quantities (see discussion below) or state fields.

* `"Log Data`" is intended to represent runtime diagnostics to indicate the status of the simulation in progress.  This data is typically written by the simulation code to the screen or some other stream or file pipe.  The volume of `"Log Data`" generated is a function of the `"Verbosity`" setting under `"Execution Control`".

* `"Walkabout Data`" is intended to be used as input to the particle tracking software Walkabout.

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
 * Gravimetric water content [volumetric water content * water density / bulk density, in kg/m^3]
 * Hydraulic Head [ (aqueous pressure - atmospheric pressure)/(rho * gravity) + z ]
 * Drawdown [m] (difference of hydraulic heads at different time moments) 
 * Aqueous mass flow rate [ kg/s ] (must use integral functional in the observation)
 * Aqueous volumetric flow rate [ m^3/s ] (must use integral functional in the observation)

Note that MaterialID will be treated as a double that is unique to each defined material.  Its value will be generated internal to Amanzi.  The log file will be appended with the (material name)->(integer) mapping used.  Also note that this list tacitly assumes the presence of Aqueous Water as one of the transported components.  Presently, it is an error if the `"Phase Definition`" above does not sufficiently define this component.


Time macros specify a rule to generate a list of time values.  They are defined in the parameter list `"Time Macros`":

* [SU] `"Time Macros`" [list] can accept multiple lists for user-named macros TMACRO

 * [SU] TMACRO [list] can accept either `"Values`" or `"Start_Period_Stop`"

  * [SU] `"Values`" [Array(double)] values of time, or 

  * [SU] `"Start_Period_Stop`" [Array(double)] values of start time (ts), period (dt) and (optionally) end time (te) to generate times, t=ts + dt*i, for any integer i. If stop time is less than start time, the time intervals have no ending.


Cycle macros specify a rule to generate or list cycle values.  They are defined in the parameter list `"Cycle Macros`":

* [SU] `"Cycle Macros`" [list] can accept multiple lists for user-named macros CMACRO

 * [SU] CMACRO [list] can accept either `"Values`" or `"Start_Period_Stop`"

  * [SU] `"Values`" [Array(int)] values of cycle number, or 

  * [SU] `"Start_Period_Stop`" [Array(int)] values of start cycle (cs), period (dc) and (optionally) end cycle (ce) to generate cycle numbers, c=cs + dc*i, for any integer i. If stop cycle < 0, the cycle intervals will not end.



Observation Data
----------------

A user may request any number of specific observations from Amanzi.  Each labeled Observation Data quantity involves a field quantity, a model, a region from which it will extract its source data, and a list of discrete times 
for its evaluation.  The observations are evaluated during the simulation and returned to the calling process through one of Amanzi arguments.

* [SU] `"Observation Data`" [list] can accept multiple lists for named observations (OBSERVATION)

  * [SU] `"Observation Output Filename`" [string] user-defined name for the file that the observations are written to.

  * [SU] OBSERVATION [list] user-defined label, can accept values for `"Variable`", `"Functional`", `"Region`", `"Time Macro`", and `"Cycle Macro`".

    * [SU] `"Variable`" [string] name of field quantities taken from the list of "Available field quantities" defined under "Time and Cycle specification"

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
      <Parameter name="Start_Period_Stop" type="Array(double)" value="{0, 3.1536e7,-1}"/>
    </ParameterList>
  </ParameterList>

  <ParameterList name="Observation Data">
    <Parameter name="Observation Output Filename" type="string" value="obs_output.out"/>
    <ParameterList name="Integrated Mass">
      <Parameter name="Region" type="string" value="All"/>
      <Parameter name="Functional" type="string" value="Observation Data: Integral"/>
      <Parameter name="Variable" type="string" value="Volumetric Water Content"/>
      <Parameter name="Time Macro" type="string" value="Annual"/>
    </ParameterList>
  </ParameterList>

In this example, the user requests an annual report of the integrated volume of water over the entire domain.



Checkpoint Data
---------------------------------

A user may request periodic dumps of Amanzi Checkpoint Data.  The user has no explicit control over the content of these files, but has the guarantee that the Amanzi run will be reproducible (with accuracies determined
by machine round errors and randomness due to execution in a parallel computing environment).  Therefore, output controls for Checkpoint Data are limited to file name generation and writing frequency, by numerical cycle number.

* [SU] `"Checkpoint Data`" [list] can accept a file name base [string] and cycle data [list] used to generate the file base name or directory base name that is used in writing Checkpoint Data. 

  * [SU] `"File Name Base`" [string]

  * [SU] `"File Name Digits`" [int] specify the number of digits that should be appended to the file name for the cycle number.

  * [SU] `"Cycle Macro`" [string] can accept label of user-defined Cycle Macro (see above)

Notes:

[SU] if only `"File Name Base`" and optionally `"File Name Digits`" are specified then only one checkpoint is
written at the end of the simulation.



Example:

.. code-block:: xml

  <ParameterList name="Cycle Macros">
    <ParameterList name="Every-5">
      <Parameter name="Start_Period_Stop" type="Array(int)" value="{0, 5, -1}"/>
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
  
  * [SU] `"Cycle Macros`" [Array(string)] can accept a list of of user-defined Cycle Macro (see above)
  
  * [SU] `"Time Macros`" [Array(string)] can accept a list of the labeled time macros (see above)

  * [S] `"Variables`" [Array(string)] can accept a list of field quantities to include in the file.  At present the unstructured code dumps all of the dependent variables in the system state.

  * [U] `"Write Regions`" [Array(string)] (empty array) write an array into the visualization file that can be used to identify a region or regions. The first entry in the regions array is marked with the value 1.0 in the array, the second with the value 2.0, and so forth. The code ignores entries in the regions array that are not valid regions that contain cells.

  * [U] `"Write Partitions`" [bool] (false) if this parameter is true, then write an array into the visualization file that contains the rank number of the processor that owns a mesh cell. 



Notes:

[SU] If only `"File Name Base`" is specified, an initial and a final visualization dump are written.



Example:

.. code-block:: xml

  <ParameterList name="Cycle Macros">
    <ParameterList name="Every-10">
      <Parameter name="Start_Period_Stop" type="Array(int)" value="{0, 10,-1}"/>
    </ParameterList>
  </ParameterList>

  <ParameterList name="Visualization Data">
    <Parameter name="File Name Base" type="string" value="chk"/>
    <Parameter name="File Name Digits" type="int" value="5"/>
    <Parameter name="Cycle Macros" type="Array(string)" value="Every-10">
  </ParameterList>

In this example, the liquid pressure and moisture content are written when the cycle number is evenly divisble by 5.




Walkabout Data
--------------

A user may request periodic dumps of Walkabout Data.  Output controls for Walkabout Data are limited to file name generation and writing frequency, by numerical cycle number.

* [U] `"Walkabout Data`" [list] can accept a file name base [string] and cycle data [list] used to generate the file base name or directory base name that is used in writing Walkabout Data. 

  * [U] `"File Name Base`" [string]

  * [U] `"File Name Digits`" [int] specify the number of digits that should be appended to the file name for the cycle number.

  * [U] `"Cycle Macro`" [string] can accept label of user-defined Cycle Macro (see above)

Notes:

[U] if only `"File Name Base`" and optionally `"File Name Digits`" are specified then only one walkabout dump is
written at the end of the simulation.



Example:

.. code-block:: xml

  <ParameterList name="Cycle Macros">
    <ParameterList name="Every-5">
      <Parameter name="Start_Period_Stop" type="Array(int)" value="{0, 5, -1}"/>
    </ParameterList>
  </ParameterList>

  <ParameterList name="Walkabout Data">
    <Parameter name="File Name Base" type="string" value="walk"/>
    <Parameter name="File Name Digits" type="int" value="5"/>
    <Parameter name="Cycle Macro" type="string" value="Every-5"/>
  </ParameterList>

In this example, Checkpoint Data files are written when the cycle number is evenly divisible by 5.




Restart
-------

A user may request a restart from a Checkpoint Data file by including the sublist 
`"Restart`" in the Execution Control list. This mode of restarting
will overwrite all other initializations of data that are called out in the input file.
The purpose of restarting Amanzi in this fashion is mostly to continue a run that has been 
terminated because its allocation of time ran out.


* [S] `"Restart`" [list]

  * [S] `"File Name`" [string] file name of the specific Checkpoint Data file to restart from

Example:

.. code-block:: xml

  <ParameterList name="Restart">
     <Parameter name="File Name" type="string" value="chk00123.h5"/>
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


Complete Example
=================

Presented below is a complete example of an Amanzi input file.  It does not exercise all the options provided for in this specification, but rather provides a concrete example of a set of self-consistent definitions
required to specify a real simulation with Amanzi envisioned functional for the Phase 2 demo deadline.

.. code-block:: xml


       <?xml version="1.0" encoding="utf-8"?>
       <!-- The input example below conforms to the current (10/4/11) input specification on 
            the Amanzi wiki.  We believe that this specification will need to be revised to better
            reflect terminology and organization used by domain scientists. GEH, VLF, MLR -->
       <!--
          BC Cribs input spec for 3D PE Analysis which involves variably saturated flow, 
          solute transport, and homogeneous property distributions within material types.
          This input file depends on a checkpoint (restart) file from a previous steady-state
          flow simulation. 
          
          Submitted by Vicky Freedman, PNNL, September 26, 2011.
          Revised by Glenn Hammond, PNNL, September 28, 2011
          
         -->
       
       <!-- GEH: Assumptions:
       
          GEH: Glenn Hammond
          VLF: Vicky Freedman
          MLR: Mark Rockhold
       
          Units: All units assumed to be SI.
          Domain: 400. x 1. x 110. meters in x, y, z (the domain width could decrease!)
       
              Bear in mind that this is hypothetical and does not reflect that actual lithofacies at
              HDVZ.  Characters { ,*,#} indicate materials.
       
                  <-                    400m                 ->
                               Crib 1        Crib 2  
                   ___________xxxxxxxx______xxxxxxxx___________ <-BC: flow = Neumann, transport = inflow 
                  |                                            |   ^
                  |                        *                   |   |
                  |                        **                  |        
                  |          *             ***       ****      |
                  |   ***  ****       #  ********* ****        |
                  |    ***** ****      ##       *****          |
                  |      ***********    ####      **           |  110m
                  |         ******       ######                |
        BC: flow->|           ***          ########            |<-BC: flow = zero flow, transport = zero flux 
                  |            *****         ##########        |
                  |              **            ############    |
                  |                               ###########  |   |
                  |____________________________________________|   v
                                                                <-BC: flow = Dirichlet pressure, transport = outflow 
       
          Simplifications:
            1. Dispersion removed
            2. Diffusion removed
            3. Zero gradient boundary condition will be assumed for solute outflow
            4. Assuming no geochemistry, solely solute transport
            5. Facies are layered instead of irregular/heterogeneous
       
       -->
       
       <ParameterList name="Main">
       
         <Parameter name="Amanzi input format version" type="string" value="0.9.2"/>
       
         <ParameterList name="General Description">
           <Parameter name="Model ID" type="string" value="Transient Richards"/>
           <Parameter name="Model name" type="string" value="BC Cribs PE Template"/>
           <Parameter name="Description" type="string" value="Unsat flow and transport"/>
           <Parameter name="Purpose" type="string" value="Provide input req. for Phase II Demo"/>
           <Parameter name="Creation date" type="string" value="09.25.11 01:28"/>
           <Parameter name="Last modified" type="string" value="09.25.11 01:28"/>
         </ParameterList>
       

         <ParameterList name="Execution Control">

           <Parameter name="Flow Mode" type="string" value="Richards"/>
           <Parameter name="Transport Mode" type="string" value="On"/>
           <Parameter name="Chemistry Mode" type="string" value="Off"/>

           <ParameterList name="Time">
             <Parameter name="Integration Mode" type="string" value="Transient"/>
             <!-- 1956 -->
             <Parameter name="Start" type="double" value="0"/>
             <!-- 2006 -->
             <Parameter name="End" type="double" value="1.5768e9"/>
           </ParameterList>

           <Parameter name="Verbosity" type="string" type="High"/>

           <ParameterList name="Numerical Control Parameters">
             <ParameterList name="Unstructured Algorithm">
               <ParameterList name="Adaptive Mesh Refinement">
                 <Parameter name="Transport Integration Algorithm" type="string" value="None"/>
               </ParameterList>
             </ParameterList>
           </ParameterList>
         </ParameterList>


         <ParameterList name="Domain">
           <Parameter name="Spatial Dimension" type="int" value="2"/>
         </ParameterList>
         <ParameterList name="Mesh">
         <!-- Uncomment this block for unstructured with an internally generated mesh
           <ParameterList name="Unstructured">
             <ParameterList name="Generate Mesh">
               <ParameterList name="Uniform Structured">
                 <Parameter name="Number of Cells" type="Array(int)" value="{800, 1, 220}"/>
                 <Parameter name="Domain Low Coordinate" type="Array(double)" value="{0.0, 0.0, 0.0}" />
                 <Parameter name="Domain High Coordinate" type="Array(double)" value="{400., 1.0, 110.}" />
               </ParameterList>
             </ParameterList>
           </ParameterList>
         -->
         <!-- Uncomment this block for unstructured with an externally generated mesh
           <ParameterList name="Unstructured">
             <ParameterList name="Read Mesh File">
               <Parameter name="File" type="string" value="mesh_filename"/>
               <Parameter name="Format" type="string" value="Exodus II"/>
             </ParameterList>
           </ParameterList>
         -->
         <!-- Uncomment this block for structured
           <ParameterList name="Structured">
             <Parameter name="Number of Cells" type="Array(int)" value="{800, 1, 220}"/>
             <Parameter name="Domain Low Coordinate" type="Array(double)" value="{0.0, 0.0, 0.0}" />
             <Parameter name="Domain High Coordinate" type="Array(double)" value="{400., 1.0, 110.}" />
           </ParameterList>
         -->
         </ParameterList>
       
         <ParameterList name="Regions">
       
         <!--
                   ____________________________________________  110.
                  |                                            |
                  |     Material 3 Region                      |
                  |                                            |
                  |____________________________________________| 60.
                  |                                            |
                  |     Material 2 Region                      |
                  |____________________________________________| 30.
                  |                                            |
                  |     Material 1 Region                      |
                  |____________________________________________| 0.
       
           -->
       
           <ParameterList name="Material 1 Region">
             <ParameterList name="Region: Box">
               <Parameter name="Low Coordinate" type="Array(double)" value="{0.0, 0.0, 0.0}"/>
               <Parameter name="High Coordinate" type="Array(double)" value="{400.0, 1.0, 30.0}"/>
             </ParameterList>
           </ParameterList>
       
           <ParameterList name="Material 2 Region">
             <ParameterList name="Region: Box">
               <Parameter name="Low Coordinate" type="Array(double)" value="{0.0, 0.0, 30.0}"/>
               <Parameter name="High Coordinate" type="Array(double)" value="{400.0, 1.0, 60.0}"/>
             </ParameterList>
           </ParameterList>
       
           <ParameterList name="Material 3 Region">
             <ParameterList name="Region: Box">
               <Parameter name="Low Coordinate" type="Array(double)" value="{0.0, 0.0, 60.0}"/>
               <Parameter name="High Coordinate" type="Array(double)" value="{400.0, 1.0, 110.0}"/>
             </ParameterList>
           </ParameterList>
       
             <!-- GEH:  
                               Crib 1        Crib 2  
                   ___________xxxxxxxx______xxxxxxxx___________ 
                  |                                            |
       
           -->
           <ParameterList name="Top Surface Outside Cribs Region A">
             <ParameterList name="Box">
               <!-- GEH: These are approximate as placeholders for how.  Vicky will provide more
                         accurate values soon. -->
               <Parameter name="Low Coordinate" type="Array(double)" value="{0.0, 0.0, 110.0}"/>
               <Parameter name="High Coordinate" type="Array(double)" value="{170.0, 1.0, 110.0}"/>
             </ParameterList>
           </ParameterList>
           <ParameterList name="Top Surface Outside Cribs Region B">
             <ParameterList name="Region: Box">
               <Parameter name="Low Coordinate" type="Array(double)" value="{173.0, 0.0, 110.0}"/>
               <Parameter name="High Coordinate" type="Array(double)" value="{190.0, 1.0, 110.0}"/>
             </ParameterList>
           </ParameterList>
           <ParameterList name="Top Surface Outside Cribs Region C">
             <ParameterList name="Region: Box">
               <Parameter name="Low Coordinate" type="Array(double)" value="{193.0, 0.0, 110.0}"/>
               <Parameter name="High Coordinate" type="Array(double)" value="{400.0, 1.0, 110.0}"/>
             </ParameterList>
           </ParameterList>
       
           <ParameterList name="90 Meter Plane Region">
             <!-- GEH: Note that we could use a 2D box for these regions too. -->
             <ParameterList name="Region: Plane">
               <Parameter name="Coordinate"  type="Array(double)" value="{0., 0., 90.}"/>
               <!-- GEH: Note the downward unit vector -->
               <Parameter name="Direction"  type="Array(double)" value="{0., 0., 1.}"/>
             </ParameterList>
           </ParameterList>
       
           <ParameterList name="Bottom Surface Region">
             <!-- GEH: Note that we could use a 2D box for these regions too. -->
             <ParameterList name="Region: Plane">
               <Parameter name="Coordinate"  type="Array(double)" value="{0., 0., 0.}"/>
               <!-- GEH: Note the downward unit vector -->
               <Parameter name="Direction"  type="Array(double)" value="{0., 0., -1.}"/>
             </ParameterList>
           </ParameterList>
       
           <ParameterList name="West Surface Region">
             <ParameterList name="Region: Plane">
               <Parameter name="Coordinate"  type="Array(double)" value="{0., 0., 0.}"/>
               <Parameter name="Direction"  type="Array(double)" value="{-1., 0., 0.}"/>
             </ParameterList>
           </ParameterList>
       
           <ParameterList name="East Surface Region">
             <ParameterList name="Region: Plane">
               <Parameter name="Coordinate"  type="Array(double)" value="{400., 0., 0.}"/>
               <Parameter name="Direction"  type="Array(double)" value="{1., 0., 0.}"/>
             </ParameterList>
           </ParameterList>
       
           <ParameterList name="South Surface Region">
             <ParameterList name="Region: Plane">
               <Parameter name="Coordinate"  type="Array(double)" value="{0., 0., 0.}"/>
               <Parameter name="Direction"  type="Array(double)" value="{0., -1., 0.}"/>
             </ParameterList>
           </ParameterList>
       
           <ParameterList name="North Surface Region">
             <ParameterList name="Region: Plane">
               <Parameter name="Coordinate"  type="Array(double)" value="{0., 1., 0.}"/>
               <Parameter name="Direction"  type="Array(double)" value="{0., 1., 0.}"/>
             </ParameterList>
           </ParameterList>
       
           <ParameterList name="Crib 1 Region">
             <ParameterList name="Region: Box">
               <!-- GEH: Assuming unit cell width in Y -->
               <Parameter name="Low Coordinate" type="Array(double)" value="{170.0, 0.0, 110.0}"/>
               <Parameter name="High Coordinate" type="Array(double)" value="{173.0, 1.0, 110.0}"/>
             </ParameterList>
           </ParameterList>
       
           <ParameterList name="Crib 2 Region">
             <ParameterList name="Region: Box">
               <!-- GEH: Assuming unit cell width in Y -->
               <Parameter name="Low Coordinate" type="Array(double)" value="{190.0, 0.0, 110.0}"/>
               <Parameter name="High Coordinate" type="Array(double)" value="{193.0, 1.0, 110.0}"/>
             </ParameterList>
           </ParameterList>
       
           <ParameterList name="Sample Point 1 Region">
             <ParameterList name="Region: Point">
               <Parameter name="Coordinate"  type="Array(double)" value="{171.5, 0.5, 50.0}"/>
             </ParameterList>
           </ParameterList>
       
           <ParameterList name="Sample Point 2 Region">
             <ParameterList name="Region: Point">
               <Parameter name="Coordinate"  type="Array(double)" value="{191.5, 0.5, 50.0}"/>
             </ParameterList>
           </ParameterList>
       
           <ParameterList name="Sample Point 3 Region">
             <ParameterList name="Region: Point">
               <Parameter name="Coordinate"  type="Array(double)" value="{181.5, 0.5, 50.0}"/>
             </ParameterList>
           </ParameterList>
       
         </ParameterList>
       
         <ParameterList name="Material Properties">
       
           <ParameterList name="Material 1">
       
             <ParameterList name="Porosity: Uniform">
               <Parameter name="Value" type="double" value="0.38"/>
             </ParameterList>
       
             <ParameterList name="Intrinsic Permeability: Anisotropic Uniform">
               <Parameter name="Horizontal" type="double" value="2.05e-8"/>
               <Parameter name="Vertical" type="double" value="2.05e-9"/>
             </ParameterList>
       
             <!-- GEH: Pressure-saturation function: van Genuchten-->
             <ParameterList name="Capillary Pressure: van Genuchten">
               <Parameter name="alpha" type="double" value="2.14e-4"/> <!-- 0.021 cm^-1 -> Pa^-1 -->
               <Parameter name="Sr" type="double" value="0.0"/>
               <Parameter name="m" type="double" value="0.601"/>
               <Parameter name="Relative Permeability" type="string" value="Mualem"/>
             </ParameterList>
       
             <Parameter name="Assigned Regions" type="Array(string)" value="{Material 1 Region}"/>

           </ParameterList>
       
           <ParameterList name="Material 2">
       
             <ParameterList name="Porosity: Uniform">
               <Parameter name="Value" type="double" value="0.36"/>
             </ParameterList>
       
             <ParameterList name="Intrinsic Permeability: Anisotropic Uniform">
               <Parameter name="Horizontal" type="double" value="4.84e-8"/>
               <Parameter name="Vertical" type="double" value="4.84e-9"/>
             </ParameterList>
       
             <ParameterList name="Capillary Pressure: van Genuchten">
               <Parameter name="alpha" type="double" value="7.35e-4"/>
               <Parameter name="Sr" type="double" value="0.0"/>
               <Parameter name="m" type="double" value="0.511"/>
               <Parameter name="Relative Permeability" type="string" value="Mualem"/>
             </ParameterList>
       
             <Parameter name="Assigned Regions" type="Array(string)" value="{Material 2 Region}"/>

           </ParameterList>
       
           <ParameterList name="Material 3">
       
             <ParameterList name="Porosity: Uniform">
               <Parameter name="Value" type="double" value="0.23"/>
             </ParameterList>
       
             <ParameterList name="Intrinsic Permeability: Anisotropic Uniform">
               <Parameter name="Horizontal" type="double" value="3.00e-9"/>
               <Parameter name="Vertical" type="double" value="3.00e-10"/>
             </ParameterList>

             <ParameterList name="Capillary Pressure: van Genuchten">
               <Parameter name="alpha" type="double" value="1.74e-4"/>
               <Parameter name="Sr" type="double" value="0.0"/>
               <Parameter name="m" type="double" value="0.420"/>
               <Parameter name="Relative Permeability" type="string" value="Mualem"/>
             </ParameterList>
       
             <Parameter name="Assigned Regions" type="Array(string)" value="{Material 3 Region}"/>

           </ParameterList>
       
         </ParameterList>
       
         <ParameterList name="Phase Definitions">
           <ParameterList name="Aqueous">
             <ParameterList name="Phase Properties">
               <ParameterList name="Viscosity: Uniform">
                 <Parameter name="Viscosity" type="double" value="8.9e-4"/>
               </ParameterList>
               <ParameterList name="Density: Uniform">
                 <Parameter name="Density" type="double" value="998.32"/>
               </ParameterList>
             </ParameterList>
             <ParameterList name="Phase Components">
               <!-- GEH: Note sure if this is what we want.  Water component with solutes.  The input spec
                         reflects this, although it refers to "Aqueous Water" instead of "Water". -->
               <ParameterList name="Water">
                 <Parameter name="Component Solutes" type="Array(string)" value="{Tc-99}"/>           
               </ParameterList>
             </ParameterList>
           </ParameterList>
         </ParameterList>


         <ParameterList name="Initial Conditions">
           <ParameterList name="IC For Domain">
             <Parameter name="Assigned Regions" type="Array(string)" value="{All}"/>
             <ParameterList name="IC: Linear Pressure">
               <Parameter name="Phase" type="string" value="Aqueous"/>
               <Parameter name="Reference Value" type="double" value="101325."/>
               <Parameter name="Reference Point" type="Array(double)" value="{0., 0., 0.}"/>
               <!-- GEH: Units of gradient are Pa/m = rho*g = 998.32 kg/m^3 * 9.81 m/s^2-->
               <Parameter name="Gradient Value" type="Array(double)" value="{0., 0., -9793.5192}"/>
             </ParameterList>
             <ParameterList name="Solute IC">
               <ParameterList name="Aqueous">
                 <ParameterList name="Water">
                   <ParameterList name="Tc-99">
                     <ParameterList name="IC: Uniform">
                       <Parameter name="Value" type="double" value="0.0"/>
                     </ParameterList>
                   </ParameterList>     
                 </ParameterList>
               </ParameterList>
             </ParameterList>
           </ParameterList>
         </ParameterList>



         <ParameterList name="Boundary Conditions">
           <ParameterList name="BC For Top Surface Outside Cribs Region">
             <Parameter name="Assigned Regions" type="Array(string)" value="{Top Surface Outside Cribs Region A, Top Surface Outside Cribs Region B, Top Surface Outside Cribs Region C}"/>
             <ParameterList name="BC: Flux">
	       <!-- GEH/VLF: These recharge intervals/rates will change. -->
               <!-- 1956, 1984 in seconds-->
               <Parameter name="Times" type="Array(double)" value="{0., 883008000.}"/>
               <Parameter name="Time Functions" type="Array(string)" value="{Constant, Constant}"/>
               <!-- Recharge = 77 mm/yr, 25mm/yr -->
               <Parameter name="Extensive Flux" type="Array(double)" value="{2.44e-9, 7.93e-10}"/>
             </ParameterList>
             <ParameterList name="Solute BC">
               <ParameterList name="Aqueous">
                 <ParameterList name="Water">
                   <ParameterList name="Tc-99">
                     <ParameterList name="BC: Inflow">
                         <!-- GEH: Throughout entire simulation, no solute enters through top surface -->
                       <Parameter name="Times" type="Array(double)" value="{0.}"/>
                       <Parameter name="Time functions" type="Array(string)" value="{Constant}"/>
                       <Parameter name="Values" type="Array(double)" value="{0.}"/>
                     </ParameterList>
                   </ParameterList>
                 </ParameterList>
               </ParameterList>
             </ParameterList>
           </ParameterList>

           <ParameterList name="BC For Crib 1 Region">
             <Parameter name="Assigned Regions" type="Array(string)" value="{Crib 1 Region}"/>
             <ParameterList name="BC: Flux">
                 <!-- GEH/VLF: These recharge intervals/rates will change. -->
                 <!-- 1956, 1956.25 in seconds-->
               <Parameter name="Times" type="Array(double)" value="{0., 7884000.}"/>
               <Parameter name="Time functions" type="Array(string)" value="{Constant, Constant}"/>
                 <!-- 11.25, 0. m/d-->
               <Parameter name="Extensive Flux" type="Array(double)" value="{1.302e-4, 0.}"/>
             </ParameterList>
             <ParameterList name="Solute BC">
               <ParameterList name="Aqueous">
                 <ParameterList name="Water">
                   <ParameterList name="Tc-99">
                     <ParameterList name="BC: Inflow">
                        <!-- 1956, 1956.25 in seconds-->
                       <Parameter name="Times" type="Array(double)" value="{0., 7884000.}"/>
                       <Parameter name="Time functions" type="Array(string)" value="{Constant, Constant}"/>
                       <Parameter name="Values" type="Array(double)" value="{1000., 0.}"/>
                     </ParameterList>
                   </ParameterList>
                 </ParameterList>
               </ParameterList>
             </ParameterList>
           </ParameterList>

           <ParameterList name="BC For Crib 2 Region">
             <Parameter name="Assigned Regions" type="Array(string)" value="{Crib 2 Region}"/>
             <ParameterList name="BC: Flux">
                 <!-- GEH/VLF: These recharge intervals/rates will change. -->
                 <!-- 1956, 1956.33, 1956.66 in seconds-->
               <Parameter name="Times" type="Array(double)" value="{0., 10406880., 20813760.}"/>
               <Parameter name="Time functions" type="Array(string)" value="{Constant, Constant, Constant}"/>
                   <!-- 0., 8.75, 0. m/d-->
               <Parameter name="Extensive Flux" type="Array(double)" value="{0., 1.013e-4, 0.}"/>
             </ParameterList>
             <ParameterList name="Solute BC">
               <ParameterList name="Aqueous">
                 <ParameterList name="Water">
                   <ParameterList name="Tc-99">
                     <ParameterList name="BC: Inflow">
                       <!-- 1956, 1956.33, 1956.66 in seconds-->
                       <Parameter name="Times" type="Array(double)" value="{0., 10406880., 20813760.}"/>
                       <Parameter name="Time functions" type="Array(string)" value="{Constant, Constant, Constant}"/>
                       <Parameter name="Values" type="Array(double)" value="{0., 900., 0.}"/>
                     </ParameterList>
                   </ParameterList>
                 </ParameterList>
               </ParameterList>
             </ParameterList>
           </ParameterList>
       
           <ParameterList name="BC For Bottom Surface Region">
             <Parameter name="Assigned Regions" type="Array(string)" value="{Bottom Surface Region}"/>
             <ParameterList name="BC: Uniform Pressure">
               <Parameter name="Times" type="Array(double)" value="{0.}"/>
               <Parameter name="Time functions" type="Array(string)" value="{Constant}"/>
               <Parameter name="Values" type="Array(double)" value="{101325.}"/>
             </ParameterList>
             <ParameterList name="Solute BC">
               <ParameterList name="Aqueous">
                 <ParameterList name="Water">
                   <ParameterList name="Tc-99">
                     <ParameterList name="BC: Outflow">
                       <Parameter name="Times" type="Array(double)" value="{0.}"/>
                       <Parameter name="Time functions" type="Array(string)" value="{Constant}"/>
                     </ParameterList>
                   </ParameterList>
                 </ParameterList>
               </ParameterList>
             </ParameterList>
           </ParameterList>
       
           <ParameterList name="BC For West Surface Region">
             <Parameter name="Assigned Regions" type="Array(string)" value="{West Surface Region}"/>
             <ParameterList name="BC: No Flow">
               <Parameter name="Times" type="Array(double)" value="{0.}"/>
               <Parameter name="Time functions" type="Array(string)" value="{Constant}"/>
             </ParameterList>
             <ParameterList name="Solute BC">
               <ParameterList name="Aqueous">
                 <ParameterList name="Water">
                   <ParameterList name="Tc-99">
                     <ParameterList name="BC: Zero Flux">
                       <Parameter name="Times" type="Array(double)" value="{0.}"/>
                       <Parameter name="Time functions" type="Array(string)" value="{Constant}"/>
                     </ParameterList>
                   </ParameterList>
                 </ParameterList>
               </ParameterList>
             </ParameterList>
           </ParameterList>

           <ParameterList name="BC For East Surface Region">
             <Parameter name="Assigned Regions" type="Array(string)" value="{East Surface Region}"/>
             <ParameterList name="BC: No Flow">
               <Parameter name="Times" type="Array(double)" value="{0.}"/>
               <Parameter name="Time functions" type="Array(string)" value="{Constant}"/>
             </ParameterList>
             <ParameterList name="Solute BC">
               <ParameterList name="Aqueous">
                 <ParameterList name="Water">
                   <ParameterList name="Tc-99">
                     <ParameterList name="BC: Zero Flux">
                       <Parameter name="Times" type="Array(double)" value="{0.}"/>
                       <Parameter name="Time functions" type="Array(string)" value="{Constant}"/>
                     </ParameterList>
                   </ParameterList>
                 </ParameterList>
               </ParameterList>
             </ParameterList>
           </ParameterList>

         </ParameterList>
       
         <ParameterList name="Output">

           <!-- GEH: The following are desired for output:
             1. Integrated water and Tc-99 mass over time (yearly) (mass balance)
             2. Water saturation, water pressure and Tc-99 concentration throughout space at specified times (plot file)
             3. Water saturation, water pressure and Tc-99 concentration over time at points in space (breakthrough)
             4. Integrate Tc-99 mass crossing the cribs and bottom boundaries over time (flux)
             5. Checkpoint files every N time steps
       
             I will attempt these calculations based on the "Observation Data" section of the input spec.
             -->


           <!-- define some handy cycle macros -->
           <ParameterList name="Cycle Macros">
             <ParameterList name="Every-5-steps">
               <Parameter name="Start_Period_Stop" type="Array(int)" value="{0, 5, -1}"/>
             </ParameterList>
           </ParameterList>

           <!-- define some handy time macros -->
           <ParameterList name="Time Macros">
             <ParameterList name="Annual">
               <Parameter name="Start_Period_Stop" type="Array(double)" value="{0, 3.1536e7, -1}"/>
             </ParameterList>

             <ParameterList name="My_times">
               <!-- 1956, 1956.1, 1956.2, 1956.3, 1956.4, 1956.4, 1956.5, 1956.6, 1956.7, 1956.8, 1956.9, 1957, 1958, 1960, 1970, 1980, 1990, 2000, 2006 -->
               <Parameter name="Values" type="Array(double)" value="{0., 3153600., 6307200., 9460800., 12614400., 1576800., 18921600., 22075200., 25228800., 28382400., 31536000., 63072000., 126144000., 441504000., 756864000., 1072224000., 1387584000., 1576800000. }"/>
             </ParameterList>

             <ParameterList name="Daily_1957-1967">
               <Parameter name="Start_Period_Stop" type="Array(double)" value="{3.1536e7, 86400., 3.46896e8}"/>
             </ParameterList>

             <ParameterList name="Daily_1957-2006">
               <Parameter name="Start_Period_Stop" type="Array(double)" value="{3.1536e7, 86400., 1.5768e9}"/>
             </ParameterList>

           </ParameterList>


           <ParameterList name="Observation Data">

             <!-- Global water and Tc-99 mass -->
             <ParameterList name="Integrated Mass">
               <Parameter name="Region" type="string" value="All"/>
               <Parameter name="Functional" type="string" value="Observation Data: Integral"/>
               <Parameter name="Variables" type="Array(string)" value="{Water Mass Density, Tc-99 Aqueous Concentration}"/>
               <Parameter name="Time Macro" type="string" value="Annual"/>
             </ParameterList>

             <!-- Point samples of water and Tc-99 -->
             <ParameterList name="Point Sample 1">
               <Parameter name="Region" type="string" value="Sample Point 1 Region"/>
               <Parameter name="Functional" type="string" value="Observation Data: Point"/>
               <Parameter name="Variables" type="Array(string)" value="{Water Mass Density, Tc-99 Aqueous Concentration}"/>
               <Parameter name="Time Macro" type="string" value="Daily_1957-1967"/>
             </ParameterList>

             <ParameterList name="Point Sample 2">
               <Parameter name="Region" type="string" value="Sample Point 2 Region"/>
               <Parameter name="Functional" type="string" value="Observation Data: Point"/>
               <Parameter name="Variables" type="Array(string)" value="{Water Mass Density, Tc-99 Aqueous Concentration}"/>
               <Parameter name="Time Macro" type="string" value="Daily_1957-1967"/>
             </ParameterList>

             <ParameterList name="Point Sample 3">
               <Parameter name="Region" type="string" value="Sample Point 3 Region"/>
               <Parameter name="Functional" type="string" value="Observation Data: Point"/>
               <Parameter name="Variables" type="Array(string)" value="{Water Mass Density, Tc-99 Aqueous Concentration}"/>
               <Parameter name="Time Macro" type="string" value="Daily_1957-1967"/>
             </ParameterList>

             <!-- cummulative flux of Tc-99 -->
             <ParameterList name="Cummulative Tc-99 Flux Integral - Bottom">
               <Parameter name="Region" type="string" value="Bottom Surface Region"/>
               <Parameter name="Functional" type="string" value="Observation Data: Cummulative Integral"/>
               <Parameter name="Variables" type="Array(string)" value="{Tc-99 Aqueous Concentration}"/>
               <Parameter name="Time Macro" type="string" value="Daily_1957-2006"/>
             </ParameterList>

             <ParameterList name="Cummulative Tc-99 Flux Integral - Crib 1">
               <Parameter name="Region" type="string" value="Crib 1 Region"/>
               <Parameter name="Functional" type="string" value="Observation Data: Cummulative Integral"/>
               <Parameter name="Variables" type="Array(string)" value="{Tc-99 Aqueous Concentration}"/>
               <Parameter name="Time Macro" type="string" value="Daily_1957-2006"/>
             </ParameterList>

             <ParameterList name="Cummulative Tc-99 Flux Integral - Crib 2">
               <Parameter name="Region" type="string" value="Crib 2 Region"/>
               <Parameter name="Functional" type="string" value="Observation Data: Cummulative Integral"/>
               <Parameter name="Variables" type="Array(string)" value="{Tc-99 Aqueous Concentration}"/>
               <Parameter name="Time Macro" type="string" value="Daily_1957-2006"/>
             </ParameterList>

             <ParameterList name="Cummulative Tc-99 Flux Integral - 90m">
               <Parameter name="Region" type="string" value="90 Meter Plane Region"/>
               <Parameter name="Functional" type="string" value="Observation Data: Cummulative Integral"/>
               <Parameter name="Variables" type="Array(string)" value="{Tc-99 Aqueous Concentration}"/>
               <Parameter name="Time Macro" type="string" value="Daily_1957-2006"/>
             </ParameterList>

           </ParameterList>
       

           <ParameterList name="Visualization Data">
             <Parameter name="File Name Base" type="string" value="viz-"/>
             <Parameter name="Cycle Macro" type="string" value="Every-10-steps"/>
             <Parameter name="Variables" type="Array(string)" value="{Aqueous Pressure, Tc-99 Aqueous Concentration}"/>
           </ParameterList>

           <ParameterList name="Checkpoint Data">
             <Parameter name="File Name Base" type="string" value="dump-"/>
             <Parameter name="Cycle Macro" type="string" value="Every-100-steps"/>
           </ParameterList>

         </ParameterList> <!-- End of Output -->
       
       </ParameterList> <!-- End of Main -->
       
