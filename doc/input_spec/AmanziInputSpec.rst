========================================
Amanzi XML Input Specification
========================================

The ASCEM simulator, Amanzi, evolves a system of conservation
equations for reacting flow in porous media.  The model is detailed in
the ASCEM report entitled "Mathematical Formulation Requirements and
Specifications for the Process Models`", which we refer to herein as
the 'Model Requirements Document (MRD)'). The purpose of the present
document is to specify the type and format of data required to execute
Amanzi on a problem within its intended scope.  This specification
should be regarded as a companion to the MRD; it relies heavily on
the detailed formulations of the models.  Where applicable, the
relevant sections of the MRD are indicated.


.. contents:: **Table of Contents**


========================================
Overview
========================================

For each Amanzi simulation, a user must specify a set of phase and
tracer components they wish to follow, and provide initial data and
boundary conditions for all affected states in terms of a set of
parameterized functionals.  Conservation of mass for each of the
specified components is given by equation 2.11 of the MRD, where the
volumetric flow rate has been specified via Darcy's law (equation
2.10).  For Darcy flow, the properties of the (rock) medium must be identified
throughout the domain, including the phase component permeabilities,
etc.  To complete the mathematical specification, any sources/sinks
for the evolved field quantities must be communicated as well.

A number of primitives are supported by Amanzi to aid in communicating this
data to the executable:

 * Domain: This is the definition of the computational domain. For structured AMR simulations, the domain is defined by opposite corners of the domain. For unstructured mesh based simulations, the domain is defined as a union of regions (defined below).

 * Region: Generically, this is a subset of the computational domain. The construct can be used also to communicate a particular boundary or computational volume for a variety of purposes, including boundary conditions, line-integrals or volumetric diagnostics. Regions are defined as a box specified by opposite corners or by a combination of surfaces (defined) below bounding the volume in a watertight manner. Although a computational domain is defined by a list of regions, some regions may be defined for other purposes such as viz or uncertainty quantification. Regions may be assigned names and/or numbers.

 * Surface: Surfaces may be defined analytically or as a triangular mesh embedded in 3D. Not all surfaces are used in the definition of regions; some may be defined for the purpose of visualization or uncertainty quantification. However, urfaces that are used to define a mesh region have additional constraints - they must not extend beyond the extents of the region or be overlapping with each other. In other words, the union of the surfaces defining a region must strictly be the boundary of the region and no post-processing of the surfaces must be required. If a surface is defined as a triangular mesh, it's definition points to a mesh file containing the mesh. Surfaces may be assigned names and/or numbers. Boundary conditions are then assigned to named or numbered surfaces.

 * Rock: This contains a characterization of the relevant properties of the flow substrate.  Because of the massive range of length scales involved in typical groundwater flow models, the domain may be comprised of several subregions with sharply distinct rock properties.

 * State: The names of each of the state quantities and their ordering in the state, along with a consistent set of initial data and boundary conditions.

 * Mesh: The mesh is defined only in unstructured mesh based simulations and the definition in the input file is merely a link to an input file. Eventually, Amanzi will require that mesh elements be linked to the regions (see above) they reside in and mesh faces be linked to the surfaces (see above) they lie on. This will enable a trivial application of boundary conditions and material properties to any mesh.


Given these primitives, a user can complete the problem specification by configuring the source terms for each component (if applicable), and instructions for generating the observation data array to output.  The sections below give a catalogue of parameters that are expected for each of the major sections of communication between Amanzi and the user.  In particular, instructions are provided for generating named sug-regions, which then are used to define other parts of the specification.

Notes:

 * Currently, the data is passed from the user into the Amanzi executable via an XML-capable parameter list object.  Thus, the input data set is described here from the point of view of a user that is constructing such an XML file for this purpose.

 * This specification incorporates functionality that the HPC developers believe to be sufficient to set up and execute a broad range of simple model problems -- the Platform group has yet to comment on the suitablity or completeness of this specification. Also, as of the writing of this document, not all functionality discussed in this specification is completely implemented in Amanzi. As a final stage in the production of this document, features that are not fully implemented will be clearly identified.

 * The Amanzi code has a dual execution path, catering to the special requirements, and exploiting many of the unique advantages of structured versus unstructured mesh implementations, respectively.  Completeness of the implementation of the models discussed here will vary between mesh schemes as well, and will evolve over time.

 * On occasion below, we refer to an array type as if it were an XML intrinsic.  Such a construct does not yet exist.  We expect that an `"array double`" value might be something like `"2.3 4.0 6`".  Alternatively, such an array can be easily written as an XML list.



Creating a Parameter List
=======================================

The XML input file must be framed at the beginning and end by the following statements:


.. code-block:: xml

  <ParameterList name="Main">

  </ParameterList>

where the name can be anything. Here I have named the main list "Main".

We have several sublists that are read and interpreted in various parts of Amanzi.  The sublist names head each of the
following sections

1. State
=======================================

The state specifies the distribution of phases, species and pressure in the system.  Generally, there
can be multiple phases (e.g. gaseous, aqueous, etc), and each is comprised of a number
of components.  Additionally, each phase may carry a number of trace species.  The tracers are assumed to
have no impact on the thermodynamic and transport properties of the phases, but may be involved in chemical
processes.  In this section, each state component is labeled and defined, in terms of physical properties.
Tracers are defined, and tracer groups are constructed in order to support requirements of the chemistry solver.
There is an implied consistency between the state identifiers defined here and 
those referenced elsewhere, such as in the chemistry specification.  To the extent possible, Amanzi will
recognize when this expected consistency is violated.

In this section, "region" labels are assumed to have been defined elsewhere.  These labels are used to
identify portions of the interior for setting initial state data, and portions of the boundary for setting
boundary conditions.  We summarize the acceptable data for "state" in a hierarchical list.
(here, and in the remainder of this document, reserved keywords and labels are `"quoted`"; representative
user-defined labels are indicated with ALL-CAPS).

* "state" (list) can accept lists for named components (COMP), `"add tracer`" (list) to add a tracer, or boundary conditions (BC), a value for the name of the dominant component, and a group (string array) identifying a group of tracers

  * COMP (list) can accept values for phase name (string), mass density (double), viscosity (double) and diffusivity (double), and a list for initial data (IC)

    * IC (list) can accept a list of region names (REGION) 

      * REGION (list) can accept a list (IC-PARAM) to specify the parameters of a named functional

        * IC-PARAM (list) can accept a set of parameter values for the functional

    * `"mass density`" (double) the density of this component

    * `"viscosity`" (double) the viscosity of this component

    * `"diffusivity`" (double) the diffusivity of this component

    * `"phase name`" (string) the phase that this component is in

  * `"add tracer`" (list) can accept a list for initial data (IC), name (string), and name of parent component (string)

    * `"name`" (string) name of tracer

    * IC (list) can accept a list of region names (REGION) 

      * REGION (list) can accept a list (IC-PARAM) to specify the parameters of a named functional

        * IC-PARAM (list) can accept a set of parameter values for the functional

    * `"parent phase component`" (string) name of carrying phase component (must be one of COMP defined in this file)

  * BC (list) can accept a list (BC-REG) named after the "region" that defines a surface that bounds the computational domain
 
    * BC-REG (list) can accept a list (BC-PARAM) to specify the parameters of a named functional

      * BC-PARAM (list) can accept a set of parameter values for the functional

  * `"dominant component`" (string) must be the name of one of the COMP lists defined above

  * `"add group`" (string array) an array of strings of names, taken from the `"add tracer`" list defined above

Initial conditions are required for each phase component over the entire computational domain.
Boundary conditions are required on all domain boundaries.  Both are constructed using a limited number
of explicitly parameterized functional forms.  If the simulation is to be intialized using a restart file,
the phase and component definitions are taken from the restart file, and initial condition instructions specified 
here are quietly ignored.  Boundary conditions are required regardless of the initial data, and must be defined
consistently.

The following parameterized distrubution functional are supported:
 * `"ic: constant`" requires `"value`" (see note below)
 * `"ic: coordinate-aligned linear`"  requires direction `"dir`" (string) of variation, `"x0_y0_slope`" (array double) specifying x0, y0, m for function of the form: `y-y0 = m*(x-x0)`.  Here y is state value, x is coordinate in `"dir`" direction.  For state values however, see note below.
 * `"ic: quadratic`" similar to linear
 * `"ic: exponential`" similar to linear

The following parameterized boundary conditions are supported:
 * `"bc: inflow`" requires `"bc: distribution`" (string) to set the distribution of the state upstream of the boundary (outside domain)
 * `"bc: outflow`"  requires `"bc: distribution`" (string) to set the distribution of the state downstream of the boundary (outside domain)
 * `"bc: seepage`" requires location `"water table height`" (double) of the water table.  If a more complex specification is needed, this should be changed to require list to define it appropriately.
 * `"bc:  noflow`" requires no parameter data

NOTES:

In initial data and boundary conditions, the user must be able to specify the desired data to implement for each state
component.  It is unclear how to create a manageable interface for this.  Additionally, there may be a desire
to initialize via file when the file is not a restart file.  To add such a method, one requires an implicit
interpolation operator to map the stored solution to the current one, and that the data files are sufficiently
self-describing and contain enough information to support this interpolation.  Amanzi does not support such an
option at this time.

Due to various physical constraints (e.g. component saturations sum to unity), initial and boundary functionals
not explicitly specified will be derived, if possible.  If insufficient or contradicting information is detected,
an error will be thrown.


Example:

.. code-block:: xml

  <ParameterList name="state">
    <Parameter name="dominant component" type="string" value="air"/>    
    <ParameterList name="air">
      <Parameter name="phase" type="string" value="gaseous"/>
      <Parameter name="mass density" type="double" value="1.2"/>
      <Parameter name="viscosity" type="double" value="0.018"/>
      <Parameter name="diffusivity" type="double" value="0."/>
      <ParameterList name="top">
        <ParameterList name="ic: constant">
          <Parameter name="value" type="double" value=".5"/>
        </ParameterList>   
      </ParameterList>   
      <ParameterList name="middle">
        <ParameterList name="ic: constant">
          <Parameter name="value" type="double" value=".4"/>
        </ParameterList>   
      </ParameterList>   
      <ParameterList name="bottom">
        <ParameterList name="ic: coordinate-aligned linear"/>
          <Parameter name="direction" type="string" value="x"/>
          <Parameter name="x0_y0_slope" type="array double" value=".4 .9 3"/>
        </ParameterList>   
      </ParameterList>   
    </ParameterList> 
    <ParameterList name="water">
      <Parameter name="phase" type="string" value="aqueous"/>
      <Parameter name="density" type="double" value="1.e3"/>
      <Parameter name="viscosity" type="double" value="1.0"/>
      <Parameter name="diffusivity" type="double" value="0."/>
    </ParameterList>   
    <ParameterList name="boundary conditions">
      <ParameterList name="XLOBC">
        <ParameterList name="inflow">
          <ParameterList name="bc: constant">
          </ParameterList> 
        </ParameterList> 
      </ParameterList> 
      <ParameterList name="XHIBC">
        <ParameterList name="outflow">
        </ParameterList> 
      </ParameterList> 
    </ParameterList>
    <ParameterList name="add tracer">
      <Parameter name="name" type="string" value="Uranium"/>
      <Parameter name="parent phase component" type="string" value="water"/>
      <ParameterList name="all">
        <ParameterList name="ic: constant">
          <Parameter name="value" type="double" value=".004"/>
        </ParameterList> 
      </ParameterList> 
    </ParameterList>

In this example, there are 2 phases (water, air).  Each phase consists of a single component.  Three
volumetric regions ("top", "middle" and "bottom"), and two boundary regions (XLOBC and XHIBC)
have been defined elsewhere.  The initial data for the fields are set using a combination of linear and
constant profile functions over the two volumetric regions.  The boundary conditions are Dirichlet inflow
on the low side and outflow on the high side.


2. Regions
=======================================

Regions are used in Amanzi to define the physical extent of the simulation domain and its bounding surfaces.
Regions are also used to specify initial data and boundary conditions, and to define output data expected
upon return from the simulator.  Currently, it is assumed that the simulation domain (region = "all") is a
parallelepiped and must be defined explicitly. 6 additional regions are automatically constructed for every
simulation: XLOBC, XHIBC, YLOBC, YHIBC, ZLOBC and ZHIBC, represent each rectangular side of "all" (aligned
with the coordinate axis).


* "regions" (list) can accept lists for named regions (REGION)

  * REGION (list) can accept lists (SHAPE) that specify a functional for its shape.  Though Amanzi currently supports only a single shape specifier per region, this limitation may be removed in the future.

    * SHAPE (list) can accept lists of shape parameters (SHAPE-PARAMS) 

      * SHAPE-PARAMS (double array or string) parameters to specify shape

Currently, Amanzi supports parameterized forms for a number of analytic shapes, as well as more complex
definitions based on triangulated surface files: point, box, arbitrary, layer.  Depending on the functional, SHAPE requires
a number of parameters:

+------------------------+-------------------------+------------------------------+---------------------------------------------------------------------------------------------+
|  shape functional name | parameters              | type(s)                      | Comment                                                                                     |
+========================+=========================+==============================+=============================================================================================+
| `"point"`              | `"loc`"                 | double array                 | Location of point in space                                                                  |
+------------------------+-------------------------+------------------------------+---------------------------------------------------------------------------------------------+
| `"box"`                | `"lo`", `"hi`"          | double array, double array   | Location of boundary points of box                                                          |
+------------------------+-------------------------+------------------------------+---------------------------------------------------------------------------------------------+
| `"arbitrary"`          | `"file`"                | string                       | Region enveloped by surface described in specified file (see note below for format of file) |
+------------------------+-------------------------+------------------------------+---------------------------------------------------------------------------------------------+
| `"layer"`            | `"file_lo`" `"file_hi`" | string, string               | Region between surfaces described in specified files (see note below for format of file)    |
+------------------------+-------------------------+------------------------------+---------------------------------------------------------------------------------------------+
| `"surface"`            | `"id1`" `"name2`" ... `"idN`" | string, string ,..., string               | Region between surfaces described in specified files (see note below for format of file)    |
+------------------------+-------------------------+------------------------------+---------------------------------------------------------------------------------------------+

Note: surface file format TBD.


Example:

.. code-block:: xml

  <ParameterList name="regions">
    <ParameterList name="all">
      <ParameterList name="box">
        <Parameter name="lo" type="double array" value="2 3 4"/>
        <Parameter name="hi" type="double array" value="4 5 8"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="top">
      <ParameterList name="box">
        <Parameter name="lo" type="double array" value="2 3 6"/>
        <Parameter name="hi" type="double array" value="4 5 8"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="middle">
      <ParameterList name="box">
        <Parameter name="lo" type="double array" value="2 3 6"/>
        <Parameter name="hi" type="double array" value="4 5 8"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="bottom">
      <ParameterList name="box">
        <Parameter name="lo" type="double array" value="2 3 4"/>
        <Parameter name="lo" type="double array" value="4 5 6"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>

In this example, a simulation domain is defined to be 2x2x4 with its lower bound at the point (2,3,4).  Three box-shaped sub-regions are defined for an unspecified
purpose.


4. Rock
=======================================

Rock properties must be specified over the entire simulation domain ("all") defined in the Region section.  This can be implemented using any combination of regions
defined above, provided that the entire domain is covered.  Currently, the regions used should be disjoint.  Amanzi may eventually support verifying this condition,
and/or specifying a precedence order for overalapping regions.

Each rock type is given a label (string) and assigned a density (double) and models (string) for porosity, permeability and capillary pressure.  Each rock is assigned to
regions (string array), a list of regions.

* "rock" (list) can accept multiple lists for named rock types (ROCK)

  * ROCK (list) can accept lists to specify a model (MODEL) for porosity, relative permeability and capillary pressure, and values for the `"density`" (double) and `"permeability`" (double array) - values in the three principal axes (currently assumed to align with the coordinate axes and grid).  It can also accept a string array `"regions`" to specify where these properties apply.

    * MODEL (list) can accept model parameters (MODEL-PARAMS) 

      * MODEL-PARAMS (double, double array) parameters to specify model (see notes below for details of each model available)

    * `"regions`" (string array) a set of regions defined above

The following models are currently supported for porosity:
 * `"porosity: file`" requires a single string "filename" specifying the name of a file.  This file must be written in a self-describing format that is consistent with that of the current meshing option (sturctured_grid or unstructured_grid).  In particular, the physical domain of the input data must completely cover the current "all" region, and the data must exist on discrete cells that are consistent with the current meshing configuration.  This option is not currently supported under the unstructured option.
 * `"porosity: uniform`" requires the porosity (double)
 * `"porosity: random`" requires the mean value of porosity and the percentage fluctuation, "porosity and fluctuation" (double array) to generate
 * `"porosity: gslib`" requires the name of a gslib-formatted file "gslib filename" to generate porosity (plus other data?)

The following models are currently supported for relative permeability:
 * `"perm: perfect`" requires no parameters, krl=krg=1
 * `"perm: linear`" requires no parameters, krl=sl and krg=sg
 * `"perm: quadratic`" requires slr, sgr (double array), krl=sc^2, krg=1-se^2, se=(sl-sg)/(1-slr-sgr)
 * `"perm: vGM`" (van Genuchten-Mualem) requires m, slr, sgr (double array), krl=sqrt(se)(1-(1-se^-m)^m)^2, krg=(1-sekg)^1/3 (1-sekg^-m)^(2m), se=(sl-slr)/(1-slr-sgr), sekg=sl/(1/sgr)

The following models are currently supported for capillary pressure:
 * `"pc: none`" requires no parameters, pc = 0
 * `"pc: linear`" requires no parameters, pc = sl
 * `"pc: vG`" requires m, sigma, slr, sgr (double array), pc=(1/sigma)(se^-m - 1)^-n, se=(sl-slr)/(1-slr-sgr)

Example:

.. code-block:: xml

  <ParameterList name="rock">
    <ParameterList name="backfill">
      <Parameter name="density" type="double" value="2.8e3"/>
      <Parameter name="permeability" type="double array" value="1240 1240 1240"/>
      <ParameterList name="porosity: uniform">
        <Parameter name="porosity" type="double" value="0.2585"/>
      </ParameterList>
      <ParameterList name="perm: vGM">
        <Parameter name="m_slr_sgr" type="double array" value="0.6585 0.0774 0"/>
      </ParameterList>
      <ParameterList name="pc: vG">
        <Parameter name="m_sigma_slr_sgr" type="double array" value="0.6585 102.1 0.0774 0"/>
      </ParameterList>
      <Parameter name="regions" type="string array" value="top bottom"/>
    </ParameterList>
    <ParameterList name="fine sand">
      <Parameter name="density" type="double" value="2.8e3"/>
      <Parameter name="permeability" type="double array" value="337.0 337.0 337.0"/>
      <ParameterList name="porosity: uniform">
        <Parameter name="porosity" type="double" value="0.3586"/>
      </ParameterList>
      <ParameterList name="perm: vGM">
        <Parameter name="m_slr_sgr" type="double array" value="0.4694 0.0837 0"/>
      </ParameterList>
      <ParameterList name="pc: vG">
        <Parameter name="m_sigma_slr_sgr" type="double array" value="0.4694 9.533 0.0837 0"/>
      </ParameterList>
      <Parameter name="regions" type="string array" value="middle"/>
    </ParameterList>
  </ParameterList>

In this example, there are two types of rock, `"backfill`" (which fills bottom and top regions) and `"fine sand`" (which fills middle region).  Both have
van Genuchten models for relative permeability and capillary pressure.


5. Source
=======================================

Specification of volumetric source terms for the various state quantities depends on the definition of regions and state labels.

Each source is given a label (string), state id (string), integrated source strength (double), distribution functional (list) and region (string).

* "source" (list) can accept multiple lists for named sources (SOURCE)

  * SOURCE (list) can accept a list to specify a distribution (DIST), and values for `"state id`", `"region`" and `"strength`".

    * DIST (list) can accept shape parameters (DIST-PARAMS) 

      * DIST-PARAMS (double, double array) parameters to specify shape model (see notes below for details of each model available)

    * `"region`" (string) a region defined above

    * `"state id`" (string) a state quantity defined above

    * `"strength`" (double) integrated source strength

The following models are currently supported for source distribution:
 * `"source: uniform`" requires no parameters
 * `"source: linear`" requires location `"loc`" (double array) of a point or two locations, `"lo`", `"hi`" specifying a line or a rectangular plane
 * `"source: quadratic`" requires location `"loc`" (double array) of a point or two locations, `"lo`", `"hi`" specifying a line or a rectangular plane
 * `"source: exponential`" requires the exponent, `"exp`" and location `"loc`" (double array) of a point or two locations, `"lo`", `"hi`" specifying a line or a rectangular plane


Example:

.. code-block:: xml

  <ParameterList name="source">
    <ParameterList name="infiltration">
      <Parameter name="state id" type="string" value="water"/>
      <Parameter name="region" type="string" value="top"/>
      <Parameter name="strength" type="double" value="7.6e-6"/>
      <ParameterList name="source: uniform">
      </ParameterList>
    </ParameterList>
    <ParameterList name="tracer discharge">
      <Parameter name="state id" type="string" value="all tracers"/>
      <Parameter name="region" type="string" value="bottom"/>
      <Parameter name="strength" type="double" value="3.6e-7"/>
      <ParameterList name="source: uniform">
      </ParameterList>
    </ParameterList>
  </ParameterList>

In this example, there is an infiltration source of water in the top region, and a discharge of all the tracers through the bottom.

6. Observation Data
=======================================

Observation data generally refers to a small list of diagnostic quantities to extract from a simulation.  This is to be contrasted with the
detailed field data typically associated with analysis and restart data.  Examples of these quantities include various volume and surface integrals
used to characterize the response of the system to variations of input data.  Computation of observation data involves applying a parameterized 
functional on a specified state quantity or flux value at specific values of the simulation time.

Each observation is given a label (string), state id (string), evaluation functional (list), region (string) and a list of times for evaluation.
The results are passed to the calling process in a structure that is consistent with the following specification.

* "observation" (list) can accept multiple lists for named observations (OBSERVATION)

  * OBSERVATION (list) can accept values for `"state id`", `"region`", `"functional`" and `"times`"

    * `"region`" (string) a region defined above

    * `"state id`" (string) a state quantity defined above

    * `"functional`" (string) choses which funcitional to apply (see below)

    * `"times`" (double array) values of time where this quantity is desired

The following observation functionals are currently supported
 * `"observation: average`" 
 * `"observation: integral`" 
 * `"observation: squared integral`" 
 * `"observation: peak value`" 

Example:

.. code-block:: xml

  <ParameterList name="observation">
    <ParameterList name="mass of water">
      <Parameter name="state id" type="string" value="water"/>
      <Parameter name="region" type="string" value="all"/>
      <Parameter name="functional" type="string" value="integral"/>
      <Parameter name="times" type="double array" value="1.e3 2.e3 2.5e3"/>
    </ParameterList>
  </ParameterList>

In this example, the user requests the volume integral of the water density over the entire domain at three different times.
Amanzi will also support integrals and point samples of phase fluxes.  Note that times specified may not necessarily fall within
the time interval of the present simuluation.  The format of the data structure used to convey the observation data includes
a flag for each time to indicate whether the quantity was successfully filled.

7. Chemistry
=======================================

This section is completely unintelligible, and needs to be re-written.  In the structured_grid implementation, the following are the only chemistry-related 
inputs currently allowed:

+---------------------+--------+----------------------------------------------------------------------+
| Name                | Type   | Description                                                          |
+=====================+========+======================================================================+
| `"do chemistry`"    | int    | If 0, disable chemistry                                              |
+---------------------+--------+----------------------------------------------------------------------+
| `"chemistry file`"  | string | Amanzi-formatted chemistry input file                                |
+---------------------+--------+----------------------------------------------------------------------+
| `"interval`"        | int    | Number of coarse-grid time steps between chemistry solver invocation |
+---------------------+--------+----------------------------------------------------------------------+
| `"splitting order`" | int    | Accuracy order of chemistry evolution (1, 2)                         |
+---------------------+--------+----------------------------------------------------------------------+

....original text...

Example:

.. code-block:: xml

  <ParameterList name="Chemistry">
    <Parameter name="Thermodynamic Database Format" type="string" value="simple" />
    <Parameter name="Thermodynamic Database File" type="string" value="fbasin-uo2-5-component.bgd" />
    <Parameter name="Verbosity" type="int" value="1" />
    <Parameter name="Activity Model" type="string" value="debye-huckel" />
    <Parameter name="Tolerance" type="double" value="1.5e-12"/>
    <Parameter name="Max Time Step (s)" type="double" value="86400.0"/>
    <Parameter name="Maximum Newton Iterations" type="int" value="150"/>
    <Parameter name="Using sorption" type="string" value="yes"/>
    <Parameter name="Free ion concentrations provided" type="string" value="yes"/>
    <ParameterList name="Initial Conditions">
      <Parameter name="Number of minerals" type="int" value="3"/>
      <Parameter name="Number of ion exchange sites" type="int" value="0"/>
      <Parameter name="Number of sorption sites" type="int" value="0"/>
      <Parameter name="Number of mesh blocks" type="int" value="1"/>
      <ParameterList name="Mesh block 1"> 
        <Parameter name="Mesh block ID" type="int" value="0"/>
        <ParameterList name="Free Ion Species">
	  <Parameter name="Free Ion Species 0" type="double" value="4.36476e-16"/>  <!-- Al+++ -->
	  <Parameter name="Free Ion Species 1" type="double" value="3.16772e-08"/>  <!-- H+ -->
	  <Parameter name="Free Ion Species 2" type="double" value="1.00000e-06"/>  <!-- HPO4-2 -->
	  <Parameter name="Free Ion Species 3" type="double" value="1.87000e-04"/>  <!-- SiO2(aq) -->
	  <Parameter name="Free Ion Species 4" type="double" value="1.84374e-20"/>  <!-- UO2++ -->
        </ParameterList>
        <ParameterList name="Minerals">
          <Parameter name="Mineral 0" type="double" value="0.15"/>  <!-- Kaolinite -->
          <Parameter name="Mineral 1" type="double" value="0.21"/>  <!-- Quartz -->
          <Parameter name="Mineral 2" type="double" value="0.0"/>   <!-- (UO2)3(PO4)2.4H2O -->
        </ParameterList>
      </ParameterList>
    </ParameterList>
  </ParameterList>


'''Note: all chemistry names and values are case sensitive.'''

+------------------------------------+---------------+------------------+-----------------------------+----------------------------------------------------------------------------------------+
|  Parameter name                    |  Type         |  Default Value   | Optional Values             | Purpose                                                                                |
+====================================+===============+==================+=============================+========================================================================================+
| `"Thermodynamic Database Format"`  | string        | `"simple`"       | `"simple"`                  | set the format of the database                                                         |
+------------------------------------+---------------+------------------+-----------------------------+----------------------------------------------------------------------------------------+
| `"Thermodynamic Database File"`    | string        | `"dummy.dbs"`    |  ---                        | path name to the chemistry database file, relative to the program execution directory. |
+------------------------------------+---------------+------------------+-----------------------------+----------------------------------------------------------------------------------------+

The following parameters are optional in the Chemistry parameter list:

+---------------------------------------+---------------+------------------+-----------------------------+-------------------------------------------------------------------------------------+
|  Parameter name                       | Type          | Default Value    | Optional Values             | Purpose                                                                             |
+=======================================+===============+==================+=============================+=====================================================================================+
| `"Verbosity"`                         | int           | 0                | 0, 1, 2, 3, 4, 5, 6, ...    | set the verbosity level of chemistry: 0=silent, 1=terse warnings, 2=verbose details,|
|                                       |               |                  |                             |  3=debug, 4=debug beaker, 5=debug database file, ....                               | 
+---------------------------------------+---------------+------------------+-----------------------------+-------------------------------------------------------------------------------------+
| `"Activity Model"`                    | string        | `"unit`"         | `"unit"`, `"debye-huckel"`  | set the model used for activity corrections                                         |
+---------------------------------------+---------------+------------------+-----------------------------+-------------------------------------------------------------------------------------+
| `"Tolerance"`                         | double        | 1.0e-12          |  ---                        | set the tolerance for newton iterations within chemistry                            |
+---------------------------------------+---------------+------------------+-----------------------------+-------------------------------------------------------------------------------------+
| `"Maximum Newton Iterations"`         | int           | 200              | ---                         | set the maximum number of newton iterations for chemistry.                          |
+---------------------------------------+---------------+------------------+-----------------------------+-------------------------------------------------------------------------------------+
| `"Max Time Step (s)"`                 | double        | 9.9e9            | ---                         | set the maximum time step allowed for chemistry.                                    |
+---------------------------------------+---------------+------------------+-----------------------------+-------------------------------------------------------------------------------------+
| `"Using sorption"`                    | string        | `"no"`           | `"yes"`                     | Tells the chemistry module whether to allocate memory for sorption.                 |
+---------------------------------------+---------------+------------------+-----------------------------+-------------------------------------------------------------------------------------+
| `"Free ion concentrations provided"`  | string        | `"no"`           | `"yes"`                     | Tells chemistry that in initial guess for free ion concentrations is provided in    |
|                                       |               |                  |                             | the xml file.                                                                       |
+---------------------------------------+---------------+------------------+-----------------------------+-------------------------------------------------------------------------------------+

The initial condition list must have a `"Mesh Block"` parameter list for each mesh block, mesh block numbering should correspond to the other mesh block lists. Each mesh block list will have parameter lists for the non-zero elements of the chemistry. Valid parameter list names are: `"Free Ion Species"` `"Minerals"` `"Ion Exchange Sites"` `"Sorption Sites"`.

Each initial condition list should contain a parameter name constructed like `"Type #"` where `"Type"` is `"Mineral"`, `"Free Ion Species"`, `"Ion Exchange Site"` `"Sorption Site"` and `"#"` in the integer identifier, starting with zero.  

'''Note: it is recommended that you include an xml comment with the species or mineral name after each initial condition. The xml parser expects every instance of `"--"` to mark a comment, so species names with negative charges should be written as `"SO4-2"` rather than `"SO4--"`.'''

7.1 Chemistry Thermodynamic data file formats 
-------------------------------------------------

7.1.1 XML Database format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

not yet implemented

7.1.2 Simple Database format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Importing thermodynamic data into the chemistry module using the `"simple"` (file extension `"bgd"`) format requires the user to explicitly specify all the species and reactions for the problem. There is no basis switching or automatic species and reaction selection. Below is an example of a `"simple"` database file for a five component uranium problem with mineral dissolution and surface complexation:

::

 <Primary Species
 # name               ; debye-huckel a0 ; charge ; GMW     

 Al+++                ;   9.0 ;   3.0 ;  26.9815
 H+                   ;   9.0 ;   1.0 ;   1.0079
 HPO4--               ;   4.0 ;  -2.0 ;  95.9793
 SiO2(aq)             ;   3.0 ;   0.0 ;  60.0843
 UO2++                ;   4.5 ;   2.0 ;  270.028

 <Aqueous Equilibrium Complexes
 # name               =  coeff primary_name  coeff primary_name  ; log10(Keq) 25C ; debye-huckel a0 ; charge ; GMW      

 OH-                  =  1.0 H2O  -1.0 H+                ;    13.9951 ;   3.5 ;  -1.0 ;  17.0073 
 AlOH++               =  1.0 H2O  1.0 Al+++  -1.0 H+     ;     4.9571 ;   4.5 ;   2.0 ;  43.9889 
 Al(OH)2+             =  2.0 H2O  1.0 Al+++  -2.0 H+     ;    10.5945 ;   4.0 ;   1.0 ;  60.9962 
 Al(OH)3(aq)          =  3.0 H2O  1.0 Al+++  -3.0 H+     ;    16.1577 ;   3.0 ;   0.0 ;  78.0034 
 Al(OH)4-             =  4.0 H2O  1.0 Al+++  -4.0 H+     ;    22.8833 ;   4.0 ;  -1.0 ;  95.0107 
 UO2OH+               =  1.0 H2O  -1.0 H+  1.0 UO2++     ;     5.2073 ;   4.0 ;   1.0 ;  287.035 
 UO2(OH)2(aq)         =  2.0 H2O  -2.0 H+  1.0 UO2++     ;    10.3146 ;   3.0 ;   0.0 ;  304.042 
 UO2(OH)3-            =  3.0 H2O  -3.0 H+  1.0 UO2++     ;    19.2218 ;   4.0 ;  -1.0 ;   321.05 
 UO2(OH)4--           =  4.0 H2O  -4.0 H+  1.0 UO2++     ;    33.0291 ;   4.0 ;  -2.0 ;  338.057 
 (UO2)2OH+++          =  1.0 H2O  -1.0 H+  2.0 UO2++     ;     2.7072 ;   5.0 ;   3.0 ;  557.063 
 (UO2)2(OH)2++        =  2.0 H2O  -2.0 H+  2.0 UO2++     ;     5.6346 ;   4.5 ;   2.0 ;   574.07 
 (UO2)3(OH)4++        =  4.0 H2O  -4.0 H+  3.0 UO2++     ;     11.929 ;   4.5 ;   2.0 ;  878.112 
 (UO2)3(OH)5+         =  5.0 H2O  -5.0 H+  3.0 UO2++     ;    15.5862 ;   4.0 ;   1.0 ;   895.12 
 (UO2)3(OH)7-         =  7.0 H2O  -7.0 H+  3.0 UO2++     ;    31.0508 ;   4.0 ;  -1.0 ;  929.135 
 (UO2)4(OH)7+         =  7.0 H2O  -7.0 H+  4.0 UO2++     ;    21.9508 ;   4.0 ;   1.0 ;  1199.16 
 UO2(H2PO4)(H3PO4)+   =  3.0 H+  2.0 HPO4--  1.0 UO2++   ;   -22.7537 ;   4.0 ;   1.0 ;   465.01 
 UO2(H2PO4)2(aq)      =  2.0 H+  2.0 HPO4--  1.0 UO2++   ;   -21.7437 ;   3.0 ;   0.0 ;  464.002 
 UO2HPO4(aq)          =  1.0 HPO4--  1.0 UO2++           ;    -8.4398 ;   3.0 ;   0.0 ;  366.007 
 UO2H2PO4+            =  1.0 H+  1.0 HPO4--  1.0 UO2++   ;   -11.6719 ;   4.0 ;   1.0 ;  367.015 
 UO2H3PO4++           =  2.0 H+  1.0 HPO4--  1.0 UO2++   ;   -11.3119 ;   4.5 ;   2.0 ;  368.023 
 UO2PO4-              =  -1.0 H+  1.0 HPO4--  1.0 UO2++  ;    -2.0798 ;   4.0 ;  -1.0 ;  364.999 

 <Minerals
 # name               =  coeff primary_name  coeff primary_name  ; log10(Keq) 25C ; GMW      ; molar volume [cm^2/mol] ; SSA [m^2/g] 

 Kaolinite            =  5.00 H2O  2.00 Al+++  -6.00 H+  2.00 SiO2(aq)  ;     6.8101 ;   258.16 ;    99.52 ;   1.0 
 Quartz               =  1.00 SiO2(aq)  ;    -3.9993 ;  60.0843 ;   22.688 ;   1.0 
 (UO2)3(PO4)2.4H2O    =  4.00 H2O  -2.00 H+  2.00 HPO4--  3.00 UO2++  ;   -27.0349 ;  1072.09 ;    500.0 ;   1.0 

 <Mineral Kinetics
 # name               ; TST ; log10_rate_constant double     moles_m2_sec 

 Kaolinite            ; TST ; log10_rate_constant    -16.699 moles_m2_sec 
 Quartz               ; TST ; log10_rate_constant      -18.0 moles_m2_sec 
 (UO2)3(PO4)2.4H2O    ; TST ; log10_rate_constant      -10.0 moles_m2_sec 

 <Surface Complex Sites
 # name               ; surface_density

 >FeOH                ; 6.3600E-03
 >AlOH                ; 6.3600E-03
 >SiOH                ; 6.3600E-03

 <Surface Complexes
 # name               =  coeff surface site  coeff primary_name  ; log10(Keq) 25C ; charge 

 >SiOUO3H3++          =  1.0 >SiOH  1.0 H2O  1.0 UO2++  ;       5.18 ;   2.0 
 >SiOUO3H2+           =  1.0 >SiOH  1.0 H2O  -1.0 H+  1.0 UO2++  ;       5.18 ;   1.0 
 >SiOUO3H             =  1.0 >SiOH  1.0 H2O  -2.0 H+  1.0 UO2++  ;       5.18 ;   0.0 
 >SiOUO3-             =  1.0 >SiOH  1.0 H2O  -3.0 H+  1.0 UO2++  ;      12.35 ;  -1.0 
 >SiOUO2(OH)2-        =  1.0 >SiOH  2.0 H2O  -3.0 H+  1.0 UO2++  ;      12.35 ;  -1.0 
 >FeOHUO3             =  1.0 >FeOH  1.0 H2O  -2.0 H+  1.0 UO2++  ;       3.05 ;   0.0 
 >FeOHUO2++           =  1.0 >FeOH  1.0 UO2++  ;      -6.63 ;   2.0 
 >AlOUO2+             =  1.0 >AlOH  -1.0 H+  1.0 UO2++  ;      -3.13 ;   1.0 

Note the following about this format:

 * Any line starting with a `"#"` or space character is a comment. 

 * Data is separated into sections, where each section of the file is starts with a line containing `"<Section Name"`. The valid section names are: `"Primary Species"`, `"Aqueous Equilibrium Complexes"`, `"Minerals"`, `"Mineral Kinetics"`, `"General Kinetics"`, `"Ion Exchange Sites"`, `"Ion Exchange Complexes"`, `"Surface Complex Sites"`, `"Surface Complexes"`. The less than character, `"<"`, should be the first character on the line and there is no space between the character and the section name.

 * Sections should be ordered so that the primary species, minerals, and exchange sites come before any reactions using those species.

 * Each line within a section is a semi-colon delimited

 * A primary species line contains the primary species must contain:

   ::

     # name               ; debye-huckel a0 ; charge ; GMW [grams/mole]    
     Al+++                ;   9.0 ;   3.0 ;  26.9815

 * An aqueous equilibrium complex line contains a reaction and data for the reaction on a single line:

   ::

     # name               =  coeff primary_name  coeff primary_name ... ; log10(Keq) 25C ; debye-huckel a0 ; charge ; GMW [grams/mole]     
     OH-                  =  1.0 H2O  -1.0 H+  ;    13.9951 ;   3.5 ;  -1.0 ;  17.0073 

   The reaction is written as product species = reactants.... The coefficient of the product aqueous complex is assumed to be 1.0, and the reactants must be primary species. The equilibrium constant is for a fixed temperature of 25C.

 * Minerals and other complexes follow the same convention as aqueous equilibrium complexes, with additional data as needed.

   ::

     <Minerals
     # name               =  coeff primary_name  coeff primary_name ... ; log10(Keq) 25C ; GMW      ; molar volume [cm^2/mol] ; SSA [m^2/g] 

     <Surface Complexes
     # name               =  coeff surface site  coeff primary_name ... ; log10(Keq) 25C ; charge 

     These are all minerals present in the system during the simulation, including those that may precipitate later. They are used for calculating saturation states, but not equilibrium or kinetic calculations.

 * The mineral kinetics section lists the name of a mineral found in the mineral section, the type of rate law, and rate parameters for that law.

   :: 

     <Mineral Kinetics
     # name               ; TST ; log10_rate_constant double     moles_m2_sec ; primary_name coeff ....
 
   Currently only the `"TST"` rate law is implemented. The keywords "log10_rate_constant" and "moles_m2_sec" must be present in the line, but no unit conversion are currently preformed. The    modifying primary species terms follow the rate constant, along with their exponent coefficients.

 * Surface complex sites are listed by name and surface density:

   ::

     <Surface Complex Sites 
     # name               ; surface_density [moles sites / m^2 mineral]


8. MPC
=======================================

Example:

.. code-block:: xml

  <ParameterList name="MPC">
    <Parameter name="Start Time" type="double" value="0.0"/>
    <Parameter name="End Time" type="double" value="0.1"/>
    <Parameter name="End Cycle" type="int" value="10000"/>
    <Parameter name="Flow model" type="string" value="Darcy"/>
    <Parameter name="disable Flow_PK" type="string" value="no"/>
    <Parameter name="disable Transport_PK" type="string" value="no"/>
    <Parameter name="disable Chemistry_PK" type="string" value="yes"/>
    <Parameter name="Viz dump cycle frequency" type="int" value="10"/>
    <Parameter name="Viz dump time frequency" type="double" value="0.05"/>
    <ParameterList name="CGNS">
      <Parameter name="File name" type="string" value="test1.cgns"/>
    </ParameterList>
  </ParameterList> 

In the MPC parameter list, the user specifies the following parameters:

 * "Start Time" the start time of the simulation
 * "End Time" the end time of the simulation
 * "End Cycle" the end cycle of the simulation 
 * "Flow model" specifies the choice of flow model.  The choices are currently `"Darcy"` for saturated flow and  `"Richards"` for unsaturated flow.
 * "disable Flow_PK" is used to disable flow in the simulation. In this case the user should specify a mesh block wise constant darcy flow field in the State namelist.
 * "disable Transport_PK" is used to disable transport in the simulation.
 * "disable Chemistry_PK" is used to disable chemistry in the simulation.
 * "Viz dump cycle frequency" is used to generate visualization dumps every so many cycles.
 * "Viz dump time frequency" is used to generate visualization dumps every so many time increments.

The sublist named "CGNS" is used to specify the filename for the CGNS visualization dumps. 



9. Transport
=======================================

Example:


.. code-block:: xml

  <ParameterList name="Transport">
    <Parameter name="CFL" type="double" value="0.5"/>   
    <!-- debug and developers options -->
    <Parameter name="enable internal tests" type="string" value="no"/>   
    <Parameter name="internal tests tolerance" type="double" value="1e-6"/>   
    <Parameter name="verbosity level" type="int" value="0"/>  
    <Parameter name="maximal time step" type="double" value="10"/>
    <!-- end of debug and developers options -->
    
    <ParameterList name="Transport BCs">
      <Parameter name="number of BCs" type="int" value="1"/>
      <ParameterList name="BC 0">
	<Parameter name="Side set ID" type="int" value="3"/>
	<Parameter name="Type" type="string" value="Constant"/>
	<Parameter name="Component 0" type="double" value="1.0"/>
	<Parameter name="Component 1" type="double" value="0.6"/>
	<Parameter name="Component 2" type="double" value="0.2"/>
      </ParameterList>  
      <ParameterList name="BC 1"> 
        <Parameter name="Side set ID" type="int" value="50001"/> 
        <Parameter name="Type" type="string" value="Constant"/> 
        <Parameter name="Component 0" type="double" value="1.0"/> 
        <Parameter name="Component 2" type="double" value="1.0e-4"/>     
      </ParameterList>
    </ParameterList>
  </ParameterList>

In the Transport sublist the user specifies the following parameters: 
 
 * "CFL" is the Courant–Friedrichs–Lewy number. It must be strictly bigger than 0 and less or equal to 1. Default value is 1. 
 * "enable internal tests" turns on/off build-in tests. This option is useful for code development; therefore its default value is "no". 
 * "verbosity level" sets up the volume of information printed out by the transport. It must be any non-negative integer. This option is useful for code development; therefore, its default value is 0.
 * "internal tests tolerance" is the relative tolerance for internal tests. This is the developers option. Default value is 1e-6.
 * "maximal time step" overwrites the calculated time step. This is the developers option.  
 	 
The boundary conditions sublist consists of a few similar sublists related to boundary side sets. The number of these sublists can be both bigger or smaller than the number of defined side sets. Each of the sublists may contain only a few components. The other components will be automatically set to zero. Note that the boundary conditions have to be set up mathematically only on influx boundaries. If it is not done, the default boundary condition is always zero.   

 * "number of BCs" is the total number of boundary conditions (i.e. subsequent sublists). 
 * "Side set ID" is the side set id in the mesh model. 
 * "Type" specifies the boundary condition. At the moment only constant boundary conditions are available. Put a ticket if you need a different type of boundary condition. 
 * "Component X" specified the value of component X on this boundary. 


10. Flow
=======================================

The flow parameter list

.. code-block:: xml

  <ParameterList name="Flow">
    ...
  </ParameterList>



specifies the parameters required by the flow process kernel.  This includes
numerical solver parameters and the specification of flow boundary conditions.
This parameter list is required if flow is enabled in the MPC parameter list
with the "disable Flow_PK" parameter.
[*This is ugly and ought to be changed to 'enabling'.*]


The following parameters must be specified in the  Flow parameter list:

* `"Max Iterations"` (int) defines the maximum number of iterations the
  flow solver is allowed to take.
  
* `"Error Tolerance"` (double) defines the error tolerance to which the
  flow solver will attempt to solve the flow equation.

* `"Nonlinear Solver"` (string) defines the choice of nonlinear solver when
  using the Richards flow model.  The valid  choices are `"JFNK"` for
  Jacobian-free Newton-Krylov and `"NLK"` for the Nonlinear Krylov method.
  This parameter is unused for the Darcy flow model.
  
* `"Preconditioner Update Frequency"` (int) sets how frequently the
  preconditioner will be recomputed during the iterative nonlinear solve of
  the Richards flow model.  With the value 1 it will be recomputed every
  iteration, with 2 every other iteration, and so forth.
  This parameter is unused for the Darcy flow model.

* `"Flow BC"` (list) defines the boundary conditions for the flow equations.
  This list consists of 0 or more primitive BC sublists and a parameter that
  gives the number of sublists to expect.
  * `"number of BCs"`
  The number of these conditions that are specified is defined by the parameter named "number of BCs". The boundary condition sublists must be named "BC00", "BC01" and so forth. Each of these boundary condition sublists must contain the following paramters:

 * "Type" defines the boundary condition type, allowed values are "Darcy Constant", "Pressure Constant", "Static Head", or "No Flow".
 * "BC value" is the value that should be applied, its interpretation depends on the parameter "Type" above.
 * "Side set ID" is the ID number of the side set in the mesh where the boundary condition should be applied.

The default boundary condition is "No Flow". It is applied to all boundary faces that are in side sets that do not have a corresponding boundary condition sublist.


Example
-------

.. code-block:: xml

  <ParameterList name="Flow">

    <Parameter name="Max Iterations" type="int" value="100"/>
    <Parameter name="Error Tolerance" type="double" value="1.0e-13"/>
    <Parameter name="Nonlinear Solver" type="string" value="NLK"/>
    <Parameter name="Preconditioner Update Frequency" type="int" value="1"/>
    
    <ParameterList name="Flow BC">
      <Parameter name="number of BCs" type="int" value="2"/>
      <ParameterList name="BC00">
	<Parameter name="Type" type="string" value="Darcy Constant"/>
	<Parameter name="BC value" type="double" value="-1.0" />
	<Parameter name="Side set ID" type="int" value="3" />
      </ParameterList>  
      <ParameterList name="BC01">
	<Parameter name="Type" type="string" value="Pressure Constant"/>
	<Parameter name="BC value" type="double" value="0.0" />
	<Parameter name="Side set ID" type="int" value="1" />
      </ParameterList>  
     </ParameterList>
  </ParameterList>



10. Control
=======================================

Here's a list of remaining parameters under the general category of "control".

 * `"maximum time step`"
 * `"maximum simulation time`"
 * `"CFL`"
 * `"initial time step`"
 * `"max step size change fraction`"
 * `"fixed time step size`"
 * `"small time step size cutoff`"
 * `"gravity vector`"
 * `"number of coarse cells across domain`"
 * `"maximum refinement level`"
 * `"refinement ratio between AMR levels`"
 * `"interval between regrid`"
 * `"regrid on restart`"
 * `"grid efficiency`"
 * `"number of error buffer cells`"
 * `"maximum grid size`"
 * `"grid blocking factor`"
 * `"fixed grid file`"
 * `"checkpoint file prefix`"
 * `"checkpoint file interval`"
 * `"restart file`"
 * `"write checkpoint files`"
 * `"number CPUs used to write checkpoint files`"
 * `"plotfile prefix`"
 * `"plotfile interval`"
 * `"write plotfiles`"
 * `"number CPUs used to write checkpoint files`"
 * `"state ids in plotfile`"
 * `"derived variables in plotfile`"


