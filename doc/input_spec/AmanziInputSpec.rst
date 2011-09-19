========================================
Amanzi XML Input Specification
========================================

.. contents:: **Table of Contents**


Overview
========

The Amanzi simulator evolves a system of conservation
equations for reacting flow in porous media, as detailed in
the ASCEM report entitled "Mathematical Formulation Requirements and
Specifications for the Process Models`" (hereafter referred to
as the 'Model Requirements Document (MRD)'). The purpose of the present
document is to specify the data required to execute Amanzi.  This specification
should be regarded as a companion to the MRD, and parameterizations of
the individual submodels are consistent between Amanzi, the MRD and this
document. Where applicable, the
relevant sections of the MRD are indicated.


Preliminary Concepts
--------------------

Amanzi solves a set of parameterized models for multiphase flow in porous media.  An Amanzi simulation is specified by providing:

* values for a parameterized PDE-based transport model, including boundary and initial conditions, constituitive laws, and parameterized/phenomenological models for fluid and chemical sources and characterizations of the porous medium,

* parameters controlling the selection of key algorithmic options and output, 

* a description of the (discrete) state of the computational sytem, including a list of the independent variables and instructions for obtaining or generating the discrete mesh, and a characterization of the (parallel) computing environment.

The primary software interface to Amanzi is a compiled C++ function, and much of the input data required is communicated through a single `Teuchos::ParameterList <http://trilinos.sandia.gov/packages/docs/r7.0/packages/teuchos/doc/html/index.html>`_.
A ParameterList consist of a simple hierarchy of parameters and lists of parameters, and is constructed directly from a similarly structured XML file.  The Amanzi input specification is defined in terms of the XML file format
used to construct a `Teuchos::ParameterList <http://trilinos.sandia.gov/packages/docs/r7.0/packages/teuchos/doc/html/index.html>`_.

In practice, Amanzi is called by a "simulation coordinator" which manages the simulation instructions and orchestrates the flow of data.  A basic simulation coordinator is
provided with the Amanzi source code distribution.  This simple stand-alone coordinator can be used to drive a simple sequence of Amanzi runs, or can serve as a template for user-generated extensions supporting more intricate workflows.  


Model Characterization
~~~~~~~~~~~~~~~~~~~~~~

For each phase in the model system, Amanzi requires the specification of component fluids and the chemical solutes they contain.  If a component exists in multiple phases, a relationship is required to compute its phase distribution as a function of the state of the stystem.
Equation 2.11 of the MRD governs the conservation and transport of each component, where the volumetric flow rate has been specified via Darcy's law (equation 2.10).  For Darcy flow, the properties of the porous medium must be specified over the entire simulation domain.  All transported phases
require an appropriate set of boundary conditions along the edge of the simulation domain.  Source/sink terms and initial data are provided for each phase component and any solutes they contain.  Additionally, the extent of physical domain
is specified, along with its discrete representation (i.e., the mesh).

Categories of Output Data
~~~~~~~~~~~~~~~~~~~~~~~~~

Output data from Amanzi is currently organized into four specific groups:

* `"Observation Data`": values generated during a run that characterize a particular feature of interest in the solution.  Such data would be used, for example, when quantifying the "response" of Amanzi to a set of targeted parameter variations.  `"Observation data`" is evaluated according to instructions given via the input parameter list, and is assembled into a special data structure that is returned to the calling routine through Amanzi's simulation driver function arguments.

* `"Visualization Data`": field-data intended for post-processing, and may include various subsets and transformations of the simulation's state data or material properties.

* `"Checkpoint Data`": the complete set of information that is required to duplicate a given execution of Amanzi.  Typically `"Checkpoint Data`" is created at periodic intervals during a long run in order to facilitate a repeatable simulation restart capability and archiving procedures. Unlike the first two data groups, there are no user-settable parameters to control the contents of a `"Checkpoint`" file - the implementation will determine the quantity, format and precision of this data in order to guarantee reproducibility.

* `"Log Data`": running commentary on the performance and status of Amanzi's execution.  Typically such data is written to a C++ stream which may be directed to a pipe or file.  The amount and detail of log data is determined by a range of verbosity controls.

Generally, `"Visualization Data`" and `"Checkpoint Data`" consists of high-dimensional field data representing snapshots of the evolving discrete variables.  These are large datasets, relative to the other types, and are most often written to disk in a file format that allows a direct repesentation of the underlying discrete mesh and parallel data distribution.


ParameterList XML
-----------------

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Input specification for each ParameterList entry in the input hierarchy consists of two parts.  First, a bulleted list defines the complete set of options available.  This is followed by example snipets of XML code to demonstrate usage.

* In many cases, Amanzi supports multiple parameterized models for a particular process.  This will be indicated in the specification by using the keyword `"MODEL(<prefix>)`".  A list of supported models is provided at the end of the section.  Each model will be given in the XML as a sublist labeled "<prefix>: STR" and the sublist will contain values for each of the specified parameters.  For example, the specification might be listed as:


 * `"Material Properties`" [list] 

  * MODEL(Porosity)

  * `"mass density`" [double]

  Here [list] indicates that this must be a ParameterList.  This specifcation will be followed by a list of valid models:

  * `"Porosity: Uniform`" requires `"Value`" [double] 
  * `"Porosity: GSLib`" requires `"Filename`" [string] 

Example:

.. code-block:: xml

    <ParameterList name="Material Properties">
      <ParameterList name="Porosity: Uniform">
        <Parameter name="Value" type="double" value="0.7"/>
      </ParameterList>   
    </ParameterList>   
 

Notation:

* Reserved keywords and labels are `"quoted and italicized`" -- these labels or values of parameters in user-generated input files must match (using XML matching rules) the specified or allowable values.  User-defined labels are indicated with ALL-CAPS, and are meant to represent a typical name given by a user - these can be names or numbers or whatever serves best the organization of the user input data.

* Where applicable, the relevant section of the MRD is referred to by section or chapter number in parentheses.



Version
=======

Each input set contains at the top level a string variable `"Amanzi Input Format Version`".  As of the most recent update of this specification, the
current version of the Amanzi input is `"1.0.0`".  If the version is unspecified, it is assumed to be earlier than `"0.9.0`".  Release notes documenting the
evolving input specification version can be found *here*.

* "Amanzi Input Format Version" [string] Three part version string

Example:

.. code-block:: xml

  <ParameterList name="Main">
    <Parameter name="Amanzi Input Format Version" type="string" value="1.0.0"/>
  </ParameterList>

Documentation
=============

The `"Documenation`" parameter list can be used to provide a brief description of the problem specified in the file.  Any number of string entries can be provided
with any label that may be useful for the user own purposes

* LABEL [string] A descriptive string

Example:

.. code-block:: xml

  <ParameterList name="Main">
    <ParameterList name="Documentation">
      <Parameter name="Simulation Objective" type="string" value="Validate workflow for parameter estimation"/>
      <Parameter name="Spatial Dimension" type="string" value="2"/>
      <Parameter name="Domain Shape" type="string" value="Rectangle: 2x1 aspect ratio"/>
      <Parameter name="Author" type="string" value="M. Day"/>
    </ParameterList>
  </ParameterList>



Mesh
=======================================

Amanzi supports a number of mesh "frameworks" used to discretize the simulation domain, including support for structured and unstructured grids.  The structured-grid option supports dynamic solution-adaptive grid generation.  Amanzi's unstructured grid options include variants that generate meshes internally "on-the-fly", and others that require the user to specify an externally-generated mesh.

Generally, the set of options for the mesh frameworks depend on whether the grid is to be generated or read in from a file.


Amanzi-generated grids:

* FRAMEWORK [list] labeled after mesh framework, accepts the following types: `"Structured-grid`", `"SimpleMesh`", `"stk::mesh`"

 * `"Domain Low Corner`" [Array double] Location of low corner of box

 * `"Domain High Corner`" [Array double] Location of high corner of box

 * `"Number Of Cells`" [Array int] the number of uniform cells in each coordinate direction

.. code-block:: xml

   <ParameterList name="Mesh">
     <ParameterList name="Structured"/>
       <Parameter name="Number of Cells" type="Array int" value="{100, 1, 100}"/>
       <Parameter name="Domain Low Corner" type="Array double" value="{0.0, 0.0, 0.0}" />
       <Parameter name="Domain High Corner" type="Array double" value="{103.2, 1.0, 103.2}" />
     </ParameterList>   
   </ParameterList>

Pre-generated grids:

* `"Framework`" [string] labeled after mesh framework, accepts the following types: `"MOAB`", `"Exodus`"

 * `"File`" [string] name of pre-generated mesh file

 * `"Format`" [string] format of pre-generated mesh file

Example

.. code-block:: xml

  <ParameterList name="Mesh">
    <ParameterList name="MOAB">
      <Parameter name="File" type="string" value="moab_filename"/>
      <Parameter name="Format" type="string" value="moab_default"/>
    </ParameterList>   
  </ParameterList>


Regions
=======================================

Regions are used in Amanzi to define subsets of the computational domain in order to specify the problem
to be solved, and the output desired.  Amanzi automatically defines the special region labeled `"All`", which is the 
entire simulation domain.  The user must additionally define the boundary surface(s) which enclose the domain.
Amanzi assumes that the union of the boundary surfaces envelopes the entire computational domain
(*i.e.* is "water-tight").  The special regions (`"All`" and the boundaries) may also serve as generic
regions (see the dicussion below for how these regions are labeled) and
can thus be used to specify other components of the problem (source terms, initial conditions, etc).

For the mesh framework options that support the `"Generate`" keyword, Amanzi implicitly defines the planes bounding the domain as regions that
automatically available to the input specification, using the following labels: `"XLOBC`", `"XHIBC`", `"YLOBC`", `"YHIBC`", `"ZLOBC`", `"ZHIBC`"

Regions specifications take the following form

 * "Regions" (list) can accept a number of lists for named regions (REGION)

   * MODEL(Region)

Amanzi supports parameterized forms for a number of analytic shapes, as well as more complex
definitions based on triangulated surface files.  

+------------------------+-----------------------------------------+------------------------------+--------------------------------------------------------------+
|  shape functional name | parameters                              | type(s)                      | Comment                                                      |
+========================+=========================================+==============================+==============================================================+
| `"Point"`              | `"Coordinate`"                          | Array double                 | Location of point in space                                   |
+------------------------+-----------------------------------------+------------------------------+--------------------------------------------------------------+
| `"Box"`                | `"Low Coordinate`", `"High Coordinate`" | Array double, Array double   | Location of boundary points of box                           |
+------------------------+-----------------------------------------+------------------------------+--------------------------------------------------------------+
| `"Plane"`              | `"Direction`", `"Location`"             | Array double, Array double   | Location of boundary points of box                           |
+------------------------+-----------------------------------------+------------------------------+--------------------------------------------------------------+
| `"Labeled Set"`        | `"label`", `"file`",                    | string, string,              | Set per label defined in mesh file (see below)               |
|                        | `"mesh framework`", `"entity`"          | string, string               |  (available for frameworks supporting the `"File`" keyword)  |
+------------------------+-----------------------------------------+------------------------------+--------------------------------------------------------------+
| `"Layer"`              | `"file#`", `"label#`"                   | (#=1,2) string, string       | Region between two surfaces                                  |
+------------------------+-----------------------------------------+------------------------------+--------------------------------------------------------------+
| `"Surface"`            | `"file`" `"label`"                      | string, string               | Labeled triangulated face set in file                        |
+------------------------+-----------------------------------------+------------------------------+--------------------------------------------------------------+

Notes

* `"Box`" and "Plane" must be bounded by coordinate-aligned lines and planes.

* The "Labeled Set" region is defined by a label that was given to sets generated in a preprocessing step and stored in a mesh-dependent data file.  For example, an "exodus::mesh" type mesh file can be processed to tag cells, faces and/or nodes with specific labels, using a variety of external tools.  Regions based on such sets are assigned a user-defined label for Amanzi, which may or may not correspond to the original label in the exodus file.  Note that the file used to express this labeled set may be in any Amanzi-supported mesh framework (the mesh framework is specified in the parameters for this option).  The `"entity`" parameter may be necessary to specify a unique set.  For example, an exodus file requires `"Cell`", `"Face`" or `"Node`" as well as a label (which is an integer).  When the mesh framework for the region is different from the current mesh framework (defined in `"Mesh`" above), the intersection of the specified region and the computational domain defines the region.  This latter option is not yet supported, but will likely be implemented as a special (piecewise-constant) case of a generalized interpolation operator.

* Surface files contain labeled triangulated face sets.  The user is responsible for ensuring that the intersections with other surfaces in the problem, including the boundaries, are `"exact`" (*i.e.* that surface intersections are `"watertight`" where applicable), and that the surfaces are contained within the computational domain.  If nodes in the surface fall outside the domain, the elements they define are ignored.

* Eventually, Amanzi will support a "geometric modeling" syntax such that complex regions can be assembled by composition with logical operators.  The next step toward this capability will likely be to allow the definition of a single region as a concatentation of a number of basic shapes.  A more general capability might include the name of an instruction file (and a label to identify a particular region in the file) to interface to a scripted modeler.

Example:

.. code-block:: xml

  <ParameterList name="Regions">
    <ParameterList name="Top Section">
      <ParameterList name="Box">
        <Parameter name="Low Coordinate" type="Array double" value="{2, 3, 5}"/>
        <Parameter name="High Coordinate" type="Array double" value="{4, 5, 8}"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="Middle Section">
      <ParameterList name="Box">
        <Parameter name="Low Coordinate" type="Array double" value="{2, 3, 3}"/>
        <Parameter name="High Coordinate" type="Array double" value="{4, 5, 5}"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="Bottom Section">
      <ParameterList name="Box">
        <Parameter name="Low Coordinate" type="Array double" value="{2, 3, 0}"/>
        <Parameter name="High Coordinate" type="Array double" value="{4, 5, 3}"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>

In this example, "Top Section", "Middle Section" and "Bottom Section" are three box-shaped regions.



Soil Properties
====================

Soil properties must be specified over the entire simulation domain (`"All`") defined in the Region section.  This can be implemented using any combination of regions
defined above, provided that the entire domain is covered.  The regions used should be disjoint.  Each soil type (Section 2.6) is given a label (string) and assigned
physical properties using from a selection of models.  A model can be as simple as `"Porosity: Uniform`", which sets the porosity in every cell to a single value, or it may be extremely 
complex.  The available models for each property are listed below the specification.  Each soil that is defined is assigned to a list of regions.

* "Soil Properties" (list) can accept multiple lists for named soil types (SOIL)

 * SOIL (list) can accept lists to specify models, and `"Assigned Regions`" to specify where this model applies

  * MODEL(Porosity)

  * MODEL(Mass Density)

  * MODEL(Intrinsic Permeability)

  * MODEL(Capillary Pressure)

  * `"Assigned Regions`" (Array string) a set of labels corresponding to defined regions

The following models are currently supported for porosity:

* `"Porosity: File`" requires the following strings: `"File`" (name of a file), `"Label`" (the label of the scalar field in the file to associate with the values of porosity).  Optionally `"Interpolation`" (the interpolation strategy: : `"Constant`" [default] or `"Linear`").  Optionally `"Framework`" (if the mesh framework with which the file was written is different from current) will indicate the format of the file.  Note that the physical domain of this input data must completely cover the union of the regions over which this property is to be evaluated.

* `"Porosity: Uniform`" requires `"Value`" [double] to specify the constant value of porosity.

* `"Porosity: Random`" requires the `"Mean And RMS Values`" [Array double]

* `"Porosity: GSLib`" requires `"File`" [string], the name of a gslib input file 


The following models are currently supported for mass density:

* `"Mass Density: File`" requires the following strings: `"File`" (name of a file), `"Label`" (the label of the scalar field in the file to associate with the values of mass density).  Optionally `"Interpolation`" (the interpolation strategy: : `"Constant`" [default] or `"Linear`").  Optionally `"Framework`" (if the mesh framework with which the file was written is different from current) will indicate the format of the file.  Note that the physical domain of this input data must completely cover the union of the regions over which this property is to be evaluated.

* `"Mass Density: Uniform`" requires `"Value`" [double] to specify the constant value of mass density of the soil.


The following models are currently supported for the absolute (soil) permeability:

* `"Intrinsic Permeability: File`" requires the following strings: `"File`" (name of a file), `"Label`" (the label of the scalar field in the file to associate with the values of intrinsic permeability).  Optionally `"Interpolation`" (the interpolation strategy: : `"Constant`" [default] or `"Linear`").  Optionally `"Framework`" (if the mesh framework with which the file was written is different from current) will indicate the format of the file.  Note that the physical domain of this input data must completely cover the union of the regions over which this property is to be evaluated.

* `"Intrinsic Permeability: Uniform`" requires `"Value`" [double] to specify the constant value of the intrinsic permeability

* `"Intrinsic Permeability: Random`" requires the `"Mean And RMS Values`" [Array double]

* `"Intrinsic Permeability: GSLib`" requires `"File`" [string], the name of a gslib input file 

* Additionally, all intrinsic permeability models optionally accept the following parameter:

  * `"Intrinsic Permeability Anisotropy`" [Array double] - (optional) indicates that the intrinsic permeability is a diagonal tensor, an the XX, YY, and ZZ are given by the specifed X value and scaled by these values.


The following models are currently supported for relative permeability (Section 2.6):

* `"Relative Permeability: Perfect`" requires no parameters, krl=krg=1

* `"Relative Permeability: Linear`" requires no parameters, krl=sl and krg=sg

* `"Relative Permeability: Quadratic`" requires slr, sgr [Array double]

* `"Relative Permeability: vGM`" (van Genuchten-Mualem) requires m, slr, sgr [Array double]

The following models are currently supported for capillary pressure (Section 3.3.2):

* `"Capillary Pressure: None`" requires no parameters, pc = 0

* `"Capillary Pressure: Linear`" requires no parameters, pc = sl

* `"Capillary Pressure: vG`" requires m, sigma, slr, sgr [Array double]

Example:

.. code-block:: xml

  <ParameterList name="Soil Properties">
    <ParameterList name="Backfill">
      <ParameterList name="Mass Density: Uniform">
        <Parameter name="Value" type="double" value="2.8e3"/>
      </ParameterList>
      <ParameterList name="Intrinsic Permeability: Uniform">
        <Parameter name="Value" type="double" value="1240"/>
        <Parameter name="Permeability Anisotropy" type="Array double" value="{1., 0.001, 0.001}"/>
      </ParameterList>
      <ParameterList name="Porosity: Uniform">
        <Parameter name="Value" type="double" value="0.2585"/>
      </ParameterList>
      <ParameterList name="Relative permeability: vGM">
        <Parameter name="m_slr_sgr" type="Array double" value="{0.6585, 0.0774, 0}"/>
      </ParameterList>
      <ParameterList name="Capillary Pressure: vG">
        <Parameter name="m_sigma_slr_sgr" type="Array double" value="{0.6585, 102.1, 0.0774, 0}"/>
      </ParameterList>
      <Parameter name="Assigned regions" type="string array" value="{Top Region, Bottom Region}"/>
    </ParameterList>
    <ParameterList name="Fine Sand">
      <ParameterList name="Mass Density: Uniform">
        <Parameter name="Value" type="double" value="2.8e3"/>
      </ParameterList>
      <ParameterList name="Intrinsic Permeability: Uniform">
        <Parameter name="Value" type="double" value="337"/>
      </ParameterList>
      <ParameterList name="Porosity: Uniform">
        <Parameter name="Value" type="double" value="0.3586"/>
      </ParameterList>
      <ParameterList name="Relative Permeability: vGM">
        <Parameter name="m_slr_sgr" type="Array double" value="{0.4694, 0.0837, 0}"/>
      </ParameterList>
      <ParameterList name="Capillary Pressure: vG">
        <Parameter name="m_sigma_slr_sgr" type="Array double" value="{0.4694, 9.533, 0.0837, 0}"/>
      </ParameterList>
      <Parameter name="Assigned Regions" type="string array" value="{middle}"/>
    </ParameterList>
  </ParameterList>

In this example, there are two types of soil, `"Backfill`" (which fills `"Bottom Region`" and `"Top Region`") and `"Fine Sand`" (which fills `"Middle Region`").  Both have
van Genuchten models for relative permeability and capillary pressure.  `"Backfill`" has an anisotropic permeability, where the vertical value is 1000 times
the horizontal values.




Mobile Phases
=======================================

The `"Mobile Phases`" parameter list is used to specify components of each of the phases that are mobile, and solutes that are contained within them.  For each such 
phase, the list identifies the set of all independent variables that are to be stored on each discrete mesh cell.
For organizational convenience, the `"Mobile Phases`" parameter list is also where the initial conditions, boundary data and source
terms are defined for each phase component.  Future versions of Amanzi will support mass transfer between phases, and this is also where
the phase distribution models will be specified.

Phases, components and solutes
------------------------------

In the general problem, multiple phases may coexist in the domain (e.g. gaseous, aqueous, etc), and each is
comprised of a number of components (section 2.2).  In turn, each component may carry a number of solutes and some of these may participate
in chemical reactions.  As a result of reactions, a chemical source or sink term may appear for the solutes involved in the reaction, including solutes in other mobile phases or in the soil matrix.  
Additionally, certain reactions such as precipitation may affect the flow properties of the soil itself during the simulation, and 
some might affect the properties of the fluid (e.g. brines affect the liquid density). While Amanzi does not currently support chemical reactions and thermal processes, the specification here allows for the existence of
the necessary data structures and input data framework.

Currently in Amanzi, inert solutes are transported in the various phase components and are treated in "complexes".  Each complex is in chemical equilibrium with itself and does not undergo phase change.
Under these conditions, knowledge of the local concentration of the "basis" or "primary" species (the terms are used here interchangeably) in a chemical complex is sufficient to determine the concentrations of all related secondary species
in the phase. Each basis species has a total component concentration and a free ion concentration. The total component concentration for each basis species is a sum of the
free ion concentrations in the phase components and its stoichiometric contribution to all secondary species. Amanzi splits the total component concentration into a set of totals for each of the transported phases
and total sorbed concentration. Given the free ion concentration of each basis species (and if there is more than one phase, a specification of the 
equilibrium phase distribution of components that appear in more than one phase), we can reconstruct the concentration of the secondary species in each phase. As a result only the basis species are maintained in the state
data structures for each phases component.

In addition to solutes in the transported phases, there may be various immobile chemical constituents within the
porous media (soil) matrix, such as "minerals" and "surface complexes". Bookkeeping for these constituents is managed in Amanzi
data structures by generalizing the "solute" concept - a slot in the state is allocated for each of these immobile species, but their concentrations are
not included in the transport/flow components of the numerical integration.  To allow selective transport of the various solutes, Amanzi
uses the concept of solute groups.   The aqueous solute concentrations are typically treated together as a group, for example, and often represent the only 
chemical constituents that are mobile.

Specification of Amanzi's numerical state is organized fundamentally around the list of phases that are present.  Each phase consists of multiple components.  For each of these,
Amanzi requires a label, a set of models that specify its physical properties (Section 4.6), and a list of solutes.  For each solute, a group membership is specified.
Note that Amanzi will eventually support the use of a master chemistry database, where the solute complexes and their chemical activity are defined.  In that case, inclusion of a particular solute in the
Amanzi input file will be conditioned on its presence in the appropriate section of the master list.

Sources and Initial and Boundary Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Mobile phase components, and solutes contained in them, require boundary conditions along the entire surface bounding the computational domain (Sections 3.3, 3.6, 3.10 and 4.3).  Generally, phase component boundary conditions are
specified in porous media systems by giving either the component pressure or Darcy velocity on the boundary, along with the phase saturation on the bounding surface.  Since mobile solutes are carried with the resulting flow,
inflowing boundary conditions for solutes are typically specified using Dirichlet conditions that define the effective solute concentration in the incoming flow.  On outflow boundaries,
no solute information is carried into the domain so no data is required. For simplicity here, any boundary conditions not explicitly set in the input are defaulted to outflow.

Volumetric source terms, used to model infiltration (Section 3.7) and a wide variety of production and loss processes, are defined for each phase component, if applicable, and include the distribution of any solutes that are carried into the domain with the phase component.

Boundary conditions and source terms may be time-dependent, in general.

The generalized specification is as follows:

* `"Mobile Phases`" (list) can accept lists named phases (PHASE).

 * PHASE (list) can accept the following lists: `"Phase Properties`", `"Phase Components`"

  * `"Phase Properties`" can accept models for viscosity and density

   * MODEL(Density)

   * MODEL(Viscosity)

  * `"Phase Components`" can accept COMP [list] named after a user-defined phase component.

   * COMP (list) can accept `"Solute Properties`" [list] to define solutes carried by the component.  Also, accepts`"Component Initial Conditions`" [list], `"Component Boundary Conditions`" [list], `"Component Sources`" [list]

    * `"Component Initial Conditions`" (list) accepts lists IC-REGION named after the user-defined region that IC function will apply over

     * IC-REGION (list) can accept a model for initial conditions, list for solute initial conditions

      * MODEL(Initial Conditions)

      * `"Solute Initial Conditions`" can accept lists SOLUTE named after indivual solutes

       * SOLUTE can accept a model for initial conditions, and a flag for the units of Dirichlet values in the model

        * MODEL(Initial Conditions)

        * `"Concentration Units`" [string] can accept `"Molar Concentration`" (moles/volume), `"Molal Concentration`" (moles/volume of water) , `"Specific Concentration`" (mass/volume of water)

    * `"Component Boundary Conditions`" (list) accepts lists BC-REGION named after the user-defined region that BC function will apply over

     * BC-REGION (list) can accept a model for boundary conditions, and list for solute booundary conditions

      * MODEL(Boundary Conditions)

      * `"Solute Boundary Conditions`" can accept lists SOLUTE named after indivual solutes

       * SOLUTE can accept a model for boundary conditions, and a flag for the units of Dirichlet values in the model

        * MODEL(Boundary Conditions)

        * `"Concentration Units`" [string] can accept `"Molar Concentration`" (moles/volume), `"Molal Concentration`" (moles/volume of water) , `"Specific Concentration`" (mass/volume of water)

    * `"Component Boundary Sources`" (list) accepts lists S-REGION named after the user-defined region that source function will apply over

     * S-REGION (list) can accept a model for a source, and list for solute sources

      * MODEL(Source)

      * `"Solute Source`" can accept lists SOLUTE named after indivual solutes

       * SOLUTE can accept a model for a source, and a flag for the units of Dirichlet values in the model

        * MODEL(Source)

        * `"Concentration Units`" [string] can accept `"Molar Concentration`" (moles/volume), `"Molal Concentration`" (moles/volume of water) , `"Specific Concentration`" (mass/volume of water)

Initial conditions are required for each phase component, and the solutes contained in them, over the entire computational domain.
Boundary conditions are required on all domain boundaries (see Sections 3.3, 4.3).  Source terms for all are optional.  All are constructed using a limited number
of explicitly parameterized model are supported for communicating initial conditions:

* `"Initial Conditions: Uniform`" requires `"Value`" [double]

* `"Initial Conditions: Gradient`" requires `"Reference Coordinate`" (Array double), `"Reference Value`" [double], and  `"Gradient Value`" (Array double)

* `"Initial Conditions: File`" requires `"File`" [string] and `"Label`" [string] - the label of the field to use.  If the file format is not compatible with the current mesh framework, `"Format`" [string] is also required.

The following parameterized boundary conditions are supported for communicating boundary conditions:

* `"Boundary Conditions: Flux`" requires `"Times`" [Array double], `"Time Functions`" [Array string] (see the note below) and one of the following: `"Extensive Volumetric Flux`" [double] or `"Extensive Mass Flux`" [double], `"Intensive Volumetric Flux`" [double] or `"Intensive Mass Flux`" [double]

* `"Boundary Conditions: Uniform Pressure`" requires `"Times`" [Array double], `"Time Functions`" [Array string] and `"Values`" [Array double]

* `"Boundary Conditions: Seepage`" requires `"Times`" [Array double], `"Time Functions`" [Array string] and `"Water Table Height`" [double] (see below)

* `"Boundary Conditions: Hydrostatic`" requires `"Times`" [Array double], `"Time Functions`" [Array string] and `"Water Table Height`" [double] (see below)

* `"Boundary Conditions: Impermeable`" requires no data

* `"Boundary Conditions: Outflow`" requires no data

The following models are currently supported for communicating source distribution:

* `"Source: Uniform Volumetric Rate`" requires `"Times`" [Array double], `"Time Functions`" [Array string], and `"Values`" [Array double].  

* `"Source: Uniform Mass Rate`" requires `"Times`" [Array double], `"Time Functions`" [Array string],  `"Values`" [Array double].  

Time Functions
~~~~~~~~~~~~~~

Boundary data and source models utilize a parameterized model for time variations that is either piecewise constant or piecewise linear.  For example:

.. code-block:: xml

      <Parameter name="Times" type="Array double" value="{1, 2, 3}"/>
      <Parameter name="Time Values" type="Array double" value="{10, 20, 30}"/>
      <Parameter name="Time Functions" type="Array string" value="{Constant, Linear}"/>    


This define four time intervals: (-inf,1), (1,2), (2,3), (3,+inf).  By assumption the function is constant over the first and last intervals.  The remaining 
two intervals are speicified by the `"Time Functions`" parameter.  Thus, the value here is 10 anytime prior to t=2. The value increases linearly from 10 to 
20 over the interval t=2 to t=3, and then is constant at 30 for t>3.


Example Phase Definition
~~~~~~~~~~~~~~~~~~~~~~~~
Due to its length, an XML example of the `"Mobile Phases`" parameter list appears in the attached file:XXX.


Output
======

Output data from Amanzi is currently organized into four specific groups: `"Observation Data`", `"Visualization Data`", `"Checkpoint Data`" and `"Log Data`".  
Each of these is controlled in different ways, reflecting their intended use.

* `"Checkpoint Data`" is intended to represent all that is necesary to repeat or continue an Amanzi run.  The specific data contained in a checkpoint dump is specific to the algorithm optoins and mesh framework selected.  Checkpoint data is special in that no interpolation is perfomed prior to writing the data files; the raw binary state is necessary.  As a result, the user is allowed to only write checkpoint information at the discrete intervals of the simulation.

* `"Visualization Data`" is intended to represent spatially complete snapshots of the solution at defined instances during the simulation.  Dependeing on the control parameters provided here, visualizatoin files may include only a fraction of the state data, and may contiain auxiliary "derived" information (see below for more discussion).

* `"Observation Data`" is intended to represent diagnostic values to be returned to the calling routine from Amanzi's simulation driver.  Observations are typically generated at arbitrary times, and frequently involve various point samplings and volumetric reductions that are interpolated in time to the desired instant.  Observations may involve derived quantities (see discussion below) or state fields.

* `"Log Data`" is intended to represent runtime diagnostics to indicate the status of the simulation in progress.  This data is typically written by the simulation code to the screen or some other stream or file pipe.  The volume of `"Log Data`" generated is typically a function of various verbosity settings for a given run.

"`Log Data`" is not explicitly controlled in this section, since it is easier to control in the context of specifying details of the algorithms.  The remaining data types are discussed in the section below.


Observation Data
----------------

A user may request any number of specific observations from Amanzi.  Each labeled observation involves a state or derived component, a model, a region from which it will extract its source data, and a list of discrete times 
for its evaluation.  The observations are evaluated during the simulation and returned to the calling process through one of Amanzi arguments.

* `"Observation Data`" [list] can accept multiple lists for named observations (OBSERVATION)

  * OBSERVATION [list] user-defined label, can accept values for `"Phase`", `"Component`", `"Solute`", `"Region`", `"Times`" and a model.

    * `"Phase`" [string] the label of a phase defined above

    * `"Component`" [string] the label of one of the components defined for this phase

    * `"Region`" [string] the label of a user-defined region

    * `"Solute`" [string] (optional) the label of one of the solutes defined for this phase component

    * `"Times`" [Array double] values of time where this quantity is desired

    * MODEL(Observation)

The following observation models are currently supported.  All of them operate on the state quantity identified.
* `"Observation: Mean`" returns the mean value of the phase or component saturation, or the solute concentration over the region
* `"Observation: Integral`" returns the integral of the phase or component saturation, or the solute concentration over the region
* `"Observation: Flux Integral`" returns the integral of the flux of the phase, component, or solute over the region
* `"Observation: Peak Value`" returns the peak value of the phase or component saturation, or the solute concentration over the region
* `"Observation: Distance to Center of Mass`" returns the distance from a given location of the center of mass of the phase or component saturation, or the solute concentration over the region.  Requires a single parameter, "Reference Location" [Array double] specifying the refnerece location.

Example:

.. code-block:: xml

  <ParameterList name="Observation">
    <ParameterList name="Center of UO+2 Mass">
      <Parameter name="Phase" type="string" value="Aqueous"/>
      <Parameter name="Component" type="string" value="Water"/>
      <Parameter name="Solute" type="string" value="UO+2"/>
      <Parameter name="Region" type="string" value="All"/>
      <ParameterList name="Observation: Distance to Center of Mass">
        <Parameter name="Reference Location" type="Array double" value="{0, 0, 100}"/>
      </ParameterList>
      <Parameter name="Times" type="Array double" value="{10, 30 , 50}">
    </ParameterList>
  </ParameterList>

In this example, the user requests that the center of mass for the solute UO+2 be computed, and that the distance from that location to the point (0, 0, 100) be returned at t={10, 30 and 50}.
The format of the data structure used to communicate the observation data back to the calling function includes a flag for each requested time to indicate whether the quantity was successfully filled.


Checkpoint Data
---------------------------------

A user may request periodic dumps of Amanzi checkpoint data.  The user has not explicit control over the content of these files, but has the guarantee that the Amanzi run will be reproducible (with accuracies determined
by machine round errors and randomness due to execution in a parallel computing environment.  Therefore, output controls for checkpoint data are limited to file name generation and writing frequency, by numerical cycle number.

* `"Checkpoint Data`" [list] can accept a file name base [string], cycle data [list] and number of digits [int] used to generate the file name (e.g. file00050 if 5 digits and time step 50)

  * `"File Name Base`" [string]

  * `"Cycle Data`" [string] can accept start, end and interval data for cycle number

    * `"Start`" [int] step number of first file

    * `"End`" [int] step number of last file, if < 0 or not present then value is not used (no stopping condition)

    * `"Interval`" [int] number of steps per file write

    * `"Steps`" [Array int] specific step numbers to write (if parameter present, the (Start, Step, Interval) ignored

Example:

.. code-block:: xml

  <ParameterList name="Checkpoint Data">
    <Parameter name="File Name Base" type="string" value="chk"/>
    <Parameter name="File Name Digits" type="int" value="5"/>
    <ParameterList name="Cycle Data">
      <Parameter name="Start" type="int" value="0"/>
      <Parameter name="End" type="int" value="-1"/>
      <Parameter name="Interval" type="int" value="5"/>
    </ParameterList>
  </ParameterList>

In this example, checkpoint data is written when the cycle number is evenly divisble by 5.


Visualization Data
---------------------------------

A user may request periodic writes of field data for the purposes of vizualization.  The user will specify explicitly what is to be included in the file at each snapshot.  Visualization files can only be written 
at intervals corresponding to the numerical time step values; writes are controlled by timestep cycle number.

* `"Visualization Data`" [list] can accept a file name base [string], cycle data [list] and number of digits [int] used to generate the file name (e.g. file00050 if 5 digits and time step 50).  It can also accept a set of lists to specify which state variables to write.

  * `"File Name Base`" [string]

  * `"Cycle Data`" [string] can accept start, end and interval data for cycle number

    * `"Start`" [int] step number of first file

    * `"End`" [int] step number of last file, if < 0 or not present then value is not used (no stopping condition)

    * `"Interval`" [int] number of steps per file write

    * `"Steps`" [Array int] specific step numbers to write (if parameter present, the (Start, Step, Interval) ignored

  * `"Variable`" [list] can accept `"Phase`" [string], `"Component`" [string] (optional), `"Solute`" [string]

    * `"Phase`" [string] the label of a phase defined above, or "All" to write all phases

    * `"Component`" [string] the label of one of the components defined for this phase, or "All" to write all components of the selected phase(s)

    * `"Solute`" [string] the label of a solute defined above, or "All" to write all solutes of the component, or "None" to write none of them.



Example:

.. code-block:: xml

  <ParameterList name="Visualization Data">
    <Parameter name="File Name Base" type="string" value="chk"/>
    <Parameter name="File Name Digits" type="int" value="5"/>
    <ParameterList name="Cycle Data">
      <Parameter name="Start" type="int" value="0"/>
      <Parameter name="End" type="int" value="-1"/>
      <Parameter name="Interval" type="int" value="5"/>
    </ParameterList>
    <ParameterList name="Variable">
      <Parameter name="Phase" type="string" value="Aqueous"/>
      <Parameter name="Component" type="string" value="Water"/>
      <Parameter name="Solute" type="string" value="UO+2"/>
    </ParameterList>
    <ParameterList name="Variable">
      <Parameter name="Phase" type="string" value="Gas"/>
      <Parameter name="Component" type="string" value="All"/>
      <Parameter name="Solute" type="string" value="All"/>
    </ParameterList>
  </ParameterList>

In this example, visalization data is written when the cycle number is evenly divisble by 5.  The files will include the concentration of UO+2 in the Aqueous Water component, and all the solues in the Gas Phase.



Execution Control
=================

       This section is highly specific to the numerical algorithm details, which
       will be a sensitive function of the mesh framework, the type of problem 
       selected, the mode requested for time integration, whether the mesh
       is dynamically adaptive, and a host of more detailed algorithm and model
       decisions.  

       The parameter set below represents a fictional calculation and depicts 
       an organization of the numerical parameters that might be appropriate.       
       The main ParameterList here is named after a labeled "type" of solve
       one might like to do.  Had this been an unsteady simulation, many of the
       linear and nonlinear solver parameters may not be applicable at all.

       It is unclear whether the inputs for this section can or should be orgainized
       at any finer a level of granularity.

       See the example XML file for a typical set of control parameters.



Complete Example
=================

Presented below is a complete example of an Amanzi input file.  It does not exercise all the options provided for in this specification, but rather provides a concrete example of a set of self-consistent definitions
required to specify a real simulation with Amanzi envisioned functional for the Phase 2 demo deadline.

.. code-block:: xml

       <!--
          Simple test problem for variably saturated flow, solute transport, spatially-
          variable properties, and time-dependent boundary conditions. This example
          could represent a basic setup for a Hanford deep vadose zone problem. 
          Note, however, that the listed parameter values are not necessarily accurate for
          such an application. 
          
          Submitted by Mark Rockhold and Vicky Freedman, PNNL, September 6, 2011.
          Rearranged and generalized by Marc Day, 9/9/11, LBNL
          Futher rearranged and generalized by Marc Day to incorporate ascem-phi comments, 9/16/11, LBNL
         -->
       <ParameterList name="Main">  
       
         <Parameter name="Amanzi input format version" type="string" value="1.0.0"/>
       
         <ParameterList name="Documentation">
           <Parameter name="Model Name" type="string" value="Steady Richards"/>
           <Parameter name="Description" type="string" value="BC Cribs" />
         </ParameterList>
       
         <ParameterList name="Mesh">
           <ParameterList name="STK::mesh">
             <Parameter name="File" type="string" value="mesh.par" />
           </ParameterList>
         </ParameterList>
         
         <!-- Region Definitions -->
         <ParameterList name="Regions">    
           <ParameterList name="Ringold Region">
             <ParameterList name="Block">
               <Parameter name="Low Coordinate" type="Array double" value="{0.0, 0.0, 0.0"/>
               <Parameter name="High Coordinate" type="Array double" value="{100.0, 100.0, 20.0}"/>
             </ParameterList>
           </ParameterList>
           
           <ParameterList name="Caliche Region">
             <ParameterList name="Block">
               <Parameter name="Low Coordinate" type="Array double" value="{0.0, 0.0, 20.0}"/>
               <Parameter name="High Coordinate" type="Array double" value="{100.0, 100.0, 25.0}"/>
             </ParameterList>
           </ParameterList>
           
           <ParameterList name="Hanford Region">
             <ParameterList name="Block">
               <Parameter name="Low Coordinate" type="Array double" value="{0.0, 0.0, 25.0}"/>
               <Parameter name="High Coordinate" type="Array double" value="{100.0, 100.0, 100.0}"/>
             </ParameterList>
           </ParameterList>
           
           <ParameterList name="Top Surface region">
             <ParameterList name="Coordinate Plane">
               <Parameter name="Coordinate Direction" type="string" value="Z"/>
               <Parameter name="Coordinate Location" type="double" value="100.0"/>
             </ParameterList>
           </ParameterList>
           
           <ParameterList name="Bottom Surface Region">
             <ParameterList name="Coordinate Plane">
               <Parameter name="Coordinate Direction" type="string" value="Z"/>
               <Parameter name="Coordinate Location" type="double" value="0.0"/>
             </ParameterList>
           </ParameterList>
           
           <ParameterList name="Crib 1 Region">
             <ParameterList name="Block">
               <Parameter name="Low Coordinate" type="Array double" value="{20.0, 20.0, 97.0}"/>
               <Parameter name="High Coordinate" type="Array double" value="{23.0, 23.0, 100.0}"/>
             </ParameterList>
           </ParameterList>
       
           <ParameterList name="Crib 2 Region">
             <ParameterList name="Block">
               <Parameter name="Low Coordiante" type="Array double" value="{40.0, 40.0, 97.0}"/>
               <Parameter name="High Coordinate" type="Array double" value="{43.0, 43.0, 100.0}"/>
             </ParameterList>
           </ParameterList>
       
           <ParameterList name="Sample Point Region">
             <ParameterList name="point">
               <Parameter name="Coord" type="Array double" value="{50.0, 50.0, 0.0}"/>
             </ParameterList>
           </ParameterList>
         </ParameterList> <!-- End of Region Definitions -->
       
       
         <ParameterList name="Soil Properties">
       
           <ParameterList name="Ringold Material">
             <Parameter name="Assigned Regions" type="Array string" value="{Ringold Region}" />
             <ParameterList name="Porosity: Uniform">
               <Parameter name="Porosity" type="double" value="0.38" />
               <Parameter name="Porosity Units" type="string" value="Null" />
             </ParameterList>
             <ParameterList name="Intrinsic Permeability: Uniform">
               <Parameter name="Intrinsic Permeability" type="double" value="200" />
               <Parameter name="Intrinsic Permeability units" type="string" value="mD" />
             </ParameterList>
             <Parameter name="Intrinsic Permeability Anisotropy" type="Array double" value="{1., 1., 0.1}" />      
             <ParameterList name="Relative Permeability: Mualem">
               <Parameter name="Pore Interaction Term ???" type="Array double" value="{0.5, 0.5, 0.5}" />
               <Parameter name="Pore Interaction Term ??? Units" type="Array string" value="{Null, Null, Null}" />
             </ParameterList>
             <ParameterList name="Capillary Pressure: vG">
               <Parameter name="vG_alpha_n_Sr" type="Array double" value="{0.03, 2.7, 0.0234}" />
               <Parameter name="vG alpha units" type="Array string" value="{cm^-1, Null,  Null}" />
             </ParameterList>      
             <ParameterList name="Dispersivity: Uniform">
               <Parameter name="Dispersivity" type="double" value="0.5" />
               <Parameter name="Dispersivity units" type="string" value="cm"/>
             </ParameterList>    
             <Parameter name="Dispersivity Anisotropy" type="Array double" value="{10., 10., 1.0}" />
           </ParameterList>
           
           <ParameterList name="Caliche Material">
             <Parameter name="Assigned Regions" type="Array string" value="{Caliche Region}" />
             <ParameterList name="Porosity: Uniform">
               <Parameter name="Porosity" type="double" value="0.40" />
             </ParameterList>
             <ParameterList name="Intrinsic Permeability: Uniform">
               <Parameter name="Intrinsic Permeability" type="double" value="500" />
               <Parameter name="Intrinsic Permeability units" type="string" value="mD" />
             </ParameterList>      
             <Parameter name="Intrinsic Permeability Anisotropy" type="Array double" value="{1., 1., 0.1}" />
             <ParameterList name="Relative Permeability: vG">
               <Parameter name="Pore Interaction Term ???" type="Array double" value="{0.5, 0.5, 0.5}" />
               <Parameter name="Pore Interaction Term ??? Units" type="Array string" value="{}" />
             </ParameterList>
             <ParameterList name="Capillary Pressure: vG">
               <Parameter name="vG_alpha_n_Sr" type="Array double" value="{0.03, 2.7, 0.0234}" />
               <Parameter name="vG alpha units" type="Array string" value="{cm^-1, Null, Null}" />
             </ParameterList>
             <ParameterList name="Dispersivity: Uniform">
               <Parameter name="Dispersivity" type="double" value="0.5" />
               <Parameter name="Dispersivity units" type="string" value="cm"/>
             </ParameterList>    
             <Parameter name="Dispersivity Anisotropy" type="Array double" value="{10., 10., 1.0}" />
           </ParameterList>
       
           <ParameterList name="Hanford Material">
             <Parameter name="Assigned Regions" type="Array string" value="{Hanford Region}" />
             <Parameter name="Particle Density" type="double" value="2.65" />
             <Parameter name="Particle Density Units" type="string" value="g cm^-3"/>
             <ParameterList name="Porosity: Uniform">
               <Parameter name="Porosity" type="double" value="0.25" />
             </ParameterList>
             <Parameter name="Compressibility" type="double" value="1.e-6" />
             <Parameter name="Compressibility units" type="string" value="psi" />
             <ParameterList name="Intrinsic Permeability: Uniform">
               <Parameter name="Intrinsic Permeability" type="double" value="500" />
               <Parameter name="Intrinsic Permeability units" type="string" value="mD" />
             </ParameterList>      
             <Parameter name="Intrinsic Permeability Anisotropy" type="Array double" value="{1., 1., 0.1}" />
             <ParameterList name="Relative Permeability: Mualem">
               <Parameter name="Pore Interaction Term ???" type="Array double" value="{0.5, 0.5, 0.5}" />
               <Parameter name="Pore Interaction Term ??? Units" type="Array string" value="{}" />
             </ParameterList>
             <ParameterList name="Capillary Pressure: vG">
               <Parameter name="vG_alpha_n_Sr" type="Array double" value="{0.03 2.7 0.0234}" />
               <Parameter name="vG alpha units" type="Array string" value="{cm^-1, Null, Null}" />
             </ParameterList>
             <ParameterList name="Dispersivity: Uniform">
               <Parameter name="Dispersivity" type="double" value="0.5" />
               <Parameter name="Dispersivity units" type="string" value="cm"/>
             </ParameterList>    
             <Parameter name="Dispersivity Anisotropy" type="Array double" value="{10., 10., 1.0}" />
           </ParameterList>
       
         </ParameterList> <!-- End of soil specification -->
         
         <ParameterList name="Phase Definitions" >
           
           <!-- Definitions for Aqueous Phase -->
           <ParameterList name="Aqueous" >
       
             <!-- Definitions for Aqueous Phase Properties -->
             <ParameterList name="Phase Properties" >
               <ParameterList name="Viscosity: Uniform">
                 <Parameter name="Viscosity" type="double" value="0.01" />
                 <Parameter name="Viscosity Units" type="string" value="g cm^-1 s^-1" />
               </ParameterList>        
               <ParameterList name="Density: Uniform">
                 <Parameter name="Density" type="double" value="0.998" />
                 <Parameter name="Density Units" type="double" value="g cm^-3" />
               </ParameterList>        
             </ParameterList> <!-- End of Definitions for Aqueous Phase Properties -->
               
             <!-- Definitions for Aqueous Phase -->
             <ParameterList name="Phase Components" >
               
               <!-- Definitions for Aqueous Water + solutes -->
               <ParameterList name="Aqueous Water">        
       
                 <ParameterList name="Solute Properties">        
                   <ParameterList name="Sodium-nitrate">
                     <ParameterList name="Diffusion Coefficient: Uniform">
                       <Parameter name="Diffusion Coefficient" type="double" value="1.57e-9" />
                       <Parameter name="Diffusion Coefficient Units" type="string" value="m^2 s^-1" />
                     </ParameterList>
                   </ParameterList>
                   
                   <ParameterList name="Tc-99">
                     <ParameterList name="Diffusion Coefficient: Uniform">
                       <Parameter name="Diffusion Coefficient" type="double" value="2.e-9" />
                       <Parameter name="Diffusion Coefficient Units" type="string" value="m^2 s^-1" />
                     </ParameterList>
                   </ParameterList>
                 </ParameterList> <!-- End of Aqueous Phase Water Solute Properties definitions -->
                 
                 <!-- Initial Conditions for Aqueous Water + solutes -->
                 <ParameterList name="Component Initial Conditions">
                   <ParameterList name="All">
                     <!--
                        <ParameterList name="Initial Condition: Steady Richards">
                          <Parameter name="Water Table Height" type="double" value="0." />          
                          <Parameter name="Water Table Height Units" type="string" value="m" />          
                          <Parameter name="Water Pressure" type="double" value="0." />          
                          <Parameter name="Water Pressure Units" type="string" value="m" />          
                          <Parameter name="Water Pressure Location" type="double" value="0." />          
                          <Parameter name="Water Pressure Location Units" type="string" value="m" />          
                        </ParameterList>          
                        -->
                     <ParameterList name="Initial Condition: Linear">
                       <Parameter name="Saturation" type="double" value="0.5" />
                       <Parameter name="Reference Coordinate" type="Array double" value="{500, 1000, 97}" />
                       <Parameter name="Reference Coordinate Units" type="Array string" value="{m, m, m}"/>
                       <Parameter name="Gradient Value" type="Array double" value="{0, 0, -9793.5192}" />
                       <Parameter name="Gradient Value Units" type="Array double" value="{1/m, 1/m, 1/m}" />
                     </ParameterList>          
                     <ParameterList name="Solute Initial Conditions">        
                       <ParameterList name="Sodium-Nitrate">
                         <ParameterList name="Initial Condition: Uniform">
                           <Parameter name="Molar Concentration" type="double" value="0.0" />
                           <Parameter name="Molar Concentration Units" type="string" value="mol/L" />
                         </ParameterList>
                       </ParameterList>
                       <ParameterList name="Tc-99">
                         <ParameterList name="Initial Condition: Uniform">
                           <Parameter name="Molal Concentration" type="double" value="0.0" />
                           <Parameter name="Molal Concentration Units" type="string" value="mol/L_water" />
                         </ParameterList>
                       </ParameterList>
                     </ParameterList>
                   </ParameterList>
                 </ParameterList> <!-- End of Initial Conditions for Aqueous Water + solutes --> 
                 
                 
                 <!-- Boundary Conditions for Aqueous Water + solutes -->
                 <ParameterList name="Component Boundary Conditions">
                   
                   <!-- Note: BC must specify (1) Darcy flux or pressure, and saturation, -->
                   <!-- (2) Hydrostatic + water table, (3) impermeable, (4) outflow -->
                   <!-- ...plus, requires solute concentrations, unless outflow -->
                   <!-- Also: Time functions assume piecewise constant, unless specified linear -->
                   <ParameterList name="Top Region">
                     <ParameterList name="BC: Flux">
                       <Parameter name="Times" type="Array double" value="{0, 100, 1000}" />
                       <Parameter name="Time function" type="string" value="Piecewise Linear" />
                       <Parameter name="time units" type="Array string" value="{s, yr, yr}" />
                       <Parameter name="Integrated Volumetric Flux" type="Array double" value="{-500, -100, -50}" />
                       <Parameter name="Integrated Volumetric Flux Units" type="Array string" value="{mm/yr,  mm/yr,  mm/yr}" />
                       
                       <ParameterList name="Solute Boundary Conditions">
                         <ParameterList name="Sodium Nitrate">
                           <Parameter name="Times" type="Array double" value="{0.0, 50.0, 100.0}" />
                           <Parameter name="Times units" type="Array string" value="{s, yr, yr}" />
                           <Parameter name="Molar Concentration" type="Array double" value="{0.0, 5.0e-03, 1.0e-01}" />
                           <Parameter name="Molar Concentration Units" type="Array string" value="{mol/L, mol/L, mol/L}" />
                         </ParameterList>
                         <ParameterList name="Tc-99">
                           <Parameter name="Times" type="Array double" value="{0.0, 50.0, 100.0}" />
                           <Parameter name="Times units" type="Array string" value="{s, yr, yr}" />
                           <Parameter name="Molar Concentration" type="Array double" value="{0.0, 5.0e-03, 1.0e-01}" />
                           <Parameter name="Molar Concentration Units" type="Array string" value="{mol/L, mol/L, mol/L}" />
                         </ParameterList>
                       </ParameterList>
                     </ParameterList>
                   </ParameterList>
                   
                   <ParameterList name="Bottom Region">
                     <ParameterList name="BC: Dirichlet Pressure">
                       <Parameter name="Times" type="Array double" value="{0.0}" />
                       <Parameter name="Time units" type="Array string" value="{s}" />
                       <Parameter name="Value" type="Array double" value="{101325}" />
                       <Parameter name="Pressure Units" type="Array string" value="{Pa}" />
                     </ParameterList>
                     <ParameterList name="BC: Dirichlet Saturation">
                       <Parameter name="Times" type="Array double" value="{0.0}" />
                       <Parameter name="Time units" type="Array string" value="{s}" />
                       <Parameter name="Value" type="Array double" value="{1.}" />
                     </ParameterList>
                     <ParameterList name="Solute Boundary Conditions">
                       <ParameterList name="BC: Outflow">
                         <Parameter name="Times" type="Array double" value="{0.0}" />
                         <Parameter name="Time Units" type="Array string" value="s" />
                       </ParameterList>
                     </ParameterList>
                   </ParameterList>
                   
                 </ParameterList>    <!-- End of Boundary Conditions for Aqueous Water + solutes -->
                 
                 
                 <!-- Sources for Aqueous Water + solutes -->
                 <ParameterList name="Component Sources">
                   
                   <ParameterList name="Crib 1 Region">
                     <ParameterList name="Source: Volumetric Rate">
                       <Parameter name="Times" type="Array double" value="{0.0, 0.25}" />
                       <Parameter name="Times Units" type="Array string" value="{s, yr}" />
                       <Parameter name="Value" type="Array double" value="{0.35, 0.35}" />
                       <Parameter name="Value Units" type="Array string" value="{m^3/d,  m^3/d}" />
                     </ParameterList>
                     <ParameterList name="Solute Sources">
                       <ParameterList name="Sodium Nitrate">
                         <ParameterList name="Source: Mass Rate">
                           <Parameter name="Times" type="Array double" value="{0.0, 0.25}" />
                           <Parameter name="Time Units" type="Array string" value="{yr, yr}" />
                           <Parameter name="Value" type="Array double" value="{0.2, 0.2}" />
                           <Parameter name="Value Units" type="Array string" value="{kg/d,  kg/d}" />
                         </ParameterList>
                       </ParameterList>
                     </ParameterList>
                   </ParameterList>
                   
                   <ParameterList name="Crib 2 Region">
                     <ParameterList name="Source: Volumetric">
                       <Parameter name="Times" type="Array double" value="{0.0, 0.25}" />
                       <Parameter name="Times Units" type="Array string" value="{s, yr}" />
                       <Parameter name="Value" type="Array double" value="{0.35, 0.35}" />
                       <Parameter name="Value Units" type="Array string" value="{m^3/d, m^3/d}" />
                     </ParameterList>
                     
                     <ParameterList name="Solute Sources">
                       <ParameterList name="Sodium Nitrate">
                         <ParameterList name="Source: Mass Rate">
                           <Parameter name="Times" type="Array double" value="{0.0, 0.25}" />
                           <Parameter name="Times Units" type="Array string" value="{yr, yr}" />
                           <Parameter name="Value" type="Array double" value="{0.2, 0.2}" />
                           <Parameter name="Value Units" type="Array string" value="{kg/d, kg/d}" />
                         </ParameterList>
                       </ParameterList>
                     </ParameterList>
                   </ParameterList>
                   
                 </ParameterList> <!-- End of Sources for Aqueous Water + solutes -->
                 
               </ParameterList> <!-- End of Definitions for Aqueous Water + solutes -->
       
             </ParameterList> <!-- End of aqueous components definitions -->
             <!--
                <ParameterList name="...some other aqueous component">
                </ParameterList>
                -->
       
           </ParameterList> <!-- End of Definitions for Aqueous Phase -->
       
           <!--
              <ParameterList name="...some other phase">
              </ParameterList>
              -->
           
           <!-- NOTE: if the same component is multiple phases, requires specify mass transfer/phase eqm model -->
           
           <!-- May want to define profiles for immobile species in the soil matrix, though not sure if -->
           <!-- it should go here or in the soil definition
           
                <ParameterList name="Immobile Solutes">
                  <ParameterList name="Initial Conditions">
                    <ParameterList name="All">
                      <ParameterList name="Solutes">        
                        <ParameterList name="Quartz">
                          <ParameterList name="Initial Condition: Uniform">
                            <Parameter name="Molar Concentration" type="double" value="1.e-5" />
                            <Parameter name="Molar Concentration Units" type="string" value="mol/L" />
                          </ParameterList>
                        </ParameterList>
                        <ParameterList name="Kaolinite">
                          <ParameterList name="Initial Condition: Uniform">
                            <Parameter name="Molal Concentration" type="double" value="1.e-3" />
                            <Parameter name="Molal Concentration Units" type="string" value="mol/L_water" />
                          </ParameterList>
                        </ParameterList>
                      </ParameterList>
                    </ParameterList>
                  </ParameterList>
                </ParameterList>       
       
                -->    
         </ParameterList> <!-- End of Phase Definitions -->
         
         <!-- Definitions for Output -->
         <ParameterList name="Output">
                  
           <!-- Definitions for Observations -->
           <ParameterList name="Observations">
             <ParameterList name="Observation 1: Volume Integrals">        
               <Parameter name="Region" type="string" value="all"/>
               <ParameterList name="Times">
                 <Parameter name="Times" type="Array double"       value="{0.0, 0.5,  1.0,  10.0,  50.0, 100.0}" />
                 <Parameter name="Times Units" type="Array string" value="{  d,   d,    d,     d,     d,     d}" />
               </ParameterList>
               <Parameter name="Functional" type="string" value="Observation: Volume Integral"/>
               <Parameter name="Phase" type="string" value="Aqueous"/>
               <Parameter name="Include Solutes" type="bool" value="True"/>
             </ParameterList>
             
             <ParameterList name="Observation 2: Point Samples">
               <Parameter name="Region" type="string" value="Sample Point"/>
               <ParameterList name="Cycles">
                 <Parameter name="Cycle Frequency" type="integer" value="5" />
                 <Parameter name="Start Cycle" type="integer" value="15" />
                 <Parameter name="End Cycle" type="integer" value="150" />
               </ParameterList>
               <Parameter name="Functional" type="string" value="Observation: Point Sample"/>
               <Parameter name="Phase" type="string" value="Aqueous"/>
               <Parameter name="Include Solutes" type="bool" value="True"/>
             </ParameterList>
             
             <ParameterList name="Observation 3: Flux Integrals">
               <Parameter name="Region" type="string" value="Bottom Surface"/>
               <ParameterList name="Sample Times">
                 <Parameter name="Cycle Frequency" type="integer" value="5" />
                 <Parameter name="Cycle Start" type="integer" value="15" />
                 <Parameter name="Cycle End" type="integer" value="150" />
               </ParameterList>
               <Parameter name="Functional" type="string" value="Observation: Mass Flux Integral"/>
               <Parameter name="Phase" type="string" value="Aqueous"/>
               <Parameter name="Include Solutes" type="bool" value="False"/>
             </ParameterList>
           </ParameterList> <!-- End of Definitions for Observations -->
           
           <!-- Definitions for Checkpoints -->
           <ParameterList name="Checkpoint Data">
             <Parameter name="File Base Name" type="string" value="dump-" />
             <ParameterList name="Cycle Data">
               <Parameter name="Start" type="int" value="0"/>
               <Parameter name="End" type="int" value="500000"/>
               <Parameter name="Interval" type="int" value="5"/>
             </ParameterList>
           </ParameterList>
       
           <!-- Definitions for Visualization -->
           <ParameterList name="Visualization Data">
             <Parameter name="File Base Name" type="string" value="viz-" />
             <ParameterList name="Cycle Data">
               <Parameter name="Start" type="int" value="0"/>
               <Parameter name="End" type="int" value="500000"/>
               <Parameter name="Interval" type="int" value="5"/>
             </ParameterList>
           </ParameterList>
       
         </ParameterList>  <!-- End of Definitions for Output -->
       
       
         <ParameterList name="Execution control">
       
           <!--
              This section is highly specific to the numerical algorithm details, which
              will be a sensitive function of the mesh framework, the type of problem 
              selected, the mode requested for time integration, whether the mesh
              is dynamically adaptive, and a host of more detailed algorithm and model
              decisions.  
       
              The parameter set below represents a fictional calculation and depicts 
              an organization of the numerical parameters that might be appropriate.       
              The main ParameterList here is named after a labeled "type" of solve
              one might like to do.  Had this been an unsteady simulation, many of the
              linear and nonlinear solver parameters may not be applicable at all.
       
              It is unclear whether the inputs for this section can or should be orgainized
              at any finer a level of granularity.
       
              -->
           <ParameterList name="Amanzi Custom Steady Unstructured Flow Solver">
       
             <Parameter name="Simulation Start Time" type="double" value="0" />
             <Parameter name="Simulation End Time" type="double" value="3e7" />
             <Parameter name="Initial Delta-t" type="double" value="1.e-7" />
             
             <Parameter name="Max iterations" type="int" value="500" />
             <Parameter name="Relative or Absolute Tolerance" type="string" value="rel"/>
             <Parameter name="Tolerance" type="double" value="1.e-6"/>
             <Parameter name="Solver" type="string" value="Bi-CGSTAB" />
             <Parameter name="Preconditioner" type="string" value="ILU" />
             
             <ParameterList name="ODE Integrator">
               <Parameter name="Nonlinear Solver Max Iterations" type="int" value="5"/>
       	<Parameter name="Nonlinear Solver Tolerance" type="double" value="0.05"/>
               <Parameter name="NKA Max vectors" type="int" value="5"/>
       	<Parameter name="NKA Drop Tolerance" type="double" value="5.0e-2"/>
       	<Parameter name="Maximum Number of BDF Tries" type="int" value="50"/>
               <ParameterList name="Verbosit">
       	  <Parameter name="Value" type="string" value="Medium"/>
       	</ParameterList>
             </ParameterList>	
             
             <ParameterList name="Diffusion Preconditioner">
       	<ParameterList name="ML Parameters">
       	  <Parameter name="ML output" type="int" value="0"/>
         	  <Parameter name="max levels" type="int" value="40"/>
       	  <Parameter name="prec type" type="string" value="MGW"/>
       	  <Parameter name="cycle applications" type="int" value="5"/>
       	  <Parameter name="aggregation: type" type="string" value="Uncoupled"/> 
       	  <Parameter name="aggregation: damping factor" type="double" value="1.333"/>
       	  <Parameter name="aggregation: threshold" type="double" value="0.0"/>
       	  <Parameter name="aggregation: nodes per aggregate" type="int" value="3"/>
       	  <Parameter name="eigen-analysis: type" type="string" value="cg"/>
       	  <Parameter name="eigen-analysis: iterations" type="int" value="30"/>
       	  <Parameter name="smoother: sweeps" type="int" value="5"/>
       	  <Parameter name="smoother: damping factor" type="double" value="1.0"/>
       	  <Parameter name="smoother: pre or post" type="string" value="both"/>
       	  <Parameter name="smoother: type" type="string" value="Gauss-Seidel"/>
       	  <Parameter name="smoother: damping factor" type="double" value="1.0"/>
       	  <Parameter name="coarse: type" type="string" value="Amesos-KLU"/>
                 <Parameter name="coarse: max size" type="int" value="100"/>	   
                 <Parameter name="coarse: damping factor" type="double" value="1.0"/>	   
       	</ParameterList>
             </ParameterList>
           </ParameterList>	
       
         </ParameterList>	
       
       
           <!-- <ParameterList name="Structured"> -->
           <!-- Mesh-framework/model-ID - specific numerical control parameters -->
           <!--   <ParameterList name="Time Step control">
                      <Parameter name="cfl" type="string" value="0.8"/>
                      <Parameter name="init_shrink" type="string" value="1"/>
                      <Parameter name="Verbosity" type="string" value="3"/>
                  </ParameterList>
       
                  <ParameterList name="AMR">
                    <Parameter name="blocking_factor" type="string" value="16"/>
                    <Parameter name="derive_plot_vars" type="string" value="gradpx gradpy gradn"/>
                    <Parameter name="grid_eff" type="string" value="0.75"/>
                    <Parameter name="max_grid_size" type="string" value="64"/>
                    <Parameter name="max_level" type="string" value="0"/>
                    <Parameter name="n_cell" type="string" value="64 64"/>
                    <Parameter name="n_error_buf" type="string" value="2 2 2"/>
                    <Parameter name="plot_file" type="string" value="temp/plt"/>
                    <Parameter name="plot_int" type="string" value="1"/>
                    <Parameter name="probin_file" type="string" value="probin"/>
                    <Parameter name="ref_ratio" type="string" value="2 2"/>
                    <Parameter name="regrid_int" type="string" value="2"/>
                  </ParameterList>
       
                  <ParameterList name="Solver Controls">
                    <ParameterList name="diffuse">
                      <Parameter name="Verbosity" type="string" value="2"/>
                      <Parameter name="diff_abs_tol" type="string" value="1.e-14"/>
                      <Parameter name="diff_tol" type="string" value="1.e-12"/>
                    </ParameterList>
                    <ParameterList name="mac">
                      <Parameter name="Verbosity" type="string" value="3"/>
                      <Parameter name="mac_abs_tol" type="string" value="1.e-14"/>
                      <Parameter name="mac_sync_tol" type="string" value="1.e-12"/>
                      <Parameter name="mac_tol" type="string" value="1.e-12"/>
                    </ParameterList>
                    <ParameterList name="mg">
                      <Parameter name="cg_solver" type="string" value="1"/>
                      <Parameter name="maxiter" type="string" value="100"/>
                      <Parameter name="smooth_on_cg_unstable" type="string" value="1"/>
                      <Parameter name="Verbosity" type="string" value="1"/>
                    </ParameterList>
                    <ParameterList name="cg">
                      <Parameter name="unstable_criterion" type="string" value="100"/>
                      <Parameter name="use_jacobi_precond" type="string" value="1"/>
                      <Parameter name="Verbosity" type="string" value="0"/>
                    </ParameterList>
                  </ParameterList>           
         </ParameterList>  -->
           -->      
       
       
       
       </ParameterList>
