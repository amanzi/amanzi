========================================
Amanzi XML Input Specification
========================================

.. contents:: **Table of Contents**


Overview
========================================

The Amanzi simulator evolves a system of conservation
equations for reacting flow in porous media, as detailed in
the ASCEM report entitled "Mathematical Formulation Requirements and
Specifications for the Process Models`" (hereafter referred to
as the 'Model Requirements Document (MRD)'). The purpose of the present
document is to specify the data required to execute Amanzi.  This specification
should be regarded as a companion to the MRD; it relies heavily on
the detailed formulations of the models.  Where applicable, the
relevant sections of the MRD are indicated.


Each Amanzi simulation requires specification of a set of phase and
tracer components, and corresponding initial data, boundary conditions and source terms.  Conservation of mass for each of the
specified components is given by equation 2.11 of the MRD, where the
volumetric flow rate has been specified via Darcy's law (equation
2.10).  For Darcy flow, the properties of the (rock) medium must be identified
throughout the domain, including the phase component permeabilities,
etc. (see `"Rock`" section below for more details).  To complete the mathematical specification, any sources/sinks
for the evolved field quantities must be communicated.

A number of primitives are provided to support communicating the problem setup to Amanzi:

 * *Region*: Used to define the computation domain and labeled sub-regions for the purposes of communicating boundary and initial conditions, material properties, solution diagnostics, and source terms. 

 * *Rock*: Structure to characterize the physical properties of the porous media.

 * *State*: List of evolved quantities and their respective initial/boundary data and source terms.

Beyond the mathematical specification of the problem, the user must provide instructions to control the details of the numerical simulation.  Some of this information tends to be generic (`e.g.`, simulation stop time, type and amount of output, etc), while some is more specific to the numerical integration options.  The control parameter database is tree-like: options at a finer level of detail are dependendent on choices made at a higher level.

Additional notes:

 * Currently, the problem setup and control data is passed from the user into the Amanzi executable via a parameter list (specifically, a `Teuchos::ParameterList <http://trilinos.sandia.gov/packages/docs/r7.0/packages/teuchos/doc/html/index.html>`_). Each entry in a ParameterList object can itself be a ParameterList, or can be data.  Supported data types include double, float, short, int, bool, string), and simple arrays of these basic types.  ParameterList objects may be initialized using an XML file; examples of the proper syntax are included below.

 * It is intended that this specification be sufficiently detailed and flexible to support set up, solution and analysis of a broad range of simple model problems.  However, it is likely that details of the specification will evolve over time, as new capabilities are implemented within the Amanzi simulator.

 * The Amanzi code has a dual execution path, catering to the special requirements, and exploiting many of the unique advantages of structured versus unstructured mesh implementations, respectively.  Completeness of the implementation of the models discussed here will vary between mesh schemes, and will evolve over time as well.


Creating a Parameter List
--------------------------------------------

The Amanzi simulator input file is an XML file in ASCII text format, and must be framed at the beginning and end by the following statements:


.. code-block:: xml

  <ParameterList name="Main">

  </ParameterList>

The value in the "name" can be anything ("Main" in this example).  Individual parameters in a ParameterList are specified as follows:


.. code-block:: xml

  <ParameterList name="Sub">
    <Parameter name="CFL" type="double" value="0.9"/>
    <Parameter name="ratio" type="Array int" value="{2, 2, 4}"/>
  </ParameterList>

In this example, the sublist "Sub" has a parameter named "CFL" that is a "double" and has the value of 0.9, and a Teuchos::Array<int>
parameter named "ratio" such that ratio[0] = 2. ratio[1]=2 and ratio[2]=4.

It is vital to recognize that the input specification has a irreducibly circular dependence; definitions in one part of the specification
must be consistent with those of another.  To the extent possible, Amanzi will
recognize when this expected consistency is violated.  For example, "region" labels referenced in the State section are assumed to
have been defined in the "regions" section; if the entire input set is scanned and the referenced region label is not defined, an
error is thrown and the simulation is aborted.

In the remainder of this document, we attempt to adhere to the following standard for presentation.  Reserved keywords and labels are
`"quoted`" (and italicized) -- these labels or values of parameters in user-generated input files must match (using XML matching rules) the specified
or allowable values.  User-defined labels are indicated with ALL-CAPS, and are meant to represent a typical name given by a user -
these can be names or numbers or whatever serves best the organization of the user input data.

Where applicable, the relevant section the MRD is referred to by section or chapter number in parentheses.


0. Version
=======================================

Each input set contains at the top level a string variable `"Amanzi input format version`".  As of the most recent update of this specification, the
current version of the Amanzi input is `"0.9.1`".  If the version is unspecified, it is assumed to be earlier than `"0.9.0`".  The only difference between 
`"0.9.0`" and `"0.9.1`" is that the "grid_option" parameter was removed, and the mesh specification was moved from the "Regions" section and into 
a new "Mesh" section (section 1 below).  Options for `"grid_option`" parameter included `"Structured`" and `"Unstructured`".  In file version
`0.9.1`", a mesh framework is specified instead (see below).


1. Mesh
=======================================

The computational mesh is specified in this section, based on the `"Mesh Framework`", which can be `"Structured`" or a set of unstructured
options, including `"Simple`", `"stk:mesh`" (+...).  The `"Generate`" sublist of Mesh takes instructions that are specific to the framework - here 'generate' 
is a generic term for actual mesh generation (by Amanzi) or ingestion (file reads) to obtain mesh data created by pacakges external to Amanzi.

Notes:

 * A number of frameworks support the generation of logically rectangular, uniformly spaced structured meshes.  Under `"Generate`", all of these take a common set of instructions through three parameters: `"Number of Cells`" (integer array), `"Domain Low Corner`" (double array) and `"Domain High Corner`" (double array).

 * For the options that assume an external package generates the mesh, the data is passed into Amanzi through a file, and the `"Generate`" parameter list includes the name of that file `"filename`".  Additionally, as discussed in the "Regions" section below, mesh files produced by external packages may contain auxiliary data that associates a tag or label with each mesh entity (cells, faces, nodes).  These labeled sets can be assigned to a named region for use here. (see below).

2. State
=======================================

The state specifies the distribution of phases, species and pressure that are stored on the discrete mesh.  Generally, there
can be multiple phases (e.g. gaseous, aqueous, etc), and each is comprised of a number
of components (section 2.2).  Additionally, each component may carry a number of trace species.  
The tracers are assumed to have no impact on the thermodynamic and transport properties of the component, but may be involved in chemical
processes (Section 2.5).

Each state component must be labeled and defined in terms of physical properties: 
mass density, viscosity, and diffusivity (Section 4.6).  Boundary conditions must be specified along the entire surface
bounding the computational domain (Sections 3.3, 3.6, 3.10 and 4.3).  
Volumetric source terms, used to model infiltration (Section 3.7) and a wide variety of source and loss processes, are defined for each state, if applicable.

Tracers are labeled and defined in terms of 
their carrier component and group membership, as necessary to support the chemistry model (Section 5).  In particular, 
the "total concentration" (Equation 5.5) is the weighted sum of all tracers in the the tracer group "Total".
This is the only group of tracers that actually moves with the phase/component.  Other tracer groups include minerals
and surface complexation; they occupy a slot in the state but do not move with the flow (see Section 5).
Tracers may have volumetric sources as well, and like component sources their specification requires a region, strength and distribution functional.
Supported functionals for initial and boundary data and source distributions are listed below.

* "state" (list) can accept lists for named components (COMP).  Also a label specifies the dominant component

  * COMP (list) can accept values for the carrying phase name (string), mass density (double), viscosity (double) and diffusivity (double). IC is a named list to specify the instructions for constructing the intial state profile, and SOURCE (string) is a list to specify a set of volumetric sources.

    * IC (list) is named after a defined REGION, or the special denotation of `"default`".  `"default`" instructions will be used to fill the complement of the sum of the named regions.

      * IC-FUNC (list) can accept a set of parameter values for the functional (see table below for parameters required for each supported functional)

    * `"mass density`" (double) the mass density of this component

    * `"viscosity`" (double) the viscosity of this component

    * `"diffusivity`" (double) the diffusivity of this component

    * `"phase`" (string) the name of the phase that carries this component

    * `"source`" (list) can accept a REGION (string)

      * REGION (string) the name of a labeled region

        * S-FUNC-COMP (list) can accept a set of parameter values for the functional (see table below for parameters required for each supported functional)

    * `"tracer`" (list) can accept a name `"name`" (string), the list named `"source`", and a list for initial data (REGION-IC) and boundary data (REGIONBC)

      * `"name`" (string) the name of the tracer

      * REGION-IC or `"default`" (list) named after a defined region can accept a list with a label corresponding to a named IC functional (IC-FUNC)

        * IC-FUNC (list) can accept a set of parameter values for the functional (see table below for parameters required for each supported functional)

      * `"source`" (list) can accept the name of a distribution functional (S-FUNC-TRACER)

        * S-FUNC-TRACER (list) can accept a set of parameter values for the functional (see table below for parameters required for each supported functional)

      * REGION-BC (list) named after a region that defines a surface bounding the computational domain can accept a list (BC-FUNC) named after a boundary data function, BC_FUNC 
 
        * BC-FUNC (list) can accept a list (BC-PARAM) to specify the parameters of a named functional

  * `"dominant component`" (string) must be the name of one of the COMP lists defined above

Note: A label that identifies a region of non-zero volume is interpreted as an initial condition or source function instruction, whereas a region defined on the boundary of the 
computational domain will be interpreted as a boundary condition instruction.

Initial conditions are required for each component and tracer over the entire computational domain.
Boundary conditions are required on all domain boundaries (see Sections 3.3, 4.3).  Source terms for all are optional.  All are constructed using a limited number
of explicitly parameterized functional forms.  If the simulation is to be intialized using a restart file,
the phase component and tracer definitions are taken from the restart file, and initial condition instructions provided
here are quietly ignored, so that restarts are possible by simply changing a single control parameter (discussed in the control section).  Boundary conditions are
required regardless of the initial data, and must be defined consistently.

The following parameterized distribution functionals are supported for communicating initial conditions:
 * `"ic: constant`" requires `"value`" (see note below)
 * `"ic: coordinate-aligned linear`"  requires direction `"dir`" (string) of variation, `"x0_y0_slope`" (array double) specifying {x0, y0, m} for function
          of the form: `y-y0 = m*(x-x0)`.  Here y is state value, x is coordinate in `"dir`" direction.  For state values however, see note below.
 * `"ic: quadratic`" similar to linear
 * `"ic: exponential`" similar to linear

The following parameterized boundary conditions are supported for communicating boundary conditions:
 * `"bc: inflow`" requires `"bc: distribution`" (string) to set the distribution of the state upstream of the boundary (outside domain)
 * `"bc: outflow`"  requires `"bc: distribution`" (string) to set the distribution of the state downstream of the boundary (outside domain)
 * `"bc: seepage`" requires location `"water table height`" (double) of the water table.  If a more complex specification is needed, this should be changed to require a list to define it appropriately.
 * `"bc:  noflow`" requires no parameter data

The following models are currently supported for communicating source distribution:
 * `"source: uniform`" requires no parameters
 * `"source: linear`" requires location `"loc`" (array double) of a point or two locations, `"lo`", `"hi`" specifying a line or a rectangular plane
 * `"source: quadratic`" requires location `"loc`" (array double) of a point or two locations, `"lo`", `"hi`" specifying a line or a rectangular plane
 * `"source: exponential`" requires the exponent, `"exp`" and location `"loc`" (array double) of a point or two locations, `"lo`", `"hi`" specifying a line or a rectangular plane


NOTE:

Because of various physical constraints (e.g. component saturations sum to unity), initial and boundary functionals
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
          <Parameter name="x0_y0_slope" type="Array double" value="{.4, .9, 3}"/>
        </ParameterList>   
      </ParameterList>   
    </ParameterList> 
    <ParameterList name="water">
      <Parameter name="phase" type="string" value="aqueous"/>
      <Parameter name="density" type="double" value="1.e3"/>
      <Parameter name="viscosity" type="double" value="1.0"/>
      <Parameter name="diffusivity" type="double" value="0."/>
      <ParameterList name="source"/>
        <ParameterList name="middle"/>
          <ParameterList name="source: uniform"/>
            <Parameter name="strength" type="double" value="20."/>
          </ParameterList>
        </ParameterList>
      </ParameterList>   
      <ParameterList name="default"/>
        <ParameterList name="ic: uniform"/>
          <Parameter name="value" type="double" value=".8"/>
        </ParameterList>
      </ParameterList>   
      <ParameterList name="middle"/>
        <ParameterList name="ic: uniform"/>
          <Parameter name="strength" type="double" value="20."/>
        </ParameterList>
      </ParameterList>   
      <ParameterList name="xlobc">
        <ParameterList name="inflow">
          <ParameterList name="bc: constant">
            <Parameter name="value" type="double" value="1."/>
          </ParameterList> 
        </ParameterList> 
      </ParameterList> 
      <ParameterList name="xhibc">
        <ParameterList name="outflow">
          <ParameterList name="bc: constant">
            <Parameter name="value" type="double" value="1."/>
          </ParameterList> 
        </ParameterList> 
      </ParameterList> 
      <ParameterList name="tracer">
        <Parameter name="name" type="string" value="Uranium"/>
        <ParameterList name="all">
          <ParameterList name="ic: constant">
            <Parameter name="value" type="double" value=".004"/>
          </ParameterList> 
        <ParameterList name="xlobc">
          <ParameterList name="inflow">
            <ParameterList name="bc: constant">
              <Parameter name="value" type="double" value=".005"/>
            </ParameterList> 
          </ParameterList> 
        </ParameterList> 
        <ParameterList name="xhibc">
          <ParameterList name="outflow">
            <ParameterList name="bc: constant">
              <Parameter name="value" type="double" value="-1"/>
            </ParameterList> 
          </ParameterList> 
        </ParameterList> 
      </ParameterList> 
    </ParameterList> 
  </ParameterList> 

In this example, there are 2 phases (water, air).  Each phase consists of a single component.  Three
volumetric regions ("top", "middle" and "bottom"), and two boundary regions (xlobc and xhibc)
have been defined elsewhere.  The initial data for the fields are set using a combination of linear and
constant profile functions over the two volumetric regions.  The boundary conditions are Dirichlet inflow
on the low side and outflow on the high side.  There is a uniform source of water over the "middle" 
region with an integrated strength of 20.



3. Regions
=======================================

Regions are used in Amanzi to define subsets of the computational domain in order to specify the problem
to be solved, and the output desired.  Amanzi automatically defines the special region labeled `"all`", which is the 
entire simulation domain.  The user must additionally define the boundary surface(s) which enclose the domain.
Amanzi assumes that the union of the boundary surfaces envelopes the entire computational domain
(*i.e.* is "water-tight").  The special regions (`"all`" and the boundaries) may also serve as generic
regions (see the dicussion below for how these regions are labeled) and
can thus be used to specify other components of the problem (source terms, initial conditions, etc).


Special note:
For the `"structured`" mesh option, the bounding surfaces are implicitly defined as the planar surfaces that surround the domain,
and are automatically generated with the following labels `"xlobc`", `"xhibc`", `"ylobc`", `"yhibc`",
`"zlobc`", `"zhibc`" that are accessible throughout the input file.

For the `"unstructured`" mesh option, Amanzi supports fixed meshes in the MOAB and MSTK formats, as well as 
a simple mesh specification that accommodates a parallelepiped domain.  In the first two cases, the domain boundaries
must be identified explicitly in the mesh file, for `"simple mesh`", the boundaries are created automatically, as
following the scheme for the `"structured`" mesh option.  In most cases, these surfaces embedded in the mesh
files will be labeled subsets of the mesh faces.

Regions specifications take the following form

 * "regions" (list) can accept lists for named regions (REGION)

   * REGION (list) can accept lists (SHAPE) that specify a functional for its shape.

     * SHAPE (list) can accept lists of shape parameters (SHAPE-PARAMS) 

       * SHAPE-PARAMS (array double or string) parameters to specify shape

Amanzi supports parameterized forms for a number of analytic shapes, as well as more complex
definitions based on triangulated surface files: point, box, arbitrary, layer.  Depending on the functional, SHAPE requires
a number of parameters:

+------------------------+----------------------------------------+------------------------------+---------------------------------------------------------------------------------------------+
|  shape functional name | parameters                             | type(s)                      | Comment                                                                                     |
+========================+========================================+==============================+=============================================================================================+
| `"point"`              | `"loc`"                                | array double                 | Location of point in space                                                                  |
+------------------------+----------------------------------------+------------------------------+---------------------------------------------------------------------------------------------+
| `"box"`                | `"lo`", `"hi`"                         | array double, array double   | Location of boundary points of box                                                          |
+------------------------+----------------------------------------+------------------------------+---------------------------------------------------------------------------------------------+
| `"labeled set"`        | `"label`", `"file`", `"mesh framework`"| string, string, string       | Set per label defined in mesh file (see below)                                              |
+------------------------+----------------------------------------+------------------------------+---------------------------------------------------------------------------------------------+
| `"arbitrary"`          | `"file`"                               | string                       | Region enveloped by surface described in specified file (see note below for format of file) |
+------------------------+----------------------------------------+------------------------------+---------------------------------------------------------------------------------------------+
| `"layer"`              | `"file_lo`" `"file_hi`"                | string, string               | Region between surfaces described in specified files (see note below for format of file)    |
+------------------------+----------------------------------------+------------------------------+---------------------------------------------------------------------------------------------+
| `"surface"`            | `"id1`" `"name2`" ... `"idN`"          | string, string ,..., string  | Region between surfaces described in specified files (see note below for format of file)    |
+------------------------+----------------------------------------+------------------------------+---------------------------------------------------------------------------------------------+

Notes

* The "labeled set" region is defined by a label that was given to sets generated in a preprocessing step and stored in a mesh-dependent data file.  For example, an "exodus" type mesh file can be processed to tag cells and/or faces with specific labels, using a variety of external tools.  Regions based on such sets are assigned a user-defined label for Amanzi, which may or may not correspond to the original label in the exodus file.  The intersection of this volume and the computational domain defines the region.  Note that the file used to express this labeled set may be in any Amanzi-supported mesh framework (the mesh framework is specified in the parameters for this option).  This option is implemented as a special (piecewise-constant) case of a generalized interpolation operator.

* Surface files contain labeled triangulated surfaces in a format that is yet to be determined.  Regardless of the format, the user is responsible for ensuring that the intersections with other surfaces in the problem, including the boundaries, are `"exact`" (*i.e.* that surface intersections are `"watertight`" where applicable), and that the surfaces are contained within the computational domain.  If nodes defining surfaces are separated by a distance *s* < `"domain_epsilon`" Amanzi will consider them coincident; if they fall outside the domain, the elements they define are ignored.


Example:

.. code-block:: xml

  <ParameterList name="regions">
    <ParameterList name="all">
      <ParameterList name="box">
        <Parameter name="lo" type="Array double" value="{2, 3, 4}"/>
        <Parameter name="hi" type="Array double" value="{4, 5, 8}"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="top">
      <ParameterList name="box">
        <Parameter name="lo" type="Array double" value="{2, 3, 6}"/>
        <Parameter name="hi" type="Array double" value="{4, 5, 8}"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="middle">
      <ParameterList name="box">
        <Parameter name="lo" type="Array double" value="{2, 3, 6}"/>
        <Parameter name="hi" type="Array double" value="{4, 5, 8}"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="bottom">
      <ParameterList name="box">
        <Parameter name="lo" type="Array double" value="{2, 3, 4}"/>
        <Parameter name="lo" type="Array double" value="{4, 5, 6}"/>
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

Each rock type (Section 2.6) is given a label (string) and assigned a density (double) and models (string) for porosity, permeability and capillary pressure.  Each rock is assigned to
regions (string array), a list of regions.

* "rock" (list) can accept multiple lists for named rock types (ROCK)

  * ROCK (list) can accept lists to specify a model (MODEL) in a way that is yet to be determined, and a string array `"regions`" to specify where these properties apply.  Generally, the complete specification of rock properties should include models for porosity, relative permeability, capillary pressure and rock permeability.  However, there appears to be a motivation to specify using porosity, permeability and "water retention".  This needs to be sorted out.

    * MODEL (list) can accept model parameters (MODEL-PARAMS) 

      * MODEL-PARAMS (double, array double) parameters to specify model (see notes below for details of each model available)

    * `"regions`" (string array) a set of labels corresponding to defined regions


The following models are currently supported for porosity:
 * `"porosity: file`" requires the following strings: `"filename`" (name of a file), `"interpolation`" (the interpolation strategy: : `"constant`" or `"linear`"), `"framework`" (the mesh framework with which the file is compatible), and `"label`" (the label of the scalar field in the file to associate with the values of porosity).  In particular, the physical domain of this input data must completely cover the union of the regions over which this property is to be evaluated.
 * `"porosity: uniform`" requires a double specifying the constant value of porosity.
 * `"porosity: random`" requires the mean value of porosity and the percentage fluctuation, "porosity and fluctuation" (array double) to generate
 * `"porosity: gslib`" requires the name of a gslib-formatted file "gslib filename" to generate porosity field

The following models are currently supported for the absolute (rock) permeability:
 * `"permeability: file`" requires the following strings: `"filename`" (name of a file), `"interpolation`" (the interpolation strategy: `"constant`" or `"linear`"), `"framework`" (the mesh framework with which the file is compatible), and `"label`" (the label of the scalar field in the file to associate with the values of permeability).  The physical domain of this input data must completely cover the union of the regions over which this property is to be evaluated.
 * `"permeability: uniform`" requires a double specifying the constant value of porosity.
 * `"permeability: random`" requires the mean value of porosity and the percentage fluctuation, "mean permeability and rms fluctuation" (array double) to generate
 * `"permeability: gslib`" requires the name of a gslib-formatted file "gslib filename" to generate permeability field
 *  NOTE: All but `"permeability: file`" may also accept the array parameter `"permeability anisotropy`" (array double) to specify that the permeability is a diagonal tensor; these values are used to scale the X, Y and Z values.

The following models are currently supported for relative permeability (Section 2.6):
 * `"relative permeability: perfect`" requires no parameters, krl=krg=1
 * `"relative permeability: linear`" requires no parameters, krl=sl and krg=sg
 * `"relative permeability: quadratic`" requires slr, sgr (array double), krl=sc^2, krg=1-se^2, se=(sl-sg)/(1-slr-sgr)
 * `"relative permeability: vGM`" (van Genuchten-Mualem) requires m, slr, sgr (array double), krl=sqrt(se)(1-(1-se^-m)^m)^2, krg=(1-sekg)^1/3 (1-sekg^-m)^(2m), se=(sl-slr)/(1-slr-sgr), sekg=sl/(1/sgr)

The following models are currently supported for capillary pressure (Section 3.3.2):
 * `"capillary pressure: none`" requires no parameters, pc = 0
 * `"capillary pressure: linear`" requires no parameters, pc = sl
 * `"capillary pressure: vG`" requires m, sigma, slr, sgr (array double), pc=(1/sigma)(se^-m - 1)^-n, se=(sl-slr)/(1-slr-sgr)

The following models are currently supported for water retention (should we support this mode of specification?):
 * `"water retention: vG`" requires m, sigma, slr (array double)

Example:

.. code-block:: xml

  <ParameterList name="rock">
    <ParameterList name="backfill">
      <Parameter name="density" type="double" value="2.8e3"/>
      <ParameterList name="permeability: uniform">
        <Parameter name="permeability" type="double" value="1240"/>
        <Parameter name="permeability anisotropy" type="Array double" value="{1., 0.001, 0.001}"/>
      </ParameterList>
      <ParameterList name="porosity: uniform">
        <Parameter name="porosity" type="double" value="0.2585"/>
      </ParameterList>
      <ParameterList name="relative permeability: vGM">
        <Parameter name="m_slr_sgr" type="Array double" value="{0.6585, 0.0774, 0}"/>
      </ParameterList>
      <ParameterList name="capillary pressure: vG">
        <Parameter name="m_sigma_slr_sgr" type="Array double" value="{0.6585, 102.1, 0.0774, 0}"/>
      </ParameterList>
      <Parameter name="regions" type="string array" value="{top, bottom}"/>
    </ParameterList>
    <ParameterList name="fine sand">
      <Parameter name="density" type="double" value="2.8e3"/>
      <ParameterList name="permeability: uniform">
        <Parameter name="permeability" type="double" value="337.0"/>
      </ParameterList>
      <ParameterList name="porosity: uniform">
        <Parameter name="porosity" type="double" value="0.3586"/>
      </ParameterList>
      <ParameterList name="relative permeability: vGM">
        <Parameter name="m_slr_sgr" type="Array double" value="{0.4694, 0.0837, 0}"/>
      </ParameterList>
      <ParameterList name="capillary pressure: vG">
        <Parameter name="m_sigma_slr_sgr" type="Array double" value="{0.4694, 9.533, 0.0837, 0}"/>
      </ParameterList>
      <Parameter name="regions" type="string array" value="{middle}"/>
    </ParameterList>
  </ParameterList>

In this example, there are two types of rock, `"backfill`" (which fills bottom and top regions) and `"fine sand`" (which fills middle region).  Both have
van Genuchten models for relative permeability and capillary pressure.  The backfill has an anisotropic permeability, where the vertical value is 1000 times
the horizontal values.


5. Observation Data
=======================================

Observation data generally refers to simple diagnostic quantities extracted from a simulation for the purposes of characterizing
the response of the system to variations of input data.  Unlike very large datasets used for post-processing and simulation
restart, observation data for any particular simulations tends to consist of only a handful of scalar values.
Examples include volume and surface integrals, such as the total water mass in the system at a specific time.
Computation of observation data involves applying a parameterized 
functional on a specified state quantity or flux value at specific simulation times.

Each observation is given a label (string), state id (string), evaluation functional (list), region (string) and a list of times for evaluation.
The resulting observations are evaluated during the simulation and returned to the calling process

* "observation" (list) can accept multiple lists for named observations (OBSERVATION)

  * OBSERVATION (list) can accept values for `"state id`", `"region`", `"functional`" and `"times`"

    * `"region`" (string) a region defined above

    * `"state id`" (string) a state quantity defined above

    * `"functional`" (string) choses which funcitional to apply (see below)

    * `"times`" (array double) values of time where this quantity is desired

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
      <Parameter name="times" type="Array double" value="{1.e3, 2.e3, 2.5e3}"/>
    </ParameterList>
  </ParameterList>

In this example, the user requests the volume integral of the water density over the entire domain at three different times.
Amanzi will also support integrals and point samples of phase fluxes.  Note that times specified may not necessarily fall within
the time interval of the present simulation.  The format of the data structure used to communicate the observation data back
to the calling function includes a flag for each requested time to indicate whether the quantity was successfully filled.




6. Chemistry
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

6.1 Chemistry Thermodynamic data file formats 
-------------------------------------------------

6.1.1 XML Database format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

not yet implemented

6.1.2 Simple Database format
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


7. MPC
=======================================

Note that this section is highly specific to the unstructured mesh options, and to running a Richards model.  The parameter list will need to be completely revamped to more generally control the structured and unstructured options simultaneously.  However, for the present time, the following is an example of supported parameters in this list:

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



8. Transport
=======================================

Note that this section is highly specific to the unstructured mesh options, and to running a Richards model.  The parameter list will need to be completely revamped to more generally control the structured and unstructured options simultaneously.  However, for the present time, the following is an example of supported parameters in this list:


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


9. Flow
=======================================

Note that this section is highly specific to the unstructured mesh options, and to running a Richards model.  The parameter list will need to be completely revamped to more generally control the structured and unstructured options simultaneously. However, for the present time, the following is an example of supported parameters in this list:

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

Here's a partial list of additional parameters under the general category of "control".  Most of these are specific to the structured grid option.  This section will require a complete revamp.

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


