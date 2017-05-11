ATS Native XML Input Specification V1
#######################################

.. contents:: **Table of Contents**

  
Syntax of the Specification
#######################################

* Input specification for each ParameterList entry consists of two parts.  
  First, a bulleted list defines the usage syntax and available options.  
  This is followed by example snipets of XML code to demonstrate usage.

* In many cases, the input specifies data for a particular parameterized model, and ATS 
  supports a number of parameterizations.  
  For example, initial data might be uniform (the value is required), or linear in y (the value 
  and its gradient are required).  
  Where ATS supports a number of parameterized models for quantity Z, the available 
  models will be listed by name, and then will be described in the subsequent section.  
  For example, the specification for an `"X`" list might begin with the following:

  * `"Y`" ``[string]`` **"default_value"**, `"other`", `"valid`", `"options`"

  * Z ``[Z-spec]`` Model for Z, choose exactly one of the following: (1) `"z1`", or (2) `"z2`" (see below) 

Here, an `"X`" is defined by a `"Y`" and a `"Z`".  
The `"Y`" is a string parameter but the `"Z`" is given by a model (which will require its own set of parameters).
The options for `"Z`" will then be described as a spec:

 * `"z1`" applies model z1.  Requires `"z1a`" ``[string]``

 * `"z2`" applies model z2.  Requires `"z2a`" ``[double]`` and `"z2b`" ``[int]``

An example of using such a specification:

.. code-block:: xml

    <ParameterList name="X">
      <Parameter name="Y" type="string" value="hello"/>
      <ParameterList name="z2">
        <Parameter name="z2a" type="double" value="0.7"/>
        <Parameter name="z2b" type="int" value="3"/>
      </ParameterList>   
    </ParameterList>   
 
Here, the user is defining X with Y="hello", and Z will be a z2 constructed with z2a=0.7 and z2b=3.

Conventions:

* Reserved keywords and labels are `"quoted and italicized`" -- these
  labels or values of parameters in user-generated input files must
  match (using XML matching rules) the specified or allowable values.

* User-defined labels are indicated with ALL-CAPS, and are meant to
  represent a typical name given by a user - these can be names or
  numbers or whatever serves best the organization of the user input
  data.

* Bold values are default values, and are used if the Parameter
  is not provided.


Symbol Index
#############

:math:`|E|` | volume of a cell :math:`[m^X]` (where :math:`X` is the dimension of the mesh)
:math:`g` | gravitational acceleration vector :math:`[m s^-2]`
:math:`h` | ponded depth, or the water head over the surface :math:`[m]`
:math:`` | alternative, in context of the subsurface, water head :math:`[m]`
:math:`h_{snow}` | snow depth :math:`[m]`
:math:`K` | absolute permeability :math:`[m^2]`
:math:`k_r` | relative permeability :math:`[-]`
:math:`n_X` | molar density of phase X :math:`[mol m^-3]`
:math:`p` | pressure of the liquid phase :math:`[Pa]`
:math:`P_{s,r}` | precipitation of rain or snow, noting that snow is always a precipitation rate in snow-water-equivalent (SWE) basis.  :math:`[m s^-1]`
:math:`Q_w` | mass source of water :math:`[mol s^-1]`
:math:`s_X` | saturation of phase X :math:`[-]`
:math:`t` | time variable :math:`[s]`
:math:`z` | elevation :math:`[m]`
:math:`\nu` | dynamic viscosity of water :math:`[Pa s]`
:math:`\phi` | porosity of the soil :math:`[-]`
:math:`\rho` | mass density of a phase :math:`[kg m^-3]`
:math:`\Theta` | extensive water content of a cell :math:`[mol]`

   

  
Main
#######################################

The `"main`" ParameterList frames the entire input spec, and must contain
one sublist for each of the following sections.

* `"mesh`" ``[mesh-spec]``  See the Mesh_ spec.

* `"regions`" ``[list]``

  List of multiple Region_ specs, each in its own sublist named uniquely by the user.

* `"coordinator`" ``[coordinator-spec]``  See the Coordinator_ spec.

* `"visualization`" ``[visualization-spec]`` A Visualization_ spec for the main mesh/domain.

* `"visualization XX`" ``[visualization-spec]``

  Potentially more than one other Visualization_ specs, one for each domain `"XX`".  e.g. `"surface`"

* `"checkpoint`" ``[checkpoint-spec]`` A Checkpoint_ spec.

* `"observations`" ``[observation-spec]`` An Observation_ spec.

* `"PKs`" ``[list]``

  A list containing exactly one sublist, a PK_ spec with the top level PK.

* `"state`" ``[list]`` A State_ spec.

  
Mesh
#####


All processes are simulated on a domain, which is discretized through a mesh.

Multiple domains and therefore meshes can be used in a single simulation, and multiple meshes can be constructed on the fly.

The base mesh represents the primary domain of simulation.  Simple, structured
meshes may be generated on the fly, or complex unstructured meshes are
provided as ``Exodus II`` files.  The base *mesh* list includes either a
GeneratedMesh_,  MeshFromFile_, or LogicalMesh_ spec, as described below.

Additionally, a SurfaceMesh_ may be formed by lifting the surface of a
provided mesh and then flattening that mesh to a 2D surface.  ColumnMeshes_
which split a base mesh into vertical columns of cells for use in 1D models
may also be generated automatically.

Finally, mesh generation is hard and error-prone.  A mesh audit is provided,
which checks for many common geometric and topologic errors in mesh
generation.  This is reasonably fast, even for big meshes, and can be done through providing a "verify mesh" option.

* `"verify mesh`" ``[bool]`` **false** Perform a mesh audit.
* `"deformable mesh`" ``[bool]`` **false** Will this mesh be deformed?

GeneratedMesh
==============

Generated mesh are by definition structured, with uniform dx, dy, and dz.
Such a mesh is specified by a bounding box high and low coordinate, and a list
of number of cells in each direction.

* `"generate mesh`" ``[list]``

  * `"domain low coordinate`" ``[Array(double)]`` Location of low corner of domain
  * `"domain high coordinate`" ``[Array(double)]`` Location of high corner of domain
  * `"number of cells`" ``[Array(int)]`` the number of uniform cells in each coordinate direction

Example:

.. code-block:: xml

   <ParameterList name="mesh">
     <ParameterList name="generate mesh"/>
       <Parameter name="number of cells" type="Array(int)" value="{{100, 1, 100}}"/>
       <Parameter name="domain low coordinate" type="Array(double)" value="{{0.0, 0.0, 0.0}}" />
       <Parameter name="domain high coordinate" type="Array(double)" value="{{100.0, 1.0, 10.0}}" />
     </ParameterList>
   </ParameterList>   


MeshFromFile
==============

Meshes can be pre-generated in a multitude of ways, then written to "Exodus
II" file format, and loaded in ATS.

* `"read mesh file`" ``[list]`` accepts name, format of pre-generated mesh file

  * `"file`" ``[string]`` name of pre-generated mesh file. Note that in the case of an
        Exodus II mesh file, the suffix of the serial mesh file must be .exo and 
        the suffix of the parallel mesh file must be .par.
        When running in serial the code will read this the indicated file directly.
        When running in parallel and the suffix is .par, the code will instead read
        the partitioned files, that have been generated with a Nemesis tool and
        named as filename.par.N.r where N is the number of processors and r is the rank.
        When running in parallel and the suffix is .exo, the code will partition automatically
        the serial file.
     
  * `"format`" ``[string]`` format of pre-generated mesh file (`"MSTK`" or `"Exodus II`")

Example:

.. code-block:: xml

    <ParameterList name="mesh">
      <ParameterList name="read mesh file">
        <Parameter name="file" type="string" value="mesh_filename.exo"/>
        <Parameter name="format" type="string" value="Exodus II"/>
      </ParameterList>   
      <Parameter name="verify mesh" type="bool" value="true" />
    </ParameterList>


LogicalMesh
==============

** Document me! **


SurfaceMesh
==============

To lift a surface off of the mesh, a side-set specifying all surface faces
must be given.  These faces are lifted locally, so the partitioning of the
surface cells will be identical to the partitioning of the subsurface faces
that correspond to these cells.  All communication and ghost cells are set up.
The mesh is flattened, so all surface faces must have non-zero area when
projected in the z-direction.  No checks for holes are performed.  Surface
meshes may similarly be audited to make sure they are reasonable for
computation.

* `"surface sideset name`" ``[string]`` The Region_ name containing all surface faces.
* `"verify mesh`" ``[bool]`` **false** Verify validity of surface mesh.
* `"export mesh to file`" ``[string]`` Export the lifted surface mesh to this filename.

Example:

.. code-block:: xml

    <ParameterList name="mesh">
      <ParameterList name="read mesh file">
        <Parameter name="file" type="string" value="mesh_filename.exo"/>
        <Parameter name="format" type="string" value="Exodus II"/>
      </ParameterList>   
      <Parameter name="verify mesh" type="bool" value="true" />
      <ParameterList name="surface mesh">
        <Parameter  name="surface sideset name" type="string" value="surface_region"/>
        <Parameter name="verify mesh" type="bool" value="true" />
        <Parameter name="export mesh to file" type="string" value="surface_mesh.exo" />
      </ParameterList>   
    </ParameterList>


ColumnMeshes
==============

** Document me! **

Example:

.. code-block:: xml

    <ParameterList name="mesh">
      <ParameterList name="read mesh file">
        <Parameter name="file" type="string" value="mesh_filename.exo"/>
        <Parameter name="format" type="string" value="Exodus II"/>
      </ParameterList>   
      <ParameterList name="column meshes">
      </ParameterList>   
    </ParameterList>





Region
##########




Regions are geometrical constructs used in Amanzi to define subsets of
the computational domain in order to specify the problem to be solved, and the
output desired. Regions may represents zero-, one-, two- or three-dimensional
subsets of physical space.  for a three-dimensional problem, the simulation
domain will be a three-dimensional region bounded by a set of two-dimensional
regions.  If the simulation domain is N-dimensional, the boundary conditions
must be specified over a set of regions are (N-1)-dimensional.

Amanzi automatically defines the special region labeled *All*, which is the 
entire simulation domain. Currently, the unstructured framework does
not support the *All* region, but it is expected to do so in the
near future.

 * `"regions`" ``[list]`` can accept a number of uniquely named lists for regions

   * ``[region-spec]`` Geometric model primitive, as described below.

Amanzi supports parameterized forms for a number of analytic shapes, as well
as more complex definitions based on triangulated surface files.


**Notes:**

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

Example:

.. code-block:: xml

   <ParameterList>  <!-- parent list -->
     <ParameterList name="regions">
       <ParameterList name="TOP SECTION">
         <ParameterList name="region: box">
           <Parameter name="low coordinate" type="Array(double)" value="{2, 3, 5}"/>
           <Parameter name="high coordinate" type="Array(double)" value="{4, 5, 8}"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="MIDDLE SECTION">
         <ParameterList name="region: box">
           <Parameter name="low coordinate" type="Array(double)" value="{2, 3, 3}"/>
           <Parameter name="high coordinate" type="Array(double)" value="{4, 5, 5}"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="BOTTOM SECTION">
         <ParameterList name="region: box">
           <Parameter name="low coordinate" type="Array(double)" value="{2, 3, 0}"/>
           <Parameter name="high coordinate" type="Array(double)" value="{4, 5, 3}"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="INFLOW SURFACE">
         <ParameterList name="region: labeled set">
           <Parameter name="label"  type="string" value="sideset_2"/>
           <Parameter name="file"   type="string" value="F_area_mesh.exo"/>
           <Parameter name="format" type="string" value="Exodus II"/>
           <Parameter name="entity" type="string" value="face"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="OUTFLOW PLANE">
         <ParameterList name="region: plane">
           <Parameter name="point" type="Array(double)" value="{0.5, 0.5, 0.5}"/>
           <Parameter name="normal" type="Array(double)" value="{0, 0, 1}"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="BLOODY SAND">
         <ParameterList name="region: color function">
           <Parameter name="file" type="string" value="F_area_col.txt"/>
           <Parameter name="value" type="int" value="25"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="FLUX PLANE">
         <ParameterList name="region: polygon">
           <Parameter name="number of points" type="int" value="5"/>
           <Parameter name="points" type="Array(double)" value="{-0.5, -0.5, -0.5, 
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





Point
======

List *region: point* defines a point in space. 
This region consists of cells containing this point.

* `"coordinate`" ``[Array(double)]`` Location of point in space.

Example:

.. code-block:: xml

   <ParameterList name="DOWN_WIND150"> <!-- parent list defining the name -->
     <ParameterList name="region: point">
       <Parameter name="coordinate" type="Array(double)" value="{-150.0, 0.0, 0.0}"/>
     </ParameterList>
   </ParameterList>




Box
======


List *region: box* defines a region bounded by coordinate-aligned
planes. Boxes are allowed to be of zero thickness in only one
direction in which case they are equivalent to planes.

* `"low coordinate`" ``[Array(double)]`` Location of the boundary point with the lowest coordinates.

* `"high coordinate`" ``[Array(double)]`` Location of the boundary points with the highest coordinates.

Example:

.. code-block:: xml

   <ParameterList name="WELL">  <!-- parent list -->
     <ParameterList name="region: box">
       <Parameter name="low coordinate" type="Array(double)" value="{-5.0,-5.0, -5.0}"/>
       <Parameter name="high coordinate" type="Array(double)" value="{5.0, 5.0,  5.0}"/>
     </ParameterList>
   </ParameterList>
  



Plane
======

List *region: plane* defines a plane using a point lying on the plane and normal to the plane.

* `"normal`" ``[Array(double)]`` Normal to the plane.

* `"point`" ``[Array(double)]`` Point in space.

Example:

.. code-block:: xml

   <ParameterList name="TOP_SECTION"> <!-- parent list -->
     <ParameterList name="region: plane">
       <Parameter name="point" type="Array(double)" value="{2, 3, 5}"/>
       <Parameter name="normal" type="Array(double)" value="{1, 1, 0}"/>
       <ParameterList name="expert parameters">
         <Parameter name="tolerance" type="double" value="1.0e-05"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>




Labeled Set
============

The list *region: labeled set* defines a named set of mesh entities
existing in an input mesh file. This is the same file that contains
the computational mesh. The name of the entity set is given
by *label*.  For example, a mesh file in the Exodus II
format can be processed to tag cells, faces and/or nodes with
specific labels, using a variety of external tools. Regions based
on such sets are assigned a user-defined label for Amanzi, which may
or may not correspond to the original label in the exodus file.
Note that the file used to express this labeled set may be in any
Amanzi-supported mesh format (the mesh format is specified in the
parameters for this option).  The *entity* parameter may be
necessary to specify a unique set.  For example, an Exodus file
requires *cell*, *face* or *node* as well as a label (which is
an integer).  The resulting region will have the dimensionality 
associated with the entities in the indicated set. 

* `"label`" ``[string]`` Set per label defined in the mesh file.

* `"file`" ``[string]`` File name.

* `"format`" ``[string]`` Currently, we only support mesh files in the "Exodus II" format.

* `"entity`" ``[string]`` Type of the mesh object (cell, face, etc).

Example:

.. code-block:: xml

   <ParameterList name="AQUIFER">
     <ParameterList name="region: labeled set">
       <Parameter name="entity" type="string" value="cell"/>
       <Parameter name="file" type="string" value="porflow4_4.exo"/>
       <Parameter name="format" type="string" value="Exodus II"/>
       <Parameter name="label" type="string" value="1"/>
     </ParameterList>
   </ParameterList>




Color Function
===============


The list *region: color function* defines a region based a specified
integer color, *value*, in a structured color function file,
*file*. The format of the color function file is given below in
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

* `"file`" ``[string]`` File name.

* `"value`" ``[int]`` Color that defines the set in a tabulated function file.

Example:

.. code-block:: xml

   <ParameterList name="SOIL_TOP">
     <ParameterList name="region: color function">
       <Parameter name="file" type="string" value="geology_resamp_2D.tf3"/>
       <Parameter name="value" type="int" value="1"/>
     </ParameterList>
   </ParameterList>





Coordinator
############



In the `"coordinator`" sublist, the user specifies global control of
the simulation, including starting and ending times and restart options.  
 
* `"start time`" ``[double]``, **0.** Specifies the start of time in model time.
 
* `"start time units`" ``[string]``, **"s"**, `"d`", `"yr`"

* `"end time`" ``[double]`` Specifies the end of the simulation in model time.
 
* `"end time units`" ``[string]``, **"s"**, `"d`", `"yr`" 

* `"end cycle`" ``[int]`` If provided, specifies the end of the simulation in timestep cycles.

* `"restart from checkpoint file`" ``[string]`` If provided, specifies a path to the checkpoint file to continue a stopped simulation.

* `"wallclock duration [hrs]`" ``[double]`` After this time, the simulation will checkpoint and end.  Not required.

* `"required times`" ``[time-control-spec]``

  A TimeControl_ spec that sets a collection of times/cycles at which the simulation is guaranteed to hit exactly.  This is useful for situations such as where data is provided at a regular interval, and interpolation error related to that data is to be minimized.
   
Note: Either `"end cycle`" or `"end time`" are required, and if
both are present, the simulation will stop with whichever arrives
first.  An `"end cycle`" is commonly used to ensure that, in the case
of a time step crash, we do not continue on forever spewing output.

Example:

.. code-block::xml

   <!-- simulation control -->
   <ParameterList name="coordinator">
     <Parameter  name="end cycle" type="int" value="6000"/>
     <Parameter  name="start time" type="double" value="0."/>
     <Parameter  name="start time units" type="string" value="s"/>
     <Parameter  name="end time" type="double" value="1"/>
     <Parameter  name="end time units" type="string" value="yr"/>
     <ParameterList name="required times">
       <Parameter name="start period stop" type="Array(double)" value="{0,-1,86400}" />
     </ParameterList>
   </ParameterList>



   

Visualization
##############

A user may request periodic writes of field data for the purposes of
visualization in the `"visualization`" sublists.  ATS accepts a visualization
list for each domain/mesh, including surface and column meshes.  These are in
separate ParameterLists, entitled `"visualization`" for the main mesh, and
`"visualization surface`" on the surface mesh.  It is expected that, for any
addition meshes, each will have a domain name and therefore admit a spec of
the form: `"visualization DOMAIN-NAME`".



Each list contains all parameters as in a IOEvent_ spec, and also:

* `"file name base`" ``[string]`` **"visdump_data"**, **"visdump_surface_data"**
  
* `"dynamic mesh`" ``[bool]`` **false**

  Write mesh data for every visualization dump, this facilitates visualizing deforming meshes.


Example:

.. code-block:: xml

  <ParameterList name="visualization">
    <Parameter name="file name base" type="string" value="visdump_data"/>
  
    <Parameter name="cycles start period stop" type="Array(int)" value="{{0, 100, -1}}" />
    <Parameter name="cycles" type="Array(int)" value="{{999, 1001}}" />

    <Parameter name="times start period stop 0" type="Array(double)" value="{{0.0, 10.0, 100.0}}"/>
    <Parameter name="times start period stop 1" type="Array(double)" value="{{100.0, 25.0, -1.0}}"/>
    <Parameter name="times" type="Array(double)" value="{{101.0, 303.0, 422.0}}"/>

    <Parameter name="dynamic mesh" type="bool" value="false"/>
  </ParameterList>




  
Checkpoint
##############

A user may request periodic dumps of ATS Checkpoint Data in the
`"checkpoint`" sublist.  The user has no explicit control over the
content of these files, but has the guarantee that the ATS run will be
reproducible (with accuracies determined by machine round errors and
randomness due to execution in a parallel computing environment).
Therefore, output controls for Checkpoint Data are limited to file
name generation and writing frequency, by numerical cycle number.
Unlike `"visualization`", there is only one `"checkpoint`" list for
all domains/meshes.



Each list contains all parameters as in a IOEvent_ spec, and also:

* `"file name base`" ``[string]`` **"checkpoint"**

* `"file name digits`" ``[int]`` **5**

  Write mesh data for every visualization dump, this facilitates visualizing deforming meshes.

Example:

.. code-block:: xml

  <ParameterList name="checkpoint">
    <Parameter name="cycles start period stop" type="Array(int)" value="{{0, 100, -1}}" />
    <Parameter name="cycles" type="Array(int)" value="{{999, 1001}}" />
    <Parameter name="times start period stop 0" type="Array(double)" value="{{0.0, 10.0, 100.0}}"/>
    <Parameter name="times start period stop 1" type="Array(double)" value="{{100.0, 25.0, -1.0}}"/>
    <Parameter name="times" type="Array(double)" value="{{101.0, 303.0, 422.0}}"/>
  </ParameterList>

In this example, checkpoint files are written when the cycle number is
a multiple of 100, every 10 seconds for the first 100 seconds, and
every 25 seconds thereafter, along with times 101, 303, and 422.  Files will be written in the form: `"checkpoint00000.h5`".


  


 
Observation
##############


Observations are a localized-in-space but frequent in time view of
data, designed to get at useful diagnostic quantities such as
hydrographs, total water content, quantities at a point, etc.  These
are designed to allow frequent collection in time without saving huge
numbers of visualization files to do postprocessing.  In fact, these
should be though of as orthogonal data queries to visualization -- vis
is pointwise in time but complete in space, while observations are
pointwise/finite in space but complete in time.

A user may request any number of specific observations from ATS.  Each
observation spec involves a field quantity, a functional reduction
operator, a region from which it will extract its source data, and a
list of discrete times for its evaluation.  The observations are
evaluated during the simulation and written to disk.

* `"observations`" [list] can accept multiple ``[observation-spec]`` entries.

An ``[observation-spec]`` consists of the following quantities:

* `"observation output filename`" [string] user-defined name for the file that the observations are written to.

* `"variable`" [string] any ATS variable used by any PK, e.g. `"pressure`" or `"surface-water_content`"

* `"region`" [string] the label of a user-defined region

* `"location name`" [string] the mesh location of the thing to be measured, i.e. `"cell`", `"face`", or `"node`"

* `"functional`" [string] the label of a function to apply to the variable across the region.  Valid functionals include:
 * `"observation data: point`" returns the value of the field quantity at a point.  The region and location name must result in a single entity being selected.
 * `"observation data: extensive integral`" returns the sum of an (extensive) variable over the region.  This should be used for extensive quantities such as `"water_content`" or `"energy`".
 * `"observation data: intensive integral`" returns the volume-weighted average of an (intensive) variable over the region.  This should be used for intensive quantities such as `"temperature`" or `"saturation_liquid`".

* Additionally, each ``[observation-spec]`` contains all parameters as in a IOEvent_ spec, which are used to specify at which times/cycles the observation is collected.

For flux observations, and additional option is available:

* `"direction normalized flux`" [bool] *false* Normalize the flux to point in the outward-normal direction.  This is important when looking at fluxes across a boundary, for instance to plot a hydrograph.


Example:

.. code-block:: xml
  
  <ParameterList name="observations" type="ParameterList">
    <!-- This measures the hydrograph out the "east" face of the surface domain -->
    <ParameterList name="surface outlet flux" type="ParameterList">
      <Parameter name="variable" type="string" value="surface-mass_flux" />
      <Parameter name="direction normalized flux" type="bool" value="true" />
      <Parameter name="region" type="string" value="east" />
      <Parameter name="functional" type="string" value="observation data: extensive integral" />
      <Parameter name="delimiter" type="string" value=" " />
      <Parameter name="location name" type="string" value="face" />
      <Parameter name="observation output filename" type="string" value="surface_outlet_flux.dat" />
      <Parameter name="times start period stop" type="Array(double)" value="{0.0,86400.0,-1.0}" />
    </ParameterList>
    <!-- This measures the total water, in mols, in the entire subsurface domain -->
    <ParameterList name="subsurface water content" type="ParameterList">
      <Parameter name="variable" type="string" value="water_content" />
      <Parameter name="region" type="string" value="computational domain" />
      <Parameter name="functional" type="string" value="observation data: extensive integral" />
      <Parameter name="delimiter" type="string" value=" " />
      <Parameter name="location name" type="string" value="cell" />
      <Parameter name="observation output filename" type="string" value="water_content.dat" />
      <Parameter name="times start period stop" type="Array(double)" value="{0.0,86400.0,-1.0}" />
    </ParameterList>
    <!-- This tracks the temperature at a point -->
    <ParameterList name="temperature_probeA" type="ParameterList">
      <Parameter name="variable" type="string" value="temperature" />
      <Parameter name="region" type="string" value="probeA" />
      <Parameter name="functional" type="string" value="observation data: point" />
      <Parameter name="delimiter" type="string" value=" " />
      <Parameter name="location name" type="string" value="cell" />
      <Parameter name="observation output filename" type="string" value="temperature_probeA.dat" />
      <Parameter name="times start period stop" type="Array(double)" value="{0.0,86400.0,-1.0}" />
    </ParameterList>
  </ParameterList>





PK
#####



A process kernel represents a single or system of partial/ordinary
differential equation(s) or conservation law(s), and is used as the
fundamental unit for coupling strategies.

Implementations of this interface typically are either an MPC
(multi-process coupler) whose job is to heirarchically couple several
other PKs and represent the system of equations, or a Physical PK,
which represents a single equation.

All PKs have the following parameters in their spec:

* `"PK type`" ``[string]``

  The PK type is a special key-word which corresponds to a given class in the PK factory.  See available PK types listed below.

* `"PK name`" ``[string]`` **LIST-NAME**

  This is automatically written as the `"name`" attribute of the containing PK sublist, and need not be included by the user.

Example:

.. code-block:: xml

  <ParameterList name="PKs">
    <ParameterList name="my cool PK">
      <Parameter name="PK type" type="string" value="my cool PK"/>
       ...
    </ParameterList>
  </ParameterList>

.. code-block:: xml

  <ParameterList name="PKs">
    <ParameterList name="Top level MPC">
      <Parameter name="PK type" type="string" value="strong MPC"/>
       ...
    </ParameterList>
  </ParameterList>

 


Base PKs
===============

There are several types of PKs, and each PK has its own valid input spec.  However, there are three main types of PKs, from which nearly all PKs derive.  Note that none of these are true PKs and cannot stand alone.


PKPhysicalBase
----------------



``PKPhysicalBase`` is a base class providing some functionality for PKs which
are defined on a single mesh, and represent a single process model.  Typically
all leaves of the PK tree will inherit from ``PKPhysicalBase``.

* `"domain`" ``[string]`` **""**, e.g. `"surface`".

  Domains and meshes are 1-to-1, and the empty string refers to the main domain or mesh.  PKs defined on other domains must specify which domain/mesh they refer to.

* `"primary variable`" ``[string]``

  The primary variable associated with this PK, i.e. `"pressure`", `"temperature`", `"surface_pressure`", etc.

* `"initial condition`" ``[initial-condition-spec]``  See InitialConditions_.

  Additionally, the following parameters are supported:

 - `"initialize faces from cell`" ``[bool]`` **false**

   Indicates that the primary variable field has both CELL and FACE objects, and the FACE values are calculated as the average of the neighboring cells.


NOTE: ``PKPhysicalBase (v)-->`` PKDefaultBase_





PKBDFBase
----------------



``PKBDFBase`` is a base class from which PKs that want to use the ``BDF``
series of implicit time integrators must derive.  It specifies both the
``BDFFnBase`` interface and implements some basic functionality for ``BDF``
PKs.

* `"initial time step`" ``[double]`` **1.**

  The initial timestep size for the PK, this ensures that the initial timestep
  will not be **larger** than this value.

* `"assemble preconditioner`" ``[bool]`` **true** 

  A flag for the PK to not assemble its preconditioner if it is not needed by
  a controlling PK.  This is usually set by the MPC, not by the user.

In the top-most (in the PK tree) PK that is meant to be integrated implicitly,
several additional specs are included.  For instance, in a strongly coupled
flow and energy problem, these specs are included in the ``StrongMPC`` that
couples the flow and energy PKs, not to the flow or energy PK itself.
  
* `"time integrator`" ``[time-integrator-spec]`` is a TimeIntegrator_.

  Note that this is only provided in the top-most ``PKBDFBase`` in the tree --
  this is often a StrongMPC_ or a class deriving from StrongMPC_.

* `"preconditioner`" ``[preconditioner-spec]`` is a Preconditioner_.

  This spec describes how to form the (approximate) inverse of the preconditioner.
  
NOTE: ``PKBDFBase  (v)-->`` PKDefaultBase_





PKPhysicalBDFBase
-------------------



A base class for all PKs that are both physical, in the sense that they
implement an equation and are not couplers, and support the implicit
integration interface.  This largely just supplies a default error norm based
on error in conservation relative to the extent of the conserved quantity.

* `"absolute error tolerance`" [double] **1.0**

  Absolute tolerance, :math:`a_tol` in the equation below.

* `"relative error tolerance`" [double] **1.0**

  Relative tolerance, :math:`r_tol` in the equation below.

By default, the error norm used by solvers is given by:

:math:`ENORM(u, du) = |du| / ( a_tol + r_tol * |u| )`

The defaults here are typically good, or else good defaults are set in the
code, so these need not be supplied.


NOTE: ``PKPhysicalBDFBase -->`` PKBDFBase_
      ``PKPhysicalBDFBase -->`` PKPhysicalBase_
      ``PKPhysicalBDFBase (v)-->`` PKDefaultBase_




Physical PKs
===============

Physical PKs are the physical capability implemented within ATS.

Flow PKs
-----------

Richards PK
^^^^^^^^^^^^^^^


Solves Richards equation:

.. math::
  \frac{\partial \Theta}{\partial t} - \nabla \frac{k_r n_l}{\mu} K ( \nabla p + \rho g \cdot \hat{z} ) = Q_w


Options:

Variable naming:

* `"domain`" ``[string]`` **""**  Defaults to the base subsurface domain.

* `"primary variable`" ``[string]`` The primary variable associated with this PK, typically `"pressure`"


Other variable names, typically not set as the default is basically always good:

* `"conserved quantity suffix`" ``[string]`` **"water_content"**  If set, changes the conserved quantity key.

* `"conserved quantity key`" ``[string]`` **"DOMAIN-CONSERVED_QUANTITY_SUFFIX"** Typically not set, default is good. ``[mol]``

* `"mass density key`" ``[string]`` **"DOMAIN-mass_density_liquid"** liquid water density ``[kg m^-3]``

* `"molar density key`" ``[string]`` **"DOMAIN-molar_density_liquid"** liquid water density ``[mol m^-3]``

* `"permeability key`" ``[string]`` **"DOMAIN-permeability"** permeability of the soil medium ``[m^2]``

* `"conductivity key`" ``[string]`` **"DOMAIN-relative_permeability"** scalar coefficient of the permeability ``[-]``

* `"upwind conductivity key`" ``[string]`` **"DOMAIN-upwind_relative_permeability"** upwinded (face-based) scalar coefficient of the permeability.  Note the units of this are strange, but this represents :math:`\frac{n_l k_r}{\mu}`  ``[mol kg^-1 s^1 m^-2]``

* `"darcy flux key`" ``[string]`` **"DOMAIN-mass_flux"** mass flux across a face ``[mol s^-1]``

* `"darcy flux direction key`" ``[string]`` **"DOMAIN-mass_flux_direction"** direction of the darcy flux (used in upwinding :math:`k_r`) ``[??]``

* `"darcy velocity key`" ``[string]`` **"DOMAIN-darcy_velocity"** darcy velocity vector, interpolated from faces to cells ``[m s^-1]``

* `"darcy flux key`" ``[string]`` **"DOMAIN-mass_flux"** mass flux across a face ``[mol s^-1]``

* `"saturation key`" ``[string]`` **"DOMAIN-saturation_liquid"** volume fraction of the liquid phase ``[-]``


Time integration and timestep control:

* `"initial time step`" ``[double]`` **1.** Max initial time step size ``[s]``.

* `"time integrator`" ``[time-integrator-spec]`` is a TimeIntegrator_.

  Note that this is only provided if this Richards PK is not strongly coupled to other PKs.

* `"initial condition`" ``[initial-condition-spec]``  See InitialConditions_.

  Additionally, the following parameter is supported:

 - `"initialize faces from cell`" ``[bool]`` **false**

   Indicates that the primary variable field has both CELL and FACE objects,
   and the FACE values are calculated as the average of the neighboring cells.

Error control:

* `"absolute error tolerance`" [double] **DERIVED** Defaults to a porosity of 0.5 * a saturation of 0.1 * n_l.  A small, but significant, amount of water.

* `"relative error tolerance`" [double] **1** Take the error relative to the amount of water present in that cell.

* `"flux tolerance`" [double] **1**

  Multiplies the error in flux (on a face) relative to the min of water in the
  neighboring cells.  Typically only changed if infiltration is very small and
  the boundary condition is not converging, at which point it can be decreased
  by an order of magnitude at a time until the boundary condition is
  satisfied.

Boundary conditions:

* `"boundary conditions`" ``[subsurface-flow-bc-spec]`` **defaults to Neuman, 0 normal flux**

Physics control:

* `"permeability rescaling`" ``[double]`` **1** Typically 1e7 or order :math:`sqrt(K)` is about right.  This rescales things to stop from multiplying by small numbers (permeability) and then by large number (:math:`\rho / \mu`).



  
 


Permafrost Flow PK
^^^^^^^^^^^^^^^^^^^^

Overland Flow, head primary variable PK
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Overland Flow, pressure primary variable, PK
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Snow Distribution PK
^^^^^^^^^^^^^^^^^^^^


Energy PKs
-----------

Advection Diffusion PK
^^^^^^^^^^^^^^^^^^^^^^^

Energy Base PK
^^^^^^^^^^^^^^^^^^^^^^^

Two-Phase subsurface Energy PK
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Three-Phase subsurface Energy PK
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Three-Phase subsurface Energy PK
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Surface Ice Energy PK
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



Surface Energy Balance PKs
------------------------------


Surface Energy Balance / Snow -- Monolithic Version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



Surface Energy Balance -- Generic Version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



Biogeochemistry
-----------------


Biogeochemistry -- Monolithic Version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



Deformation
-------------


Volumetric Deformation
^^^^^^^^^^^^^^^^^^^^^^



MPCs
===============

MPCs couple other PKs, and are the non-leaf nodes in the PK tree.

WeakMPC
----------

StrongMPC
----------

Physical MPCs
===============
 coupling is an art, and requires special off-diagonal work.  Physical MPCs can derive from default MPCs to provide special work.

Coupled Water MPC
--------------------


Subsurface MPC
--------------------

Permafrost MPC
--------------------


State
##############

State consists of two sublists, one for evaluators and the other for
atomic constants.  The latter is currently called `"initial
conditions`", which is a terrible name which must be fixed.

example:

.. code-block:: xml
                
  <ParameterList name="state">
    <ParameterList name="field evaluators">
      ...
    </ParameterList>
    <ParameterList name="initial conditions">
      ...
    </ParameterList>
  </ParameterList>
 

Field Evaluators
=================

Many field evaluators exist, but most derive from one of four base types.

Field Evaluator Base Classes
-------------------------------

PrimaryVariableEvaluator
^^^^^^^^^^^^^^^^^^^^^^^^^^^

SecondaryVariableEvaluator
^^^^^^^^^^^^^^^^^^^^^^^^^^^

SecondaryVariablesEvaluator
^^^^^^^^^^^^^^^^^^^^^^^^^^^

IndependentVariableEvaluator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

While these provide base functionality, all of the physics are in the
following derived classes.

Water Content
-----------------

Water content is the conserved quantity in most flow equations, including
Richard's equation with and without ice.  A variety of evaluators are provided
for inclusion of multiple phases.

RichardsWaterContentEvaluator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Evaluator type: `"richards water content`"

Evaluates water content in cell E.

.. math::
  \Theta = \phi n_{{liq}} s_{{liq}} |E|

* `"my key`" ``[string]`` **DOMAIN_water_content** Set by code, not user. [mol]
* `"porosity key`" ``[string]`` **DOMAIN_porosity** Names the porosity variable. [-]
* `"saturation liquid key`" ``[string]`` **DOMAIN_saturation_liquid** Names the saturation variable. [-]
* `"molar density liquid key`" ``[string]`` **DOMAIN_molar_density_liquid** Names the density variable. [mol m^-3]
* `"cell volume key`" ``[string]`` **DOMAIN_cell_volume** Names the cell volume variable. [m^3]

Note that in the defaults, DOMAIN is determined from the name of the evaluated data, which is set by the name of the list.

Example:

.. code-block:: xml

  <ParameterList name="water_content">
    <Parameter name="evaluator type" type="string" value="richards water content"/>
  </ParameterList>



RichardsWaterContentWithVaporEvaluator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Evaluator type: `"richards water content with vapor`"

Evaluates water content in cell E.

.. math::
  \Theta = \phi (n_{{liq}} s_{{liq}} + n_{{gas}} s_{{gas}} \omega) |E|

* `"my key`" ``[string]`` **DOMAIN_water_content** Set by code, not user. [mol]
* `"porosity key`" ``[string]`` **DOMAIN_porosity** Names the porosity variable. [-]
* `"saturation liquid key`" ``[string]`` **DOMAIN_saturation_liquid** Names the saturation variable. [-]
* `"saturation gas key`" ``[string]`` **DOMAIN_saturation_gas** Names the gas saturation variable. [-]
* `"molar density liquid key`" ``[string]`` **DOMAIN_molar_density_liquid** Names the density variable. [mol m^-3]
* `"molar density gas key`" ``[string]`` **DOMAIN_molar_density_gas** Names the gas density variable. [mol m^-3]
* `"mol fraction vapor in gas key`" ``[string]`` **DOMAIN_mol_frac_gas** Names the molar fraction of water vapor in the gas phase variable. [-]
* `"cell volume key`" ``[string]`` **DOMAIN_cell_volume** Names the cell volume variable. [m^3]

Note that in the defaults, DOMAIN is determined from the name of the evaluated data, which is set by the name of the list.

Example:

.. code-block:: xml

  <ParameterList name="water_content">
    <Parameter name="evaluator type" type="string" value="richards water content with vapor"/>
  </ParameterList>



PermafrostWaterContentEvaluator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Evaluator type: `"permafrost water content`"

Evaluates water content in cell E.

.. math::
  \Theta = \phi (n_{{ice}} s_{{ice}} + n_{{liq}} s_{{liq}} + n_{{gas}} s_{{gas}} \omega) |E|

* `"my key`" ``[string]`` **DOMAIN_water_content** Set by code, not user. [mol]
* `"porosity key`" ``[string]`` **DOMAIN_porosity** Names the porosity variable. [-]
* `"saturation ice key`" ``[string]`` **DOMAIN_saturation_ice** Names the ice saturation variable. [-]
* `"saturation liquid key`" ``[string]`` **DOMAIN_saturation_liquid** Names the liquid saturation variable. [-]
* `"saturation gas key`" ``[string]`` **DOMAIN_saturation_gas** Names the gas saturation variable. [-]
* `"molar density ice key`" ``[string]`` **DOMAIN_molar_density_ice** Names the ice density variable. [mol m^-3]
* `"molar density liquid key`" ``[string]`` **DOMAIN_molar_density_liquid** Names the liquid density variable. [mol m^-3]
* `"molar density gas key`" ``[string]`` **DOMAIN_molar_density_gas** Names the gas density variable. [mol m^-3]
* `"mol fraction vapor in gas key`" ``[string]`` **DOMAIN_mol_frac_gas** Names the molar fraction of water vapor in the gas phase variable. [-]
* `"cell volume key`" ``[string]`` **DOMAIN_cell_volume** Names the cell volume variable. [m^3]

Note that in the defaults, DOMAIN is determined from the name of the evaluated data, which is set by the name of the list.

Example:

.. code-block:: xml

  <ParameterList name="water_content">
    <Parameter name="evaluator type" type="string" value="permafrost water content"/>
  </ParameterList>





Surface Water potential surfaces
---------------------------------

Evaluators for 

SurfaceElevation
^^^^^^^^^^^^^^^^^^

Evaluator type: `"meshed elevation`"

Evaluates the z-coordinate and the magnitude of the slope :math:``|\nambla_h z|``

* `"elevation key`" ``[string]`` **elevation** Name the elevation variable. [m]
* `"slope magnitude key`" ``[string]`` **slope_magnitude** Name the elevation variable. [-]
* `"dynamic mesh`" ``[bool]`` **false** Lets the evaluator know that the elevation changes in time, and adds the `"deformation`" dependency.

Example:

.. code-block:: xml

  <ParameterList name="elevation">
    <Parameter name="evaluator type" type="string" value="meshed elevation"/>
  </ParameterList>




SurfacePotential
^^^^^^^^^^^^^^^^^^^

Evaluator type: ""

.. math::
  h + z

* `"my key`" ``[string]`` **pres_elev** Names the surface water potential variable, h + z [m]
* `"height key`" ``[string]`` **ponded_depth** Names the height variable. [m]
* `"elevation key`" ``[string]`` **elevation** Names the elevation variable. [m]


NOTE: This is a legacy evaluator, and is not in the factory, so need not be in
the input spec.  However, we include it here because this could easily be
abstracted for new potential surfaces, kinematic wave, etc, at which point it
would need to be added to the factory and the input spec.

NOTE: This could easily be replaced by a generic AdditiveEvaluator_




SnowSurfacePotential
^^^^^^^^^^^^^^^^^^^^^^

Evaluator type: "snow skin potential"

.. math::
  h + z + h_{{snow}} + dt * P_{{snow}}

* `"my key`" ``[string]`` **snow_skin_potential** Names the potential variable evaluated [m]
* `"ponded depth key`" ``[string]`` **ponded_depth** Names the surface water depth variable. [m]
* `"snow depth key`" ``[string]`` **snow_depth** Names the snow depth variable. [m]
* `"precipitation snow key`" ``[string]`` **precipitation_snow** Names the snow precipitation key. [m]
* `"elevation key`" ``[string]`` **elevation** Names the elevation variable. [m]
* `"dt factor`" ``[double]`` A free-parameter factor for providing a time scale for diffusion of snow precipitation into low-lying areas.  Typically on the order of 1e4-1e7. This timestep times the wave speed of snow provides an approximate length of how far snow precip can travel.  Extremely tunable! [s]

NOTE: This is equivalent to a generic AdditiveEvaluator_

Example:

.. code-block:: xml

  <ParameterList name="snow_skin_potential" type="ParameterList">
    <Parameter name="field evaluator type" type="string" value="snow skin potential" />
    <Parameter name="dt factor" type="double" value="864000.0" />
  </ParameterList>





Generic Evaluators
---------------------------------

Several generic evaluators are provided.








InitialConditions
=================

Initial condition specs are used in two places -- in the PK_ spec
which describes the initial condition of primary variables, and in the
initial conditions sublist of state, in which the value of atomic
constants are provided.  In Amanzi, this list is also used for initial
conditions of primary variables are specified here, not within the PK
list (hence the name of this sublist).  In ATS, this sublist is pretty
much only used for constant scalars and constant vectors.

This list needs to be renamed -- it has nothing to do with inital conditions anymore.

Initialization of constant scalars
------------------------------------

A constant scalar field is the global (with respect to the mesh)
constant.  At the moment, the set of such fields includes atmospheric
pressure.  The initialization requires to provide a named sublist with
a single parameter `"value`".

.. code-block:: xml

  <ParameterList name="fluid_density">
    <Parameter name="value" type="double" value="998.0"/>
  </ParameterList>


Initialization of constant vectors
------------------------------------

A constant vector field is the global (with respect to the mesh)
vector constant.  At the moment, the set of such vector constants
includes gravity.  The initialization requires to provide a named
sublist with a single parameter `"Array(double)`". In two dimensions,
is looks like

.. code-block:: xml

  <ParameterList name="gravity">
    <Parameter name="value" type="Array(double)" value="{0.0, -9.81}"/>
  </ParameterList>


Initialization of scalar fields
------------------------------------

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
------------------------------------
 
A variable tensor (or vector) field is defined similarly to a variable
scalar field.  The difference lies in the definition of the function
which is now a multi-values function.  The required parameters are
`"Number of DoFs`" and `"Function type`".

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


Initialization from a file
------------------------------------

Some data can be initialized from files. Additional sublist has to be
added to named sublist of the `"state`" list with the file name and
the name of attribute.  For a serial run, the file extension must be
`".exo`".  For a parallel run, it must be `".par`".  Here is an
example:

.. code-block:: xml

  <ParameterList name="permeability">
    <ParameterList name="exodus file initialization">
      <Parameter name="file" type="string" value="mesh_with_data.exo"/>
      <Parameter name="attribute" type="string" value="perm"/>
    </ParameterList>
  </ParameterList>



example:

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

    </ParameterList>
  </ParameterList>




Time integrators, solvers, and other mathematical specs
####################################################################################

Common specs for all solvers and time integrators, used in PKs.


TimeIntegrator
=================

Linear Solver Spec
===================

For each solver, a few parameters are used:

* `"iterative method`" ``[string]`` `"pcg`", `"gmres`", or `"nka`"

  defines which method to use.

* `"error tolerance`" ``[double]`` **1.e-6** is used in the convergence test.

* `"maximum number of iterations`" ``[int]`` **100** is used in the convergence test.

* `"convergence criteria`" ``[Array(string)]``  **{"relative rhs"}** specifies multiple convergence criteria. The list
  may include `"relative residual`", `"relative rhs`", and `"absolute residual`", and `"???? force once????`"

* `"size of Krylov space`" ``[int]`` is used in GMRES iterative method. The default value is 10.

.. code-block:: xml

     <ParameterList name="my solver">
       <Parameter name="iterative method" type="string" value="gmres"/>
       <Parameter name="error tolerance" type="double" value="1e-12"/>
       <Parameter name="maximum number of iterations" type="int" value="400"/>
       <Parameter name="convergence criteria" type="Array(string)" value="{relative residual}"/>
       <Parameter name="size of Krylov space" type="int" value="10"/>

       <ParameterList name="VerboseObject">
         <Parameter name="Verbosity Level" type="string" value="high"/>
       </ParameterList>
     </ParameterList>


Preconditioner
===================

These can be used by a process kernel lists to define a preconditioner.  The only common parameter required by all lists is the type:

 * `"preconditioner type`" ``[string]`` **"identity"**, `"boomer amg`", `"trilinos ml`", `"block ilu`" ???
 * `"PC TYPE parameters`" ``[list]`` includes a list of parameters specific to the type of PC.

Example:

.. code-block:: xml

     <ParameterList name="my preconditioner">
       <Parameter name="type" type="string" value="trilinos ml"/>
        <ParameterList name="trilinos ml parameters"> ?????? check me!
            ... 
        </ParameterList>
     </ParameterList>


Hypre's Boomer AMG
-------------------

Internal parameters for Boomer AMG include

* `"tolerance`" ``[double]`` if is not zero, the preconditioner is dynamic 
  and approximate the inverse matrix with the prescribed tolerance (in
  the energy norm ???).

* `"smoother sweeps`" ``[int]`` **3** defines the number of smoothing loops. Default is 3.

* `"cycle applications`" ``[int]`` **5** defines the number of V-cycles.

* `"strong threshold`" ``[double]`` **0.5** defines the number of V-cycles. Default is 5.

* `"relaxation type`" ``[int]`` **6** defines the smoother to be used. Default is 6 
  which specifies a symmetric hybrid Gauss-Seidel / Jacobi hybrid method. TODO: add others!

* `"coarsen type`" ``[int]`` **0** defines the coarsening strategy to be used. Default is 0 
  which specifies a Falgout method. TODO: add others!

* `"max multigrid levels`" ``[int]`` optionally defined the maximum number of multigrid levels.

* `"number of functions`" ``[int]`` **1**  Any value > 1 tells Boomer AMG to use the `"systems 
  of PDEs`" code.  Note that, to use this approach, unknowns must be ordered with 
  DoF fastest varying (i.e. not the native Epetra_MultiVector order).  By default, it
  uses the `"unknown`" approach in which each equation is coarsened and
  interpolated independently.  **Getting this correct is very helpful!**
  
  * `"nodal strength of connection norm`" ``[int]`` tells AMG to coarsen such
    that each variable has the same coarse grid - sometimes this is more
    "physical" for a particular problem. The value chosen here for nodal
    determines how strength of connection is determined between the
    coupled system.  I suggest setting nodal = 1, which uses a Frobenius
    norm.  This does NOT tell AMG to use nodal relaxation.
    Default is 0.

* `"verbosity`" ``[int]`` **0** prints a summary of run time settings and
  timing information to stdout.  `"1`" prints coarsening info, `"2`" prints
  smoothing info, and `"3`'" prints both.

Example:
  
.. code-block:: xml

  <ParameterList name="boomer amg parameters">
    <Parameter name="tolerance" type="double" value="0.0"/>
    <Parameter name="smoother sweeps" type="int" value="3"/>
    <Parameter name="cycle applications" type="int" value="5"/>
    <Parameter name="strong threshold" type="double" value="0.5"/>
    <Parameter name="coarsen type" type="int" value="0"/>
    <Parameter name="relaxation type" type="int" value="3"/>
    <Parameter name="verbosity" type="int" value="0"/>
    <Parameter name="number of functions" type="int" value="1"/>
  </ParameterList>




Trilinos ML
-------------------

Internal parameters of Trilinos ML includes

Example:

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
-------------------


The internal parameters for block ILU are as follows:

Example:

.. code-block:: xml

  <ParameterList name="block ilu parameters">
    <Parameter name="fact: relax value" type="double" value="1.0"/>
    <Parameter name="fact: absolute threshold" type="double" value="0.0"/>
    <Parameter name="fact: relative threshold" type="double" value="1.0"/>
    <Parameter name="fact: level-of-fill" type="int" value="0"/>
    <Parameter name="overlap" type="int" value="0"/>
    <Parameter name="schwarz: combine mode" type="string" value="Add"/>
    </ParameterList>
  </ParameterList>




Indentity
-------------------
The default, no PC applied.



NonlinearSolver
===================




Other Common Specs
##########################################

IOEvent
===================



The IOEvent is used for multiple objects that need to indicate simulation times or cycles on which to do something.

* `"cycles start period stop`" ``[Array(int)]`` 

    The first entry is the start cycle, the second is the cycle
    period, and the third is the stop cycle or -1, in which case there
    is no stop cycle. A visualization dump is written at such
    cycles that satisfy cycle = start + n*period, for n=0,1,2,... and
    cycle < stop if stop != -1.0.

* `"cycles start period stop N`" ``[Array(int)]`` 

    If multiple cycles start period stop parameters are needed, then
    use these parameters with N=0,1,2,...

* `"cycles`" ``[Array(int)]`` 
  
    An array of discrete cycles that at which a visualization dump is
    written.

* `"times start period stop`" ``[Array(double)]`` 

    The first entry is the start time, the second is the time period,
    and the third is the stop time or -1, in which case there is no
    stop time. A visualization dump is written at such times that
    satisfy time = start + n*period, for n=0,1,2,... and time < stop
    if stop != -1.0.  Note that all times units are in seconds.

* `"times start period stop n`" ``[Array(double)]``

    If multiple start period stop parameters are needed, then use this
    these parameters with n=0,1,2,..., and not the single `"times
    start period stop`" parameter.  Note that all times units are in
    seconds.

* `"times`" ``[Array(double)]`` 

    An array of discrete times that at which a visualization dump
    shall be written.  Note that all times units are in seconds.
 


VerboseObject
===================



This allows control of log-file verbosity for a wide variety of objects
and physics.

* `"verbosity level`" ``[string]`` **GLOBAL_VERBOSITY**, `"low`", `"medium`", `"high`", `"extreme`"

   The default is set by the global verbosity spec, (fix me!)  Typically,
   `"low`" prints out minimal information, `"medium`" prints out errors and
   overall high level information, `"high`" prints out basic debugging, and
   `"extreme`" prints out local debugging information.

Note: while there are other options, users should typically not need them.
Instead, developers can use them to control output.
   
Example:

.. code-block:: xml

  <ParameterList name="verbose object">
    <Parameter name="verbosity level" type="string" value="medium"/>
    <Parameter name="name" type="string" value="my header"/>
    <Parameter name="hide line prefix" type="bool" value="false"/>
    <Parameter name="write on rank" type="int" value="0"/>
  </ParameterList>



   

Function
===================


Analytic, algabraic functions of space and time are used for a variety of
purposes, including boundary conditions, initial conditions, and independent
variables.

For initial conditions, functions are prescribed of space only, i.e.

:math:`u = f(x,y,z)`

For boundary conditions and independent variables, functions are also a
function of time:

:math:`u = f(t,x,y,z)`

A ``[function-spec]`` is used to prescribe these functions.




It is straightforward to add new functions as needed.

Constant Function
-------------------------


Constant function is defined as :math:`f(x) = a`, for all :math:`x`. 

* `"value`" ``[double]`` The constant to be applied.

Example:

.. code-block:: xml

  <ParameterList name="function-constant">
    <Parameter name="value" type="double" value="1.0"/>
  </ParameterList>


  

Tabular Function
-------------------------


A piecewise function of one variable.

A tabular function is tabulated on a series of intervals; given values
:math:`{{x_i}}, {{y_i}},, i=0, ... n-1` and functional forms :math:`{{f_j}},,
j=0, ... n-2` a tabular function :math:`f(x)` is defined as:

.. math::
  \begin{matrix}
  f(x) &=& y_0, & x \le x_0,\\
  f(x) &=& f_{{i-1}}(x)  & x \in (x_{{i-1}}, x_i],\\
  f(x) &=& y_{{n-1}}, & x > x_{{n-1}}.
  \end{matrix}

The functional forms :math:`{f_j}` may be constant, which uses the left endpoint, i.e.

:math:`f_i(x) = y_i`,

linear, i.e.

:math:`f_i(x) = ( y_i * (x - x_i) + y_{{i+1}} * (x_{{i+1}} - x) ) / (x_{{i+1}} - x_i)`

or arbitrary, in which the :math:`f_j` must be provided.

The :math:`x_i` and :math:`y_i` may be provided in one of two ways -- explicitly in the input spec or from an HDF5 file.  The length of these must be equal, and the :math:`x_i` must be monotonically increasing.  Forms, as defined on intervals, must be of length equal to the length of the :math:`x_i` less one.

Explicitly specifying the data:

* `"x values`" ``[Array(double)]`` the :math:`x_i`
* `"y values`" ``[Array(double)]`` the :math:`y_i`
* `"forms`" ``[Array(string)]`` **linear**, `"constant`", `"USER_DEFINED`"
* `"USER_DEFINED`" ``[function-spec]`` user-provided functional forms on the interval
* `"x coordinate`" ``[string]`` **t**, `"x`", `"y`", `"z`" defines which coordinate direction the :math:`x_i` are formed, defaulting to time.

The below example defines a function that is zero on interval :math:`(-\infty,\,0]`,
linear on interval :math:`(0,\,1]`, constant (`f(x)=1`) on interval :math:`(1,\,2]`, 
square root of `t` on interval :math:`(2,\,3]`,
and constant (`f(x)=2`) on interval :math:`(3,\,\infty]`.

Example:

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
  

Loading table from file (note that `"USER_DEFINED`" is not an option here, but could be made so if requested):


* `"file`" ``[string]`` filename of the HDF5 data
* `"x header`" ``[string]`` name of the dataset for the :math:`x_i` in the file
* `"y header`" ``[string]`` name of the dataset for the :math:`y_i` in the file
* `"forms`" ``[Array(string)]`` **linear**, `"constant`"

The example below would perform linear-interpolation on the intervals provided by data within the hdf5 file `"my_data.h5`".

Example:

.. code-block:: xml
  
  <ParameterList name="function-tabular">
    <Parameter name="file" type="string" value="my_data.h5"/>
    <Parameter name="x coordinate" type="string" value="t"/>
    <Parameter name="x header" type="string" value="/time"/>
    <Parameter name="y header" type="string" value="/data"/>
  </ParameterList>




Smooth step Function
-------------------------


A smooth :math:`C^2` function `f(x)` on interval :math:`[x_0,\,x_1]` is
defined such that `f(x) = y_0` for `x < x0`, `f(x) = y_1` for `x > x_1`, and
monotonically increasing for :math:`x \in [x_0, x_1]` through cubic
interpolation.

Example:

.. code-block:: xml

  <ParameterList name="function-smooth-step">
    <Parameter name="x0" type="double" value="0.0"/>
    <Parameter name="y0" type="double" value="0.0"/>
    <Parameter name="x1" type="double" value="1.0"/>
    <Parameter name="y1" type="double" value="2.0"/>
  </ParameterList>




Polynomial Function
-------------------------


A generic polynomial function is given by the following expression:

.. math::
  f(x) = \sum_{{j=0}}^n c_j (x - x_0)^{{p_j}}

where :math:`c_j` are coefficients of monomials,
:math:`p_j` are integer exponents, and :math:`x_0` is the reference point.

Example:

.. code-block:: xml

  <ParameterList name="function-polynomial">
    <Parameter name="coefficients" type="Array(double)" value="{{1.0, 1.0}}"/>
    <Parameter name="exponents" type="Array(int)" value="{{2, 4}}"/>
    <Parameter name="reference point" type="double" value="0.0"/>
  </ParameterList>


  

Multi-variable linear Function
------------------------------


A multi-variable linear function is formally defined by
 
.. math::
  f(x) = y_0 + \sum_{{j=0}}^{{n-1}} g_j (x_j - x_{{0,j}}) 

with the constant term "math:`y_0` and gradient :math:`g_0,\, g_1\,..., g_{{n-1}}`.
If the reference point :math:`x_0` is specified, it must have the same
number of values as the gradient.  Otherwise, it defaults to zero.
Note that one of the parameters in a multi-valued linear function can be time.
Here is an example:

.. code-block:: xml

  <ParameterList name="function-linear">
    <Parameter name="y0" type="double" value="1.0"/>
    <Parameter name="gradient" type="Array(double)" value="{{1.0, 2.0, 3.0}}"/>
    <Parameter name="x0" type="Array(double)" value="{{2.0, 3.0, 1.0}}"/>
  </ParameterList>


  

Separable Function
------------------


A separable function is defined as the product of other functions such as

.. math::
  f(x_0, x_1,...,x_{{n-1}}) = f_1(x_0)\, f_2(x_1,...,x_{{n-1}})

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




Additive Function
------------------


An additive function simply adds two other function results together.

.. math::
  f(x) = f_1(x) + f_2(x)

where :math:`f_1` is defined by the `"function1`" sublist, and 
:math:`f_2` by the `"function2`" sublist:

.. code-block:: xml

  <ParameterList name="function-additive">
    <ParameterList name="function1">
      function-specification
    </ParameterList>
    <ParameterList name="function2">
      function-specification
    </ParameterList>
  </ParameterList>



Multiplicative Function
--------------------------


A multiplicative function simply multiplies two other function results together.

.. math::
  f(x) = f_1(x) * f_2(x)

where :math:`f_1` is defined by the `"function1`" sublist, and 
:math:`f_2` by the `"function2`" sublist:

.. code-block:: xml

  <ParameterList name="function-multiplicative">
    <ParameterList name="function1">
      function-specification
    </ParameterList>
    <ParameterList name="function2">
      function-specification
    </ParameterList>
  </ParameterList>



Composition Function
--------------------------


Function composition simply applies one function to the result of another.

.. math::
  f(x) = f_1( f_2(x) )

where :math:`f_1` is defined by the `"function1`" sublist, and 
:math:`f_2` by the `"function2`" sublist:

.. code-block:: xml

  <ParameterList name="function-composition">
    <ParameterList name="function1">
      function-specification
    </ParameterList>
    <ParameterList name="function2">
      function-specification
    </ParameterList>
  </ParameterList>



Piecewise Bilinear Function
---------------------------


A piecewise bilinear function extends the linear form of the tabular function to two variables.

Define :math:`i(x) = i : x_i < x <= x_{{i+1}}` and similarly :math:`j(y) = j : y_j < y <= y_{{j+1}}` for monotonically increasing :math:`x_i` and :math:`y_j`.

Given a two-dimensional array :math:`u_{{i,j}}`, :math:`f` is then defined by
bilinear interpolation on :math:`u_{{i(x),j(y)}}, u_{{i(x)+1,j(y)}},
u_{{i(x),j(y)+1}}, u_{{i(x)+1,j(y)+1}}, if :math:`(x,y)` is in
:math:`[x_0,x_n] \times [y_0,y_m]`, linear interpolation if one of :math:`x,y`
are out of those bounds, and constant at the corner value if both are out of
bounds.
 
* `"file`" ``[string]`` HDF5 filename of the data
* `"row header`" ``[string]`` name of the row dataset, the :math:`x_i`
* `"row coordinate`" ``[string]`` one of `"t`",`"x`",`"y`",`"z`"
* `"column header`" ``[string]`` name of the column dataset, the :math:`y_i`
* `"column coordinate`" ``[string]`` one of `"t`",`"x`",`"y`",`"z`"
* `"value header`" ``[string]`` name of the values dataset, the :math:`u_{{i,j}}`

Example:

.. code-block:: xml

  <ParameterList name="function-bilinear">
    <Parameter name="file" type="string" value="pressure.h5"/>
    <Parameter name="row header" type="string" value="/time"/>
    <Parameter name="row coordinate" type="string" value="t"/>
    <Parameter name="column header" type="string" value="/x"/>
    <Parameter name="column coordinate" type="string" value="x"/>
    <Parameter name="value header" type="string" value="/pressure"/>
  </ParameterList>




Distance Function
-------------------


A distance function calculates distance from reference point :math:`x_0`
using by the following expression:

.. math::
  f(x) = \sum_{j=0}^{n} m_j (x_j - x_{0,j})^2

Note that the first parameter in :math:`x` can be time.
Here is an example of a distance function using isotropic metric:

Example:
.. code-block:: xml

  <ParameterList name="function-distance">
    <Parameter name="x0" type="Array(double)" value="{1.0, 3.0, 0.0}"/>
    <Parameter name="metric" type="Array(double)" value="{1.0, 1.0, 1.0}"/>
  </ParameterList>




Monomial Function
-------------------


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




Standard Math Function
-------------------------

These functions allow to set up non-trivial time-dependent boundary conditions 
which increases a set of analytic solutions that can be used in convergence 
analysis tests.

.. math::
  f(x) = A * operator( p * (x - s) )

or

.. math::
  f(x) = A * operator(x-s, p)

Note that these operate only on the first coordinate, which is often time.
Function composition can be used to apply these to other coordinates (or
better yet a dimension could/should be added upon request).

* `"operator`" ``[string]`` specifies the name of a standard mathematical function.
  Available options are `"cos`", `"sin`", `"tan`", `"acos`", `"asin`", `"atan`", 
  `"cosh`", `"sinh`", `"tanh`", `"exp`", `"log`", `"log10`", `"sqrt`", `"ceil`",
  `"fabs`", `"floor`", `"mod`", and `"pow`".

* `"amplitude`" ``[double]`` specifies a multiplication factor `a` in formula `a f(x)`. 
  The multiplication factor is ignored by function `mod`. Default value is 1.

* `"parameter`" ``[double]`` **1.0** specifies additional parameter `p` for
  math functions with two arguments. These functions are `"a pow(x[0], p)`"
  and `"a mod(x[0], p)`".  Alternative, scales the argument before
  application, for use in changing the period of trig functions.

* `"shift`" ``[double]`` specifies a shift of the function argument. Default is 0.

Example:

.. code-block:: xml

  <ParameterList name="function-standard-math">
    <Parameter name="operator" type="string" value="sqrt"/>
    <Parameter name="amplitude" type="double" value="1e-7"/>
    <Parameter name="shift" type="double" value="0.1"/>
  </ParameterList>

This example defines function `1e-7 sqrt(t-0.1)`.
 







