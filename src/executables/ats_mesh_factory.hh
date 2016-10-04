/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! Simple wrapper that takes a ParameterList and generates all needed meshes.

/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/


/*!
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

*/

#ifndef ATS_MESH_FACTORY_HH_
#define ATS_MESH_FACTORY_HH_

#include "Teuchos_ParameterList.hpp"
#include "State.hh"


namespace ATS {

void
createMeshes(Teuchos::ParameterList& plist,
             const Teuchos::RCP<Epetra_MpiComm>& comm,
             const Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel>& gm,
             Amanzi::State& s);


} // namespace ATS

#endif
