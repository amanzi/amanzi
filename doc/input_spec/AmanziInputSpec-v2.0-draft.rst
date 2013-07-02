==============================================
Amanzi XML Input Specification (Version 2.0.x)
==============================================

.. contents:: **Table of Contents**

Overview
========

The Amanzi simulator evolves a system of conservation equations for
reacting flows in porous media, as detailed in the ASCEM report
entitled "Mathematical Formulation Requirements and Specifications for
the Process Models`" (hereafter referred to as the 'Model Requirements
Document (MRD)'). The purpose of the present document is to specify
the data required to execute Amanzi.  This specification should be
regarded as a companion to the MRD, and parameterizations of the
individual submodels are consistent between Amanzi, the MRD and this
document. Where applicable, the relevant sections of the MRD are
indicated.


Mesh
====

Amanzi supports both structured and unstructured numerical solution approaches.  This flexibility has a direct impact on the selection and design of the underlying numerical algorithms, the style of the software implementations, and, ultimately, the complexity of the user-interface.  "Mesh`" is used to select between the following options:

* `"Structured`": This instructs Amanzi to use BoxLib data structures and an associated paradigm to numerically represent the flow equations.  Data containers in the BoxLib software library, developed by CCSE at LBNL, are based on a hierarchical set of uniform Cartesian grid patches.  `"Structured`" requires that the simulation domain be a single coordinate-aligned rectangle, and that the "base mesh" consists of a logically rectangular set of uniform hexahedral cells.  This option supports a block-structured approach to dynamic mesh refinement, wherein successively refined subregions of the solution are constructed dynamically to track "interesting" features of the evolving solution.  The numerical solution approach implemented under the `"Structured`" framework is highly optimized to exploit regular data and access patterns on massively parallel computing architectures.

* `"Unstructured`": This instructs Amanzi to use data structures provided in the Trilinos software framework.  To the extent possible, the discretization algorithms implemented under this option are largely independent of the shape and connectivity of the underlying cells.  As a result, this option supports an arbitrarily complex computational mesh structure that enables users to work with numerical meshes that can be aligned with geometrically complex man-made or geostatigraphical features.  Under this option, the user typically provides a mesh file that was generated with an external software package.  The following mesh file formats are currently supported: `"Exodus 2`" (see example), `"MSTK`" (see example), `"MOAB`" (see example).  Amanzi also provides a rudmentary capability to generate unstructured meshes automatically.

.. code-block:: xml

   <mesh class=unstructured framework=["mstk"|"stk::mesh"|"moab"|"simple"]>

      <comments> May be included in the Mesh element </comments>

      <generate>
         <number_of_cells nx = "integer value"  ny = "integer value"  nz = "integer value"/>
         <box  low_coordinates = "x_low,y_low,z_low" high_coordinates = "x_high,y_high,z_high"/>
      </generate>

   </mesh>

testing.

.. code-block:: xml

  <mesh framework="mstk"> <!-- default is MSTK for unstructured -->
   <dimension>3</dimension>
   <generate>
     <number_of_cells nx = "64"  ny = "56"  nz = "107"/>
     <box  low_coordinates = "0.0,0.0,0.0" high_coordinates = "320.0,280.0,107.0"/>
   </generate>
  </mesh>


Output
======

Output data from Amanzi is currently organized into four specific groups: `"Observation Data`", `"Visualization Data`", `"Checkpoint Data`" `"Diagnostic Output`" and `"Log Data`".  
Each of these is controlled in different ways, reflecting their intended use.

* `"Checkpoint Data`" is intended to represent all that is necesary to repeat or continue an Amanzi run.  The specific data contained in a Checkpoint Data dump is specific to the algorithm optoins and mesh framework selected.  Checkpoint Data is special in that no interpolation is perfomed prior to writing the data files; the raw binary state is necessary.  As a result, the user is allowed to only write Checkpoint Data at the discrete intervals of the simulation.

* `"Visualization Data`" is intended to represent spatially complete snapshots of the solution at defined instances during the simulation.  Dependeing on the control parameters provided here, visualizatoin files may include only a fraction of the state data, and may contiain auxiliary "derived" information (see below for more discussion).

* `"Observation Data`" is intended to represent diagnostic values to be returned to the calling routine from Amanzi's simulation driver.  Observations are typically generated at arbitrary times, and frequently involve various point samplings and volumetric reductions that are interpolated in time to the desired instant.  Observations may involve derived quantities (see discussion below) or state fields.

* `"Diagnostic Output`" is intended to represent diagnostic values to be written to stdout during a simulation. The available diagnostics are for the most part analogous to what is available as observations under the Observation Data capability. 

* `"Log Data`" is intended to represent runtime diagnostics to indicate the status of the simulation in progress.  This data is typically written by the simulation code to the screen or some other stream or file pipe.  The volume of `"Log Data`" generated is a function of the `"Verbosity`" setting under `"Execution Control`".

"`Log Data`" is not explicitly controlled in this section, since it is easier to control in the context of specifying details of the algorithms.  The remaining data types are discussed in the section below.

Observations
------------

The Observations element holds all the observations that the user is
requesting from Amanzi, as well as meta data, such as the name of the
file that Amanzi will write observations to.  The observations are
collected by their phase. Thus, the ''observations'' element has the
following requirements

.. code-block:: xml

   <observations>

     Required Elements: filename, phase
     Optional Elements: NONE

   </observations>

The *filename* element contains the filename for the observation output,
and may include the full path.

.. code-block:: xml

     <filename>OptionalPath/ObservationsFileName</filename>

The *phase* element requires that the name of the phase be specified
in an attribute:

.. code-block:: xml

     <phase name="Name of Phase (Required)">

       Required Elements: NONE 
       Optional Elements: observation (one observation element block for each observation)

     </phase>

In this release the only valid phase name is ''aqueous''.  The
observation element requires a field quantity be given as an 
attribute, and elements for a region, a model (functional)
with which it will extract its source data, and a list of
discrete times for its evaluation.  The observations are evaluated
during the simulation and returned to the calling process through one
of Amanzi arguments. 

.. code-block :: xml

   <observation variable="Field Quantity (Required: see above for list of valid fields)">

     Required Elements: assigned_region, functional, one of either time_macro or cycle_macro
     Optional Elements: NONE
     
   </observation>

Here the elements are ... 



Example:

.. code-block :: xml

   <observations>
     <filename>observation.out</filename> 
       <phase name="aqueous">
         <observation variable="H+ Aqueous concentration">
           <assigned_region>Well_1</assigned_region>
           <functional>point</functional>
           <time_macro>Every year</time_macro>
         </observation>
	 <observation variable="UO2++ Aqueous concentration">
	   <assigned_region>Well_3</assigned_region>
	   <functional>point</functional>
	   <time_macro>Every year</time_macro>
	 </observation>
       </phase>
     </observations>


Observations (take 2)
---------------------

Here I'm just experimenting with a much more compressed "use case" style.  Not sure if this 
gives enough detail.

.. code-block :: xml

   <observations>
     <!-- Amanzi will write the observation to the file specified here (required) -->
     <filename>OptionalPath/ObservationsFileName</filename>
       <!-- Phases: required attribute is name -->
       <phases name="aqueous">
         <!-- Observations: List as many observations as desired,
                            Required Attribute is variable
         -->
         <observation variable="Name of field quantity from list of 'Available field quantities' defined above)">
           <assigned_region>name of region</assigned_region>
           <functional>name of observation functional (from list below)</functional>
           <time_macro>name of a time macro (from definitions)</time_macro>
         </observation>
       </phases>
   </observations>

            
