=============================================================
Change Log for Amanzi XML Input Specification (Version 2.x.x)
=============================================================

.. contents:: **Table of Contents**

Overview
========

The Amanzi simulator will continue to expand the available features and enhance current functionality.  This progress will result in modifications and additions to the XML input specification.  Specific versions of the XML input specification will only work with the corresponding Amanzi release version.  The following is a description of the changes made to the XML input specification and the Amanzi release version(s) that accept the input version.

Version 2.2.1 (Amanzi Release - devel - 8925 )
==================================================

Changeset 8925

NERSC build on  March 2, 2016.

Edison build path /project/projectdirs/m1012/amanzi/install/edison/mpich-7.3.1-gnu-5.2.0/Release-TPLs-0.92.23/default-160302/

Global Definitions
------------------

* Modified pattern for exponential so that expressions such as "1." and ".1" (no trailing or leading 0) are valid.  This primarily effects importing into Akuna.  This update is also included in the definition of coordinates.

process_kernels->chemistry
--------------------------

* Listing of the engine specific input file and database file have been moved the chemistry element as input_filename and database, respectively.

* The attribute algorithm for the element transport has been moved to unstr_transport_controls.  This is more consistent with our other process kernel numerical control structure.

phases
------

* The section solutes has been removed.  All non-reactive solutes are now listed as primaries.

execution_controls
------------------

* Moved the restart attribute up a level to a subelement of execution_controls.  This enforces the restriction of only 1 restart file being specified and read.

unstr_transport_controls
------------------------

* The subelement algorithm has been moved here.

unstr_chemistry_controls
------------------------

* Moved chemistry control options from the geochemistry section (and engine specific subelements) to the numerical_controls section.  This is more consistent with our other numerical control structures.  Note, only some of the options specified in subelements are valid for a given engine.  See the example.xml for comments.

geochemistry
------------

* Only verbosity and constraints are now specified under this section.  Numerical control options have been moved to the unstr_chemistry_control section.  The engine specific input file and database file have been moved to the process_kernels sections.

* The constraint attribute initial_guess has been renamed value.  Note, that this is the initial value and may be modified by equilibration in the chemistry engine.

materials
---------

* Non-reactive solutes are now named as primaries.  Reference to solute in the input file is being changed to primary.

initial_conditions/boundary_conditions
--------------------------------------

* The initial_condition/boundary_condition subelement solute_component has been renamed geochemistry_component.

* The initial_condition/boundary_condition subelement geochemistry has been removed.

Version 2.1.1 (Amanzi Release - devel - 8326 )
==================================================

Changeset 8326

NERSC build on  Aug 18, 2015.

Edison build path /project/projectdirs/m1012/amanzi/install/hopper/mpich-7.1.1-gnu-4.9.2/Release-TPLs-0.92.19/default-150818

Hopper build path /project/projectdirs/m1012/amanzi/install/hopper/mpich-7.1.1-gnu-4.9.2/Release-TPLs-0.92.19/default-150818

viz
---

* Modification to the write_regions element.  To make this field more useful for production the user can now specify a list of "field" subelements.  The name given to the field will appear in the list of available fields to visualize.  Each region listed for the given field will be assigned and colored by an integer id.


Version 2.1.1 (Amanzi Release - devel - 8245 )
==================================================

Changeset 8245

NERSC build on  Aug 6, 2015.

Edison build path /project/projectdirs/m1012/amanzi/install/hopper/mpich-7.1.1-gnu-4.9.2/Release-TPLs-0.92.19/default-150806

Hopper build path /project/projectdirs/m1012/amanzi/install/hopper/mpich-7.1.1-gnu-4.9.2/Release-TPLs-0.92.19/default-150806

str_amr_controls
----------------

* Several options take a series if integer values.  These values were indicated using a sequence of subelements named "int".  This has been updated to be a space separated list of integers within in the specific element.  The specific options require a minimum number of entries (either equal to the number of amr levels or the number of amr levels -1 ).  Any additional values will be ignored.


Version 2.1.1 (Amanzi Release - devel - 8217 )
==================================================

Changeset 8217

NERSC build on  Aug 3, 2015.

Edison build path /project/projectdirs/m1012/amanzi/install/hopper/mpich-7.1.1-gnu-4.9.2/Release-TPLs-0.92.18/default-150803

Hopper build path /project/projectdirs/m1012/amanzi/install/hopper/mpich-7.1.1-gnu-4.9.2/Release-TPLs-0.92.18/default-150803

Regions
-------

* Added tolerance attribute to the region types plane, polygonal_surface.  This attribute is optional.  It species a tolerance either side of the plane/surface that will be explored to find face centroids.


Version 2.1.1 (Amanzi Release - devel - 8153 )
==================================================

Changeset 8153

NERSC build on  July 20, 2015.

Edison build path /project/projectdirs/m1012/amanzi/install/hopper/mpich-7.1.1-gnu-4.9.2/Release-TPLs-0.92.18/default-150720

Hopper build path /project/projectdirs/m1012/amanzi/install/hopper/mpich-7.1.1-gnu-4.9.2/Release-TPLs-0.92.18/default-150720

* The version number of the schema has been updated and coincides with the 0.84 release of amanzi.  Note that a new link has been created on the NERSC machines called 'release-0.84'.  The new 'devel' link will continue to be updated as new builds are created and the input schema is updated.  The new version number of the schema for 'devel' will be 2.2.0.

Version 2.1.1 (Amanzi Release - devel - 8135 )
==================================================

Changeset 8135

NERSC build on  July 14, 2015.  

Edison build path /project/projectdirs/m1012/amanzi/install/hopper/mpich-7.1.1-gnu-4.9.2/Release-TPLs-0.92.18/default-150714  

Hopper build path /project/projectdirs/m1012/amanzi/install/hopper/mpich-7.1.1-gnu-4.9.2/Release-TPLs-0.92.18/default-150714

Build
-----

* The new python script UpdateSpec_210-211.py is now installed along with the executable and schema file in $INSTALL/bin.  This was added for convenience of users to update their input files as the schema continues to evolve.

Tools
-----

* Added the python script UpdateSpec_210-211.py to the repository in tools/install.  This script reads in an existing 2.1.0 input file and writes out an updated version consistent with the latest 2.1.0 spec (which will be updated to 2.1.1 soon).

Output
------

* Added `"vis`" element option `"write_regions`" to documentation.  This has been available for awhile but was not included in the documentation.

Version 2.1.0 (Amanzi Release - devel - 8061/8064)
==================================================

Changeset 8061/8064

NERSC build on  June 18, 2015.  

Edison build path /project/projectdirs/m1012/amanzi/install/hopper/mpich-7.1.1-gnu-4.9.2/Release-TPLs-0.92.18/default-150618  

Hopper build path /project/projectdirs/m1012/amanzi/install/hopper/mpich-7.1.1-gnu-4.9.2/Release-TPLs-0.92.18/default-150618

Numerical Controls
------------------

* Changed `"bdf1_integration_method`" attributes to elements.  This was done for consistency and readability.

* Changed element name from `"unstr_pseudo_time_integrator`" to `"unstr_initialization`".

* Added the parameter `"clipping_pressure`" to the renamed `"unstr_initialization`" list.

* Removed the parameter `"initialize_with_darcy`" from the section `"unstr_initialization`".  This section's parameters are used to initialize the steady time step and `"initialize_with_darcy`" is already specified under the steady-state section.

* Added the option `"darcy_solver`" the parameter `"method`" in the `"unstr_initialization`" list.


Version 2.1.0 (Amanzi Release - devel - 8005)
=============================================

Changeset 8005

NERSC build on  June 2, 2015.  

Edison build path /project/projectdirs/m1012/amanzi/install/hopper/mpich-7.1.1-gnu-4.9.2/Release-TPLs-0.92.17/default-150602  

Hopper build path /project/projectdirs/m1012/amanzi/install/hopper/mpich-7.1.1-gnu-4.9.2/Release-TPLs-0.92.17/default-150602

Execution Controls
------------------

* Added translation of execution time periods to `"Time Period Controls`" in the 1.2.3 input spec.

Output
------

* Both the 2.1.0 input schema and 1.2.3 input spec are moving towards using plural macros for time and cycle specifications.  This will affect vis, observations, checkpoint, and walkabout elements.  To help users transition the input translator (2.1.0 to 1.2.3) will read singular `"time_macro`" and `"cycle_macro`" and translate these to `"Time Macros`" and `"Cycle Macros`" with a single macro specified.  Also, the input parser for the unstructured algorithm will continue to read both singular and plural forms.  Note, these are temporary measures to ease transition.  Please updating input files to use the plural forms.


Version 2.1.0 (Amanzi Release - devel - 7926)
=============================================

Changeset 7926

NERSC build on  May 12, 2015.  

Edison build path /project/projectdirs/m1012/amanzi/install/hopper/mpich-7.1.1-gnu-4.9.2/Release-TPLs-0.92.17/default-150512  

Hopper build path /project/projectdirs/m1012/amanzi/install/hopper/mpich-7.1.1-gnu-4.9.2/Release-TPLs-0.92.17/default-150512

Numerical Controls
------------------

* Added missing preconditioner options under `"unstr_steady-state_controls`".  This update also included fixing the translation of the `"preconditioner`" option under `"unstr_linear_solver`"  which was incorrectly being mapped to `"steady preconditioner`" instead of `"linear solver preconditioner`".

* Cleaned up preconditioner specification for all locations.  For each unstructured numerical control with a `"preconditioner`" subelement, the valid options are the strings `"trilinos_ml`", `"hypre_amg`", or `"block_ilu`".  Options for each preconditioner have been consolidated in the subelement `"numerical_controls`" -> `"unstructured_controls`" -> `"preconditioners`".  The element `"preconditioners`" has a subelement for each preconditioner.  Each preconditioner has subelements for its specific options.

Version 2.1.0 (Amanzi Release - devel - 7688)
=============================================

Changeset 7688

NERSC build on  May 8, 2015.  

Edison build path /project/projectdirs/m1012/amanzi/install/hopper/mpich-7.1.1-gnu-4.9.2/Release-TPLs-0.92.17/default-150508  

Hopper build path /project/projectdirs/m1012/amanzi/install/hopper/mpich-7.1.1-gnu-4.9.2/Release-TPLs-0.92.17/default-150508

.. Model Description
.. -----------------

.. Definitions
.. -----------

Process Kernels
---------------

* Moved attributes from `"flow`" and `"transport`" elements that were only valid under the unstructured algorithm.  The `"flow`" attributes `"discretization_method`", `"rel_perm_method`", `"atmospheric_pressure`", and `"preconditioning_strategy`" are now subelements located under `"numerical_controls`" -> `"unstructured_controls`" -> `"unstr_flow_controls`".  The `"transport`" attributes `"algorithm`" and `"sub_cycling`" are now subelements located under `"numerical_controls`" -> `"unstructured_controls`" -> `"unstr_transport_controls`".

.. Phases
.. ------

.. Execution Controls
.. ------------------

Numerical Controls
------------------

* Added new sections under `"numerical_controls`" -> `"unstructured_controls`" for process kernel options that are specific to the unstructured algorithm.  The new sections are `"unstr_flow_controls`" and `"unstr_transport_controls`".  Options currently available were moved from the process kernels under `"process_kernels`".

    * `"discretization_method`" is now an element located under `"unstr_flow_controls`".  Valid options for this element are `"fv-default`", `"fv-monotone`", `"fv-multi_point_flux_approximation`", `"fv-extended_to_boundary_edges`", `"mfd-default`", `"mfd-optimized_for_sparsity`", `"mfd-support_operator`", `"mfd-optimized_for_monotonicity`", and `"mfd-two_point_flux_approximation`".


    * `"rel_perm_method`" is now an element located under `"unstr_flow_controls`".  Valid options for this element are `"upwind-darcy_velocity`", `"upwind-gravity`", `"upwind-amanzi`", `"other-arithmetic_average`", and `"other-harmonic_average`".  The default option is `"upwind-darcy_velocity`".

    * `"preconditioning_strategy`" is now an element located under `"unstr_flow_controls`".  Valid options for this element are `"diffusion_operator`" and `"linearized_operator`".  The default option is `"linearized_operator`".

    * `"algorithm`" is now an element located under `"unstr_transport_controls`".  Valid options for this element are `"explicit first-order`", `"explicit second-order`", and `"implicit upwind`".  The default option is `"explicit first-order`".

    * `"sub_cycling`" is now an element located under `"unstr_transport_controls`".  Valid options for this element are `"on`" and `"off`".  The default option is `"off`".

* Added an element for specifying a petsc options file.  By default, the file named .petsc will automatically be read.  However, if the user wishes to use a different filename this option will specify that filename.  The new element is `"petsc_options_file`" and is located under `"numerical_controls`" -> `"structured_controls`".

.. Geochemistry
.. ------------

.. Materials
.. ---------

.. Initial Conditions
.. ------------------

.. Boundary Conditions
.. -------------------

.. Sources
.. -------

Output
------

* For the observation output options, the element `"time_macro`" has been updated to `"time_macros`" to allow users to provide a list of time macros to be utilized.


Version 2.1.0 (Amanzi Release - devel - 7478)
=============================================

Changeset 7688


Process Kernels
---------------

* Added flow process options `"rel_perm_method`" and `"preconditioning_strategy`" as attributes.  These options are only valid for the unstructured algorithm.


Version 2.1.0 (Amanzi Release - devel - 7434)
=============================================

Changeset 7434

Materials
---------

* Stubbed in ability for file read for the material properties permeability, porosity, particle_Density, specific_storage, specific_yield, tortuosity, molecular_diffusion, viscosity, density.  Capability current available for only permeability.  
  
.. Made write_regions minOccurs=1 (why?)

Version 2.1.0 (Amanzi Release - devel - 7427)
=============================================

Changeset 7427

Output
------

* Added `"write_regions`" sub-element to the vis element. A list of regions can be given in this element similar to assigned_regions.  The specified regions will be written to the visualization file.  This is useful for debugging or easy visualization of regions for demonstrations. 

Version 2.1.0 (Amanzi Release - devel - 7412)
=============================================

Changeset 7412

Output
------
 
* Added a new observation called `"solute_volumetric_flow_rate`".  Subelements include `"filename`", `'assigned_regions`", `"functional`", `"time_macro`", and `"solute`".  The volumetric flow rat for the specified solute will be written out.


Version 2.1.0 (Amanzi Release - devel - 7335)
=============================================

Changeset 7335

Sources
-------

* Added `"diffusion_dominated_release`" as a solute component for liquid phase sources.


Version 2.1.0 (Amanzi Release - devel - 7298)
=============================================

Changeset 7298

Materials
---------

* Expanded dispersion tensor models.  New dispersion tensor types are now `"uniform_isotropic`", `"burnett_frind`", and `"lichtner_kelkar_robinson`".

Version 2.1.0 (Amanzi Release - devel - 7277)
=============================================

Changeset 7277

Numerical Controls
------------------

 * Added sub-element `"error_control_options`" to both `"unstr_steady-state_controls`" and `"unstr_pseudo_time_integrator`".


Version 2.1.0 (Amanzi Release - devel - 7266)
=============================================

Changeset 7266

Materials
---------

* Started added file read capability for `"permeability`".

Boundary Conditions
-------------------

* For hydrostatic boundary condition (uniform and linear) add attribute `"submodel`".


Version 2.1.0 (Amanzi Release - devel - 7256)
=============================================

Changeset 7256

Numerical Controls
------------------

* Added `"unstr_steady-state_controls`" subelements `"restart_tolerance_factor`" and `"restart_tolerance_relaxation_factor`".

.. Version 2.1.0 (Amanzi Release - devel - ####)
.. =============================================

.. Changeset 7688

.. NERSC build on  May 8, 2015.  

.. Edison build path /project/projectdirs/m1012/amanzi/install/hopper/mpich-7.1.1-gnu-4.9.2/Release-TPLs-0.92.17/default-150508  

.. Hopper build path /project/projectdirs/m1012/amanzi/install/hopper/mpich-7.1.1-gnu-4.9.2/Release-TPLs-0.92.17/default-150508

.. Model Description
.. -----------------

.. Definitions
.. -----------

.. Process Kernels
.. ---------------

.. Phases
.. ------

.. Execution Controls
.. ------------------

.. Numerical Controls
.. ------------------

.. Geochemistry
.. ------------

.. Materials
.. ---------

.. Initial Conditions
.. ------------------

.. Boundary Conditions
.. -------------------

.. Sources
.. -------

.. Output
.. ------

