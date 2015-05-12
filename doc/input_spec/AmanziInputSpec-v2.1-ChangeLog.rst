=============================================================
Change Log for Amanzi XML Input Specification (Version 2.x.x)
=============================================================

.. contents:: **Table of Contents**

Overview
========

The Amanzi simulator will continue to expand the available features and enhance current functionality.  This progress will result in modifications and additions to the XML input specificition.  Specific versions of the XML input specificiation will only work with the corresponding Amanzi release version.  The following is a description of the changes made to the XML input specification and the Amanzi release version(s) that accept the input version.

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

* Added an element for specifing a petsc options file.  By default, the file named .petsc will automatically be read.  However, if the user wishes to use a different filename this option will specify that filename.  The new element is `"petsc_options_file`" and is located under `"numerical_controls`" -> `"structured_controls`".

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

* For the obeservation output options, the element `"time_macro`" has been updated to `"time_macros`" to allow users to provide a list of time macros to be utilized.


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

