=============================================================
Change Log for Amanzi XML Input Specification (Version 2.x.x)
=============================================================

.. contents:: **Table of Contents**

Overview
========

The Amanzi simulator will continue to expand the available features and enhance current functionality.  This progress will result in modifications and additions to the XML input specificition.  Specific versions of the XML input specificiation will only work with the corresponding Amanzi release version.  The following is a description of the changes made to the XML input specification and the Amanzi release version(s) that accept the input version.

Version 2.1.0 (Amanzi Release - devel)
======================================

Model Description
-----------------

Definitions
-----------

Process Kernels
---------------

* Moved attributes from `"flow`" and `"transport`" elements that were only valid under the unstructured algorithm.  The `"flow`" attributes `"discretization_method`", `"rel_perm_method`", `"atmospheric_pressure`", and `"preconditioning_strategy`" are now located under `"numerical_controls`" -> `"unstructured_controls`" -> `"unstr_flow_controls`".  The `"transport`" attributes `"algorithm`" and `"sub_cycling`" are now located under `"numerical_controls`" -> `"unstructured_controls`" -> `"unstr_transport_controls`".

Phases
------

Execution Controls
------------------

Numerical Controls
------------------

* Added new sections under `"numerical_controls`" -> `"unstructured_controls`" for process kernel options that are specific to the unstructured algorithm.  The new sections are `"unstr_flow_controls`" and `"unstr_transport_controls`".  Options currently available were moved from the process kernels under `"process_kernels`".
* `"discretization_method`" is now an element located under `"unstr_flow_controls`".  Valid options for this element are `"fv-default`", `"fv-monotone`", `"fv-multi_point_flux_approximation`", `"fv-extended_to_boundary_edges`", `"mfd-default`", `"mfd-optimized_for_sparsity`", `"mfd-support_operator`", `"mfd-optimized_for_monotonicity`", and `"mfd-two_point_flux_approximation`".
* `"rel_perm_method`" is now an element located under `"unstr_flow_controls`".  Valid options for this element are `"upwind-darcy_velocity`", `"upwind-gravity`", `"upwind-amanzi`", `"other-arithmetic_average`", and `"other-harmonic_average`".  The default option is `"upwind-darcy_velocity`".
* `"preconditioning_strategy`" is now an element located under `"unstr_flow_controls`".  Valid options for this element are `"diffusion_operator`" and `"linearized_operator`".  The default option is `"linearized_operator`".
* `"algorithm`" is now an element located under `"unstr_transport_controls`".  Valid options for this element are `"explicit first-order`", `"explicit second-order`", and `"implicit upwind`".  The default option is `"explicit first-order`".
* `"sub_cycling`" is now an element located under `"unstr_transport_controls`".  Valid options for this element are `"on`" and `"off`".  The default option is `"off`".

Geochemistry
------------

Materials
---------

Initial Conditions
------------------

Boundary Conditions
-------------------

Sources
-------

Output
------
