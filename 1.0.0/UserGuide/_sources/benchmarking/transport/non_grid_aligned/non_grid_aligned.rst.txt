Non-grid-aligned transport with calcite precipitation along mixing interface
============================================================================

Overview and Capabilities tested
--------------------------------

This test example performs the simulation of calcite precipitation in a 2D rectangular flow domain. In this domain, two different solutions are injected at different flow rates at top and bottom half of the inlet face. As a result of the different flow rates, the directions of flow and transport are not aligned with the rectangluar grid, and the mixing zone that develops is curved. Along the curved interface, geochemical disequilibrium drives precipitation of calcite. This example tests the following capabilities: 

* 2D flow
* 2D advective transport 
* 2D dispersive transport
* Geochemical reactions

	* Mineral precipitation

For details on this test, see :ref:`about_non_grid_aligned`.
	
Background
----------

Mixing-induced precipitation is a problem that has recently received some attention in the literature, especially in pore scale studies :cite:`nga-Yoon_etal_mixing_2012`, :cite:`nga-Tartakovsky_etal_mixing_2008`. To accurately predict precipitation along mixing interfaces, accurate prediction of diffusive transport is required. This is especially true, when these transport processes are not aligned with the grid that is used to discretize the domain. Amanzi uses numerical schemes that accurately simulate non-grid-aligned advective and diffusive transport.  In this release, the advection-diffusion equation is solved using an operator splitting framework with an explicit scheme for the advection operator and an implicit scheme for diffusion. A simple operator-split approach is used; future releases will include options for higher-order schemes to couple the processes temporally.

The spatial discretization schemes available in Amanzi for advective and diffusive processes are optimized for the selected mesh framework (structured AMR or unstructured).  For unstructured meshes, Amanzi provides a second-order approximation of advective fluxes with a MUSCL-type (Monotonic Upstream-Centered Scheme for Conservation Laws) scheme :cite:`nga-Barth_lecture_1994` that combines local linear reconstruction of solute concentrations with a tensorial slope limiter.  The tensorial slope limiter adds robustness to the scheme and reduces numerical diffusion compared to the more conventional (scalar) Barth-Jesperson limiter.  The tensor diffusion fluxes are computed with either standard central differencing or a second-order MFD (Mimetic Finite Difference) method :cite:`nga-Beirao_Lipikov_Manzini_2014` optimized to guarantee solution monotonicity.  For structured AMR meshes, Amanzi implements an unsplit Godunov scheme :cite:`nga-Bell_etal_unsplit_1988`, :cite:`nga-Nonaka_etal_unsplit_3D_2011` for advection, and a standard 9-point (2D) or 27-point (3D) tensor diffusion scheme.  The Godunov scheme robustly guarantees monontonicity for non-grid-aligned advection.


Model
-----

Flow and transport 
~~~~~~~~~~~~~~~~~~

For this test, the domain is fully saturated.  The velocity field satisfies the
mass conservation equation with a Darcy approximation for volumetric fluxes:

.. math::
  \left(\frac{S_s}{g} + \frac{S_y}{Lg}\right)
    \frac{\partial p}{\partial t} 
  + \boldsymbol{\nabla}\cdot(\rho \boldsymbol{q}) = Q.

Solutes are advected with the velocity field given above, and are subject to diffusion fluxes arising 
from the properties of the solutes in the fluid medium, as well as dispersion terms due to the 
fluid flow field.  The following expresses conservation of the molar density of the i-th solute 

.. math::
  \frac{\partial (\phi s_l C_i)}{\partial t} 
  + \nabla \cdot \boldsymbol{J}_i^{\text{adv}} 
  = Q_i 
  - \nabla \cdot \boldsymbol{J}_i^{\text{disp}}

where the transport flux has the form

.. math::
  \boldsymbol{J}_i^\text{disp} = - \phi s_l \boldsymbol{D} \nabla C_i,

with the dispersion tensor :math:`boldsymbol{D}` and liquid saturation :math:`s_l`.

Calcite precipitation
~~~~~~~~~~~~~~~~~~~~~

Calcite precipitation can be described as

:math:`Ca^{2+} + CO_3^{2-} \rightarrow CaCO_3(s)`

The rate expression is 

:math:`r = S \cdot k \cdot (1 - \displaystyle\frac{Q}{K_{sp}})`

where 
:math:`S`
is the reactive surface area, 
:math:`k`
is the intrinsic rate constant, 
:math:`Q`
is the ion activity product, and
:math:`K_{sp}`
is the solubility constant of calcite. 

Problem Specification
---------------------

Flow and transport 
~~~~~~~~~~~~~~~~~~

The domain is a 2D dimensional rectangle, 60x50 cm in size. Material properties are homogeneous in the domain:

* Porosity = 0.38
* Hydraulic conductivity = 0.38 cm/s

The inlet face is divided in two equal halves. The top half is subject to the following conditions:

* flow rate = 0.50 cm/s
* :math:`[Ca^{2+}] = 0.05 M`

The bottom half is subject to the following conditions:

* flow rate = 0.26 cm/s
* :math:`[CO_3^{2-}] = 0.05 M`

The initial concentration of all species in the domain is :math:`10^{-10} M`

The dispersion coefficients for longitudinal and transverse component (relative to the local velocity vector) are:

* :math:`\alpha_{L} = \alpha_{T} = 0.0001 m`

The molecular diffusivities for all species are zeros.

Geochemistry
~~~~~~~~~~~~

The calcite mineral surface area is:

* :math:`S = 250 \text{ m}^2 \text{/m}^3`

While the instrinsic rate constant for calcite dissolution is:

* :math:`k = 10^{-11} \text{ mol/cm}^2 \text{s}`

With calcite solubility being:

* :math:`\text{log}(K_{sp}) = -8.4801`

Results and Comparison
----------------------

Expected results
~~~~~~~~~~~~~~~~

Precipitation of calcite is expected to occur in the zone where the two solutions mix. Because the flow rate in the top half is faster, the mixing zone curves downwards, into the bottom half of the domain. Since the precipitation of calcite is fast relative to transport, the mixing zone is narrow; the effective reaction rate is transport-limited.  We do not have an analytic solution for this problem. Due to the discontinuous boundary condition, we cannot expect formal convergence of either the structured AMR or unstructured algorithms to their second-order design rate.  However, we anticipate a robust, monotonic solution with minimal cross-stream diffusion/dispersion.

Simulation results
~~~~~~~~~~~~~~~~~~

Simulation results show a good agreement with expected results. Precipitation of calcite is indicated by its volume fraction at time 72 seconds (see Figure). For the structured and unstructured discretizations, the solution profiles for the precipitated calcite are similar narrow bands between the inflowing solutes.  In both case, the solute profiles are monotonic and well-behaved at all mesh resolutions.  With additional refinement (not shown), the magnitude of the peak calcite volume fraction increases (due the increased vertical gradients of precipitating solutes at the inflow boundary condition), but its concentration stays properly confined to a narrow zone at the interface. These results demonstrate that Amanzi is capable of robustly capturing non-grid-aligned processes in both the structured and unstructured mesh frameworks.

.. plot:: benchmarking/transport/non_grid_aligned/non_grid_aligned.py
   :align: left

References
----------

.. bibliography:: /bib/ascem.bib
   :filter: docname in docnames
   :style:  alpha
   :keyprefix: nga-

.. _about_non_grid_aligned:


About
-----

* Benchmark simulators: N/A
* Files

  * Amanzi input file/s (native chemistry): non_grid_aligned-u.xml, calcite_dbs.bgd
  * Amanzi input file/s (Alquimia chemistry): non_grid_aligned-u-alq.xml, calcite_dbs.bgd

* Location: testing/benchmark/chemistry/non_grid_aligned_dispersion/
* Author: K. Lipnikov, M. Day, S. Molins 
* Documentation: S. Molins
* Created on: March 10, 2014
