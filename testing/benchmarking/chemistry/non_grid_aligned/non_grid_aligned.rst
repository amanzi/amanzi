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

Mixing-induced precipitation is a problem that has recently received some attention in the literature, especially in pore scale studies [Yoon2012]_, [Tartakovsky2008]_. To accurately predict precipitation along mixing interfaces, accurate prediction of diffusive-dispersive transport. This is especially true, when these transport processes are not aligned with the grid that is used to discretize the domain. Amanzi uses numerical schemes that are capable to simulate accurately non-grid-aligned processes such as advective and dispersive transport when a rectangular grid is used.

[Add description here of Konstantin's and Marc's numerical schemes].

Model
-----

Flow and transport 
~~~~~~~~~~~~~~~~~~

Flow is solved with the single-phase saturated flow equation

.. math::
  \left(\frac{S_s}{g} + \frac{S_y}{Lg}\right)
    \frac{\partial p}{\partial t} 
  + \boldsymbol{\nabla}\cdot(\rho \boldsymbol{q}) = Q

Transport is solved with the advection-dispersion equation

.. math::
  \frac{\partial (\phi s_l C_i)}{\partial t} 
  + \nabla \cdot \boldsymbol{J}_i^{\text{adv}} 
  = Q_i 
  - \nabla \cdot \boldsymbol{J}_i^{\text{disp}}

where the dispersive flux has the form

.. math::
  \boldsymbol{J}_i^\text{disp} = - \phi s_l \boldsymbol{D} \nabla C_i

Calcite precipitation
~~~~~~~~~~~~~~~~~~~~~

Calcite precipitation can be described as

:math:`Ca^{2+} + CO_3^{2-} \rightarrow CaCO_3(s)`

The rate expression is 

:math:`r = S \cdot k \cdot (1 - \frac{Q}{K_{sp}})`

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

While the bottom half is subject to the following conditions:

* flow rate = 0.26 cm/s
* :math:`[CO_3^{2-}] = 0.05 M`

Initially, the concentration of all species in the domain is :math:`10^{-10} M`

The dispersion coefficients for longitudinal and transverse component are:

* :math:`\alpha_{L} = \alpha_{T} = 0.0001 m`

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

Precipitation of calcite is expected to occur in the zone where the two solutions mix. Because the flow rate in the top half is faster, the mixing zone curves downwards and is located in the bottom half of the domain. Because the precipitation of calcite is relatively fast relative to transport, the mixing zone is relatively narrow with the effective reaction rate being mostly transport-limited.

Simulation results
~~~~~~~~~~~~~~~~~~

Simulation results show a good agreement with expected results. Precipitation of calcite is indicated by its volume fraction at time 72 seconds (see Figure). This result demonstrates that the handling of dispersion in Amanzi is capable of capturing non-grid-aligned processes correctly.

.. plot:: prototype/chemistry/non_grid_aligned/plot_non_grid_aligned.py

..   :align: left

References
----------

.. [Yoon2012] Hongkyu Yoon, Albert J. Valocchi, Charles J. Werth, and Thomas Dewers (2012) Pore-scale simulation of mixing-induced calcium carbonate precipitation and dissolution in a microfluidic pore network, Water Resour. Res., 48, W02524, doi:10.1029/2011WR011192.
.. [Tartakovsky2008] A. M. Tartakovsky, G. Redden, P. C. Lichtner, T. D. Scheibe, and P. Meakin (2008) Mixing-induced precipitation: Experimental study and multiscale numerical analysis, Water Resour. Res., 44, W06S04, doi:10.1029/2006WR005725.

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
