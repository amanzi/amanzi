Steady-State Unconfined Flow with a Seepage Boundary Condition
==============================================================

Capabilities Tested
-------------------

This one-dimensional flow problem tests the Amanzi unsaturated (Richards) flow process kernel with spatially-varying recharge and drain boundary
conditions in an unconfined aquifer.
Capabilities tested include:

  * single-phase, one-dimensional flow
  * steady-state flow
  * variably-saturated (both saturated and unsaturated) flow
  * constant-head (Dirichlet) boundary conditions
  * no-flow (Neumann) boundary conditions
  * spatially-varying recharge/drain (Cauchy) boundary conditions
  * homogeneous porous medium
  * isotropic porous medium
  * non-uniform, non-orthogonal mesh

For details on this test, see :ref:`about_unconfined_seepage`. 


Background
----------

Unconfined aquifers are often fed by surface recharge (e.g. precipitation), which can be
represented in numerical models as a Cauchy boundary condition. Similarly, if the water table
intersects the surface topography, perhaps due to slope of the ground surface, a seepage
face develops that can likewise be represented with the Cauchy boundary condition.

Flow in an unconfined aquifer is inherently multi-dimensional, which generally precludes 
simple closed-form analytic solutions. In many practical applications however, vertical 
gradients and velocities are small relative to their horizontal components and can be 
neglected following Dupuit :cite:`us-Dupuit_1863` and Forchheimer :cite:`us-Forchheimer_1930`.
The Dupuit-Forchheimer theory of free-surface flow assumes that :cite:`us-Freeze_1979`:

#. flow is horizontal and equipotential lines are vertical, and 
#. hydraulic gradient is equal to the slope of the free surface. 

In this problem, we present an analytical solution for the elevation of the water table in an 
unconfined aquifer where the upper bounds are at once draining along a seepage face and being
fed by areal recharge. Analytical solution results are compared to those generated using *Amanzi*.


Model
-----

An analytic solution for the elevation of the water table can be
readily derived if we make Dupuit-Forchheimer assumptions of free-surface flow --- that is, 
that flow is horizontal and the hydraulic gradient is equal to the slope of the free surface 
:cite:`us-Freeze_1979`.

Let :math:`L_s` [m] denote the unknown location of the seepline, and
:math:`h_s` [m] denote the hydraulic head or height of the water table at this
location. Between the left boundary and seepline, the analytic
solution for hydraulic head in the saturated zone (:cite:`us-Aleman_PORFLOW_2007`, Equation 4.3.5) :math:`h` is
analogous to *Amanzi* unconfined aquifer test case #1 (:ref:`amanzi_unconfined_no_recharge_1D`), with
:math:`(L_s,h_s)` taking the place of :math:`(L,h_L)`:

	.. math:: h^2 = h_0^2 + (h_s^2 - h_0^2) \frac{x}{L_s} + \frac{Q_{src}L_s^2}{K}\left( \frac{x}{L_s} \right) \left(1 - \frac{x}{L_s} \right),\: 0 \leqslant x \leqslant L_s,
		:label: unconfinedLeft

where :math:`Q_{src}=-Q` the seepline is position at :math:`x=L_s` [m].

The hydraulic head, :math:`h` [m], is also the height of the water table. To
the right of the seepline, any surface water is assumed to readily
drain off such that the hydraulic head or water table elevation
coincides exactly with the ground elevation, that is,

	.. math:: h = 50 \text{ ft} \left(2 - \frac{x}{L}  \right),\: L_s \leqslant x \leqslant L
		:label: unconfinedRight

The location of the seepline is obtained by recognizing that Darcy's law and 
mass conservation across the vertical line :math:`x=L_s` requires 
(:cite:`us-Aleman_PORFLOW_2007`, Equations 4.4.3 and 4.4.4)

	.. math:: \frac{dh}{dx} \vert_{x=L_s^-} = \frac{1}{h_s} \left[ \frac{h_s^2 - h_0^2}{2L_s} - \frac{Q_{src} L_s}{2K} \right] = \frac{h_L - h_s}{L - L_s} = \frac{dh}{dx} \vert_{x=L_s^+},
		:label: massConstraint

where 
     
where
	.. math:: h_s = 50 \text{ ft} \left(2 - \frac{L_s}{L}  \right)
		:label: elevationConstraint

Simultaneous solution of Equations :eq:`massConstraint` and
:eq:`elevationConstraint` for the specific parameters defined in the
test problem schematic yields :math:`L_s = 829 ft`.


Problem Specification
---------------------


Schematic
~~~~~~~~~

Consider the following scenario involving steady-state groundwater
flow in an unconfined aquifer that discharges to a sloped ground
surface along a seepage face (:cite:`us-Aleman_PORFLOW_2007`, Section 4.4):

.. figure:: schematic/porflow_4.4.1.png
    :figclass: align-center
    :width: 5in

    **Schematic of the model domain and boundary conditions.**

The ground elevation slopes from 100 ft at :math:`x=0` to 50 ft at
:math:`x=L`, and the location of the seepline is unknown *a priori*.


Mesh
~~~~

The mesh uses constant grid discretization in the x-direction (:math:`\Delta x = 25` ft) and variable discretization in the z-direction (:math:`2.5` ft :math:`\leqslant \Delta z \leqslant 5` ft).
To conform to the physical domain depicted in the test problem
schematic, a conformal grid is used for the *Amanzi* simulation:

.. figure:: mesh/porflow_4.4.3.png
    :figclass: align-center
    :width: 5in

    **Unstructured two-dimensional mesh with upper bounds conforming to ground surface.**


Variables
~~~~~~~~~

With the Dupuit approximation the analytic solution given by Equation
:eq:`unconfinedLeft` is one-dimensional in the horizontal coordinate
and describes only the saturated zone. Because *Amanzi* does not
directly solve a reduced governing equation set based on the Dupuit
assumption, a two-dimensional :math:`(x,z)` simulation of the combined
saturated and unsaturated zones using the :cite:`us-Richards_1931` equation is
required. Thus a vertical hydraulic conductivity and parameters
defining moisture characteristic curves for the unsaturated zone are
required beyond the material properties implied by Equation
:eq:`unconfinedLeft`.  Input parameters for the numerical simulation
are summarized as:

* Domain (2D)
      
  * :math:`x_{min} = z_{min} = 0 \text{ [ft]}`
  * :math:`x_{max} = L = 1000 \text{ [ft]}`
  * :math:`z_{max} = 100` at :math:`x = 0` and :math:`50` at :math:`x = L \text{ [ft]}`

* Material properties
  
  * isotropic hydraulic conductivity:     :math:`K = 1 \text{ [ft/d]} =  3.528 \times 10^{-6} \text{ [m/s]}`

    * derived from:    :math:`K=\frac{k \rho g}{\mu}`, where permeability :math:`k = 3.6098 \times 10^{-13} \text{ [m}^2\text{]}` (1 ft = 0.3048 ft)

  * porosity:    :math:`\phi = 0.3`
  * fluid density:    :math:`\rho = 998.2 \: \text{[kg/m}^3\text{]}`
  * dynamic viscosity:    :math:`\mu = 1.002 \times 10^{-3} \: \text{[Pa} \cdot \text{s]}` 
  * gravitational acceleration:    :math:`g = 9.807 \: \text{[m/s}^2\text{]}` 
  * van Genuchten :cite:`us-vanGenuchten_1980` - Mualem :cite:`us-Mualem_1976` parameters

    * :math:`\alpha = 1.0212 \times 10^{-4} \text{ [Pa}^{-1} \text{]}`
    * :math:`S_r = 0.25`
    * :math:`m = 0.09090`

* Boundary conditions

  * no-flow (Neumann) boundary condition prescribed at :math:`z_{min}`
  * prescribed hydraulic head (Dirichlet): :math:`h(0) = 80,\: h(L) = 50 \text{ [ft]}`
  * recharge along the top surface = 1 ft/y for :math:`0 \leqslant x \leqslant L_s`


Results and Comparison
----------------------

 .. image:: figures/hydraulic_head.png
    :width: 4in
    :align: center

.. include:: table_values.txt


References
----------

.. bibliography:: /bib/ascem.bib
   :filter: docname in docnames
   :style:  alpha
   :keyprefix: us-
		

.. _about_unconfined_seepage:
	
About
-----

* Directory:  testing/verification/flow/richards/steady-state/unconfined_seepage_1d

* Authors:  Markus Berndt

* Maintainer:  David Moulton (moulton@lanl.gov)

* Input Files:

  * amanzi_unconfined_seepage_1d-u.xml

    * Spec Version 2.3, unstructured mesh framework
    * mesh:  porflow4_4.exo

      * two-dimensional mesh with conformal (non-orthogonal) grid

.. * Mesh Files:

  .. * porflow4_4.exo
 
    * two-dimensional mesh with conformal grid


.. Status
.. ------
