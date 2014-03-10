Steady-State Unconfined Flow with a Seepage Boundary Condition
==============================================================

Capabilities Tested
-------------------



Background
----------

Consider the following scenario involving steady-state groundwater
flow in an unconfined aquifer that discharges to a sloped ground
surface along a seepage face ([Aleman2007]_, Section 4.4):

.. image:: schematic/porflow_4.4.1.png
   :width: 5in
   :align: center


The ground elevation slopes from 100 ft at :math:`x=0` to 50 ft at
:math:`x=L`, and the location of the seepline is unknown *a priori*.
This test case presents an opportunity to verify *Amanzi*
implementation of spatially-varying recharge and drain boundary
conditions along the top surface of the domain, as well as, correct
simulation of two-dimensional variably-saturated flow in an unconfined
aquifer and adjoining vadose zone similar to related *Amanzi*
verification cases. The test problem motivates use of a numerical grid
that conforms to the ground surface, and the simulation can also serve
as a test of *Amanzi* mass conservation etc.  on a non-orthogonal
grid.

Model
-----

An analytic solution for the elevation of the water table can be
readily derived if vertical gradients and velocities in the saturated
zone are assumed to be negligible relative to their horizontal
counterparts following [Dupuit1863]_ and [Forchheimer1930]_.  The
Dupuit-Forchheimer theory of free-surface flow specifically assumes
that [Freeze1979]_:
1) flow is horizontal and equipotential lines are vertical, and 
2) hydraulic gradient is equal to the slope of the free surface.

Let :math:`L_s` denote the unknown location of the seepline, and
:math:`h_s` the hydraulic head or height of the water table at this
location. Between the left boundary and seepline, the analytic
solution for hydraulic head in the saturated zone :math:`h` is
analogous to *Amanzi* unconfined aquifer test case #1 ([Aleman2007]_,
Equation 4.3.5; **<link to unconfined test case>**) with
:math:`(L_s,h_s)` taking the place of :math:`(L,h_L)`:

	.. math:: h^2 = h_0^2 + (h_s^2 - h_0^2) \frac{x}{L_s} + \frac{Q_{src}L_s^2}{K}\left( \frac{x}{L_s} \right) \left(1 - \frac{x}{L_s} \right), 0 \leqslant x \leqslant L_s
		:label: unconfinedLeft

The hydraulic head :math:`h` is also the height of the water table. To
the right of the seepline, any surface water is assumed to readily
drain off such that the hydraulic head or water table elevation
coincides exactly with the ground elevation, that is,

	.. math:: h = 50 ft \left(2 - \frac{x}{L}  \right), L_s \leqslant x \leqslant L
		:label: unconfinedRight

The location of the seepline is obtained by recognizing that Darcy's law and 
mass conservation across the vertical line :math:`x=L_s` requires 
([Aleman2007]_, Equations 4.4.3 and 4.4.4)

	.. math:: \frac{dh}{dx} \vert_{x=L_s^-} = \frac{1}{h_s} \left[ \frac{h_s^2 - h_0^2}{2L_s} - \frac{Q_{src} L_s}{2K} \right] = \frac{h_L - h_s}{L - L_s} = \frac{dh}{dx} \vert_{x=L_s^+}
		:label: massConstraint

where
	.. math:: h_s = 50 ft \left(2 - \frac{L_s}{L}  \right)
		:label: elevationConstraint

Simultaneous solution of Equations :eq:`massConstraint` and
:eq:`elevationConstraint` for the specific parameters defined in the
test problem schematic yields :math:`L_s = 829 ft`.


Problem Specification
~~~~~~~~~~~~~~~~~~~~~


Mesh
----
 
To conform to the physical domain depicted in the test problem
schematic, a conformal grid is used for the *Amanzi* simulation:

.. image:: mesh/porflow_4.4.3.png
   :width: 5in
   :align: center

With the Dupuit approximation the analytic solution given by Equation
:eq:`unconfinedLeft` is one-dimensional in the horizontal coordinate
and describes only the saturated zone. Because *Amanzi* does not
directly solve a reduced governing equation set based on the Dupuit
assumption, a two-dimensional :math:`(x,z)` simulation of the combined
saturated and unsaturated zones using the [Richards1931]_ equation is
required. Thus a vertical hydraulic conductivity and parameters
defining moisture characteristic curves for the unsaturated zone are
required beyond the material properties implied by Equation
:eq:`unconfinedLeft`.  Input parameters for the numerical simulation
are summarized as:

* Domain (2D)

	* :math:`x_{min} = z_{min} = 0`
	* :math:`x_{max} = L = 1000 ft`
	* :math:`z_{max} = 100 ft` at :math:`x = 0` and :math:`50 ft` at :math:`x = L`

* Boundary conditions

	* no-flow prescribed at the :math:`z_{min}` boundary
	* prescribed hydraulic head: :math:`h(0) = 80 ft, h(L) = 50 ft`
	* recharge along the top surface = 1 ft/y for :math:`0 \leqslant x \leqslant L_s`

* Material properties

	* :math:`\rho = 998.2 \: kg/m^3, \mu = 1.002 \times 10^{-3} \: Pa\cdot s, g = 9.807 \: m/s^2` 
	* hydraulic conductivity :math:`K = 1 ft/d`
	* van Genuchten [vanGenuchten1980]_ - Mualem [Mualen1976]_ parameters
		* :math:`\alpha = 1.0212e-04 Pa^{-1}`
		* :math:`S_r = 0.25`
		* :math:`m = 0.09090`

* Model discretization

	* :math:`\Delta x = 25 ft`
	* variable: :math:`2.5 ft \leqslant \Delta z \leqslant 5 ft`

Results and Comparison
~~~~~~~~~~~~~~~~~~~~~~


 .. image:: figures/hydraulic_head.png
    :width: 4in
    :align: center


.. .. include:: table_values.txt


References
----------

.... [Aleman2007] Aleman, S. 2007. *PORFLOW Testing and Verification Document*. Savannah River National Laboratory technical report WSRC-STI-2007-00150 Rev 0. 193 p.
 
.. [Dupuit1863] Dupuit, J. 1863. *Estudes Thèoriques et Pratiques sur le mouvement des Eaux dans les canaux dècouverts et à travers les terrains permèables* (Second Edition ed.). Paris: Dunod.

.. [Forchheimer1930] Forchheimer, P. 1930. *Hydraulik*. 3d ed., Teubner, Leipzig, Berlin. 

.. [Freeze1979] Freeze, R.A. and J.A. Cherry. 1979. *Groundwater*. Prentice-Hall, Inc. Englewood Cliffs, NJ. 604 p.

.. [Mualen1976] Mualem, Y. 1976. *A new model predicting the hydraulic conductivity of unsaturated porous media*. Water Resour. Res. 12:513-522.

.. [Richards1931] Richards, L.A. 1931. *Capillary conduction of liquids through porous mediums*. Physics 1 (5): 318-333.

.. [vanGenuchten1980] van Genuchten, M. Th. 1980. *A Closed-form Equation for Predicting the Hydraulic Conductivity of Unsaturated Soils*. Soil Sci. Soc. Am. J. 44: 892-898.


About
-----

Status
------
