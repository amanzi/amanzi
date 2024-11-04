.. _amanzi_unconfined_no_recharge_1D:

Steady-State Flow in an Unconfined Aquifer: Head, and Flux Boundary Conditions
==============================================================================

Capabilities Tested
-------------------



Background
----------

Flow in an unconfined aquifer is inherently multi-dimensional, which generally
precludes simple closed-form analytic solutions. In many practical applications however,
vertical gradients and velocities are small relative to their horizontal counterparts
and can be neglected following Dupuit :cite:`ur-Dupuit_1863` and Forchheimer :cite:`ur-Forchheimer_1930`.
The Dupuit-Forchheimer theory of free-surface flow assumes that :cite:`ur-Freeze_1979`:

#. flow is horizontal and equipotential lines are vertical, and 
#. hydraulic gradient is equal to the slope of the free surface. 
   
With these assumptions and
assuming the Cartesian :math:`x`- and :math:`y`-coordinates align with the 
principal axes of the hydraulic conductivity tensor,
the general expression for saturated flow in an unconfined aquifer becomes (e.g., :cite:`ur-deMarsily_1986`)

	.. math:: \frac{\partial}{\partial x} \left[\int_\sigma^h K_{xx}dz\frac{\partial h}{\partial x}\right]
		+ \frac{\partial}{\partial y} \left[\int_\sigma^h K_{yy}dz\frac{\partial h}{\partial y}\right]
		= \omega_d \frac{\partial h}{\partial t} + Q
		:label: generalODE

where :math:`h` = hydraulic head (and elevation of the water table) [L],
:math:`t` = time [T],
:math:`K` = saturated hydraulic conductivity tensor [L/T], 
:math:`\sigma` = elevation of aquifer base [L], 
:math:`\omega_d` = specific yield or drainage porosity [-], and
:math:`Q` = volumetric flow rate per unit area withdrawn from the aquifer (sink) [L/T]. 
The hydraulic or total head (:math:`h`, [L]) is the sum of pressure head (:math:`P/\rho g`, [L]) 
and elevation (:math:`z`, [L])

	.. math:: h = \frac{P}{\rho g}+z

where :math:`\rho` = density [M/L\ :sup:`3`\ ], and :math:`g` = gravitational acceleration [L/T\ :sup:`2`\ ]. 
For constant, isotropic, hydraulic conductivity, a horizontal aquifer base at 
:math:`\sigma = z = 0`, and one-dimensional flow in the :math:`x`-direction, Equation :eq:`generalODE` becomes simply

	.. math:: \frac{d^2h^2}{dx^2} = \frac{2Q}{K}
		:label: ode

The ordinary differential equation :eq:`ode` is linear in :math:`h^2` and easily solved by 
direct integration as

	.. math:: h^2 = C_1 x + C_2
		:label: generalSoln

where the integration constants :math:`C_1` and :math:`C_2` depend on the boundary conditions
and volumetric sink :math:`Q`.  
Following :cite:`ur-Aleman_PORFLOW_2007`, Section 4.2), two particular solutions are considered
for testing *Amanzi* implementation of prescribed hydraulic head, no-flow, and recharge boundary conditions, 
Darcy's law, and mass conservation on a problem involving variably saturated conditions and quasi 2D flow
(approximately horizontal in the saturated zone and approximately vertical in the unsaturated zone).

Analytic solution #1: Prescribed head at boundaries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When hydraulic head is prescribed at both boundaries as

	.. math:: 
		h(0) &= h_0\\
		h(L) &= h_L
		:label: BCs1

the analytic solution :eq:`generalSoln` for hydraulic head becomes

	.. math:: h^2 = h_0^2 + (h_L^2 - h_0^2) \frac{x}{L} + \frac{Q_{src}L^2}{K}\left( \frac{x}{L} \right) \left(1 - \frac{x}{L} \right)
		:label: specificSoln1

where :math:`L` = domain length [L], and :math:`Q_{src} \equiv -Q`. 


Analytic solution #2: Prescribed head and no-flow boundaries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the alternative boundary conditions

	.. math:: 
		h(0) &= h_0\\
		h'(L) &= 0 \text{ (no-flow)}
		:label: BCs2

the analytic solution :eq:`generalSoln` becomes

	.. math:: h^2 = h_0^2 + \frac{Q_{src}L^2}{K}\left( \frac{x}{L} \right) \left(2 - \frac{x}{L} \right)
		:label: specificSoln2

where again :math:`Q_{src} \equiv -Q`. 


Amanzi verification test problem #1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
To generate numerical results the following specifications are considered for 
analytic solution #1 (:cite:`ur-Aleman_PORFLOW_2007`, Figure 4.2.1):

	.. image:: schematic/porflow_4.2.1.png
		:scale: 35 %
		:align: center

With the Dupuit approximation the analytic solution given by Equation :eq:`specificSoln1`
is one-dimensional in the horizontal coordinate and describes only the saturated zone. 
Because *Amanzi* does not directly solve a reduced governing equation set
equivalent to Equation :eq:`generalODE`, a two-dimensional :math:`(x,z)` simulation of the combined 
saturated and unsaturated zones using the Richards (1931) equation is required. 
Thus a vertical hydraulic conductivity and parameters defining moisture characteristic
curves for the unsaturated zone are required beyond the material properties implied by 
Equation :eq:`specificSoln1`. To minimize vertical gradients consistent with the Dupuit
assumption, the vertical hydraulic conductivity is set 10x higher than the
horizontal conductivity. To minimize non-vertical flow in the unsaturated zone
(and preserve the uniform distribution of recharge applied to the top of the model domain),
van Genuchten (1980) - Mualem (1976) parameters consistent with a gravel
are selected. Input parameters for the numerical simulation are summarized as:

* Domain (2D)

	* :math:`x_{min} = z_{min} = 0`
	* :math:`x_{max} = L = 100 ft, z_{max} = 60 ft`

* Boundary conditions

	* no-flow prescribed at the :math:`z_{min}, z_{max}` boundaries (:math:`Q_{src} = 0`)
	* prescribed hydraulic head at the x-coordinate boundaries: :math:`h(0) = 40 ft, h(L) = 20 ft`

* Material properties

	* :math:`\rho = 998.2 \: kg/m^3, \mu = 1.002e-3 \: Pa\cdot s, g = 9.807 \: m/s^2` 
	* :math:`K_{xx} = 10^{-3} ft/s`
	* :math:`K_{zz} = 10 \cdot K_{xx}` 
	* van Genuchten :cite:`ur-vanGenuchten_1980` - Mualem :cite:`ur-Mualem_1976` parameters for 
          a gravel based on :cite:`ur-Phifer_data_2006`:

		* :math:`\alpha = 0.143 cm^{-1} (1.46e-3 Pa^{-1})`
		* :math:`S_r = 0.052`
		* :math:`m = 0.314`

* Model discretization

	* :math:`\Delta x = 1 ft, \Delta z = 1 ft`

For these input specifications, *Amanzi* simulation output is expected to closely match

	.. math:: h [ft] = \sqrt{1600 - 12 x}
		:label: expectedH1

from Equation :eq:`specificSoln1`. This is demonstrated with the next figure.

.. plot:: verification/unconfined_flow/unconfined_no_recharge_1d/amanzi_unconfined_no_recharge_1d.py
   :align: center

(TBD) somehow insert table comparing analytic and Amanzi hydraulic head


Amanzi verification test problem #2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
Similarly, to generate numerical results for analytic solution #2 
the following specifications are considered (Aleman 2007, Figure 4.2.2):

.. image:: ../unconfined_recharge_1d/schematic/porflow_4.2.2.png
   :scale: 35 %
   :align: center

Input parameters for the numerical simulation are summarized as:

* Domain (2D)

	* :math:`x_{min} = z_{min} = 0`
	* :math:`x_{max} = L = 1640 ft, z_{max} = 240 ft`

* Boundary conditions

	* no-flow prescribed at the :math:`x_{max}, z_{min}` boundaries
	* recharge at :math:`z_{max}` boundary, :math:`Q_{src} = 0.0328 ft/d`
	* prescribed hydraulic head at the x-coordinate boundary: :math:`h(0) = 164 ft`

* Material properties

	* :math:`\rho = 998.2 \: kg/m^3, \mu = 1.002e-3 \: Pa\cdot s, g = 9.807 \: m/s^2` 
	* :math:`K_{xx} = 3.28 ft/d`
	* :math:`K_{zz} = 10 \cdot K_{xx}` 
	* van Genuchten :cite:`ur-vanGenuchten_1980` - Mualem :cite:`ur-Mualem_1976` parameters for 
          a gravel based on :cite:`ur-Phifer_data_2006`:

		* :math:`\alpha = 0.143 cm^{-1} (1.46e-3 Pa^{-1})`
		* :math:`S_r = 0.052`
		* :math:`m = 0.314`

* Model discretization

	* :math:`\Delta x = 1 ft, \Delta z = 2 ft`

For these input specifications, *Amanzi* simulation output is expected to closely match

	.. math:: h [ft] = 164 \sqrt{1 + \left( \frac{x}{1640} \right) \left( 2 - \frac{x}{1640} \right)}
		:label: expectedH2

from Equation :eq:`specificSoln2`. 


Amanzi verification test results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The presented results show good match of numerical data with semi-analytic solution.

.. image:: ../unconfined_recharge_1d/hydraulic_head.png
   :scale: 75 %
   :align: center


References
----------

.. bibliography:: /bib/ascem.bib
   :filter: docname in docnames
   :style:  alpha
   :keyprefix: ur-

.. _about_unconfined_no_recharge:
	    

About
-----

* Directory: testing/verification/flow/richards/steady-state/unconfined_no_recharge_1d

* Authors:  

* Maintainer(s): David Moulton (moulton@lanl.gov) 

* Input Files:

  * amanzi_unconfined_no_recharge_1d.xml, Spec 2.3
  * auto generated structured mesh



