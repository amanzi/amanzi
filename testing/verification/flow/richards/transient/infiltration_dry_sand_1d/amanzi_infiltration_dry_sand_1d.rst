Unsaturated Flow Simulations
============================

Transient, one-dimensional, unsaturated flow
--------------------------------------------

Vanderborght et al. (2005) present an analytic solution for a one-dimensional
wetting front traveling through an initially dry homogeneous sand. The test
problem is relevant to vadose zone applications with transient surface 
infiltration. Infiltration sources could include rainfall, runoff, ponding, 
surface disposal facilities (e.g., cribs), or pipe leaks. The homogeneous sand 
system is initially at low saturation under a uniform tension. A constant 
infiltration is then applied at the surface.  The wetting front moves 
downward as an intact traveling wave with the front position changing with time.

This transient one-dimensional unsaturated flow problem tests the *Amanzi* 
Richards equation process kernel. Features tested include homogeneous material 
properties, constant rate infiltration boundary condition, specified pressure 
boundary conditions, van Genuchten saturation function, and Mualem relative 
permeability function. Infiltration into a dry coarse-grained sediment such as sand
is numerically challenging because of the highly nonlinear behavior of the relative 
permeability function, which leads to a large pressure drop across a sharp wetting 
front. Computational efficiency is also being tested because small time steps are 
typically required when the infiltration begins. 

The analytical solution is based on the one-dimensional Richards (1931) equation
given by

	.. math:: \frac{d\theta(\psi)}{dt} = \frac{d}{dz} \left[ K(\psi) \left( \frac{d\psi}{dz} + 1 \right) \right]
		:label: RichardsEquation

where the :math:`z` coordinate is *positive in the upward direction* [L], 
:math:`\psi` = pressure head (:math:`P / \rho g`) [L],
:math:`K` = hydraulic conductivity [L/T], and
:math:`\theta` = volumetric water content [-].
The van Genuchten (1980) - Mualem (1976) functional forms are chosen for the 
water retention and unsaturated hydraulic conductivity curves:

	.. math:: \theta(\psi) = \theta_r + \frac{\theta_s - \theta_r}{[1+(\alpha|\psi|)^n]^{1-1/n}}
		:label: waterRetention

	.. math:: K(\psi) = K_s \frac{\lbrace1-(\alpha|\psi|)^{n-1}[1+(\alpha|\psi|)^n]^{(1-n)/n}\rbrace^2}{[1+(\alpha|\psi|)^n]^{\ell(1-1/n)}}
		:label: unsaturatedConductivity

The subscripts :math:`r` and :math:`s` denote residual and saturated conditions, 
respectively, and :math:`\alpha`, :math:`n`, and :math:`\ell` are empirical constants.


Mualem, Y. 1976. *A new model predicting the hydraulic conductivity of unsaturated porous media*. Water Resour. Res. 12:513–522.

Richards, L.A. 1931. *Capillary conduction of liquids through porous mediums*. Physics 1 (5): 318–333.

van Genuchten, M. Th. 1980. *A Closed-form Equation for Predicting the Hydraulic Conductivity of Unsaturated Soils*. Soil Sci. Soc. Am. J. 44: 892–898.

Vanderborght, J., R. Kasteel, M. Herbst, M. Javaux, D. Thiéry, M. Vanclooster, 
C. Mouvet, and H. Vereecken. 2005. *A Set of Analytical Benchmarks to Test Numerical 
Models of Flow and Transport in Soils*. Vadose Zone Journal v. 4. 206-221.
doi:10.2113/4.1.206

Amanzi verification test problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To generate numerical results the following specific conditions are considered.

	.. image:: Schematic.jpg
		:scale: 35 %
		:align: center

Input parameters for the numerical simulation are summarized as:

* Domain (3D)

	* :math:`x_{min} = y_{min} = z_{min} = 0`
	* :math:`x_{max} = y_{max} = 1 m`, :math:`z_{max} = 4 m`

* Initial conditions

	* pressure head, :math:`\psi = -4 m` (:math:`-39156 Pa` gage, :math:`62169 Pa` absolute pressure)

* Boundary conditions

	* no-flow across :math:`x` and :math:`y` coordinate boundaries
	* pressure head, :math:`\psi(0) = -4 m` (:math:`-39156 Pa` gage, :math:`62169 Pa` absolute pressure)
	* infiltration / Darcy velocity = :math:`1 m/d` at top boundary (:math:`1.11553e-2 kg/s`)

* Material properties

	* :math:`\rho = 998.2 \: kg/m^3, \mu = 1.002e-3 \: Pa\cdot s, g = 9.80665 \: m/s^2` 
	* saturated hydraulic conductivity :math:`K_s = 1 m/d`
	* equivalent *Amanzi* intrinsic permeability :math:`k = \mu K_s/\rho g = 1.18472e-11 m^2`
	* Vanderborght (2005) van Genuchten (1980) - Mualem (1976) parameters
		* :math:`\alpha = 0.15 cm^{-1}`
		* :math:`\theta_r = 0.045`
		* :math:`\theta_s = 0.43`
		* :math:`n = 3.0`
		* :math:`\ell = 0.5`
	* equivalent *Amanzi* van Genuchten (1980) - Mualem (1976) parameters
		* :math:`\alpha/\rho g = 0.0015323 Pa^{-1}`
		* residual saturation, :math:`S_r = \theta_r/\theta_s = 0.104651`
		* :math:`m = 1 - 1/n = 0.6667`
		* :math:`\ell = 0.5`

* Model discretization

	* :math:`\Delta x = \Delta y = 1.0 m` (one cell)
	* :math:`\Delta z = 0.01 m` (400 cells)

Amanzi verification test results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

(TBD) somehow insert table comparing analytic and Amanzi hydraulic head

.. include:: table_values.txt

(TBD) somehow insert plot comparing analytic and Amanzi hydraulic head

.. plot:: amanzi_steady_linear.py


Further test information
~~~~~~~~~~~~~~~~~~~~~~~~
Amanzi Test:  amanzi_infiltration_dry_sand_1d

Description:  One-Dimensional Transient Unsaturated Flow:  infiltration in an initially dry soil

Test Type: Verification

Process Kernel Tested:  Richards Equation

Test Developer(s):  Castleton/Yabusaki

Last Update:  September 23, 2013

Files:  

Unstructured: 

Input file:  amanzi_infiltration_dry_sand_1d.xml

Output file:   plot_data.h5

Structured:  none

Verification codes:  TravelingWave.cpp, Functor.hpp, Soils.hpp

Location:  /amanzi/testing/verification/flow/richards/transient/infiltration_dry_sand_1d

Hardware Environment:  poblano, New Mexico Consortium, Los Alamos, NM

		Architecture??, 

Software Environment:  ??

Parametric Studies:  No
