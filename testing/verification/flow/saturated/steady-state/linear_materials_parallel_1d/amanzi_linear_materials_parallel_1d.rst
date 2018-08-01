Steady-State One-Dimensional Flow: Materials in Parallel
========================================================

Capabilities Tested
-------------------

This one-dimensional model of steady-state flow through a saturated porous
medium with constant properties tests the Amanzi
implementation of prescribed hydraulic head boundary conditions,
Darcy's law, and mass conservation on an elementary problem with discrete heterogeneity.  
Capabilities tested include:
  
  * one-dimensional representation
  * steady-state
  * saturated flow
  * heterogenous porous medium
  * prescribed hydraulic head boundary conditions

For details on this test, see :ref:`about_linear_materials_parallel_1d`.


Background
----------
For one-dimensional, steady-state, flow through a saturated porous medium with constant properties, 
the general governing differential equation expressing mass conservation and Darcy's law :cite:`matp-Darcy_1856` becomes simply

	.. math:: \frac{d^2h}{dx^2} = 0
		:label: ode_materials_parallel

where the total head (:math:`h`, [L]) is the sum of pressure head (:math:`P/\rho g`, [L]) 
and elevation (:math:`z`, [L])

	.. math:: h = \frac{P}{\rho g}+z

:math:`\rho` = density [M/L\ :sup:`3`\ ], :math:`g` = gravitational acceleration [L/T\ :sup:`2`\ ], 
and :math:`x` = horizontal distance [L]. The ordinary differential equation :eq:`ode_materials_parallel` is easily solved by 
direct integration as

	.. math:: h = C_1 x + C_2
		:label: generalSoln_materials_parallel

where the integration constants :math:`C_1` and :math:`C_2` depend on the boundary conditions.

For a simple heterogeneous porous medium composed of two constant-property materials in parallel, 
Equation :eq:`generalSoln_materials_parallel` can be applied to each subregion separately. To analyze this 
special case, let the subscripts *1* and *2* denote the two subregions.



Model
-----
The analytic solution for prescribed inlet and outlet pressures is shown below.
When hydraulic head is prescribed at both boundaries as

	.. math:: 
		h(0) &= h_0\\
		h(L) &= h_L
		:label: bc_materials_parallel

the analytic solution :eq:`generalSoln` for hydraulic head in each subregion (:math:`h_i`, [L]) becomes

	.. math:: 
		h_i = (h_L - h_0) \frac{x}{L} + h_0, i=1,2
		:label: specificSoln_materials_parallel

where :math:`L` = domain length [L]. The volumetric flowrate per unit area through a porous medium, 
or Darcy velocity (:math:`U`, [L/T]), is defined by Darcy's law as

	.. math:: U = -\frac{k}{\mu\rho g}\frac{dh}{dx} = -K\frac{dh}{dx}
		:label: DarcyVel_materials_parallel

where :math:`k` = intrinsic permeability [L\ :sup:`2`\ ],
:math:`\mu` = viscosity [M/LT], and 
:math:`K` = hydraulic conductivity [L/T]. 
Applying Equation :eq:`DarcyVel_materials_parallel` to each subregion using Equation :eq:`specificSoln_materials_parallel` yields

	.. math:: 
		U_i = K_i\frac{h_0 - h_L}{L}, i=1,2
		:label: specificDarcyVel_materials_parallel

Note that the hydraulic head and Darcy velocity in each subregion are independent of the properties of
the other subregion.


Problem Specification
---------------------
The analytic solutions for hydraulic head and Darcy velocity can be used to test Amanzi
implementation of prescribed hydraulic head boundary conditions, Darcy's law, and mass conservation
on an elementary problem with discrete heterogeneity.


Schematic
~~~~~~~~~
The domain is shown in the following schematic.

.. figure:: schematic/schematic.png 
    :figclass: align-center
    :width: 400 px

    **One-dimensional, steady-state flow through a saturated porous medium with constant properties**


Mesh
~~~~
A steady-flow mesh is applied. The mesh consists of 400 cells: 20 grid cells in the x-direction, 2 cells in the y-direction, and 1 cell in the z-direction. Mesh discretization is as follows: :math:`\Delta x = 5 \: m, \: \Delta y = 1 \: m,` and :math:`\Delta z = 10 \: m`. 


Variables
~~~~~~~~~
To generate numerical results the following specifications are considered:

* Domain

	* :math:`x_{min} = y_{min} = z_{min} = 0`
	* :math:`x_{max} = 100 \: m, \: y_{max} = 2 \: m, \: z_{max} = 10 \: m`

* Horizontal flow in the x-coordinate direction

	* no-flow prescribed at the :math:`y_{min}, \: y_{max}, \: z_{min}, \: z_{max}` boundaries
	* prescribed hydraulic head at the x-coordinate boundaries: :math:`h(0) = 20 \: m, \: h(L) = 19 \: m`

* Material properties:

	* :math:`\rho = 998.2 \: kg/m^3, \:  \mu = 1.002 \times 10^{-3} \: Pa\cdot s, \: g = 9.807 \: m/s^2` 
	* :math:`K_1 = 1.0 \: m/d` :math:`(k = 1.1847 \times 10^{-12} \: m^2)` for :math:`0 \: m \leqslant y \leqslant 1 \: m`
	* :math:`K_2 = 10 \: m/d` :math:`(k = 1.1847 \times 10^{-11} \: m^2)` for :math:`1 \: m \leqslant y \leqslant 2 \: m`

* Model discretization

	* :math:`\Delta x = 5 \: m, \Delta y = 1 \: m, \Delta z = 10 \: m`

For these input specifications, Amanzi simulation output is expected to closely match

	.. math:: h_i = 20m -\frac{x}{100m}, \: i=1,2
		:label: expectedH_materials_parallel

and

	.. math:: 
		U_1 &= 0.01 \: m/d\\
		U_2 &= 0.1 \: m/d
		:label: expectedU_materials_parallel

following Equations :eq:`specificSoln_materials_parallel` and :eq:`specificDarcyVel_materials_parallel`.


Results and Comparison
----------------------
The discretization is exact for linear solutions, and it is clear in the figure that
Amanzi has reproduced the exact solution.

.. plot:: amanzi_linear_materials_parallel_1d.py

This is also visible in the following table.

.. include:: table_values.txt


References
----------

.. bibliography:: /bib/ascem.bib
   :filter: docname in docnames
   :style:  alpha
   :keyprefix: matp- 


.. _about_linear_materials_parallel_1d:

About
-----
* Directory: testing/verification/flow/saturated/steady-state/linear_materials_parallel_1d

* Authors:  Greg Flach

* Maintainer(s): David Moulton, moulton@lanl.gov

* Input Files:

  * amanzi_linear_materials_parallel_1d-s.xlm 

    * Spec Version 2.3.0, structured mesh framework
    * mesh:  steady-flow_mesh.h5
    * runs

  * amanzi_linear_materials_parallel_1d-u.xml

    * Spec Version 2.3.0, unstructured mesh framework
    * mesh:  generated in file
    * runs

* Mesh Files:

  * steady-flow_mesh.h5
  * unstructured mesh is generated in file

* Analytic solution computed with golden output

  * Subdirectory: golden_output

  * Input Files:
  
    * steady-flow_data.h5

Status
~~~~~~
.. todo:: 

  * Documentation:
  * keb: Is this really 1D flow or is it horizontal flow with 2 dimensions?
  * keb: List what is expected out of Amanzi simulation output.
