Steady-State One-Dimensional Flow with Head Boundary Conditions
---------------------------------------------------------------

Introduction
~~~~~~~~~~~~

For one-dimensional, steady-state, flow through a saturated porous medium with constant properties, 
the general governing differential equation expressing mass conservation and Darcy's law becomes simply

	.. math:: \frac{d^2h}{dx^2} = 0
		:label: ode_linear_head_head

where the total head (:math:`h`, [L]) is the sum of pressure head (:math:`P/\rho g`, [L]) 
and elevation (:math:`z`, [L])

	.. math:: h = \frac{P}{\rho g}+z

:math:`\rho` = density [M/L\ :sup:`3`\ ], :math:`g` = gravitational acceleration [L/T\ :sup:`2`\ ], 
and :math:`x` = horizontal distance [L]. The ordinary differential equation :eq:`ode_linear_head_head` is easily solved by 
direct integration as

	.. math:: h = C_1 x + C_2
		:label: generalSoln_linear_head_head

where the integration constants :math:`C_1` and :math:`C_2` depend on the boundary conditions.

Analytic solution for prescribed inlet and outlet pressures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When hydraulic head is prescribed at both boundaries as

	.. math:: 
		h(0) &= h_0\\
		h(L) &= h_L
		:label: bc_linear_head_head

the analytic solution :eq:`generalSoln_linear_head_head` for hydraulic head becomes

	.. math:: h = (h_L - h_0) \frac{x}{L} + h_0
		:label: specificSoln_linear_head_head

where :math:`L` = domain length [L]. For these boundary conditions the volumetric flowrate per unit area, 
or Darcy velocity (:math:`U`, [L/T]), is defined by Darcy's law as

	.. math:: U = -\frac{k}{\mu}\rho g \frac{dh}{dx} = -K\frac{dh}{dx} = K\frac{h_0 - h_L}{L}
		:label: DarcyVel_linear_head_head

where :math:`k` = intrinsic permeability [L\ :sup:`2`\ ],
:math:`\mu` = viscosity [M/LT], and 
:math:`K` = hydraulic conductivity [L/T]. 

Amanzi verification test problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The analytic solutions for hydraulic head and Darcy velocity can be used to test Amanzi
implementation of prescribed hydraulic head boundary conditions, Darcy's law, and mass conservation
on an elementary problem. To generate numerical results the following specifications are considered:

* Domain

	* :math:`x_{min} = y_{min} = z_{min} = 0`
	* :math:`x_{max} = 100 m, y_{max} = 2 m, z_{max} = 10 m`

* Horizontal flow in the x-coordinate direction

	* no-flow prescribed at the :math:`y_{min}, y_{max}, z_{min}, z_{max}` boundaries
	* prescribed hydraulic head at the x-coordinate boundaries: 
          :math:`h(0) = 20 \text{[m]}, h(L) = 19 \text{[m]}`

* Material properties:

	* :math:`\rho = 998.2 \: \text{[kg/m}^3\text{]}`, 
          :math:`\mu = 1.002e-3 \: \text{[Pa} \cdot \text{s]}`, 
          :math:`g = 9.807 \: \text{[m/s}^2\text{]}` 
	* :math:`K = 1.0 \text{[m/d]}` 
          (permeability: :math:`k = 1.1847 \times 10^{-12} \text{[m}^2\text{]})`

* Model discretization

	* :math:`\Delta x = 5 \text{[m]}, \Delta y = 2 \text{[m]}, \Delta z = 10 \text{[m]}`

For these input specifications, Amanzi simulation output is expected to closely match

	.. math:: h = 20 -\frac{x}{100} \text{[m]}
		:label: expectedH_linear_head_head

and

	.. math:: U = 1.0 \text{[m/d]}
		:label: expectedU_linear_head_head

following Equations :eq:`specificSoln_linear_head_head` and :eq:`DarcyVel_linear_head_head`.

Amanzi verification test results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The discretization is exact for linear solutions, and it is clear in the figure that
Amanzi has reproduced the exact solution. 

.. plot:: amanzi_linear_head_head_1d.py
   :align: center


This is also visible in the following table:

.. include:: table_values.txt


FIXME:
~~~~~~

  * probably do not want the points on the boundary 
  * we may want to plot flux?
  * values in the table need to be sorted.
