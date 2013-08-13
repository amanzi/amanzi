A Simple Flow Example
=====================

Example
--------

Saturated 3-D Steady-State Test Flow Problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Problem Definition
~~~~~~~~~~~~~~~~~~~

* Horizontal flow in the **x-direction**

* Model domain:

	* :math:`X_{min} = Y_{min} = Z_{min} = 0`
	* :math:`Y_{max} = Z_{max} = 50` *meters*

* Model domain: 10 meters X 5 meters X 5 meters

* Total flow applied (flux) at **x = 0 plane** is *1.95E-2* 
  :math:`kg/s/m^2`

* Constant head specified at **x = 100 meters plane**: 120 meters 

* Uniform Intrinsic Permeability, 
  :math:`\kappa =` *1.0E-12*
  :math:`m^2`

* Water Viscosity,
  :math:`\mu =` *1.002E-3 Pa - s*

* Water Density,
  :math:`\rho =` *998.2*
  :math:`kg/m^3`

Expected result, according to Darcy's Law
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

Darcy's Law is defined as: 

.. math:: Q/A = -K \frac{dh}{dx}

which, in this case, is based on isotropic soil conductivity, *K*, and the hydraulic head gradient in the x-direction.  

Calculated head gradient (-): ?

Calculated hydraulic conductivity:

.. math:: K = \kappa \rho g / \mu = 9.67*10^{-6} \frac{m}{s}

         \frac{Q}{A} = \frac{1.95*10^{-2}}{998.2} \frac{kg/s/m^2}{kg/m^3} = 9.76*10^{-6} \frac{m}{s}

**Pressure p(x,z)** = 
:math:`[ h(x) - z ]\rho g`

Tables are the one weakness in reStructured text, but they are 
reasonable.  We will autogenerate these from the results of the
Amanzi run and the analytic solution calculations.


Plot the solutions
~~~~~~~~~~~~~~~~~~

.. plot:: prototype/steady-linear/amanzi_steady_linear.py

.. include:: table_values.txt
