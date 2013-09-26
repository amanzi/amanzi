A Simple Flow Example
=====================

Example
--------

Saturated 3-D Steady-State Test Flow Problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Problem Definition
~~~~~~~~~~~~~~~~~~~

* Horizontal flow in the **x-direction**

* Model domain: :math:`10\text{[m]}\times 5\text{[m]}\times 5\text{[m]}`

* Constant flux (at the right boundary): :math:`\mathbf{u}(0,t)\cdot\mathbf{x} = 1.95\times 10^{-2} \text{ [kg/s/m}^2\text{]}`

* Constant head (at the left boundary): :math:`h(100,t)=120\text{ [m]}`

* Uniform Intrinsic Permeability, 
  :math:`\kappa = 1.0\times 10^{-12} \text{ [m}^2\text{]}`

* Water Viscosity,
  :math:`\mu = 1.002 \times 10^{-3} \text{ [Pa}\cdot\text{s]}`

* Water Density,
  :math:`\rho = 998.2 \text{[ kg/m}^3\text{]}`


Expected result, according to Darcy's Law
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

Darcy's Law is defined as: 

.. math:: 
   :label: DarcyFlux

   Q/A = -K \frac{dh}{dx}

which, in this case, is based on isotropic soil conductivity, *K*, and the hydraulic head gradient in the x-direction.  

Calculated head gradient (-): ?

Calculated hydraulic conductivity:

.. math:: K = \kappa \rho g / \mu = 9.67*10^{-6} \frac{m}{s}

         \frac{Q}{A} = \frac{1.95*10^{-2}}{998.2} \frac{kg/s/m^2}{kg/m^3} = 9.76*10^{-6} \frac{m}{s}

**Pressure:**  :math:`p(x,z) = [ h(x) - z ]\rho g`

Let's see if we can create a link to Darcy's law (Equation :eq:`DarcyFlux`).

Tables are the one weakness in reStructured text, but they are 
reasonable.  We will autogenerate these from the results of the
Amanzi run and the analytic solution calculations.



Plot the solutions
~~~~~~~~~~~~~~~~~~

.. plot:: prototype/steady-linear/amanzi_steady_linear.py
   :align: center

.. include:: table_values.txt

