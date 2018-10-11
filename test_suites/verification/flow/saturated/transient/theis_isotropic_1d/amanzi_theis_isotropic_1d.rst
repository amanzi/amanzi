Transient One-Dimensional Confined Flow to a Pumping Well (Theis)
=================================================================

Capabilities Tested
-------------------

This transient one-dimensional (radial) flow problem tests the Amanzi flow process kernel. Capabilities tested include:

  * single-phase, one-dimensional flow 
  * transient flow
  * saturated flow
  * constant-rate pumping well
  * constant-head (Dirichlet) boundary conditions
  * specified volumetric flux (Neumann) boundary conditions
  * homogeneous porous medium
  * isotropic porous medium
  * uniform mesh

For details on this test, see :ref:`about_theis`.


Background
----------

Groundwater resource evaluation is crucial to understanding the
concept of groundwater yield. This concept is vital because it helps
determine the maximum allowable pump rate out of an aquifer and
the consequences on the water table.  For example,
pumping water out of the water table (unconfined aquifer) may dry up
nearby wells due to the fall in the saturated thickness of the
aquifer. This can be prevented by modeling the fall of the water table
given a certain pumping rate.

To apply the analysis developed by Theis, three parameters must be known:

* :math:`T`, transmissivity
* :math:`S`, storativity
* :math:`Q`, well-head pumping rate

Transmissivity [L\ :sup:`2`/T] of the aquifer is defined as: 

.. math:: T = Kb

where :math:`K` is the hydraulic conductivity of the aquifer  [L/T] and :math:`b` is the
saturated thickness of the aquifer [L].  Transmissivity values greater than 0.015 
:math:`\frac{m^2}{s}` represent aquifers capable of well exploitation.
Storativity is a dimensionless parameter that describes the amount of
water released by the aquifer per unit unit decline in hydraulic head in the aquifer,
per unit area of the aquifer.
Storativity can be calculated using:

.. math:: S = S_s b

where :math:`S_s` is the *specific storage* [L\ :sup:`-1`], which describes
the volume of water an aquifer releases from storage, per unit volume of aquifer, per 
unit decline in hydraulic head. For a confined aquifer, storativity is the vertically
integrated specific storage value.

Lastly, the constant pumping rate, :math:`Q`, is the volume of water
discharged from the well per unit time [L\ :sup:`3`/T].


Model
-----

Theis (1935) developed an analytical solution for transient (non-
steady state) drawdown for a fully penetrating well by imposing the
boundary conditions: :math:`h = h_0` for :math:`t = 0` and 
:math:`h \Rightarrow h_0` as :math:`r \Rightarrow \infty` :cite:`theis-Theis1935relation`. The equation assumes an infinite and uniform confined aquifer.  When a well is
pumped the water table declines toward the well and flow is induced
toward the well from all directions. Theoretically, this flow can be
idealized by purely radial symmetric flow and can be decribed by the
equation below (this is analogous to heat flow by conduction developed
by Fourier):

.. math:: 
     \frac{\partial^2 h}{\partial r^2} 
   + \frac{1}{r} \frac{\partial h}{\partial r} 
   = \frac{S}{T} \frac{\partial h}{\partial t},

where :math:`h` is hydraulic head [m], :math:`r` is radial distance from 
the well [m], :math:`S` is the aquifer storativity [-], :math:`T` is the 
aquifer transmissivity [m\ :sup:`2`/s], and :math:`t` is time [s]. 

The analytical solution of drawdown as a function of time and distance
is found to be:

.. math:: s = h(r,0) - h(r,t) = \frac{Q }{4 \pi T} W(u),

where

.. math:: W(u) = \int_u^\infty \frac{e^{-\psi}}{\psi} d\psi,

and

.. math:: u(r,t) = \frac{r^2 S}{4 T t},

for which :math:`s` is head drawdown [m], :math:`Q` is the volumetric discharge
rate of the well [L\ :sup:`3`/T], :math:`W(u)` is the well function,
and :math:`u` is the argument of the well function.  The integral approximation 
for the well function, *W(u)*, for :math:`0 < u < 1` is

.. math::
      W(u) = -0.577 - log(u) + .99 u - 0.2499 u^2 +0.055 u^3 - 0.00976 u^4 +
      0.00108 u^5
   
and for :math:`u \geq 1`,

.. math:: W(u) =(\frac{exp[-u]}{u})  \frac{u^4 +8.57333 u^3  +18.05902
	  u^2 + 8.63476 u +  0.26777}{u^4 + 9.57332 u^3 + 25.63296
	  u^2 + 21.09965 u + 3.95850}  

Note, :math:`W(u)` can easily be determined from existing tables once
:math:`u(r, t)` is found which is a measure of aquifer response time. For a
given value of *t* one can construct a draw down curve with respect to
the distance from the pumping well, *r*. 


.. _plot_table_Theis:


Problem Specification
---------------------


Schematic
~~~~~~~~~

Note, the values in the schematic correlate to the values found in
:ref:`Variables`.

.. figure:: schematic/Theis.png 
    :figclass: align-center

    **Illustration of transient drawdown in a confined aquifer.**
		    

Mesh
~~~~

The mesh is generated in the input file. It consists of cells with size :math:`\Delta x=4` m, :math:`\Delta y=4` m, and :math:`\Delta z=10` m. It has 600 grid cells in the x-direction and 600 grid cells in the y-direction, with a single cell in the z-direction.  

.. .. figure:: figures/mesh.png
    :figclass: align-center
    :width: 600 px

    .. **Computational mesh with 360,000 cells.**

.. _Variables:

Variables
~~~~~~~~~

* Domain:

  * :math:`x_{min} = y_{min} = z_{min} = -1200` m
  * :math:`x_{max} = 1200` m, :math:`y_{max} = 1200` m, :math:`z_{max} = 10` m

* Boundary and initial conditions:

  * initial hydraulic head:   :math:`h(r,0)=20.0 \: \text{[m]}`
  * constant-head (Dirichlet) lateral boundary conditions:   :math:`h(x_{max},t)=h(y_{max},t)=20.0 \: \text{[m]}`
  * no-flow (Neumann) upper and lower boundary conditions
  * well-head pumping rate:   :math:`Q=4.0 \: \text{[m}^3\text{/s]}`

* radial distance from well:    :math:`r \: \text{[m]}`
* duration of pumping time:    :math:`t \: \text{[s]}`

* Material properties:

  * fluid density:    :math:`\rho = 1000.0 \: \text{[kg/m}^3\text{]}`
  * dynamic viscosity:    :math:`\mu = 1.002 \times 10^{-3} \: \text{[Pa} \cdot \text{s]}` 
  * gravitational acceleration:    :math:`g = 9.807 \: \text{[m/s}^2\text{]}`
  * storativity:    :math:`S=7.5 \times 10^{-4} \: \text{[-]}`

    * derived from:    :math:`S=S_s b`, where :math:`S_s=7.5 \times 10^{-5} \: \text{[m}^{-1} \text{]}` and :math:`b=10 \: \text{[m]}`
  * porosity:    :math:`\phi = 0.2`
  * transmissivity:    :math:`T=4.7 \times 10^{-4} \: \text{[m}^2\text{/s]}`

    * derived from:    :math:`T=Kb`, where :math:`K=\frac{k \rho g}{\mu}`
    * intrinsic permeability:    :math:`k = 2.35 \times 10^{-11} \: \text{[m}^2\text{]}` 


Results and Comparison
----------------------

.. plot:: amanzi_theis_isotropic_1d.py
	  :align: center


.. include:: table_values_theis.txt


References
----------

.. bibliography:: /bib/ascem.bib
   :filter: docname in docnames
   :style:  alpha
   :keyprefix: theis-


.. _about_theis:

About
-----

* Directory:  testing/verification/flow/saturated/transient/theis_isotropic_1d

* Authors:  Dylan Harp, Alec Thomas, David Moulton 

* Maintainer:  David Moulton, moulton@lanl.gov

* Input Files:

  * amanzi_theis_isotropic_1d-u.xml

    * Spec Version 2.3, unstructured mesh framework
    * mesh:  generated internally

* Analytic solution computed with model_theis_isotropic_1d.py

