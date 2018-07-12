.. _Theis:

Transient One-Dimensional Confined Flow to a Pumping Well (Theis)
=================================================================

Capabilties Tested
------------------


Background
----------

Theis (1935) developed an analytical solution for transient (non
steady state) drawdown for a fully penetrating well by imposing the
boundary conditions: :math:`h = h_0` for :math:`t = 0` and 
:math:`h \Rightarrow h_0` as :math:`r \Rightarrow \infty`.  The equation
assumes an infinite and uniform confined aquifer.  When a well is
pumped the water table declines toward the well and flow is induced
toward the well from all directions. Theoretically, this flow can be
idealized by purely radial symmetric flow and can be decribed by the
equation below (this is analogous to heat flow by conduction developed
by Fourier)

.. math:: 
     \frac{\partial^2 h}{\partial r^2} 
   + \frac{1}{r} \frac{\partial h}{\partial r} 
   = \frac{S}{T} \frac{\partial h}{\partial t}
          

The analytical solution of drawdown as a function of time and distance
is found to be:

.. math:: s = h(r,0) - h(r,t) = \frac{Q W(u)}{4 \pi T} 
   = \int_u^\infty \frac{exp[-\tau]}{\tau} d\tau = \frac{Q W(u)}{4\pi T}

where, 

.. math:: u(r,t) = \frac{r^2 S}{4 T t}

and the integral approximation for the well function, *W(u)*, for :math:`0 < u < 1` is  

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


Model
-----

Groundwater resource evaluation is crucial to understanding the
concept of groundwater yield. This concept is vital because it helps
determine the maximum allowable pump rate out of an aquifer and aid in
the understanding of the consequnces on the water table.  For example,
pumping water out of the water table (unconfined aquifer) may dry up
near by wells due to the fall in the saturated thickness of the
aquifer. This can be prevented by modeling the fall of the water table
given a certain pumping rate.

To apply the analysis developed by Theis, three parameters must be known:

* *T*, transmissivity
* *S*, storativity
* *Q*, constant pumping rate

Transimissivity of the aquifer is defined as 

.. math:: T = Kb

where :math:`K` is the hydraulic conductivity of the aquifer and :math:`b` is the
saturated thickness.  Transmissivity values greater than 0.015 
:math:`\frac{m^2}{s}` represent aquifers capable of well exploitation.
Storativity is a dimensionless parameter that describes the amount of
water released by the aquifer per unit volume of the aquifer.
Storativity can be calculated using

.. math:: S = S_s b

where :math:`S_s` is the *specific storage* and is unique to each aquifer.
Again, :math:`b` is the saturated thickness of the aquifer.  Specific
storage represents the volume of water released per unit volume of the
aquifer per unit decline in hydraulic head.

Lastly, the constant pumping rate, :math:`Q`, is the volume of water
discharged from the well per unit time.

.. _plot_table_Theis:


Problem Specification
---------------------


Schematic
~~~~~~~~~

Note, the values in the schematic correlate to the values found in
:ref:`plot_table_Theis`.

.. figure:: schematic/Theis.png 
    :figclass: align-center

    **Figure 1.2:  Illustration of transient drawdown**
		    
.. _Variables:

Mesh
~~~~


Variables
~~~~~~~~~

* *Q* constant pumping rate
* :math:`h(r,0)` initial water table table height
* *T* transmissivity 
* *W(u)* well function
* *r* radial distnace measured outward from well
* *S* storage coefficient 
* *t* duration of pumping time

Results and Comparison
----------------------

.. plot:: amanzi_theis_isotropic_1d.py
	  :align: center


.. include:: table_values_theis.txt


References
----------


About
-----

* Directory: testing/verification/transport/saturated/steady-state/dispersion_aligned_point_2d

* Authors: Dylan Harp, Alec Thomas, David Moulton 

* Maintainer(s): David Moulton (moulton@lanl.gov)

* Input Files:

  * 

    * Spec Version 2.2, unstructured mesh framework
    * mesh:  
    * runs
 
  * amanzi_dispersion_aligned_point_2d-s.xml

    * Spec Version 1.2, structured AMR framework
    * runs

* Mesh Files:

  * amanzi_dispersion_aligned_point_2d.exo

    * two-dimensional statically refined mesh
    * treated as an unstructured polygonal mesh

  * amanzi_dispersion_aligned_point_2d-1layer.exo

    * three-dimensional statically refined mesh
    * one layer of cells in the z-direction

* Analytic solution computed with AT123D-AT

  * Subdirectory: at123d-at

  * Input Files: 

    * at123d-at_centerline.list, at123d-at_centerline.in
    * at123d-at_slice_x=0.list, at123d-at_slice_x=0.in,  
    * at123d-at_slice_x=420.list, at123d-at_slice_x=420.in

