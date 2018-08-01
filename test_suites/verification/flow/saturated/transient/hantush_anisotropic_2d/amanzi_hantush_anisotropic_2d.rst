Transient Flow to a Pumping Well in an Anisotropic Aquifer (Hantush)
====================================================================

Capabilities Tested
-------------------

This two-dimensional transient flow problem 

Background
----------

Hantush and Thomas :cite:`han-Hantush_1966_method` developed an analytical solution for measuring drawdown
during constant discharge for a completely penetrating well in a
homogenous anisotropic nonleaky infinite confined aquifer. 

Model
-----

The soil permeabilty varies in the *x* and *y* directions; thus, anisotropic
condictions exist. Hantush and Thomas simplifies to the Theis solution :cite:`han-Theis1935relation`, which
assumes isotropic soil profile, when the permeabilty is equal in all
directions for all times. This solution will yield contour lines of
equal draw down by an elliptical shape.     

.. math::
    \frac{\partial }{\partial x} (T_x \frac{\partial h}{\partial x})+\frac{\partial }{\partial y} (T_y \frac{\partial h}{\partial y})
    = S \frac{\partial h}{\partial t} + Q \delta(x) \delta(y)

The initial conditions is

.. math::  h(x,y,0)=h_0

Due to the varying hydraulic conductivities two transmissivities are
defined as 

.. math:: T_x = K_xb \; \; and \;\; T_y=K_yb

Hantush and Thomas found the solution to the governing equation as

.. math:: s=h_0-h(r,t)=\frac{Q}{4 \pi \sqrt{T_x T_y}} \int_\phi^\infty
	  \frac{exp[-\tau]}{\tau} d\tau = \frac{Q}{4 \pi \sqrt{T_x T_y}} \; W(\phi)

where

.. math:: \phi = \frac{(x^2T_t + y^2T_x)\;S}{4T_xT_yt}

Notice when :math:`T_x=T_y`, :math:`\phi` is now equal to *u* and the
problem simplifies to the Theis solution.  The variables in
the equations above are defined in :ref:`Variables` with subtle
differences.  We have now defined transmissivity in two directions and
redefined the well function, *W*, to apply to :math:`\phi` instead of
*u*.  The integral is still solved as a negative exponential integral.  


Problem Specification
---------------------


Schematic
~~~~~~~~~

.. figure:: schematic/ellipse.png
    :figclass: align-center
    :width: 600 px

    **Schematic of an equal drawdown curve around a well in an anisotropic aquifer.**


Mesh
~~~~


Variables
~~~~~~~~~


Results and Comparison
----------------------

.. plot:: amanzi_hantush_anisotropic_2d.py
          :align: center

       
References
----------

.. bibliography:: /bib/ascem.bib
   :filter: docname in docnames
   :style:  alpha
   :keyprefix: han-

	    
.. _about_aligned_dispersion:

About
-----


