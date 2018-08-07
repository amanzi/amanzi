Transient Flow to a Pumping Well in an Anisotropic Aquifer (Hantush)
====================================================================

Capabilities Tested
-------------------

This two-dimensional transient flow problem in a non-leaky aquifer with constant pumping tests the Amanzi saturated flow process kernel in anistropic media.
Capabilities tested include:

  * transient, single-phase flow
  * saturated flow conditions 
  * drawdown prediction in anistropic media
  * statically refined (non-uniform) mesh

For details on this test, see :ref:`about_hantush_anisotropic`.


Background
----------

An aquifer is homogeneous and anisotropic when the hydraulic conductivity is a function of
direction only. Because of this, the transmissivity of such an aquifer varies with direction.
The transmissivity in the major direction of anisotropy may be 2 - 10 times greater than 
in the minor direction of anisotripy. Estimates of groundwater flow may be in significant error
if the effect of such anisotropy is not accounted for. The in-place determination of the
hydraulic properties of anisotropic aquifers, then, is of great practical importance. 

Hantush and Thomas :cite:`han-Hantush_1966_method` developed an analytical solution for measuring
drawdown during constant discharge for a completely penetrating well in a
homogenous anisotropic, nonleaky confined aquifer of infinite extent. Their work follows 
Hantush's previously published work :cite:`han-Hantush_1966_wells`, extending it by describing
a method for obtaining the transmissivity in any direction and storativity of a homogeneous
anisotropic non-leaky aquifer by delineating the expected elliptical shape of an equal residual
drawdown curve. 


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
problem simplifies to the Theis solution :cite:`han-Theis1935relation`.  The variables in
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

The mesh consists of 12,208 cells. There is a single cell in the z-direction, which is uniform :math:`\Delta z=5.0` m everywhere.

.. figure:: figures/mesh.png
    :figclass: align-center

    **Unstructured computational mesh with 12208 cells.**


Variables
~~~~~~~~~

* :math:`Q=2.0` constant pumping rate, [m\ :sup:`3`/s]
* :math:`\phi=0.3` constant porosity
* :math:`S_s=7.5 \times 10^{-5}` specific storage (:math:`S=S_s b`), [m\ :sup:`-1`]
* :math:`k_x=2.3543 \times 10^{-11}, \: k_y=k_z=2.3543 \times 10^{-12}` permeability tensor (:math:`k_i \rho g /\mu=K_i=T_i/b)`, [m\ :sup:`2`] 
* :math:`\rho=998.20` water density, [kg/m\ :sup:`3`]
* :math:`g=9.81` gravity constant, [m/s\ :sup:`2`]
* :math:`\mu=1.002 \times 10^{-3}` dynamic viscosity, [kg/s/m]
* :math:`t=86400` total simulation time (= 1 day), [s]

Initial condition: zero drawdown everywhere in domain

Boundary conditions: zero drawdown on four lateral boundaries


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

	    
.. _about_hantush_anisotropic:

About
-----

* Directory: testing/verification/flow/saturated/transient/hantush_anisotropic_2d

* Authors: Alec Thomas, Konstantin Lipnikov

* Maintainer: David Moulton (moulton@lanl.gov)

* Input Files:

  * amanzi_hantush_anisotropic_2d-u.xml

    * Spec Version 2.3, unstructured mesh framework
    * mesh:  porflow4_6.exo

* Mesh Files:

  * porflow4_6.exo


.. todo::

  * Documentation:

    * Decide whether to keep structured run
    * Include info about analytic solution calculation?
    * convert units in Variables to be same as in Model?
