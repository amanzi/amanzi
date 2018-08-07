2D Dispersive Transport in a Uniform Flow Field (diagonal)
==========================================================

Capabilities Tested
-------------------

This two-dimensional transport problem --- with a constant rate mass
loading point source in a steady-state flow field --- tests the Amanzi
advection and dispersion features of the transport process kernel. It is a slight modification
of the previous (aligned) test: the direction of the velocity field is rotated by
:math:`45^\circ` counter-clockwise with respect to the mesh.
Capabilities tested include:
  
  * single-phase, steady-state flow induced by prescribed constant velocity field 
  * saturated flow conditions  
  * solute source with constant mass loading rate 
  * advection in an isotropic medium
  * longitudinal and transverse dispersivities
  * statically refined (non-uniform) mesh

For details on this test, see :ref:`about_diagonal_dispersion`.


Background
----------

This test problem is a slight modification of the previous test,
the direction of the velocity field is rotated by :math:`45^\circ` 
counter clockwise with respect to the mesh.


Model
-----

The analytical solution addresses the advection-dispersion equation 
(saturation :math:`s_l = 1`):

.. math::
  \phi \frac{\partial C}{\partial t} 
  = Q + {\rm div}(\boldsymbol{D} \nabla C),

where :math:`\boldsymbol{D}` is the dispersion tensor

.. math::
  \boldsymbol{D} = \begin{pmatrix}
  D_{xx} & D_{xy} \\[0.5ex]
  D_{yx} & D_{yy}
  \end{pmatrix}.

Let :math:`\boldsymbol{v} = (v_x,\,v_y)` denote the pore velocity,
:math:`\tau` the torsuosity, and :math:`D_m` the molecular diffusion.
Then the diagonal entries in the dispersion tensor are

.. math::
  D_{xx} = \alpha_L \frac{v_x^2}{\| \boldsymbol{v}\|}
  + \alpha_T \frac{v_y^2}{\| \boldsymbol{v}\|}
  + \phi \tau D_m, 
  \qquad
  D_{yy} = \alpha_L \frac{v_y^2}{\| \boldsymbol{v}\|}
  + \alpha_L \frac{v_x^2}{\| \boldsymbol{v}\|}
  + \phi \tau D_m.,

The off-diagonal entries are

.. math::
  D_{xy} = D_{yx} 
  = (\alpha_L - \alpha_T) \frac{v_x\, v_y}{\| \boldsymbol{v}\|}.


Problem Specification
---------------------


Schematic
~~~~~~~~~


Mesh
~~~~

The background mesh consists of square cells with size :math:`H=15` m.
It has 83 grid cells in the x-direction and 37 grid cells in the y-direction. 
The mesh is gradually refined toward the source such that the well is
represented by a square cell of :math:`h=0.46875` [m] (:math:`h = H/32`).
The mesh refinement adds 8.4% more cells.

.. figure:: figures/mesh.png 
    :figclass: align-center
    
    **Computational mesh with 7468 cells and static refinement around the well.**


Variables
~~~~~~~~~

* :math:`Q=8.1483 \times 10^{-8}` constant pumping rate [kg/s/m]
* :math:`\boldsymbol{q}=(1.3176 \times 10^{-6},\,0.0)` constant Darcy velocity [m/s]
* :math:`\phi=0.35` constant porosity
* :math:`k=1.0 \times 10^{-10}` intrinsic permeability [m\:sup:`2`]
* :math:`\alpha_L=21.3` longitudinal dispersivity [m]
* :math:`\alpha_T=4.3` transverse dispersivity [m]
* :math:`D_m=0.0` molecular diffusion coefficient [m\ :sup:`2`\/s]
* :math:`T=1400` simulation time [d]

Initial condition: :math:`C(x,0)=0` [kg/m\ :sup:`3`\]

Boundary conditions: :math:`C(x,t)=0` [kg/m\ :sup:`3`\]


Results and Comparison
----------------------

The plume structure is characterized by three line cuts.
The first cut is given by line :math:`y=0` that goes through the well.
The two other cuts are given by lines :math:`x=0` and :math:`x=424`.

.. plot:: verification/transport/dispersion_45_point_2d/amanzi_dispersion_45_point_2d-a.py
   :align: center

.. include:: table_centerline.txt

The analytic data were computed with the AT123DAT software package.
A difference observed near downstream boundary for the Amanzi with the 
first-order transport scheme (Amanzi(1st)).
The second-order transport scheme provides excellent match.

.. plot:: verification/transport/dispersion_45_point_2d/amanzi_dispersion_45_point_2d-b.py
   :align: center


References
----------

.. bibliography:: /bib/ascem.bib
   :filter: docname in docnames
   :style:  alpha
   :keyprefix: dd-

	    
.. _about_diagonal_dispersion:

About
-----

* Directory:  testing/verification/transport/saturated/steady-state/dispersion_45_point_2d

* Author:  Konstantin Lipnikov 

* Maintainer:  David Moulton (moulton@lanl.gov)








* Input Files:

  * amanzi_dispersion_45_point_2d.xml

    * Spec Version 2.3, unstructured mesh framework
    * mesh:  amanzi_dispersion_45_point_2d.exo

* Mesh Files:

  * amanzi_dispersion_45_point_2d.exo
