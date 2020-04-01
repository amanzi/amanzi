.. _amanzi_dispersion_45_point_2d:

2D Dispersive Transport in a Uniform Flow Field (diagonal)
==========================================================

Capabilities Tested
-------------------

This two-dimensional transport problem --- with a constant rate mass
loading point source in a steady-state flow field --- tests the Amanzi
flow process kernel as well as advection and dispersion features of the transport process kernel. It is a slight modification
of the previous test :ref:`amanzi_dispersion_aligned_point_2d`: the direction 
of the velocity field is rotated by :math:`45^\circ` counter-clockwise with respect to the 
mesh. Thus, the test simulates non-grid-aligned advective and dispersive transport. 
Capabilities tested include:
  
  * single-phase, one-dimensional flow
  * two-dimensional transport
  * steady-state saturated flow
  * constant-rate solute mass injection well 
  * advective transport
  * dispersive transport
  * longitudinal and transverse dispersivity
  * non-reactive (conservative) solute
  * homogeneous porous medium
  * isotropic porous medium
  * statically refined (non-uniform) mesh
  
For details on this test, see :ref:`about_diagonal_dispersion`.


Background
----------

Numerical dispersion caused by flow and transport processes not aligned with the 
rectangular grid can cause numerical solutions to be innacurate. Such flow and transport
regimes can arise if different solutes are injected at different rates, or if the model domain
is characterized by geologic media with non-uniform anisotropy, to name a couple examples. It
is therefore desirable to minimize this numerical dispersion in numerical flow and transport
schemes.

This test problem is a slight modification of the previous test :ref:`amanzi_dispersion_aligned_point_2d`: the direction of the velocity field is rotated by :math:`45^\circ` 
counter-clockwise with respect to the mesh. 


Model
-----

The analytical solution addresses the advection-dispersion equation 
(saturation :math:`s_l = 1`):

.. math::
  \phi \frac{\partial C}{\partial t} 
  = Q + {\nabla \cdot }(\phi \boldsymbol{D} \nabla C),

where :math:`\phi` is porosity, :math:`C` is the solute concentration [kg/m\ :sup:`3`], :math:`t` is time [s], :math:`Q` is solute mass injection rate [kg of tracer/m\ :sup:`3`/s], and :math:`\boldsymbol{D}` is the dispersion tensor:

.. math::
  \boldsymbol{D} = \begin{pmatrix}
  D_{xx} & D_{xy} \\[0.5ex]
  D_{yx} & D_{yy}
  \end{pmatrix}.

Let :math:`\boldsymbol{v} = (v_x,\,v_y)` denote the pore velocity.
Then, the diagonal entries in the dispersion tensor are

.. math::
  D_{xx} = \alpha_L \frac{v_x^2}{\| \boldsymbol{v}\|}
  + \alpha_T \frac{v_y^2}{\| \boldsymbol{v}\|},
  \qquad
  D_{yy} = \alpha_L \frac{v_y^2}{\| \boldsymbol{v}\|}
  + \alpha_L \frac{v_x^2}{\| \boldsymbol{v}\|},

where :math:`\alpha_L` is longitudinal dispersivity [m] and :math:`\alpha_T` is transverse dispersivity [m].
The off-diagonal entries are:

.. math::
  D_{xy} = D_{yx} 
  = (\alpha_L - \alpha_T) \frac{v_x\, v_y}{\| \boldsymbol{v}\|}.


Problem Specification
---------------------


Schematic
~~~~~~~~~

The flow direction and source location are shown in the following schematic.

.. figure:: schematic/schematic.png 
    :figclass: align-center
    :width: 600 px

    **Plume from point source loading in a constant uniform groundwater flow field.**


Mesh
~~~~

The background mesh consists of square cells with size :math:`H=15` m.
It has 83 grid cells in the x-direction and 83 grid cells in the y-direction. 
The mesh is gradually refined toward the source such that the well is
represented by a square cell of :math:`h=0.46875` [m] (:math:`h = H/32`).
The mesh refinement adds 8.4% more cells.

.. figure:: figures/mesh.png 
    :figclass: align-center
    
    **Computational mesh with 7468 cells and static refinement around the well.**


Variables
~~~~~~~~~

* initial concentration condition:    :math:`C(x,0)=0 \: \text{[kg/m}^3\text{]}`
* constant solute mass injection rate:    :math:`Q=8.1483 \times 10^{-8} \: \text{[kg/m}^3\text{/s]}`
* volumetric fluid flux:    :math:`\boldsymbol{q}=(1.3176 \times 10^{-6},\,1.3176 \times 10^{-6}) \: \text{[m}^3 \text{/(m}^2 \text{ s)]}` or :math:`\text{[m/s]}`

* Material properties:

  * isotropic hydraulic conductivity:    :math:`K = 84.41 \: \text{[m/d]}`
    derived from :math:`K=\frac{k \rho g}{\mu}`, where permeability :math:`k = 1.0 \times 10^{-10} \text{ [m}^2\text{]}`
  
  * porosity:    :math:`\phi=0.35` 
  * longitudinal dispersivity:    :math:`\alpha_L=21.3 \: \text{[m]}` 
  * transverse dispersivity:    :math:`\alpha_T=4.3 \: \text{[m]}` 
  * fluid density:    :math:`\rho = 998.2 \: \text{[kg/m}^3\text{]}`
  * dynamic viscosity:    :math:`\mu = 1.002 \times 10^{-3} \: \text{[Pa} \cdot \text{s]}` 
  * gravitational acceleration:    :math:`g = 9.807 \: \text{[m/s}^2\text{]}` 

* total simulation time:    :math:`t=1400 \: \text{[d]}`

.. Boundary conditions: :math:`C(x,t)=0 \: \text{[kg/m}^3\text{]}`


Results and Comparison
----------------------

The plume structure is characterized by three line cuts.
The first cut is given by line :math:`y=x` that goes through the well.
The two other cuts are given by tranverse cuts :math:`y=-x` and :math:`y=\sqrt{848}-x`.

.. plot:: verification/transport/dispersion_45_point_2d/amanzi_dispersion_45_point_2d-a.py
   :align: center

.. include:: table_centerline.txt

The analytic data were computed with the AT123DAT software package.
A difference observed near downstream boundary for the Amanzi with the 
first-order transport scheme (Amanzi(1st)).
The second-order transport scheme provides excellent match.

.. plot:: verification/transport/dispersion_45_point_2d/amanzi_dispersion_45_point_2d-b.py
   :align: center

.. include:: table_cross-section-b.txt


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

  * amanzi_dispersion_45_point_2d-u.xml, Spec Version 2.3
  * mesh amanzi_dispersion_45_point_2d.exo

.. * Mesh Files:

  .. * amanzi_dispersion_45_point_2d.exo

* Analytic solution computed with AT123D-AT

  * Subdirectory: at123d-at

  * Input Files: 

    * at123d-at_centerline.list, at123d-at_centerline.in
    * at123d-at_slice_x=0.in

.. todo:: 

  * write script to generate the tranverse cut at distance 420 m from well [jpo]
  * need \*.list file for slice x=0? [jpo] 
