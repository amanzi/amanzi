2D Dispersive Transport in a Uniform Flow Field (aligned)
=========================================================

Capabilities Tested
-------------------

This two-dimensional transport problem --- with a constant rate mass
loading point source in a steady-state flow field --- tests the Amanzi flow process kernel as well as advection and dispersion features of the 
transport process kernel.  
Capabilities tested include:
  
  * single-phase, one-dimensional flow
  * two-dimensional transport
  * steady-state flow
  * saturated flow
  * constant-rate solute mass injection well 
  * advective transport
  * dispersive transport
  * longitudinal and transverse dispersivity
  * non-reactive (conservative) solute
  * homogeneous porous medium
  * isotropic porous medium
  * statically refined (non-uniform) mesh

For details on this test, see :ref:`about_aligned_dispersion`.


Background
----------

This test problem is motivated by confined flow applications with
transport from a relatively small source.  A typical source would be a
well that is fully screened over the confined flow interval.  In this
case, the steady-state flow field is uniform and the injection well tracer 
mass source rate
does not affect the flow field.  The dominant transport processes are
advection and dispersion.  Dispersion results in the spreading of the
released mass in all directions, including upstream of the source
location.  Longitudinal dispersivity (in the axis of flow) is
typically up to an order of magnitude larger than the transverse dispersivity.
Accuracy of the simulated concentrations in the vicinity of the source
is sensitive to the resolution of the grid and the numerical transport
scheme.

The analytical solution uses Green's
functions to integrate solutions for point and line sources in each of
the principle coordinate directions to generate advective-dispersive
transport from a point, line, planar, or rectangular sources.  The
verification data used in this test is generated from AT123D, a
program of these analytical solutions for one-, two-, or
three-dimensional transport of heat, dissolved chemicals, or dissolved
radionuclides in a homogeneous aquifer subject to a uniform,
stationary regional flow field :cite:`da-Yeh_AT123D_1981`.  The assumption is that the
modeled aquifer is infinite in the direction of flow but can be
considered finite or infinite in the vertical and transverse
horizontal directions.  The test problem specification is described in
Problem 5.2 of :cite:`da-Aleman_PORFLOW_2007`.  The verification data from AT123D that
is used to compare against Amanzi output is listed in Tables 5.2.3,
5.2.4, and 5.2.5 of that document.


Model
-----

The analytical solution addresses the advection-dispersion equation 
(saturation :math:`s_l = 1`):

.. math::
  \phi \frac{\partial C}{\partial t} 
  = Q + {\nabla \cdot }(\boldsymbol{D} \nabla C),

where :math:`\phi` is porosity, :math:`C` is the solute concentration [kg/m\ :sup:`3`], :math:`t` is time [s], :math:`Q` is solute mass injection rate [kg of tracer/m\ :sup:`3`/s], and :math:`\boldsymbol{D}` is the dispersion tensor:

.. math::
  \boldsymbol{D} = \begin{pmatrix}
  D_{xx} & D_{xy} \\[0.5ex]
  D_{yx} & D_{yy}
  \end{pmatrix}.

Let :math:`\boldsymbol{v} = (v_x,\,v_y)` denote the pore velocity,
:math:`\tau` the tortuosityi [-], and :math:`D_m` the molecular diffusion [m\ :sup:`2`/s].
Then the diagonal entries in the dispersion tensor are

.. math::
  D_{xx} = \alpha_L \frac{v_x^2}{\| \boldsymbol{v}\|}
  + \alpha_T \frac{v_y^2}{\| \boldsymbol{v}\|}
  + \phi \tau D_m, 
  \qquad
  D_{yy} = \alpha_L \frac{v_y^2}{\| \boldsymbol{v}\|}
  + \alpha_L \frac{v_x^2}{\| \boldsymbol{v}\|}
  + \phi \tau D_m,

where :math:`\alpha_L` is longitudinal dispersivity [m] and :math:`\alpha_T` is transverse dispersivity [m]. The off-diagonal entries are:

.. math::
  D_{xy} = D_{yx} 
  = (\alpha_L - \alpha_T) \frac{v_x\, v_y}{\| \boldsymbol{v}\|}.


Problem Specification
---------------------


Schematic
~~~~~~~~~

The domain, flow direction, and source location are shown in the following schematic.

.. figure:: schematic/schematic.png 
    :figclass: align-center
    :width: 600 px

    **a) Plume from point source loading in a constant uniform groundwater flow field. b) Problem domain.**
                    

Mesh
~~~~

The background mesh consists of square cells with size :math:`H=15` m.
It has 83 grid cells in the x-direction and 37 grid cells in the y-direction. 
The mesh is gradually refined toward the source such that the well is
represented by a square cell of :math:`h=0.46875` [m] (:math:`h = H/32`).
The mesh refinement adds 17% more cells.

.. figure:: figures/mesh.png 
    :figclass: align-center

    **Computational mesh with 3650 cells and static refinement around the well.**


Variables
~~~~~~~~~

* constant solute mass injection rate:    :math:`Q=8.1483 \times 10^{-8} \: \text{[kg/m}^3\text{/s]}`
* volumetric fluid flux:    :math:`\boldsymbol{q}=(1.8634 \times 10^{-6},\,0.0) \: \text{[m}^3 \text{/(m}^2 \text{ s)]}` or :math:`\text{[m/s]}`

* Material properties:

  * isotropic hydraulic conductivity:    :math:`K = 84.41 \: \text{[m/d]}`

    * derived from:    :math:`K=\frac{k \rho g}{\mu}`, where permeability, :math:`k = 1.0 \times 10^{-10} \text{ [m}^2\text{]}`
  
  * porosity:    :math:`\phi=0.35` 
  * longitudinal dispersivity:    :math:`\alpha_L=21.3 \: \text{[m]}` 
  * transverse dispersivity:    :math:`\alpha_T=4.3 \: \text{[m]}` 
  * molecular diffusion coefficient:    :math:`D_m=0.0 \: \text{[m}^2\text{/s]}` 
  * fluid density:    :math:`\rho = 998.2 \: \text{[kg/m}^3\text{]}`
  * dynamic viscosity:    :math:`\mu = 1.002 \times 10^{-3} \: \text{[Pa} \cdot \text{s]}` 
  * gravitational acceleration:    :math:`g = 9.807 \: \text{[m/s}^2\text{]}` 

* total simulation time:    :math:`t=1400 \: \text{[d]}`

Initial concentration condition: :math:`C(x,0)=0 \: \text{[kg/m}^3\text{]}`

.. Boundary conditions: :math:`C(x,t)=0 \: \text{[kg/m}^3\text{]}`


Results and Comparison
----------------------

The plume structure is characterized by three line cuts.
The first cut is given by line :math:`y=0` that goes through the well.
The two other cuts are given by lines :math:`x=0` and :math:`x=424`.

.. plot:: amanzi_dispersion_aligned_point_2d-a.py
   :align: center

.. include:: table_centerline.txt

The analytic data were computed with the AT123DAT software package.
A difference is observed near the downstream boundary for Amanzi with 
the first-order transport scheme (boxes), while the second-order transport 
scheme provides excellent match (circles).

.. plot:: amanzi_dispersion_aligned_point_2d-b.py
   :align: center

.. include:: table_cross-section-b.txt

.. plot:: amanzi_dispersion_aligned_point_2d-c.py
   :align: center

.. include:: table_cross-section-c.txt


References
----------

.. bibliography:: /bib/ascem.bib
   :filter: docname in docnames
   :style:  alpha
   :keyprefix: da-

	    
.. _about_aligned_dispersion:

About
-----

* Directory: testing/verification/transport/saturated/steady-state/dispersion_aligned_point_2d

* Authors:  Konstantin Lipnikov, Marc Day, David Moulton, Steve Yabusaki 

* Maintainer(s): Konstantin Lipnikov, Marc Day 

* Input Files:

  * amanzi_dispersion_aligned_point_2d-u.xml 

    * Spec Version 2.2, unstructured mesh framework
    * mesh:  amanzi_dispersion_aligned_point_2d.exo
 
  * amanzi_dispersion_aligned_point_2d-s.xml

    * Spec Version 1.2, structured AMR framework

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


Status
~~~~~~

  * Input Files:

    * Version 2.2 - unstructured: runs 2D problem, results are in excellent agreement
    * Version 1.2 - structured AMR: runs

  * Documentation:

    * Complete for unstructured mesh framework, including line plots and tables.

.. todo:: 

  * Documentation:

    * Decide whether to run this as a 2D or 3D problem
    * Do we need a short discussion on numerical methods (i.e., discretization, splitting, solvers)?
    * Store *Gold Standard* simulation results (need name and location)?
    * Add plots for structured AMR results to these plots or make subplots
    * What about colormap plot of the field?
    * Schematics should include citation
    * Schematics should be redone with tikz.
    * Citations may/should be improved with the natbib extension for sphinx?
    * Tables need to be centered and captions added (more sphinx tuning/extensions)
    * Results from Structured AMR mesh infrastructure runs need to be added.
    * Should we separate authors/maintainers into code and documentation?
