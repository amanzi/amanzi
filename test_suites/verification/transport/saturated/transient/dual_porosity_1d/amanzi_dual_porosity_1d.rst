1D Dispersive Transport with Dual Porosity Model
================================================

Capabilties Tested
------------------

This one-dimensional transport problem, with a constant contaminant 
boundary condition and a steady-state flow field, tests the Amanzi
transport dual porosity model.  
Capabilties tested include,
  
  * prescribed constant velocity field 
  * constant boundary condition for tracer
  * advection in an isotropic medium, longitudinal dispersivity
  * generalized dual porosity

For details on this test, see :ref:`about_dual_porosity`.


Background
----------

This test problem is motivated by confined flow applications with
transport from a relatively small source. A typical source would be a
well that is fully screened over the confined flow interval. In this
case, the steady-state flow field is uniform and the tracer influx 
does not affect the flow field. The dominant transport processes are
advection, dispersion in fractured rock, and solute exchange with 
immobile fluid in matrix. Matrix play the role of storage and delays
propagation of tracer down the flow stream.  Dispersion results in the 
spreading of the released mass in all directions including upstream 
of the source location. 
Accuracy of the break-through curve at a distant point down the stream
is sensitive to the time step in Amanzi's operator split approach.

Model
-----

The analytical solution addresses the advection-dispersion equation in the fracture

.. math::
  \frac{\partial (\phi_f\, C_f)}{\partial t} 
  =
  - \boldsymbol{\nabla} \cdot (\boldsymbol{q} C_f) 
  + \boldsymbol{\nabla} \cdot (\phi_f\, (\boldsymbol{D} + \tau_f M) \boldsymbol{\nabla} C_f) 
  - \frac{\phi_m\,\tau_m}{L_m}\, M \nabla C_m,

and the linear ODE in the matrix:

.. math::
  \frac{\partial (\phi_m\, C_m)}{\partial t} = -\nabla\cdot (\phi_m\, \tau_m\, M \nabla C_m).

Here
:math:`\phi` is porosity,
:math:`\boldsymbol{q}` is the Darcy velocity,
:math:`\boldsymbol{D}` is dispersion tensor, 
:math:`M` is molecular diffusion coefficient, and
:math:`L_m` is the characteristic matrix depth defined typically as the ratio of a matrix block
volume to its surface area.
Subscripts :math:`f` and :math:`m` indicate fracture and matrix media, respectively. 

We neglect transverse dispersion, so that the dispersion tensor is diagonal:

.. math::
  \boldsymbol{D} = \begin{pmatrix}
  D_{xx} & 0 \\[0.5ex]
  0      & D_{yy}
  \end{pmatrix},
  \qquad
  D_{xx} = \alpha_L v_x, \quad D_{yy} = 0.
  

Problem Specification
---------------------

Schematic
~~~~~~~~~

The domain, flow direction, and source location are shown in the schematic below.
Note that :math:`B = L_m` here.

.. figure:: schematic/schematic.pdf
    :figclass: align-center
    :width: 400 px

    *Periodic structure of fractures and matrix.*
                    

Mesh
~~~~

The background mesh consists of square cells with size :math:`H=1` m.
It has 100 grid cells in the x-direction and 1 grid cell in the y-direction. 


Variables
~~~~~~~~~

* :math:`\boldsymbol{q}=(7.922 \cdot 10^{-9},\,0.0)` constant Darcy velocity [m/s]
* :math:`\phi_f=0.0001` constant fracture porosity
* :math:`\phi_m=0.25` constant matrix porosity
* :math:`\alpha_L=1.0` longitudinal dispersivity [m]
* :math:`\alpha_T=0.0` transverse dispersivity [m]
* :math:`M=2.65 \cdot 10^{-9}` molecular diffusion coefficient [m^2/s]
* :math:`\tau_m = 1.0` matrix tortuosity [-]
* :math:`\tau_f = 1.0` fracture tortuosity [-]
* :math:`L_m = 0.8` characteristic matrix depth [m]
* :math:`T=200` simulation time [y]

Initial condition: :math:`C_m(x,0)=C_f(x,0) = 0` [kg/m^3]

Boundary conditions: :math:`C_f(x,t)=1` [kg/m^3] at :math:`x=0.0` of fracture.


Results and Comparison
----------------------

The concentrantaion at the fracture end point as the function of time
(results are in a good agreement with the analytic solution):

.. plot:: amanzi_dual_porosity_1d.py
   :align: center

The verification data used in this test is generated from a Fortran code
shared by the authors of :cite:`Sudicky_et_al_1982`. 


References
----------

.. bibliography:: /bib/ascem.bib
   :filter: docname in docnames
   :style:  alpha
   :keyprefix: da-

	    
.. _about_dual_porosity:

About
-----

* Directory: testing/verification/transport/saturated/transient/dual_porosity_1d

* Authors:  Konstantin Lipnikov, David Moulton

* Maintainer(s): Konstantin Lipnikov

* Input Files:

  * amanzi_dual_porosity_1d-u.xml 

  * Spec Version 2.3, unstructured mesh framework
 

* Analytic solution computed with Sudicky's code:

  * Subdirectory: sudicky

  * Input Files: 

    * tracer_conc.txt


.. todo:: 

  * Add references
  * Complete for unstructured mesh framework, including line plots. Tables will be added.
  * Decide whether to add simple dual porosity model with Warren-Root parameter.
  * Do we need a short discussion on numerical methods (i.e., discretization, splitting, solvers)?
