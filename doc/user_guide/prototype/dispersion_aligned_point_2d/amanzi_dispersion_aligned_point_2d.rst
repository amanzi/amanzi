2D Dispersive Transport in a Uniform Flow Field
===============================================

Introduction
~~~~~~~~~~~~
This test problem is motivated by confined flow applications with transport 
from a relatively small source.  
A typical source would be a well that is fully screened over the confined flow interval.  
In this case, the steady-state flow field is uniform and the mass source rate 
does not affect the flow field.  
The dominant transport processes are advection and dispersion.  
Dispersion results in the spreading of the released mass in all directions 
including upstream of the source location.   
Longitudinal dispersivity (in the axis of flow) is typically several 
times larger than the transverse dispersivity.  
Accuracy of the simulated concentrations in the vicinity of the source 
is sensitive to the resolution of the grid and the numerical transport scheme.  

This two-dimensional transport problem for a constant rate mass loading point 
source in a steady-state flow field tests the Amanzi advection and dispersion process kernels.  
Options tested include specified velocity field, constant mass loading rate 
source term, advection, longitudinal and transverse dispersivities, and nonuniform gridding.
The analytical solution uses Green's functions to integrate solutions for point and 
line sources in each of the principle coordinate directions to generate 
advective-dispersive transport from a point, line, planar, or rectangular sources.  
The verification data used in this test is generated from AT123D, a program of 
these analytical solutions for one-, two-, or three-dimensional transport of 
heat, dissolved chemicals, or dissolved radionuclides in a homogeneous 
aquifer subject to a uniform, stationary regional flow field (Yeh, 1981). 
The assumption is that the modeled aquifer is infinite in the direction of 
flow but can be considered finite or infinite in the vertical and 
transverse horizontal directions. 
The test problem specification is described in (Aleman, 2007) as Problem 5.2.  
The verification data from AT123D that is used to compare against Amanzi 
output is listed in Tables 5.2.3, 5.2.4, and 5.2.5 of that document.

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
  + \phi \tau D_m.

The off-diagonal entries are

.. math::
  D_{xy} = D_{yx} 
  = (\alpha_L - \alpha_T) \frac{v_x\, v_y}{\| \boldsymbol{v}\|}.


Schematic
~~~~~~~~~

Note, the values in the schematic correlate to the values found in
:ref:`Plot-Table`.

.. figure:: schematic/schematic.png 
    :figclass: align-center

    **a) Plume from point source loading in a constant uniform groundwater flow field. b) Problem domain.**
                    
.. _Variables:
        
Defining Variables
~~~~~~~~~~~~~~~~~~

* *Q* constant pumping rate
* :math:`h(r,0)` initial water table table height
* *r* radial distnace measured outward from well
* *t* duration of pumping time


The Amanzi grid resolution is nominally 15 m x 15 m but with 0.2 m x 0.2 m source cell.  This results in 84 grid cells in the x-direction and 37 grid cells in the y-direction.  The internal Amanzi calculational units are kilograms, meters, and seconds.  The nominal parameters specified in Aleman (2007) and the Amanzi conversions are presented below.

