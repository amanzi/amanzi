Transient Flow in a 2D Confined Aquifer with a Linear Strip
===========================================================

Capabilities Tested
-------------------

This two-dimensional flow problem --- with a constant pumping rate in a heterogeneous confined aquifer --- tests the Amanzi saturated flow process kernel. 
Capabilities tested include:
 
  * single-phase, two-dimensional flow
  * transient flow
  * saturated flow
  * constant-rate pumping well
  * constant-head (Dirichlet) boundary conditions
  * heterogeneous porous medium
  * isotropic porous medium
  * uniform mesh

For details on this test, see :ref:`about_butler_strip_2d`.


Background
----------

Butler and Liu (1991) :cite:`strip-Butler_Liu_linear_strip_1991` developed a semi-analytical solution for calculating drawdown in an aquifer system, in which an infinite linear strip of one material is embedded in a matrix of different hydraulic properties. The problem is interested in the drawdown as a function of location and time due to pumping from a fully penetrating well located in any of three zones. The problem is solved analytically in the Fourier-Laplace space and the drawdown is solved numerically by inversion from the Fourier-Laplace space to the real space.

The solution reveals several interesting features of flow in this configuration dependening on the relative contrast in material transmissivity. If the transmissivity of the strip is much higher than that of the matrix, linear and bilinear flow regimes dominate during the pumping test. If the contrast between matrix and strip properties is not as extreme, radial flow dominates.


Model
-----

Flow within zones that do not contain the pumping well can be described mathematically in terms of drawdown :math:`s` as :eq:`flow_nowell`

.. math:: \frac{\partial ^2 s_i}{\partial x^2} 
   + \frac{\partial ^2 s_i}{\partial y^2} 
   = \frac{S_i}{T_i} \frac{\partial s_i}{\partial t},
  :label: flow_nowell

where 
:math:`s_i` is the drawdown [L] in material :math:`i`,
:math:`t` is time [T],
:math:`T_i` is the transmissivity [L\ :sup:`2`\/T] of material :math:`i`, and
:math:`S_i` is the storage coefficient [-] of material :math:`i`.

Flow within zones that contain the pumping well can be represented as

.. math:: \frac{\partial ^2 s_i}{\partial x^2} 
   + \frac{\partial ^2 s_i}{\partial y^2} 
   + \frac{Q}{T_i} \delta(x-a)\delta(y-b)
   = \frac{S_i}{T_i} \frac{\partial s_i}{\partial t},
  :label: flow_well

where
:math:`Q` is the pumping rate [L\ :sup:`3`\/T] from well located at :math:`(a,b)`,
:math:`\delta(x)` is the Direc delta function, being 1 for :math:`x = 0` and :math:`0 \text{ otherwise}`.

The initial conditions are the same for all three zones:

.. math:: s_i(x,y,0) = 0.
  :label: ic_ButlerLiu_strip

The boundary conditions are:

.. math:: 
.. math::    s_1(-\infty, y, t) =  0,\\
.. math::      s_i(x,\pm\infty, t) =  0, \\
.. math::     s_3(\infty, y, t) =  0,\\
.. math::     s_1(-d, y, t) =  s_2(-d, y, t),\\
.. math::     s_2(0, y, t) =  s_3(0, y, t),\\
.. math::      T_1\frac{\partial s_1(-d,y,t)}{\partial x} = T_2\frac{\partial s_2(-d,y,t)}{\partial x},\\
.. math::      T_2\frac{\partial s_2(0,y,t)}{\partial x} = T_3\frac{\partial s_3(0,y,t)}{\partial x}.
  :label: bc_ButlerLiu_strip


Problem Specification
---------------------


Schematic
~~~~~~~~~

The domain configuration and well locations are indicated in the following schematic. The origin of the coordinate system is shown in the figure as 'o'.

.. figure:: schematic/butler_strip_schematic.png
    :figclass: align-center
    :width: 600 px

    **Schematic of the Butler and Liu's Linear Strip verification problem.**


Mesh
~~~~

The background mesh is :math:`1202 \: m \times 1202 \: m \times 1 \: m` and consists of 361,201 cells. There are 601 cells in the x-direction, 601 cells in the y-direction, and 1 cell in the z-direction.  


Variables
~~~~~~~~~

* Domain:

  * :math:`x_{min} = y_{min} = 1202, z_{min} = 0 \text{ [m]}`
  * :math:`x_{max} = y_{max} = 1202, z_{max} = 1 \text{ [m]}`
  * aquifer thickness:    :math:`b=z_{max}-z_{min} = 1 \text{ [m]}`
  * Zone 1 (left zone):   
    
    * :math:`-1202 \leq x \leq -10`
    * :math:`-1202 \leq y \leq 1202`
    * :math:`0 \leq z \leq 1`

  * Zone 2 (strip):   
    
    * :math:`-10 \leq x \leq 10`
    * :math:`-1202 \leq y \leq 1202`
    * :math:`0 \leq z \leq 1`
  
  * Zone 3 (right zone):   
    
    * :math:`10 \leq x \leq 1202`
    * :math:`-1202 \leq y \leq 1202`
    * :math:`0 \leq z \leq 1`

  * pumping well location:    :math:`(a,b) = (0,0) \text{ [m]}`
  * observation well locations:

    * :math:`(x_{obs24},y_{obs24},z_{obs24}) = (24,0,1) \text{ [m]}`
    * :math:`(x_{obs100},y_{obs100},z_{obs100}) = (100,0,1) \text{ [m]}`

* Boundary and initial conditions:

  * initial condition:    :math:`s(x,y,0)=0 \text{ [m]}`
  * constant-head (Dirichlet) boundary conditions:    :math:`s(x_{min,max},y_{min,max},t) = 0 \text{ [m]}`
  * well-head pumping rate:    :math:`Q = -11.5485 \text{ [m}^3 \text{/s]} = 1000 \text{ [m}^3 \text{/d]}`
  * duration of pumping:    :math:`t_{max} = 31.7 \text{ [days]}`

* Material properties:

  * transmissivity (all isotropic):

    * :math:`T_1 = 0.11574 \text{ [m}^2 \text{/s]}`
    * :math:`T_2 = 0.011574 \text{ [m}^2 \text{/s]}`
    * :math:`T_3 = 0.0011574 \text{ [m}^2 \text{/s]}`
    
      * derived from:    :math:`T=Kb`, where :math:`K=\frac{k \rho g}{\mu}`

      * intrinsic permeability:    :math:`k_1 = 1.187 \times 10^{-8}, k_2 = 1.187 \times 10^{-9}, k_2 = 1.187 \times 10^{-10} \text{ [m}^2 \text{]}`

  * storativity:   
    
    * :math:`S_1=5.0 \times 10^{-3} \: \text{[-]}`
    * :math:`S_2=2.0 \times 10^{-3} \: \text{[-]}`
    * :math:`S_3=2.0 \times 10^{-4} \: \text{[-]}`

      * derived from:    :math:`S=S_s b`, where :math:`b=10 \: \text{[m]}`

  * porosity:    :math:`\phi_{1,2,3} = 0.25`

.. * Width of the strip: :math:`\;\; d = 18 \;m`

.. * Pumping well location :math:`\;\; (-9\; m, 0\; m)`

.. The boundary conditions are given as: constant pressure of 1.07785 MPa (i.e., head = 100 m) at all four boundaries and initially the pressure is 1.07785 MPa (head = 100 m) everywhere in the domain. 

.. Observation well locations :math:`(15\; m, 0\; m)` and :math:`(91\; m, 0\; m)`, which gives the distance between the pumping well and observation wells :math:`r = 24 \;m` and :math:`r = 100 \;m`.


Results and Comparison
----------------------

.. _plot_ButlerLiu_strip:

.. plot:: amanzi_butler_strip_2d.py
             :align: center


The comparison shows that the results from the Amanzi model match the analytical solution very well at early time, and that they deviate when the effect of pumping hits the constant head boundary of the domain. Note that the analytical solution was developed for unbounded domain, so it is therefore expected that the two solutions will deviate from each other at late time.  To show that such a deviation is indeed caused by the boundary effect, we also conducted numerical simulations using 
FEHM, a widely used numerical simulator for simulating heat and mass flow in subsurface environment :cite:`strip-Zyvoloski_FEHM_summary_1997`. It is showed that the results from Amanzi are almost the same as those from FEHM, see :cite:`strip-Lu_Harp_Birdsell_benchmarking_2014` for detailed comparison.


References
----------

.. bibliography:: /bib/ascem.bib
   :filter: docname in docnames
   :style:  alpha
   :keyprefix: strip-


.. _about_butler_strip_2d:

About
-----

* Directory: testing/verification/flow/saturated/transient/butler_strip_2d

* Authors:  Zhiming Lu (zhiming@lanl.gov),  Dylan Harp (dharp@lanl.gov)

* Maintainer(s):  Zhiming Lu,  Dylan Harp

* Input Files: 
  
  * amanzi_butler_strip_2d-u.xml

    * Spec: Version 2.3, unstructured mesh framework
    * Mesh: generated internally 

* Analytical Solutions

  * Directory: analytic/

  * Executable: butler_strip.x, compiled from FORTRAN code under the Linux environment.

  * Input Files:

    * now.dat

  * Output Files:

    * drdn.dat,  drawdown as a function of time for all observation wells.


Status
~~~~~~

The analytical solution was solved using a FORTRAN code modified from the original code from Greg Ruskauf.
We may need to implement the algorithm by ourselves or get permission from Greg Ruskauf for using the code.
As the flow problem was solved analytically in the Laplace-Fourier transformed space, one needs to implement
numerical inversion from the Laplace-Fourier transformed space back to the real space.

