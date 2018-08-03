Transient Drawdown due to Pumping Wells in a Bounded Domain
===========================================================

Capabilities Tested
-------------------


Background
----------

As we see from the comparison of results from Amanzi and Butler's solutions, it is inevitable that the numerical solutions from the Amanzi model will not match the results from these analytical solutions at later time, because these analytical solutions were derived for the unbounded domains while in numerical simulations the domain is always bounded.  In this example, we will verify the Amanzi model using the solution for the bounded domain with uniform hydraulic properties.


Problem Specification
---------------------

Transient flow in saturated uniform porous media can be represented by

.. math:: \frac{\partial ^2 h}{\partial x^2} 
   + \frac{\partial ^2 h}{\partial y^2} 
   + \sum_{i=1}^{N_w} \frac{Q_i}{T} \delta(x-x_i)\delta(y-y_i) H(t-t_i)
   = \frac{S}{T} \frac{\partial h}{\partial t}
  :label: flow

where 
:math:`h` [L] is the head,
:math:`t` [T]is the time,
:math:`T` [L\ :sup:`2`\/T] is the transmissivity, 
:math:`S` [-] is the storage coefficient,
:math:`Q_i` [L\ :sup:`3`\/T]is the pumpage at a  well located at :math:`(x_i,y_i)` that starts pumping at :math:`t_i`,
:math:`N_w` is the number of pumping wells,
:math:`\delta(x)` is the Direc delta function, being 1 for :math:`x = 0` and 0 otherwise, and
:math:`H(x)` is the Heaviside function, being 1 for :math:`x \ge 0` and 0 otherwise.

Initially the head is a constant everywhere in the domain:

.. math:: h(x,y,0) = H_0.
  :label: ic_bounded_domain_2D

The boundary conditions are:

.. math:: 
.. math::    h({\bf x}, t) =  0, \text{   for } {\bf x} \in \Gamma_D\\
.. math::    T \nabla h({\bf x,t}) \cdot {\bf n}({\bf x})  = q({\bf x}, t)  \text{  for } {\bf x} \in \Gamma_N\\
  :label: bc_bounded_domain_2D

where :math:`\Gamma_D` is the Dirichlet boundary and :math:`\Gamma_D` is the Neumann boundary.

The drawdown solution can be written as

.. math:: s(x,y,t) = -\frac{4}{DT} \sum_{n=0}^{\infty} \sum_{m=1}^{\infty}
   a_n \sin(\alpha_m x) \cos(\beta_n y) 
   \sum_{i =1}^{N_w} \frac{Q_i \sin(\alpha_n x_i) \cos(\beta_n y_i)} {\alpha_m^2 +\beta_n^2}
   H(t-t_i) \left [1 - e^{-\frac{T}{S}(\alpha_m^2 + \beta_n^2)(t-t_i)} \;\right ]
  :label: solution

where :math:`\alpha_m = m \pi/L_1, m=1,2,\cdots`, 
:math:`\beta_n = n \pi/L_2, n=0,2,\cdots`, 
:math:`L_1` and :math:`L_2` are the domain size in the x and y directions, respectively,
:math:`D = L_1L_2` is the area of the domain,
:math:`a_0 =1/2`, and :math:`a_n =1` for :math:`n \ge 1`,


Problem Specification
---------------------


Schematic
~~~~~~~~~

The domain configuration and well locations are indicated in the following schematic.

.. figure:: schematic/config_bounded.jpg
    :figclass: align-center
    :width: 600 px

    **Schematic of verification problem for bounded domains.**

    
Mesh
~~~~

The model domain is 2400 m :math:`\times` 2400 m. It has 3600 grid cells: 600 cells in the x-direction, 600 cells in y-direction, and 1 cell in the z-direction. 


Variables
~~~~~~~~~

* Transmissivity: :math:`\;\; T = 0.011574 \; m^2/s` 
* Storativity: :math:`\;\; S = 2\times 10^{-4}` 
* Pumping rate: :math:`\;\; Q = 1000 \;m^3/day \:(= 0.011574 \;m^3/s)`
* Pumping well location (1200 m, 1200 m) and pumping starts at :math:`t = 0`

Observation well locations (1224 m, 1200 m) and (1300 m, 1200 m), so  their distance  to the pumping well is 24 m and 100 m, respectively.

The boundary conditions are given as: constant pressure head of 1.07785 MPa (i.e., 100 m) at the left and the right  boundaries, and a no-flow condition was imposed on the upper and lower boundaries. Initially the pressure head is 1.07785 MPa everywhere in the domain.


Results and Comparison
----------------------

.. _Plot_BoundedDomain2D:


Comparison of  Analytic Solution and Amanzi Results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. plot:: amanzi_boundedDomain_2d.py
   :align: center

The comparison shows that the results from the Amanzi model are nearly identical to those from the analytical solution.
Detailed comparison can be found in :cite:`bd-Lu_Harp_Birdsell_benchmarking_2014`.


References
----------

.. bibliography:: /bib/ascem.bib
   :filter: docname in docnames
   :style:  alpha
   :keyprefix: bd-


.. _about_bounded_domain:

About
-----

* Directory: testing/verification/flow/saturated/transient/boundedDomain

* Authors:  Zhiming Lu (zhiming@lanl.gov),  Dylan Harp (dharp@lanl.gov)

* Maintainer(s):  Zhiming Lu,  Dylan Harp

* Input Files: 
  
  * amanzi_boundedDomain_2d.xml
 
    * Spec: Version 2.0
    * mesh: Generated in running time
    * runs

* Analytical Solutions

  * Directory: analytic/

  * Executable: boundedDomain.x, compiled from FORTRAN code under Linux environment.

  * Input Files:

    * input

  * Output Files:
   
    * test_h_tr.dat,  drawdown as a function of time for all observation wells


Status
~~~~~~


