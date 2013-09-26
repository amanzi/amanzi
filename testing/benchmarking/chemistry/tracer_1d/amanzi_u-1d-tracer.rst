1D Conservative tracer transport
================================

Overview
--------

This test example performs the simulation of advective transport of a single conservative tracer component.

Capabilities tested
~~~~~~~~~~~~~~~~~~~

* 1D flow
* 1D advective (single-component) transport 

About
~~~~~

* Test case ID: 1SSConTran-tracer
* Test type: Benchmark
* Benchmark simulator: PFlotran
* Files: 
  
  * Amanzi input file: amanzi-1d-tracer.xml
  * Benchmark simulator input file: 1d-tracer.in

* Location: amanzi/examples/examples/phase2/chemistry/1d-tracer
* Author: B. Andre, G. Hammond
* Testing and Documentation: S. Molins
* Last tested on: Aug 31, 2013

Introduction
------------

When running a reactive transport problem, it is good practice to include a non-reactive component or tracer. Results obtained for this conservative tracer can be compared to results for reactive components. This comparison can provide insights into the effects of reactions on the fate of the reactive species, e.g. retardation of species subject to sorption. The problem presented here simulates the conservative (advective) transport of a single component in a 1D domain. The flow and transport components of this test problem are used as basis to develope the following reactive transport test problems: :doc:`../1d-tritium/amanzi_u-1d-tritium`, :doc:`../1d-calcite/amanzi_u-1d-calcite`, :doc:`../1d-ion-exchange/amanzi_u-1d-ion-exchange`, :doc:`../1d-surface-complexation/amanzi_u-1d-surface-complexation`, :doc:`../1d-farea/amanzi_u-1d-farea`.

Model
-----

Flow and Transport
~~~~~~~~~~~~~~~~~~

The flow equation solved is the saturated flow equation:

  :math:`\frac{\partial}{\partial t} (\phi \rho) + \nabla(\rho \mathbf{q}) = 0`

The transport equation solved is the transient advective transport equation:

  :math:`\frac{\partial}{\partial t} (c)+ \nabla(\mathbf{q} c) = 0`

Problem Specification
---------------------

* 1D model domain length: 100 meters,  
  :math:`[x_{min},x_{max}] = [0, 100]`

* The domain is discretized with 100 cells in the x-direction, 1 meter long each. 

* Porosity is 0.25

* The flow is horizontal in the x-direction. The uniform flow velocity in x-direction is
  :math:`7.913 \cdot 10^{-6} kg/s`
  (or 
  :math:`7.927 \cdot 10{-9} m^3/s`
  ) specified as mass flux boundary conditions ("BC: Flux") at the inlet end of the domain, i.e. "Inward Mass Flux" at 
  :math:`x_{min} = 0 m`
  ). The oulet end of the domain (
  :math:`x_{max} = 100 m`
  ) is assigned a "BC: Uniform Pressure" boundary condition with *p = 0*. The flow velocity value is chosen so that it takes 1 year for water to move 1 cell of the discretization. Flow is solved to steady state before starting the transport solve ("Initialize To Steady").

* Simulation time = 50 years
 
* The concentration of the tracer is set to 
  :math:`10^{-20} mol/L`
  initially everywhere in the domain. From time 0, a solution with a tracer concentration of
  :math:`10^{-3} mol/L` 
  is injected during the entire simulation time (50 years).

Results and Comparison
---------------------- 

Expected results
~~~~~~~~~~~~~~~~

The flow solution is trivial and is independent of the choice of permeability. The flow velocity everywhere in the domain should be
:math:`7.927 \cdot 10{-9} m^3/s`

The advective front moves with the flow velocity. In the simulation time (50 years), the front moves half way through the domain. In the absence of diffusion, spreading of the front can be attributed to numerical dispersion added by the numberical scheme employed. Both Amanzi and PFlotran (the benchmark simulator used in the example) add a moderate amount of numerical dispersion to the solution. In the figure below, the solution by Amanzi at time 50 years is compared to results obtained with PFlotran (see below) along the length of the domain. 

Simulation results
~~~~~~~~~~~~~~~~~~

This figure shows the tracer front at time 50 years as a function of length along the flow path. Amanzi and PFlotran solutions are compared.

.. plot:: prototype/chemistry/tracer_1d/tracer_1d.py
   :align: center
