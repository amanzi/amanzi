1D Ion exchange
===============

Overview
--------

This test example performs the simulation of cation exchange on a single exchange site in a 1D flow domain. 

Features tested
~~~~~~~~~~~~~~~

* 1D flow
* 1D advective transport 
* Geochemical reactions

	* Ion exchange (equilibrium)

About
~~~~~

* Test case ID: 1SSConTran-ion-exchange
* Test type: Benchmark
* Benchmark simulator: PFlotran 
* Files:

  * Amanzi input file: amanzi-u-1d-ion-exchange.xml
  * Benchmark simulator input file: 1d-ion-exchange.in

* Location: amanzi/examples/examples/phase2/chemistry/1d-ion-exchange
* Author: B. Andre, G. Hammond
* Testing and Documentation: S. Molins
* Last tested on Aug 31 2014

Introduction
------------

Ion exchange description...
The simulation is run to 50 years.

Model
-----

Flow and transport 
~~~~~~~~~~~~~~~~~~

See the :doc:`../1d-tracer/amanzi_u-1d-tracer` example.

Primary species
~~~~~~~~~~~~~~~

Four primary species are used: 3 cations (
:math:`Na^+`
,
:math:`Ca^{2+}`
,
:math:`Mg^{2+}`
)
and 1 anion (
:math:`Cl^-`
) to charge-balance the solution.

Ion exchange 
~~~~~~~~~~~~

:math:`Na^+`
,
:math:`Ca^{2+}`
,
:math:`Mg^{2+}`
exchange on the single bulk site
:math:`X-`
with a cation exchange capacity (CEC) of 750.0 :math:`eq/m^3`. Note that the CEC is a material property, defined in 'Material Properties'.

The three equilbrium exchange reactions are (by convention, secondary species are given in the left hand side, while primary species are in the right hand side):

* :math:`NaX = Na^+ + X^-\;log(K)=1.0`
* :math:`CaX_2 = Ca^{2+} + 2 X^-\;log(K)=0.2953`
* :math:`MgX_2 = Mg^{2+} + 2 X^-\;log(K)=0.1666`

Problem specifications
----------------------

See the :doc:`../1d-tracer/amanzi_u-1d-tracer` example.

Results and Comparison
----------------------

Expected results
~~~~~~~~~~~~~~~~

These are the expected results.

Simulation results
~~~~~~~~~~~~~~~~~~

Here go the figure and table.

