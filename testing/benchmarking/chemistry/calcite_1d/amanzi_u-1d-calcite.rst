1D Calcite dissolution
1D Calcite dissolution
======================

Overview
--------

This test example performs the simulation of calcite dissolution in a 1D flow domain. 

Capabilities tested
~~~~~~~~~~~~~~~~~~~

* 1D flow
* 1D advective transport 
* Geochemical reactions

	* Aqueous complexation reactions (equilibrium)
	* Mineral dissolution

About
~~~~~

* Test case ID: 1SSConTran-calcite
* Test type: Benchmark
* Benchmark simulator: PFlotran
* Files

  * Amanzi input file/s (native chemistry): amanzi-1d-calcite.xml, calcite.bgd
  * Amanzi input file/s (Alquimia chemistry): amanzi-1d-calcite-alq.xml, 1d-calcite.in, calcite.dat 
  * Benchmark simulator input and output file/s: pflotran/1d-calcite.in, pflotran/calcite.dat, pflotran/1d-calcite.h5

* Location: testing/benchmark/chemistry/calcite_1d/
* Author: B. Andre, G. Hammond
* Testing and Documentation: S. Molins
* Last tested on: Sept 30, 2013
	
Introduction
------------

Carbonate minerals are present in many subsurface environments and contribute to their buffering capacity. Calcite dissolution is represented here with a kinetic rate expression based on the transition state theory. In this test example, a solution under saturated with calcite is injected at x=0 into a 100-m porous domain containing calcite; as a result, calcite dissolves raising the pH and the concentration of Ca in the effluent end of the domain. The simulation is run to 50 years.

Model
-----

Flow and transport 
~~~~~~~~~~~~~~~~~~

See the :doc:`../tracer_1d/amanzi_u-1d-tracer` example.

Aqueous complexation
~~~~~~~~~~~~~~~~~~~~

Six reactions equilibrium reactions are considered in aqueous phase (by convention, secondary species are given in the left hand side, while primary species are in the right hand side):

 * :math:`\ce{OH^- <=> H2O - H^+}`,  
   :math:`\text{ } \log(K)=13.9951`              
 * :math:`\ce{CO3^{--} <=>  -H^+ + HCO3^-}`, 
   :math:`\text{ } \log(K)=10.3288`
 * :math:`\ce{CO2_{(aq)} <=> - H2O + H^+ + HCO3^-}`, 
   :math:`\text{ } \log(K)=-6.3447`
 * :math:`\ce{CaOH_+ <=> H2O - H^+ + Ca^{++}}`, 
   :math:`\text{ } log(K)=12.85000`, 
 * :math:`\ce{CaHCO3^+ <=> HCO3^- + Ca^{++}}`, 
   :math:`\text{ } \log(K)=-1.0467`
 * :math:`\ce{CaCO3_{(aq)} <=> - H^+ + HCO3^- + Ca^{++}}`,
   :math:`\text{ } \log(K)=7.0017`

Calcite dissolution
~~~~~~~~~~~~~~~~~~~

Calcite dissolution can be expressed as

:math:`CaCO_3(s) \rightarrow Ca^{2+} + CO_3^{2-}`

The rate expression is 

:math:`r = S \cdot k \cdot (1 - \frac{Q}{K_{sp}})`

where 
:math:`S`
is the reactive surface area, 
:math:`k`
is the intrinsic rate constant, 
:math:`Q`
is the ion activity product, 
:math:`Q`
is the ion activity product, and
:math:`K_{sp}`
is the solubility constant of calcite. 

Problem Specification
---------------------

Flow and transport 
~~~~~~~~~~~~~~~~~~

See the :doc:`../tracer_1d/amanzi_u-1d-tracer` example.

Geochemistry
~~~~~~~~~~~~

To do.

Results and Comparison
----------------------

Expected results
~~~~~~~~~~~~~~~~

These are the expected results.

Simulation results
~~~~~~~~~~~~~~~~~~

The figure shows the concentration of total calcium along the length of the column at the end of the simulation. 

.. plot:: prototype/chemistry/calcite_1d/calcite_1d.py
   :align: center


