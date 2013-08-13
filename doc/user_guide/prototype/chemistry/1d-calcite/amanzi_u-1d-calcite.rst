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
* Benchmark simulator: PFlotran (*input file:* 1d-calcite.in)
* Files

  * Amanzi input file: amanzi-1d-calcite.xml
  * Benchmark simulator input file: 1d-calcite.in

* Location: amanzi/examples/examples/phase2/chemistry/1d-calcite
* Author: B. Andre, G. Hammond
* Testing and Documentation: S. Molins
* Last tested on: Aug 31, 2013
	
Introduction
------------

Carbonate minerals are present in many subsurface environments and contribute to their buffering capacity. Under common subsurface flow conditons, calcite dissolution is a relatively fast geochemical reaction leading often to local geochemical equilibrium and sharp dissolution fronts. Calcite dissolution is represented with a kinetic rate expression based on the transition state theory. In this test example, a solution under saturated with calcite is injected at x=0 into a 100-m porous domain containing calcite; as a result, calcite dissolves raising the pH and concentration of cations at the effluent of the domain. The simulation is run to 50 years.

Model
-----

Flow and transport 
~~~~~~~~~~~~~~~~~~

See the :doc:`../1d-tracer/amanzi_u-1d-tracer` example.

Aqueous complexation
~~~~~~~~~~~~~~~~~~~~

Six reactions equilibrium reactions are considered in aqueous phase (by convention, secondary species are given in the left hand side, while primary species are in the right hand side):

* :math:`OH^- = H_2O - H^+\;log(K)=13.9951`
* :math:`CO_3^{--} =  - H^+ + HCO_3^-\;log(K)=10.3288`
* :math:`CO_2(aq) =  - H_2O + H^+ + HCO_3^-\;log(K)=-6.3447`
* :math:`CaOH_+ = H_2O - H^+ + Ca^{++}\;log(K)=12.85000`
* :math:`CaHCO_3^+ = HCO_3^- + Ca^{++}\;log(K)=-1.0467`
* :math:`CaCO_3(aq) =  - H^+ + HCO_3^- + Ca^{++}\;log(K)=7.0017`

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

See the :doc:`../1d-tracer/amanzi_u-1d-tracer` example.

Results and Comparison
----------------------

Expected results
~~~~~~~~~~~~~~~~

These are the expected results.

Simulation results
~~~~~~~~~~~~~~~~~~

Here go the figure and table.

