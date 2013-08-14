1D Tritium first-order decay
============================

Overview
--------

This test example performs the simulation of decay of tritium, a radioactive isotope of hydrogen, in a 1D flow domain. 

Capabilities tested
~~~~~~~~~~~~~~~~~~~

* 1D flow
* 1D advective transport 
* Geochemical reactions

	* First order decay (aqueous kinetics)

About
~~~~~

* Test case ID: 1SSConTran-tritium
* Test type: Benchmark
* Benchmark simulator: PFlotran
* Files

  * Amanzi input file: amanzi-u-1d-tritium.xml
  * Benchmark simulator input file: 1d-tritium.in

* Location: amanzi/examples/examples/phase2/chemistry/1d-tritium
* Author: B. Andre, G. Hammond
* Testing and Documentation: S. Molins
* Last tested on Aug 31, 2013

Introduction
------------

Tritium is a radioactive isotope of hydrogen that was added to the atmosphere during the nuclear bomb testing during the 1950's and 1960's at well above background concentrations. After the end of widespread nuclear bomb testing, tritium concentration dropped. Concentrations of tritium in groundwater reflect when the water was in contact with the atmosphere, and they are used as an estimate of the rate of groundwater recharge. Tritium concentrations decay with time according to a first order expression. This first order decay is typically characterized by the half life (
:math:`t_{1/2}`
). In this example, tritium injected at constant concentration at the left boundary decays as it travels down gradient a 1D flow domain. The simulation is run to 50 years.

Model
-----

Flow and transport 
~~~~~~~~~~~~~~~~~~

See the :doc:`../1d-tracer/amanzi_u-1d-tracer` example.

First order (kinetic) decay
~~~~~~~~~~~~~~~~~~~~~~~~~~~

One single reaction is considered: tritium decay. This reaction is homogeneous and is modeled with a first order kinetic rate expression:

:math:`r = - \lambda c` 

where 
:math:`c`
is the concentration of tritium, and 
:math:`\lambda`
is the first order constant. This is related to the half life of tritium according to:

:math:`\lambda = \frac{ln(2)}{t_{1/2}}`

In this example, the half life of tririum is taken as 
:math:`t_{1/2} = 13.31 s`
. Thus, 
:math:`\lambda = 1.78577 \cdot 10^{-9} s^{-1}`

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

