1D Surface complexation
=======================

Overview
--------

This test example performs the simulation of 
:math:`Zn^{2+}`
surface complexation on weak and strong iron hydroxide sites in a 1D flow domain. 

Capabilities tested
~~~~~~~~~~~~~~~~~~~

* 1D flow
* 1D advective transport 
* Geochemical reactions

	* Aqueous complexation reactions (equilibrium)
	* Surface complexation reactions (equilibrium)

Information about this test problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Test case ID: 1SSConTran-surface-complexation
* Test type: Benchmark testing
* Benchmark simulator: PFlotran
* Files:

  * Amanzi input file: amanzi-u-1d-surface-complexation.xml
  * Benchmark simulator input file: 1d-surface-complexation.in

* Location: amanzi/examples/examples/phase2/chemistry/1d-surface-complexation
* Author: B. Andre, G. Hammond
* Testing and Documentation: S. Molins
* Last tested on Aug 31 2013	

Introduction
------------

Description of surface complexation... The simulation is run to 50 years.

Model
-----

Flow and transport 
~~~~~~~~~~~~~~~~~~

See the :doc:`../1d-tracer/amanzi_u-1d-tracer` example.

Aqueous complexation
~~~~~~~~~~~~~~~~~~~~

Five reactions equilibrium reactions are considered in aqueous phase that include the aqueous speciation of 
:math:`Zn^{2+}` (by convention, secondary species are given in the left hand side, while primary species are in the right hand side):

* :math:`OH^- = H_2O - H^+\;log(K)=13.9951`
* :math:`Zn(OH)_2(aq) =   2 H_2O  -2 H^+ + Zn^{2+}\;log(K)=17.3282`
* :math:`Zn(OH)_3^- =  3 H_2O  -3 H^+ + Zn^{2+}\;log(K)=28.8369`
* :math:`Zn(OH)_4^{2-} =  4 H_2O  -4 H^+ + Zn^{2+}\;log(K)=41.6052`
* :math:`ZnOH^+ =   H_2O -1 H^+ + 1 Zn^{2+}\;log(K)=8.96`

Surface complexation
~~~~~~~~~~~~~~~~~~~~

Two distinct iron hydroxide (FeOH) surface site are available for sorption of 
:math:`Zn^{2+}`
: weak (>FeOH_w) and strong (>FeOH_s). Surface site concentrations for >FeOH_w and >FeOH_w are defined in 'Material Properties'. No explicit definition of the mineral phase on which 
:math:`Zn^{2+}`
is required. 

In total, six surface complexation reactions are considered:

* :math:`>FeOH_{2w}^+ = >FeOH_w + H^+\;log(K)=-7.18`
* :math:`>FeO^-_w =   >FeOH_w - H^+\;log(K)=8.82`
* :math:`>FeOHZn^+_w = >FeOH_w - H^+ + Zn^{2+}\;log(K)=2.32`
* :math:`>FeOH_{2s}^+ = >FeOH_s + H^+\;log(K)=-7.18`
* :math:`>FeO^-_s =   >FeOH_s - H^+\;log(K)=8.82`
* :math:`>FeOHZn^+_s =  >FeOH_s - H^+\;log(K)=-0.66`

Problem Specification
---------------------

Flow and transport 
~~~~~~~~~~~~~~~~~~

See the :doc:`../1d-tracer/amanzi_u-1d-tracer` example.

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

Here go the figure and table.
