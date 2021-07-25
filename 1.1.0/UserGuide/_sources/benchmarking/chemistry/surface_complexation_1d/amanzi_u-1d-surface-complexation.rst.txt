.. raw:: latex
	 
   \clearpage
   
1D Surface complexation
=======================

Overview and Capabilities tested
--------------------------------

This test example performs the simulation of 
:math:`\ce{Zn^{2+}}`
surface complexation on weak and strong iron hydroxide sites in a 1D flow domain. 

* 1D flow
* 1D advective transport 
* Geochemical reactions

	* Aqueous complexation reactions (equilibrium)
	* Surface complexation reactions (equilibrium)

For details on this test, see :ref:`about_surface_complexation`.

Background
----------

Introduction
------------

Description of surface complexation... The simulation is run to 50 years.

Model
-----

Flow and transport 
~~~~~~~~~~~~~~~~~~

See the :doc:`../tracer_1d/amanzi_u-1d-tracer` example.

Aqueous complexation
~~~~~~~~~~~~~~~~~~~~

Five equilibrium reactions are considered in aqueous phase that include the aqueous speciation of 
:math:`\ce{Zn^{2+}}` (by convention, secondary species are given in the left hand side, while primary species are in the right hand side):

* :math:`\ce{OH^- <=> H_2O - H^+}`,
  :math:`\; log(K)=13.9951`
* :math:`\ce{Zn(OH)_{2(aq)} <=>  2 H_2O - 2 H^+ + Zn^{2+}}`,
  :math:`\; log(K)=17.3282`
* :math:`\ce{Zn(OH)_3^- <=> 3 H_2O - 3 H^+ + Zn^{2+}}`,
  :math:`\; log(K)=28.8369`
* :math:`\ce{Zn(OH)_4^{2-} <=>  4 H_2O  - 4 H^+ + Zn^{2+}}`,
  :math:`\; log(K)=41.6052`
* :math:`\ce{ZnOH^+ = H_2O -1 H^+ + 1 Zn^{2+}}`,
  :math:`\; log(K)=8.96`

Surface complexation
~~~~~~~~~~~~~~~~~~~~

Two distinct iron hydroxide (:math:`\ce{FeOH}`) surface sites are available for sorption of 
:math:`\ce{Zn^{2+}}`: weak :math:`\ce{({>}FeOH_{w})}` and strong :math:`\ce{({>}FeOH_{s})}`.
Surface site concentrations for :math:`\ce{>FeOH_{w}}` and :math:`\ce{>FeOH_{s}}` are
defined in 'Material Properties'. No explicit definition of the mineral phase on which 
:math:`\ce{Zn^{2+}}` is required. 

In total, six surface complexation reactions are considered:

* :math:`\ce{{>}FeOH_{2w}^+ <=> {>}FeOH_{w} + H^+}`,
  :math:`\; log(K)=-7.18`
* :math:`\ce{{>}FeO^{-}_{w} <=> {>}FeOH_{w} - H^+}`,
  :math:`\; log(K)=8.82`
* :math:`\ce{{>}FeOHZn^{+}_{w} <=> {>}FeOH_{w} - H^+ + Zn^{2+}}`,
  :math:`\; log(K)=2.32`
* :math:`\ce{{>}FeOH_{2s}^+ <=> {>}FeOH_{s} + H^+}`,
  :math:`\; log(K)=-7.18`
* :math:`\ce{{>}FeO^{-}_{s} <=> {>}FeOH_{s} - H^+}`,
  :math:`\; log(K)=8.82`
* :math:`\ce{{>}FeOHZn^{+}_{s} <=> {>}FeOH_{s} - H^+}`
  :math:`\; log(K)=-0.66`

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

.. plot:: benchmarking/chemistry/surface_complexation_1d/surface_complexation_1d.py


.. _about_surface_complexation:	  
	  
About
-----

* Benchmark simulator: PFlotran

* Files:

  * Amanzi input file: amanzi-u-1d-surface-complexation.xml
  * Benchmark simulator input file: 1d-surface-complexation.in

* Location: amanzi/examples/examples/phase2/chemistry/1d-surface-complexation
* Author: B. Andre, G. Hammond
* Testing and Documentation: S. Molins
* Last tested on Aug 31 2013	
