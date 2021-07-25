.. raw:: latex
	 
   \clearpage

1D Ion exchange
===============

Overview and Capabilities tested
--------------------------------

This test example performs the simulation of cation exchange on a single exchange site in a 1D flow domain, testing the following capabilities:

* 1D flow
* 1D advective transport 
* Geochemical reactions

	* Ion exchange (equilibrium)

For details on this test, see :ref:`about_ion_exchange`.

Background
----------

Ion exchange description to come.

Model
-----

Flow and transport 
~~~~~~~~~~~~~~~~~~

See the :doc:`../tracer_1d/amanzi_u-1d-tracer` example.

Primary species
~~~~~~~~~~~~~~~

Four primary species are used: 3 cations (
:math:`\ce{Na^+}`
,
:math:`\ce{Ca^{2+}}`
,
:math:`\ce{Mg^{2+}}`
)
and 1 anion (
:math:`\ce{Cl^-}`
) to charge-balance the solution.

Ion exchange 
~~~~~~~~~~~~

:math:`\ce{Na^+}`
,
:math:`\ce{Ca^{2+}}`
,
:math:`\ce{Mg^{2+}}`
exchange on the single bulk site
:math:`\ce{X^-}`
. The three equilbrium exchange reactions are (by convention, secondary species are given in the left hand side, while primary species are in the right hand side):

* :math:`\ce{NaX} = \ce{Na^+} + \ce{X^-}`,
  :math:`\text{ } \log(K)=1.0`
* :math:`\ce{CaX_2} = \ce{Ca^{2+}} + 2 \ce{X^-}`,
  :math:`\text{ } \log(K)=0.2953`
* :math:`\ce{MgX_2} = \ce{Mg^{2+}} + 2 \ce{X^-}`,
  :math:`\text{ } \log(K)=0.1666`

Problem specifications
----------------------

Flow and transport 
~~~~~~~~~~~~~~~~~~

See the :doc:`../tracer_1d/amanzi_u-1d-tracer` example. The simulation is run to 50 years.

Geochemistry 
~~~~~~~~~~~~

The cation exchange capacity (CEC) of the bulk site is 750.0 :math:`\text{ eq/m}^3`.

Results and Comparison
----------------------

Expected results
~~~~~~~~~~~~~~~~

Expected results to come.

Simulation results
~~~~~~~~~~~~~~~~~~

The figure below shows both the aqueous and sorbed concentrations of :math:`\ce{Na^+}, \ce{Ca^+}, \ce{Mg^{2+}}, \ce{Cl^-}` along the flow direction at 50 years.

.. plot:: benchmarking/chemistry/ion_exchange_1d/ion_exchange_1d.py

..   :align: left

.. _about_ion_exchange:

About
-----

* Benchmark simulator: PFlotran 
* Files:

  * Amanzi input file/s (native chemistry):  amanzi-u-1d-ion-exchange.xml
  * Amanzi input file/s (Alquimia chemistry): amanzi-1d-ion-exchange-alq.xml, 1d-ion-exchange.in, ion-exchange.dat 
  * Benchmark simulator input file: 1d-ion-exchange.in, ion-exchange.dat

* Location: testing/benchmarking/chemistry/ion_exchange_1d
* Author: B. Andre, G. Hammond
* Testing and Documentation: S. Molins
* Last tested on Oct 4, 2013
