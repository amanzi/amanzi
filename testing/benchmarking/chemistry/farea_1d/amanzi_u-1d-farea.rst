1D F-Area Full Geochemistry
===========================

Overview
--------

This test example performs the simulation of a complex geochemical network in a simple 1D flow domain, combining reaction types described separately in :doc:`../tritium_1d/amanzi_u-1d-tritium`, :doc:`../calcite_1d/amanzi_u-1d-calcite`, and :doc:`../ion_exchange_1d/amanzi_u-1d-ion-exchange`. The reaction network is based on the Savannah River Site F-Area geochemistry. 

..  comment out for now
    , and :doc:`../surface_complexation_1d/amanzi_u-1d-surface-complexation`


Capabilities tested
~~~~~~~~~~~~~~~~~~~

* 1D flow
* 1D advective transport 
* Geochemical reactions

	* Aqueous complexation reactions (equilibrium)
	* Radioactive decay
	* Mineral dissolution
	* Ion exchange
	* Surface complexation

About
~~~~~

* Test case ID: 1SSConTran-farea-full
* Test type: Benchmark
* Benchmark simulator: PFlotran
* Files

  * Amanzi input file: amanzi-u-1d-farea.xml
  * Benchmark simulator input file: 1d-farea.in

* Location: amanzi/examples/examples/phase2/chemistry/1d-farea
* Author: B. Andre, G. Hammond
* Testing and Documentation: S. Molins, H.M. Wainright
* Last tested on: November 1, 2013
	
Introduction
------------

This problem intends to demonstrate the ASCEM capability of modeling complex geochemical reactions and contaminant transport. The model is based on the one developed at the Savannah River Site F-Area for predicting the uranium and acidic plume migration in the groundwater [Bea2012]_. 

The SRS is located in south-central South Carolina, near Aiken, approximately 100 miles from the Atlantic Coast. The SRS F-Area seepage basins received approximately 7.1 billion litters of acidic, low-level waste solutions from 1950 through 1989 [Millings2012]_. After the basins were closed and capped in 1991, the site has been under active remediation, including the pump-and-treat and hybrid funnel-and-gate system. 

The highly acidic plume including many radionuclides (e.g., uranium isotopes, strontium-90, iodine-129, technetium, and tritium) is developed from the basins to downgradient in groundwater. The concentration of U-238 and pH are of the current main interest at the site. Understanding and predicting the plume mobility and its fate requires modeling complex geochemistry, since pH is influenced by mineral dissolution and precipitation, and the uranium mobility is greatly influenced by pH.

Model
-----

Flow and transport 
~~~~~~~~~~~~~~~~~~

See the :doc:`../tracer_1d/amanzi_u-1d-tracer` example.

Geochemistry
~~~~~~~~~~~~

The primary geochemical system consists of 13 reactive chemical components and 8 minerals [Bea2012]_. A detailed list of reactions and geochemical parameters is included in the table below. Further detail can be found in [Freshley2012]_.

The geochemical processes include:

* Equilibrium aqueous complexation
* Kinetically controlled mineral dissolution and precipitation
* Mineral adsorption/desorption (nonelectrostatic surface complexation)
* Equilibrium ion exchange

The sorption model is based on a single-site equilibrium, pH-dependent surface complexation model, which provides the principal control on the uranium migration rate. The nonelectrostatic model used here is applied to the bulk sediment rather than to specific pure mineral phases that serve as sorbents in the conventional electrostatic model. Additionally, the nonelectrostatic model does not assume the presence of well-defined mineral phases [Davis1998]_. An ion exchange model includes reactions involving the major cations (:math:`Ca^{2+}`, :math:`Na{+}`, and :math:`Al^{3+}`, along with :math:`H^+`) and provides primary pH buffering along with the mineral reactions.

Problem Specification
---------------------

Flow and transport 
~~~~~~~~~~~~~~~~~~

See the :doc:`../tracer_1d/amanzi_u-1d-tracer` example.

Geochemistry 
~~~~~~~~~~~~

Add tables.


Results and Comparison
----------------------

Expected results
~~~~~~~~~~~~~~~~

These are the expected results.

Simulation results
~~~~~~~~~~~~~~~~~~

These are the simulation results.

.. plot:: prototype/chemistry/farea_1d/farea_1d.py

..   :align: left

References
----------

.. [Bea2012] Sergio A. Bea, Haruko Wainwright, Nicolas Spycher, Boris Faybishenko, Susan S. Hubbard, Miles E. Denham (2013) Identifying key controls on the behavior of an acidic-U(VI) plume in the Savannah River Site using reactive transport modeling, Journal of Contaminant Hydrology, 151:34-54. 

.. [Davis1998] Davis, J.A., J.A. Coston, D.B. Kent, and C. C. Fuller (1998). Application of the Surface Complexation Concept to Complex Mineral Assemblages, Environmental Science & Technology, Vo. 32, No. 19, 2820-2828.

.. [Freshley2012] Freshley, M.S. Hubbard, et al. (2012) Advanced Simulation Capability for Environmental Management (ASCEM) Phase II Demonstration. ASCEM-SITE-2012-01; DOI 10.2172/1055500 `link <http://www.osti.gov/bridge/product.biblio.jsp?query_id=1&page=0&osti_id=1055500>`_ 

.. [Millings2012] M.R. Millings, B.B. Looney, M.E. Denham (2012) Geochemical modeling of F Area seepage basin composition and variability. `pdf file <http://sti.srs.gov/fulltext/SRNL-STI-2012-00269.pdf>`_

