.. raw:: latex
	 
   \clearpage

1D F-Area Full Geochemistry
===========================

Overview and Capabilities tested
--------------------------------

This test example performs the simulation of a complex geochemical network in a simple 1D flow domain, combining reaction types described separately in :doc:`../tritium_1d/amanzi_u-1d-tritium`, :doc:`../calcite_1d/amanzi_u-1d-calcite`, and :doc:`../ion_exchange_1d/amanzi_u-1d-ion-exchange`. The reaction network is based on the Savannah River Site F-Area geochemistry. This example tests the following capabilities: 

..  comment out for now
    , and :doc:`../surface_complexation_1d/amanzi_u-1d-surface-complexation`

* 1D flow
* 1D advective transport 
* Geochemical reactions

	* Aqueous complexation reactions (equilibrium)
	* Radioactive decay
	* Mineral dissolution
	* Ion exchange
	* Surface complexation

For details on this test, see :ref:`about_farea`.

Background
----------

This problem intends to demonstrate the ASCEM capability of modeling complex geochemical reactions and contaminant transport. The model is based on the one developed at the Savannah River Site F-Area for predicting the uranium and acidic plume migration in the groundwater [Bea2013]_. 

The SRS is located in south-central South Carolina, near Aiken, approximately 100 miles from the Atlantic Coast. The SRS F-Area seepage basins received approximately 7.1 billion litters of acidic, low-level waste solutions from 1950 through 1989 [Millings2012]_. After the basins were closed and capped in 1991, the site has been under active remediation, including the pump-and-treat and hybrid funnel-and-gate system. 

The highly acidic plume including many radionuclides (e.g., uranium isotopes, strontium-90, iodine-129, technetium, and tritium) is developed from the basins to downgradient in groundwater. The concentration of U-238 and pH are of the current main interest at the site. Understanding and predicting the plume mobility and its fate requires modeling complex geochemistry, since pH is influenced by mineral dissolution and precipitation, and the uranium mobility is greatly influenced by pH.

Model
-----

Flow and transport 
~~~~~~~~~~~~~~~~~~

See the :doc:`../tracer_1d/amanzi_u-1d-tracer` example.

Geochemistry
~~~~~~~~~~~~

The primary geochemical system consists of 13 reactive chemical components and 8 minerals [Bea2013]_. A detailed list of reactions and geochemical parameters is included in the tables below. Further detail can be found in [Freshley2012]_.

The geochemical processes include:

* Equilibrium aqueous complexation (Table 1)
* Mineral adsorption/desorption (nonelectrostatic surface complexation) (Table 2)
* Equilibrium ion exchange (Table 2)
* Kinetically controlled mineral dissolution and precipitation (Table 3)

The sorption model is based on a single-site equilibrium, pH-dependent surface complexation model, which provides the principal control on the uranium migration rate. The nonelectrostatic model used here is applied to the bulk sediment rather than to specific pure mineral phases that serve as sorbents in the conventional electrostatic model. Additionally, the nonelectrostatic model does not assume the presence of well-defined mineral phases [Davis1998]_. An ion exchange model includes reactions involving the major cations (:math:`Ca^{2+}`, :math:`Na{+}`, and :math:`Al^{3+}`, along with :math:`H^+`) and provides primary pH buffering along with the mineral reactions.

**Table 1. Aqueous complexes.**

+----------------------------------------------------------------------------------+------------------------------------------+
| Reaction                                                                         | logK (25ºC)                              |
+==================================================================================+==========================================+
| :math:`\ce{OH^- <=> H2O - H^+}`                                                  | 13.99                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{AlOH^{2+} <=> H2O - H^+ + Al^{3+}}`                                   |  4.96                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{Al(OH)2^+ <=> 2 H2O -2 H^+ + Al^{3+}}`                                | 10.59                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{Al(OH)3(aq) <=> 3 H2O -3 H^+ + Al^{3+}}`                              | 16.16                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{Al(OH)4^- <=> 4 H2O -4 H^+ + Al^{3+}}`                                | 22.88                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{CaOH^+ <=> H2O - H^+ + Ca^{2+}}`                                      | 12.85                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{CaHCO3^+ <=> H2O - H^+ + Ca^{2+} + CO2(aq)}`                          |  5.3                                     |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{CaCO3(aq) <=> H2O -2 H^+ + Ca^{2+} + CO2(aq)}`                        | 13.35                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{CO3^{2-} <=> H2O -2 H^+ + CO2(aq)}`                                   | 16.67                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{HCO3^- <=> H2O - H+ + CO2(aq)}`                                       |  6.34                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{NaCO3^- <=> H2O -2 H^+ + CO2(aq) + Na^+}`                             | 16.16                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{NaHCO3(aq) <=>  H2O - H^+ + CO2(aq) + Na^+}`                          |  6.19                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{NaOH(aq) <=> H2O - H^+ + Na^+}`                                       | 14.78                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{MgCO3(aq) <=> H2O -2 H^+ + CO2(aq) + Mg^{2+}}`                        | 13.69                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{MgOH^+ <=> H2O - H^+ + Mg^{2+}}`                                      | 11.79                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{MgHCO3^+ <=>  H2O - H^+ + CO2(aq) + Mg^{2+}}`                         |  5.31                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{(UO2)2(OH)2^{2+} <=> 2 H2O -2 H^+ +2 UO2^{2+}}`                       |  5.63                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{(UO2)2CO3(OH)3^- <=> 4 H2O -5 H^+ + CO2(aq) + 2 UO2^{2+}}`            | 17.57                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{(UO2)2OH^{3+} <=> H2O -1 H^+ + 2 UO2^{2+}}`                           |  2.71                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{(UO2)3(CO3)6^{6-} <=> 6 H2O -12 H^+ + 6 CO2(aq) + 3 UO2^{2+}}`        | 46.13                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{(UO2)3(OH)4^{2+} <=> 4 H2O -4 H^+ + 3 UO2^{2+}}`                      | 11.93                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{UO2(OH)4^{2-} <=> 4 H2O -4 H^+ + UO2^{2+}}`                           | 33.03                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{(UO2)3(OH)5^+ <=> 5 H2O -5 H^+ + 3 UO2^{2+}}`                         | 15.59                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{(UO2)3(OH)7^- <=> 7 H2O -7 H^+ + 3 UO2^{2}}`                          | 31.05                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{(UO2)3O(OH)2(HCO3)^+ <=> 4 H2O -5 H^+ + CO2(aq) + 3 UO2^{2+}}`        | 16.06                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{(UO2)4(OH)7^+ <=> 7 H2O -7 H^+ + 4 UO2^{2+}}`                         | 21.95                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{UO2NO3^+ <=> NO3^- + UO2^{2+}}`                                       | -0.28                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{UO2OH^+ <=> H2O -1 H^+ + 1 UO2^{2+}}`                                 |  5.21                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{UO2(OH)2(aq) <=> 2 H2O -2 H^+ + UO2^{2+}}`                            | 10.31                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{UO2(OH)3^- <=> 3 H2O -3 H^+ + UO2^{2+}}`                              | 19.22                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{UO2CO3(aq) <=> H2O -2 H^+ + CO2(aq) + UO2^{2+}}`                      |  7.01                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{UO2(CO3)2^{2-} <=> 2 H2O -4 H^+ + 2 CO2(aq) + UO2^{2+}}`              | 16.44                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{UO2(CO3)3^{4-} <=> 3 H2O -6 H^+ + 3 CO2(aq) + UO2^{2+}}`              | 28.46                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{CaUO2(CO3)3^{2-} <=> 3 H2O -6 H^+ + Ca^{2+} + 3 CO2(aq) + UO2^{2+}}`  | 22.84                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{Ca2UO2(CO3)3(aq) <=> 3 H2O -6 H^+ + 2 Ca^{2+} + 3 CO2(aq) + UO2^{2}}` | 19.32                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{MgUO2(CO3)3^{2-} <=> 3 H2O -6 H^+ + 3 CO2(aq) + Mg^{2+} + UO2^{2+}}`  | 23.91                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{UO2OSi(OH)3^+ <=> 2 H2O - H^+ + SiO2(aq) + UO2^{2+}}`                 |  2.48                                    |
+----------------------------------------------------------------------------------+------------------------------------------+

**Table 2. Surface complexation and cation-exchange reactions**

+----------------------------------------------------------------------------------+------------------------------------------+
| Reaction                                                                         | logK (25ºC)                              |
+==================================================================================+==========================================+
| Surface Complexation (*)                                                         |                                          |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{({>}SO)UO2^{+} <=> {>}SOH^{-} - H^{+} + UO2^{2+}}`                    | -0.44                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| Cation Exchange (Gaines-Thomas convention)                                       |                                          |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{Na^+X <=> 1.0 Na^+ + X^-}`                                            |  1.0                                     |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{Ca^{2+}X <=> Ca^{2+} + 2 X^-}`                                        |  0.316                                   |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{Al^{3+}X <=> Al^{3+} + 3 X^-}`                                        |  1.71                                    |
+----------------------------------------------------------------------------------+------------------------------------------+
| :math:`\ce{H^+X <=> H^+ + X^-}`                                                  |  0.025                                   |
+----------------------------------------------------------------------------------+------------------------------------------+


(*) Bulk site concentration = 0.1801 :math:`\text{moles sites m}^{-3}`

**Table 3. Mineral dissolution/precipitation reactions**

+----------------------------------------------------------------------+---------------------+----------------------------------+
| Reaction                                                             | logK (25ºC)         |  Reference                       |
+======================================================================+=====================+==================================+
| :math:`\ce{Quartz <=> SiO2(aq)}`                                     | -3.7501             | [SNL2007]_, [Guillaumont2003]_   |
+----------------------------------------------------------------------+---------------------+----------------------------------+
| :math:`\ce{Kaolinite <=> 5 H2O -6 H^+ +2 Al^{3+} + 2 SiO2(aq)}`      |  7.57               | [YangSteefel2008]_               |
+----------------------------------------------------------------------+---------------------+----------------------------------+
| :math:`\ce{Goethite <=> 2 H2O -3 H^+ + Fe^{3+}}`                     |  0.1758             |                                  |
+----------------------------------------------------------------------+---------------------+----------------------------------+
| :math:`\ce{Schoepite <=> 3 H2O -2 H^+ + UO2^{2+}}`                   |  4.8443             | [SNL2007]_, [Guillaumont2003]_   |
+----------------------------------------------------------------------+---------------------+----------------------------------+
| :math:`\ce{Gibbsite <=> 3 H2O -3 H^+ + Al^{3+}}`                     |  7.738              | [Pokrovskii-Helgeson1995]_       |
+----------------------------------------------------------------------+---------------------+----------------------------------+ 
| :math:`\ce{Jurbanite <=> 6 H2O -1 H^+ + Al^{3+} + SO4^{2-}}`         | -3.8                | [Nordstrom1982]_                 |
+----------------------------------------------------------------------+---------------------+----------------------------------+ 
| :math:`\ce{Basaluminite <=> 15 H2O -10 H^+ + 4 Al^{3+} + SO4^{2-}}`  | 22.251              | [Nordstrom1982]_                 |
+----------------------------------------------------------------------+---------------------+----------------------------------+ 
| :math:`\ce{Opal <=>  SiO2(aq)}`                                      | -3.005              | [Sonnenthal-Spycher2000]_        |
+----------------------------------------------------------------------+---------------------+----------------------------------+ 


Problem Specification
---------------------

Flow and transport 
~~~~~~~~~~~~~~~~~~

See the :doc:`../tracer_1d/amanzi_u-1d-tracer` example.

Geochemistry 
~~~~~~~~~~~~

Rate expression for mineral dissolution-precipitation reactions

.. math::
   r_j= A_j \times k_j \times a_{H^+}^p \times (1 - \Omega)
  :label: mineralrate

**Table 4. Kinetic parameters for mineral reactions**

+---------------------------------+---------------------+-----------+
| Mineral reaction                | :math:`log(k_j)`    |  p        |
+=================================+=====================+===========+
| :math:`\ce{Quartz}`             | -13.345             |  0        |
+---------------------------------+---------------------+-----------+
| :math:`\ce{Kaolinite}`          | -12.967             |  0.777    |
+---------------------------------+---------------------+-----------+
| :math:`\ce{Goethite}`           | -7.94               |  0        |
+---------------------------------+---------------------+-----------+
| :math:`\ce{Schoepite}`          |  0.301              |  0        |
+---------------------------------+---------------------+-----------+
| :math:`\ce{Gibbsite}`           | -11.5               |  0        |
+---------------------------------+---------------------+-----------+ 
| :math:`\ce{Jurbanite}`          | -8                  |  0        |
+---------------------------------+---------------------+-----------+ 
| :math:`\ce{Basaluminite}`       | -8                  |  0        |
+---------------------------------+---------------------+-----------+ 
| :math:`\ce{Opal}`               | -12.135             |  0        |
+---------------------------------+---------------------+-----------+ 

**Table 5. Chemical composition for the background (initial) and seepage (left boundary) solutions**

+-------------------------+-----------------------------+-----------------------------+-----------------------------+
| Component               | Background                  |  Seepage                    | Units                       |
+=========================+=============================+=============================+=============================+
| :math:`\ce{pH}`         |  5.4                        |   2.5                       | pH units                    |
+-------------------------+-----------------------------+-----------------------------+-----------------------------+
| :math:`\ce{Na}`         | :math:`2.78\times10^{-4}`   | :math:`3.05\times10^{-4}`   | :math:`\text{mol kgw}^{-1}` |
+-------------------------+-----------------------------+-----------------------------+-----------------------------+
| :math:`\ce{Cl}`         | :math:`9.98\times10^{-3}`   | :math:`3.39\times10^{-5}`   | :math:`\text{mol kgw}^{-1}` |
+-------------------------+-----------------------------+-----------------------------+-----------------------------+
| :math:`\ce{CO2(aq)}`    | :math:`1.23\times10^{-5}`   | :math:`1.07\times10^{-5}`   | :math:`\text{mol kgw}^{-1}` |
+-------------------------+-----------------------------+-----------------------------+-----------------------------+
| :math:`\ce{Al}`         | :math:`2.2\times10^{-8}`    | :math:`10^{-8}`             | :math:`\text{mol kgw}^{-1}` |
+-------------------------+-----------------------------+-----------------------------+-----------------------------+
| :math:`\ce{Fe(III)}`    | :math:`2.5\times10^{-16}`   | :math:`2.41\times10^{-6}`   | :math:`\text{mol kgw}^{-1}` |
+-------------------------+-----------------------------+-----------------------------+-----------------------------+
| :math:`\ce{K}`          | :math:`3.32\times10^{-5}`   | :math:`1.72\times10^{-6}`   | :math:`\text{mol kgw}^{-1}` |
+-------------------------+-----------------------------+-----------------------------+-----------------------------+
| :math:`\ce{Ca}`         | :math:`10^{-5}`             | :math:`10^{-5}`             | :math:`\text{mol kgw}^{-1}` |
+-------------------------+-----------------------------+-----------------------------+-----------------------------+
| :math:`\ce{Mg}`         | :math:`5.35\times10^{-3}`   | :math:`2.47\times10^{-6}`   | :math:`\text{mol kgw}^{-1}` |
+-------------------------+-----------------------------+-----------------------------+-----------------------------+
| :math:`\ce{U(VI)}`      | :math:`1.25\times10^{-10}`  | :math:`3.01\times10^{-5}`   | :math:`\text{mol kgw}^{-1}` |
+-------------------------+-----------------------------+-----------------------------+-----------------------------+
| :math:`\ce{Nitrate}`    | :math:`10^{-3}`             | :math:`10^{-2}`             | :math:`\text{mol kgw}^{-1}` |
+-------------------------+-----------------------------+-----------------------------+-----------------------------+
| :math:`\ce{SO4}`        | :math:`2.25\times10^{-5}`   | :math:`4.8\times10^{-5}`    | :math:`\text{mol kgw}^{-1}` |
+-------------------------+-----------------------------+-----------------------------+-----------------------------+
| :math:`\ce{SiO2(aq)}`   | :math:`1.77\times10^{-4}`   | :math:`1.18\times10^{-4}`   | :math:`\text{mol kgw}^{-1}` |
+-------------------------+-----------------------------+-----------------------------+-----------------------------+
| :math:`\ce{^3H}`        | :math:`10^{-15}`            | :math:`2.17\times10^{-9}`   | :math:`\text{mol kgw}^{-1}` |
+-------------------------+-----------------------------+-----------------------------+-----------------------------+
| :math:`\ce{pCO2(g)}`    | :math:`10^{-3.5}`           | :math:`10^{-3.5}`           | :math:`\text{mol kgw}^{-1}` |
+-------------------------+-----------------------------+-----------------------------+-----------------------------+

**Table 6. Mineral volume fractions (intial), reactive surface areas and molar volumes**

+---------------------------------+---------------------+------------------------+------------------------+
| Mineral reaction                | Volume fraction     |  Surface area          | Molar Volume           |
+                                 +---------------------+------------------------+------------------------+
|                                 | :math:`m^3/m^3`     | :math:`cm^2/cm^3`      | :math:`cm^3/mol`       |
+=================================+=====================+========================+========================+
| :math:`\ce{Quartz}`             |  0.88               |  3262.3                | 22.68                  |
+---------------------------------+---------------------+------------------------+------------------------+
| :math:`\ce{Kaolinite}`          |  0.11               |  59093.9               | 99.52                  |
+---------------------------------+---------------------+------------------------+------------------------+
| :math:`\ce{Goethite}`           |  0.016              |  11076.3               | 20.82                  |
+---------------------------------+---------------------+------------------------+------------------------+
| :math:`\ce{Schoepite}`          |  0.0                |  0.1                   | 66.08                  |
+---------------------------------+---------------------+------------------------+------------------------+
| :math:`\ce{Gibbsite}`           |  0.0                |  0.1                   | 31.95                  |
+---------------------------------+---------------------+------------------------+------------------------+
| :math:`\ce{Jurbanite}`          |  0.0                |  0.1                   | 218.93                 |
+---------------------------------+---------------------+------------------------+------------------------+ 
| :math:`\ce{Basaluminite}`       |  0.0                |  0.1                   | 29.0                   |
+---------------------------------+---------------------+------------------------+------------------------+
| :math:`\ce{Opal}`               |  0.0                |  0.1                   | 126.0                  |
+---------------------------------+---------------------+------------------------+------------------------+ 

Results and Comparison
----------------------

.. Expected results
   ~~~~~~~~~~~~~~~~
   These are the expected results.

Simulation results
~~~~~~~~~~~~~~~~~~

.. plot:: prototype/chemistry/farea_1d/farea_1d.py

..   :align: left

References
----------

.. [Bea2013] Sergio A. Bea, Haruko Wainwright, Nicolas Spycher, Boris Faybishenko, Susan S. Hubbard, Miles E. Denham (2013) Identifying key controls on the behavior of an acidic-U(VI) plume in the Savannah River Site using reactive transport modeling, Journal of Contaminant Hydrology, 151:34-54. 

.. [Davis1998] Davis, J.A., J.A. Coston, D.B. Kent, and C. C. Fuller (1998). Application of the Surface Complexation Concept to Complex Mineral Assemblages, Environmental Science & Technology, Vo. 32, No. 19, 2820-2828.

.. [Freshley2012] Freshley, M.S. Hubbard, et al. (2012) Advanced Simulation Capability for Environmental Management (ASCEM) Phase II Demonstration. ASCEM-SITE-2012-01; DOI 10.2172/1055500 `link <http://www.osti.gov/bridge/product.biblio.jsp?query_id=1&page=0&osti_id=1055500>`_ 

.. [Guillaumont2003] Guillaumont, R., T. Fanghanel, J. Fuger, I. Grenthe, V. Neck, D.A. Palmer, and M. H. Rand.  (2003)  Update on the Chemical Thermodynamics of Uranium, Neptunium, Plutonium, Americium, and Technetium.  Chemical Thermodynamics 5, OECD Nuclear Energy Agency, etd., Elsevier, Amsterdam.

.. [Pokrovskii-Helgeson1995] Pokrovskii, V.A. and Helgeson, H.C. (1995) Thermodynamic properties of aqueous species and the solubilities of minerals at high pressures and temperatures: the system Al2O3 –H2O–NaCl., Am. J. Sci. 295, 1255–1342.

.. [Millings2012] M.R. Millings, B.B. Looney, M.E. Denham (2012) Geochemical modeling of F Area seepage basin composition and variability. `pdf file <http://sti.srs.gov/fulltext/SRNL-STI-2012-00269.pdf>`_

.. [Nordstrom1982] Nordstrom, D.K. (1982) Aqueous pyrite oxidation and the consequent formation of secondary iron minerals. In Kittrick, J.A., Fanning, D.S. and Hossner, L.R. (eds.): Acid Sulfate Weathering. Soil Sci. Soc. Am. Publ.: 37—56.

.. [SNL2007] SNL (Sandia National Laboratories) (2007)  Qualification of thermodynamic data for geochemical modeling of mineral–water interactions in dilute systems (eds. Wolery T. J. and Jove Colon C. F.), Report ANL-WIS-GS-000003 REV 01. Las Vegas, Nevada: Sandia National Laboratories. ACC: DOC.20070619.0007.

.. [Sonnenthal-Spycher2000] Sonnenthal, E., and Spycher, N.  (2000) Drift-scale coupled processes model, Analysis and model report (AMR) N0120/U0110, Yucca Mountain Nuclear Waste Disposal Project, Lawrence Berkeley National Laboratory, Berkeley, California.

.. [YangSteefel2008] Yang, L. and C.I. Steefel. (2008) Kaolinite dissolution and precipitation kinetics at 22°C and pH 4. Geochimica Cosmochimica Acta 72(1), 99-116. 


.. _about_farea:

About
-----

* Benchmark simulator: PFlotran
* Files

  * Amanzi input file: amanzi-u-1d-farea.xml
  * Benchmark simulator input file: 1d-farea.in

* Location: testing/benchmarking/chemistry/farea_1d
* Author: B. Andre, G. Hammond
* Testing and Documentation: S. Molins, H.M. Wainright
* Last tested on: November 1, 2013


