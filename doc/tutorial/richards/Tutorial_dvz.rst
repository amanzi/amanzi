Transient vadose zone flow and Tc-99 transport
----------------------------------------------

This tutorial illustrates *Amanzi* simulation of transient vadose zone flow and 
non-reactive (tracer) transport in the context of Tc-99 migration at the DOE 
Hanford BC Cribs and Trenches site. 
The tutorial involves a two-dimensional heterogeneous three-layer subsurface
system with time- and space-varying infiltration at the ground surface. 

BC Cribs and Trenches site
~~~~~~~~~~~~~~~~~~~~~~~~~~

The two-dimensional geometry and salient features of the Hanford BC Cribs and 
Trenches site, as represented in this simplified tutorial, are defined in the 
following schematic diagram:

	.. image:: Schematic.png
		:scale: 50 %
		:align: center

General recharge at the site (:math:`U`) increased following development of the Hanford
facility in 1956, and Cribs B-17 and B-18 represent localized sources of water 
infiltration and Tc-99 contamination (:math:`C_{Tc-99}`) to the subsurface:

+------------------------+-------------------+------------------------------------------+
|  ID                    | :math:`U` (mm/yr) |   :math:`C_{Tc-99}` (mol/m\ :sup:`3`\ )  |
+========================+===================+==========================================+
| General site, pre-1956 | 3.5               |  N/A                                     |
+------------------------+-------------------+------------------------------------------+
| General site, post-1956| 47                |  N/A                                     |
+------------------------+-------------------+------------------------------------------+
| B-17, Jan 1956         |  8025             |  1.88e-6                                 |
+------------------------+-------------------+------------------------------------------+
| B-18, Feb-Mar 1956     |  10439            |  2.27e-6                                 |
+------------------------+-------------------+------------------------------------------+

The three facies underlying the site have the following material properties:

+-------------------------------------------------+---------------+---------------+---------------+
| Material property                               | Facies 1      | Facies 2      | Facies 3      |
+=================================================+===============+===============+===============+
| Porosity                                        | 0.4082        | 0.2206        | 0.2340        |
+-------------------------------------------------+---------------+---------------+---------------+
| Particle density (kg/m\ :sup:`3`\ )             | 2720.0        | 2720.0        | 2720.0        |
+-------------------------------------------------+---------------+---------------+---------------+
| Horizontal permeability (m\ :sup:`2`\ )         | 1.9976e-12    | 6.9365e-11    | 2.0706e-09    |
+-------------------------------------------------+---------------+---------------+---------------+
| Vertical permeability (m\ :sup:`2`\ )           | 1.9976e-13    | 6.9365e-12    | 2.0706e-09    |
+-------------------------------------------------+---------------+---------------+---------------+
| Capillary pressure model                        | van Genuchten | van Genuchten | Brooks-Corey  |
+-------------------------------------------------+---------------+---------------+---------------+
| Relative permeability model                     | Mualem        | Mualem        | Mualem        |
+-------------------------------------------------+---------------+---------------+---------------+
| :math:`\alpha` (Pa\ :sup:`-1`\ )                | 1.9467e-04    | 2.0260e-03    | 2.0674e-03    |
+-------------------------------------------------+---------------+---------------+---------------+
| :math:`S_r`                                     | 0.0           | 0.0           | 0.0           |
+-------------------------------------------------+---------------+---------------+---------------+
| van Genuchten :math:`m`                         | 0.2294        | 0.2136        | N/A           |
+-------------------------------------------------+---------------+---------------+---------------+
| Brooks-Corey :math:`\lambda`                    | N/A           | N/A           | 0.3006        |
+-------------------------------------------------+---------------+---------------+---------------+


*Amanzi* input specifications
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A working *Amanzi* input file for this tutorial simulation is listed in a 
separate section below. 
*Amanzi* input file structure and XML element names are intended to be reasonably 
self-explanatory; however, selective commentary is provided here for further guidance.
Currently *Amanzi* requires the SI units indicated in the ``model_description/units`` 
element. The corresponding units of density, viscosity, pressure, 
hydraulic conductivity, and permeability  are thus 
kg/m\ :sup:`3`\ , Pa :math:`\cdot` s, Pa, m/s, and m\ :sup:`2`\ , respectively. 
Simulation outputs also use these dimensional units. 
The use of named constants in the  ``definitions`` element, e.g., 

``<constant name="pre_1956_recharge" type="area_mass_flux" value="1.1071E-7"/>``

allows the user to use a fixed name throughout the input file, but change its 
numerical value in only one location. Here the time- and space-varying infiltration 
parameters are input as defined constants. ``time_macro name="Macro 1"`` defines the
output times for visualization in seconds; for example, the last time is
``6.33834396E10`` seconds or year 2008.5 for a 365.25 day year. 
``Macro 2`` is defined but not used in this
example. Under ``process_kernels`` the chemistry kernel is turned off because the
simulation involves only non-reactive transport. In the ``phases``
element, equation-of-state (eos) computations are turned off in favor of direct 
specification of fluid density and viscosity. Under the ``mesh`` element
a three-dimensional grid is specified, however, only one cell is specified
for the :math:`y`-coordinate (``ny = "1"``) making the simulation effectively
2D in the :math:`(x,z)` coordinates. The grid spacing is 0.5 meters 
horizontally and 0.42 meters in the vertical dimension. In the 
``initial_conditions`` element, a hydrostatic profile is specified through a
linear pressure profile with slope/gradient = :math:`- \rho g` = -9789 Pa/m. 
The reference point for :math:`P = P_{atm}` = 101325 Pa is set to 
:math:`z` = 0.5 meters rather than the water table elevation of 
:math:`z` = 0.0 meters to be compatible with a prior benchmarking comparison to
a STOMP code model. The boundary conditions
at the top surface are defined in the problem statement as a volumetric fluxes 
[m\ :sup:`3`\ /m\ :sup:`2`\ d = m/d]. *Amanzi* currently requires boundary 
conditions of this type to be specified as mass fluxes, 
:math:`\rho U` [kg\ :sup:`3`\ /m\ :sup:`2`\ s]. The ``inward_mass_flux`` values
in the XML input have these units.


*Amanzi* execution
~~~~~~~~~~~~~~~~~~

*Amanzi* is executed from the command-line using this or analogous command for
the user's specific installation:

``mpirun -n 4 /ascem/amanzi/install/current/bin/amanzi --xml_schema=/ascem/amanzi/install/current/bin/amanzi.xsd --xml_file=dvz_3_layer_2d-isv2.xml``

Here execution using 4 processor cores is specified by ``-n 4``. Successful completion
is marked by ``SIMULATION_SUCCESSFUL`` followed by a timing summary ::

	Amanzi::SIMULATION_SUCCESSFUL

	**********************************************************
	***                   Timing Summary                   ***
	**********************************************************



*Amanzi* results
~~~~~~~~~~~~~~~~

The figure shows saturation dynamics. Saturation grows underneath the 
left and right cribs during their operational cycles. Note it takes time 
for two moving plumes to penetrate into the middle soil which leads to saturation 
increase on the interface between soils. 

  .. image:: saturation.png
	:scale: 50 %
	:align: center


*Amanzi* XML input file
~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: dvz_3_layer_2d-isv2.xml
    :language: xml
