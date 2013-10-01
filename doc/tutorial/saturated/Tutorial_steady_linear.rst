Tutorials
=========

One-dimensional, steady-state, saturated flow
---------------------------------------------

This tutorial uses a simple saturated flow simulation problem to demonstrate the
basic workflow required to create an input file, run *Amanzi*, and view results.

Problem definition
~~~~~~~~~~~~~~~~~~

Consider the following steady-state saturated flow scenario:

	.. image:: ProblemDefinition.png
		:scale: 35 %
		:align: center

Although the domain is three-dimensional, the imposed boundary conditions create a 
one-dimensional flow in the :math:`x`-direction. The Darcy velocity is expected
to be :math:`U = 1.688` m/d = :math:`1.954 \times 10^{-5}` m/s throughout the domain 
considering mass conservation. From Darcy's law

	.. math:: U = \frac{\rho g k}{\mu} \cdot \frac{h_{in} - h_{out}}{L} = K \frac{h_{in} - h_{out}}{L}
		:label: DarcysLaw

where :math:`K = \rho g k / \mu` is hydraulic conductivity, 
the inlet hydraulic head corresponding to the prescribed inlet flux :math:`U` is 

	.. math:: h_{in} = \frac{U L}{K} + h_{out}
		:label: inletHead

Assuming :math:`\rho = 998.2` kg/m\ :sup:`3`\  and
:math:`\mu = 1.002 \times 10^{-3}` Pa :math:`\cdot` s for water at 20C, and
:math:`g = 9.807` m/s\ :sup:`2`\ ,
the hydraulic conductivity of the porous medium is 
:math:`K = 9.7698 \times 10^{-6}` m/s = :math:`0.844` m/d
and the computed hydraulic head at the inlet should be

	.. math:: h_{in} = \frac{1.688 \text{ m/d} \cdot 100\text{ m}}{0.844 \text{ m/d}} + 120\text{ m} = 320 \text{ m}

Analogous to Equation :eq:`inletHead` the expected hydraulic head at an 
arbitary location :math:`x` is

	.. math:: h(x) = \frac{U (L - x)}{K} + h_{out}
		:label: observationHead

where :math:`(L - x)` is the distance upgradient of the outlet.

*Amanzi* input specifications
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Amanzi* input takes the form of an eXtensible Markup Language (XML) file 
with the extension ``*.xml``. Although any text editor may be used to create or 
edit the contents of the file, an editor that provides XML syntax-highlighting
improves readability and is helpful toward identifying syntax errors.
The following image is an example of XML syntax-highlighting in the *Notepad++* 
editor for Microsoft Windows (http://notepad-plus-plus.org/):

	.. image:: SyntaxHighlightingExample.png
		:scale: 50 %
		:align: center

A working *Amanzi* input file for this tutorial simulation is listed in a 
separate section below. The *Amanzi* input specification provides 
more general information on required content, options, and formatting. 
At the top level is an XML element named ``amanzi_input``

.. literalinclude:: amanzi_steady_linear-isv2.xml
    :language: xml
    :lines: 1-2,190

with an attribute ``type`` that defines the model grid type, either
``unstructured`` or ``structured``. A second attribute defines the *Amanzi* version.
Nested within this element are additional XML 
elements with associated content and attributes to define other model inputs. 
Model parameters are defined through a combination of element *attributes*, 
e.g. ``water`` in the ``liquid_phase`` element,

.. literalinclude:: amanzi_steady_linear-isv2.xml
    :language: xml
    :lines: 23

and element *content*, e.g. ``998.2`` kg/m\ :sup:`3`\  in the ``density`` element,

.. literalinclude:: amanzi_steady_linear-isv2.xml
    :language: xml
    :lines: 25

according to the *Amanzi* input specification. 
Although not required by the XML standard or *Amanzi*, 
indentation using tabs and/or spaces is commonly used to indicate
the hierarchy of elements and improve readibility. Tabs are used in this tutorial 
example to show the file hierarchy. Although *Amanzi* XML element names and the 
file structure are intended to be reasonably self-explanatory, selective commentary 
is provided here for further guidance. 

Currently *Amanzi* requires the SI units indicated in the ``model_description/units`` 
element. The units of density, viscosity, pressure, hydraulic conductivity, and 
permeability  are thus kg/m\ :sup:`3`\ , Pa :math:`\cdot` s, Pa, m/s, 
and m\ :sup:`2`\ , respectively. 
Simulation outputs also use these dimensional units. In the ``phases``
element, equation-of-state (eos) computations are turned off in favor of direct 
specification of fluid density and viscosity. 
Under ``execution_controls`` the ``start`` and ``end``
may both be zero for a steady-state calculation (``mode = "steady"``). Coordinates,
such as those specified under ``mesh/generate/box``, must be comma-separated. 
Whitespace is allowed and parenthesis may be used to bracket coordinates in the form 
``(x,y,z)`` if desired, e.g.,

.. literalinclude:: amanzi_steady_linear-isv2.xml
    :language: xml
    :lines: 66

Location points are defined in the ``regions`` element in anticipation of 
hydraulic head and pressure observations. 
To avoid ambiguity and the need for interpolation, 
the observation points are defined so that each
location coincides exactly with a grid node (computational cell center). 
The grid spacing in the :math:`x`-direction is 1 meter, so cell centers are positioned at
:math:`x_i = 0.5, 1.5, 2.5, ..., 99.5`. Similarly, the :math:`y`-coordinates of grid 
nodes are :math:`y_i = 5, 10, 15, 20,` and :math:`25`, and finally 
:math:`z_i = 0.5, 1.5, 2.5, ..., 49.5`. Point ``Well 1`` is located just inside the 
inlet boundary near the center of the face. Point ``Well 3`` is similarly located 
one-half cell inside the downstream boundary. Point ``Well 2`` is located near the 
center of the domain (within one-half cell). Points ``Well 2t`` and ``Well 2b`` are
positioned above and below ``Well 2`` at the top and bottom of the domain. 
Although the simulation problem is steady-state, an initial condition is required
**?why?**. The *Amanzi* ``pressure`` variable is absolute pressure. 
Here atmospheric pressure is specified (``101325.0`` Pa). The boundary condition
at :math:`x=0` is defined in the problem statement as a volumetric flux 
[m\ :sup:`3`\ /m\ :sup:`2`\ d = m/d]. *Amanzi* currently requires boundary conditions
of this type to be specified as a mass flux, 
:math:`\rho U` [kg\ :sup:`3`\ /m\ :sup:`2`\ s]. Performing this calculation for 
:math:`\rho = 998.2` kg/m\ :sup:`3`\  and 
:math:`U = 1.688` m/d = :math:`1.954 \times 10^{-5}` m/s
yields ``1.95e-2`` kg\ :sup:`3`\ /m\ :sup:`2`\ s for ``inward_mass_flux``. 
For this steady-state saturated flow-only simulation porosity does not affect
the solution and its value is arbitrarily set to 0.25. 
In the ``definitions`` element, the ``time_macro`` element defines the time at which
model observations will be written to output according to the 
``output/observations`` element that follows. Because the simulation is steady-state,
the ``time`` is set to zero. Observations are requested for hydraulic head and 
(absolute) pressure at the five well locations, and to be written to 
``observations.out``.

*Amanzi* execution
~~~~~~~~~~~~~~~~~~

*Amanzi* is executed from the command-line using this or analogous command for
the user's specific installation:

``mpirun -n 4 /ascem/amanzi/install/current/bin/amanzi --xml_schema=/ascem/amanzi/install/current/bin/amanzi.xsd --xml_file=amanzi_steady_linear-isv2.xml``

Here execution using 4 processor cores is specified by ``-n 4``. Successful completion
is marked by ``SIMULATION_SUCCESSFUL`` followed by a timing summary ::

	Amanzi::SIMULATION_SUCCESSFUL

	**********************************************************
	***                   Timing Summary                   ***
	**********************************************************



*Amanzi* results
~~~~~~~~~~~~~~~~

For the requested observation points, the following results are expected in terms
of hydraulic head (Equation :eq:`observationHead`) and pressure 
(:math:`P = (h - z)\rho g + P_{atm}`):

+---------+-------+----------+-----------+
|  ID     | x (m) |   h (m)  |   P (Pa)  |
+=========+=======+==========+===========+
| Well 1  |  0.5  |  319.00  |  2974498  |
+---------+-------+----------+-----------+
| Well 2  | 50.5  |  221.00  |  1995564  |
+---------+-------+----------+-----------+
| Well 2t | 50.5  |  221.00  |  1760619  |
+---------+-------+----------+-----------+
| Well 2b | 50.5  |  221.00  |  2240297  |
+---------+-------+----------+-----------+
| Well 3  | 99.5  |  121.00  |  1036208  |
+---------+-------+----------+-----------+

In output file ``observations.out``, *Amanzi* produces the following results:

+---------+-------+----------+-----------+
|  ID     | x (m) |   h (m)  |   P (Pa)  |
+=========+=======+==========+===========+
| Well 1  |  0.5  |  318.96  |  2974022  |
+---------+-------+----------+-----------+
| Well 2  | 50.5  |  218.98  |  1995310  |
+---------+-------+----------+-----------+
| Well 2t | 50.5  |  218.98  |  1760374  |
+---------+-------+----------+-----------+
| Well 2b | 50.5  |  218.98  |  2240035  |
+---------+-------+----------+-----------+
| Well 3  | 99.5  |  121.00  |  1036172  |
+---------+-------+----------+-----------+

The output from *Amanzi* is observed to agree closely with the expected hydraulic
head and pressure results from the analytic solution.

*Amanzi* XML input file
~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: amanzi_steady_linear-isv2.xml
    :language: xml
