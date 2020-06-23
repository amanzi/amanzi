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
		      
The *Amanzi* input specification document (**insert link**) defines 
the content and format of an XML input file. A complete *Amanzi* input file
for this tutorial is listed in a separate section below. 
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

.. amanzi_xml_include:: amanzi_steady_linear-isv2.xml
    :string: <liquid_phase name="water">
    :indent: 2

and element *content*, e.g. ``998.2`` kg/m\ :sup:`3`\  in the ``density`` element,

.. amanzi_xml_include:: amanzi_steady_linear-isv2.xml
    :string: <density>
    :indent: 2

as defined by the *Amanzi* input specification. 
Although not required by the XML standard or *Amanzi*, 
indentation using tabs and/or spaces is commonly used to indicate
the hierarchy of elements and improve readibility. Tabs are used in this tutorial 
example to show the file hierarchy. Although *Amanzi* XML element names and the 
file structure are intended to be reasonably self-explanatory, selective commentary 
is provided here for further guidance. 

*Amanzi* offers the user two numerical gridding approaches. 
The ``structured`` grid option refers to an orthogonal grid of rectangular (2D) or 
brick (3D) elements that may be further subdivided through 
Adaptive Mesh Refinement (AMR). The ``unstructured`` grid option can accommodate a
network of non-orthogonal cells that may not be connected on a regular
pattern defined by :math:`(i,j,k)` coordinate indices. The ``unstructured``
option also accommodates structured grids, but not multi-level AMR as with the
``structured`` grid option. A simple orthogonal structured grid with no refinement
will be used for this tutorial, and either gridding option would suffice. 
However, the ``unstructured`` option is selected in this example through the 
``type`` attribute in the ``amanzi_input`` XML element:

.. amanzi_xml_include:: amanzi_steady_linear-isv2.xml
    :string: <amanzi_input

Nested within the ``amanzi_input`` element are additional XML elements defining
various attributes of the numerical simulation to be performed. Although the order
of these elements is arbitrary, users will typically want to sequence elements to 
mirror their conceptualization of the physical problem. For this tutorial we begin
with the ``model_description`` element:

.. amanzi_xml_include:: amanzi_steady_linear-isv2.xml
    :start: <model_description
    :end: </model_description>

Besides model identification, this block defines the dimensional units to be used 
for the numerical simulation. Currently *Amanzi* requires the SI units indicated 
in the ``model_description/units`` element.
The units of density, viscosity, pressure, hydraulic conductivity, and 
permeability  are thus kg/m\ :sup:`3`\ , Pa :math:`\cdot` s, Pa, m/s, 
and m\ :sup:`2`\ , respectively. 
Simulation outputs also use these dimensional units. 

The next input block is the ``process_kernels`` element:

.. amanzi_xml_include:: amanzi_steady_linear-isv2.xml
    :start: <process_kernels>
    :end: </process_kernels>

The tutorial problem involves only fully-saturated flow, so the ``saturated`` model
is selected here, and transport and chemistry are explicitly turned ``off``. 
The latter may be implicitly turned off by omitting those two entries.

In the ``phases`` element

.. amanzi_xml_include:: amanzi_steady_linear-isv2.xml
    :start: <phases>
    :end: </phases>

``water`` is the only phase present (as a liquid), and 
equation-of-state (eos) computations are turned off in favor of direct 
specification of fluid density and viscosity as 998.2 kg/m\ :sup:`3`\ 
and 1.002e-03 Pa :math:`\cdot` s, respectively.

Under ``execution_controls`` the ``steady`` mode and first-order backward-difference
(backward Euler) differencing scheme (``bdf1``) are chosen:

.. amanzi_xml_include:: amanzi_steady_linear-isv2.xml
    :start: <execution_controls>
    :end: </execution_controls>

The ``start`` and ``end`` times may both be zero for a steady-state calculation. 
Controls for the selected numerical scheme are defined in the block

.. amanzi_xml_include:: amanzi_steady_linear-isv2.xml
    :start: <numerical_controls>
    :end: </numerical_controls>

Here recommended default values have been selected.

The computational mesh may be imported (read from file) or, for simple geometries, 
generated by *Amanzi*. A simple mesh is adequate for this simulation so the
latter option is selected for this tutorial example in the ``mesh`` element:

.. amanzi_xml_include:: amanzi_steady_linear-isv2.xml
    :start: <mesh framework="mstk">
    :end: </mesh>

Note that a 3D domain is specified through the ``dimension`` element. The
coordinates under ``mesh/generate/box`` refer to the 
:math:`(x_{min}, y_{min}, z_{min})` and :math:`(x_{max}, y_{max}, z_{max})` 
corners of the computational domain, and must be 
comma-separated. Whitespace is allowed and parentheses may be used to bracket 
coordinates in the form ``(x,y,z)`` if desired. The ``nx``, ``ny``, and ``nz``
values define the number of computational cells in each coordinate direction. 
Discretization is specified for the :math:`y` and :math:`z` coordinates,
although no variability is expected in these directions.
The user could effectively create a one-dimensional simulation by setting
``ny = 1`` and ``nz = 1``.

The ``regions`` element is used to define volumes for material property assignments

.. amanzi_xml_include:: amanzi_steady_linear-isv2.xml
    :start: <regions>
    :end: </region>

surfaces for boundary condition assignments

.. amanzi_xml_include:: amanzi_steady_linear-isv2.xml
    :start: <region name ="Upstream">
    :end: </region>

and discrete points

.. amanzi_xml_include:: amanzi_steady_linear-isv2.xml
    :lines: 66-71

The location points are defined in anticipation of 
hydraulic head and pressure observations to be specified below. 
To avoid ambiguity and the need for interpolation, 
the observation points are defined so that each
location coincides exactly with a grid node (computational cell center). 
The grid spacing in the :math:`x`-direction is :math:`1.0` meter, 
so cell centers are positioned at :math:`x_i = 0.5, 1.5, 2.5, ..., 99.5`. 
Similarly, the :math:`y`-coordinates of grid 
nodes are :math:`y_i = 5, 10, 15, 20,` and :math:`25`, and finally 
:math:`z_i = 0.5, 1.5, 2.5, ..., 49.5`. Point ``Well 1`` is located just inside the 
inlet boundary near the center of the face. Point ``Well 3`` is similarly located 
one-half cell inside the downstream boundary. Point ``Well 2`` is located near the 
center of the domain (within one-half cell). Points ``Well 2t`` and ``Well 2b`` are
positioned above and below ``Well 2`` at the top and bottom of the domain. 

Although the simulation problem is steady-state, an initial condition is required
**? why ?**:

.. amanzi_xml_include:: amanzi_steady_linear-isv2.xml
    :lines: 73-83

The *Amanzi* ``pressure`` variable is absolute pressure. 
Here atmospheric pressure(``101325.0`` Pa) is specified  over the modeling domain 
using the ``Entire Domain`` region specified above.

Next boundary conditions are specified as

.. amanzi_xml_include:: amanzi_steady_linear-isv2.xml
    :lines: 85-105

The boundary condition at :math:`x=0` is defined in the problem statement as a 
volumetric flux :math:`U` [m\ :sup:`3`\ /m\ :sup:`2`\ d = m/d]. 
*Amanzi* currently requires boundary conditions of this type to be specified as 
a mass flux, :math:`\rho U` [kg\ :sup:`3`\ /m\ :sup:`2`\ s]. 
Performing this calculation for :math:`\rho = 998.2` kg/m\ :sup:`3`\  and 
:math:`U = 1.688` m/d = :math:`1.954 \times 10^{-5}` m/s
yields ``1.95e-2`` kg\ :sup:`3`\ /m\ :sup:`2`\ s for ``inward_mass_flux``. 

Material properties are defined in the ``materials`` element

.. amanzi_xml_include:: amanzi_steady_linear-isv2.xml
    :lines: 107-116

For this steady-state saturated flow-only simulation porosity does not affect
the solution and its value is arbitrarily set to 0.25. The ``permeability`` values
are intrinsic permeability :math:`k` [m\ :sup:`2`\] rather than hydraulic 
conductivity :math:`K = \rho g k / \mu` [m/s], which is sometimes loosely referred
to as "permeability".
 
In the ``definitions`` element

.. amanzi_xml_include:: amanzi_steady_linear-isv2.xml
    :lines: 118-124

the ``time_macro`` element defines the time at which
model observations will be written to output according to the 
``output/observations`` element that follows. Because the simulation is steady-state,
the ``time`` is set to zero. 

The final input is the ``output`` block

.. amanzi_xml_include:: amanzi_steady_linear-isv2.xml
    :lines: 126-139

:math:`\vdots`

.. amanzi_xml_include:: amanzi_steady_linear-isv2.xml
    :lines: 156-165

:math:`\vdots`

.. amanzi_xml_include:: amanzi_steady_linear-isv2.xml
    :lines: 183-188

Observations, to be written to ``observations.out``, are requested for hydraulic head 
and (absolute) pressure at the five well locations previously defined under 
``regions/point``. The ``vis`` element specifies 5 significant figures for the
visualization file with prefix ``steady-flow``.

*Amanzi* execution
~~~~~~~~~~~~~~~~~~

*Amanzi* is executed from the command-line using this or analogous command for
the user's specific installation:

``mpirun -n 4 install/current/bin/amanzi --xml_schema=install/current/bin/amanzi.xsd --xml_file=amanzi_steady_linear-isv2.xml``

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
| Well 1  |  0.5  |  319.01  |  2974465  |
+---------+-------+----------+-----------+
| Well 2  | 50.5  |  219.00  |  1995306  |
+---------+-------+----------+-----------+
| Well 2t | 50.5  |  219.00  |  1760595  |
+---------+-------+----------+-----------+
| Well 2b | 50.5  |  219.00  |  2240256  |
+---------+-------+----------+-----------+
| Well 3  | 99.5  |  121.00  |  1036175  |
+---------+-------+----------+-----------+

The output from *Amanzi* is observed to agree closely with the expected hydraulic
head and pressure results from the analytic solution.
Minor mismatch is due to rounding of input data.

*Amanzi* XML input file
~~~~~~~~~~~~~~~~~~~~~~~

.. amanzi_xml_include:: amanzi_steady_linear-isv2.xml


