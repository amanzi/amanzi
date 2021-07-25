Units in Amanzi
---------------

Internally, *Amanzi* uses SI units for length, time, mass, and
concentration as indicated in the ``model_description/units``
element. The only exception is the mol/L unit for the concentration.
Although, it will eventually support the automatic
translation of popular units in an input file to SI units internally,
at present, this translation is only supported by the unstructured component of *Amanzi*. 
For example, the units for density, viscosity, pressure, hydraulic
conductivity, and permeability are kg/m\ :sup:`3`\ , Pa :math:`\cdot`
s, Pa, m/s, and m\ :sup:`2`\ , respectively.  Simulation outputs also
use these units.

Beyond dimensional units, the physical entities that are tracked in
internal *Amanzi* calculations, and are often preferred for *Amanzi*
input and output, may be less familiar to those accustomed to
hydrologic applications with constant material properties.
Specifically, to address the full range of multiphase-multicomponent
reactive transport simulation capabilities specified in the Amanzi
process models design document, *Amanzi*'s numerical calculations
use the following physical parameters:

      *	Absolute pressure :math:`P` (Pa) rather than
	total/hydraulic head :math:`h = (P-P_{atm})/\rho g + z` (m) or 
	pressure head :math:`(P-P_{atm})/\rho g` (m)

      *	Intrinsic permeability :math:`k` (m\ :sup:`2`\ ) instead of
	hydraulic conductivity :math:`K = \rho g k/\mu` (m/s)

      *	Mass flux :math:`\dot{m} = \rho \dot{q}` (kg/m\ :sup:`2`\ s) instead of 
	volumetric flux (Darcy velocity) :math:`\dot{q}` 
	(m\ :sup:`3`\ /m\ :sup:`2`\ s = m/s)

Pressure is preferred over head for variable fluid density because
:math:`\rho` is no longer a constant in the head expression
:math:`(P-P_{atm})/\rho g`.  Hydraulic conductivity :math:`K = \rho g
k/\mu` (m/s) is a property of both the porous matrix material and the
fluid. Permeability, which is strictly a material property, is a more
fundamental entity. When the fluid properties are constant, hydraulic
conductivity is proportional to intrinsic permeability and a
convenient surrogate for the latter.  However, permeability is
preferred when fluid properties are variable to separate the effects
of variations in the porous medium versus the fluid filling its pore
space. Similarly, when fluid properties are constant, volumetric flux
is proportional to mass flux and a convenient surrogate for the
latter, but direct specification of mass flux is more convenient for
varying fluid properties.
