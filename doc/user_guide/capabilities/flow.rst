Flow Models
-----------



Single-Phase Richards
~~~~~~~~~~~~~~~~~~~~~

The Richards equation provides the necessary physics to represent the
flow of a liquid phase (typically water) under partially saturated
conditions, with the assumption made that the gas phase is
infinitely mobile and has no effect on the flow of the aqueous phase.  
The Richards equation is derived from the conservation of
liquid mass continuity equation:

.. math::
  \frac{\partial (\phi\, s_l \rho)}{\partial t} 
  =
  \boldsymbol{\nabla} \cdot (\rho \boldsymbol{q}) + Q,

where :math:`s_l` is the liquid saturation, 
:math:`\phi` is the porosity,
:math:`\rho` is the liquid density, 
:math:`Q` is a source or sink,
and :math:`\boldsymbol{q}` is the liquid velocity given by Darcy's Law:

.. math::
  \boldsymbol{q} 
  = -\frac{\boldsymbol{K} k_r}{\mu} 
  (\boldsymbol{\nabla} p - \rho \boldsymbol{g}),
where :math:`\boldsymbol{K}` is the tensor of absolute perbeability of rock,
:math:`k_r(s_l)` is the relative permeability,
:math:`\mu` is the liquid viscosity,
:math:`p` is the liquid pressure,
and :math:`\boldsymbol{g}` is the gravity vector.

Typical models of relative permeability 
are the van Genuchten-Mualem relations and the Brooks-Corey-Burdine relations.
The primary variable is liquid pressure :math:`p`, so that
we write these models as functions of capillary pressure :math:`p_c`.


Van Genuchten Relative Permeability
...................................

The Mualem relative permeability function for the liquid phase derived
from the van Genuchten saturation function is given by

.. math::
  k_r = s_e^{\ell} \left\{1 - \left[1 - s_e^{1/m} \right]^m \right\}^2.

The Burdine relative permeability function has the form

.. math::
  k_r = s_e^{\ell} \left\{1 - \left[1 - s_e^{1/m} \right]^m \right\}.


Brooks-Corey Relative Permeability
..................................

Combined with the Brooks-Corey saturation function, the Mualem
relative permeability function is given by

.. math::
  k_r = \big(s_e\big)^{\ell+2+2/\lambda}
      = \left(\alpha | p_c | \right)^{-((\ell+2)\lambda+2)}.

The Burdine form originally considered by (Brooks, 1964) is given by

.. math::
  k_r = \left(s_e \right)^{ \ell+1+2/\lambda}
      = \left(\alpha |p_c| \right)^{-((\ell+1)\lambda+2)}.


Single-Phase Saturated
~~~~~~~~~~~~~~~~~~~~~~

The most basic case of flow is that of a single phase in a porous medium.  
Notwithstanding its simplicity, it has a wide application to describing 
subsurface processes.
The governing equations are the mass balance and Darcy's Law:

.. math::
  \frac{\partial (\phi \rho)}{\partial t} 
  + \boldsymbol{\nabla}\cdot(\rho \boldsymbol{q}) = Q,

where :math:`\rho` denotes the fluid density (sometimes considered to be a
function of :math:`p`), :math:`Q` represents a source/sink term, 
and :math:`\boldsymbol{q}` is the Darcy velocity:

.. math::
  \boldsymbol{q} 
  = -\frac{\boldsymbol{K}}{\mu} (\boldsymbol{\nabla} p - \rho \boldsymbol{g}).


Boundary Conditions
~~~~~~~~~~~~~~~~~~~









