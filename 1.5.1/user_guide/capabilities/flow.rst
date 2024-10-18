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
  -\boldsymbol{\nabla} \cdot (\rho \boldsymbol{q}) + Q,

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
 
When the Richards process kernel is a part of a complex multi-process non-isothermal simulation,
fluid density and viscosity are functions of temperature and pressure.



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
  \left(\frac{S_s}{g} + \frac{S_y}{Lg}\right)
    \frac{\partial p}{\partial t} 
  + \boldsymbol{\nabla}\cdot(\rho \boldsymbol{q}) = Q,

where :math:`S_s` is the specific storage, :math:`S_y` is the specific yield,
:math:`g = \| \boldsymbol{g}\|` is the gravity value, 
:math:`p` is pressure, 
:math:`Q` is a source/sink term, 
and :math:`\boldsymbol{q}` is the Darcy velocity:

.. math::
  \boldsymbol{q} 
  = -\frac{\boldsymbol{K}}{\mu} (\boldsymbol{\nabla} p - \rho \boldsymbol{g}).

The specific storage controls the amount of water that the aquifer releases from
storage while remaining fully saturated.
For cells at the water table where confining conditions do not exist, specific 
yield is used instead of specific storage.
Specific yield is the drainable porosity when the water table moves (e.g.
during pumping), indicating the fraction of the aquifer volume that will be drained 
under the force of gravity.
It is implemented in *Amanzi* using the characteristic vertical size 
:math:`L` of a yield cell. 
In the case of a structured mesh, :math:`L` is the cell height. 

The specific yield is an approximation of a more complex process of movement
of the water table, e.g. using the Richards equation.
Without :math:`S_y`, the linear model overpredicts the pressure response at 
observation wells because the water released from drainage is orders of magnitude
greater than the amount of water released from elastic deformation, as 
simulated with :math:`S_s`.


Boundary Conditions
~~~~~~~~~~~~~~~~~~~
For a single phase model, three types of boundary conditions are supported:
(a) prescribed pressure or head; 
(b) prescribed flux; 
(c) seepage face. 

The prescribed pressure or head (Dirichlet) boundary condition
occurs whenever the flow domain is adjacent to a body of open water.
Mathematically, it is given by 

.. math::
  p(\boldsymbol{x},t) = p_0(\boldsymbol{x},t),
  \qquad \boldsymbol{x} \in \Gamma_D, \  t > 0,

and :math:`p_0` is the prescribed pressure function.

The flux (Neumann) boundary condition specifies 
the flux normal to the boundary surface or, equivalently, the Darcy velocity 
normal to the boundary:

.. math::
  \boldsymbol{q}(\boldsymbol{x},t) \cdot \boldsymbol{n} = q_0(\boldsymbol{x},t),
  \qquad \boldsymbol{x} \in \Gamma_N, \  t > 0,

where :math:`q_0` is the given function. 
This condition is often used for surface infiltration and 
and :math:`q_0` is refered to as the infiltration velocity.

The seepage boundary condition is defined as the boundary where water 
leaves the ground surface and continues to flow in a thin film along its surface.
The seepage face boundary condition becomes the atmospheric pressure condition there.
This is a nonlinear condition since 
the geometry of the seepage face is defined by the solution and may vary with time. 
In *Amanzi*, this boundary condition is simulated using a dynamic switch from the
prescribed pressure boundary condition to a prescribed flux boundary condition 
representing the recharge.


Initial Conditions
~~~~~~~~~~~~~~~~~~

An initial conditon specifies pressure at time :math:`T=0` inside the
computational domain:

.. math::
     p(\boldsymbol{x}, 0) = p_{0}(\boldsymbol{x})
  \qquad \boldsymbol{x} \in \Omega.


Source Terms
~~~~~~~~~~~~

The source term :math:`Q` is a given function specifying usually location of
wells inside the computational domain:

.. math::
     Q(\boldsymbol{x}, t) = Q_{0}(\boldsymbol{x},t),
  \qquad \boldsymbol{x} \in \Omega,\quad t > 0.











