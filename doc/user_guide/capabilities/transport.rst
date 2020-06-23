Transport Models
----------------

The transport component of *Amanzi* models a set of physical processes 
that lead to movement of dissolved and solid contaminants in the subsurface, 
treating the chemical and microbiological reactions that can affect the transport 
rate through a retardation effect as a separate set of processes.  
The governing PDE for the component :math:`C_i` reads

.. math::
  \frac{\partial (\phi s_l C_i)}{\partial t} 
  + \nabla \cdot \boldsymbol{J}_i^{\text{adv}} 
  = Q_i 
  - \nabla \cdot \boldsymbol{J}_i^{\text{disp}} 
  - \nabla \cdot \boldsymbol{J}_i^{\text{diff}}

where :math:`\phi` is the porosity, :math:`s_l` is the liquid saturation, 
:math:`\boldsymbol{J}_i^\text{adv}` is the advective flux, 
:math:`\boldsymbol{J}_i^\text{disp}` is the dispersive flux, 
:math:`\boldsymbol{J}_i^\text{diff}` is the diffusive flux (often grouped 
with the dispersive flux), and  :math:`Q_i` is the 
summation of the various source terms.


Advective Flux
~~~~~~~~~~~~~~

Advection involves the translation in space of dissolved or suspended material at the 
rate of movement of the bulk fluid phase.  
The advective flux of a dissolved species :math:`C_i` in porous 
media can be described mathematically as

.. math::
  \boldsymbol{J}_i^\text{adv} = \boldsymbol{q}_l C_{i},  

where :math:`\boldsymbol{q}_l` is the Darcy velocity.


Dispersive Flux
~~~~~~~~~~~~~~~

Dispersion of a dissolved species involves its spreading along tortuous pathways 
in a porous medium caused by mixing effects.  
Dispersion takes place in the direction of the flow (longitudinal) and normal to the 
flow (transverse).  
A conventional Eulerian Fickian representation of dispersion is assumed, which may be 
taken as the asymptotic limiting form of the dispersion tensor. 
The dispersive flux has the form

.. math::
  \boldsymbol{J}_i^\text{disp} = - \phi s_l \boldsymbol{D} \nabla C_i,

where :math:`\boldsymbol{D}` denotes the dispersion tensor.
The dispersion tensor takes different forms depending on whether the media 
is isotropic or anisotropic. 
For an isotropic medium with no preferred axis of symmetry the dispersion 
tensor has the well-known form:

.. math::
  \boldsymbol{D} 
  = \alpha_T \|\boldsymbol{v}\| \boldsymbol{I} 
  + \left(\alpha_L-\alpha_T \right) 
    \frac{\boldsymbol{v} \boldsymbol{v}}{\|\boldsymbol{v}\|},

characterized by the two parameters :math:`\alpha_L` [m] and :math:`\alpha_T` [m] 
referred to as the longitudinal and transverse dispersivity, respectively. 
The vector :math:`\boldsymbol{v}` [m/s] denotes the average pore velocity,
and :math:`\boldsymbol{I}` is the identity matrix.  


Diffusive Flux
~~~~~~~~~~~~~~

Molecular diffusion is often indistinguishable from mechanical dispersion 
as a process operating in porous media, and thus the two are often lumped 
together to form a hydrodynamic dispersion term.  
Molecular diffusion is an entropy-producing process in which the random motion 
of molecules causes spreading or homogenization of a concentration field.  
Atomistic representations of molecular diffusion capture this random motion, 
but continuum models of the kind considered here typically represent only 
the average behavior of the molecules.  

Molecular diffusion is usually described in terms of Fick's First Law, which 
states that the diffusive flux is proportional to the concentration gradient.
Since water-rock interaction commonly takes place in porous materials, it is 
important to account for the effect of tortuosity, which is defined as the ratio 
of the path length the solute would follow in water alone, :math:`L`, relative 
to the tortuous path length it would follow in porous media, :math:`L_e`:

.. math::
  \tau_{L} = (L/L_e)^2.

The diffusive flux, then, is given by

.. math::
  \boldsymbol{J}_{i}^\text{diff} 
  = - \phi s_l D_i \tau_{L} \nabla C_i.

where :math:`D_i` is referred to as the diffusion coefficient and 
is specific to the chemical component considered as indicated by the subscript :math:`i`. 


Boundary Conditions
~~~~~~~~~~~~~~~~~~~

A first-type or Dirichlet condition involves specification of a fixed value of the concentration, 
:math:`C_i` at the boundary location:

.. math::
  C_i(\boldsymbol{x}, t) = C_{i,0}(\boldsymbol{x}, t)
  \qquad \boldsymbol{x} \in \Gamma^{in},\quad t > 0,

where :math:`C_0` is a given function.
Pure advective transport requires to set up the Dirichlet boundary conditions only on
the inflow boundary of the computational domain. 

A second-type or Neumann boundary condition involves specification of the flux

.. math::
  J_i(\boldsymbol{x}, t) = J_{i,0}(\boldsymbol{x}, t),
  \qquad \boldsymbol{x} \in \partial \Omega,\quad t > 0,

where :math:`J_{i,0}` is a given flux function. 
The Dirichlet boundary condition on the outflow part of the computational domain may 
result in parabolic and/or exponential boundary layers; therefore it should be used
with caution.


Initial Conditions
~~~~~~~~~~~~~~~~~~

An initial condition specifies concentration at time :math:`T=0` inside the
computational domain:

.. math::
  C_i(\boldsymbol{x}, 0) = C_{i,0}(\boldsymbol{x})
  \qquad \boldsymbol{x} \in \Omega.


Source Terms
~~~~~~~~~~~~

The source term :math:`Q_i` is a given function specifying usually location of
wells inside the computational domain:

.. math::
  Q_i(\boldsymbol{x}, t) = Q_{i,0}(\boldsymbol{x},t),
  \qquad \boldsymbol{x} \in \Omega,\quad t > 0.


