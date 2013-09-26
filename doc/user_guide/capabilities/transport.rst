Transport Models
----------------

We use "Transport" to refer to the set of physical processes 
that lead to movement of dissolved and solid contaminants in the subsurface, 
treating the chemical and microbiological reactions that can affect the transport 
rate through a retardation effect as a separate set of processes.  
The principal transport processes to be considered are advection, mechanical dispersion,
and molecular diffusion.  

Recall the equation for mass conservation can be written as

.. math::
  \frac{\partial (\phi \sum_\alpha [s_{\alpha} C_\alpha])}{\partial t} 
  +
  \boldsymbol{\nabla} \cdot \boldsymbol{J}_{\text{adv}} 
  = Q
  - \boldsymbol{\nabla} \cdot \boldsymbol{J}_{\text{disp}} 
  - \boldsymbol{\nabla} \cdot \boldsymbol{J}_{\text{diff}}

where :math:`\boldsymbol{J_{\text{adv}}}` refers to the advective flux, 
:math:`\boldsymbol{J}_{\text{disp}}` is the dispersive flux, 
:math:`\boldsymbol{J}_{\text{diff}}` is the diffusive flux (often grouped 
with the dispersive flux), and  :math:`Q` is the 
summation of the various source terms (which may include reactions).


Advective Flux
~~~~~~~~~~~~~~

Advection involves the translation in space of dissolved or suspended material at the 
rate of movement of the bulk fluid phase.  
No modification of the shape of a front and no dilution occurs when transport 
is purely via advection---a sharp front remains so when undergoing purely advective transport.  
The advective flux, :math:`\boldsymbol{J}_{adv}`, of a dissolved species in porous 
media can be described mathematically as

.. math::
  \boldsymbol{J}_{adv} =\phi s_{\alpha} \boldsymbol{v}_{\alpha} C_{i},  

where :math:`\phi` is the porosity, :math:`s_{\alpha}` is the saturation of 
phase :math:`\alpha`, and  :math:`\boldsymbol{v}_{\alpha}` is the average linear 
velocity of the phase, and :math:`C_i` is the concentration of the :math:`i`-th 
species.


Dispersive Flux
~~~~~~~~~~~~~~~



Diffusive Flux
~~~~~~~~~~~~~~



Boundary Conditions
~~~~~~~~~~~~~~~~~~~



Source Terms
~~~~~~~~~~~~
