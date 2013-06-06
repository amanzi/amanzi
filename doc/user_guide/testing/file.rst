Theis Analysis
~~~~~~~~~~~~~~~

Theis (1935) developed an analytical solution for transient (non steady state) drawdown for a fully penetrating well by imposing the boundary conditions:
:math:`h = h_0` for 
:math:`t = 0` and 
:math:`h \Rightarrow h_0` as 
:math:`r \Rightarrow \infty`.  The equation assumes an infinite and uniform unconfined aquifer.  When a well is pumped the water table declines toward the well and flow is induced toward the well from all directions. Theoretically, this flow can be idealized by purely radial symmetric flow and can be decribed by the equation below (this is analogous to heat flow by conduction developed by Fourier).

:math:`\frac{\partial^2 h}{\partial r^2} + \frac{1}{r} \frac{\partial h}{\partial r} = \frac{S}{T} \frac{\partial h}{\partial t}`

The analytical solution of drawdown as a function of time and distance is found to be:

:math:`s = h_o - h(r,t) = \frac{Q W(u)}{4 \pi T} = \int_u^\infty \frac{exp[-u]}{u} du = \frac{Q W(u)}{4\pi T}`

where, 
:math:`u(r,t) = \frac{r^2 S}{4 T t}`

Note,
:math:`W(u)` can easily be determined from existing tables once 
:math:`u(r, t)` is found which is a measure of aquifer response time. For a given value of *t* one can construct a draw down curve with respect to the distance from the pumping well, *r*.  

* *Q* constant pump rate
* :math:`h_o` initial water table table height
* *T* transmissivity (product of conductivity and saturated thickness)
* *W(u)* well function
* *r* radial distnace measured outward from well
* *S* storage coefficient 
* *t* duration of pumping time
  
