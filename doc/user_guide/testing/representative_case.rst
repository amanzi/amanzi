Representative Case
===================

Example
--------

Saturated 3-D Steady-State Test Flow Problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Problem Definition
~~~~~~~~~~~~~~~~~~~

* Horizontal flow in the **x-direction**

* Model domain:

	* :math:`X_{min} = Y_{min} = Z_{min} = 0`
	* :math:`Y_{max} = Z_{max} = 50` *meters*

* Model domain: 10 meters X 5 meters X 5 meters

* Total flow applied (flux) at **x = 0 plane** is *1.95E-2* 
  :math:`kg/s/m^2`

* Constant head specified at **x = 100 meters plane**: 120 meters 

* Uniform Intrinsic Permeability, 
  :math:`\kappa =` *1.0E-12*
  :math:`m^2`

* Water Viscosity,
  :math:`\mu =` *1.002E-3 Pa - s*

* Water Density,
  :math:`\rho =` *998.2*
  :math:`kg/m^3`

Expected result, according to Darcy's Law
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

Darcy's Law is defined as: 

.. math:: Q/A = -K \frac{dh}{dx}

which, in this case, is based on isotropic soil conductivity, *K*, and the hydraulic head gradient in the x-direction.  

Calculated head gradient (-): ?

Calculated hydraulic conductivity:

.. math:: K = \kappa \rho g / \mu = 9.67*10^{-6} \frac{m}{s}

         \frac{Q}{A} = \frac{1.95*10^{-2}}{998.2} \frac{kg/s/m^2}{kg/m^3} = 9.76*10^{-6} \frac{m}{s}

**Pressure p(x,z)** = 
:math:`[ h(x) - z ]\rho g`

+------+------+------+-----------------------+
| x [m]| z [m]| h [m]|       P [MPa]         |
+======+======+======+========+==============+
|      |      |      |Expected|Amanzi Output |                        
+------+------+------+--------+--------------+
|0.5   |49.5  |319   |2.63    |2.74          |
+------+------+------+--------+--------------+
|50.5  |25.5  |219   |1.90    |2.0           |
+------+------+------+--------+--------------+
|50.5  |49.5  |219   |1.65    |1.76          |
+------+------+------+--------+--------------+
|50.5  |0.5   |219   |2.14    |2.24          |
+------+------+------+--------+--------------+
|99.5  |49.5  |121   |0.70    |0.80          |
+------+------+------+--------+--------------+

Code
-----

The code below is an example from the same values from above.  It is a
simple case with 1D flow.  The code presented will illustrate how
users can set up their input file. 

Structure
~~~~~~~~~

The strucutre of the Amanzi code is presented below.  The code
begins with Execution Control by setting up a mesh.  Here the user will
identify three types of models: flow, transport, and chemistry.  Next,
the user must chose between transient or steady-state time integration
mode.  
  
.. literalinclude:: f.xml
    :language: xml
    :lines: 14-44
   

Then the mesh is generated from the user given size.  The number of cells
per plane is specified to utilize the finite difference method to
solve for unkowns.  The domain is further specified by defining the
lower corner and the upper corner.  

.. literalinclude:: f.xml
     :language: xml
     :lines: 46-64

The next major portion is defining the Region which consisting of
domains, boundaries, and well locations.  [talk about domain again,
really what is this defining].

.. literalinclude:: f.xml
     :language: xml
     :lines: 66-79

The boundaries are specified to define the 'upstream' domain and the
'downstream' domain.  

.. literalinclude:: f.xml
     :language: xml
     :lines: 80-94

Then the coordinates of each well are specified
to determine the pressure at that location within the well. Note, finite
difference method solves linear differential equations exactly at the
center of each cell and approximates nonlinear equation.  If a well
location is not in center of a cell increasing the amount of cells
will yield a closer approximation effectively reducing the distance of
the where the solution is derived and the well location. 

.. literalinclude:: f.xml
     :language: xml 
     :lines: 95-126

Lastly, the properties of the subsurface and the fluid are defined.
The subsurface properties are permeability and characterization of the
porosity i.e. uniform or nonuniform.  The fluid, typically water, is
defined by its density and viscosity. 

.. literalinclude:: f.xml
     :language: xml
     :lines: 128-163

Amanzi can then run once initial coniditions and boundary conditions
are applied.     

 
