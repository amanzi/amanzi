Model
======

Setup
-------

Confined Aquifer
~~~~~~~~~~~~~~~~

For Theis analysis to applicable three parameters must be known:

* *T*, transmissivity
* *S*, storativity
* *Q*, constant pumping rate

Transimissivity is defined as 

.. math:: T = Kb

where *K* is the hydraulic conductivity of the aquifer and *b* is the saturated thickness.  Transmissivity values greater than 0.015
:math:`\frac{m^2}{s}` represent aquifers capable of well exploitation.  
Storatvity is a dimensionless parameter that describes the amount of water released by the aquifer per unit volume of the aquifer.  Storativity can be calculated using 

.. math:: S = S_s b

where
:math:`S_s` is the *specific storage* unique to each aquifer.  Again, *b* is the saturated thickness of the aquifer.  Specific storage represents the volume of water released per unit volume of the aquifer per unit decline in hydraulic head.  

Lastly, the constant pumping rate, *Q*, is the volume of water discharged from the well per unit time.  

Unconfined Aquifer
~~~~~~~~~~~~~~~~~~

Note, there are multiple models to estimate the drawdown in unconfined aquifers with varying degrees of accuracy.  An important difference between a confined and an unconfined aquifer is the hydraulic gradients that are created.  This desrcibes how the water is flowing in the aquifer to the well. 

For example, the flow in a confined aquifer is purely horizontal and no dewatering of the geologic system occurs.  The mechanisms for water production in these wells are 1) expansion of the water and 2) compaction of the aquifer.  However, in an unconfined aquifer there is both an induced horizontal and vertical flow.  The water produced in these wells come from confined delivery and the dewatering of the unconfined aquifer.  Hence, there is a new coefficent of storage called *specific yield*.  This paramter is always far greater than the storativities of confined aquifers.    

There are three approaches to estimate the cone of depression for an unconfined aquifer.  The three approaches are: 

1. saturated-unsaturated flow system in which water table drawdowns are accompanied by a change in the unsaturated moisture content above the water table
2. Delayed water table response
3. Confined aquifer equation (simplest)

It is important to note that when using the third option the specific yield term is used rather than the storativity coefficient, *S*.  Also the saturated thickness, *b* is considered to be the initial height of the water table.  




Schematic
~~~~~~~~~

