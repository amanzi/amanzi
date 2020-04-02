Chemistry Models
----------------

The chemistry component of Amanzi models a set of geochemical processes that lead to the transformation of reactant species to product species. Amanzi can simulate a comprehensive set of reaction types.  These reactions can be classified into *homogeneous* if they occur within the aqueous phase or *heterogeneous* if, in addition to the aqueous phase, they involve the solid phase. These reactions can also be classified into equilibrium or kinetic depending on whether local equilibrium can be assumed. The following table summarizes the reactions available by type.

+--------------------+----------------------------+---------------------------------------+
| *Reaction Types*   | Homogeneous                | Heterogeneous                         |
+--------------------+----------------------------+---------------------------------------+
| Equilibrium        | **Aqueous Complexation**   | **Surface Complexation**,             |
|                    |                            | **Ion Exchange**                      |
+--------------------+----------------------------+---------------------------------------+
| Kinetic            | **Aqueous Kinetics**       | **Mineral Dissolution-Precipitation** |
+--------------------+----------------------------+---------------------------------------+

Stoichiometry
~~~~~~~~~~~~~

The transformation from reactant to product species occurs according to a stoichiometric relationship. The stoichiometry of a generic reaction involving :math:`N_c` reactants (:math:`A_j`) and a product (:math:`A_i`) can be written as

.. math::
  \sum_{j=1}^{N_c} \nu_{ij}~A_j = A_i
  :label: stoichiometry

where :math:`\nu_{ij}` is the stoichiometric coefficients of reactant *j* in reaction *i*. 

Geochemical Equilibrium
~~~~~~~~~~~~~~~~~~~~~~~

The local equilibrium assumption is applicable when reaction rates are much faster than transport rates, and therefore the relationship between the concentrations of products and reactants can be expressed with the law of mass action:

.. math::
  \log(K_i) = \sum_{j=1}^{N_c}{\nu_{ij}~log(\gamma_j \cdot C_j)} - log(\gamma_i \cdot C_i)
  :label: equilibrium

where :math:`K_i` is the equilibrium constant, :math:`C_i` and :math:`C_j` are the concentrations of reactant and product species, and :math:`\gamma_i` and :math:`\gamma_j` are the corresponding activity coefficients. 

Total Component Concentrations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Geochemical equilibrium between aqueous species makes it possible to calculate the concentration of a set of secondary species as a function of a set of primary species. Further, one can define a set of components, one for each primary species, the mass of which are not affected by geochemical aqueous equilibrium reactions. The total mass of a component  (:math:`U_j`) is defined as the sum of the mass of a primary species and the mass of a number of secondary species multiplied by the corresponding stoichiometric coefficient:

.. math::
  U_j = C_j + \sum_{i=1}^{N_r}{\nu_{ij}~C_i}
  :label: totalcomponent

where :math:`N_r` is the number of aqueous equilibrium reactions.

This approach is widely used in geochemical models (:cite:`chem-Lichtner_1985`, :cite:`chem-Yeh_rt-evaluation_1989`, :cite:`chem-Steefel_1994`, :cite:`chem-Steefel_1996`) to reduce the size of the system to solve, and more specifically the number of transport equations.

Transport and Chemistry Coupling: Operator Splitting Approach
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Amanzi uses the operator spliting aproach to couple transport and reaction processes :cite:`chem-Yeh_rt-evaluation_1989`. In the operator splitting approach, transport and chemistry are solved sequentially over the same time step, first transport and then chemistry. One advantage of this approach is that transport equations can be solved using linear solvers. Additionally, the non-linear chemical equations can be solved on a cell-by-cell basis in what is an embarassingly parallel problem. To access these advantages, the transport solver does not operate on the concentrations of each reactant (:math:`C_j`) and product species (:math:`C_i`). Rather, transport operates on the total component concentrations (:math:`U_j`) as defined in eq. :eq:`totalcomponent`. Please see :doc:`transport` for details of the transport capabilities of Amanzi. Following transport, the chemical processes affecting the mass balance of all components are solved over the same timestep:

.. math::
   \frac{\partial (\phi s_l U_j)}{\partial t} 
  = r_j
  :label: operatorsplitting

This set of equations is solved together with eqs. :eq:`totalcomponent` and :eq:`equilibrium` to know the concentration of the individual species, as well as with the equations that describe the kinetic reaction rates :math:`r_j` (see eq. :eq:`kinetics` below).

Geochemical Kinetics
~~~~~~~~~~~~~~~~~~~~

When local equilibrium cannot be assumed, reaction rates need to be calculated explicitly. The reaction rates are in general a non-linear function of concentrations of the geochemical species. Additionally, for heterogeneous reactions, the rates can be a function of material properties (symbolically: :math:`\Psi`), e.g. reactive surface area in the case of mineral dissolution-precipitation. The particular expression depends on the reaction type. For example, mathematical formulations include the transition state theory rate law for mineral dissolution-precipitation or a first order dependence on concentration for radioactive decay. For the sake of brevitiy, the rate expression for reaction *k* is provided here only as a generic function:

.. math::
  r_k = f(C_i,\Psi)
  :label: kinetics

The specific mathematical formulations employed for the rate expression depend on the geochemical engine. The choices for geochemical engines are discussed below.

Geochemical Engines
~~~~~~~~~~~~~~~~~~~

Taking advantage of the operator splitting approach, whereby transport and chemistry are solved sequentially, Amanzi offers a flexible approach to use geochemical engines. The objective of these geochemical engines is to solve the chemistry problem, i.e. equations :eq:`operatorsplitting`, :eq:`totalcomponent`, :eq:`equilibrium` and :eq:`kinetics`. These geochemical engines can be grouped into two groups:

* Amanzi's native geochemical engine: a set of basic geochemical capabilites included in Amanzi.
* External geochemical engines: any existing geochemical code that is coupled to Amanzi through using the Alquimia API.

Native Geochemistry
~~~~~~~~~~~~~~~~~~~

The native geochemical engine implements in Amanzi the methods to solve the equations eqs. :eq:`operatorsplitting`, :eq:`totalcomponent`, :eq:`equilibrium` and :eq:`kinetics` for the types of reactions outlined above. Specifically for kinetic reactions, the following mathematical formulation are available for reaction rate expressions (:math:`r_k`)

+-----------------------------------+----------------------------+-----------------------------------------------------------+
|  Kinetic Reaction Types           | Rate Expression Type       | Mathematical Formulation                                  |
+===================================+============================+===========================================================+
| Mineral Dissolution-Precipitation | Transition State Theory    | :math:`r_k = k \times A_s \times (1 -Q/K_s)`              |
+-----------------------------------+----------------------------+-----------------------------------------------------------+
| Aqueous Kinetics                  | First order dependence     | :math:`r_k = \lambda \times C_i`                          |
+-----------------------------------+----------------------------+-----------------------------------------------------------+

where :math:`k` and :math:`\lambda` are rate constants; :math:`A_s` is the reactive surface area of the mineral, a material property; :math:`Q` is the ion activity product; and :math:`K_s` is the solubility or equlibrium constant of the mineral reaction.

The reaction network is specified through a its own geochemical database file (typically with extension .bgd) that is specific to the problem at hand. The total concentrations of all components (:math:`U_i`) are specified in the Amanzi input file.  

The reader is kindly directed to the :ref:`Benchmark Testing/Chemistry <sec-benchmarks-chemistry>` section for examples.


Alquimia API
~~~~~~~~~~~~

Alquimia is an Application Programming Interface (API) that exposes the functionality of a geochemical engine to Amanzi. Alquimia does not perform any geochemical calculations itself. The geochemical engine is responsible for all geochemical calculations, and must provide a wrapper library that exactly conforms to the Alquimia API. Thus, the geochemical capabilities of Amanzi when using the Alquimia interface will depend on the geochemical engine of choice. That means that they can provide Amanzi with those capabilities or specific formulation not available in the native geochemical engine. 

Currently, the geochemical capabilities of the reactive transport code PFloTran 
(http://ees.lanl.gov/pflotran/ and https://bitbucket.org/pflotran) 
and CrunchFlow (https://bitbucket.org/crunchflow/crunchtope/downloads) are available within Amanzi through the Alquimia interface. 
These capabilities are described in the documentation for these packages. 
Some examples are available in the Amanzi documentation (:ref:`Benchmark Testing/Chemistry <sec-benchmarks-chemistry>`).

References
~~~~~~~~~~

.. bibliography:: /bib/ascem.bib
   :filter: docname in docnames
   :style:  alpha
   :keyprefix: chem-
