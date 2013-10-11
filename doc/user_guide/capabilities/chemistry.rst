Chemistry Models
----------------

The chemistry component of Amanzi models a set of geochemical processes that lead to the transformation of reactant to products according to a stoichiometric relationship. 

Stoichiometry
~~~~~~~~~~~~~

The stoichiometric relation for a generic reaction involving :math:`N_c` reactants (:math:`A_j`) and a product (:math:`A_i`) can be written as

.. math::
  \sum_{j=1}^{N_c} \nu_{ij}~A_j = A_i
  :label: stoichiometry

where :math:`\nu_{ij}` is the stoichiometric coefficients of the :math:`j-th` reactant in the :math:`i-th` reaction. 

Total component concentrations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the operator splitting approach employed by Amanzi (see below), the transport solver does not operate on the concentrations of each reactant (:math:`C_j`) and product (:math:`C_i`) species. Rather, transport operates on total concentrations of geochemical components (:math:`U_j`), defined as geochemical entities whose mass is not affected by geochemical equilibrium reactions. For a set of :math:`N_r` equilibrium reactions, the total mass of a component is the sum of the mass of a primary species (or component species) and the secondary species multiplied by the corresponding stoichiometric coefficients.

.. math::
  U_j = C_j + \sum_{i=1}^{N_r}{\nu_{ij}~C_i}
  :label: totalcomponent

Geochemical equilibrium
~~~~~~~~~~~~~~~~~~~~~~~

The local equilibrium assumption is applicable when reaction rates are much faster than transport rates, and therefore the concentrations of products and reactants can be expressed with the mass action law:

.. math::
  \log(K_i) = \sum_{j=1}^{N_c}{\nu_{ij}~log(\gamma_j \cdot C_j)} - log(\gamma_i \cdot C_i)
  :label: equilibrium

where :math:`K_i` is the equilibrium constant, and :math:`gamma_i` are the activity coefficients. 

Geochemical kinetics
~~~~~~~~~~~~~~~~~~~~

When local equilibrium cannot be assumed, reaction rates are calculated explicitly. E.g. for reaction :math::`i`:

.. math::
  r_i = f(C_i,C_j)
  :label: kinetics

Here, we choose to write the rate expression as a generic function but several forms of this function are used depending on the reaction of interest. 

Operator Splitting Approach
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The operator splitting approach involves solving the mass balance equation of each component for a given timestep in two steps. The first step involves the transport operators (see :doc:`transport`). Then, for the same timestep geochemical processes are solved:

.. math::
   \frac{\partial (\phi s_l U_i)}{\partial t} 
  = r_i
  :label: operatorsplitting

This equation is solved together with eqs. :eq:`totalcomponent`, :eq:`equilibrium` and :eq:`kinetics`.

Geochemical Engines
~~~~~~~~~~~~~~~~~~~

Taking advantage of the operator splitting approach, Amanzi offers a flexible approach to use geochemical engines that solve eqs. :eq:`operatorsplitting`, :eq:`totalcomponent`, :eq:`equilibrium` and :eq:`kinetics`. Amanzi has its own geochemical engine (native geochemistry) but it also uses the Alquimia API to couple existing geochemical engines to Amanzi. 

Native Geochemistry
~~~~~~~~~~~~~~~~~~~

The native geochemical engine can simulate a comprehensive set of types of reactions. These include:

* Equilibrium: 

  * Aqueous complexation
  * Surface complexation
  * Ion exchange

* Kinetics:

  * Mineral dissolution/preciptation
  * First-order decay (kinetic)

A detailed description of these reaction types and their formulation is provided in this user guide as a set of benchmarking examples (see :doc:`../testing/index`).


Alquimia Geochemistry
~~~~~~~~~~~~~~~~~~~~~

The geochemical capabilities when using the Alquimia interface depend on the simulation capabilites of the geochemical engine of choice. Currently, the geochemical capabilities of the reactive transport code PFloTran can be accessed in Amanzi through the Alquimia interface. For full details of PFloTran, see http://ees.lanl.gov/pflotran/ and https://bitbucket.org/pflotran/pflotran-dev/wiki/Home.
