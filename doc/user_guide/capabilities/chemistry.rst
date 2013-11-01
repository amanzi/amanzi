Chemistry Models
----------------

The chemistry component of Amanzi models a set of geochemical processes that lead to the transformation of reactant species to product species. Amanzi can simulate a comprehensive set of reaction types.  These reactions can be classified into homogeneous if they occur within the aqueous phase or heterogeneous if in addition to the aqueous phase they involve the solid phase. These reactions can also be classified into equilibrium or kinetic depending on whether local equilibrium can be assumed. The following table summarizes the reactions available by type.

+--------------------+----------------------------+---------------------------------------+
| *Reaction Types*   | Homogeneous                | Heterogeneous                         |
+--------------------+----------------------------+---------------------------------------+
| Equilibrium        | **Aqueous Complexation**   | **Surface Complexation**,             |
|                    |                            | **Ion Exchange**                      |
+--------------------+----------------------------+---------------------------------------+
| Kinetic            | **First-order decay**      | **Mineral Dissolution-Precipitation** |
+--------------------+----------------------------+---------------------------------------+

Stoichiometry
~~~~~~~~~~~~~

The transformation from reactant to product species occurs according to a stoichiometric relationship. The stoichiometry of a generic reaction involving :math:`N_c` reactants (:math:`A_j`) and a product (:math:`A_i`) can be written as

.. math::
  \sum_{j=1}^{N_c} \nu_{ij}~A_j = A_i
  :label: stoichiometry

where :math:`\nu_{ij}` is the stoichiometric coefficients of reactant *j* in reaction *i*. 

Geochemical equilibrium
~~~~~~~~~~~~~~~~~~~~~~~

The local equilibrium assumption is applicable when reaction rates are much faster than transport rates, and therefore the relationship between the concentrations of products and reactants can be expressed with the law of mass action:

.. math::
  \log(K_i) = \sum_{j=1}^{N_c}{\nu_{ij}~log(\gamma_j \cdot C_j)} - log(\gamma_i \cdot C_i)
  :label: equilibrium

where :math:`K_i` is the equilibrium constant, :math:`C_i` and :math:`C_j` are the concentrations of reactant and product species, and :math:`\gamma_i` and :math:`\gamma_j` are the corresponding activity coefficients. 

Total component concentrations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Geochemical equilibrium between aqueous species makes it possible to calculate the concentration of a set of secondary species as a function of a set of primary species. Further, one can define a set of components, one for each primary species, the mass of which are not affected by geochemical aqueous equilibrium reactions. The total mass of a component  (:math:`U_j`) is defined as the sum of the mass of a primary species and the mass of a number of secondary species multiplied by the corresponding stoichiometric coefficient:

.. math::
  U_j = C_j + \sum_{i=1}^{N_r}{\nu_{ij}~C_i}
  :label: totalcomponent

where :math:`N_r` is the number of aqueous equilibrium reactions.

This approach is widely used in geochemical models ([Lichtner1985]_, [Yeh-Tripathi1989]_, [Steefel-Lasaga1994]_, [Steefel-MacQuarrie1996]_) to reduce the size of the system to solve, and more specifically the number of transport equations.

Transport and Chemistry Coupling: Operator Splitting Approach
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Amanzi uses the so-called Operator Spliting Aproach to couple transport and reaction processes [Yeh-Tripathi1989]_. In the operator splitting approach, transport and chemistry are solved sequentially over the same time step, first transport and then chemistry. One advantage of this approach is that transport equations can be solved using linear solvers. Additionally, the non-linear chemical equations can be solved on a cell-by-cell basis in what is an embarassingly parallel problem. To access these advantages, the transport solver does not operate on the concentrations of each reactant (:math:`C_j`) and product species (:math:`C_i`). Rather, transport operates on the total component concentrations (:math:`U_j`) as defined above :eq:`totalcomponent`. Please see :doc:`transport` for details of the transport capabilities of Amanzi. Following transport, the chemical processes affecting the mass balance of all components are solved over the same timestep:

.. math::
   \frac{\partial (\phi s_l U_j)}{\partial t} 
  = r_j
  :label: operatorsplitting

This set of equations is solved together with eqs. :eq:`totalcomponent` and :eq:`equilibrium` to know the concentration of the individual species, as well as with the equations that describe the kinetic reaction rates :math:`r_j`.

Geochemical kinetics
~~~~~~~~~~~~~~~~~~~~

When local equilibrium cannot be assumed, reaction rates need to be calculated explicitly. The reaction rates are in general a non-linear function of concentrations of the geochemical species. Additionally, for heterogeneous reactions, the rates can be a function of material properties (symbolically: :math:`\Psi`), e.g. reactive surface area in the case of mineral dissolution-precipitation. The particular expression depends on the reaction type. For example, mathematical formulations include the transition state theory rate law for mineral dissolution-precipitation or a first order dependence on concentration for radioactive decay. For the sake of brevitiy, the rate expression for reaction *k* is provided here only as a generic function and the reader is directed to the chemistry benchmarking examples for specific mathematical formulations (see :doc:`../benchmarking/chemistry/index`):

.. math::
  r_k = f(C_i,\Psi)
  :label: kinetics

Geochemical Engines
~~~~~~~~~~~~~~~~~~~

Taking advantage of the operator splitting approach, Amanzi offers a flexible approach to use geochemical engines that solve eqs. :eq:`operatorsplitting`, :eq:`totalcomponent`, :eq:`equilibrium` and :eq:`kinetics`. Amanzi has its own geochemical engine (native geochemistry) but it also uses the Alquimia API to couple existing geochemical engines to Amanzi. 

Native Geochemistry
~~~~~~~~~~~~~~~~~~~

The native geochemical engine implements in Amanzi the methods to solve the equations eqs. :eq:`operatorsplitting`, :eq:`totalcomponent`, :eq:`equilibrium` and :eq:`kinetics` for the types of reactions outlined above. The reaction network is specified through a its own geochemical database file (typically with extension .bgd) that is specific to the problem at hand. The total concentrations of all components (:math:`U_i`) are specified in the Amanzi input file.

Alquimia Geochemistry
~~~~~~~~~~~~~~~~~~~~~

The geochemical capabilities and specific implementations for chemistry when using the Alquimia interface depend on the geochemical engine of choice. Currently, the geochemical capabilities of the reactive transport code PFloTran can be accessed in Amanzi through the Alquimia interface. For full details of PFloTran, see http://ees.lanl.gov/pflotran/ and https://bitbucket.org/pflotran/pflotran-dev/wiki/Home.


References
~~~~~~~~~~

.. [Lichtner1985] Lichtner, P. (1985), Continuum model for simultaneous chemical-reactions and mass-transport in hydrothermal systems, Geochim. Cosmochim. Acta, 49(3), 779–800, doi:10.1016/0016-7037(85)90172-3.

.. [Steefel-Lasaga1994] Steefel, C. I., and A. C. Lasaga (1994), A coupled model for transport of Hydrothermal fluxes of major elements, Juan de Fuca flank 1755 multiple chemical species and kinetic precipitation/dissolution reactions with application to reactive flow in single phase hydrothermal systems, Am. J. Sci., 294, 529–592.

.. [Steefel-MacQuarrie1996] Steefel, C. I., and K. MacQuarrie (1996), Approaches to modeling of reactive transport in porous media, in Reactive Transport In Porous Media, Rev. in Min., vol. 34, edited by P. C. Lichtner, C. I. Steefel, and E. H. Oelkers, pp. 83–129, Min. Soc. Am., Washington, D.C.

.. [Yeh-Tripathi1989] Yeh, G. T., and V. S. Tripathi (1989), A critical evaluation of recent developments in hydrogeochemical transport models of reactive multichemical components, Water Resour. Res., 25, 93–108.
