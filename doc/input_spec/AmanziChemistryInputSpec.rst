=============================================
Amanzi Chemistry Database Input Specification
=============================================

.. contents:: **Table of Contents**


Overview
========

Documentation of the `"simple`" database format used by the amanzi chemistry library.

Importing thermodynamic data into the chemistry module using the `"simple"` (file extension `"bgd"`) format requires the user to explicitly specify all the species and reactions for the problem. There is no basis switching or automatic species and reaction selection. 

In this format, any line starting with a `"#"` or space character is a comment. 

Data is separated into sections, where each section of the file is starts with a line containing `<Section Name`. The less than character, `<`, should be the first character on the line and there is no space between the character and the section name. NOTE(bandre): It would be a trivial code change to get rid of the spacing requirements....

  Sections should be ordered so that the primary species, minerals, and exchange sites come before any reactions using those species.

  The valid section names are:

  * `Primary Species`
  * `Aqueous Equilibrium Complexes`
  * `Minerals`
  * `Mineral Kinetics`
  * `General Kinetics`
  * `Ion Exchange Sites`
  * `Ion Exchange Complexes`
  * `Surface Complex Sites`
  * `Surface Complexes`
  * `Sorption Isotherms`

Each line within a section is semi-colon delimited. White space around values is removed.

Database Format
===============

Primary Species
~~~~~~~~~~~~~~~

A primary species line contains, in order: the primary species name, debye-huckel a0 parameter, charge and gram molecular weight.

Example::

  <Primary Species
  # name      ; debye-huckel a0 ; charge ; GMW [grams/mole]
  Al+++       ; 9.0             ;    3.0 ; 26.9815



Aqueous Equilibrium Complexes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An aqueous equilibrium complex, or secondary species, line contains a reaction and data for the reaction on a single line. The reaction is written as `"product species = reactants....`" The coefficient of the product aqueous complex is assumed to be 1.0, and the reactants must be primary species. The equilibrium constant is for a fixed temperature of 25C.
  
Example::

  <Aqueous Equilibrium Complexes
  # name               =  coeff primary_name  coeff primary_name ... ; log10(Keq) 25C ; debye-huckel a0 ; charge ; GMW [grams/mole]
  OH-                  =  1.0 H2O  -1.0 H+                           ;    13.9951     ;   3.5           ;  -1.0  ;  17.0073


Aqueous Kinetics
~~~~~~~~~~~~~~~~

Example::

  <General Kinetics
    1.00 Pu_238 <->   1.00 U_234 ;   1.00 Pu_238 ;   2.50000E-10 ;   1.00 U_234 ;   0.00000E+00
    1.00 U_234 <->   1.00 Th_230 ;   1.00 U_234 ;   8.90000E-14 ;   1.00 Th_230 ;   0.00000E+00
    1.00 Th_230 <->   1.00 Ra_226 ;   1.00 Th_230 ;   2.90000E-13 ;   1.00 Ra_226 ;   0.00000E+00
    1.00 Ra_226 <->   1.00 Pb_210 ;   1.00 Ra_226 ;   1.40000E-11 ;   1.00 Pb_210 ;   0.00000E+00


Minerals
~~~~~~~~
Minerals and other complexes follow the same convention as aqueous equilibrium complexes, with additional data as needed. This section contains all minerals present in the system during the simulation, including those initiall not present but that may precipitate later. They are used for calculating saturation states, but not equilibrium or kinetic calculations.

Example::

  <Minerals
  # name               =  coeff primary_name  coeff primary_name ... ; log10(Keq) 25C ; GMW      ; molar volume [cm^2/mol] ; SSA [m^2/g] 


Mineral Kinetics
~~~~~~~~~~~~~~~~

The mineral kinetics section lists the name of a mineral found in the mineral section, the type of rate law, and rate parameters for that law. Currently only the `TST` rate law is implemented. The keywords "log10_rate_constant" and "moles_m2_sec" must be present in the line, but no unit conversion are currently preformed. The modifying primary species terms follow the rate constant, along with their exponent coefficients.

Example::

  <Mineral Kinetics
  # name               ; TST ; log10_rate_constant double     moles_m2_sec ; primary_name coeff ....

Equilibrium Surface Complexation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Surface complex sites are listed by name and surface density. At this time all surface sites are assumed to occur on the bulk material rather than a specific mineral.

Example::

  <Surface Complex Sites 
  # name  ; site_density [moles sites / m^3 bulk]
  >FeOH_w ;   3.70000E-04

Surface Complexes....

Example::

  <Surface Complexes
  # name     =  coeff surface site  coeff primary_name ... ; log10(Keq) 25C ; charge 
  >FeOH2+_w  =   1.00 >FeOH_w  1.00 H+ ;  -7.18000E+00 ;   1.00
  >FeO-_w    =   1.00 >FeOH_w -1.00 H+ ;   8.82000E+00 ;  -1.00
  >FeOHZn+_w =   1.00 >FeOH_w -1.00 H+  1.00 Zn++ ;   2.32000E+00 ;   1.00


Ion Exchange
~~~~~~~~~~~~

Example::

  <Ion Exchange Sites
  # name ; charge ; mineral (ignored)
  X- ; -1.0 ; Kaolinite

  <Ion Exchange Complexes
  # name   = coeff primar_name  coeff site_name  ; exchange coeff
  NaX = 1.0 Na+  1.00 X-  ;   1.00000E+00
  UO2X2 = 1.0 UO2++  2.00 X-  ;   2.23872E-01
  CaX2 = 1.0 Ca++  2.00 X-  ;   3.16228E-01
  AlX3 = 1.0 Al+++  3.00 X-  ;   1.71133E+00
  HX = 1.0 H+  1.00 X-  ;   2.51189E-02



Sorption Isotherms
~~~~~~~~~~~~~~~~~~

Example::

  <Isotherms
  # Primary Species Name ; linear ; Kd
  Pu_238 ; linear ;   2.00000E+07
  # Primary Species Name ; langmuir ; Kd  langmuir_b
  U_234 ; langmuir ;   5.00000E+06    1.0
  # Primary Species Name ; langmuir ; Kd  freundlich_n
  Th_230 ; freundlich ;   1.00000E+07  1.0



Example
=======

Below is an example of a `"simple"` database file for a five component uranium problem with mineral dissolution and surface complexation:

  
::

  <Primary Species
  # name               ; debye-huckel a0 ; charge ; GMW
  Al+++                ;   9.0 ;   3.0 ;  26.9815
  H+                   ;   9.0 ;   1.0 ;   1.0079
  HPO4--               ;   4.0 ;  -2.0 ;  95.9793
  SiO2(aq)             ;   3.0 ;   0.0 ;  60.0843
  UO2++                ;   4.5 ;   2.0 ;  270.028

  <Aqueous Equilibrium Complexes
  # name               =  coeff primary_name  coeff primary_name  ; log10(Keq) 25C ; debye-huckel a0 ; charge ; GMW      
  OH-                  =  1.0 H2O  -1.0 H+  ;    13.9951 ;   3.5 ;  -1.0 ;  17.0073 
  AlOH++               =  1.0 H2O  1.0 Al+++  -1.0 H+  ;     4.9571 ;   4.5 ;   2.0 ;  43.9889 
  Al(OH)2+             =  2.0 H2O  1.0 Al+++  -2.0 H+  ;    10.5945 ;   4.0 ;   1.0 ;  60.9962 
  Al(OH)3(aq)          =  3.0 H2O  1.0 Al+++  -3.0 H+  ;    16.1577 ;   3.0 ;   0.0 ;  78.0034 
  Al(OH)4-             =  4.0 H2O  1.0 Al+++  -4.0 H+  ;    22.8833 ;   4.0 ;  -1.0 ;  95.0107 
  UO2OH+               =  1.0 H2O  -1.0 H+  1.0 UO2++  ;     5.2073 ;   4.0 ;   1.0 ;  287.035 
  UO2(OH)2(aq)         =  2.0 H2O  -2.0 H+  1.0 UO2++  ;    10.3146 ;   3.0 ;   0.0 ;  304.042 
  UO2(OH)3-            =  3.0 H2O  -3.0 H+  1.0 UO2++  ;    19.2218 ;   4.0 ;  -1.0 ;   321.05 
  UO2(OH)4--           =  4.0 H2O  -4.0 H+  1.0 UO2++  ;    33.0291 ;   4.0 ;  -2.0 ;  338.057 
  (UO2)2OH+++          =  1.0 H2O  -1.0 H+  2.0 UO2++  ;     2.7072 ;   5.0 ;   3.0 ;  557.063 
  (UO2)2(OH)2++        =  2.0 H2O  -2.0 H+  2.0 UO2++  ;     5.6346 ;   4.5 ;   2.0 ;   574.07 
  (UO2)3(OH)4++        =  4.0 H2O  -4.0 H+  3.0 UO2++  ;     11.929 ;   4.5 ;   2.0 ;  878.112 
  (UO2)3(OH)5+         =  5.0 H2O  -5.0 H+  3.0 UO2++  ;    15.5862 ;   4.0 ;   1.0 ;   895.12 
  (UO2)3(OH)7-         =  7.0 H2O  -7.0 H+  3.0 UO2++  ;    31.0508 ;   4.0 ;  -1.0 ;  929.135 
  (UO2)4(OH)7+         =  7.0 H2O  -7.0 H+  4.0 UO2++  ;    21.9508 ;   4.0 ;   1.0 ;  1199.16 
  UO2(H2PO4)(H3PO4)+   =  3.0 H+  2.0 HPO4--  1.0 UO2++  ;   -22.7537 ;   4.0 ;   1.0 ;   465.01 
  UO2(H2PO4)2(aq)      =  2.0 H+  2.0 HPO4--  1.0 UO2++  ;   -21.7437 ;   3.0 ;   0.0 ;  464.002 
  UO2HPO4(aq)          =  1.0 HPO4--  1.0 UO2++  ;    -8.4398 ;   3.0 ;   0.0 ;  366.007 
  UO2H2PO4+            =  1.0 H+  1.0 HPO4--  1.0 UO2++  ;   -11.6719 ;   4.0 ;   1.0 ;  367.015 
  UO2H3PO4++           =  2.0 H+  1.0 HPO4--  1.0 UO2++  ;   -11.3119 ;   4.5 ;   2.0 ;  368.023 
  UO2PO4-              =  -1.0 H+  1.0 HPO4--  1.0 UO2++  ;    -2.0798 ;   4.0 ;  -1.0 ;  364.999 

  <Minerals
  # name               =  coeff primary_name  coeff primary_name  ; log10(Keq) 25C ; GMW      ; molar volume [cm^2/mol] ; SSA [m^2/g] 

  Kaolinite            =  5.00 H2O  2.00 Al+++  -6.00 H+  2.00 SiO2(aq)  ;     6.8101 ;   258.16 ;    99.52 ;   1.0 
  Quartz               =  1.00 SiO2(aq)  ;    -3.9993 ;  60.0843 ;   22.688 ;   1.0 
  (UO2)3(PO4)2.4H2O    =  4.00 H2O  -2.00 H+  2.00 HPO4--  3.00 UO2++  ;   -27.0349 ;  1072.09 ;    500.0 ;   1.0 

  <Mineral Kinetics
  # name               ; TST ; log10_rate_constant double     moles_m2_sec 

  Kaolinite            ; TST ; log10_rate_constant    -16.699 moles_m2_sec 
  Quartz               ; TST ; log10_rate_constant      -18.0 moles_m2_sec 
  (UO2)3(PO4)2.4H2O    ; TST ; log10_rate_constant      -10.0 moles_m2_sec 

  <Surface Complex Sites
  # name               ; surface_density

  >FeOH                ; 6.3600E-03
  >AlOH                ; 6.3600E-03
  >SiOH                ; 6.3600E-03

  <Surface Complexes
  # name               =  coeff surface site  coeff primary_name  ; log10(Keq) 25C ; charge 

  >SiOUO3H3++          =  1.0 >SiOH  1.0 H2O  1.0 UO2++  ;       5.18 ;   2.0 
  >SiOUO3H2+           =  1.0 >SiOH  1.0 H2O  -1.0 H+  1.0 UO2++  ;       5.18 ;   1.0 
  >SiOUO3H             =  1.0 >SiOH  1.0 H2O  -2.0 H+  1.0 UO2++  ;       5.18 ;   0.0 
  >SiOUO3-             =  1.0 >SiOH  1.0 H2O  -3.0 H+  1.0 UO2++  ;      12.35 ;  -1.0 
  >SiOUO2(OH)2-        =  1.0 >SiOH  2.0 H2O  -3.0 H+  1.0 UO2++  ;      12.35 ;  -1.0 
  >FeOHUO3             =  1.0 >FeOH  1.0 H2O  -2.0 H+  1.0 UO2++  ;       3.05 ;   0.0 
  >FeOHUO2++           =  1.0 >FeOH  1.0 UO2++  ;      -6.63 ;   2.0 
  >AlOUO2+             =  1.0 >AlOH  -1.0 H+  1.0 UO2++  ;      -3.13 ;   1.0 


 
