/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*!

The `"aqueous equilibrium complexes`" section is a list of aqueous equiliblium reactions.
Each sublist is named after a secondary species and contains the following parameters:

* `"ion size parameter`" [double] is an empirical parameter that provides agreement 
  between measured activity coefficients and ionic strength. In theory, it is the diameter
  of the hydrated ion.

* `"charge`" [int] is the ion charge. The net charge of an ion is non-zero since the
  total number of electrons is unequal to the total number of protons.

* `"gram molecular weight`" [double] is amount of a molecular substance whose weight, 
  in grams, is numerically equal to the molecular weight of that substance. 

* `"reaction`" [string] is the equilibrium reaction equation.

* `"equilibrium constant`" [double] is the value of reaction quotient at chemical equilibrium.
  The reaction quotient is defined mathematically as the ratio of the activities of the product 
  species to that of the reactant species. The stoichiometric coefficients are taken into account 
  as exponents of activities.

In a non-isothermal simulation, the equlibrium constant is a function of temperature.
In such a case we have one additional parameter `"temperature`" and  parameter
`"equilibrium constant`" becomes a list that defines a piecewice linear function logK(T).
This is the `tabular function`_ with linear forms.

* `"T`" [Array(double)] is the array of temperature points.

* `"Keq`" [Array(double)] is the matching array of equilibroum constant values.

.. code-block:: xml

  <ParameterList name="thermodynamic database">
    <ParameterList name="aqueous equilibrium complexes">
      <ParameterList name="OH-">
        <Parameter name="ion size parameter" type="double" value="3.5"/>
        <Parameter name="charge" type="int" value="-1"/>
        <Parameter name="gram molecular weight" type="double" value="17.0073"/>
        <Parameter name="reaction" type="string" value="1.0 H2O  -1.0 H+"/>
        <Parameter name="equilibrium constant" type="double" value="13.9951"/>
        <Parameter name="temperature" type="double" value="298.15"/>
        <ParameterList name="equilibrium constant">
          <Parameter name="T" type="Array(double)" value="{273.15,  298.15,  333.15,  373.15,  423.15,  473.15,  523.15,  573.15}"/>
          <Parameter name="Keq" type="Array(double)" value="{14.9398, 13.9951, 13.0272, 12.2551, 11.6308, 11.2836, 11.1675, 11.3002}"/>
        </ParameterList>
      <ParameterList name="HCO3-">
        <Parameter name="ion size parameter" type="double" value="4.0"/>
        <Parameter name="charge" type="int" value="-1"/>
        <Parameter name="gram molecular weight" type="double" value="61.0171"/>
        <Parameter name="reaction" type="string" value="1.0 H2O  -1.0 H+  1.0 CO2(aq)"/>
        <Parameter name="equilibrium constant" type="double" value="6.3447"/>
        <Parameter name="temperature" type="double" value="298.15"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>

Below are a few examples of aqueous reaction.
Each line has five fields: reaction equation, logarithm of the equlibrium constant at 25C, 
ion size parameter, ion charge, and the gram molecular weight.

.. code-block:: txt

   AlHPO4+    =  1.0 Al+++ 1.0 HPO4--    -7.4       4.0     1.0   122.961
   CaCl+      =  1.0 Ca++  1.0 Cl-        0.6956    4.0     1.0    75.5307
   CaCl2(aq)  =  1.0 Ca++  2.0 Cl-        0.6436    3.0     0.0   110.9834
   CaCO3(aq)  =  1.0 Ca+2  1.0 CO3-2     -3.151     3.5     0.0    17.0073
   CaHCO3+    =  1.0 Ca++  1.0 HCO3-     -1.0467    4.0     1.0   101.0951
   CaHPO4(aq) =  1.0 Ca++  1.0 HPO4--    -2.7400    3.0     0.0   136.0573
   CO3--      = -1.0 H+    1.0 HCO3-     10.3288    4.5    -2.0    60.0092
   FeCl+      =  1.0 Fe++  1.0 Cl-        0.1605    4.0     1.0    91.2997 
   FeCl2(aq)  =  1.0 Fe++  2.0 Cl-        2.4541    3.0     0.0   126.752 
   FeCl4--    =  1.0 Fe++  4.0 Cl-       -1.9       4.0    -2.0   197.658 
   FeHCO3+    =  1.0 Fe++  1.0 HCO3-     -2.72      4.0     1.0   116.864 
   FeHPO4(aq) =  1.0 Fe++  1.0 HPO4--    -3.6       3.0     0.0   151.826 
   FeF+       =  1.0 Fe++  1.0 F-        -1.36      4.0     1.0    74.8454 
   FeSO4(aq)  =  1.0 Fe++  1.0 SO4--     -2.2       3.0     0.0   151.911 
   H2PO4-     =  1.0 H+    1.0 HPO4--    -7.2054    4.0    -1.0    96.9872
   H3PO4(aq)  =  2.0 H+    1.0 HPO4--    -9.3751    3.0     0.0    97.9952
   H2SO4(aq)  =  2.0 H+    1.0 SO4--      1.0209    3.0     0.0    98.0795 
   HCl(aq)    =  1.0 H+    1.0 Cl-       -0.67      3.0     0.0    36.4606 
   HNO3(aq)   =  1.0 H+    1.0 NO3-       1.3025    3.0     0.0    63.0129 
   HSO4-      =  1.0 H+    1.0 SO4--     -1.9791    4.0    -1.0    97.0715 
   KCl(aq)    =  1.0 K+    1.0 Cl-        1.4946    3.0     0.0    74.551
   KHPO4-     =  1.0 K+    1.0 HPO4--    -0.78      4.0    -1.0   135.078
   KSO4-      =  1.0 K+    1.0 SO4--     -0.8796    4.0    -1.0   135.162
   NaHCO3(aq) =  1.0 Na+   1.0 HCO3-     -0.1541    3.0     0.0    84.0069 
   NaCl(aq)   =  1.0 Na+   1.0 Cl-        0.777     3.0     0.0    58.4425 
   NaF(aq)    =  1.0 Na+   1.0 F-         0.9976    3.0     0.0    41.9882  
   NaHPO4-    =  1.0 Na+   1.0 HPO4--    -0.92      4.0    -1.0   118.969
   NaNO3(aq)  =  1.0 Na+   1.0 NO3-       1.044     3.0     0.0    84.9947
   NaSO4-     =  1.0 Na+   1.0 SO4--     -0.82      4.0    -1.0   119.053
   MgCO3(aq)  =  1.0 Mg+2  1.0 CO3-2     -2.928     3.5     0.0    17.0073
   OH-        =  1.0 H2O  -1.0 H+        13.9951    3.5    -1.0    17.00730
   P2O7----   = -1.0 H2O   2.0 HPO4--     3.7463    4.0    -4.0   173.9433
   PO4---     = -1.0 H+    1.0 HPO4--    12.3218    4.0    -3.0    94.9714
   UO2Cl+     =  1.0 Cl-   1.0 UO2++     -0.1572    4.0     1.0   305.48 
   UO2Cl2(aq) =  2.0 Cl-   1.0 UO2++      1.1253    3.0     0.0   340.933 
   UO2F+      =  1.0 F-    1.0 UO2++     -5.0502    4.0     1.0   289.026 
   UO2F2(aq)  =  2.0 F-    1.0 UO2++     -8.5403    3.0     0.0   308.024 
   UO2F3-     =  3.0 F-    1.0 UO2++    -10.7806    4.0    -1.0   327.023 
   UO2F4--    =  4.0 F-    1.0 UO2++    -11.5407    4.0    -2.0   346.021 
   UO2HPO4(aq)= 1.0 HPO4-- 1.0 UO2++     -8.4398    3.0     0.0   366.007 
   UO2NO3+    =  1.0 NO3-  1.0 UO2++     -0.2805    4.0     1.0   332.033 
   UO2SO4(aq) =  1.0 SO4-- 1.0 UO2++     -3.0703    3.0     0.0   366.091 
   UO2(SO4)2-- = 2.0 SO4-- 1.0 UO2++     -3.9806    4.0    -2.0   462.155 

   Al2(OH)2++++	 = -2.0 H+    2.0 Al+++   2.0 H2O         7.6902    5.5	   4.0    87.9778
   Al3(OH)4(5+)	 = -4.0 H+    3.0 Al+++   4.0 H2O        13.8803    6.0    5.0   148.9740
   Al(OH)2+      =  2.0 H2O   1.0 Al+++  -2.0 H+         10.5945    4.0    1.0    60.9962 
   Al(OH)3(aq)   =  3.0 H2O   1.0 Al+++  -3.0 H+         16.1577    3.0    0.0    78.0034 
   AlH2PO4++     =  1.0 Al+++ 1.0 H+      1.0 HPO4--     -3.1       4.5    2.0   123.969
   AlO2-         =  2.0 H2O   1.0 Al+++  -4.0 H+         22.8833    4.0   -1.0    58.9803
   AlOH++        =  1.0 H2O   1.0 Al+++  -1.0 H+          4.9571    4.5    2.0    43.9889 
   CaCO3(aq)     = -1.0 H+    1.0 Ca++    1.0 HCO3-       7.0017    3.0    0.0   100.0872
   CaH2PO4+      =  1.0 Ca++  1.0 H+      1.0 HPO4--     -1.4000    4.0    1.0   137.0652
   CaP2O7--      = -1.0 H2O   1.0 Ca++    2.0 HPO4--     -3.0537    4.0   -2.0   214.0213
   CaPO4-        = -1.0 H+    1.0 Ca++    1.0 HPO4--      5.8618    4.0   -1.0   135.0494
   CaOH+         = -1.0 H+    1.0 Ca++    1.0 H2O        12.8500    4.0    1.0    57.0853
   CO2(aq)       = -1.0 H2O   1.0 H+      1.0 HCO3-      -6.3447    3.0    0.0    44.0098
   H2P2O7--      = -1.0 H2O   2.0 H+      2.0 HPO4--    -12.0709    4.0   -2.0   175.9592
   H2S(aq)       =  2.0 H+    1.0 SO4--  -2.0 O2(aq)    131.329     3.0    0.0    34.0819 
   H3P2O7-       = -1.0 H2O   2.0 HPO4--  3.0 H+        -14.4165    4.0   -1.0   176.9671
   H4P2O7(aq)    = -1.0 H2O   2.0 HPO4--  4.0 H+        -15.9263    3.0    0.0   177.9751
   HAlO2(aq)     =  2.0 H2O   1.0 Al+++  -3.0 H+         16.4329    3.0    0.0    59.9883
   HCO3-         =  1.0 H2O  -1.0 H+      1.0 CO2(aq)    6.34470    4.0   -1.0    61.01710
   HO2-          =  1.0 H2O  -1.0 H+      0.5 O2(aq)     28.302     4.0   -1.0    33.0067 
   HP2O7---      = -1.0 H2O   1.0 H+      2.0 HPO4--     -5.4498    4.0   -3.0   174.9513
   HS-           =  1.0 H+    1.0 SO4--  -2.0 O2(aq)    138.317     3.5   -1.0    33.0739 
   Fe2(OH)2++++  =  1.0 H2O   2.0 Fe++    0.5 O2(aq)    -14.0299    5.5    4.0   145.709 
   FeCO3(aq)     =  1.0 Fe++ -1.0 H+      1.0 HCO3-       5.5988    3.0    0.0   115.856 
   FeH2PO4+      =  1.0 Fe++  1.0 H+      1.0 HPO4--     -2.7       4.0    1.0   152.834 
   Fe(OH)2(aq)   =  2.0 H2O   1.0 Fe++   -2.0 H+         20.6       3.0    0.0    89.8617 
   Fe(OH)3-      =  3.0 H2O   1.0 Fe++   -3.0 H+         31.0       4.0   -1.0   106.869 
   Fe(OH)4--     =  4.0 H2O   1.0 Fe++   -4.0 H+         46.0       4.0   -2.0   123.876 
   FeOH+         =  1.0 H2O   1.0 Fe++   -1.0 H+          9.5       4.0    1.0    72.8543 
   FeOH++        =  0.5 H2O   1.0 Fe++    0.25 O2(aq)    -6.3       4.5    2.0    72.8543 
   FePO4-        =  1.0 Fe++ -1.0 H+      1.0 HPO4--      4.3918    4.0   -1.0   150.818 
   KHSO4(aq)     =  1.0 K+    1.0 H+      1.0 SO4--      -0.8136    3.0    0.0   136.17
   KOH(aq)       =  1.0 H2O   1.0 K+     -1.0 H+         14.46      3.0    0.0    56.1056
   KP2O7---      = -1.0 H2O   1.0 K+      2.0 HPO4--      1.4286    4.0   -3.0   213.042
   MgOH+         =  1.0 H2O  -1.0 H+      1.0 Mg++      11.78510    4.0    1.0    41.3123
   NaCO3-        = -1.0 H+    1.0 HCO3-   1.0 Na+         9.8144    4.0   -1.0    82.9990
   NaOH(aq)      =  1.0 H2O   1.0 Na+    -1.0 H+         14.7948    3.0    0.0    39.9971
   NH3(aq)       =  1.5 H2O   0.5 N2(aq) -0.75 O2(aq)    58.2305    3.0    0.0    17.0306 
   UO2CO3(aq)    = -1.0 H+    1.0 HCO3-   1.0 UO2++       0.6634    3.0    0.0   330.037 
   UO2(CO3)2--   = -2.0 H+    2.0 HCO3-   1.0 UO2++       3.7467    4.0   -2.0   390.046 
   UO2(CO3)3---- = -3.0 H+    3.0 HCO3-   1.0 UO2++       9.4302    4.0   -4.0   450.055 
   UO2H2PO4+     =  1.0 H+    1.0 HPO4--  1.0 UO2++     -11.6719    4.0    1.0   367.015 
   UO2H3PO4++    =  2.0 H+    1.0 HPO4--  1.0 UO2++     -11.3119    4.5    2.0   368.023 
   UO2OH+        =  1.0 H2O  -1.0 H+      1.0 UO2++       5.2073    4.0    1.0   287.035 
   UO2PO4-       = -1.0 H+    1.0 HPO4--  1.0 UO2++      -2.0798    4.0   -1.0   364.999 
   UO2(OH)2(aq)  =  2.0 H2O  -2.0 H+      1.0 UO2++      10.3146    3.0    0.0   304.042 
   UO2(OH)3-     =  3.0 H2O  -3.0 H+      1.0 UO2++      19.2218    4.0   -1.0   321.05 
   UO2(OH)4--    =  4.0 H2O  -4.0 H+      1.0 UO2++      33.0291    4.0   -2.0   338.057 
   (UO2)2OH+++   =  1.0 H2O  -1.0 H+      2.0 UO2++       2.7072    5.0    3.0   557.063 
   (UO2)2(OH)2++ =  2.0 H2O  -2.0 H+      2.0 UO2++       5.6346    4.5    2.0   574.07 
   (UO2)3(OH)4++ =  4.0 H2O  -4.0 H+      3.0 UO2++      11.929     4.5    2.0   878.112 
   (UO2)3(OH)5+  =  5.0 H2O  -5.0 H+      3.0 UO2++      15.5862    4.0    1.0   895.12 
   (UO2)3(OH)7-  =  7.0 H2O  -7.0 H+      3.0 UO2++      31.0508    4.0   -1.0   929.135 
   (UO2)4(OH)7+  =  7.0 H2O  -7.0 H+      4.0 UO2++      21.9508    4.0    1.0   1199.16 
   UO2(H2PO4)(H3PO4)+ = 3.0 H+ 2.0 HPO4-- 1.0 UO2++     -22.7537    4.0    1.0   465.01 
   UO2(H2PO4)2(aq) =    2.0 H+ 2.0 HPO4-- 1.0 UO2++     -21.7437    3.0    0.0   464.002 
   Zn(OH)2(aq)   =  2.0 H2O  -2.0 H+      1.0 Zn++       17.3282    3.0    0.0    99.4047
   Zn(OH)3-      =  3.0 H2O  -3.0 H+      1.0 Zn++       28.8369    4.0   -1.0   116.41200
   Zn(OH)4--     =  4.0 H2O  -4.0 H+      1.0 Zn++       41.6052    4.0   -2.0   133.41940
   ZnOH+         =  1.0 H2O  -1.0 H+      1.0 Zn++        8.9600    4.0    1.0    82.39730

   Ca2UO2(CO3)3(aq) =  2.0 Ca++ -3.0 H+     3.0 HCO3-   1.0 UO2++        0.2864   4.0    0.0   530.215 
   CaUO2(CO3)3--    =  1.0 Ca++ -3.0 H+     3.0 HCO3-   1.0 UO2++        3.8064   4.0   -2.0   530.215 
   CH4(aq)          =  1.0 H2O   1.0 H+     1.0 HCO3-  -2.0 O2(aq)     144.141    3.0    0.0     0.0
   NaAlO2(aq)       =  2.0 H2O   1.0 Na+    1.0 Al+++  -4.0 H+          23.6266   3.0    0.0    81.9701
   NaHP2O7--        = -1.0 H2O   1.0 Na+    1.0 H+      2.0 HPO4--      -6.8498   4.0   -2.0   197.941
   NaHSiO3(aq)      =  1.0 H2O   1.0 Na+   -1.0 H+      1.0 SiO2(aq)     8.304    3.0    0.0   100.081
   Fe+++            = -0.5 H2O   1.0 Fe++   1.0 H+      0.25 O2(aq)     -8.49     9.0    3.0    55.847 
   Fe3(OH)4(5+)     =  2.5 H2O   3.0 Fe++  -1.0 H+      0.75 O2(aq)    -19.1699   6.0    5.0   235.57 
   Fe(OH)2+         =  1.5 H2O   1.0 Fe++  -1.0 H+      0.25 O2(aq)     -2.82     4.0    1.0    89.8617 
   Fe(OH)3(aq)      =  2.5 H2O   1.0 Fe++  -2.0 H+      0.25 O2(aq)      3.51     3.0    0.0   106.869 
   Fe(OH)4-         =  3.5 H2O   1.0 Fe++  -3.0 H+      0.25 O2(aq)     13.11     4.0   -1.0   123.876 
   FeCO3+           = -0.5 H2O   1.0 Fe++   1.0 HCO3-   0.25 O2(aq)     -7.8812   4.0    1.0   115.856 
   MgHCO3+          =  1.0 H2O  -1.0 H+     1.0 CO2(aq) 1.0 Mg++         5.309    4.0    1.0    85.3221
   N3-              =  0.5 H2O  -1.0 H+     1.5 N2(aq) -0.25 O2(aq)     77.7234   4.0   -1.0    42.0202 
   NH4+             =  1.5 H2O   1.0 H+     0.5 N2(aq) -0.75 O2(aq)     48.9895   2.5    1.0    18.0385 
   U+++             = -0.5 H2O   1.0 H+     1.0 UO2++  -0.75 O2(aq)     64.8028   5.0    3.0   238.029 
   U++++            = -1.0 H2O   2.0 H+     1.0 UO2++  -0.5 O2(aq)      33.949    5.5    4.0   238.029 
   UO2+             =  0.5 H2O  -1.0 H+     1.0 UO2++  -0.25 O2(aq)     20.0169   4.0    1.0   270.028 
   UO2OSi(OH)3+     =  2.0 H2O  -1.0 H+     1.0 SiO2(aq) 1.0 UO2++       2.4810   9.0    1.0   365.135
   (UO2)2CO3(OH)3-  =  3.0 H2O  -4.0 H+     1.0 HCO3-   2.0 UO2++       11.2229   4.0   -1.0   651.087 

   Fe(SO4)2- = -0.5 H2O  1.0 Fe++  1.0 H+      2.0 SO4--   0.25 O2(aq)   -11.7037   4.0   -1.0   247.974 
   FeCl++    = -0.5 H2O  1.0 Fe++  1.0 H+      1.0 Cl-     0.25 O2(aq)    -7.6792   4.5    2.0    91.2997 
   FeCl2+    = -0.5 H2O  1.0 Fe++  1.0 H+      2.0 Cl-     0.25 O2(aq)   -10.62     4.0    1.0   126.752 
   FeCl4-    = -0.5 H2O  1.0 Fe++  1.0 H+      4.0 Cl-     0.25 O2(aq)    -7.7      4.0   -1.0   197.658 
   FeF++     = -0.5 H2O  1.0 Fe++  1.0 H+      1.0 F-      0.25 O2(aq)   -12.6265   4.5    2.0    74.8454 
   FeF2+     = -0.5 H2O  1.0 Fe++  1.0 H+      2.0 F-      0.25 O2(aq)   -16.8398   4.0    1.0    93.8438 
   FeH2PO4++ = -0.5 H2O  1.0 Fe++  2.0 H+      1.0 HPO4--  0.25 O2(aq)   -12.66     4.5    2.0   152.834 
   FeHPO4+   = -0.5 H2O  1.0 Fe++  1.0 H+      1.0 HPO4--  0.25 O2(aq)   -18.67     4.0    1.0   151.826 
   FeNO3++   = -0.5 H2O  1.0 Fe++  1.0 H+      1.0 NO3-    0.25 O2(aq)    -9.49     4.5    2.0   117.852 
   FeSO4+    = -0.5 H2O  1.0 Fe++  1.0 H+      1.0 SO4--   0.25 O2(aq)   -10.4176   4.0    1.0   151.911 
   MgUO2(CO3)3-- = 3.0 H2O -6.0 H+ 3.0 CO2(aq) 1.0 Mg++    1.00 UO2++     23.9105   3.0   -2.0   500.0
   NH4SO4-   =  1.5 H2O  1.0 H+    0.5 N2(aq)  1.0 SO4--  -0.75 O2(aq)    57.2905   4.0   -1.0   114.102

*/

#ifndef AMANZI_CHEMISTRY_AQUEOUS_EQUILIBRIUM_COMPLEX_HH_
#define AMANZI_CHEMISTRY_AQUEOUS_EQUILIBRIUM_COMPLEX_HH_

#include <vector>

#include "SecondarySpecies.hh"

namespace Amanzi {
namespace AmanziChemistry {

// forward declarations from chemistry
class MatrixBlock;

class AqueousEquilibriumComplex : public SecondarySpecies {
 public:
  AqueousEquilibriumComplex() : SecondarySpecies(){};
  AqueousEquilibriumComplex(int id,
                            const std::string& name,
                            const Teuchos::ParameterList& plist,
                            const std::vector<Species>& primary_species);
  ~AqueousEquilibriumComplex(){};

  // update molalities
  virtual void Update(const std::vector<Species>& primary_species, const Species& water_species);

  // add stoichiometric contribution of complex to total
  virtual void AddContributionToTotal(std::vector<double>* total);

  // add derivative of total with respect to free-ion to dtotal
  virtual void
  AddContributionToDTotal(const std::vector<Species>& primary_species, MatrixBlock* dtotal);

  void Display(const Teuchos::Ptr<VerboseObject> vo) const;
  void DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const;
};

} // namespace AmanziChemistry
} // namespace Amanzi

#endif
