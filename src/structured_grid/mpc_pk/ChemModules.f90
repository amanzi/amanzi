!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>                                                                     >
!>                                                                     >
!>       Modules related to reactive biogeochemical transport          >
!>                                                                     >
!>                                                                     >
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Batch_PhysiChem_Conditions
! 
         SAVE
! 
! ----------
! ...... Double precision variables     
! ----------
! 
         REAL(KIND = 8) :: tc2          ! Temperature, oC
         REAL(KIND = 8) :: tk2          ! Temperature, K
         REAL(KIND = 8) :: tmpmin       ! Maximum temperature of log K fit range
         REAL(KIND = 8) :: tmpmax       ! Minimum temperature of log K fit range
         REAL(KIND = 8) :: phi2         ! Porosity
         REAL(KIND = 8) :: sg2          ! Current gas saturation          
         REAL(KIND = 8) :: sl2          ! Current liquid saturation
         REAL(KIND = 8) :: phisl2       ! = phi2*sl2
         REAL(KIND = 8) :: phisl2_old   ! phisl2 at the previous time step
         REAL(KIND = 8) :: ph2          ! pH value
         REAL(KIND = 8) :: aw2          ! Water activity
         REAL(KIND = 8) :: str          ! ionic strength
         REAL(KIND = 8) :: Pt2          ! Total pressure, bar
         REAL(KIND = 8) :: densw2       ! Density of water, kg/m^3 
         REAL(KIND = 8) :: a_fmr2       ! Modified active fracture area for reaction  
         REAL(KIND = 8) :: fugcoeCO22   ! CO2 fugacity coefficient from flow module  
         REAL(KIND = 8) :: Rh2          ! Relative humidity obtained from EOS
!
         REAL(KIND = 8) :: sumsalts, &  ! Sum of salts weights in kg per kg water
                           dliq,     &  ! Water density in g/cc (kg/l) = densw2E-3
                           vliq,     &  ! Volume liquid (liters) 
                           factw        ! Conversion factor = kg h2o liq/liter liquid
! 
      END MODULE Batch_PhysiChem_Conditions
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE MaxChem_Dimensions
! 
! 
         SAVE
! 
! ----------
! ...... Integer Maximum dimension variables for geochemical system     
! ----------
! 
         INTEGER(KIND=8)  :: MaxNpri     ! Maximum number of primary species
         INTEGER(KIND=8)  :: MaxNtrx     ! Maximum Number of kinetic reactions among primary species
         INTEGER(KIND=8)  :: MaxNaqx     ! Maximum number of secondary species
         INTEGER(KIND=8)  :: MaxNaqt     ! = MaxNpri + MaxNaqx
         INTEGER(KIND=8)  :: MaxNmin     ! Maximum number of equilibrium and kinetic minerals
         INTEGER(KIND=8)  :: MaxNss      ! Maximum number of solid solutions 
         INTEGER(KIND=8)  :: MaxCpss     ! Maximum number of solid solution endmembers  
         INTEGER(KIND=8)  :: MaxNgas     ! Maximum number of gasesous species
         INTEGER(KIND=8)  :: MaxNads     ! Maximum number of surface complexes
         INTEGER(KIND=8)  :: MaxNexc     ! Maximum number of exchangeable cations 
         INTEGER(KIND=8)  :: MXsites     ! Maximun number of exchangeable sites
         INTEGER(KIND=8)  :: MaxNkdd     ! Maximum number of species with linear Kd adsorption or/and decay
!
         INTEGER(KIND=8)  :: MaxNnr      ! = MaxNpri + MaxNmin + MaxNgas + 2
!
         INTEGER(KIND=8)  :: MaxNiwtype  !  Maximum number of initial water types 
         INTEGER(KIND=8)  :: MaxNbwtype  !  Maximum number of boundary (injection) water types
!
         INTEGER(KIND=8)  :: MaxNmtype   !  Maximum number of initial mineral zones
         INTEGER(KIND=8)  :: MaxNgtype   !  Maximum number of initial gas zones
         INTEGER(KIND=8)  :: MaxNppzon   !  Maximum number of prosity-permeability zones
         INTEGER(KIND=8)  :: MaxNdtype   !  Maximum number of surface complex zones 
         INTEGER(KIND=8)  :: MaxNkdtype  !  Maximum number of lineral Kd adsorption zones 
         INTEGER(KIND=8)  :: MaxNxtype   !  Maximum number of exchange zones     
!
         INTEGER(KIND=8)  :: MaxNtmp     ! Maximum number of temperature points in the database
         INTEGER(KIND=8)  :: MaxNpair    ! Maximum of ion pairs for using Pitzer ion-interaction model
!
         INTEGER(KIND=8)  :: MaxNum_Elem        !                                 !!!!!! Tempor..!!!! Delete leter
         INTEGER(KIND=8)  :: MaxNum_Conx       !                                 !!!!!! Tempor..!!!! Delete leter
!
! 
      END MODULE MaxChem_Dimensions
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************!
!
!
      MODULE Chem_Dimensions
!
! 
         SAVE
! 
! ----------
! ...... Integer dimension variables for geochemical system     
! ----------
! 
         INTEGER(KIND=8)  :: npri     ! Number of primary species
         INTEGER(KIND=8)  :: ntrx     ! Number of kinetic reactions among primary species
         INTEGER(KIND=8)  :: naqx     ! Number of secondary species
         INTEGER(KIND=8)  :: naqt     ! =npri + naqx
         INTEGER(KIND=8)  :: nmequ    ! Number of equlibrium minerals
         INTEGER(KIND=8)  :: nmkin    ! Number of kinetic minerals
         INTEGER(KIND=8)  :: nmin     ! Total number of minerals
         INTEGER(KIND=8)  :: nss      ! Number of solid solutions 
         INTEGER(KIND=8)  :: ngas     ! Number of gasesous species
         INTEGER(KIND=8)  :: nads     ! Number of surface complexes
         INTEGER(KIND=8)  :: nexc     ! Number of exchangeable cations 
         INTEGER(KIND=8)  :: NXsites  ! Number of exchangeable sites
         INTEGER(KIND=8)  :: nkdd     ! number of species with linear Kd adsorption or/and decay
         INTEGER(KIND=8)  :: ntot     ! = npri + nmin +ngas    
!
         INTEGER(KIND=8)  :: ii       !  Counter for water zones 
         INTEGER(KIND=8)  :: iwtype   !  Initial water type index (for CHDUMP subroutine)
         INTEGER(KIND=8)  :: niwtype  !  Number of initial water types 
         INTEGER(KIND=8)  :: nbwtype  !  Number of boundary (injection) water types
!
         INTEGER(KIND=8)  :: nmtype   !  Number of initial mineral zones
         INTEGER(KIND=8)  :: ngtype   !  Number of initial gas zones
         INTEGER(KIND=8)  :: nppzon   !  Number of prosity-permeability zones
         INTEGER(KIND=8)  :: ndtype   !  Number of surface complex zones 
         INTEGER(KIND=8)  :: nkdtype  !  Number of lineral Kd adsorption zones 
         INTEGER(KIND=8)  :: nxtype   !  Number of exchange zones     
!
! 
      END MODULE Chem_Dimensions
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Species_Index
!
! 
         SAVE
! 
! ----------
! ...... Integer variables for special species     
! ----------
! 
         INTEGER(KIND=2)  :: ne        ! Index of e (electron)
         INTEGER(KIND=2)  :: nw        ! Index of species H2O
         INTEGER(KIND=2)  :: nh        ! Index of species H+
         INTEGER(KIND=2)  :: noh       ! Index of species OH-
         INTEGER(KIND=2)  :: no2aq     ! Index of species O2(aq)
         INTEGER(KIND=2)  :: nd        ! Index of primary surface species
         INTEGER(KIND=2)  :: nx        ! Index of master exchangeable species
         INTEGER(KIND=2)  :: nco2      ! Index of aqueous CO2 species
         INTEGER(KIND=2)  :: nco2g     ! Index of gaseous CO2 species
         INTEGER(KIND=2)  :: nh2       ! Index of H2(aq) among all aqueous species 
         INTEGER(KIND=2)  :: nh2g      ! Index of gaseous H2 species
!
! ----------
! ...... Integer control variables     
! ----------
!
          INTEGER(KIND=1)  :: ico2gt0   ! =1: initial Pco2>0         
          INTEGER(KIND=1)  :: ih2gt0    ! =1: initial Ph2>0         
          INTEGER(KIND=1)  :: iaqxs     ! =1: define secondary species in CHEMICAL.INP
! 
! 
      END MODULE Species_Index
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Chem_Constants
!
! 
         SAVE
! 
! ----------
! ...... Double precision parameters
! ----------
! 
         REAL(KIND = 8), PARAMETER :: gc      = 0.0831456d0  ! Gas constant in liter*bar/K/mole
         REAL(KIND = 8), PARAMETER :: rmh2o   = 55.50868d0   ! Moles of H2O per kg water
         REAL(KIND = 8), PARAMETER :: faraday = 96485.d0     ! Faraday constant
         REAL(KIND = 8), PARAMETER :: rgas    = 8.3145d0     ! Gas constant
         REAL(KIND = 8), PARAMETER :: rgasm   = 8.3144d-3    ! Gas constant
         REAL(KIND = 8), PARAMETER :: epsi    = 7.08d-11     ! epsi used for double-layer adsorption model
!
         REAL(KIND = 8), PARAMETER :: & ! Areas and permeability changes calculations                  
              ppi     = 3.1415926536d0,       &
              pi3     = 9.42477796076938d0,   &
              pi4     = 12.5663706143592d0,   &
              pi6     = 18.8495559215388d0,   &
              pid128  = 2.45436926061703d-2,  &
              tfdpi   = 0.238732414637843d0 
! 
! 
      END MODULE Chem_Constants
! 
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Aqueous_Species
!
! 
         SAVE
!
         REAL(KIND = 8) :: xh2o    ! H2O mole fraction
         REAL(KIND = 8) :: denswa  ! Liquid (Water+solute) density (kg/l)
! 
! ----------
! ...... Derived-type variables     
! ----------
! 
         TYPE Primary_Species
!
            CHARACTER(LEN = 20) :: napri  ! Name of primary species such as 'H+'
!  
            REAL(KIND = 8) :: zp   ! The ion electric charge of primary species
            REAL(KIND = 8) :: cp   ! Concentrations of primary species          
            REAL(KIND = 8) :: u2   ! Total dissolved concentrations of the primary species
!
            INTEGER(KIND = 1) :: NoTrans   ! NoTrans=1: not subject to transport, =0 for transport
            INTEGER(KIND = 1) :: icon   ! Flag indicating the type of constraint controlling the solute content 
!                   that is given under u2:
!                      =1 the concentration of the species is constrained by the total dissolved
!                         concentration U2, except for H2O which is assumed to have a value of 1.0.
!                      =2 the concentration of the species calculated through charge balance.
!                      =3 the activity of the species is fixed during the inialization. 
!                         In this case, the following two variables: CGUESS = CTOT = the 
!                         fixed activity. 
!                   For example, if the pH of an initial water is fixed at 7, 
!                   we set CGUESS = CTOT = 10-7 for H+ activity. 
!
            INTEGER(KIND = 2):: kddp ! Kd-decay species indexing
!    
            REAL(KIND = 8) :: cguess ! Initial guess for the concentration of the primary species, 
!                   in moles/kg H2O (molal) for species other than H2O and in kg for H2O.
            REAL(KIND = 8) :: tt     ! Total concentrations of the primary species (aqueous+solid+gas) 
            REAL(KIND = 8) :: gamp   ! Acitivity coefficient of primary species
            REAL(KIND = 8) :: cr          ! Rate in term of primary species due to mineral dissolution/precipitation
            REAL(KIND = 8) :: cr0         ! Previous iteration cr
!
            REAL(KIND = 8), DIMENSION(:), POINTER :: du2   ! Derivatives of total dissolved concentrations with 
!                       respect to concentrations of primary species
            REAL(KIND = 8), DIMENSION(:), POINTER :: dgamp ! Derivative of activity coefficient of primary species
            REAL(KIND = 8), DIMENSION(:), POINTER :: dr    ! Derivative of rate   (needs change, a same name appears in MULTI)    !!!!!!!!!!
!
            REAL(KIND = 8) :: crx        ! Generation rate of the primary species from kinetic reactions among primary species
            REAL(KIND = 8), DIMENSION(:), POINTER :: drx   ! Derivative of rate from kinetic reactions among primary species 
            REAL(KIND = 8), DIMENSION(:), POINTER :: ds2   ! Sum derivatives from other phases including ion-exchange,
!                       surface complexation, aqueous and mineral kinetics 
! 
         END TYPE Primary_Species
!
! 
         TYPE Secondary_Species
!
            CHARACTER(LEN = 20) :: naaqx  ! Name of secondary species
!
            INTEGER(KIND = 1) :: ncps     ! Number of primary species involved in the secondary species
!
            REAL(KIND = 8) :: cs         ! Concentrations of secondary species 
            REAL(KIND = 8) :: zs         ! The ion electric charge of secondary species
            REAL(KIND = 8) :: aks        ! Log10 K
            REAL(KIND = 8) :: gams       ! Acitivity coefficient of secondary species
!
            INTEGER(KIND = 1),DIMENSION(:), POINTER :: icps  ! Index of primary species in reaction stoichiometry
!
            REAL(KIND = 8), DIMENSION(:), POINTER :: stqs    ! Stoichiometric coefficient of primary species in secondary species 
            REAL(KIND = 8), DIMENSION(:), POINTER :: dgams   ! Derivative of activity coefficient of secondary species
            REAL(KIND = 8), DIMENSION(5) :: akcoes           ! Regression coefficients 
!
         END TYPE Secondary_Species
! 
! 
         TYPE PriSec_Species
!
            CHARACTER(LEN = 20) :: naaqt  ! Name of aqueous (primary+secondary) species
!
            REAL(KIND = 8) :: a0     ! Debye-Huckel a0 parameter (see Appendix H in the manual for details)
            REAL(KIND = 8) :: z      ! The ion electric charge of all aqueous species
            REAL(KIND = 8) :: zz2    ! z**2
            REAL(KIND = 8) :: ct     ! Concentrations of all aqueous species 
            REAL(KIND = 8) :: gamt   ! Acitivity coefficient of all aqueous species
! 
            REAL(KIND = 8), DIMENSION(:), POINTER :: dgamt ! Derivative of activity coefficient of all aqueous species
!
         END TYPE PriSec_Species
!
!
! ----------
! ...... Derived-type arrays     
! ----------
! 
         TYPE(Primary_Species),   ALLOCATABLE, DIMENSION(:) :: PSpe
!
         TYPE(Secondary_Species), ALLOCATABLE, DIMENSION(:) :: SSpe
!
         TYPE(PriSec_Species),    ALLOCATABLE, DIMENSION(:) :: ASpe
! 
! 
      END MODULE Aqueous_Species
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Kinetics_Among_Primary_Species  ! Including aqueous kinetics 
!                     and sorption kinetics and biodegradation                                  
!
! 
         SAVE
!
         INTEGER(KIND = 2) :: nrx    ! Counter of kinetic reactions       ! Pass through argument  !!
! 
! ----------
! ...... Derived-type variables     
! ---------- 
! 
         TYPE Primary_Kinetics
!
!
            INTEGER(KIND = 2)  :: ncp_rx     
!                  Number of primary species involved in the kinetic reaction
!
            INTEGER(KIND = 2)  :: i_mod     
!                  Rate law index: i_mod=1 for multiple monod kinetics
!
            INTEGER(KIND = 2)  :: n_mech    
!                  Number of mechanisms involved in the kinetic reaction
!

            INTEGER(KIND = 2), DIMENSION(5) :: ncp_rx1  ! (Maximum 5 mechiansms)   
!                  Number of species in the product term for each mechanisms
!
            INTEGER(KIND = 2), DIMENSION(5) :: ncp_rx2  ! (Maximum 5 mechiansms)   
!                  Number of species in the Monod term for each mechanisms
!
            INTEGER(KIND = 2), DIMENSION(5) :: ncp_rx3  ! (Maximum 5 mechiansms)   
!                  Number of species in the Inhibition term for each mechanisms
!
!
            INTEGER(KIND = 2), DIMENSION(5,5) :: ia1    ! Product, (5 mechanisms, 5 species)
!                  Indicator using activity (=1), concentration (=2), or total conc (=3)
!
            INTEGER(KIND = 2), DIMENSION(5,5) :: ia2    ! Monod, (5 mechanisms, 5 Monods)
!                  Indicator using activity (=1), concentration (=2), or total conc (=3)
!
            INTEGER(KIND = 2), DIMENSION(5,8) :: ia3    ! Inhib, (5 mechanisms, 8 Inhibitions)
!                  Indicator using activity (=1), concentration (=2), or total conc (=3)
!
!
            INTEGER(KIND = 2), DIMENSION(10) :: icprx      ! (Maximum 10 species) 
!                                  Species index for the kinetic reaction formulation
!
            INTEGER(KIND = 2), DIMENSION(5,5) :: icprx1     ! (5 mechanisms, 5 species)
!                                  Species index for the Product term of the rate law
!
            INTEGER(KIND = 2), DIMENSION(5,5) :: icprx2     ! (5 mechanisms, 5 Monod terms)
!                                  Species index for the Monod terms of the rate law
!
            INTEGER(KIND = 2), DIMENSION(5,8) :: icprx3     ! (5 mechanisms, 8 Inhibition terms)
!                                  Species index for the Inhibition terms of the rate law
!
!
            REAL(KIND = 8), DIMENSION(5) :: rkaq       !  (Maximum 5 mechanisms)
!                  Rate constant or maximum specific growth rate for each mechanism
!
            REAL(KIND = 8), DIMENSION(10) :: s_rx      !  (Maximum 10 species) 
!                  Stoichiometric coefficients of primary species in the kinetic reaction
!
            REAL(KIND = 8), DIMENSION(10) :: stqrx     !  (Maximum 10 species) 
!                  Same as the previous s_rx
!
            REAL(KIND = 8), DIMENSION(5,5) :: s_rx1    ! (5 mechanisms, 5 species)
!                  Stoichiometric coefficients of species in the Product term of the rate law
!
            REAL(KIND = 8), DIMENSION(5,5) :: s_rx2    ! (5 mechanisms, 5 Monod terms)
!                  Stoichiometric coefficients of species in the Monod term
!
            REAL(KIND = 8), DIMENSION(5,8) :: s_rx3    ! (5 mechanisms, 8 Inhibition terms)
!                  Stoichiometric coefficients of species in the Inhibition term
!
!
            REAL(KIND = 8) :: rkin2rx    ! Over all rate of the kinetic reaction 
!                                         (in terms of unit product speices)
!
!
            CHARACTER(LEN = 20), DIMENSION(10)  :: nam_rx      ! (Maximum 10 species) 
!                                  Name of primary species involved in the kinetic reaction
!
            CHARACTER(LEN = 20), DIMENSION(5,5) :: nam_rx1     ! (5 mechanisms, 5 species)
!                                  Name of species involved in the Product term of the rate law
!
            CHARACTER(LEN = 20), DIMENSION(5,5) :: nam_rx2     ! (5 mechanisms, 5 Monod terms)
!                                  Name of species involved in the Monod terms of the rate law
!
            CHARACTER(LEN = 20), DIMENSION(5,8) :: nam_rx3     ! (5 mechanisms, 8 Inhibition terms)
!                                  Name of species involved in the Inhibition terms of the rate law
!
!
         END TYPE Primary_Kinetics
!
!
! ----------
! ...... Derived-type array     
! ----------
! 
         TYPE(Primary_Kinetics),   ALLOCATABLE, DIMENSION(:) :: PKin
! 
! 
      END MODULE Kinetics_Among_Primary_Species
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Mineral_Phases
!
! 
         SAVE
!
         REAL(KIND = 8) :: sumvol2          ! Sum of mineral volume fractions (<=1)
!
! ----------
! ...... Derived-type variables for general mineral properties    
! ----------
! 
         TYPE Mineral_General
!
            CHARACTER(LEN = 20) :: namin    ! Name of mineral phase such as 'calcite'
!  
            INTEGER(KIND = 1) :: ncpm       ! Number of primary species involved in the mineral reaction
            INTEGER(KIND = 2) :: kineq      ! Pointer for minerals specified under equilibrium
            INTEGER(KIND = 1) :: ikin       ! Equilibrium or kinetic indicator; =1:kinetic, =0: equilibrium
            INTEGER(KIND = 1) :: idispre    ! A flag for the type of kinetic constraint: 
!                                               1 for dissolution only, 
!                                               2 for precipitation only 
!                                               3 for both (mineral can either precipitate or dissolve).  
!                                             Always set IDISPRE = 0 if IKIN = 0 and 
!                                                        IDISPRE > 0 if IKIN = 1.
            INTEGER(KIND = 1) :: iss        ! Solid solution index for mineral m
!
            REAL(KIND = 8) :: vmin       ! mole volume, cm**3/mol
            REAL(KIND = 8) :: dmolwm     ! molecular weight
            REAL(KIND = 8) :: cm         ! Mineral concentration, mol/kg H2O
            REAL(KIND = 8) :: pre2       ! Mineral volume fraction, ===> pre(:,:)
            REAL(KIND = 8) :: si2        ! Mineral saturation index (log10 Q/K)
            REAL(KIND = 8) :: si2_old    ! Mineral saturation index at the previous time step
            REAL(KIND = 8) :: akm        ! Log10 K
            REAL(KIND = 8) :: vol2       ! Mineral volume fraction in terms of solid
!
            REAL(KIND = 8) :: akin        ! Log10 K                                         !!!!!!!!!!
            REAL(KIND = 8) :: si2k        ! Working array for saturation index (log10 Q/K, or Q/K)
!
            REAL(KIND = 8) :: ssq0        ! Full supersaturation "window" in +log(K) units
            REAL(KIND = 8) :: sst1        ! Temperature (C) above which that window starts to
!                                      decrease exponentially with temperature (typically 25 C)
            REAL(KIND = 8) :: sst2        ! Temperature (C) at which ssq0 becomes one hundreth of
!                                      starting value (essentially no window anymore) 
            REAL(KIND = 8) :: ssq
!
            INTEGER(KIND = 1), DIMENSION(:), POINTER :: icpm  ! Index of primary species in reaction stoichiometry
!
            REAL(KIND = 8), DIMENSION(5)          :: akcoem   ! Regression coefficients 
            REAL(KIND = 8), DIMENSION(:), POINTER :: stqm     ! Stoichiometric coefficient of primary species in the mineral
! 
         END TYPE Mineral_General
!
!
! ----------
! ...... Solid solutions     
! ----------
! 
         INTEGER(KIND = 1), ALLOCATABLE, DIMENSION(:)   :: ncpss   ! Mumber of endmembers in solid solution index n (n=iss(m))
         INTEGER(KIND = 1), ALLOCATABLE, DIMENSION(:,:) :: icpss   ! Mineral index of endmember k, in solid solution index n
!
!
! ----------
! ...... Derived-type variables for mineral kinetic parameters    
! ----------
! 
         TYPE Min_Kine_Para
!
            INTEGER(KIND = 1) :: idep      ! A flag for dissolution rate constant dependence on pH or multiple mechanisms.  
!                        If idep = 0: pH dependent rate constants and multiple mechanisms 
!                           are not considered. 
!                        If idep = 1: the rate dependence on pH. 
!                        If idep = 2: need to include data for multiple mechanisms
!
            INTEGER(KIND = 1) :: imflg2     ! A flag for surface area conversion:
!                                               IMFLG = 0 for cm2/g mineral
!                                               IMFLG = 1 for m2 rock area/m3 medium
!                                               IMFLG = 2 for m2/m3 mineral
!
            REAL(KIND = 8) :: amin2      ! Mineral surface area, unit depends on the imflg2 above
            REAL(KIND = 8) :: rad2       ! Radius of mineral grain (in m) used to calculate surface area for 
!                         initial formation of secondary phase. If RAD = 0.0, the initial 
!                         surface area is calculated from RNUCL 
            REAL(KIND = 8) :: rkf        ! Rate constant at 25°C, mol/m2/sec 
            REAL(KIND = 8) :: ea         ! The activation energy, kJ/mol
            REAL(KIND = 8) :: ck1        ! The exponents respectively in rate law (Eq.B.5 in the manual).
            REAL(KIND = 8) :: ck2        ! Similar to ck1
            REAL(KIND = 8) :: acfdiss
            REAL(KIND = 8) :: bcfdiss
            REAL(KIND = 8) :: ccfdiss
!       acfdiss, acfdiss, and acfdiss above should be set to zero, unless a different form of rate 
!           constant dependence with temperature is desired.  
!       This alternate form is: log(k) = a + b·T + c/T, where T is absolute temperature in K 
!       and log is in base 10.  To enable this option, rkf must be set to 1.0, EA must be set 
!       to 0.0, ck1 and ck2 can be set to any value, and acfdiss, acfdiss, and acfdiss must 
!       be specified as the coefficients a, b, and c, respectively, in the above expression.
! 
!      
            INTEGER(KIND = 1) :: ideprec     ! A flag for precipitation rate constant, similar meaning to idep above
            INTEGER(KIND = 1) :: ndis        ! Number of additional mechanisms for dissolution 
            INTEGER(KIND = 1) :: npre        ! Number of additional mechanisms for precipitation
            INTEGER(KIND = 1) :: nplaw       ! Precipitation law index. nplaw = 0 for Eqs. B.5 and B.12) 
!                         in Appendix B of the manual; NPLAW = 1 for Eq. B.8.  
!
            REAL(KIND = 8) :: rkprec         ! Similar to rkf above for dissolution
            REAL(KIND = 8) :: eaprec         !
            REAL(KIND = 8) :: ck1prec        !
            REAL(KIND = 8) :: ck2prec        !
            REAL(KIND = 8) :: acfprec        !
            REAL(KIND = 8) :: bcfprec        !
            REAL(KIND = 8) :: ccfprec        !
            REAL(KIND = 8) :: rnucl       ! The initial volume fraction (Vmineral/Vsolid) to be assumed for 
!                         calculating initial effective surface area if the mineral is not 
!                         present at the start of a simulation but precipitates as a new 
!                         reaction product.  If zero, RNUCL is assumed to be 10e-5. 

            REAL(KIND = 8) :: ssqk0       ! Log (Q/K) gap (supersaturation window, see Eq, B.13 in Appendix B). 
!                         A zero value represents no gap. 
            REAL(KIND = 8) :: sstk1       ! Temperature (in °C) at which to begin reducing gap
            REAL(KIND = 8) :: sstk2       ! Temperature (in °C) endpoint at which the gap has diminished to 
!                         nearly zero (1% of original value).  The gap decreases exponentially 
!                         from the first (SSTK1) to the second (SSTK2) temperature, and SSTK2 
!                         must always be greater than SSTK1.
            REAL(KIND = 8) :: ssqk        ! A working array calculated from above three variables    
!      
            REAL(KIND = 8) :: aH1         ! Activity of H+ near which pH influence starts
            REAL(KIND = 8) :: aH2         ! Activity of H+ near which pOH influence starts
            REAL(KIND = 8) :: aH1p        ! Similar to aH1 above but here for precipitation
            REAL(KIND = 8) :: aH2p        !
            REAL(KIND = 8) :: aHexp       ! aHexp is aH+ exponent (slope of log(rate) with pH)
            REAL(KIND = 8) :: aHexpp      !
            REAL(KIND = 8) :: aOHexp      ! aOH exponent (slope of log(rate) with pOH)
            REAL(KIND = 8) :: aOHexpp     !
!
            REAL(KIND = 8) :: rkf_ds      ! Dissolution rate constants from additional mechanisms 
            REAL(KIND = 8) :: rkprec_ds   ! Precipitation rate constants from additional mechanisms 
            REAL(KIND = 8) :: rkin2       ! Dissolution/precipitation rate
!
            INTEGER(KIND = 1), DIMENSION(5) :: nspds    ! Number of speciess involved in one mechanism
            INTEGER(KIND = 1), DIMENSION(5) :: nsppr    ! Similar to nspds but for precipitation
            INTEGER(KIND = 1), DIMENSION(5,5) :: ids    ! Dependent species pointer in ct for dissoltion
            INTEGER(KIND = 1), DIMENSION(5,5) :: idsp   ! Dependent species pointer in ct for precipitation
!
            REAL(KIND = 8), DIMENSION(5)   :: rkds      ! Dissolution rate constant for the mechanism
            REAL(KIND = 8), DIMENSION(5)   :: eads      ! Activation energy for the mechanism
            REAL(KIND = 8), DIMENSION(5)   :: rkdsp     ! Similar to rkds but for precipitation
            REAL(KIND = 8), DIMENSION(5)   :: eadsp     ! The activation energy, kJ/mol
            REAL(KIND = 8), DIMENSION(5,5) :: expds     ! Exponatial term for the mechanism
            REAL(KIND = 8), DIMENSION(5,5) :: expdsp    ! Similar to expds but for precipitation
!
            CHARACTER(LEN = 20), DIMENSION(5,5) :: nadis  ! Name of species involved in dissolution rate constant 
            CHARACTER(LEN = 20), DIMENSION(5,5) :: napre  ! Name of species involved in precipitation rate constant 
!
!
         END TYPE Min_Kine_Para
!
! 
         INTEGER(KIND = 1) :: ndep    ! Number of minerals with kinetic rates from mechanisms other than neutral
!
         INTEGER(KIND = 1) :: nsalt   ! Number of precipitates at dry grid blocks
!
         INTEGER(KIND = 2), DIMENSION(:), POINTER :: isalt   ! A flag to indicate that the mineral may precipitate in a dry grid block 
!                       as a result of complete evaporation (when liquid saturation < sl1min 
!                       specified in the solute.inp file), or if there is water flux into 
!                       the grid block that dries out during the flow step (and therefore 
!                       liquid saturation is zero). The mineral with M1 = 1 precipitates first, 
!                       with isalt( ) = 2 second, and so on. If this flag is set to zero, 
!                       then the mineral will not be formed in the dry grid block.
!
! ----------
! ...... Derived-type arrays     
! ----------
! 
         TYPE(Mineral_General),   ALLOCATABLE, DIMENSION(:) :: MGen
!
         TYPE(Min_Kine_Para),     ALLOCATABLE, DIMENSION(:) :: MKin
!
! 
      END MODULE Mineral_Phases
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Gaseous_Species
!
! 
         SAVE
!
! ----------
! ...... Derived-type variables for gaseous species    
! ----------
! 
         TYPE Gas_Species
!
            CHARACTER(LEN = 20) :: nagas    ! Name of gaseous species such as 'co2(g)'
!
            INTEGER(KIND = 1) :: ncpg       ! Number of primary species involved in the gas reaction
!
            REAL(KIND = 8) :: dmwgas     ! Molecular weight, g/mol
            REAL(KIND = 8) :: diamol     ! Molecular diameter, m
            REAL(KIND = 8) :: akg        ! Log10 K
            REAL(KIND = 8) :: cmg        ! Amount of gaseous species, mol/kg H2O
            REAL(KIND = 8) :: cg         ! Gas partial pressure, bar
            REAL(KIND = 8) :: gamg       ! Gas fugacity coefficient
            REAL(KIND = 8) :: sig2       ! Saturation index
!
            INTEGER(KIND = 1), DIMENSION(:), POINTER :: icpg  ! Index of primary species in reaction stoichiometry
!
            REAL(KIND = 8), DIMENSION(5)          :: akcoeg   ! Regression coefficients 
            REAL(KIND = 8), DIMENSION(:), POINTER :: stqg     ! Stoichiometric coefficient of primary species in the gas
!
!
         END TYPE Gas_Species
!
!
! ----------
! ...... Derived-type arrays for gaseous species     
! ----------
! 
         TYPE(Gas_Species),   ALLOCATABLE, DIMENSION(:) :: GSpe
!
! 
      END MODULE Gaseous_Species
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Ion_Exchange
!
! 
         SAVE
!
!
         INTEGER(KIND = 1) :: iex       ! Convention index: 1 = Gaines-Thomas; 
!                                                           2 = Vanselow; 
!                                                           3 = Gapon 
!                                         The value of iex must be the same for 
!                                              all the exchanged cations.
!
         REAL(KIND = 8) :: cec2  ! Cation exchange capacity (meq/100 g of solid).
!
!
! ----------
! ...... Derived-type variables for exchanged species    
! ----------
! 
         TYPE Exchanged_Species
!
            CHARACTER(LEN = 20) :: naexc  ! Name of exchanged species such as "Na+'
!      
            INTEGER(KIND = 1) :: nbx ! Exchange species index to primary species
!
            REAL(KIND = 8) :: ekx    ! Slectivity of cations
            REAL(KIND = 8) :: cx     ! Exchanged species concentrations
!
            REAL(KIND = 8), DIMENSION(:), POINTER :: stqx  ! Stoichiometric coefficient of primary species in the exchange reaction
            REAL(KIND = 8), DIMENSION(:), POINTER :: dcx   ! Derivates of exchanged species concentrations
!
!
         END TYPE Exchanged_Species
!
!
! ----------
! ...... Derived-type arrays for exchanged species     
! ----------
!
         TYPE(Exchanged_Species),   ALLOCATABLE, DIMENSION(:) :: XSpe
!
! 
! ----------
! ...... Arrays for multi-site exchange   
! ----------
! 
!
         REAL(KIND = 8), DIMENSION(:), POINTER ::      cecM   ! Batch CEC for multi-sites
!
         REAL(KIND = 8), DIMENSION(:,:), POINTER ::    ekxM, &! Selectivity
                                                       cxM    ! Exchanged species concentrations
!
         REAL(KIND = 8), DIMENSION(:,:,:), POINTER ::  dcxM   ! Devivatives

!
! 
      END MODULE Ion_Exchange
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Surface_Complexes
!
! 
         SAVE
!
!
         INTEGER(KIND = 1) :: iadmod   ! Indicator of adsorption model (0, 1, 2), 
!                                          =0: not considering potential term
!                                          =1: constant capatance model
!                                          =2: double layer model
!
         REAL(KIND = 8) :: tads2       ! Total surface sites (concentrations)
         REAL(KIND = 8) :: supadn2     !
         REAL(KIND = 8) :: adfactor    !
         REAL(KIND = 8) :: phip2       ! Potential term
         REAL(KIND = 8) :: phip2new    !
!
!
! ----------
! ...... Derived-type variables for exchanged species    
! ----------
! 
         TYPE Surface_Species
!
            CHARACTER(LEN = 20) :: naads  ! Name of surface species
!
            INTEGER(KIND = 1) :: ncpad    ! Number of primary species involved in the surface reactions
            INTEGER(KIND = 1) :: n_nd     ! Index of master ad. species in primary list
!
            REAL(KIND = 8) :: zd          ! Charge of surface species
            REAL(KIND = 8) :: akd         ! Log10 K
            REAL(KIND = 8) :: cd          ! Surface species comcentration, mol/kg H2O
!
            INTEGER(KIND = 1), DIMENSION(:), POINTER :: icpad ! Index of primary species in reaction stoichiometry
!
            REAL(KIND = 8), DIMENSION(5)          :: akcoead  ! Regression coefficients 
            REAL(KIND = 8), DIMENSION(:), POINTER :: stqd  ! Stoichiometric coefficient of primary species in the surface reaction
            REAL(KIND = 8), DIMENSION(:), POINTER :: dcd   ! Derivative of the concentrationb
!
!
         END TYPE Surface_Species
!
!
! ----------
! ...... Derived-type arrays for surface species     
! ----------
! 
         TYPE(Surface_Species),   ALLOCATABLE, DIMENSION(:) :: DSpe
!
! 
      END MODULE Surface_Complexes
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Kd_Adsorption_Decay
!
! 
         SAVE
!
!
! ----------
! ...... Derived-type variables for species with linear Kd adsorption or/and decay    
! ----------
! 
         TYPE KdDecay_Species

            CHARACTER(LEN = 20) :: nakdd     ! Name of species with Kd or/and decay 
!
            REAL(KIND = 8)      :: decayc,  &! Decay constants, 1/s
!               If decayc<0 calculate the following two parameters: ln(decayc) = a - b/T
        	                   a_TDecay,&! Thermal decay parameter, a
                                   b_TDecay  ! Thermal decay parameter, b
            REAL(KIND = 8)      :: sden2     ! Solid density, kg/dm3        
            REAL(KIND = 8)      :: vkd2      ! Kd value(in (l/kg which is mass/kg 
!                                                solid divided by mass/l solution). 
!                                              If SDEN2 = 0.0, VKD2 automatically represents 
!                                                retardation factor (>= 1).    
!
         END TYPE KdDecay_Species
!
!
! ----------
! ...... Derived-type arrays for with linear Kd adsorption and decay    
! ----------
! 
         TYPE(KdDecay_Species),   ALLOCATABLE, DIMENSION(:) :: KSpe
!
! 
      END MODULE Kd_Adsorption_Decay
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Batch_Chem_Jacobian
!
! 
         SAVE
!
!
         INTEGER(KIND = 8) :: nmat   ! Dimension of the Jacobian matrix for solving the geochemical system

         INTEGER(KIND = 8) :: nsat   ! Number of minerals in the Jacobian. Minerals get kicked off 
!                                if exhausted (if abundance=0 and log(q/k) < 0)
!
         INTEGER(KIND = 8), DIMENSION(:), POINTER :: isat  ! Index of mineral reminded in the Jacobian
         INTEGER(KIND = 8), DIMENSION(:), POINTER :: mout  ! Indicator
!
         INTEGER(KIND = 8) :: nsatg  ! Similar to nsat above but here for gas
!
         INTEGER(KIND = 8), DIMENSION(:), POINTER :: isatg   ! Index of gas reminded in the Jacobian
         INTEGER(KIND = 8), DIMENSION(:), POINTER :: moutg   ! Indicator
	    
         REAL(KIND = 8), DIMENSION(:,:), POINTER :: amat   ! Jacobian Matrix
!
         REAL(KIND = 8), DIMENSION(:), POINTER :: bmat     ! Right-hand side residual term
!
         INTEGER(KIND = 8), DIMENSION(:), POINTER :: indx  ! Working array
!
         REAL(KIND = 8), DIMENSION(:), POINTER :: fp       ! Convergence parameters  
!
! 
      END MODULE Batch_Chem_Jacobian
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!  All arraies in the following module will be dellocated 
!         after assigning to each grid block except sden and vkd.
!
      MODULE IniChemZoneProperties
!
! 
         SAVE
!
!
         INTEGER(KIND = 1), DIMENSION(:,:), POINTER :: icon_iniZ
!
         REAL(KIND = 8), DIMENSION(:,:), POINTER :: c_iniZ    ! Species concentration, mol/kg 
         REAL(KIND = 8), DIMENSION(:,:), POINTER :: U_iniZ    ! Total dissolved concentrations of initial waters
         REAL(KIND = 8), DIMENSION(:),   POINTER :: Tc_iniZ   ! Temperature   
         REAL(KIND = 8), DIMENSION(:),   POINTER :: pH_iniZ   ! pH
         REAL(KIND = 8), DIMENSION(:),   POINTER :: aw_iniZ   ! Water activity 
         REAL(KIND = 8), DIMENSION(:,:), POINTER :: pfug_iniZ ! Gas partial pressure, bar    
!
         REAL(KIND = 8), DIMENSION(:,:), POINTER :: Uboud     ! Total dissolved con. of boundary waters
!
!
         INTEGER(KIND = 1), DIMENSION(:,:), POINTER :: imflg_iniZ   ! Flag for surface area unit
!
         REAL(KIND = 8), DIMENSION(:,:), POINTER :: vol_iniZ   ! Mineral volume fraction  
         REAL(KIND = 8), DIMENSION(:,:), POINTER :: rad_iniZ   ! Grain radius  
         REAL(KIND = 8), DIMENSION(:,:), POINTER :: amin_iniZ  ! Mineral surface area   
         REAL(KIND = 8), DIMENSION(:), POINTER :: sumvol_iniZ  ! Sum of all mineral volumes   
!
!
         REAL(KIND = 8), DIMENSION(:,:), POINTER :: cg_iniZ    ! Gas partial pressure     
!
!
         INTEGER(KIND = 4), DIMENSION(:), POINTER :: ipptyp_iniZ  ! Permeability-porosity laws, 
!                                      1: Carman-Kozeny, 
!                                      3: cubic law, 
!                                      4: modified Cubic Law, 
!                                      5: Verma-Pruess 
!
         REAL(KIND = 8), DIMENSION(:), POINTER :: apppar_iniZ   ! Fracture aperture for modified Cubic Law    
         REAL(KIND = 8), DIMENSION(:), POINTER :: bpppar_iniZ   ! Fracture spacing  for modified Cubic Law  
!
!
         REAL(KIND = 8), DIMENSION(:), POINTER :: sadsdum_iniZ     
         REAL(KIND = 8), DIMENSION(:), POINTER :: tads_iniZ   ! Total surface sites  
!
         REAL(KIND = 8), DIMENSION(:,:), POINTER :: sden  ! Solid density, kg/dm**3
         REAL(KIND = 8), DIMENSION(:,:), POINTER :: vkd   ! Kd value (in (l/kg which is mass/kg solid divided by mass/l solution). 
!                               If SDEN =0.0, VKD automatically represents retardation factor (³ 1).  
!

         REAL(KIND = 8), DIMENSION(:,:), POINTER :: cec_iniZ   ! Cation exchange capacity for multi-sites    

!
!
      END MODULE IniChemZoneProperties
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Para_DH_ActivityCoefficients
! 
      SAVE
!
!
         REAL(KIND=8) :: adh, bdh, bi, bil
!
!
      END MODULE Para_DH_ActivityCoefficients
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!.....Variables and Arrays in the Pitzer Model
!
!
      MODULE IonPairs_Index
!
! 
         SAVE
!
!
         INTEGER (KIND = 4) :: &
             mpair,  &  ! Maximum number of ionic pairs and triplets
             npair,  &  ! Number of cation-anoin ionic pairs
             npair1, &  ! Number of cation-anoin ionic pairs plus number of neutral-cation pairs
             npair2, &  ! Number of cation-anoin ionic pairs plus number of neutral-cation pairs plus number of neutral-anion pairs
             npos,   &  ! Number of cations in the system
             nneg,   &  ! Number of anions in the system
             nneu,   &  ! Number of neutral species in the aqueous system
             ncc,    &  ! Number of cation-cation pairs.
             naa,    &  ! Number of anion-anion pairs.
             ncca,   &  ! Number of cation1-cation2-anion triplets.
             ncaa,   &  ! Number of cation-anion1-anion2 triplets.
             nnca       ! Number of neutral-cation-cation triplets.
!
         INTEGER (KIND = 4), DIMENSION(:), POINTER :: &
             ippos(:), &! Reserved for later development, not currently used.
             ipneg(:), &! Reserved for later development, not currently used.
             ipnue(:)   ! Reserved for later development, not currently used.
!
         INTEGER (KIND = 4), DIMENSION(:), POINTER :: &
             itpos(:), &! Index of cations in the entire chemical system
             itneg(:), &! Index of anions in the entire chemical system
             itneu(:)   ! Index of neutral species in the entire chemical system
!
         INTEGER (KIND = 4), DIMENSION(:,:), POINTER :: &
             ipair,    &! Reserved for later development, not currently used.
             icc,      &! Index of the cations in cation-cation pair icc().
             iaa,      &! Index of the anions in anion-anion pair iaa().
             icca,     &! Index of the ions in cation1-cation2-anion triplets icca().
             icaa,     &! Index of the ions in cation-anion1-anion2 triplets icaa().
             inca,     &! Index of the neutral-cation-cation triplets inca().
             ipp        ! Index of pairs, i, j are the indexes of ions or neutral species in the ionic pair ipp(i,j).
!
!
      END MODULE IonPairs_Index
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!.....Pitzer Ionic Interaction Parameters
!
!
      MODULE IonInteract_Parameters
!
! 
         SAVE
!
!
         REAL (KIND = 8) ::  &
             aphi,     &  ! a-phi parameter
             segmaz,   &  ! Intermediate variable
             otherstr, &  ! Intermediate variable
             otherstr2    ! Intermediate variable
!
         REAL (KIND = 8), DIMENSION(:), POINTER :: &
             beta0,      &  ! beta0
             beta1,      &  ! beta1
             beta2,      &  ! beta2
             cmx,        &  ! C
             ramdapos,   &  ! Lambda of neutral and cation
             ramdaneg,   &  ! Lambda of neutral and anion
             theta_cc,   &  ! Mixing term parameters 
             theta_aa,   &  ! Mixing term parameters
             psi_cca,    &  ! Triplets coefficient
             psi_caa,    &  ! Triplets coefficient
             zeta_nca,   &  ! Triplets coefficient
             alfamx1,    &  ! alpha1
             alfamx2,    &  ! alpha2
             dbmx_dstr,  &  ! Intermediate array
             ddbmx_dstr, &  ! Intermediate array
             dbphimx_dstr   ! Intermediate array
!
         REAL (KIND = 8), DIMENSION(:,:), POINTER :: &     
             betat0,     &  ! Temperature dependent coefficients of beta0
             betat1,     &  ! Temperature dependent coefficients of beta1
             betat2,     &  ! Temperature dependent coefficients of beta2
             cmxt,       &  ! Temperature dependent coefficients of C
             ramdapost,  &  ! Temperature dependent coefficients of lambda 
             ramdanegt,  &  ! Temperature dependent coefficients of lambda
             thetat_cc,  &  ! Temperature dependent coefficients of mixing term
             thetat_aa,  &  ! Temperature dependent coefficients of mixing term
             psit_cca,   &  ! Temperature dependent coefficient of triplets
             psit_caa,   &  ! Temperature dependent coefficient of triplets
             zetat_nca      ! Temperature dependent coefficient of triplets
!
!
      END MODULE IonInteract_Parameters
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!.....Derivative terms and intermediate derivative terms for the Pitzer
!
!
      MODULE DerivativeTerms_Pitzer
!
! 
         SAVE
!
!
         REAL (KIND = 8), DIMENSION(:), POINTER :: &     
             bphimx,    &   ! B phi term
             bmx            ! B term
!
         REAL (KIND = 8), DIMENSION(:,:), POINTER :: &    
             dbphimx,   &   ! B phi Derivative term
             dbmx,      &   ! B Derivative term
             ddbmx,     &   ! Second order Derivative term of B phi.
             dgamt,     &   ! Derivative term of activity coefficients of total aqueous species
             dgamp,     &   ! Derivative term of activity coefficients of primary species
             dgams          ! Derivative term of activity coefficients of secondary species
!
!
      END MODULE DerivativeTerms_Pitzer
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE CutOff_Threshold
!
! 
         SAVE
! 
! ----------
! ...... Double precision variables     
! ----------
! 
!........Cutoff for skipping geochemical calculations 
!
         REAL(KIND = 8) :: sl1min      ! Minimum liquid saturation
         REAL(KIND = 8) :: d1min       ! Minimum inter-block distance
         REAL(KIND = 8) :: stimax      ! Maximum stoichiometric ionic strength 
         REAL(KIND = 8) :: phimin      ! Porosity cutoff (Default: 1.0D-10), 
!                   if porosity <= phimin, only mineral dissolution is allowed 
         REAL(KIND = 8) :: cnfact      ! Factor for reaction source terms 
!                     =1 fully implicit, =0 fully explicit (for kinetic minerals only)
!
         REAL(KIND = 8) :: str_threshold  ! Ionic activity models switch point 
!                                             between DH and Pitzer models
!
      END MODULE CutOff_Threshold
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE ChemOption_Variables
!
!  
         SAVE
! 
! ----------
! ...... Integer variables     
! ----------
! 
         INTEGER(KIND=1)  :: ispia     ! Iteration scheme between transport and reaction; 
!                                          0: sequential iteration 
!                                          2: sequential no iteration
!
         INTEGER(KIND=1)  :: inibound  ! Identifying boundary water solution, 0: No, 1: Yes 
!
         INTEGER(KIND=1)  :: isolvc    ! Flag for the linear equation solver for transport
!                                          2 - DSLUBC, a bi-conjugate gradient solver
!                                          3 - DSLUCS, a Lanczos-type bi-conjugate gradient solver
!                                          4 - DSLUGM, a general minimum residual solver
!                                          5 - DLUSTB, a stabilized bi-conjugate gradient solver            
!
         INTEGER(KIND=1)  :: ngas1     ! Inclusion of gaseous chemical species transport 
!                                          0 - not included
!                                          1 - included
!                   If gas partial pressure remains constant with time, set NGAS1=0. 
!                   For using EOS2 and ECO2 modules, always set NGAS1=0 because CO2 pressure is 
!                       handled in the flow simulation.
!
         INTEGER(KIND=1)  :: ichdump   ! Flag to enable printout of chemical speciation results 
!                       0 - disabled
!                       1 - printout of chemical speciation at each grid block and each time step
!                       2 - printout of chemical speciation at times specified by NWTI in the 
!                           Record_7 and grid blocks specified in Record_8 of solute.inp. 
!                   If this option is enabled, the program will abort after the output of 
!                       speciation results for the first 1000 grid blocks and/or time steps, 
!                       to avoid accidentally filling up disk space.
!
         INTEGER(KIND=1)  :: kcpl      ! Flag to consider feedback effects of changes of porosity, permeability, 
!                       and capillary pressure due to mineral dissolution and precipitation on fluid flow.
!                       0 - disabled 
!                       1 - enabled
!                       2 - only monitor the changes (printout), but without feedback on fluid flow.
!
         INTEGER(KIND=1)  :: Ico2h2o   ! Flag to consider effects of CO2 and H2O reaction source/sink terms on fluid flow 
!                   calculations. ICO2H2O is only used for the EOS2 and ECO2 flow modules. 
!                   For other flow modules, set ICO2H2O = 0.
!                       0 - effects ignored 
!                       1 - only effects of CO2 reaction source/sink terms
!                       2 - effects of both CO2 and H2O reaction source/sink terms
!
         INTEGER(KIND=1)  :: numdr     ! Flag for calculation of derivatives of mineral kinetic rates with respect to 
!                                       concentrations of primary species.
!                       0 - Analytical method (normally used) 
!                       1 - Numerical method
! 
! ----------
! ...... Double precision variables     
! ----------
! 
         REAL(KIND = 8) :: rcour       ! Both a variable and a flag to limit the time step size
!                       =0.0 - No Courant number limitation on transport
!                       >0.0 - Transport time steps for aq. species and gas are limited 
!                              by Courant number    
!                       <0.0 - Transport time steps for aq. species limited by Courant number
!

         REAL(KIND = 8) :: wtime       ! Time weighting factor, ranging from 0.0 to 1.0.  
!                                          WTIME = 1.0 (implicit) is suggested.
         REAL(KIND = 8) :: wupc        ! Upstream weighting factor, ranging from 0.0 to 1.0.  
!                                          WUPC = 1.0 (fully upstream) is suggested.
!
      END MODULE ChemOption_Variables
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Convergence_Variables
!
!  
         SAVE
! 
! ----------
! ...... Integer variables     
! ----------
! 
         INTEGER(KIND=8)  :: maxitpfl
         INTEGER(KIND=8)  :: maxitptr    ! Maximum allowed iterations for reactive transport equations
         INTEGER(KIND=8)  :: maxitpch    ! Maximum allowed iterations for solving geochemical system
         INTEGER(KIND=8)  :: maxitpad    ! maximum number of  iterations allowed for solving sorption via surface complexation. 
         INTEGER(KIND=8)  :: iterfl      ! Iterations to solve flow
         INTEGER(KIND=8)  :: itermod     !
         INTEGER(KIND=8)  :: itertr      ! Iteration number for solving reactive transport equations
         INTEGER(KIND=8)  :: iterch      ! Iteration number for solving batch geochemistry
         INTEGER(KIND=8)  :: maxitch     ! Maximum iterations to solve chemistry
         INTEGER(KIND=8)  :: iterad      ! Iteration number for solving adsorption
         INTEGER(KIND=8)  :: maxitad     ! Maximum iterations to solve adsorption
         INTEGER(KIND=8)  :: ireturn     ! = 1, indicate to return from subroutine due to convergence
!
         INTEGER(KIND=8)  :: max_chem_it
! 
! ----------
! ...... Double precision variables     
! ----------
! 
         REAL(KIND = 8) :: tolfl       ! Maximum tolerance for flow 
         REAL(KIND = 8) :: toltr       ! Maximum tolerance for reactive transport equations
         REAL(KIND = 8) :: tolch       ! Maximum tolerance for solving geochemical system
         REAL(KIND = 8) :: tolad       ! Maximum tolerance for solving adsorption
         REAL(KIND = 8) :: toldc       ! Concentration tolerance for QSS (chemical steady-state)
         REAL(KIND = 8) :: toldr       ! Mineral kinetic rate tolerance for QSS (chemical steady-state
         REAL(KIND = 8) :: averitch    ! Average iterations to solve chemistry
         REAL(KIND = 8) :: countch     ! Counts number of times going through iteration process   
         REAL(KIND = 8) :: averitad    ! Average iterations to solve adsorption
         REAL(KIND = 8) :: countad     ! Iteration number counts for adsorption

         REAL(KIND = 8) :: delt_ch
! 
! ----------
! ...... Character variables     
! ----------
! 
        CHARACTER (LEN = 16) :: delt_conne
        CHARACTER (LEN =  5) :: id_chem
!
!
      END MODULE Convergence_Variables
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE WriteChemControl
!
!  
         SAVE
! 
! ----------
! ...... Integer variables     
! ----------
! 
         INTEGER(KIND=8)  :: nwti        ! Writing frequency in time
         INTEGER(KIND=8)  :: nwnod       ! Number of grid blocks for writing in time
         INTEGER(KIND=8)  :: nwcom       ! Number of components for writing
         INTEGER(KIND=8)  :: nwmin       ! Number of minerals for writing
         INTEGER(KIND=8)  :: iwcomt      ! =1: write total component concentrations; 
!                                          =0: write species concentrations
         INTEGER(KIND=8)  :: iconflag    ! Flag for concentration unit 
!                                          =1:mol/l; 
!                                          =0: mol/kg H2O)
         INTEGER(KIND=8)  :: minflag     ! Flag for mineral unit (=0: change in mol/m3; 
!                                          =1; change in volume fraction;  
!                                          =2; current volume fraction.
!
! ----------
! ...... Integer arrays     
! ----------
! 
         INTEGER(KIND = 8), DIMENSION(:), POINTER :: iwnod   ! Grid block numbers for writing
         INTEGER(KIND = 8), DIMENSION(:), POINTER :: iwcom   ! Species numbers for writing
         INTEGER(KIND = 8), DIMENSION(:), POINTER :: iwmin   ! Mineral numbers for writing
!
! ----------
! ...... Character array     
! ----------
! 
         CHARACTER(LEN =  5), DIMENSION(:), POINTER  :: elemw ! Name of grid blocks for time evolution

!
      END MODULE WriteChemControl
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Grid_HydroTherm
!
!  
         SAVE
! 
! ----------
! .......Real arrays for hydrological and thermal properties     
! ----------
!
         REAL(KIND = 8), DIMENSION(:), POINTER :: &
              SLOLD,      &  ! Liquid saturation at the previou time step
              SGOLD,      &  ! Gas saturation at the previou time step
              SL1,        &  ! Liquid saturation at the previous time step
              SG1,        &  ! Gas saturation at the current time step
              PHI0,       &  ! Porosity at t = 0
              PHIOLD,     &  ! Porosity at the previous time step
              tc,         &  ! Temperature at the current time step, oC
              tcold,      &  ! Temperature at the previous time step, oC
              Tc_EOS9,    &  ! Initial temperature using EOS9, oC
              densw,      &  ! Density of water, kg/m^3 
              fugcoeCO2,  &  ! CO2 fugacity coefficient from ECO2 and ECO2N
              Rh             ! Relative humidity obtained from EOS
!
!
      END MODULE Grid_HydroTherm
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Grid_AqueousPhase
!
!  
         SAVE
! 
! ----------
! .......Real 1-D arrays for chemical conditions     
! ----------
!
         REAL(KIND = 8), DIMENSION(:), POINTER :: &
              ph,        &      ! pH values
              str_node,  &      ! Ionic strength
              aw                ! Water activity
! 
! ----------
! .......Real 2-D arrays for aqueous Concentrations     
! ----------
!
         REAL(KIND = 8), DIMENSION(:,:), POINTER :: & 
              c,    &   ! Concentrations of all aqueous species, mol/kg H2O
              ut,   &   ! Total dissolved concentrations of component 
              utem, &   ! Tempo. store ut 
              utold     ! Total dissolved concentrations of component 
!                           at the previous time step, mol/l
!                           at the current time step, mol/l
!
!
      END MODULE Grid_AqueousPhase
!
!  
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Grid_MineralPhases
!
!  
         SAVE
! 
! ----------
! .......Real 2-D arrays for general mineral phases     
! ----------
!
         REAL(KIND = 8), DIMENSION(:,:), POINTER :: &
              pinit, &  ! Mineral abundance at t = 0, mol/dm3 medium
              pre0,  &  ! Mineral abundance at the previous time step, mol/dm**3 medium
              pre,   &  ! Mineral abundance at the current  time step, mol/dm**3 medium
              amin,  &  ! Reactive surface area, m**2/m**3 mineral
              rad,   &  ! Radius of mineral grain, m
              SIM,   &  ! Mineral saturation index (log10(Q/K))
              RKIN0, &  ! Mineral dissolution rate (negative values for precipitation) 
!                          at the previous time step
              rkin      ! Mineral dissolution rate at the current time step
! 
! ----------
! .......Real 1-D arrays for reaction source/sink term feed back to the flow simulation     
!........For modules EOS2, or ECO2  
! ----------
!
         REAL(KIND = 8), DIMENSION(:), POINTER :: &
           Rh2o,  &     ! H2O generated from reactions
           Rco2,  &     ! CO2 generated from reactions
           SMco2        ! CO2 trapped in solid phase
!
!
      END MODULE Grid_MineralPhases
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Grid_ChemSteadyState
!
!  
         SAVE
! 
! ----------
! ...... Integer variables     
! ----------
! 
         INTEGER(KIND = 4) :: IFLOWSS, &   ! Indicator for flow steady-state
                              JSTEADY, &   ! Indicator for chemical steady-state
                              NBLOCK,  &   ! Block where mineral exhausted 
                              NMINERAL    ! Exhausted Mineral numbers
! 
! ----------
! ...... Double precision variables     
! ----------
! 
         REAL(KIND = 8)    :: TIMESTEA    ! Time of chemical steady-state, s

! 
! ----------
! .......Real 2-D arrays for general mineral phases     
! ----------
!
         REAL(KIND = 8), DIMENSION(:,:), POINTER :: &
              DELTAP, & ! Mineral abundance increment for chemical steady-state
              DRATE     ! Rate for chemical qusi-steady state
!
!
      END MODULE Grid_ChemSteadyState
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Grid_ExchangeAdsorption
!
!  
         SAVE
! 
! ----------
! .......Real 2-D array for Cation exchange     
! ----------
!
         REAL(KIND = 8), DIMENSION(:,:), POINTER :: cec ! Cation exchange capacity, meq/100 g of solid
!
! ----------
! .......Real 3-D arrays for ion exchange     
! ----------
!
         REAL(KIND = 8), DIMENSION(:,:,:), POINTER :: xcads ! Exchanged species concentrations for grid blocks
!
! ----------
! .......Real 1-D arrays for surface complexation (adsorption)     
! ----------
!
         REAL(KIND = 8), DIMENSION(:), POINTER :: &
              sads, &
              psi,  &         
              tads, &        ! Total surface site
              supadn, &
              phip          ! Potential term
!
! ----------
! .......Real 2-D arrays for surface complexation (adsorption)     
! ----------
!
         REAL(KIND = 8), DIMENSION(:,:), POINTER :: d ! Concentrations of surface complexes, mol/kg H2O
!
!
      END MODULE Grid_ExchangeAdsorption
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Grid_ChemZones
!
!
         SAVE
! 
! ----------
! .......Integer 1-D arrays for Chemical zone number assigned to grid block     
! ----------
!
         INTEGER(KIND = 4), DIMENSION(:), POINTER :: &
                 izoneiw,&  ! Number of initial water zones given in chemical.inp
                 izonebw,&  ! Number of boundary (injection) water zones
                 izonem, &  ! Number of initial mineral zones
                 izoneg, &  ! Number of initial gas zones
                 izoned, &  ! Number of initial surface complex zones
                 izonex, &  ! Number of initial cation exchange zones
                 izonpp, &  ! Number of permeability-porosity relationship zones
                 izonekd    ! Number of initial Kd zone
!
!
      END MODULE Grid_ChemZones
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Grid_CoupledTerms
!
!
         SAVE
! 
! ----------
! .......Integer 1-D arrays     
! ----------
!
         INTEGER(KIND = 4), DIMENSION(:), POINTER :: &
           iskip1,   &
           ICALLCH
!
! ----------
! .......Real 2-D arrays for solving coupled reaction-transport equations     
! ----------
!
         REAL(KIND = 8), DIMENSION(:,:), POINTER :: &
           rhand,   &   ! Right-hand side term for transport of aqueous component
           rsource, &   ! Source/sink term for transport of aqueous component
           prev_c,  &   ! Concentration at the previous iteration
           prev_ut, &   ! Concentration from transport at previous iteration
           rsourceq    ! Eq. minerals/gases source terms
!
         REAL(KIND = 8), DIMENSION(:,:), POINTER :: &
           RSOURCEG    ! Gas source terms due to reactions
!
!
      END MODULE Grid_CoupledTerms
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Grid_GaseousSpecies
!
!
         SAVE
!
! ----------
! .......Real 1-D arrays for gaseous species     
! ----------
!
         REAL(KIND = 8), DIMENSION(:), POINTER :: &
           PFUGOLD2,&  ! Partial pressure for a gaseous species, bar
           rhandg      ! Right-hand side term for transport of gaseous species
!
! ----------
! .......Real 2-D arrays for gaseous species     
! ----------
!
         REAL(KIND = 8), DIMENSION(:,:), POINTER :: &
           GP0,    & ! Gaseous species concentrations at the previous time step, mol/dm3 medium
           GP,     & ! Gaseous species concentrations at the current  time step, mol/dm3 medium   
           QG,     & ! Gas influx from the boundary              
           PFUGOLD,& ! Partial pressure at the previous time step, bar
           PFUG      ! Partial pressure at the current time step, bar
! 
! ----------
! .......variable and arrays for grid blocks connected to constant pressure 
!..........boundary of gaseous species (such as oxygen diffusion from the atmosphere)     
! ----------
!
         INTEGER(KIND = 4) :: NELG  ! Total number of blocks connected to atmosphere
!
         INTEGER(KIND = 4), DIMENSION(:), POINTER :: IELG    ! Pointer
!
         REAL(KIND = 8), DIMENSION(:), POINTER :: DISG,  &! Distance, m
              PFUGB2   ! Partial pressure pressure of one gas, bar
!
         REAL(KIND = 8), DIMENSION(:,:), POINTER ::  PFUGB   ! Partial pressure, bar
!
!
      END MODULE Grid_GaseousSpecies
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Grid_ProsityPermChanges
!
!
         SAVE
!
! ----------
! .......Arrays for changes in porosity and permeability 
!..........due to mineral alteration   
! ----------
!
!........Integer 1-D array
!
         INTEGER(KIND = 4), DIMENSION(:), POINTER :: ikplaw  ! Permeability-porosity laws; 
!                             1: Carman-Kozeny, 
!                             3: cubic law, 
!                             4: modified Cubic Law, 
!                             5: Verma-Pruess 
!
!........Real 1-D arrays
!
         REAL(KIND = 8), DIMENSION(:), POINTER :: &
           aparpp, &    ! Fracture aperture for modified Cubic Law
           bparpp, &    ! Fracture spacing  for modified Cubic Law
           pcfact, &    ! Capillary pressure scaling factor
           phim         ! Store porosity for printout
!
!........Real 2-D arrays
!
         REAL(KIND = 8), DIMENSION(:,:), POINTER :: &
           perm0,  &    ! Previous permeability
           perm,   &    ! Current  permeability
           permm       ! Store permeability for printout
!
! ----------
! .......Clay swelling (simple model for geothermal reaearch) 
! ----------
!
!........Integer variable
!
         INTEGER(KIND = 4) :: iswell   ! Indicator of clay swelling
!
!........Real 1-D array
!
         REAL(KIND = 8), DIMENSION(:), POINTER :: vmin0      ! Initial mole volume
!
!........Real 2-D array
!
         REAL(KIND = 8), DIMENSION(:,:), POINTER :: vmin_old   ! Previous mole volume for all node
!
!
      END MODULE Grid_ProsityPermChanges
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Grid_SaltPrecipitates
!
!............Salt precipitates at dry grid blocks
!
!
         SAVE
!
! ----------
!........Integer 1-D array
! ----------
!
         INTEGER(KIND = 4), DIMENSION(:), POINTER :: IDRY ! Indicator of dry grid block for calculating salt precipitates
!
! ----------
!........Real 2-D array
! ----------
!
         REAL(KIND = 8), DIMENSION(:,:), POINTER :: & 
           adry, &
           adryr,&
           adryr0
!
         REAL(KIND = 8), DIMENSION(:,:), POINTER :: drypre               
!
!
      END MODULE Grid_SaltPrecipitates
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Grid_AreaReduction
!
!.................Area reduction from active fracture model
!
         SAVE
!
! ----------
!........Real 1-D arrays
! ----------
!
         REAL(KIND = 8), DIMENSION(:), POINTER :: &
           a_fm2,&    ! Advection area reduction (flow from F to M) calculated   
!                         from active fracture model
           a_fmd     ! Diffusion area reduction (Both sides)
!
         REAL(KIND = 8), DIMENSION(:), POINTER :: &
           a_fmr     ! Save modified active fracture area for reaction
!
!
      END MODULE Grid_AreaReduction
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Grid_ChemDiffusion
!
!
         SAVE
!
! ----------
!........Double precision variables for diffusion coefficients
! ----------
!
         REAL(KIND = 8) :: DIFUN       ! Diffusion coefficient (m2/s) for aqueous  species. 
!                   DIFUN is multiplied by the tortuosity, defined in 
!                   rock property block of the flow input, and liquid saturation. 
!                   Notice that if tortuosity in flow input is zero, the program 
!                   computes it from Millington and Quirk (1961).
         REAL(KIND = 8) :: DIFUNG      ! Diffusion coefficients (m2/s) of the medium for gaseous species. 
!                   If DIFUNG  <  0.0, the program computes gaseous diffusion 
!                   coefficients as function of temperature and pressure 
!                   according to Eq. A.1 in Appendix A of the manual.
!
! ----------
!........Real 1-D arrays for time stepping use to calculate diffusion time 
! ----------
!
         REAL(KIND = 8), DIMENSION(:), POINTER :: &
                         DIFUNT_L,  &        ! For liquid phase     
                         DIFUNT_G            ! For gas phase 
!
!
      END MODULE Grid_ChemDiffusion
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      MODULE Global_MassBalance
!
!
         SAVE
!
! ----------
!........Real 1-D arrays
! ----------
!
         REAL(KIND = 8), DIMENSION(:), POINTER :: &
           SOLUTINP, &  ! Input of chemical component to the flow system, mole
           SOLUTOUT, &  ! Output from the system, mole
           SOLUTINI, &  ! Initial amount of component in the aqueous phase, mole 
           SOLIDINI, &  ! Initial amount of component in the solid phase, mole 
           SOLUTNOW, &  ! Current amount of component in the aqueous phase, mole 
           SOLIDNOW     ! Current amount of component in the solid phase, mole 
!
!
      END MODULE Global_MassBalance
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
