!
!
!
      SUBROUTINE ReadDimension
! 
! ... Modules to be used 
! 
      USE MaxChem_Dimensions
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*           Read array dimension variables from a input file          *
!*                                                                     *
!*                   Version 1.0 - July 19, 2006                       *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE                                      
!
! -------
! ... Integer variables
! -------
! 
      INTEGER :: i
      INTEGER :: icall = 0
      SAVE       icall

! -------
! ... Character variables
! -------
! 
      CHARACTER(LEN = 100) :: label
      CHARACTER(LEN =  20) :: name
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>>=>=>=>=>>=>=>=>=>>=>=>=>=>=>=>=>
!
! 
! ----------
! ....Open units for input and files     
! ----------
!
      OPEN (UNIT=25,FILE='dimension.inp',STATUS='UNKNOWN')  
      OPEN (UNIT=26,FILE='dimension.out',STATUS='UNKNOWN') 
!
      WRITE(*,*) '   --> reading array dimension variables'
!
! 
! ----------
! ....Read title part     
! ----------
!
      do i=1,6
         read  (25,'(a100)')  label     
         write (26,'(a100)')  label     
      end do
!
      do i=1,2
         read  (25,'(a100)')  label     
      end do
!
!
! ----------
! ....Chemical system     
! ----------
!
!.....Maximum number of primary species
!
      read  (25,'(a20,i10,   a100)')  Name, MaxNpri, label     
      write (26,'(a20,i10,1x,a100)')  Name, MaxNpri, label     
!
!.....Maximum number of kinetic reactions among primary species
!.....Including aqueous and sorption kinetics and biodegradation
!
      read  (25,'(a20,i10,   a100)')  Name, MaxNtrx, label     
      write (26,'(a20,i10,1x,a100)')  Name, MaxNtrx, label     
!
!.....Maximum number of secondary species
!
      read  (25,'(a20,i10,   a100)')  Name, MaxNaqx, label     
      write (26,'(a20,i10,1x,a100)')  Name, MaxNaqx, label     
!
!.....Maximum number of aqueous species
!
      MaxNaqt     = MaxNpri + MaxNaqx
!
!.....Maximum number of minerals (equilibrium + kinetics) 
!
      read  (25,'(a20,i10,   a100)')  Name, MaxNmin, label     
      write (26,'(a20,i10,1x,a100)')  Name, MaxNmin, label     
!
!.....Maximum number of solid solutions 
!
      read  (25,'(a20,i10,   a100)')  Name, MaxNss,  label     
      write (26,'(a20,i10,1x,a100)')  Name, MaxNss,  label     
!
!.....Maximum number of components in each solid solution
!
      read  (25,'(a20,i10,   a100)')  Name, MaxCpss, label     
      write (26,'(a20,i10,1x,a100)')  Name, MaxCpss, label     
!
!.....Maximum number of gasesous species
!
      read  (25,'(a20,i10,   a100)')  Name, MaxNgas, label     
      write (26,'(a20,i10,1x,a100)')  Name, MaxNgas, label     
!
!.....Maximum number of surface complexes
!
      read  (25,'(a20,i10,   a100)')  Name, MaxNads, label     
      write (26,'(a20,i10,1x,a100)')  Name, MaxNads, label     
!
!.....Maximum number of exchangeable ions
!
      read  (25,'(a20,i10,   a100)')  Name, MaxNexc, label     
      write (26,'(a20,i10,1x,a100)')  Name, MaxNexc, label     
!
!.....Maximum number of exchangeable sites
!
      read  (25,'(a20,i10,   a100)')  Name, MXsites, label     
      write (26,'(a20,i10,1x,a100)')  Name, MXsites, label     
!
!.....Maximum number of species with linear Kd adsorption or/and decay
!
      read  (25,'(a20,i10,   a100)')  Name, MaxNkdd, label     
      write (26,'(a20,i10,1x,a100)')  Name, MaxNkdd, label     
!
!.....The dimension used for Chemical Jacobian Matrix
!
      MaxNnr      = MaxNpri + MaxNmin + MaxNgas + 2
!
! 
! ----------
! ....Chemical zones     
! ----------
!
!.....Maximum number of initial water types 
!
      read  (25,'(a20,i10,   a100)')  Name, MaxNiwtype, label     
      write (26,'(a20,i10,1x,a100)')  Name, MaxNiwtype, label     
!
!.....Maximum number of boundary (injection) water types
!
      read  (25,'(a20,i10,   a100)')  Name, MaxNbwtype, label     
      write (26,'(a20,i10,1x,a100)')  Name, MaxNbwtype, label     
!
!.....Maximum number of initial mineral zones
!
      read  (25,'(a20,i10,   a100)')  Name, MaxNmtype,  label     
      write (26,'(a20,i10,1x,a100)')  Name, MaxNmtype,  label     
!
!.....Maximum number of initial gas zones
!
      read  (25,'(a20,i10,   a100)')  Name, MaxNgtype,  label     
      write (26,'(a20,i10,1x,a100)')  Name, MaxNgtype,  label     
!
!.....Maximum number of prosity-permeability zones
!
      read  (25,'(a20,i10,   a100)')  Name, MaxNppzon,  label     
      write (26,'(a20,i10,1x,a100)')  Name, MaxNppzon,  label     
!
!.....Maximum number of surface complex zones
!
      read  (25,'(a20,i10,   a100)')  Name, MaxNdtype,  label     
      write (26,'(a20,i10,1x,a100)')  Name, MaxNdtype,  label     
!
!.....Maximum number of lineral Kd adsorption zones
!
      read  (25,'(a20,i10,   a100)')  Name, MaxNkdtype, label     
      write (26,'(a20,i10,1x,a100)')  Name, MaxNkdtype, label     
!
!.....Maximum number of exchange zones  
!
      read  (25,'(a20,i10,   a100)')  Name, MaxNxtype,  label     
      write (26,'(a20,i10,1x,a100)')  Name, MaxNxtype,  label     
!
! 
! ----------
! ....Database    
! ----------
!
!.....Maximum number of temperature points in the database  
!
      read  (25,'(a20,i10,   a100)')  Name, MaxNtmp,  label     
      write (26,'(a20,i10,1x,a100)')  Name, MaxNtmp,  label     
!
!.....Maximum number of ion pairs for using Pitzer ion-interaction model  
!
      read  (25,'(a20,i10,   a100)')  Name, MaxNpair,  label     
      write (26,'(a20,i10,1x,a100)')  Name, MaxNpair,  label     
!
! 
! ----------
! ....Mesh                                        ! Later deleted, read from TOUGH+    
! ----------
!
      MaxNum_Elem =  8000                         !!!!!! Tempor..!!!! Delete leter
      MaxNum_Conx = 20000  !                      !!!!!! Tempor..!!!! Delete leter
!
! 
! ----------
! ....End read     
! ----------
!
      do i=1,5
         read  (25,'(a100)')  label     
         write (26,'(a100)')  label     
      end do
!
! 
! ----------
! ....Close units     
! ----------
!
      CLOSE (unit = 25)    ! Input  file
      CLOSE (unit = 26)    ! Output file
!
!
      return
!
      end  SUBROUTINE ReadDimension
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE AllocateChemArrayMem
! 
! ... Modules to be used 
! 
      USE MaxChem_Dimensions
      USE Aqueous_Species
      USE Kinetics_Among_Primary_Species 
      USE Mineral_Phases
      USE Gaseous_Species
      USE Ion_Exchange
      USE Surface_Complexes
      USE Kd_Adsorption_Decay
      USE Batch_Chem_Jacobian
      USE WriteChemControl
      USE IonPairs_Index
      USE IonInteract_Parameters
      USE DerivativeTerms_Pitzer
      USE IniChemZoneProperties
!c      USE Grid_HydroTherm
!c      USE Grid_AqueousPhase
!c      USE Grid_MineralPhases
!c      USE Grid_ChemSteadyState
!c      USE Grid_ExchangeAdsorption
!c      USE Grid_ChemZones
!c      USE Grid_CoupledTerms
!c      USE Grid_GaseousSpecies
!c      USE Grid_ProsityPermChanges
!c      USE Grid_SaltPrecipitates
!c      USE Grid_AreaReduction
!c      USE Grid_ChemDiffusion
!c      USE Global_MassBalance
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*      Allocate array memory for reactive geochemical transport       *
!*                                                                     *
!*                  Version 1.0 - October 13, 2005                     *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE                                      
! 
! -------
! ... Integer variables
! -------
! 
      INTEGER :: i, MaxNnn
      INTEGER :: icall = 0
      SAVE       icall
! 
! -------
! ... Integer array
! -------
! 
      INTEGER, DIMENSION(1:100) :: ierror = 99999
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>>=>=>=>=>>=>=>=>=>>=>=>=>=>=>=>=>
!
!
!
!***********************************************************************
!*                                                                     *
!*                  Allocate primary species                           *
!*                                                                     *
!***********************************************************************
!
!
! ----------
! ....Allocate memory to derived-type arrays  
! ----------
!
      ALLOCATE(PSpe(MaxNpri), STAT=ierror(1)) 
!              PSpe = Primary Species array
      ALLOCATE(ASpe(MaxNaqt), STAT=ierror(2)) 
!              ASpe = Aqueous SPEcies array
!
! ----------
! ....Allocate memory to arrays within the derived type: PSpe  
! ----------
! 
      DO_PSpe: &
      DO i = 1,MaxNpri 
         ALLOCATE(PSpe(i)%du2  (MaxNpri), &
                  PSpe(i)%ds2  (MaxNpri), &
                  PSpe(i)%dgamp(MaxNpri),  STAT = ierror(3))
         ALLOCATE(PSpe(i)%dr   (MaxNpri), &
                  PSpe(i)%drx  (MaxNpri),  STAT = ierror(4))
      END DO DO_PSpe 
!
!
!***********************************************************************
!*                                                                     *
!*         Allocate kinetic reactions among primary species            *
!*                                                                     *
!***********************************************************************
!
!
! ----------
! ....Allocate memory to derived-type arrays  
! ----------
!
      ALLOCATE(PKin(MaxNtrx), STAT=ierror(5)) 
!              PKin = Primary Kinetics
!
!
!***********************************************************************
!*                                                                     *
!*                Allocate secondary species                           *
!*                                                                     *
!***********************************************************************
!
!
      if (MaxNaqx > 0)   then
!
! ----------
! ...... Allocate memory to derived-type arrays for secondary species  
! ----------
!
         ALLOCATE(SSpe(MaxNaqx), STAT=ierror(6)) 
!                 SSpe = Secondary SPEcies array
!
! ----------
! ...... Allocate memory to arrays within the derived type: SSpe and ASpe  
! ----------
! 
         DO_SSpe: &
         DO i = 1,MaxNaqx 
            ALLOCATE(SSpe(i)%icps (MaxNpri), STAT = ierror(7))
            ALLOCATE(SSpe(i)%stqs (MaxNpri), STAT = ierror(8))
            ALLOCATE(SSpe(i)%dgams(MaxNpri), STAT = ierror(9))
         END DO DO_SSpe  
!
         DO_ASpe: &
         DO i = 1,MaxNaqt 
            ALLOCATE(ASpe(i)%dgamt(MaxNpri), STAT = ierror(10))
         END DO DO_ASpe 
!
      end if
!
!    
!***********************************************************************
!*                                                                     *
!*                 Allocate minerals in the system                     *
!*                                                                     *
!***********************************************************************
!
!
      if (MaxNmin > 0)   then
!
! ----------
! ...... Allocate memory to derived-type arrays for minerals 
! ----------
!
         ALLOCATE(MGen(MaxNmin), STAT=ierror(11)) 
!                 MGEN = General data related to mineral phases
!
         ALLOCATE(MKin(MaxNmin), STAT=ierror(12)) 
!                 MKin = Mineral dissolution and precipitation kinetic parameters
!
         ALLOCATE(isalt(0:MaxNmin), STAT=ierror(13)) 
!
         ALLOCATE(ncpss(MaxNss),       STAT=ierror(14)) 
!
! ----------
! ...... Allocate memory to arrays within the derived type: MGen  
! ----------
!  
         DO_MGen: &
         DO i = 1,MaxNmin 
            ALLOCATE(MGen(i)%icpm (MaxNpri), &
                     MGen(i)%stqm (MaxNpri), STAT = ierror(15))
!            ALLOCATE(MGen(i)%stqm (MaxNpri), STAT = ierror(16))
         END DO DO_MGen  
!
!
         ALLOCATE(icpss(MaxNss,MaxCpss), STAT = ierror(16))
!
      end if
!
! 
!***********************************************************************
!*                                                                     *
!*               Allocate gaseous species in the system                *
!*                                                                     *
!***********************************************************************
!
!
      if (MaxNgas > 0) then
!
! ----------
! ...... Allocate memory to derived-type arrays  
! ----------
!
         ALLOCATE(GSpe(MaxNgas), STAT=ierror(17)) 
!                 GSpe = Gaseous SPEcies array
!
! ----------
! ...... Allocate memory to arrays within the derived type: GSpe  
! ----------
! 
         DO_GSpe: &
         DO i = 1,MaxNgas 
            ALLOCATE(GSpe(i)%icpg (MaxNpri), &
                     GSpe(i)%stqg (MaxNpri), STAT = ierror(18))
!            ALLOCATE(GSpe(i)%stqg (MaxNpri), STAT = ierror(20))
         END DO DO_GSpe  
!
      end if
!
! 
!***********************************************************************
!*                                                                     *
!*               Allocate surface complexes in the system              *
!*                                                                     *
!***********************************************************************
!
!
      if (MaxNads > 0)   then 
!
! ----------
! ...... Allocate memory to derived-type arrays  
! ----------
!
         ALLOCATE(DSpe(MaxNads), STAT=ierror(19)) 
!                 DSpe = Surface SPEcies array
!
! ----------
! ...... Allocate memory to arrays within the derived type: DSpe 
! ----------
! 
         DO_DSpe: &
         DO i = 1,MaxNads 
            ALLOCATE(DSpe(i)%stqd (MaxNpri), &
                     DSpe(i)%dcd  (MaxNpri), STAT = ierror(20))
!            ALLOCATE(DSpe(i)%dcd  (MaxNpri), STAT = ierror(23))
         END DO DO_DSpe  
!
      end if
!
! 
!***********************************************************************
!*                                                                     *
!*           Allocate species (primary) with linear Kd and decay       *
!*                                                                     *
!***********************************************************************
!
!
      if (MaxNkdd > 0)               then
!
! ----------
! ...... Allocate memory to derived-type arrays  
! ----------
!
         ALLOCATE(KSpe(MaxNkdd), STAT=ierror(21)) 
!                 KSpe = Kd and decay SPEcies array
!
      end if
!
! 
!***********************************************************************
!*                                                                     *
!*               Allocate exchange cations of the system               *
!*                                                                     *
!***********************************************************************
!
!
      if (MaxNexc > 0)               then
! 
! ----------
! ...... Allocate memory to derived-type arrays  
! ----------
! 
         ALLOCATE(XSpe(MaxNexc), STAT=ierror(22)) 
!              XSpe = eXahanged SPEcies array
!
! ----------
! ...... Allocate memory to arrays within the derived type: XSpe  
! ----------
! 
         DO_XSpe: &
         DO i = 1,MaxNexc 
            ALLOCATE(XSpe(i)%stqx(MaxNpri), STAT = ierror(23))
            ALLOCATE(XSpe(i)%dcx (MaxNpri), STAT = ierror(24))
         END DO DO_XSpe 
!
!
! ----------
! ...... Allocate memory to arrays for multi-site exchange 
! ----------
! 
         ALLOCATE(cecM(MXsites), STAT=ierror(25))   ! Batch CEC for multi-sites
!
         ALLOCATE(ekxM(MXsites, MaxNexc), &         ! Selectivity
                  cxM (MXsites, MaxNexc), &         ! Exchanged species concentrations
                  STAT=ierror(26)) 
!
         ALLOCATE(dcxM(MXsites, MaxNexc, MaxNpri),& ! Devivatives
                  STAT=ierror(27)) 
!
      end if
!
! 
!***********************************************************************
!*                                                                     *
!*        Allocate memory to arrays for the chemical Jacobian          *
!*                                                                     *
!***********************************************************************
!
!
      ALLOCATE(isat (MaxNmin), STAT=ierror(28)) 
      ALLOCATE(mout (MaxNmin), STAT=ierror(29)) 
!
      ALLOCATE(isatg(MaxNgas), STAT=ierror(30)) 
      ALLOCATE(moutg(MaxNgas), STAT=ierror(31)) 
!
      ALLOCATE(amat(MaxNnr,MaxNnr), STAT=ierror(32)) 
      ALLOCATE(bmat(MaxNnr),     STAT=ierror(33)) 
      ALLOCATE(indx(MaxNnr),     STAT=ierror(34))                     !!! 
!
      ALLOCATE(fp(MaxNpri),      STAT=ierror(35)) 
!
! 
!***********************************************************************
!*                                                                     *
!*        Allocate chemical composition of initial and                 *
!*                      boundary (injection) water zones               *
!*                                                                     *
!***********************************************************************
!
!
      MaxNnn =   MaxNiwtype &  !  Maximum number of initial water types 
               + MaxNbwtype    !  Maximum number of boundary (injection) water types
!
      ALLOCATE(icon_iniZ(MaxNnn, MaxNpri), STAT=ierror(36)) 
      ALLOCATE(U_iniZ   (MaxNnn, MaxNpri), STAT=ierror(37)) 
      ALLOCATE(c_iniZ   (MaxNnn, MaxNaqt), STAT=ierror(38)) 
      ALLOCATE(pH_iniZ(MaxNnn), aw_iniZ(MaxNnn), STAT=ierror(39)) 
      ALLOCATE(Tc_iniZ  (MaxNnn), STAT=ierror(40)) 
!
      ALLOCATE(pfug_iniZ(MaxNiwtype, MaxNgas), STAT=ierror(41)) 
!
      if (MaxNbwtype > 0)   then
         ALLOCATE(Uboud(MaxNbwtype, MaxNpri), STAT=ierror(42)) 
      end if
!
!
!***********************************************************************
!*                                                                     *
!*                   Allocate initial mineral zones                    *
!*                                                                     *
!***********************************************************************
!
!
      if (MaxNmtype > 0)   then
         ALLOCATE(vol_iniZ   (MaxNmtype, MaxNmin), STAT=ierror(43)) 
         ALLOCATE(rad_iniZ   (MaxNmtype, MaxNmin), STAT=ierror(44)) 
         ALLOCATE(amin_iniZ  (MaxNmtype, MaxNmin), STAT=ierror(45)) 
         ALLOCATE(imflg_iniZ (MaxNmtype, MaxNmin), STAT=ierror(46)) 
         ALLOCATE(sumvol_iniZ(MaxNmtype),          STAT=ierror(47)) 
      end if
!
! 
!***********************************************************************
!*                                                                     *
!*                    Allocate initial gas zones                       *
!*                                                                     *
!***********************************************************************
!
!
      if (MaxNgtype > 0)   then
         ALLOCATE(cg_iniZ(MaxNgtype,MaxNgas), STAT=ierror(48)) 
      end if
!
! 
!***********************************************************************
!*                                                                     *
!*               Allocate permeability-porosity law zones              *
!*                                                                     *
!***********************************************************************
!
!
      if (MaxNppzon > 0)   then
         ALLOCATE(ipptyp_iniZ(MaxNppzon), STAT=ierror(49)) 
         ALLOCATE(apppar_iniZ(MaxNppzon), STAT=ierror(50)) 
         ALLOCATE(bpppar_iniZ(MaxNppzon), STAT=ierror(51)) 
      end if
!
! 
!***********************************************************************
!*                                                                     *
!*                Allocate surface complexation zones                  *
!*                                                                     *
!***********************************************************************
!
!
      if (MaxNdtype > 0)   then
         ALLOCATE(sadsdum_iniZ(MaxNdtype), STAT=ierror(52)) 
         ALLOCATE(tads_iniZ(MaxNdtype),    STAT=ierror(53)) 
      end if
!
! 
!***********************************************************************
!*                                                                     *
!*                 Allocate linear Kd adsorption zones                 *
!*                                                                     *
!***********************************************************************
!
!
      if (MaxNkdtype > 0)   then
         ALLOCATE(sden(MaxNkdtype,MaxNpri), STAT = ierror(54))
         ALLOCATE(vkd (MaxNkdtype,MaxNpri), STAT = ierror(55))
      end if
!
! 
!***********************************************************************
!*                                                                     *
!*                  Allocate cation exchange zones                     *
!*                                                                     *
!***********************************************************************
!
!
      if (MaxNxtype > 0)   then
         ALLOCATE(cec_iniZ(MaxNxtype, MXsites), STAT=ierror(56)) 
      end if
!
!
!***********************************************************************
!*                                                                     *
!*                 Allocate writing control arrays                     *
!*                                                                     *
!***********************************************************************
!
! 
      ALLOCATE(IWNOD (MaxNaqt), STAT=ierror(56)) 
      ALLOCATE(IWCOM (MaxNaqt), STAT=ierror(57)) 
!
      ALLOCATE(IWMIN (MaxNmin), STAT=ierror(58)) 
!
      ALLOCATE(elemw (100)   ) ! Name of grid blocks for time evolution

!
!***********************************************************************
!*                                                                     *
!*    Allocate hydrological and thermal properties for grib blocks     *
!*                                                                     *
!***********************************************************************
!
! 
!      ALLOCATE(SLOLD (MaxNum_Elem),    ! Liquid saturation at previou time step
!     &         SGOLD (MaxNum_Elem),    ! Gas saturation at previou time step
!     &         SL1   (MaxNum_Elem),    ! Liquid saturation at previous time step
!     &         SG1   (MaxNum_Elem),    ! Gas saturation at the current time step
!     &         PHI0  (MaxNum_Elem),    ! Porosity at t = 0
!     &         PHIOLD(MaxNum_Elem),    ! Porosity at the previous time step
!     &         tc    (MaxNum_Elem),    ! Temperature at current time step, oC
!     &         tcold (MaxNum_Elem),    ! Temperature at previous time step, oC
!     &         Tc_EOS9(MaxNum_Elem),   ! Initial temperature using EOS9, oC
!     &         densw  (MaxNum_Elem),   ! Density of water, kg/m^3 
!     &         fugcoeCO2(MaxNum_Elem), ! CO2 fugacity coefficient from ECO2 and ECO2N
!     &         Rh(MaxNum_Elem),        ! Relative humidity obtained from EOS
!     &         STAT=ierror(59))
!
!
!***********************************************************************
!*                                                                     *
!*              Allocate aqueous phase for grib blocks                 *
!*                                                                     *
!***********************************************************************
!
! 
!c      ALLOCATE(ph       (MaxNum_Elem),        ! pH values
!c     &         str_node (MaxNum_Elem),        ! Ionic strength
!c     &         aw       (MaxNum_Elem),        ! Water activity
!c     &         STAT=ierror(60))
!
!c      ALLOCATE(
!c     &   c (MaxNum_Elem, MaxNaqt), ! Concentrations of all aqueous species, mol/kg
!c     &   STAT=ierror(61))    
!
!c      ALLOCATE(
!c     &   ut (MaxNum_Elem, MaxNpri),    ! Current dissolved concent. of component, mol/l 
!c     &   utem (MaxNum_Elem, MaxNpri),  ! Tempo. store ut 
!c     &   utold (MaxNum_Elem, MaxNpri), ! Previous dissolved concent. of component 
!c     &   STAT=ierror(62))    
!
!
!***********************************************************************
!*                                                                     *
!*              Allocate mineral phases for grib blocks                *
!*                                                                     *
!***********************************************************************
!
! 
      if (MaxNmin > 0)   then
!
!c         ALLOCATE(
!c     &   pinit (MaxNum_Elem,MaxNmin+1), ! Mineral abundance at t = 0, mol/dm3 medium
!c     &   pre0  (MaxNum_Elem,MaxNmin),   ! Mineral abundance at the previous time step, mol/dm**3 medium
!c     &   pre   (MaxNum_Elem,MaxNmin),   ! Mineral abundance at the current  time step, mol/dm**3 medium
!c     &   amin  (MaxNum_Elem,MaxNmin),   ! Reactive surface area, m**2/m**3 mineral
!c     &   rad   (MaxNum_Elem,MaxNmin),   ! Radius of mineral grain, m
!c     &   sim   (MaxNum_Elem,MaxNmin),   ! Mineral saturation index (log10(Q/K))
!c     &   rkin0 (MaxNum_Elem,MaxNmin),   ! Mineral dissolution rate (negative values for precipitation) 
!c!                                         at the previous time step
!c     &   rkin  (MaxNum_Elem,MaxNmin),   ! Mineral dissolution rate at the current time step
!c     &   STAT=ierror(63))    
!c!
!c         ALLOCATE( 
!c     &   Rh2o  (MaxNum_Elem),           ! H2O generated from reactions
!c     &   Rco2  (MaxNum_Elem),           ! CO2 generated from reactions
!c     &   SMco2 (MaxNum_Elem),           ! CO2 trapped in solid phase
!c     &   STAT=ierror(64))    
!c!
      end if

!
!***********************************************************************
!*                                                                     *
!*             Allocate arrays for chemical steady-state               *
!*                                                                     *
!***********************************************************************
!
! 
      if (MaxNmin > 0)   then
!c         ALLOCATE( 
!c     &   DELTAP (MaxNum_Elem,MaxNmin),  ! Mineral abundance increment for chemical steady-state
!c     &   DRATE  (MaxNum_Elem,MaxNmin),  ! Rate for chemical qusi-steady state
!c     &   STAT=ierror(65))    
      end if
!
!
!***********************************************************************
!*                                                                     *
!*           Ion exchange and adsorption for grid blocks               *
!*                                                                     *
!***********************************************************************
!
! 
      if (MaxNexc > 0)   then
!
!c         ALLOCATE( 
!c     &   cec (MaxNum_Elem, MXsites),     ! Cation exchange capacity, meq/100 g of solid
!c     &   STAT=ierror(66))    
!c!
!c         ALLOCATE(
!ccc     &   xcads (MaxNum_Elem, MaxNexc), ! Concentrations of exchangeable cations, mol/kg H2O
!ccc     &   STAT=ierror(67))    
!c     &   xcads (MaxNum_Elem, MXsites, MaxNexc), ! Concentrations of exchangeable cations, mol/kg H2O
!c     &   STAT=ierror(67))    
!c!
!c      end if
!c!
!c!
!c      if (MaxNads > 0)   then
!c!
!c         ALLOCATE( 
!c     &   sads   (MaxNum_Elem),
!c     &   psi    (MaxNum_Elem),            
!c     &   tads   (MaxNum_Elem),         ! Total surface site
!c     &   supadn (MaxNum_Elem),
!c     &   phip   (MaxNum_Elem),         ! Potential term
!c     &   STAT=ierror(68))    
!c!
!c         ALLOCATE(
!c     &   d (MaxNum_Elem, MaxNads),    ! Concentrations of surface complexes, mol/kg H2O
!c     &   STAT=ierror(69))    
!c!
      end if
!
!
!***********************************************************************
!*                                                                     *
!*                   Chemical zones for grid blocks                    *
!*                                                                     *
!***********************************************************************
!
! 
!c      ALLOCATE( 
!c     &   izoneiw (MaxNum_Elem),  ! Number of initial water zones given in chemical.inp
!c     &   izonebw (MaxNum_Elem),  ! Number of boundary (injection) water zones
!c     &   izonem  (MaxNum_Elem),  ! Number of initial mineral zones
!c     &   izoneg  (MaxNum_Elem),  ! Number of initial gas zones
!c     &   izoned  (MaxNum_Elem),  ! Number of initial surface complex zones
!c     &   izonex  (MaxNum_Elem),  ! Number of initial cation exchange zones
!c     &   izonpp  (MaxNum_Elem),  ! Number of permeability-porosity relationship zones
!c     &   izonekd (MaxNum_Elem),  ! Number of initial Kd zone
!c     &   STAT=ierror(70))    
!
!
!***********************************************************************
!*                                                                     *
!*                 Coupled reaction-transport terms                    *
!*                                                                     *
!***********************************************************************
!
! 
!c      ALLOCATE( 
!c     &   iskip1  (MaxNum_Elem),   
!c     &   icallch (MaxNum_Elem),
!c     &   STAT=ierror(71))    
!c!
!c      ALLOCATE(
!c     &   rhand   (MaxNum_Elem,MaxNpri),    ! Right-hand side term for transport of aqueous component
!c     &   rsource (MaxNum_Elem,MaxNpri),    ! Source/sink term for transport of aqueous component
!c     &   prev_c  (MaxNum_Elem,MaxNpri),    ! Concentration at the previous iteration
!c     &   prev_ut (MaxNum_Elem,MaxNpri),    ! Concentration from transport at previous iteration
!c     &   rsourceq(MaxNum_Elem,MaxNpri),    ! Eq. minerals/gases source terms
!c     &   STAT=ierror(72))    
!c!
!c      if (MaxNgas > 0)   then
!c         ALLOCATE( 
!c     &   RSOURCEG (MaxNum_Elem,MaxNgas),   ! Gas source terms due to reactions
!c     &   STAT=ierror(73))    
!c      end if
!
!
!***********************************************************************
!*                                                                     *
!*                 Gaseous species for grid blocks                     *
!*                                                                     *
!***********************************************************************
!
!
      if (MaxNgas > 0)   then
!	 
!c         ALLOCATE( 
!c     &   PFUGOLD2 (MaxNum_Elem),    ! Partial pressure (bar) for one gas    
!c     &   rhandg   (MaxNum_Elem),    ! Right-hand side term for transport of gaseous species
!c     &   STAT=ierror(74))    
!c!
!c         ALLOCATE( 
!c     &   GP0 (MaxNum_Elem,MaxNgas),     ! Gaseous species concentrations at the previous time step, mol/dm3 medium
!c     &   GP  (MaxNum_Elem,MaxNgas),     ! Gaseous species concentrations at the current  time step, mol/dm3 medium   
!c     &   QG  (MaxNum_Elem,MaxNgas),     ! Gas influx from the boundary  
!c     &   PFUGOLD(MaxNum_Elem,MaxNgas),  ! Partial pressure at the previous time step, bar
!c     &   PFUG (MaxNum_Elem,MaxNgas),    ! Partial pressure at the current time step, bar
!c     &   STAT=ierror(75))    
!c!
!c!
!c         ALLOCATE(
!c     &   IELG (100),               ! Pointer
!c     &   STAT=ierror(76))    
!c!
!c         ALLOCATE( 
!c     &   DISG   (100),             ! Distance, m
!c     &   PFUGB2 (100),             ! Partial pressure pressure of one gas, bar
!c     &   STAT=ierror(77))    
!c!
!c         ALLOCATE( 
!c     &   PFUGB (100,MaxNgas),           ! Partial pressure, bar
!c     &   STAT=ierror(78))    
!c!
      end if
!
!
!***********************************************************************
!*                                                                     *
!*              Porosity and permeability changes                      *
!*                                                                     *
!***********************************************************************
!
!                                                        Check only allocate MaxNmin > 0  !! 
!c      ALLOCATE( 
!c     &   ikplaw (MaxNum_Elem),       ! Permeability-porosity laws; 
!c     &   STAT=ierror(79))    
!c!
!c      ALLOCATE( 
!c     &   aparpp (MaxNum_Elem),       ! Fracture aperture for modified Cubic Law
!c     &   bparpp (MaxNum_Elem),       ! Fracture spacing  for modified Cubic Law
!c     &   pcfact (MaxNum_Elem),       ! Capillary pressure scaling factor
!c     &   phim   (MaxNum_Elem),       ! Store porosity for printout
!c     &   STAT=ierror(80))    
!c!
!c      ALLOCATE( 
!c     &   perm0 (3, MaxNum_Elem),     ! Previous porosity
!c     &   perm  (3, MaxNum_Elem),     ! Current porosity
!c     &   permm (3, MaxNum_Elem),     ! Store porosity for printout
!c     &   STAT=ierror(81))    
!c!
!c      if (MaxNmin > 0)    then
!c!
!c         ALLOCATE( 
!c     &   vmin0 (MaxNmin),            ! Initial mole volume
!c     &   STAT=ierror(82))    
!c!
!c         ALLOCATE( 
!c     &   vmin_old (MaxNum_Elem, MaxNmin),  ! Previous mole volume for all grid blocks
!c     &   STAT=ierror(83))    
!c!
!c      end if
!
!
!***********************************************************************
!*                                                                     *
!*              Dalt precipitates at dry grid blocks                   *
!*                                                                     *
!***********************************************************************
!
!
!c      ALLOCATE( 
!c     &   IDRY (MaxNum_Elem), ! Indicator of dry grid block for calculating salt precipitates
!c     &   STAT=ierror(84))    
!c!
!c      ALLOCATE( 
!c     &   adry   (MaxNum_Elem, MaxNpri),
!c     &   adryr  (MaxNum_Elem, MaxNpri),
!c     &   adryr0 (MaxNum_Elem, MaxNpri),
!c     &   STAT=ierror(85))    
!c!
!c      if (MaxNmin > 0)    then
!c         ALLOCATE( 
!c     &   drypre (MaxNum_Elem, MaxNmin),    
!c     &   STAT=ierror(86))    
!c      end if
!
!
!***********************************************************************
!*                                                                     *
!*            Area reduction due to active fracture model              *
!*                                                                     *
!***********************************************************************
!
!
!c      ALLOCATE( 
!c     &   a_fm2 (MaxNum_Conx),    ! Advection area reduction (flow from F to M) calculated
!c!                               from active fracture model
!c     &   a_fmd (MaxNum_Conx),    ! Diffusion area reduction (Both sides)
!c     &   STAT=ierror(87))    
!c!
!c      ALLOCATE( 
!c     &   a_fmr(MaxNum_Elem),      ! Save modified active fracture area for reaction
!c     &   STAT=ierror(88))    
!c!
!
!***********************************************************************
!*                                                                     *
!*                Allocate chemical diffusion terms                    *
!*                                                                     *
!***********************************************************************
!
!c!
!c      ALLOCATE( 
!c     &   DIFUNT_L (MaxNum_Conx),         ! For liquid phase            
!c     &   DIFUNT_G (MaxNum_Conx),         ! For gas phase 
!c     &   STAT=ierror(89))    
!c!
!
!***********************************************************************
!*                                                                     *
!*                       Global mass balance                           *
!*                                                                     *
!***********************************************************************
!
!
!c      ALLOCATE ( 
!c     &   SOLUTINP (MaxNpri),   ! Input of chemical component to the flow system, mole
!c     &   SOLUTOUT (MaxNpri),   ! Output from the system, mole
!c     &   SOLUTINI (MaxNpri),   ! Initial amount of component in the aqueous phase, mole 
!c     &   SOLIDINI (MaxNpri),   ! Initial amount of component in the solid phase, mole 
!c     &   SOLUTNOW (MaxNpri),   ! Current amount of component in the aqueous phase, mole 
!c     &   SOLIDNOW (MaxNpri),   ! Current amount of component in the solid phase, mole 
!c     &   STAT=ierror(90))    
!
!
!***********************************************************************
!*                                                                     *
!*       Allocate arrays for the Pitzer ion-interaction model          *
!*                                                                     *
!***********************************************************************
!
!
!c      IF_PitzerModel: 
!c     &if (mopr(9) > 0)   then
!c!
!c!------------
!c!........Integer Arrays in the Pitzer Model
!c!------------
!c!
!c         ALLOCATE (
!c     &       ippos (MaxNaqt),       
!c     &       ipneg (MaxNaqt),         
!c     &       ipnue (MaxNaqt),         
!c     &       itpos (MaxNpair),        
!c     &       itneg (MaxNpair),        
!c     &       itneu (MaxNpair),
!c     &                      STAT=ierror(91))        
!c!
!c         ALLOCATE (
!c     &       ipair (MaxNpair,2),      
!c     &       icc   (MaxNpair,2),        
!c     &       iaa   (MaxNpair,2),        
!c     &       icca  (MaxNpair,3),       
!c     &       icaa  (MaxNpair,3),       
!c     &       inca  (MaxNpair,3),       
!c     &       ipp   (MaxNaqt, MaxNaqt),      
!c     &                      STAT=ierror(92))        
!c!
!c!------------
!c!........Pitzer Ionic Interaction Parameters
!c!------------
!c!
!c         ALLOCATE (
!c     &       beta0        (MaxNpair),         
!c     &       beta1        (MaxNpair),         
!c     &       beta2        (MaxNpair),         
!c     &       cmx          (MaxNpair),           
!c     &       ramdapos     (MaxNpair),      
!c     &       ramdaneg     (MaxNpair),      
!c     &       theta_cc     (MaxNpair),      
!c     &       theta_aa     (MaxNpair),      
!c     &       psi_cca      (MaxNpair),       
!c     &       psi_caa      (MaxNpair),       
!c     &       zeta_nca     (MaxNpair),      
!c     &       alfamx1      (MaxNpair),       
!c     &       alfamx2      (MaxNpair),       
!c     &       dbmx_dstr    (MaxNpair),     
!c     &       ddbmx_dstr   (MaxNpair),    
!c     &       dbphimx_dstr (MaxNpair),  
!c     &                      STAT=ierror(93))        
!c!
!c         ALLOCATE (
!c     &       betat0    (MaxNpair,4),         
!c     &       betat1    (MaxNpair,4),         
!c     &       betat2    (MaxNpair,4),         
!c     &       cmxt      (MaxNpair,4),           
!c     &       ramdapost (MaxNpair,4),      
!c     &       ramdanegt (MaxNpair,4),      
!     &       thetat_cc (MaxNpair,4),      
!c     &       thetat_aa (MaxNpair,4),      
!c     &       psit_cca  (MaxNpair,4),       
!c     &       psit_caa  (MaxNpair,4),       
!c     &       zetat_nca (MaxNpair,4),      
!c     &                      STAT=ierror(94))        
!
!------------
!........Derivative terms and intermediate derivative terms
!------------
!
!c         ALLOCATE (
!c     &       bphimx (MaxNpair),          
!c     &       bmx    (MaxNpair),             
!c     &                      STAT=ierror(95))        
!c!
!c         ALLOCATE (
!c     &       dbphimx (MaxNpair,MaxNpri),   
!c     &       dbmx    (MaxNpair,MaxNpri),      
!c     &       ddbmx   (MaxNpair,MaxNpri),     
!c     &       dgamt   (MaxNaqt, MaxNpri),      
!c     &       dgamp   (MaxNpri, MaxNpri),      
!c     &       dgams   (MaxNaqx, MaxNpri),      
!c     &                      STAT=ierror(96))        
!c!
!c      end if  IF_PitzerModel
! 
!
!***********************************************************************
!*                                                                     *
!*                      Check memory alocation                         *
!*                                                                     *
!***********************************************************************
!
!
      WRITE(38,*)      
      WRITE(38,*) '   --- Memory Aallocation ---  '
      WRITE(38,*)
!
      DO_CheckChemInput: &
      DO I=1,96
! 
         IF(ierror(I) == 0) THEN
            WRITE(38,6015) I
         ELSE IF(ierror(I) == 99999) THEN
            CONTINUE
         ELSE
            WRITE(38,6020) I
            STOP
         END IF
! 
      END DO DO_CheckChemInput
!
!
 6015 FORMAT(T2,' Memory allocation at point ',i3, &                   
                ' in subroutine "ReadChemInput" was successful')
!
 6020 FORMAT(//,20('ERROR-'),//, &                                         
             T2,' Memory allocation at point ',i3, &                       
                ' in subroutine "ReadChemInput" was unsuccessful',//, &  
             T2,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!', &       
             //,20('ERROR-'))
!
      RETURN

      END SUBROUTINE AllocateChemArrayMem
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE DellocateChemZoneArrays
! 
! ... Modules to be used 
! 
      USE MaxChem_Dimensions
      USE IniChemZoneProperties
!c      USE Grid_ChemZones
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*        Dellocate array memory related initial chemical zones        *
!*                                                                     *
!*                  Version 1.0 - October 28, 2005                     *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE                                      
! 
! -------
! ... Integer variables
! -------
! 
      INTEGER :: i
      INTEGER :: icall = 0
      SAVE       icall
! 
! -------
! ... Integer array
! -------
! 
      INTEGER, DIMENSION(1:100) :: ierror = 99999
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>>=>=>=>=>>=>=>=>=>>=>=>=>=>=>=>=>
!
! 
!***********************************************************************
!*                                                                     *
!*       Dellocate chemical composition of initial water zones         *
!*                                                                     *
!***********************************************************************
!
!
      DEALLOCATE(icon_iniZ, STAT=ierror(1)) 
      DEALLOCATE(U_iniZ,    STAT=ierror(2)) 
      DEALLOCATE(c_iniZ,    STAT=ierror(3)) 
      DEALLOCATE(pH_iniZ, aw_iniZ,  STAT=ierror(4)) 
      DEALLOCATE(Tc_iniZ,   STAT=ierror(5)) 
!
      DEALLOCATE(pfug_iniZ, STAT=ierror(6)) 
!
!
!***********************************************************************
!*                                                                     *
!*                   DEALLOCATE initial mineral zones                  *
!*                                                                     *
!***********************************************************************
!
!
      if (MaxNmtype > 0)   then
         DEALLOCATE(vol_iniZ,    STAT=ierror(7)) 
         DEALLOCATE(rad_iniZ,    STAT=ierror(8)) 
         DEALLOCATE(amin_iniZ,   STAT=ierror(9)) 
         DEALLOCATE(imflg_iniZ,  STAT=ierror(10)) 
         DEALLOCATE(sumvol_iniZ, STAT=ierror(11)) 
      end if
!
! 
!***********************************************************************
!*                                                                     *
!*                    DEALLOCATE initial gas zones                     *
!*                                                                     *
!***********************************************************************
!
!
      if (MaxNgtype > 0)   then
         DEALLOCATE(cg_iniZ, STAT=ierror(12)) 
      end if
!
! 
!***********************************************************************
!*                                                                     *
!*               DEALLOCATE permeability-porosity law zones            *
!*                                                                     *
!***********************************************************************
!
!
      if (MaxNppzon > 0)   then
         DEALLOCATE(ipptyp_iniZ, STAT=ierror(13)) 
         DEALLOCATE(apppar_iniZ, STAT=ierror(14)) 
         DEALLOCATE(bpppar_iniZ, STAT=ierror(15)) 
      end if
!
! 
!***********************************************************************
!*                                                                     *
!*                DEALLOCATE surface complexation zones                *
!*                                                                     *
!***********************************************************************
!
!
      if (MaxNdtype > 0)   then
         DEALLOCATE(sadsdum_iniZ, STAT=ierror(16)) 
         DEALLOCATE(tads_iniZ,    STAT=ierror(17)) 
      end if
!
! 
!***********************************************************************
!*                                                                     *
!*                  DEALLOCATE cation exchange zones                   *
!*                                                                     *
!***********************************************************************
!
!
      if (MaxNxtype > 0)   then
         DEALLOCATE(cec_iniZ, STAT=ierror(18)) 
      end if
!
!
!***********************************************************************
!*                                                                     *
!*                   Chemical zones for grid blocks                    *
!*                                                                     *
!***********************************************************************
!
! 
!c      DEALLOCATE( 
!c     &   izoneiw,  ! Number of initial water zones given in chemical.inp
!c     &   izonem,   ! Number of initial mineral zones
!c     &   izoneg,   ! Number of initial gas zones
!c     &   izoned,   ! Number of initial surface complex zones
!c     &   izonex,   ! Number of initial cation exchange zones
!c     &   STAT=ierror(19))    
!
!
!***********************************************************************
!*                                                                     *
!*        ENSURE PROPER MEMORY DEALLOCATION - STOP IF PROBEMS          *
!*                                                                     *
!***********************************************************************
! 
!
      WRITE(38,*)      
      WRITE(38,*) '   --- Memory Deallocation ---  '
      WRITE(38,*)
!
      DO_DeChemCheck: &
      DO i=1,19 
!
         IF (ierror(i) == 0) THEN
            WRITE(38,6001) i
         ELSE IF(ierror(i) == 99999) THEN
            CONTINUE
         ELSE
            WRITE(38,6002) i   ! Unsuccesful deallocation
            STOP
         END IF
! 
      END DO DO_DeChemCheck
!
!
 6001 FORMAT(T2,' Memory deallocation at point ',i3, &
                ' in subroutine "DellocateChemZoneArrays" was', &
                ' successful')
!
 6002 FORMAT(//,20('ERROR-'),//, &
            T2,' Memory deallocation at point ',i3,&
               ' in subroutine "DellocateChemZoneArrays" was',&
               ' unsuccessful',&
           //,T2,'!!!!!!!    THE EXECUTION WAS STOPPED',&
                  '    !!!!!!!', &
            //,20('ERROR-'))
!
!
      RETURN
!
      END SUBROUTINE DellocateChemZoneArrays
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
