!
!
!
      SUBROUTINE ludcmp(a,       & ! Squared matrix
                        n,       & ! Actual matrix dimension
                        np,      & ! Maximum matrix dimension
                        indx,    & ! Working array
                        d)         !
!     &                  timetot,  ! Simulation time, s
!     &                  ielem,    ! Grid block number 
!     &                  iterch)   ! Chemistry iteration number

!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*         LU decomposition for the system of linear equations         *
!*                                                                     *
!*                  Version 1.0 - September 22, 2005                   *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE                                      
! 
! ----------
! ....Parameters
! ----------
! 
      INTEGER(KIND=8), PARAMETER :: nmax = 100
!
      REAL(KIND = 8), PARAMETER :: tiny = 1.0d-20
! 
! -------
! ... Integer variables
! -------
! 
      INTEGER(KIND=8)  :: np  
      INTEGER(KIND=8)  :: n, i, j, k, ielem, iterch, imax
! 
! -------
! ... Integer array
! -------
! 
      INTEGER(KIND=8), DIMENSION(np) :: indx
!
! -------
! ... Double precision variables
! -------
! 
      REAL(KIND=8) :: d, timetot, aamax, sum, dum
!
! -------
! ... Double precision arrays
! -------
! 

      REAL(KIND=8), DIMENSION(nmax)  :: vv         
      REAL(KIND=8), DIMENSION(np,np) :: a
!
!	 
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>
!
!
      d = 1.d0
!
      do i=1,n
!
         aamax = 0.d0

         do j=1,n
            if (dabs(a(i,j)) > aamax) aamax = dabs(a(i,j))
         end do
!
!........ Added to stop program and write to output file
!
         if (aamax == 0.d0)then
!
            write(32,*)'Singular Matrix in Chemical Solver, STOP'
            write(34,*)'Singular Matrix in Chemical Solver, STOP'
            print *, 'Singular Matrix in Chemical Solver, STOP'
!
            stop
!
         end if
!
         vv(i) = 1.d0/aamax
!
      end do
!
!
      DO_j: &
      do j=1,n
!
         if (j > 1) then
!
            do i=1,j-1
!
               sum = a(i,j)
!
               if (i > 1)then
!
                  do k=1,i-1
                     sum = sum - a(i,k)*a(k,j)
                  end do
!
                  a(i,j) = sum
!
               end if
!
            end do
!
         end if
!
         aamax = 0.d0
!
         do i=j,n
!
            sum=a(i,j)
!
            if (j > 1)then
!
               do k=1,j-1
                  sum=sum-a(i,k)*a(k,j)
               end do
!
               a(i,j) = sum
!
            end if
!
            dum = vv(i)*dabs(sum)

            if (dum >= aamax) then
               imax  = i
               aamax = dum
            end if
!
         end do
!
         if (j /= imax) then

            do k=1,n
               dum       = a(imax,k)
               a(imax,k) = a(j,k)
               a(j,k)    = dum
            end do
!
            d = -d
            vv(imax) = vv(j)
!
         end if
!
         indx(j) = imax
!
         if (j /= n) then
!
            if (a(j,j) == 0.d0) a(j,j) = tiny
!
            dum = 1.d0/a(j,j)
!
            do i=j+1,n
               a(i,j) = a(i,j)*dum
            end do
!
         end if
!
      end do  DO_j
!
!
      if (a(n,n) == 0.d0) a(n,n) = tiny
!
!
      return
!
      end  SUBROUTINE ludcmp
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE lubksb(a,n,np,indx,b)
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*         LU substitution for the system of linear equations          *
!*                                                                     *
!*                  Version 1.0 - September 22, 2005                   *     
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
      INTEGER(KIND=8)  :: np  
      INTEGER(KIND=8)  :: n, ii, i, j, ll
!
! -------
! ... Integer array
! -------
! 
      INTEGER(KIND=8), DIMENSION(np) :: indx
!
! -------
! ... Double precision variables
! -------
! 
      REAL(KIND=8) :: sum
!
! -------
! ... Double precision arrays
! -------
! 
      REAL(KIND=8), DIMENSION(np)    :: b         
      REAL(KIND=8), DIMENSION(np,np) :: a
!
!	 
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>
!
!
      ii = 0
!
      do i=1,n

         ll    = indx(i)
         sum   = b(ll)
         b(ll) = b(i)
!
         if (ii /= 0) then
!
            do j=ii,i-1
               sum = sum - a(i,j)*b(j)
            end do
!
         else if (sum /= 0.d0) then
!
            ii=i
!
         end if
!
         b(i) = sum
!
      end do
!
!
      do i=n,1,-1
!
         sum = b(i)
!
         if (i < n) then
!
            do j=i+1,n
               sum = sum - a(i,j)*b(j)
            end do
!
         end if
!
         b(i) = sum/a(i,i)
!
      end do
!
!
      return

      end  SUBROUTINE lubksb
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE assign
!
! ... Modules to be used 
! 
      USE Chem_Dimensions
      USE Batch_PhysiChem_Conditions
      USE Aqueous_Species
      USE Mineral_Phases
      USE Gaseous_Species
      USE Surface_Complexes
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*               Assigns log10K values to each reaction                *
!*               Also assigns temperature at current grid block        *
!*                                                                     *
!*                  Version 1.0 - September 26, 2005                   *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE                                      
! 
! -------
! ....Integer parameters
! -------
! 
      INTEGER (KIND=4), PARAMETER :: MaxD = 5
! 
! -------
! ... Integer variables
! -------
! 
      INTEGER(KIND = 8) :: j, k
      INTEGER(KIND = 8) :: icall = 0
      SAVE       icall 
!
! -------
! ... Double precision variables
! -------
! 
      REAL(KIND=8) :: LogK
!
! -------
! ... Double precision array
! -------
! 
      REAL(KIND=8), DIMENSION(MaxD) :: coef                      
!
!	 
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>
!
!
      icall = icall + 1
!
!
! -------
!.....Call the following subroutine to get parameters for DH activity coefficents
! -------
!
      CALL hkfpar(tc2)
!
!
! -------
!.....Get log10K for aqueous secondary species
! -------
!
      do j=1,naqx
!
         do k=1,5
            coef(k) = SSpe(j)%akcoes(k)
         end do
!
         SSpe(j)%aks = LogK(coef,tc2,MaxD)
!
      end do
!
!
! -------
!.....Get log10K for equlibrium minerals
! -------
!
      do j=1,nmequ
!
         do k=1,5
            coef(k) = MGen(j)%akcoem(k)
         end do
         MGen(j)%akm = LogK(coef,tc2,MaxD)
!
!........Add temperature dependent supersaturation window ssq0 (in log(K) units)
!..........exponential decrease from temp = sst1 to temp = sst2 (1/100 of initial value)
!
         MGen(j)%ssq = MGen(j)%ssq0
!
         if (tc2 > MGen(j)%sst1 .and. &
            MGen(j)%sst1 /= 0.d0 .and. MGen(j)%sst2 /= 0.d0)   then
!
            MGen(j)%ssq = MGen(j)%ssq0* &
                 dexp( -4.61d0/(MGen(j)%sst2 - MGen(j)%sst1) &
                 * (tc2 - MGen(j)%sst1) )        
         end if
!	          
      end do
!
!
! -------
!......Get log10K for kinetic minerals
! -------
!
      do j=1,nmkin
!
         do k=1,5
            coef(k) = MGen(nmequ+j)%akcoem(k)
         end do
         MGen(j)%akin = LogK(coef,tc2,MaxD)                              !!!! akin can be removed, use akm   !!! 
         MGen(nmequ+j)%akm = LogK(coef,tc2,MaxD) 
!
!........Add temperature dependent supersaturation window ssqk0 (in log(K) units)
!......... exponential decrease from temp = sstk1 to temp = sstk2 (1/100 of initial value)
!
         MKin(j)%ssqk = MKin(j)%ssqk0
!
         if (tc2 > MKin(j)%sstk1 .and. & 
            MKin(j)%sstk1 /= 0.d0 .and. MKin(j)%sstk2 /= 0.d0) then
!
            MKin(j)%ssqk = MKin(j)%ssqk0* &
              dexp( -4.61d0/(MKin(j)%sstk2 - MKin(j)%sstk1) &
              * (tc2 - MKin(j)%sstk1) )          
!
        end if
!
      end do
!
!
! -------
!......Get log10K for gaseous species
! -------
!
      do j=1,ngas
!
         do k=1,5
            coef(k) = GSpe(j)%akcoeg(k)
         end do
         GSpe(j)%akg = LogK(coef,tc2,MaxD)
!
      end do
!
      return
!
      end   SUBROUTINE assign
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      FUNCTION LogK(coef,temp,MaxDim)   
!
! ... Modules to be used 
! 
      USE Batch_PhysiChem_Conditions  
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*     Calculate log10K according to regression coefficients and T     *
!*                                                                     *
!*                  Version 1.0 - September 14, 2005                   *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE                                      
! 
! -------
! ... Integer variable
! -------
! 
      INTEGER(KIND=4)  :: MaxDim
!
! -------
! ... Real variables
! -------
! 
      REAL(KIND=8) :: temp, LogK, tk
!
! -------
! ... Real array
! -------
! 
      REAL(KIND=8), DIMENSION(MaxDim) :: coef                      
!
!	 
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>
!
!
      temp=max(temp,tmpmin)
      temp=min(temp,tmpmax)
!
      tk  =   temp + 273.15d0              
!
      LogK = coef(1)*dlog(tk) &
          + coef(2) &
          + coef(3)*tk & 
          + coef(4)/tk &
          + coef(5)/(tk**2)
!
      return
!
      end
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE hkfpar(T)
!
! ... Modules to be used 
! 
      USE Para_DH_ActivityCoefficients 
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*       Calculate parameters as a function of temperature             *
!*         to calculate activity coefficents with equations in         *
!*         Helgeson, Kirkham and Flowers, 1981, A.J.S. p.1249-1516     *
!*         (see subroutine dh_hkf81).                                  *
!*                                                                     *
!*                  Version 1.0 - September 23, 2005                   *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE                                      
! 
! -------
! ....Double precision parameters
! -------
! 
      REAL(KIND=8), PARAMETER :: aft1 =  0.49276542d+00,&
                                 aft2 =  0.31857945d-03,&
                                 aft3 =  0.11628933d-04,&
                                 aft4 = -0.52038832d-07,&
                                 aft5 =  0.12633045d-09
!
      REAL(KIND=8), PARAMETER :: bft1 =  0.32476341d+00,&
                                 bft2 =  0.12018502d-03,&
                                 bft3 =  0.79530646d-06,&
                                 bft4 = -0.30410531d-08,&
                                 bft5 =  0.56931304d-11
!
      REAL(KIND=8), PARAMETER :: bift1 =  0.26538636d+01,&
                                 bift2 = -0.52889569d-02,& 
                                 bift3 = -0.11009615d-03,& 
                                 bift4 =  0.46820513d-06,&
                                 bift5 = -0.11104895d-08 
!
      REAL(KIND=8), PARAMETER :: bilft1 = -0.14769091d+02,&
                                 bilft2 =  0.22563951d+00,&
                                 bilft3 = -0.10225385d-02,&
                                 bilft4 =  0.30349650d-05,&
                                 bilft5 = -0.35468531d-08
!
! -------
! ... Double precision variables
! -------
! 
      REAL(KIND=8)   :: T    ! Temperature in oC
      REAL(KIND=8)   :: bihat
!
!	 
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>
!
!   bhat NaCl (b NaCl = bhat/(2.303RT)) from Table 29 (bi here)
!   b Na+Cl- from Table 30 (bil here)
!   A and B Debye-Huckel (adh and bdh) from Table 1 (cols. 3 and 4)
!
!  Polynomial regression coefficients a, b, c, d, e, and f are
!  stored in data statements below to calculate parameters as:
!     parameter(T) = a + b*T + c*T + d*T**2 + e*T**3 + f*T**4
!  where T is temperature in C.
!
!  These are good for the range 0 to 300 C, BUT(!) the data
!  available to fit bi and bil did not go down to T = 0 C 
!  (first point at 25 C) so the extrapolation may not be too good
!  below 25 C for bi and bil (A and B data cover the entire range),
!  but I checked that the extrapolated bi and bil values below 
!  25 C vary smoothly down to 0 C.
!
!  Note the following units:
!    aft    yields adh in kg**0.5 mol**(-0.5)
!    bft    yields bdh in kg**0.5 mol**(-0.5) Anstrom**(-1)
!    bift   yields bihat(NaCl) in kg/mol * 1e+3
!    bilft  yields bil (Na+Cl-) in kg/mol * 1e+2
!
!	 
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>
!
!
      adh = func(aft1, aft2, aft3, aft4, aft5, T)
!
      bdh = func(bft1, bft2, bft3, bft4, bft5, T)
!
      bihat = func(bift1, bift2, bift3, bift4, bift5, T)*1.d-3
!
      bi = bihat/(2.303d0*1.987d0*(T+273.15d0))    ! [kg cal/mol]
!
      bil = func(bilft1, bilft2, bilft3, bilft4, bilft5, T)*1.d-2
!
!
      return
!
!
!***********************************************************************
!*                                                                     *
!*                        INTERNAL PROCEDURES                          *
!*                                                                     *
!***********************************************************************
!
!
      CONTAINS
! 
! ----------
! .......Double precision functions
! ----------
! 
         REAL(KIND = 8) FUNCTION func(a, b, c, d, e, T)
! 
            REAL(KIND = 8) :: a, b, c, d, e, T
!
            func = a + (b+(c+(d+e*T)*T)*T)*T
!
         END FUNCTION func
! 
!
      end  SUBROUTINE hkfpar
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Get_ActivityCoefficient(no_ch)
!
! ... Modules to be used 
! 
      USE MaxChem_Dimensions
      USE Chem_Dimensions
      USE Batch_PhysiChem_Conditions,    only : str                    
      USE CutOff_Threshold 
      USE Species_Index
      USE Aqueous_Species
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*      Computes activity coefficients of aqueous species and the      *
!*        activity of water using pitzer model or an extended          *
!*        Debye-Huckel (DH) model according to a user-specified        *
!*        ionic strength threshold and a user option                   *
!*                                                                     *
!*                  Version 1.0 - September 30, 2005                   *     
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
      INTEGER(KIND = 8) :: no_ch, i, j, n, m, ncp, nprsec, mopr9
      INTEGER(KIND = 8) :: icall = 0
      SAVE       icall 
!
! -------
! ... Real variables
! -------
! 
      REAL(KIND=8) :: mstar, mchr, sum, sum2
!      common/str_thres/str_threshold  ! Ionic strength threshold for switch between pitzer and DH
      REAL(KIND=8) :: str2, stroot, stroo2, capgam, xh2rmh2, csjxh
!
! -------
! ... Real array
! -------
! 
      REAL(KIND=8), DIMENSION(MaxNaqt) :: cpion
!
!	 
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>
!
!
      icall = icall + 1
!     
!
      no_ch = 0
!
! ----------
!.....Calculate ionic strength and other global concentrations.
!.... All concentrations must be in molal scale (mol/kgH2O) 
! ----------
!
      sum   = 0.d0             ! Compute true ionic strength
      sum2  = 0.d0             ! Compute stoichiometric ionic str
      mstar = 0.d0             ! Total solute in solution
      mchr  = 0.d0             ! Total solute excluding neutral species
!
!.....Primary species contribution 
!
      do i=1,npri
         ASpe(i)%ct = PSpe(i)%cp

         if (      i /= nw  &     ! Not H2O
            .and. i /= ne   &     ! Not e-
            .and. i /= nd) then   ! Not primary surface species
!
            sum   = sum   + ASpe(i)%zz2*PSpe(i)%cp
            sum2  = sum2  + ASpe(i)%zz2*PSpe(i)%u2/xh2o
            mstar = mstar + PSpe(i)%cp
!
            if (ASpe(i)%z /= 0.d0) then 
		     mchr = mchr + PSpe(i)%cp
            end if
!
            cpion(i) = PSpe(i)%cp
!		 
         end if
!
      end do
!
!
!.....Secondary species contribution
!
      do i=1,naqx
!
         nprsec = npri + i
         ASpe(nprsec)%ct = SSpe(i)%cs
!
         sum   = sum   + ASpe(nprsec)%zz2 *SSpe(i)%cs
         mstar = mstar + SSpe(i)%cs
! 
         if (ASpe(nprsec)%z /= 0.d0) then
!
            mchr = mchr + SSpe(i)%cs
!
!...........Total concentrations excluding neutral species
!
            ncp = SSpe(i)%ncps
!
            do n=1,ncp
               j = SSpe(i)%icps(n)
               cpion(j) = cpion(j) + SSpe(i)%stqs(n)*SSpe(i)%cs
            end do
!
         end if
!	    
      end do
!
!
      if (sum < 0.d0)  sum = 0.d0
      str = 0.5d0*sum             ! Ionic strength
!
!
      if (str > stimax) then
         no_ch = 1
         return
      end if  
!
!
      str2 = str
!
      if (sum2 > 0.d0) str2 = 0.5d0*sum2
      stroot = dsqrt(str)
      stroo2 = dsqrt(str2)
!
!
!.....Conversion factor for molality scale (from eq. 122, 169)
!.....  Note, this factor turns out equal to log(H2O mole fraction)
!.....  D-H yields gamma for mole fraction scale convention, we
!       use molality scale convention, therefore needs conversion 
!
      capgam = -dlog10(1.d0 + 0.01801528d0*mstar)
!
!
         CALL Extended_DH_hkf81 &
          (str2, stroot, stroo2, capgam, cpion, mstar, mchr)
!
      return
      
!
      end  SUBROUTINE Get_ActivityCoefficient
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Extended_DH_hkf81 &
                 (str2, stroot, stroo2, capgam, cpion, mstar, mchr)
!
! ... Modules to be used 
! 
      USE MaxChem_Dimensions
      USE Chem_Dimensions
      USE Batch_PhysiChem_Conditions
      USE Species_Index
      USE Aqueous_Species
      USE Para_DH_ActivityCoefficients 
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*      Computes activity coefficients of aqueous species and the      *
!*        activity of water using the extended DH                      *   
!*                                                                     *
!*                  Version 1.0 - September 30, 2005                   *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
!
!  Flowers, 1981, A.J.S. p.1249-1516. Also computes neutral species activities 
!  from Setchenov equation (Langmuir 1997, Aqueous Environmental 
!  Geochemistry, Prentice Hall, p. 144). 
!  Charged species: individual ion activity coefficients calculated from:
!     equation 298 (which includes 121, 122, 129, 130, 169, 170 and 297) 
!      - uses true ionic strenght. 
!     use bhat NaCl (b NaCl = bhat/(2.303RT)) from Table 29 (bi here)
!     use b Na+Cl- from Table 30 (bil here)
!     use Rej from Table 3 (input in thermodynamic database as 
!       a0 variable, but note that values are NOT a0 values)
!     use A and B Debye-Huckel (adh and bdh) from Table 1 (cols. 3 and 4)
!     caclulate a0 from input Rej values and eq. 125,
!       assuming other dominant anion (for cations) is Cl- (rej=1.81 A)
!       and cation (for anions) is Na+ (rej=1.91 A)
!
!  Activity of water calculated from:
!     osmotic coefficient using equation 190 and same parameters
!       as above, assuming similar simplifications as done for 
!       the calculation of activity coefficients, but using
!       stoichiometric ionic strength. 
!     equation 106 relating activity of water to the osm. coef.
!
!  Regression coefficients a, b, c, d, e to obtain A, B, bi abd bil
!  parameter as a function of temperature were calculated using
!  4th order polynomials as follow:
!        f(T) = a + b*T + c*T**2 + d*T**3 + e*T**4 
!
!  Neutral species: assume gamma = 1 or,  
!  for weak acids and dissolved gases:
!    log(activity coef)=sltout*ionic strength  
!    where sltout is input in a0 variable as 100+sltout 
!    For now, no temperature dependence on sltout is assumed as
!    the effect is small compared to variation of solubility (K) with temp. 
!   
!
!***********************************************************************
!***********************************************************************
!
!
      IMPLICIT NONE                                      
! 
! -------
! ... Integer variables
! -------
! 
      INTEGER(KIND = 8) :: i, j, npripj
      INTEGER(KIND = 8) :: icall = 0
      SAVE       icall 
!
! -------
! ... Real variables
! -------
! 
      REAL(KIND=8) :: mstar, mchr, sum, sum2
      REAL(KIND=8) :: summt, rex, azero, zabsi, lambda, lambd2 
      REAL(KIND=8) :: omega, gamlog, sltout, osmo, gsalt, gamln
      REAL(KIND=8) :: str2, stroot, stroo2, capgam
      REAL(KIND=8) :: CC, FF, GG, EE, HH
!
! -------
! ... Real arrays
! -------
! 
      REAL(KIND=8), DIMENSION(MaxNpri) :: dstr
      REAL(KIND=8), DIMENSION(MaxNaqt) :: cpion
!
!	 
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>
!
!
      summt = 0.d0 
!
      do i=1,naqt
!
         if (i <= npri) then 
            PSpe(i)%gamp      = 1.0d0
         else
            SSpe(i-npri)%gams = 1.d0
         end if
!
! ----------
!........The following hkf parameters at given temperature are calculated 
!........in SUBROUTINE hkfpar(T) and called by SUBROUTINE assign:
!............. adh in kg**0.5 mol**(-0.5)  (this is A Debye-Huckel)
!..............bdh in kg**0.5 mol**(-0.5) Anstrom**(-1) (this is B Debye-Huckel)
!..............bi  in kg/mol * 1e+3 (this is bhat (NaCl) )
!..............bil in kg/mol * 1e+2 (this is bil (Na+ Cl-) )
! ----------
!
         rex = 1.81d0
         if (ASpe(i)%z < 0.d0) rex = 1.91d0
         zabsi = dabs(ASpe(i)%z)
         azero = 2.d0*(ASpe(i)%a0 + zabsi*rex)/(zabsi + 1.d0)
         lambda = 1.d0 + bdh*azero*stroot
         lambd2 = 1.d0 + bdh*azero*stroo2
!
! ----------
!........Charged species
! ----------
!
         if (ASpe(i)%z /= 0.d0) then
!
            omega  = 1.66027d05*ASpe(i)%zz2/ASpe(i)%a0
!
            gamlog = - adh*ASpe(i)%zz2*stroot/lambda  &
                     + capgam  &
                     + (omega*bi*str + (bil - 0.19d0*(zabsi-1.d0))*str)
!
!...........Primary charged species
!
            if (i <= npri) then
!
               PSpe(i)%gamp = 10.d0**gamlog     ! Activity coef. of primary species
!
               summt = summt + cpion(i) * &
                      (adh*ASpe(i)%zz2 &
                      / (azero*azero*azero*bdh*bdh*bdh*str2)&
                      * (lambd2- 1.d0/lambd2 - 2.d0*dlog(lambd2))&
                     + capgam/(0.0180153d0*mstar)&
                      - 0.5d0*(omega*bi*str2 + &
                      (bil-0.19d0*(zabsi-1.d0))*mchr*0.5d0))
!
!...........Secondary charged species
!
            else
               SSpe(i-npri)%gams = 10.d0**gamlog    !activity coef. of secondary species
            end if
!
! ----------
!........Neutral species
! ----------
!
         else
!
            if (ASpe(i)%a0 > 100.d0) then
!
               sltout = ASpe(i)%a0 - 100.d0     ! Read in as sltout + 100 (flag)
               gamlog = sltout*str
!
               if (i <= npri .and. i /= nw) then
                  PSpe(i)%gamp      = 10.d0**gamlog   
               else
                  SSpe(i-npri)%gams = 10.d0**gamlog  
               end if
!
            end if
!
        end if
!
      end do
!
!
      osmo = -2.303d0*summt/mstar
      PSpe(nw)%gamp = dexp(-osmo*mstar/55.50837d0)    ! Activity coef. of water
      if(ne > 0) PSpe(ne)%gamp = 1.d0            ! Activity coef. of electron
      if(nd > 0) PSpe(nd)%gamp = 1.d0            ! Activity coef. of surface species
!
!
! ----------
!.....Salting out effect for gaseous speices, Wolery's EQ3NR Eq. (91)
! ----------
!
!      if (str > 0.3d0) then
      if (str > 0.0d0) then
!
!
         cc = -1.0312d0
         ff =  0.0012806d0
         gg =  255.9d0
         ee =  0.4445d0
         hh = -0.001606d0
!
         gamln =   (cc+ff*tk2+gg/tk2)*str &
                 - (ee+hh*tk2)*(str/(str+1.0d0))
         gsalt = dexp(gamln)

         do i=1,npri
!
            if (    PSpe(i)%napri == 'co2(aq)' & 
              .or. PSpe(i)%napri == 'CO2(aq)' &
!
              .or. PSpe(i)%napri == 'ch4(aq)' &
              .or. PSpe(i)%napri == 'CH4(aq)' &
!     
              .or. PSpe(i)%napri == 'h2(aq)'  &
              .or. PSpe(i)%napri == 'H2(aq)'  &
!
              .or. PSpe(i)%napri == 'h2s(aq)' &
              .or. PSpe(i)%napri == 'H2S(aq)' &
!
              .or. PSpe(i)%napri == 'o2(aq)'  &
              .or. PSpe(i)%napri == 'O2(aq)'  &
!
              .or. PSpe(i)%napri == 'so2(aq)' &
              .or. PSpe(i)%napri == 'SO2(aq)')   then 
!
                     PSpe(i)%gamp = gsalt
!              
            end if
!
         end do
!
         do i=1,naqx
            if (    SSpe(i)%naaqx == 'co2(aq)'& 
              .or. SSpe(i)%naaqx == 'CO2(aq)' &
!      
              .or. SSpe(i)%naaqx == 'ch4(aq)' &
              .or. SSpe(i)%naaqx == 'CH4(aq)' &
!
              .or. SSpe(i)%naaqx == 'h2(aq)'  &
              .or. SSpe(i)%naaqx == 'H2(aq)'  &
!
              .or. SSpe(i)%naaqx == 'h2s(aq)' &
              .or. SSpe(i)%naaqx == 'H2S(aq)' &
!
              .or. SSpe(i)%naaqx == 'o2(aq)'  &
              .or. SSpe(i)%naaqx == 'O2(aq)'  &
!
              .or. SSpe(i)%naaqx == 'so2(aq)' &
              .or. SSpe(i)%naaqx == 'SO2(aq)')   then
!
                     SSpe(i)%gams = gsalt
!
            end if
!
         end do
!
      end if
!
!
      return
!
!
      end  SUBROUTINE Extended_DH_hkf81
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE cs_cp
!
! ... Modules to be used 
! 
      USE Chem_Dimensions
      USE Chem_Constants
      USE CutOff_Threshold 
      USE Species_Index
      USE Aqueous_Species
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*      Calculates total solute in solution and partial derivatives    *
!*        with respect to component species (mass balance in moles)    *
!*        derivative terms of gamma are included for using the Pitzer  *
!*        ionic activity model                                         *   
!*                                                                     *
!*                  Version 1.0 - October 05, 2005                     *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
!
      IMPLICIT NONE                                      
! 
! -------
! ... Integer variables
! -------
! 
      INTEGER(KIND = 8) :: i, j, k, l, ncp, m, n
      INTEGER(KIND = 8) :: icall = 0
      SAVE       icall 
!
! -------
! ... Double precision variables
! -------
! 
      REAL(KIND=8) :: dlstmx, xh2rmh2, cpixh2o, csjxh, cplxh
!
!	 
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>
!
!
      dlstmx = dlog10(stimax)      ! Ionic strength cutoff
!
! ----------
!.....Calculate the concentrations of secondary (derived) species.
!.......cp is concentration of primary species, of secondary species
!...   (in moles per kg water)
! ----------
!
       icall=icall+1
      do j=1,naqx
!
         SSpe(j)%cs =           &   ! SSPe = secondary species 
!                                  ! cs  will be concentration of secondary species (mol/kg H2O) 	    
          - dlog10(SSpe(j)%gams)&  ! gams = activity coef. of secondary species 
          - SSpe(j)%aks           ! aks  = log10K of secondary species dissociation reaction
!
         ncp = SSpe(j)%ncps        ! Number of primary species involved
!
         do n=1,ncp
!
            i = SSpe(j)%icps(n)    ! Get primary species index
!
            SSpe(j)%cs =  SSpe(j)%cs &
                        + SSpe(j)%stqs(n) &      ! stqs = stoichometric coefficient
                         *dlog10(PSpe(i)%cp &    ! cp   = concentration of primary species (mol/kg)
                                *PSpe(i)%gamp) ! gamp = activity coefficient of primary species
! 
         end do
!
!........Add if/else cut-off to avoid convergence problems in rare cases
!..........note: cs(j) above is first calculated in log form then converted below
!
         if (SSpe(j)%cs > dlstmx) then                    !!!!  ?????   should stop simulation
            SSpe(j)%cs = SSpe(j)%cs/10.d0                         !!!!  ?????
         else                      
            SSpe(j)%cs = 10.d0**SSpe(j)%cs     
         end if
!
      end do                       
!
!
! ----------
!     Total solute in solution is u2 = total moles (primary + secondary species).
!       First get u2 contributions from primary species and their derivatives
!       temp variables to limit calcs inside loop
! ----------
!
      xh2rmh2 = xh2o*rmh2o
!
      do i=1,npri
!
         PSpe(i)%u2 = PSpe(i)%cp*xh2o 
!
         do k=1,npri
            PSpe(i)%du2(k) = 0.0d0   
         end do 
!
         PSpe(i)%du2(i) = PSpe(i)%u2  ! Diagonal terms, using relative increment                   
!
      end do
!
      PSpe(nw)%u2      = xh2rmh2      ! For water
      PSpe(nw)%du2(nw) = xh2rmh2      ! Diagonal term, also use relative increment 
!
!
      DO_j: &
      do j=1,naqx
!
         ncp   = SSpe(j)%ncps
         csjxh = SSpe(j)%cs*xh2o
!
         DO_n: &
         do n=1,ncp
!
            i = SSpe(j)%icps(n)
            PSpe(i)%u2 = PSpe(i)%u2 + SSpe(j)%stqs(n)*csjxh
            cplxh = 0.0d0
!
            DO_l: &
            do l=1,ncp
!
               k = SSpe(j)%icps(l)
!
               if (k /= nw)   then
                  cplxh = cplxh + &
                          SSpe(j)%stqs(l)*PSpe(k)%dgamp(i)/PSpe(i)%gamp 
               else
                  cplxh = cplxh + &
                          SSpe(j)%stqs(l)*PSpe(k)%dgamp(i)
               end if
!
            end do  DO_l
!
            DO_m: &
            do m=1,ncp
!
               k = SSpe(j)%icps(m)
               cpixh2o = PSpe(k)%cp*xh2o
!
               if (k /= nw)   then
!
                  PSpe(i)%du2(k) = PSpe(i)%du2(k) + &
                     SSpe(j)%stqs(m)*(csjxh*SSpe(j)%stqs(n) &
                     - cpixh2o*SSpe(j)%cs*SSpe(j)%dgams(i)/SSpe(j)%gams &    ! Derivative of gam
                     + SSpe(j)%cs*cpixh2o*cplxh)
!
               else
!
                  PSpe(i)%du2(k) = PSpe(i)%du2(k) + &
!
                    csjxh*(SSpe(j)%stqs(n)&
!
                    - SSpe(j)%dgams(i)/SSpe(j)%gams &  ! Derivative of gam
!
                    + cplxh)
!
               end if
!
            end do  DO_m
!
         end do  DO_n
!
      end do  DO_j
!
!

      return
!
      end  SUBROUTINE cs_cp
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!

      SUBROUTINE MineralSaturationIndex
! 
! ... Modules to be used 
! 
      USE Chem_Dimensions
      USE Aqueous_Species
      USE Mineral_Phases
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*    This routine computes saturation indexes of all minerals         *
!*                                                                     *
!*                  Version 1.0 - July 20, 2005                        *     
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
      INTEGER(KIND = 8) :: m, n, i, ncp 
      INTEGER(KIND = 8) :: icall = 0
      SAVE       icall 
! 
! -------
! ... Double precision variable
! -------
! 
      REAL(KIND = 8) :: paim
! 
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of SUBROUTINE ineralSaturationIndex
!
!
      icall = icall + 1
!
!
       do m=1,nmin
!
         paim = 0.0d0
         ncp = MGen(m)%ncpm  ! ! Number of primary species involved in the reaction
!
         do n=1,ncp
!
            i = MGen(m)%icpm(n)
!
            paim =   paim &
                  + MGen(m)%stqm(n)  &       ! Stoichiometric coefficient of primary species in the mineral
                    *(dlog10(PSpe(i)%cp) &   ! Concentration of primary species
                  + dlog10(PSpe(i)%gamp))    ! Activity coefficient of primary species
!
         end do
!       
         MGen(m)%si2_old = MGen(m)%si2        ! Save the old Q/K                !!!!!
         MGen(m)%si2 = paim - MGen(m)%akm     ! log10 Q/K
      end do
!
!
        return
!
        end  SUBROUTINE MineralSaturationIndex
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE cmgas_cp(iinit)
!
!                                          
! ... Modules to be used 
! 
      USE Chem_Dimensions
      USE ChemOption_Variables,   only : ngas1
      USE Aqueous_Species
      USE Gaseous_Species
      USE Batch_Chem_Jacobian
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*    This routine computes saturation indexes of gases at             *
!*         equilibrium and fugacities of gases at equilibrium.         *     
!*                                                                     *
!*                  Version 1.0 - July 20, 2005                        *     
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
      INTEGER(KIND = 8) :: i, j, n, ncp, iinit
      INTEGER(KIND = 8) :: icall = 0
      SAVE       icall
! 
! -------
! ... Double precision variable
! -------
!
! 
      REAL(KIND = 8) :: paig
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>>=>  Main body of SUBROUTINE cmgas_cp
!
!
      icall = icall + 1
!
!
      do j=1,ngas
!
         paig = 0.0d0
         ncp  = GSpe(j)% ncpg   ! Number of primary species involved in the reaction
!
         do n=1,ncp
            i = GSpe(j)%icpg(n) ! Index of primary species in reaction stoichiometry
            paig =  paig  &
                 + GSpe(j)%stqg(n)   &     ! Stoichiometric coefficient of primary species in the gas
                   *(dlog10(PSpe(i)%cp) &  ! Concentration of primary species
                 + dlog10(PSpe(i)%gamp))  ! Activity coefficient of primary species
         end do
!
         if (nsatg > 0)          then  
!
! ----------
! ..........Saturation index of gases 
! ---------- 
!         
!           Call subroutine to obtain gas fugacity coefficient, gamg subroutine        
!
            CALL GasFugacityCoefficients
!
!...................sig2 = log10 Q/K/fgas
            GSpe(j)%sig2 = paig &
                        - dlog10(GSpe(j)%cg) &   ! Gas partial pressure, bar
                        - dlog10(GSpe(j)%gamg)  &! Gas fugacity coefficient
                        - GSpe(j)%akg           ! log10K  
!
                                 else                      
!
! ----------
! ..........Gas fugacity for undersaturated gas.
! ..........Note there are cases where nsatg=0 is initially 0 but ngas1=0
!...........iinit=1 for initialization (calls from inchm.f), iinit=0 for regular calls
! ---------- 
!
            CALL GasFugacityCoefficients
!
            if (ngas1 > 0  .or.  iinit == 1)   then
                GSpe(j)%cg = paig - GSpe(j)%akg ! log f = log(q/k) always, even if no gas exsolved
                GSpe(j)%cg = 10.d0**GSpe(j)%cg  ! Gas partial pressure, bar
                GSpe(j)%cg = GSpe(j)%cg / GSpe(j)%gamg ! Gas partial pressure, bar
            end if
!
         end if
!
      end do
!
      return
!
      end  SUBROUTINE cmgas_cp
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE GasFugacityCoefficients
! 
! ... Modules to be used 
! 
      USE Chem_Dimensions
      USE Batch_PhysiChem_Conditions
      USE Species_Index
      USE Gaseous_Species
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*        This routine computes gas fugucity coefficient               *
!*                                                                     *
!*                  Version 1.0 - July 20, 2005                        *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE                                      
! 
! 
! -------
! ... Integer variables
! -------
! 
      INTEGER(KIND = 8) :: ig       
      INTEGER(KIND = 8) :: icall = 0
      SAVE       icall
!
! -------
! ... Double precision variable
! -------
! 
      REAL(KIND = 8) :: paig
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>>=>  Main body of SUBROUTINE cmgas_cp
!
!
      icall = icall + 1
!
!
      do ig=1,ngas
!
         GSpe(ig)%gamg = 1.0d0         ! Default gas fugacity coefficients
!
!         if (ig == nco2g) then
!
!...........For ! For EOS2, ECO2 and ECO2N modules
!
!c            if (ico2m == 1)         then   
!c!
!c               if (ieos == 14)   then
!c                  GSpe(ig)%gamg = fugcoeCO22
!c	                           else
!c                  CALL fuga_coe_CO2(Pt2,Tk2,ig)
!c               end if
!c!
!c            end if
!c!
!c         end if
!c!
      end do
!
!
      return

      end SUBROUTINE GasFugacityCoefficients
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE fuga_coe_CO2(Pt,Ta,ig)
!
! ... Modules to be used 
! 
      USE Gaseous_Species
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*            Calculate CO2 gas fugacity coefficients                  *     
!*                                                                     *
!*                  Version 1.0 - October 5, 2005                      *     
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
      INTEGER(KIND = 8) :: ig 
      INTEGER(KIND = 8) :: icall = 0
      SAVE       icall 
! 
! -------
! ....Double precision parameters
! -------
!  
      REAL(KIND = 8), PARAMETER ::  aaa = -1430.87d+0
      REAL(KIND = 8), PARAMETER ::  bbb =  3.598d+0
      REAL(KIND = 8), PARAMETER ::  ccc = -227.376d-5
      REAL(KIND = 8), PARAMETER ::  ddd =  347.644d-2
      REAL(KIND = 8), PARAMETER ::  eee = -1042.47d-5
      REAL(KIND = 8), PARAMETER ::  fff =  846.271d-8
! 
! -------
! ... Double precision variables
! -------
! 
      REAL(KIND = 8) :: Pt, Ta, Fefec
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>>=>=>=>=>=>=>=>=>=>=>=>=>=>=>>=>
!
!
      icall = icall + 1
!
!
!.....Assume real gases and ideal mixture
!........CO2 fugacity coefficient is calculated from eq.(14) 
!........of the paper of Nicolas Spycher et al., 1988
!........For 50-350 0C and 0-500 bar
!
      GSpe(ig)%gamg =   (aaa/(Ta*Ta) + bbb/Ta + ccc)*Pt
!
      GSpe(ig)%gamg =   GSpe(ig)%gamg + (ddd/(Ta*Ta) &
                       + eee/Ta + fff)*Pt*Pt/2.0d0
!
      GSpe(ig)%gamg =   dexp(GSpe(ig)%gamg)
!
      if (GSpe(ig)%gamg <= 0.0D0   .or. &
          GSpe(ig)%gamg >  1.0D0)  then
!
         write (32,*) 'There is a problem with fugacity coefficient'
         print *, 'There is a problem with fugacity coefficient'
!
         stop

      end if
!
      Fefec = 1.00d0        ! Effective factor
      GSpe(ig)%gamg = Fefec*GSpe(ig)%gamg
!
!
      return
!
      end  SUBROUTINE fuga_coe_CO2
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
      SUBROUTINE cx_ct
! 
! ... Modules to be used 
! 
      USE MaxChem_Dimensions
      USE Chem_Dimensions
      USE Species_Index
      USE Batch_PhysiChem_Conditions,         only : phi2, phisl2
      USE Aqueous_Species
      USE Ion_Exchange
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*          Calculate concentrations of exchanged cations              *     
!*                                                                     *
!*                   Version 1.0 - July 05, 2006                       *     
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
      INTEGER(KIND = 8) :: i, j, n_x, n_j
      INTEGER(KIND = 8) :: icall = 0
      SAVE       icall 
! 
! -------
! ... Real variables
! -------
! 
      REAL(KIND = 8) :: z_x, z_j, a1, a2, a3, z1, z2, z3 
      REAL(KIND = 8) :: p0, p1, p2, cecmol
! 
! -------
! ... Real arrays
! -------
! 
      REAL(KIND=8), DIMENSION(MaxNexc) :: dum, bx, zr
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>>=>=>=>=>=>=>=>=>=>>=>=>=>=>
!
!
      icall = icall + 1
!
!
! ----------
!.....The terms of addition of total eq. fraction as a funtion of bx(nx)
! ----------
!
      a1 = 0.d0
      a2 = 0.d0
      a3 = 0.d0
!
      n_x = XSpe(nx)%nbx    ! Primary exchanged species index in the aqueous primary list
      z_x = ASpe(n_x)%z     ! The electric charge of the primary exchanged species
!
! ----------
!.....Gaines-Thomas and vanselow conventions
! ----------
!
!
      IF_GainesThomas_or_Vanselow: &
      if (iex == 1 .or. iex == 2) then
!
!.........iex     is convention index: 1 = Gaines-Thomas, 2 = Vanselow; 3 = Gapon 
!
         do j=1,nexc
!
            n_j = XSpe(j)%nbx     ! Index in the list of priamry species
            z_j = ASpe(n_j)%z     ! Electric charge
!
            zr(j)  = z_j/z_x     
!
            dum(j) =  XSpe(j)%ekx**(-z_j)&
                    * PSpe(n_j)%cp*PSpe(n_j)%gamp&
                    *(PSpe(n_x)%cp*PSpe(n_x)%gamp)&
                    **(-zr(j))
!
            if (zr(j) == 1.d0)  a1 = a1 + dum(j)
!
            if (zr(j) == 2.d0)  a2 = a2 + dum(j)
!
            if (zr(j) == 3.d0)  a3 = a3 + dum(j)
!
         end do
!
!........Resolution of: a1*bx(nx) + a2*bx(nx)**2 + a3*bx(nx)**3 - 1.0 = 0
!
         IF_a3EQ0:&
         if (a3 == 0.d0)                           then
!
            if (a2 == 0.d0)   then
               if (a1 /= 0.d0)  bx(nx) = 1.d0/a1
                              else
               bx(nx) = (-a1 + dsqrt(a1*a1 + 4.d0*a2))*0.5d0/a2
            end if
!
                                                   else
!
            IF_a1EQ0_a2EQ0:&
            if (a1 == 0.d0 .and. a2 == 0.d0) then
               bx(nx) = (1.d0/a3)**(1.d0/3.d0)
!
                                           else
!
               p1 =  a1/a3
               p2 =  a2/a3
               p0 = -1.d0/a3
!
               CALL cubic (p2,p1,p0,z1,z2,z3)
!
               if (z1 > 0.d0 .and. z1 < 1.d0 .and. &
                   z2 > 0.d0 .and. z2 < 1.d0 .and. &
                   dabs(z1-z2)/z1 > 1.d-3)            go to 2000
!
               if (z1 > 0.d0 .and. z1 < 1.d0 .and. &
                   z3 > 0.d0 .and. z3 < 1.d0 .and. &
                   dabs(z1-z3)/z1 > 1.d-3)            go to 2000
!
               if (z3 > 0.d0 .and. z3 < 1.d0 .and. &
                   z2 > 0.d0 .and. z2 < 1.d0 .and. &
                   dabs(z3-z2)/z3 > 1.d-3)            go to 2000
!
               bx(nx) = z1
!
               if (z2 > 0.d0 .and. z2 < 1.d0)  bx(nx) = z2
!
               if (z3 > 0.d0 .and. z3 < 1.d0)  bx(nx) = z3
!
            end if  IF_a1EQ0_a2EQ0
!
         end if  IF_a3EQ0
!
      end if  IF_GainesThomas_or_Vanselow
!
!
! ----------
!.....Gapon convention
! ----------
!
      IF_GaponConvention:&
      if (iex == 3) then
!
         do j=1,nexc
!
            n_j = XSpe(j)%nbx     ! Index in the list of priamry species
            z_j = ASpe(n_j)%z     ! Electric charge
!
            dum(j) = (1.d0/XSpe(j)%ekx)&
                   *(PSpe(n_j)%cp*PSpe(n_j)%gamp)**(1.d0/z_j)&
                   *(PSpe(n_x)%cp*PSpe(n_x)%gamp)**(-1.d0/z_x)
!
            a1 = a1 + dum(j)
!
         end do
!
         bx(nx) = 1.d0/a1
!
      end if  IF_GaponConvention
!
!
      if (bx(nx) > 1.d0 .or. bx(nx) < 0.d0)    go to 2000
!
! ----------
!.....The rest of eq. fractions as function of bx(nx)
! ----------
!
      do j=1,nexc
!
         if (iex == 1 .or. iex == 2)   then
!
            bx(j) = dum(j)*bx(nx)**(zr(j))
!
                                       else
!
            bx(j) = dum(j)*bx(nx)
!
         end if
!
      end do
!
!
! ----------
!.....Conversion of bx (eq. fraction) into cx (mol solute ads/dm3 sol)
! ----------
!
      cecmol = cec2*2.65d0*(1.d0-phi2)*1.d-2/phisl2
!
      do j=1,nexc
!
!........Gaines-thomas and Gapon conventions
!
         if (iex == 1 .or. iex == 3)   then 
            n_j = XSpe(j)%nbx     ! Index in the list of priamry species
            z_j = ASpe(n_j)%z     ! Electric charge
	      XSpe(j)%cx = bx(j)*cecmol/z_j
         end if
!
!........Vanselow convention
!
         if (iex == 2)  XSpe(j)%cx = bx(j)*cecmol
!
      end do
!
!
      return
!
!
2000  continue
!
      write (*,2100)  z1, z2, z3
      write (32,2100) z1, z2, z3
2100  format(///,1x,3e14.3,'error or ambiguity in ion exchange calc.')
!
      stop
!
!
      end  SUBROUTINE cx_ct
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
      SUBROUTINE cubic (a,b,c,z1,z2,z3)
! 
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*     Calculate real roots of a cubic equation (newton-raphson)       *     
!*                                                                     *
!*                   Version 1.0 - July 05, 2006                       *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE                                      
! 
! -------
! ... Double precision variables
! -------
! 
      REAL(KIND = 8) :: a, b, c, z, z1, z2, z3, f, h 
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>>=>=>=>=>=>=>=>=>=>>=>=>=>=>
!
!       
      z = 1.d0  ! Trial value
!
40    f = z*z*z + a*z*z + b*z + c
!
      h = f/(3.d0*z*z + 2.d0*a*z + b)
!
      if (dabs(h/z) <= 1.0d-4)   go to 90
!
      z = z - h
!
      go to 40
!
90    z1 = z
!
      z  = 0.01d0  ! Second root trial
!
110   f = z*z*z + a*z*z + b*z + c
!
      h = f/(3.d0*z*z + 2.d0*a*z + b)
!
      if (dabs(h/z) <= 1.0d-4)   go to 160
!
      z = z - h
!
      go to 110
!
160   z2 = z
!
170   z3 = -a - z1 - z2
!
!
      return
!
      end  SUBROUTINE cubic
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE dcx_dcp
! 
! ... Modules to be used 
! 
      USE MaxChem_Dimensions
      USE Chem_Dimensions
      USE Aqueous_Species
      USE Ion_Exchange
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*          Evaluate derivatives for cation exchange                   *     
!*                                                                     *
!*                   Version 1.0 - July 05, 2006                       *     
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
      INTEGER(KIND = 8) :: i, j
      INTEGER(KIND = 8) :: icall = 0
      SAVE       icall 
! 
! -------
! ... Real variable
! -------
! 
      REAL(KIND = 8) :: dd
!
! -------
! ... Real array
! -------
! 
      REAL(KIND=8), DIMENSION(MaxNexc) :: cxold
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>>=>=>=>=>=>=>=>=>=>>=>=>=>=>
!
!
      icall = icall + 1
!      
! ----------
!.....Numerical derivative of cx respecte a cp
! ----------
!       
      do j=1,nexc
         cxold(j) = XSpe(j)%cx
      end do
!
      do i=1,npri
!
         dd = PSpe(i)%cp*1.0d-07
         PSpe(i)%cp = PSpe(i)%cp + dd
!
         CALL cx_ct
!
         do j=1,nexc
            XSpe(j)%dcx(i) = (XSpe(j)%cx - cxold(j))/dd
            XSpe(j)%cx     = cxold(j)
         end do
!
         PSpe(i)%cp = PSpe(i)%cp - dd
!
      end do
!
!
      return
!
      end  SUBROUTINE dcx_dcp
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE ADSORPTION(IT)
! 
! ... Modules to be used 
! 
      USE Chem_Dimensions
      USE Chem_Constants
      USE Species_Index
      USE Batch_PhysiChem_Conditions
      USE Convergence_Variables
      USE Aqueous_Species
      USE Surface_Complexes
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*          Solve adsorption by Newton-Rapson iteration method         *     
!*                                                                     *
!*                  Version 1.0 - October 12, 2005                     *     
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
      INTEGER(KIND = 8) :: k, it, nND
      INTEGER(KIND = 8) :: icall = 0
      SAVE       icall 
! 
! -------
! ... Real variables
! -------
! 
      REAL(KIND = 8) :: CAPPAINV, S, SUMY, SUMZY, SUMNUY, DPHIP2, &
                        FNA1, FNA2, DS, S2, SUMZHZY, SUMZNUY, DSRE
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>>=>=>=>=>=>=>=>=>=>>=>=>=>=>
!
!
!.....Double layer thickness (dm)
!
      CAPPAINV = 3.05d0*1.0D-09/DSQRT(STR)
!
!.....Calculate alpha term
!
      ADFACTOR = CAPPAINV*FARADAY/(EPSI*SUPADN2*RGAS*tk2)     !  ALFA
!
!.....Start iteration scheme for Delta(PHI) and Delta(S)
!
      S = PSpe(nd)%CP
!
      ITERAD = 0
!
299   ITERAD = ITERAD + 1
!
!.....Call CD_CP to obtain the amount of surface complexes and potencial term
!
      CALL  CD_CP
!
      IF (ITERAD == 1)   DPHIP2 = PHIP2NEW - PHIP2
!
!.....Evaluate the terms for calculating new Delta(S)
!
      SUMY = 0.d0
!
      DO K=1,NADS
         SUMY = SUMY + DSpe(K)%CD
      END DO
!
      FNA2 = SUMY + S - TADS2
!
      SUMZY = 0.d0
!
      DO K=1,NADS
         SUMZY = SUMZY + DSpe(K)%ZD*DSpe(K)%CD
      END DO
!
      SUMNUY = 0.d0
!
      DO K=1,NADS
         nND = DSpe(k)%n_ND     ! Index of master ad. species in the primary list
         SUMNUY = SUMNUY + DSpe(K)%STQD(nND)*DSpe(K)%CD
      END DO
!
!.....Calculate the new Delta(S) and restrict S to be positive 
!
      DS = -(SUMZY*DPHIP2 + FNA2)/(SUMNUY/S + 1.d0)        ! Calculate  Delta(S)
      S2 = S/2.d0
!
      IF(DS < 0.d0 .AND. DABS(DS) >= S2)  DS = -S2

      S    = S + DS                                      ! Determine  S
      DSRE = DS/S
!
!.....Evaluate the terms for calculating new Delta(PHI)
!
      FNA1 = PHIP2NEW + ADFACTOR*SUMZY
!
      SUMZHZY=0.d0
!
      DO K=1,NADS
         SUMZHZY = SUMZHZY + DSpe(K)%ZD*DSpe(K)%ZD*DSpe(K)%CD
      END DO
!
      SUMZNUY = 0.d0
!
      DO K=1,NADS
         nND = DSpe(k)%n_ND     ! Index of master ad. species in the primary list
         SUMZNUY = SUMZNUY + DSpe(K)%ZD*DSpe(K)%STQD(nND)*DSpe(K)%CD
      END DO
!
!.....Determine new Delta(PHI) and PHI
!
      DPHIP2 = -(ADFACTOR*SUMZNUY*DS/S + FNA1)/&
                (ADFACTOR*SUMZHZY + 1.d0)           ! Caculate Delta (PHI)
!
      PHIP2NEW    = PHIP2NEW + DPHIP2                  ! New PHI
      PSpe(ND)%CP = S
      PHIP2       = PHIP2NEW         
!
!.....Convergence judgement
!
      IF (ITERAD > MAXITPAD)   THEN
!
         IF (IT > 0)      THEN  
            WRITE (32,500)
500         FORMAT(/2X,'ERROR (convergence problem): Reduce time ', &
            'increment', /2X,'  or adjust convergence criteria',&
             ' regarding solving adsorption')
!
            IRETURN = 1
!
            RETURN 
!
                          ELSE
            WRITE(32,505)
505         FORMAT(/2X,'ERROR (convergence problem) in initilization',&      
            ' of adsorption', /2X,'  adjust convergence criteria',&
             ' regarding solving adsorption')
!
            STOP
!
         END IF
!
      END IF
!
      IF (DABS(DPHIP2) <= TOLAD .AND. DABS(DSRE) <= TOLAD) GOTO 499       
!
      GOTO 299
!	     
!.....Complete iteration scheme and call CD_CP and DCD_DCP
!
499   CALL CD_CP
      CALL DCD_DCP
!
!           
      RETURN
!
      END  SUBROUTINE ADSORPTION
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE cd_cp
! 
! ... Modules to be used 
! 
      USE Chem_Dimensions
      USE Aqueous_Species
      USE Surface_Complexes
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*              Calculate the amount of surface complexs               *     
!*                                                                     *
!*                  Version 1.0 - October 12, 2005                     *     
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
      INTEGER(KIND = 8) :: i, k, n, ncp
      INTEGER(KIND = 8) :: icall = 0
      SAVE       icall 
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>>=>=>=>=>=>=>=>=>=>>=>=>=>=>
!
!
! ----------
!.....Calculate the amount of surface complexs
! ----------
!
      do k=1,nads
!
         DSpe(k)%cd = -DSpe(k)%akd
         ncp  = DSpe(k)%ncpad   ! Number of primary species involved in the reaction
!
         do n=1,ncp
!
            i = DSpe(k)%icpad(n) ! Index of primary species in reaction stoichiometry
!	   	    
            DSpe(k)%cd = DSpe(k)%cd + DSpe(k)%stqd(n)&
                         *(dlog10(PSpe(i)%cp)+dlog10(PSpe(i)%gamp)) ! Changed ln to log10    !!!
         end do
!
         DSpe(k)%cd = DSpe(k)%cd + DSpe(k)%zd*phip2
         DSpe(k)%cd = 10.d0**(DSpe(k)%cd)                           ! Changed ln to log10
!
      end do
!
! ----------
!.....Calculate PHI depending on surface potencial 
! ----------
!
      PHIP2NEW = 0.d0
      DO K=1,NADS
         PHIP2NEW = PHIP2NEW - ADFACTOR*DSpe(k)%ZD*DSpe(k)%CD
      END DO
!
!
      return
!
      end  SUBROUTINE cd_cp
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE cd_cp0
! 
! ... Modules to be used 
! 
      USE Chem_Dimensions
      USE Aqueous_Species
      USE Surface_Complexes
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*              Calculate the amount of surface complexs               *     
!*             Similar to CD_CP but without potential term             *
!*                                                                     *
!*                  Version 1.0 - October 12, 2005                     *     
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
      INTEGER(KIND = 8) :: i, k, n, ncp
      INTEGER(KIND = 8) :: icall = 0
      SAVE       icall 
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>>=>=>=>=>=>=>=>=>=>>=>=>=>=>
!
!
      do k=1,nads
!
         DSpe(k)%cd = -DSpe(k)%akd
         ncp  = DSpe(k)%ncpad   ! Number of primary species involved in the reaction
!
         do n=1,ncp 
!
            i = DSpe(k)%icpad(n) ! Index of primary species in reaction stoichiometry
!
            DSpe(k)%cd = DSpe(k)%cd + DSpe(k)%stqd(n)&
                         *(dlog10(PSpe(i)%cp) + dlog10(PSpe(i)%gamp))  
!
         end do

         DSpe(k)%cd = DSpe(k)%cd + DSpe(k)%zd*phip2
         DSpe(k)%cd = 10.d0**(DSpe(k)%cd)       

      end do
!
!
      return
!
      end  SUBROUTINE cd_cp0
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE ADMODEL0
! 
! ... Modules to be used 
! 
      USE Chem_Dimensions
      USE Species_Index
      USE Aqueous_Species
      USE Surface_Complexes
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*             Surface complexs without potential term                 *
!*                                                                     *
!*                  Version 1.0 - October 12, 2005                     *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE                                      
! 
! -------
! ... Double precision variables
! -------
! 
      REAL(KIND = 8) :: SUMCDD
!
! -------
! ... Integer variables
! -------
! 
      INTEGER(KIND = 8) :: i, k, n, ncp
      INTEGER(KIND = 8) :: icall = 0
      SAVE       icall 
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>>=>=>=>=>=>=>=>=>=>>=>=>=>=>
!
!
! ----------
!.....Calculate initial guess of XOH
! ----------
!
      do k=1,nads
!
         DSpe(k)%cd = -DSpe(k)%akd
         ncp  = DSpe(k)%ncpad   ! Number of primary species involved in the reaction
!
         do n=1,ncp
!
            i = DSpe(k)%icpad(n) ! Index of primary species in reaction stoichiometry
!
            if (i == nd)  go to 350
!
            DSpe(k)%cd = DSpe(k)%cd + DSpe(k)%stqd(n)&
                        *(dlog10(PSpe(i)%cp) + dlog10(PSpe(i)%gamp)) 
!
         end do
!
350      continue
!
         DSpe(k)%cd = 10.d0**(DSpe(k)%cd)   
!
      end do
!
!
      SUMCDD = 0.d0
!
      DO k=1,NADS
         SUMCDD = SUMCDD + DSpe(k)%CD
      END DO
!
      PSpe(nd)%CP = TADS2/(1.d0 + SUMCDD)
!
! ----------
!.....Directly calculate CD
! ----------
!
      IADMOD=0             ! Adsorption model 
!
      IF (IADMOD == 0)          THEN
!
         DO k=1,NADS
            DSpe(k)%CD = DSpe(k)%CD*PSpe(nd)%CP
         END DO
!
      END IF
!
      CALL DCD_DCP
!
!              
      RETURN
!
      END  SUBROUTINE ADMODEL0
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE dcd_dcp
! 
! ... Modules to be used 
! 
      USE Chem_Dimensions
      USE Aqueous_Species
      USE Surface_Complexes
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*        The derivative of CD with regards to primary species         *
!*                                                                     *
!*                  Version 1.0 - October 12, 2005                     *     
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
      INTEGER(KIND = 8) :: i, k, n, ncp
      INTEGER(KIND = 8) :: icall = 0
      SAVE       icall 
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>>=>=>=>=>=>=>=>=>=>>=>=>=>=>
!
!       
      do k=1,nads
!
         ncp  = DSpe(k)%ncpad   ! Number of primary species involved in the reaction

         do n=1,ncp
!
            i = DSpe(k)%icpad(n) ! Index of primary species in reaction stoichiometry
!
            DSpe(k)%dcd(i) = DSpe(k)%cd*DSpe(k)%stqd(n)/PSpe(i)%cp
!
         end do
!
      end do
!
!
      return
!
      end  SUBROUTINE dcd_dcp
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE CR_CP(deltex)      ! deltex is tiem step in s
! 
! ... Modules to be used 
! 
      USE MaxChem_Dimensions
      USE Chem_Dimensions
      USE Species_Index
      USE Chem_Constants
      USE Batch_PhysiChem_Conditions
      USE Convergence_Variables
      USE Aqueous_Species
      USE Mineral_Phases
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*          This routine calculates mineral reaction rates             *
!*              and their analytical derivatives                       *     
!*                                                                     *
!*                   Version 1.0 - June 23, 2006                       *     
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
      INTEGER(KIND = 8) :: i, j, k, l, m, n, ik, ncp, ikq, iph
      INTEGER(KIND = 8) :: icall = 0
      SAVE icall
! 
! -------
! ... Double precision variables
! -------
! 
      REAL(KIND = 8) :: eadum, deltex, rkfdum_dis, aH, ratetmp, &
                        sumqk, xmi, rk2xh, deriv1, fdrv1, ratet
! 
! -------
! ... Double precision arrays
! -------
! 
      REAL(KIND = 8), DIMENSION(MaxNmin) :: sgn, rkfdum, fdeltag,&
                       rkhdum, rkohdum, fdeltagh, fdeltagoh,&
                       dumh, dumoh, phterm, deriv, qkterm1,&
                       qkterm2, ssqk10
!
! -------
!.....Still use F77, (does not with F90)                !!!
! -------
!
!      REAL(KIND = 8), DIMENSION(MaxNmin) :: skold
!      SAVE skold
!
!      INTEGER, DIMENSION(MaxNmin) :: nflip
!      SAVE nflip
!
      double precision skold(MaxNmin)
      integer nflip(MaxNmin)
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>>=>>=>=>=>=>=>=>=>=>>=>>=>=>=>=>=>=>=>=>
!
!
      icall = icall + 1
!
! ----------
!.....For flip-flops convergence problems
! ----------
!
      if (iterch <= 1) then                                  !!! Check later  !!!
          do i=1,nmin
             nflip(i) = 0
             skold(i) = 1.d0
          end do
      end if 


!
! ----------
!... .Factor at given temperature tk2 activation energy Ea
! ----------
!
      eadum = ((1.0d0/tk2) - (1.0d0/298.15d0))/rgasm
!
! ----------
!.....Part of rate constant contributed from dependent species
! ----------
!
      if (ndep > 0)   then
         CALL RateAdditionalMechanisms(eadum)
      end if
!
! ----------
!.....Conversion factor
! ----------
!
      sumsalts = 0.d0        ! Sum of salts weights in kg (assume zero for now)
      dliq  = densw2*1.d-3   ! Liquid density in g/cc (kg/l) (assume 1. for now)
!cc      vliq  = (xh2o + sumsalts)/dliq   ! Liquid volume
      vliq  = 1.0d0
      factw = xh2o/vliq      ! Conversion factor = kg h2o liq/liter liquid
      ratet = deltex*phi2*sl2*factw
!
! ----------
!.....The saturation index for kinetic minerals
! ----------
!
      DO_KineticMinerals:&         
      do i=1,nmkin
!
         m = nmequ + i                       ! Stoichiometries ordered 1 to nmin
         MGen(m)%si2k = 10.d0**MGen(m)%si2   ! Convert from log10Q/K to Q/K      May be directly from Index subroutine
         print *, m, MGen(m)%idispre, MGen(m)%ikin
!
! ------------
!........Add skold and nflip stuff to avoid flip-flops with
!........small amounts of minerals if mineral not allowed to ppt.
!........and mineral with saturation window
! ------------
!
         if (MGen(m)%idispre == 1 .or. MKin(i)%ssqk /= 0.d0)   then
!
            if (MGen(i)%si2k >= 1.0d0 .and. skold(i) < 1.0d0)   then 
               nflip(i) = nflip(i) + 1
	      end if
!
            if (nflip(i) > 5)   then
               MGen(i)%si2k = skold(i)+(1.d0-skold(i))/1.1d0
            end if
!
            skold(i) = MGen(i)%si2k
!
         end if 
!
! ------------
!........To avoid overshoot of precipitation rate, give cutoff of Q/K
! ------------
!
!........................................MGen(m)%ikin == 1 for Lasaga TST rate law
         if (MGen(i)%si2k  > 1.0d4 .and. MGen(m)%ikin == 1)   then
	      MGen(i)%si2k = 1.0d4                                      !!!!!
         end if
!
!
!***********************************************************************
!*                                                                     *
!*            Dissolution rates (mol/s, positive values)               *                            
!*                                                                     *
!***********************************************************************
!
! ------------
!........Rate pH dependence - defaults for no dependence
! ------------
!
         phterm(i) = 1.d0     ! Dependence R = R*phterm
         deriv(i)  = 0.d0     ! Deriv of phterm w/resp H+ mol (cp(nh))
!
! ------------
!........Arrays to limit use of power function
! ------------
!
         qkterm1(i) = 0.d0    ! ((q/k)**ck2 - 1)**ck1
         qkterm2(i) = 0.d0    ! (q/k)**ck2
!
!
         IF_Dissolution:&
         if (MGen(i)%si2k < 1.0d0)  then
!
            if (MGen(m)%idispre /= 2)   then
!
!..............First calculate rate constant (neutral mechanish)
!
               rkfdum(i) = MKin(i)%rkf*dexp(-MKin(i)%ea*eadum)&
                         *(10.d0**(MKin(i)%acfdiss      + &
                                   MKin(i)%bcfdiss*tk2  + &
                                   MKin(i)%ccfdiss/tk2))
!
               rkfdum_dis = rkfdum(i)
!
!..............Rate constant contributed from additional mechanisms 
!
               if (ndep > 0)   then 
                  rkfdum(i) = rkfdum(i) + MKin(i)%rkf_ds 
               end if
!
!...........The mineral is only allowed for precipitation
!
            else if (MGen(m)%idispre == 2)  then
!
               rkfdum(i) = 0.d0
!
            end if
!
!...........If the mineral has been exhausted
!
            if (MGen(m)%pre2 < 1.0d-20)  then 
               rkfdum(i) = 0.0d0             
            end if
!
!.......... Functional dependence of rate on delta G-reaction  
!...........Use arrays to reduce use of power function later
!
!
!...........Lasaga TST rate law
!
            if (MGen(m)%ikin == 1)   then
!
               qkterm2(i) = MGen(i)%si2k**MKin(i)%ck2
               fdeltag(i) = 1.0d0 - qkterm2(i)
               qkterm1(i) = fdeltag(i)**MKin(i)%ck1
!
!...........Hellmann-Tisserand rate law 
!...........(Hellmann and Tisserand, Geochimica et Cosmochimica Acta, 70, 364-383, 2006)
!
            else if (MGen(m)%ikin == 2)   then
!
               qkterm2(i) = dlog(MGen(i)%si2k)
               qkterm2(i) = dabs(qkterm2(i))             ! g, take absolute value
               fdeltag(i) = dexp(-1.0d0*MKin(i)%ck1&
                           *qkterm2(i)**MKin(i)%ck2)     !     exp(-ng**m)
               qkterm1(i) = 1.0d0 - fdeltag(i)           ! 1 - exp(-ng**m)
!
            end if 
!
            MKin(i)%rkin2 = rkfdum(i)*MKin(i)%amin2*qkterm1(i)
!
!...........Rate dependence on aH+/aOH-
!
            IF_pH_Dpeendence:&
            if (MKin(i)%idep == 1)   then
!
               aH = PSpe(nh)%gamp*PSpe(nh)%cp
!
               if (aH > MKin(i)%aH1)   then
                  phterm(i) = (aH/MKin(i)%aH1)**MKin(i)%aHexp
                  deriv(i)  = MKin(i)%aHexp*phterm(i)     ! Relative increment
               else if (aH < MKin(i)%aH2)   then
                  phterm(i) = (aH/MKin(i)%aH2)**(-MKin(i)%aOHexp)
                  deriv(i)  = -MKin(i)%aOHexp*phterm(i)   ! Relative increment
               end if
!
               MKin(i)%rkin2 = MKin(i)%rkin2*phterm(i)
!
            end if  IF_pH_Dpeendence
!
         end if  IF_Dissolution
!
!
!***********************************************************************
!*                                                                     *
!*             Precipitation rates (mol/s, negative values)            *                            
!*                                                                     *
!***********************************************************************
!
         ssqk10(i) = 10.d0**MKin(i)%ssqk
!
         IF_Precipitation:&
         if (MGen(i)%si2k > ssqk10(i)) then
!
!...........will precipitate under kin. up to the supersaturation
!...........window specified for equilibrium, then equilibrium takes over
!
            
            ikq = MGen(i)%kineq
            if (ikq > 0) then
               if (MGen(ikq)%si2 > MGen(ikq)%ssq)   then
                  fdeltag(i)    = 1.d0
                  rkfdum(i)     = 0.d0
                  MKin(i)%rkin2 = 0.d0
               end if
!
!
!...........The mineral is only allowed for precipitation
!           idispre= 2 for precipitation; idispre= 3 allowed for both dis. and pre.
!
            else if (MGen(m)%idispre /= 1) then
!
!..............First calculate rate constant (neutral mechanism)
!
               if (MKin(m)%ideprec == 5)   then 
                  rkfdum(i) = rkfdum_dis/10.0d0**MGen(i)%akin  ! for reversible (kpre=kdis/Keq)
!                                          akin(i) was obtained in SUBROUTINE ASSIGN  !!!
!                                          Use akm(  )                !!!!!!!!!!!!
                  go to 199
               end if
!
               rkfdum(i) = MKin(i)%rkprec*dexp(-MKin(i)%eaprec*eadum)&
                         *(10.d0**(MKin(i)%acfprec      + &
                                   MKin(i)%bcfprec*tk2  + &
                                   MKin(i)%ccfprec/tk2))
!
!..............Rate constant contributed from additional mechanisms 
!
               if (ndep > 0)   then 
                  rkfdum(i) = rkfdum(i) + MKin(i)%rkprec_ds  
               end if
!
199            continue
!
!..............Calculate difference between Q and K for bounding precipitation rate
!..............Allow for different precipitation rate laws
!..............Make sure precipitation gives a negative rate (negative sign)
!..............fdeltag is negative.  We make it positive to avoid bombs at
!..............odd powers and for consistency with derivative further below
!
               if (MKin(i)%nplaw == 0)   then
!
!.................Subtract supersaturation constant
!.................Use arrays to reduce use of power function later
!
!
!.................Lasaga TST rate law
!
                  if (MGen(m)%ikin == 1)   then
!
                     qkterm2(i) = MGen(i)%si2k**MKin(i)%ck2prec
                     fdeltag(i) = qkterm2(i)-ssqk10(i)
                     qkterm1(i) = fdeltag(i)**MKin(i)%ck1prec
!
!.................Hellmann-Tisserand rate law 
!
                  else if (MGen(m)%ikin == 2)   then
!
                     qkterm2(i) = dlog(MGen(i)%si2k)            ! g, Positive value for precipitation
                     fdeltag(i) = dexp(-1.0d0*MKin(i)%ck1prec&
                                 *qkterm2(i)**MKin(i)%ck2prec)  !     exp(-ng**m)
                     qkterm1(i) = 1.0d0 - fdeltag(i)            ! 1 - exp(-ng**m)
!
                  end if 
!
                  MKin(i)%rkin2 = -rkfdum(i)*MKin(i)%amin2*qkterm1(i)
!
!
               else if (MKin(i)%nplaw == 1)   then
!
                  fdeltag(i) = MGen(i)%si2k**MKin(i)%ck2prec
!
!.................Try something that drops off quickly as q/k -> 1
!
                  MKin(i)%rkin2 = -(rkfdum(i)*MKin(i)%amin2*&
                                  (fdeltag(i) - (1.d0/fdeltag(i)**2)))
!
               end if
!
!...........The mineral is only allowed for dissolution
!
            else if (MGen(m)%idispre == 1)   then
!
               fdeltag(i)    = 1.d0
               rkfdum(i)     = 0.d0
               MKin(i)%rkin2 = 0.d0
!
            end if
!
!...........Rate dependence on aH+/aOH-
!
            IF_pH_Dependence:&
            if (MKin(i)%ideprec == 1)   then
!
               aH = PSpe(nh)%gamp*PSpe(nh)%cp
!
               if (aH > MKin(i)%aH1p)   then
                  phterm(i) = (aH/MKin(i)%aH1p)**MKin(i)%aHexpp
                  deriv(i)  = MKin(i)%aHexpp*phterm(i)     ! Relative increment - mult by cp(nh)
               else if (aH < MKin(i)%aH2p)   then
                  phterm(i) = (aH/MKin(i)%aH2p)**(-MKin(i)%aOHexpp)
                  deriv(i)  = -MKin(i)%aOHexpp*phterm(i)   !relative increm
               end if
!
               MKin(i)%rkin2 = MKin(i)%rkin2*phterm(i)
!
            end if  IF_pH_Dependence
!
         end if  IF_Precipitation
!
!
!***********************************************************************
!*                                                                     *
!*             Dis/Pre rates at special cases                          *                            
!*                                                                     *
!***********************************************************************
!
!........Mineral at Equilibrium
!
         if (MGen(i)%si2k >= 1.d0 .and. &
             MGen(i)%si2k <= ssqk10(i))   then
!
            rkfdum(i)     = 0.d0
            MKin(i)%rkin2 = 0.d0
!
         end if
!
!........Set maximum rate so not so great to overshoot for dissolution
!
         if (MGen(i)%si2k < 1.0d0 .and. MGen(m)%pre2 .lt. 1e-15) then
            MKin(i)%rkin2 = 0.d0
         else
            if (ratet > 1.d-15) then
               ratetmp = MGen(m)%pre2/ratet
               MKin(i)%rkin2 = min(MKin(i)%rkin2,ratetmp)
            end if
            
         end if
!
!
      end do  DO_KineticMinerals
!
!
!***********************************************************************
!*                                                                     *
!*               Ideal solid solutions (ai=xi)                         *                            
!*                                                                     *
!***********************************************************************
!
!.....Works only for kin minerals and without exponents on q/k terms for now
!.....Calculate mole fractions as fractions of q/k
!.....Note: no additional derivatives needed for this scheme
!
      do n=1,nss
!
         ncp   = ncpss(n)    ! No. of components in solid solution n
         sumqk = 0.d0
!
         do k=1,ncp
            m = icpss(n,k)   ! Mineral index of endmember k in solid sol n 
            i = m-nmequ
            sumqk = sumqk + MGen(i)%si2k
         end do
!
         do k=1,ncp
!
            m   = icpss(n,k)   ! mineral index of endmember k in solid sol n 
            i   = m - nmequ
            xmi = MGen(i)%si2k/sumqk  !mole fraction of endmember
!
!...........correct the reaction rate to reflect solid solution
!
            MKin(i)%rkin2 = MKin(i)%rkin2 + &
                            dabs(MKin(i)%rkin2/qkterm1(i))*(xmi-1.d0)
!
         end do
!
      end do
!
!
!***********************************************************************
!*                                                                     *
!*  Calculate rates in terms of primary species and their derivatives  *                            
!*                                                                     *
!***********************************************************************
!
! ----------
!.....Initialize reates cr and derivative dr
! ----------
!
      do i=1,npri
!
         PSpe(i)%cr = 0.d0         ! cr is moles tied up in kinetic minerals
!
         do j=1,npri
            PSpe(i)%dr(j) = 0.d0   ! dr is derivative of cr i with resprect to j
         end do
!
      end do
!
! ----------
!.....Compute rates cr and their derivetives dr
! ----------
!
      DO_EachKineticMinerals:&
      do i=1,nmkin
!     
         m     = nmequ + i
         ncp   = MGen(m)%ncpm
         rk2xh = MKin(i)%rkin2*xh2o
!
!........For precipitation
! 
         if (MGen(i)%si2k > ssqk10(i) .and. &
             MGen(m)%idispre /= 1)          then
!
            if (MKin(i)%nplaw /= 1)    then
!
!..............Lasaga TST rate law
!
               if (MGen(m)%ikin == 1)   then
!
                  fdrv1 = -rkfdum(i)*MKin(i)%amin2* &
                      MGen(i)%si2k*xh2o&
                     * (MKin(i)%ck1prec*qkterm1(i)/fdeltag(i))&
                     * (MKin(i)%ck2prec*qkterm2(i)/MGen(i)%si2k)
!
!..............Hellmann-Tisserand rate law 
!
               else if (MGen(m)%ikin == 2)   then
                  fdrv1 = - rkfdum(i)*MKin(i)%amin2*xh2o &               ! - K*A
                          *fdeltag(i)*MKin(i)%ck1prec*MKin(i)%ck2prec & ! *exp(-ng**m)*n*m
                          *qkterm2(i)**(MKin(i)%ck2prec - 1.0d0)       ! *g**(m-1)
! 
               end if
!
            else
               fdrv1 = - rkfdum(i)*MKin(i)%amin2* &
                        MGen(i)%si2k*xh2o * &
                        (MKin(i)%ck2prec*fdeltag(i)/MGen(i)%si2k + &
                        2.d0*MKin(i)%ck2prec/(fdeltag(i)**2*MGen(i)%si2k))
            end if
!
!
!........For dissolution
! 
         else if (MGen(i)%si2k     < 1.d0 .and. &
                  MGen(m)%idispre /= 2)        then
!
!...........Lasaga TST rate law
!
            if (MGen(m)%ikin == 1)   then
!
                fdrv1 = - rkfdum(i)*MKin(i)%amin2* &
                     MGen(i)%si2k*xh2o &
                     * (MKin(i)%ck1*qkterm1(i)/fdeltag(i)) &
                     * (MKin(i)%ck2*qkterm2(i)/MGen(i)%si2k)
!
!...........Hellmann-Tisserand rate law 
!
            else if (MGen(m)%ikin == 2)   then
!
               fdrv1 = - rkfdum(i)*MKin(i)%amin2*xh2o &                  
                      * fdeltag(i)*MKin(i)%ck1*MKin(i)%ck2 &
                      * qkterm2(i)**(MKin(i)%ck2 - 1.0d0)
!
              end if
!
!
         end if
!
!
         DO_n: &
         do n=1,ncp
!
            j = MGen(m)%icpm(n)
!
!...........Changes below to account for rates in moles/kgh2o/sec, 
!...........not just moles/sec.  Will not change much unless xh2o varies
!...........much from 1.
!
            PSpe(j)%cr = PSpe(j)%cr + MGen(m)%stqm(n)*rk2xh  ! Rates in terms of Primary species
            iph = 0
!
            DO_k: &
            do k=1,ncp
!
               l = MGen(m)%icpm(k)
               deriv1 = 0.d0       ! Store parts of dr(j,l)
!
!..............Derivative with respect to primary species except water
!
               if (l /= nw)   then
!
!.................Derivatives for precipitation or dissolution, 
!...................skipping suppressed phases
!.................Note that earlier, fdeltag was made always positive
!
                  if (MGen(i)%si2k > ssqk10(i) .and. &
                      MGen(m)%idispre /= 1)           then
!
                      deriv1 = MGen(m)%stqm(k)*MGen(m)%stqm(n)*fdrv1
!
!.....................pH dependent rates - addition for deriv. with respect to H+
!
                      if (l == nh) then  
                         deriv1 = deriv1*phterm(i) + MGen(m)%stqm(n)* &
                                  rk2xh*deriv(i)/phterm(i)                    
                         iph = 1
                      end if
!
                   else if (MGen(i)%si2k < 1.d0 .and. &
                            MGen(m)%idispre /= 2)     then
!
                      deriv1 = MGen(m)%stqm(k)*MGen(m)%stqm(n)*fdrv1
!
!.....................pH dependent rates - addition for deriv. with respect to H+
!
                      if (l == nh)   then 					 
                         deriv1 = deriv1*phterm(i) + MGen(m)%stqm(n) &
                                  *rk2xh*deriv(i)/phterm(i)                    
                         iph = 1 
                      end if
!
                   end if
!
!..................Now add the derivative parts to the total derivative so far 
!
                   PSpe(j)%dr(l) = PSpe(j)%dr(l) + deriv1
! 
!..............Next 3 statements for derivative w/ respect to water
!
               else
!
                  PSpe(j)%dr(l) = PSpe(j)%dr(l) + MGen(m)%stqm(n)*rk2xh    !relative increment scheme
!
               end if
!
            end do  DO_k
!
!...........pH dependent rates - in case H+ is not in mineral stoichiometry
!.......... pH dependence introduces pH in equation
!
            if (iph == 0) then
               PSpe(j)%dr(nh) = PSpe(j)%dr(nh) + MGen(m)%stqm(n)* &
                                rk2xh*deriv(i)/phterm(i)
            end if
!
         end do  DO_n
!
      end do  DO_EachKineticMinerals
!
!
      return
!
      end  SUBROUTINE CR_CP
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE RateAdditionalMechanisms(eadum)
! 
! ... Modules to be used 
! 
      USE Chem_Dimensions
      USE Aqueous_Species
      USE Mineral_Phases
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*        Calculate part of rate constant contributed from             *
!*                  additional mechanisms                              *     
!*                                                                     *
!*                  Version 1.0 - October 19, 2005                     *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      IMPLICIT NONE                                      
! 
! -------
! ... Double precision variables
! -------
! 
      REAL(KIND = 8) :: eadum, sum_an, cj, aj, gamj
!
! -------
! ... Integer variables
! -------
! 
      INTEGER(KIND = 8) :: i, j, isp, nsp, js
      INTEGER(KIND = 8) :: icall = 0
      SAVE       icall 
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>>=>>=>=>=>=>=>=>=>=>>=>>=>=>=>=>=>=>=>=>
!
!
      icall = icall + 1
!
! ----------
!.....For species dependent rate law: H+,  k(h+), expo(h+). term.
!.................................... oh-, k(oh-), expo(oh-). term.
!.................................... nadum2, rkds,  expds
!.................................... k=k(h+)*[H+]**expo +....
! ----------
!
      DO_KineticMinerals: &
      do i=1,nmkin
!
! -------------
!........For dissolution
! -------------
!
         IF_Dissolution: &
         if (MKin(i)%idep == 2)  then
!............MKin represents class of kinetic properties of minerals
!
            MKin(i)%rkf_ds = 0.0d0   
!
            DO_DisMechanisms: &
            do j=1,MKin(i)%ndis        
!
               sum_an = 0.0d0
               nsp = MKin(i)%nspds(j)     ! Number of species involved in one mechanism
!
               do isp=1,nsp
                  js     = MKin(i)%ids(j,isp)
                  cj     = ASpe(js)%ct
                  gamj   = ASpe(js)%gamt
                  aj     = cj*gamj
                  sum_an = sum_an+aj**MKin(i)%expds(j,isp)
               end do
!
               MKin(i)%rkf_ds = MKin(i)%rkf_ds + &
                                MKin(i)%rkds(j)  &
                                *dexp(-MKin(i)%eads(j)*eadum)*sum_an
!
            end do  DO_DisMechanisms
!
         end if  IF_Dissolution
!
! -------------
!........For precipitation
! -------------
!
         IF_Precipitation: &
         if (MKin(i)%ideprec == 2)   then
!
            MKin(i)%rkprec_ds = 0.0d0   

            DO_PreMechanisms: &
            do j=1,MKin(i)%npre
!
               sum_an = 0.0d0
               nsp    = MKin(i)%nsppr(j)  ! Number of species involved in one mechanism
!
               do isp=1,nsp
                  js     = MKin(i)%idsp(j,isp)
                  cj     = ASpe(js)%ct
                  gamj   = ASpe(js)%gamt
                  aj     = cj*gamj
                  sum_an = sum_an + aj**MKin(i)%expdsp(j,isp)
               end do
!
               MKin(i)%rkprec_ds = MKin(i)%rkprec_ds + &
                                   MKin(i)%rkdsp(j) &
                                  *dexp(-MKin(i)%eadsp(j)*eadum)*sum_an
!
            end do  DO_PreMechanisms
!
         end if  IF_Precipitation
!
      end do  DO_KineticMinerals
!
!
      return
!
      end  SUBROUTINE RateAdditionalMechanisms
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
