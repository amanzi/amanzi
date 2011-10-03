! :: ----------------------------------------------------------
! :: COREREACT VODE-based solver
! :: ----------------------------------------------------------
      subroutine FORT_CHEM_2ND(ielem,pth,dim_pth,ut,dim_ut,c,dim_c, &
                   pre, dim_pre, amin, dim_amin, rad, dim_rad, &
                   pfug, dim_pfug,nc,deltex,funcCount)

      USE MaxChem_Dimensions
      USE Chem_Dimensions
      USE Species_Index
      USE Chem_Constants
      USE Batch_PhysiChem_Conditions
      USE Convergence_Variables
      USE CutOff_Threshold
      USE ChemOption_Variables
      USE Aqueous_Species
      USE Mineral_Phases
      USE Gaseous_Species

      integer nc,n
      integer(8) ielem
      integer(8) dim_pth,dim_ut, dim_c, dim_pre
      integer(8) dim_amin, dim_rad, dim_pfug
      REAL(KIND=8)  pth(20)
      REAL(KIND=8)  ut(dim_ut), c(dim_c), pre(dim_pre)
      REAL(KIND=8)  amin(dim_amin), rad(dim_rad)
      REAL(KIND=8)  pfug(dim_pfug)
      REAL(KIND=8)  deltex
      REAL(KIND=8)  funcCount

      external FEX_2ND, JEX_2ND
      integer  NEQ, ITOL,ITASK,ISTATE,IOPT,LRW,LIW,MF
      integer  IPAR(1)
      REAL(KIND=8) T, TOUT, RTOL
      REAL(KIND=8) Y(dim_ut+dim_pre), YOLD(dim_ut+dim_pre)
      REAL(KIND=8) ATOL(dim_ut+dim_pre), RPAR(1)

      integer, POINTER, SAVE :: IWORK(:)
      REAL(KIND=8), POINTER, SAVE :: RWORK(:)

      NEQ    = npri+1
      LRW    = 22+9*NEQ+2*NEQ**2
      LIW    = 30+NEQ
      IPAR   = 0
      T      = 0
      TOUT   = deltex
      ITOL   = 1
      ITASK  = 1
      ISTATE = 1
      IOPT   = 0
      MF     = 22

      RTOL   = 1.d-12
      ATOL   = 1.d-12

      if (.NOT. ASSOCIATED(RWORK)) then
         allocate(RWORK(22+9*nc+2*nc**2))
         allocate(IWORK(30+nc))
      end if
!     assign the inputs to appropriate structure in modules
      call Pre_Chem(pth,ut,c,pre,rad,amin,pfug)
      do n = 1,npri
         Y(n) = PSpe(n)%tt
      end do
      Y(NEQ) = pre(1)

      call FEX_2ND (NEQ, T, Y, YOLD, RPAR, IPAR)
      funcCount = 0.0
      if (dabs(YOLD(NEQ)) .gt. 0.d0) then
 99       call DVODE(FEX_2ND,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE, &
                     IOPT,RWORK,LRW,IWORK,LIW,JEX_2ND,MF,RPAR,IPAR) 
         funcCount = funcCount + dble(IWORK(12))
!         print *, "funcCount =", IWORK(12)
         if (ISTATE .ne. 2) then
            print *, "ISTATE = ", ISTATE
            ISTATE = 1
            stop
            goto 99
         end if
      end if
      call AQEOUS_EQB_2ND(Y,NEQ,100)
!      call Batch_Chem_2ND(Y,NEQ,100)

      call Post_Chem(pth,ut,c,pre,rad,amin,pfug)    
      pre(1) = Y(NEQ)

      end

! :: ----------------------------------------------------------
! :: Compute rate of reaction for each components.
! :: ----------------------------------------------------------
      subroutine FEX_2ND (NEQ, T, Y, YDOT, RPAR, IPAR)

      USE Chem_Dimensions
      USE Species_Index
      USE Chem_Constants
      USE Batch_PhysiChem_Conditions
      USE Aqueous_Species
      USE Mineral_Phases
      USE Gaseous_Species

      integer NEQ, IPAR(7),m,imk
      REAL(KIND=8)  T, Y(NEQ), YDOT(NEQ)
      REAL(KIND=8)  RPAR(IPAR(1))
      REAL(KIND=8)  actfrc, vfrtmp, anucl, rkmol, slf, amin2(nmkin)
      integer(8) i
      REAL(KIND=8) deltex

      call AQEOUS_EQB_2ND(Y,NEQ,100)
!      call Batch_Chem_2ND(Y,NEQ,100)

      if((a_fmr2 .lt. sl2) .and. (sl2 .gt. 0.d0) .and.  &
     	 (a_fmr2 .gt. 0.d0)) then
         actfrc = a_fmr2/(phi2*densw2*sl2)
      else
         actfrc = 1.d0/(phi2*densw2)
      end if

      vfrtmp = 1.0d0 - phi2

      do m=1,nmkin
         amin2(m) = MKin(m)%amin2
         imk = nmequ + m
         if (MKin(m)%rad2 .gt. 0.d0 .and. Y(NEQ) .lt. 0.d0)  then
            anucl = (ppi*0.125d0/MKin(m)%rad2)
            MKin(m)%amin2 = anucl*actfrc
         else
!           Assign mol/L medium based on mineral amount or by preset 
!           minimum volume fraction (rnucl)
            rkmol = dmax1(Y(NEQ),(MKin(m)%rnucl*vfrtmp/MGen(imk)%vmin))

!           Calculate area in m^2/kg H2O
            MKin(m)%amin2 = amin2(m)*rkmol*MGen(imk)%vmin*actfrc

         end if
      end do

      call REACTION_RATE(deltex) 

      do m = 1,nmkin
         MKin(m)%amin2 = amin2(m)
      end do

      if (nmkin .gt. 0) then
         do i=1,npri
            YDOT(i) = PSpe(i)%cr            
        end do
        YDOT(NEQ) = -Mkin(1)%rkin2*phisl2*PSpe(nw)%cp

      end if

      end

! :: ----------------------------------------------------------
! :: Construct Jacobian.  Currently a dummy function.
! :: ----------------------------------------------------------
      subroutine JEX_2ND(NEQ, T, Y, ML, MU, PD, NRPD, RPAR, IPAR)

      integer NEQ, ML,MU, NRPD, IPAR(7)
      REAL(KIND=8) PD(NRPD,NEQ), RPAR(IPAR(1)), T, Y(NEQ)

      end


! :: ----------------------------------------------------------
! :: Solve the nonlinear species equilibrium problem
! ::   (1) A full Newton approach would be too tedious due to 
! ::       various definitions of activity coefficients.
! :: ----------------------------------------------------------
      subroutine Batch_Chem_2ND(Y,NEQ,it)

      USE MaxChem_Dimensions
      USE Chem_Dimensions
      USE Species_Index
      USE Chem_Constants
      USE Batch_PhysiChem_Conditions
      USE Convergence_Variables
      USE CutOff_Threshold
      USE Aqueous_Species
      USE Mineral_Phases
      USE Gaseous_Species
      USE ChemOption_Variables
      USE Batch_Chem_Jacobian

      integer it, NEQ
      REAL(KIND=8)  Y(NEQ)

      integer i, j, k, ncp, n, iter, iter2
      integer(8) no_ch
      REAL(8)  csjxh, err_cs

      logical plusminus
      integer ns_nwt
      integer INFO
      integer IPIV(npri)
      REAL(8) y_nwt(npri), tmp(npri)
      REAL(8) J_nwt(npri,npri)
      REAL(8) dcsdcp(naqx,npri), newt_alpha
      REAL(8) cptmp(npri),tmp2

      ns_nwt = npri

      no_ch  = 0
      call Get_ActivityCoefficient(no_ch)

      xh2o = Y(nw)/rmh2o
      MGen(1)%pre2 = Y(npri+1)

      do iter = 1,it  
         J_nwt = 0.d0
         y_nwt = 0.d0
         do i = 1,npri
            J_nwt(i,i) = 1.d0
!            if (dabs(PSpe(i)%cp) .lt. 1.d-14) then
!               PSpe(i)%cp = 0.d0;
!            end if
            y_nwt(i)   = PSpe(i)%cp - Y(i)/xh2o
         end do
         y_nwt(nw) = PSpe(nw)%cp - xh2o
      
         dcsdcp = 0
         do j=1,naqx
            ncp   = SSpe(j)%ncps
            csjxh = SSpe(j)%cs
            tmp2 = 0.d0
            do n = 1,ncp
               tmp2 = tmp2 + SSpe(j)%stqs(n)
            end do
            do n = 1,ncp
               i = SSpe(j)%icps(n)
	       dcsdcp(j,i) = tmp2*csjxh*PSpe(i)%dgamp(i)/PSpe(i)%gamp - &
      	          csjxh*SSpe(j)%dgams(i)/SSpe(j)%gams
               if (PSpe(i)%cp .gt. 1.d-14 ) then
                  dcsdcp(j,i)=dcsdcp(j,i)+ &
                   csjxh*SSpe(j)%stqs(n)/PSpe(i)%cp
!               else
!                  dcsdcp(j,i) = csjxh*(SSpe(j)%stqs(n)*
!     &                 (PSpe(i)%dgamp(i)/PSpe(i)%gamp + 1.d0/PSpe(i)%cp) 
!     &                 -SSpe(j)%dgams(i)/SSpe(j)%gams)
!                  dcsdcp(j,i) = (csjxh/PSpe(i)%cp)*(SSpe(j)%stqs(n)*
!     &                 (PSpe(i)%cp*PSpe(i)%dgamp(i)/PSpe(i)%gamp + 1.d0) 
!     &                 -PSpe(i)%cp*SSpe(j)%dgams(i)/SSpe(j)%gams)
               end if

!               if (dabs(dcsdcp(j,i)) .gt. 1e10) then
!		  print *, i,j,dcsdcp(j,i),PSpe(i)%cp
!		  print *, csjxh, SSpe(j)%stqs(n)
!		  print *, PSpe(i)%dgamp(i),PSpe(i)%gamp
!		  print *, SSpe(j)%dgams(i),SSpe(j)%gams
!	       end if
               y_nwt(i) = y_nwt(i) + SSpe(j)%stqs(n)*csjxh
            end do
         end do

         do i = 1,npri
            do j = 1,naqx
               ncp = SSpe(j)%ncps
               do n = 1,ncp
                  k = SSpe(j)%icps(n)
                  J_nwt(k,i) = J_nwt(k,i) + SSpe(j)%stqs(n)*dcsdcp(j,i)
               end do
            end do
         end do
            
         y_nwt = - y_nwt

         err_cs =  0.d0
         do i = 1,ns_nwt
            err_cs = err_cs + y_nwt(i)*y_nwt(i)
         end do 
         err_cs = err_cs/dble(ns_nwt)
         err_cs = sqrt(err_cs)
         do i = 1,npri
            cptmp(i) = PSpe(i)%cp
         end do
!         print *, iter, err_cs
        
         if (err_cs .le. 1d-16) then 
            goto 199
         end if

!        solve the linear system of equations.
         call DGETRF(ns_nwt,ns_nwt,J_nwt,ns_nwt,IPIV,INFO)
         if (INFO .ne. 0) then
            print*, INFO
            print*, J_nwt
            J_nwt = 0.0
            do i = 1,npri
               do j = 1,naqx
                  ncp = SSpe(j)%ncps
                  do n = 1,ncp
                     k = SSpe(j)%icps(n)
                     J_nwt(k,i)=J_nwt(k,i)+SSpe(j)%stqs(n)*dcsdcp(j,i)
                     print*,k,i,J_nwt(k,i)
                     print*,n,j,SSpe(j)%stqs(n),dcsdcp(j,i)
                  end do
               end do
            end do
            stop
         end if
         call DGETRS('N',ns_nwt,1,J_nwt,ns_nwt,IPIV,y_nwt,ns_nwt,INFO)


!        adjust Newton step, if needed.
         newt_alpha = 2.0
 119     plusminus = .false.         
         newt_alpha = 0.5d0*newt_alpha
         do i = 1,npri
            cptmp(i) = PSpe(i)%cp + newt_alpha*y_nwt(i)
            if (cptmp(i) < 0.d0) then
               plusminus = .true.               
            end if
         end do
         if (plusminus) then
            goto 119
         else
            do i = 1,npri
               PSpe(i)%cp = cptmp(i)
            end do
         end if
         
!        loop twice since the activity coefficients are functions of SSpe()%cs
         do iter2 = 1,2
            call Get_ActivityCoefficient(no_ch)
            do j = 1,naqx
               SSpe(j)%cs = 1.d0
               ncp = SSpe(j)%ncps
               do n = 1,ncp
                  i = SSpe(j)%icps(n)
                  SSpe(j)%cs = SSpe(j)%cs* &
                      (PSpe(i)%gamp*PSpe(i)%cp)**SSpe(j)%stqs(n)
               end do
               SSpe(j)%cs = SSpe(j)%cs/(SSpe(j)%gams*10**SSpe(j)%aks)
            end do
         end do
      end do 
      
      if (iter .ge. it) then
	 print *, "Equilibrium solve has not converged:  ", it, err_cs
      end if	
 199  do i=1,npri   
         ASpe(i)%gamt      = PSpe(i)%gamp  
      end do
      do i=1,naqx                  
         ASpe(npri+i)%gamt = SSpe(i)%gams                   
      end do

      PSpe(nw)%cp = xh2o

      end

! :: ----------------------------------------------------------
! :: Solve the nonlinear species equilibrium problem
! ::   (1) A full Newton approach would be too tedious due to 
! ::       various definitions of activity coefficients.
! ::   (2) Simplified the algorithm, and removed fluff.
! :: ----------------------------------------------------------
      subroutine AQEOUS_EQB_2ND(Y,NEQ,it)

      USE Chem_Dimensions
      USE Species_Index
      USE Chem_Constants
      USE Aqueous_Species
      USE Mineral_Phases

      integer it, NEQ
      REAL(8)  Y(NEQ)

      integer i, j, k, ncp, n, iter, iter2
      integer(8) no_ch
      REAL(8) err_cs

      logical plusminus
      integer ns_nwt
      integer INFO
      integer IPIV(npri)
      REAL(8) y_nwt(npri)
      REAL(8) J_nwt(npri,npri)
      REAL(8) dcsdcp(naqx,npri), newt_alpha
      REAL(8) cptmp(npri), cstmp(naqx)

      ns_nwt = npri

      no_ch  = 0
      call Get_ActivityCoefficient(no_ch)

      xh2o = Y(nw)/rmh2o
      MGen(1)%pre2 = Y(npri+1)

      do iter = 1,it  

         J_nwt = 0.d0
         y_nwt = 0.d0
         dcsdcp = 0.d0
!        Determine residual
         do i = 1,npri
            y_nwt(i) = PSpe(i)%cp - Y(i)/xh2o
         end do
         y_nwt(nw) = PSpe(nw)%cp - xh2o
         do j =1,naqx
            ncp = SSpe(j)%ncps
            cstmp(j) = 1.d0
            do n=1,ncp
               i = SSpe(j)%icps(n)
               cstmp(j) = cstmp(j)* &
                 (PSpe(i)%gamp*PSpe(i)%cp)**SSpe(j)%stqs(n)
            end do
            cstmp(j) = cstmp(j)/(SSpe(j)%gams*10**SSpe(j)%aks)
         end do
         do j =1,naqx
            ncp = SSpe(j)%ncps
            do n=1,ncp
               i = SSpe(j)%icps(n)
!               y_nwt(i) = y_nwt(i) + SSpe(j)%stqs(n)*SSpe(j)%cs
               y_nwt(i) = y_nwt(i) + SSpe(j)%stqs(n)*cstmp(j)
            end do
         end do
         

!        Check residual norm
         err_cs =  0.d0
         do i = 1,ns_nwt
            err_cs = err_cs + y_nwt(i)*y_nwt(i)
         end do 
         err_cs = sqrt(err_cs/dble(ns_nwt))
!         print *, iter, err_cs

         if (err_cs .le. 1d-14) then 
            goto 199
         else
            y_nwt = - y_nwt
         end if
         dcsdcp = 0
         do j = 1,naqx
            ncp = SSpe(j)%ncps
            do n = 1,ncp
               i = SSpe(j)%icps(n)
               if (PSpe(i)%cp .gt. 1.d-14) then
                  dcsdcp(j,i) = dcsdcp(j,i) + &
                    SSpe(j)%stqs(n)*SSpe(j)%cs/PSpe(i)%cp
               end if
            end do
         end do

         do i = 1,npri
            J_nwt(i,i) = J_nwt(i,i) + 1.d0
            do j = 1,naqx
               ncp = SSpe(j)%ncps
               do n = 1,ncp
                  k = SSpe(j)%icps(n)
                  J_nwt(k,i) = J_nwt(k,i)+ SSpe(j)%stqs(n)*dcsdcp(j,i)
               end do
            end do
         end do
        
!        solve the linear system of equations.
         call DGETRF(ns_nwt,ns_nwt,J_nwt,ns_nwt,IPIV,INFO)
         if (INFO .ne. 0) then
            print*, INFO
            print*, J_nwt
            J_nwt = 0.0
            do i = 1,npri
               do j = 1,naqx
                  ncp = SSpe(j)%ncps
                  do n = 1,ncp
                     k = SSpe(j)%icps(n)
                     J_nwt(k,i)=J_nwt(k,i)+SSpe(j)%stqs(n)*dcsdcp(j,i)
                     print*,k,i,J_nwt(k,i)
                     print*,n,j,SSpe(j)%stqs(n),dcsdcp(j,i)
                  end do
               end do
            end do
            stop
         end if
         call DGETRS('N',ns_nwt,1,J_nwt,ns_nwt,IPIV,y_nwt,ns_nwt,INFO)

!        adjust Newton step, if needed.
         newt_alpha = 2.0
 119     plusminus = .false.         
         newt_alpha = 0.5d0*newt_alpha
         do i = 1,npri
            cptmp(i) = PSpe(i)%cp + newt_alpha*y_nwt(i)
            if (cptmp(i) < 0.d0) then
               plusminus = .true.               
            end if
         end do
         if (plusminus) then
            goto 119
         else
            do i = 1,npri
               PSpe(i)%cp = cptmp(i)
            end do
         end if
         
!        loop twice since the activity coefficients are functions of SSpe()%cs
         do iter2 = 1,2
            call Get_ActivityCoefficient(no_ch)
            do j = 1,naqx
               SSpe(j)%cs = 1.d0
               ncp = SSpe(j)%ncps
               do n = 1,ncp
                  i = SSpe(j)%icps(n)
                  SSpe(j)%cs = SSpe(j)%cs* &
                       (PSpe(i)%gamp*PSpe(i)%cp)**SSpe(j)%stqs(n)
               end do
               SSpe(j)%cs = SSpe(j)%cs/(SSpe(j)%gams*10**SSpe(j)%aks)
            end do
         end do
      end do 
      
      if (iter .ge. it) then
	 print *, "Equilibrium solve has not converged. ", it, err_cs
      end if	

 199  do i=1,npri   
         ASpe(i)%gamt      = PSpe(i)%gamp  
      end do
      do i=1,naqx                  
         ASpe(npri+i)%gamt = SSpe(i)%gams                   
      end do

      PSpe(nw)%cp = xh2o

      end
     
! :: -----------------------------------------------------
! :: Determine the reaction rate. Adapted from ChemBatch.f
! ::   (1) omitted everything on derivatives of reaction rate
! ::   (2) pH dependence included.
! ::   (3) no equilibrium reaction involving mineral.
! ::   (4) ignore the "additional mechanisms"
! :: -----------------------------------------------------
      SUBROUTINE REACTION_RATE(deltex) 
 
      USE Chem_Dimensions
      USE Species_Index
      USE Chem_Constants
      USE Batch_PhysiChem_Conditions
      USE Aqueous_Species
      USE Mineral_Phases

      IMPLICIT NONE                                      
 
      INTEGER(8) i, j,  m, n, ncp 
      REAL(8) deltex, ratetmp, ratet
      REAL(8) eadum, aH,  phterm, rk2xh, st
      REAL(8) fdeltag, qkterm1, qkterm2, rkfdum


! ::: remove temperature dependence for now.
!     eadum = 1.d0/tk2 - 1.d0/298.15/rgasm
      eadum = 0.d0

! ::: conversion factor, not sure if any use.
      sumsalts = 0.d0        
      dliq  = densw2*1.d-3   
      vliq  = 1.0d0
      factw = xh2o/vliq      
      ratet = deltex*phi2*sl2*factw

! ::: Loop over minerals
      do i=1,nmkin

! :::::: Compute mineral saturation ratio
         st = 1.d0/(10**MGen(i)%akm)
         ncp = MGen(i)%ncpm

         do n=1,ncp
            j = MGen(i)%icpm(n)
            st = st*(PSpe(j)%cp*PSpe(j)%gamp)**MGen(i)%stqm(n)
         end do
         
         MGen(i)%si2k = st   
         MKin(i)%rkin2 = 0.d0

! :::::: Case 1: solution is not saturated and there are minerals
!                => dissolution
         if (MGen(i)%si2k < 1.0d0 .and. MGen(i)%pre2 > 1.0d-16 .and. &
             MGen(i)%idispre /=2)  then

            rkfdum = MKin(i)%rkf

!           Lasaga TST rate law
            if (MGen(i)%ikin == 1)   then
               qkterm2 = MGen(i)%si2k**MKin(i)%ck2
               fdeltag = dabs(1.0d0 - qkterm2)
               qkterm1 = fdeltag**MKin(i)%ck1

!           Hellmann-Tisserand rate law 
            else if (MGen(i)%ikin == 2)   then
               qkterm2 = dlog(MGen(i)%si2k)
               qkterm2 = dabs(qkterm2)             
               fdeltag = dexp(-1.0d0*MKin(i)%ck1 &
                            *qkterm2**MKin(i)%ck2)     
               qkterm1 = 1.0d0 - fdeltag     
            end if

            MKin(i)%rkin2 = rkfdum*MKin(i)%amin2*qkterm1

!           Dependence on pH
            if (MKin(i)%idep == 1) then
               aH = PSpe(nh)%gamp*PSpe(nh)%cp
               if (aH > MKin(i)%aH1)   then
                  phterm = (aH/MKin(i)%aH1)**MKin(i)%aHexp
               else if (aH < MKin(i)%aH2)   then
                  phterm = (aH/MKin(i)%aH2)**(-MKin(i)%aOHexp)
               end if
               MKin(i)%rkin2 = MKin(i)%rkin2*phterm
            end if


! :::::: Case 2: solution is saturated => precipitation
         elseif (MGen(i)%si2k .ge. 1.0 .and. MGen(i)%idispre /= 1) then

            rkfdum = MKin(i)%rkprec

            if (MKin(i)%nplaw == 1) then
               fdeltag = MGen(i)%si2k**MKin(i)%ck2prec
               MKin(i)%rkin2 = -rkfdum*MKin(i)%amin2* &
                               (fdeltag-1.d0/fdeltag**2)
            else
!              Lasaga TST rate law
               if (MGen(i)%ikin == 1)   then
                  qkterm2 = MGen(i)%si2k**MKin(i)%ck2prec
                  fdeltag = dabs(1.0-qkterm2)
                  qkterm1 = fdeltag**MKin(i)%ck1prec
                  
!              Hellmann-Tisserand rate law 
               else if (MGen(i)%ikin == 2)   then
                  qkterm2 = dlog(MGen(i)%si2k)     
                  qkterm2 = dabs(qkterm2)
                  fdeltag = dexp(-1.0d0*MKin(i)%ck1prec &
                              *qkterm2**MKin(i)%ck2prec) 
                  qkterm1 = 1.0d0 - fdeltag  
               end if
               MKin(i)%rkin2 = -rkfdum*MKin(i)%amin2*qkterm1
            end if

!           Dependence on pH
            if (MKin(i)%ideprec == 1)   then
               aH = PSpe(nh)%gamp*PSpe(nh)%cp
               if (aH > MKin(i)%aH1p) then
                  phterm = (aH/MKin(i)%aH1p)**MKin(i)%aHexpp
               else if (aH < MKin(i)%aH2p)   then
                  phterm = (aH/MKin(i)%aH2p)**(-MKin(i)%aOHexpp)
               end if
               MKin(i)%rkin2 = MKin(i)%rkin2*phterm
            end if

         end if

! :::::: This check should not be necessary => delete in next cleanup
!        Set maximum rate so not so great to overshoot for dissolution
!         if (MGen(i)%idispre /= 2 .and. MGen(i)%si2k < 1.0d0) then
!            if (ratet > 1.d-15) then
!               ratetmp = MGen(i)%pre2/ratet
!               MKin(i)%rkin2 = min(MKin(i)%rkin2,ratetmp)
!               print *, "dissolution here... exiting ", m, 
!     &              MGen(i)%idispre, MGen(i)%ikin
!               stop
!            end if
!         end if

      end do 

! ::: Calculate rates in terms of primary species and their derivatives  *  
      do i=1,npri
         PSpe(i)%cr = 0.d0
      end do

!     Compute rates cr 
      do i=1,nmkin
         m     = nmequ + i
         ncp   = MGen(i)%ncpm
         rk2xh = MKin(i)%rkin2*xh2o

         do n=1,ncp
            j = MGen(i)%icpm(n)
            PSpe(j)%cr = PSpe(j)%cr + MGen(i)%stqm(n)*rk2xh  
         end do 

      end do 

      end  

! :: ----------------------------------------------------------
! :: Initialize Chem Modules
! :: ----------------------------------------------------------
      subroutine INITIALIZE_CHEM_MODULES

      USE MaxChem_Dimensions
      USE Chem_Dimensions
      USE Species_Index
      USE Chem_Constants
      USE Batch_PhysiChem_Conditions
      USE Convergence_Variables
      USE CutOff_Threshold
      USE ChemOption_Variables
      USE Aqueous_Species
      USE Mineral_Phases
      USE Gaseous_Species
      USE Batch_Chem_Jacobian
      
      call Read_Chem
      MGen(1)%idispre = 2 
      end

! :: ---------------------------------------------------------
! :: Assign value to data structure. Adapted from ChemCouple.f
! :: ---------------------------------------------------------
      SUBROUTINE PRE_CHEM(pth,ut,conc,pre,rad,amin,pfug)

      USE Chem_Dimensions
      USE Chem_Constants
      USE CutOff_Threshold
      USE Convergence_Variables
      USE Batch_PhysiChem_Conditions
      USE Species_Index
      USE Aqueous_Species
      USE Gaseous_Species
      USE Mineral_Phases

      IMPLICIT NONE                                      

      INTEGER(KIND = 8) :: i,j, k, m, n, ncp, nnn, imk
      INTEGER(KIND = 8) :: icall = 0
      REAL(KIND = 8), DIMENSION(20)    :: pth                
      REAL(KIND = 8), DIMENSION(npri)  :: ut              
      REAL(KIND = 8), DIMENSION(naqt)  :: conc             
      REAL(KIND = 8), DIMENSION(nmin)  :: Pre
      REAL(KIND = 8), DIMENSION(nmin)  :: Pre0
      REAL(KIND = 8), DIMENSION(nmkin) :: Rad  
      REAL(KIND = 8), DIMENSION(nmkin) :: Amin
      REAL(KIND = 8), DIMENSION(ngas)  :: pfug 

!     Get thermo-physical-chemical condition variables for the batch system
      Pt2    = pth(1)    
      tc2    = pth(2) 
      sl2    = pth(3)
      phi2   = pth(8)
      denswa = pth(9)
      tk2    = tc2 + 273.15d0
      sg2    = 1.0d0-sl2  
      densw2 = denswa*1.0d3 
      a_fmr2 = 1.0d0  
      phisl2 = phi2*sl2
      xh2o   = denswa
      CALL assign
      do i=1,npri
         PSpe(i)%tt = ut(i)           
         PSpe(i)%u2 = conc(i)*xh2o      
      end do
      PSpe(nw)%u2= xh2o*rmh2o   
      
      do j=1,naqx
         ncp = SSpe(j)%ncps
         m=npri+j
         do n=1,ncp
            k = SSpe(j)%icps(n)
            PSpe(k)%u2 = PSpe(k)%u2 + SSpe(j)%stqs(n)*conc(m)*xh2o
         enddo
      enddo

      do n=1,npri
         PSpe(n)%cp = conc(n)
      end do
      do n=1,naqx
         SSpe(n)%cs = conc(npri+n)
      end do
      PSpe(nw)%cp = denswa

      do m=1,nmin
         MGen(m)%cm = 0.0d0
      end do

      do m=1,nmin
         MGen(m)%pre2  = pre (m) 
      end do

      do m=1,nmkin
         MKin(m)%rad2  = rad (m) 
         MKin(m)%amin2 = amin(m) 
      end do

      end 

! ::: -----------------------------------------------------------
! ::: Assign value from data structure. Adapted from ChemCouple.f
! ::: -----------------------------------------------------------
      SUBROUTINE Post_Chem(pth,ut,C,pre,rad,amin,pfug) 

      USE Chem_Dimensions
      USE Chem_Constants,              only : rmh2o
      USE Batch_PhysiChem_Conditions
      USE Species_Index
      USE Aqueous_Species
      USE Gaseous_Species
      USE Mineral_Phases
      USE Ion_Exchange
      USE Surface_Complexes

      IMPLICIT NONE                                      

      INTEGER(KIND = 8) :: i, j, k, m, n, ncp, nkkn, ip_yes, nhs
      REAL(KIND = 8) :: phislv, phisvd, pretmp, utemp

      REAL(KIND = 8), DIMENSION(20)   :: PTH
      REAL(KIND = 8), DIMENSION(Npri) :: UT
      REAL(KIND = 8), DIMENSION(Naqt) :: C 
      REAL(KIND = 8), DIMENSION(Nmin) :: Pre
      REAL(KIND = 8), DIMENSION(Nmin) :: Pre0
      REAL(KIND = 8), DIMENSION(Nmkin) :: Rad
      REAL(KIND = 8), DIMENSION(Nmkin) :: Amin 
      REAL(KIND = 8), DIMENSION(Ngas) :: GP
      REAL(KIND = 8), DIMENSION(Ngas) :: Pfug 

      dliq     = densw2/1000.d0                  
      pth(9)   = dliq
      sumsalts = 0.d0   
      vliq     = 1.0d0  
      factw    = PSpe(nw)%cp/vliq 
      phislv = phisl2/vliq

      do n=1,npri
         c(n) = PSpe(n)%cp
         ut(n) = c(n) 
      end do
      c(nw) = rmh2o
      ut(nw) = rmh2o

      do n=1,naqx
         c(npri+n) = SSpe(n)%cs
      end do

      do j=1,naqx
         ncp = SSpe(j)%ncps
         do k=1,ncp
            n = SSpe(j)%icps(k)
            utemp = ut(n) + SSpe(j)%stqs(k)*c(npri+j)
            ut(n) = utemp
         end do
      end do

      do n=1,npri
         ut(n) = ut(n)*factw              
      end do


!    Calculate pH value
      if (nh > 0 .and. nh <= npri)  then
	   ph2 = -dlog10(PSpe(nh)%gamp*PSpe(nh)%cp)    
      end if
      if (nh > npri) then
         nhs = nh - npri
	   ph2 = -dlog10(SSpe(nhs)%gams*SSpe(nhs)%cs)                
      end if
      pth(11) = ph2      
!     Save ionic strength and water activity
      pth(12) = str            
      pth(13) = PSpe(nw)%gamp  

      do m=1,nmequ
         pre (m) = pre(m) + MGen(m)%cm*phislv
         if (pre(m) <= 1.0d-30)  pre(m) = 0.0d0
      end do

      do m=1,nmkin      
         nkkn   = nmequ + m
         pretmp = pre(nkkn) 
         pre(nkkn) = pretmp
         if (pretmp <= 1.0d-30)  pre(nkkn) = 0.0d0

!        For minerals specified at equilibrium precip and kinetic
!        dissol move the amount calculated under equil to the amount 
!        calculated under kinetics and reset pre of equil. mineral to 0.d0
         if (MGen(m)%kineq /= 0) then
            pre(nkkn) = pre(nkkn) + pre(MGen(m)%kineq)
            pre(MGen(m)%kineq) = 0.d0
         end if 
      end do
      do m=1,nmkin
         rad (m) = MKin(m)%rad2 
         amin(m) = MKin(m)%amin2 
      end do
      
      end 

! :: ----------------------------------------------------------
! :: Initialize Secondary Species
! :: ----------------------------------------------------------
      subroutine INIT_SECONDARY(pth,dim_pth,ut,dim_ut,c,dim_c, &
                   pre , dim_pre, amin, dim_amin, rad, dim_rad,& 
                   pfug, dim_pfug)

      USE MaxChem_Dimensions
      USE Chem_Dimensions
      USE Species_Index
      USE Chem_Constants
      USE Batch_PhysiChem_Conditions
      USE Convergence_Variables
      USE CutOff_Threshold
      USE ChemOption_Variables
      USE Aqueous_Species
      USE Mineral_Phases
      USE Gaseous_Species

      integer n
      integer(8) dim_pth,dim_ut, dim_c, dim_pre
      integer(8) dim_amin, dim_rad, dim_pfug
      REAL(KIND=8)  pth(20)
      REAL(KIND=8)  ut(dim_ut), c(dim_c), pre(dim_pre)
      REAL(KIND=8)  amin(dim_amin), rad(dim_rad)
      REAL(KIND=8)  pfug(dim_pfug)
      REAL(KIND=8)  Y(dim_ut+dim_pre)


      call Pre_Chem(pth,ut,c,pre,rad,amin,pfug)

      do n = 1,npri
         Y(n) = PSpe(n)%tt
      end do

      Y(npri+1) = pre(1)

      call  AQEOUS_EQB_2ND(Y,npri+1,100)
!      call Batch_Chem_2ND(Y,npri+1,100)    
      call Post_Chem(pth,ut,c,pre,rad,amin,pfug)    

      pre(1) = Y(npri+1)

      end

! :: ----------------------------------------------------------
! :: Determine if mineral is saturated.
! :: ----------------------------------------------------------
      subroutine MINERAL_SAT(st,Y,NEQ)

      USE Chem_Dimensions
      USE Aqueous_Species
      USE Mineral_Phases
      USE Gaseous_Species

      integer NEQ
      integer(8) i,j,n,ncp
      REAL(KIND=8) Y(NEQ), st
      call  AQEOUS_EQB_2ND(Y,npri+1,100)
!      call Batch_Chem_2ND(Y,NEQ,100)

      do i=1,nmin

         st = 1.d0/(10**MGen(i)%akm)
         ncp = MGen(i)%ncpm

         do n=1,ncp
            j = MGen(i)%icpm(n)
            st = st*(PSpe(j)%cp*PSpe(j)%gamp)**MGen(i)%stqm(n)
         end do
      end do

      end 
