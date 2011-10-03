!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      subroutine Read_Chem 
!
!
! ... Modules to be used 
! 
      USE MaxChem_Dimensions
      USE Chem_Dimensions
      USE Batch_PhysiChem_Conditions
      USE Species_Index
      USE Aqueous_Species
      USE Mineral_Phases
      USE Gaseous_Species
      USE Convergence_Variables
      USE CutOff_Threshold
      USE ChemOption_Variables
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*          Read TOUGHREACT chemical system setup for                  *
!*           reactive chemical transport with XXXXX                    *
!*                  Version 1.0 - May 06, 2008                         *     
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
      INTEGER(KIND = 8)   :: i,j,isp,n,nel,nsp,nmp
      INTEGER(KIND = 8)   :: icall = 0
      CHARACTER(LEN = 39) :: dum
      CHARACTER(LEN =120) :: dum0
      CHARACTER(LEN =200) :: line
      CHARACTER(LEN =100) :: form
      SAVE       icall 
!
! -------
!
      icall = icall + 1
!
      tmpmax=300.0d0
      tmpmin=0.0d0
!
      OPEN(unit=27,file='chemsystem.dat',status='old')
      write(*,*) '   --> Reading Chemical System Data' 
      read(27,*,err=100,end=110)
      read(27,*,err=100,end=110)
      read(27,*,err=100,end=110)
      read(27,*,err=100,end=110)
      OPEN(unit=28,file='chemsystem.out',status='unknown')
      write(28,*)
      write(28,"('    Chemical System Setup:')")
      write(28,*)
      write(28,*)

      read(27,"(a39,I8)",err=100,end=110) dum,nel
      write(28,"('       Number of Initial Data Points = ',I8)") nel

      read(27,"(a39,I8)",err=100,end=110) dum,npri
      write(28,"(' Number of Component Aqueous Species = ',I8)") npri

      read(27,"(a39,I8)",err=100,end=110) dum,ntrx
      write(28,"('          Number of Aqueous Kinetics = ',I8)") ntrx

      read(27,"(a39,I8)",err=100,end=110) dum,naqx
      write(28,"('  Number of Derived Aquesous Species = ',I8)") naqx

      read(27,"(a39,I8)",err=100,end=110) dum,nmequ
      write(28,"('   Number of Minerals in Equilibrium = ',I8)") nmequ

      read(27,"(a39,I8)",err=100,end=110) dum,nmkin
      write(28,"('      Number of Minerals in Kinetics = ',I8)") nmkin

      read(27,"(a39,I8)",err=100,end=110) dum,nss
      write(28,"('           Number of Solid Solutions = ',I8)") nss

      read(27,"(a39,I8)",err=100,end=110) dum,ngas
      write(28,"('                     Number of Gases = ',I8)") ngas

      read(27,"(a39,I8)",err=100,end=110) dum,nads
      write(28,"('       Number of Surface Comlexation = ',I8)") nads

      read(27,"(a39,I8)",err=100,end=110) dum,nexc
      write(28,"('      Number of Exchangeable Cations = ',I8)") nexc

      read(27,"(a39,I8)",err=100,end=110) dum,nxsites
      write(28,"('        Number of Exchangeable Sites = ',I8)") nxsites

      read(27,"(a39,I8)",err=100,end=110) dum,nw
      write(28,"('                      Index of Water = ',I8)") nw

      read(27,"(a39,I8)",err=100,end=110) dum,nh
      write(28,"('                         Index of H+ = ',I8)") nh

      read(27,"(a39,I8)",err=100,end=110) dum,noh
      write(28,"('                        Index of OH- = ',I8)") noh

      read(27,"(a39,I8)",err=100,end=110) dum,ne
      write(28,"('                         Index of e- = ',I8)") ne

      read(27,"(a39,I8)",err=100,end=110) dum,no2aq
      write(28,"('                     Index of O2(aq) = ',I8)") no2aq

      read(27,"(a39,I8)",err=100,end=110) dum,nd
      write(28,"('                       Index of xoh- = ',I8)") nd

      read(27,"(a39,I8)",err=100,end=110) dum,nh2
      write(28,"('                     Index of H2(aq) = ',I8)") nh2

      read(27,"(a39,I8)",err=100,end=110) dum,nco2
      write(28,"('        Index of Aqueous Carbonates) = ',I8)") nco2

      read(27,"(a39,I8)",err=100,end=110) dum,nco2g
      write(28,"('                    Index of CO2 gas = ',I8)") nco2g

!      read(27,"(a39,I8)",err=100,end=110) dum,nch4
!      write(28,"('                    Index of CH4(aq) = ',I8)") nch4

!      read(27,"(a39,I8)",err=100,end=110) dum,nh2s
!      write(28,"('     Index of Aqueous Surfer Spesies = ',I8)") nh2s

      read(27,"(a39,I8)",err=100,end=110) dum,nh2
      write(28,"('                     Index of H2(aq) = ',I8)") nh2

      read(27,"(a39,I8)",err=100,end=110) dum,nx
      write(28,"('Index of Master Exchangeable Species = ',I8)") nx

      read(27,*)
      write(28,*)
      read (27,'(a39,I8,E12.4)')  dum, maxitpch,tolch   ! Maximum allowed iterations for solving geochemical system
      write(28,'(a39,I8,E12.4)')  dum, maxitpch,tolch
      read (27,'(a39,I8,E12.4)')  dum, maxitpad,tolad   ! maximum number of  iterations allowed for solving sorption via surface complexation. 
      write(28,'(a39,I8,E12.4)')  dum, maxitpad,tolad
      read (27,'(a39,E12.4)')  dum, stimax           ! maximum allowed ionic strength for HKF Model. 
      write(28,'(a39,E12.4)')  dum, stimax

      naqt=npri+naqx
      nmin=nmequ+nmkin
      ngas1=ngas

      read(27,*,err=200,end=210)
      read(27,*,err=200,end=210)
      read(27,*,err=200,end=210)
      read(27,*,err=200,end=210)


      MaxNpri     =         npri+1 ! Maximum number of primary species                                                                 
      MaxNtrx     =         ntrx+1 ! Maximum number of kinetic reactions among primary species including aqueous and sorption kinetics 
      MaxNaqx     =         naqx+1 ! Maximum number of secondary species                                                               
      MaxNmin     =         nmin+1 ! Total number of minerals (equilibrium + kinetics)                                                 
      MaxNss      =          1 ! Maximum number of solid solutions                                                                 
      MaxCpss     =          1 ! maximum number of components in each solid solution                                               
      MaxNgas     =         ngas+1 ! Maximum number of gasesous species                                                                
      MaxNads     =         nads+1 ! Maximum number of surface complexes                                                               
      MaxNexc     =         nexc+1 ! Maximum number of exchangeable cations                                                            
      MXsites     =         nxsites+1 ! Maximum number of exchangeable sites                                                              
      MaxNkdd     =          1 ! Maximum number of species with linear Kd adsorption or/and decay                                  
      MaxNiwtype  =          1 ! Maximum number of initial water types                                                             
      MaxNbwtype  =          1 ! Maximum number of boundary (injection) water types                                                
      MaxNmtype   =          1 ! Maximum number of initial mineral zones                                                           
      MaxNgtype   =          1 ! Maximum number of initial gas zones                                                               
      MaxNppzon   =          1 ! Maximum number of prosity-permeability zones                                                      
      MaxNdtype   =          1 ! Maximum number of surface complex zones                                                           
      MaxNkdtype  =          1 ! Maximum number of lineral Kd adsorption zones                                                     
      MaxNxtype   =          1 ! Maximum number of exchange zones                                                                  
      MaxNtmp     =         10 ! Maximum number of temperature points in the database                                              
      MaxNpair    =          1 ! Maximum of ion pairs for using Pitzer ion-interaction model                                       
      MaxNum_Elem =         nel                         !!!!!! Tempor..!!!! Delete leter
      MaxNum_Conx =       4*nel                     !!!!!! Tempor..!!!! Delete leter
      MaxNaqt     = MaxNpri + MaxNaqx
      MaxNnr      = MaxNpri + MaxNmin + MaxNgas + 2
! 
!      CALL ReadDimension

      CALL AllocateChemArrayMem
      
      do i=1,npri
         read(27,'(x,a20,f6.2,2x,f8.4)',err=200,end=210) &
           PSpe(i)%napri,ASpe(i)%z,ASpe(i)%a0
      enddo
!
      write(28,*)
      write(28,*) 'Primary Species:'
      write(28,*)
      write(28,'(x,a20,2a8,a10)') 'Name                ',' Charge ', &
              '    a0  '
      do i=1,npri
         write(28,'(x,a20,f6.2,2x,f8.4)') &
           PSpe(i)%napri,ASpe(i)%z,ASpe(i)%a0
      enddo
!
!
!      read(27,*)
!      read(27,*)
!      read(27,*)
!      do i=1,ntrx
!      read(27,'(4x,3I4,<max(1,ncp_rx(i))>(f8.2,I4))')
!     +     i_mod(i),n_mech(i),ncp_rx(i),
!     +     (stqrx(i,j),icprx(i,j),j=1,ncp_rx(i))
!      do,j=1,n_mech(i)
!      read(27,'(e12.4,I3,<max(1,ncp_rx1(i,j))>(I3,f8.2,I3),
!     +           e12.4,I3,<max(1,ncp_rx2(i,j))>(I3,f8.2,I3),
!     +           e12.4,I3,<max(1,ncp_rx3(i,j))>(I3,f8.2,I3))')
!     +  rkaq(i,j),ncp_rx1(i,j),(ia1(i,j,n),s_rx1(i,j,n),
!     +  icprx1(i,j,n),n=1,ncp_rx1(i,j)),       
!     +  rkaq(i,j),ncp_rx2(i,j),(ia2(i,j,n),s_rx2(i,j,n),
!     +  icprx2(i,j,n),n=1,ncp_rx2(i,j)),       
!     +  rkaq(i,j),ncp_rx3(i,j),(ia3(i,j,n),s_rx3(i,j,n),
!     +  icprx3(i,j,n),n=1,ncp_rx3(i,j))
!      enddo       
!      if(i_mod(i)==2)
!     +  read(27,'(5e16.8)') (coef_rx(i,n),n=1,5)  
!      enddo
!
!      write(28,*)
!      if(ntrx.eq.0) then
!	write(28,'(a24)') 'Aqueous Kinetics:   None'
!      else
!      write(28,'(a17)') 'Aqueous Kinetics:'
!      endif
!      write(28,*)
!      do i=1,ntrx
!      write(28,'(4I4,<max(1,ncp_rx(i))>(f8.2,I4))')
!     +     i,i_mod(i),n_mech(i),ncp_rx(i),
!     +     (stqrx(i,j),icprx(i,j),j=1,ncp_rx(i))
!      do,j=1,n_mech(i)
!      write(28,'(e12.4,I3,<max(1,ncp_rx1(i,j))>(I3,f8.2,I3),
!     +           e12.4,I3,<max(1,ncp_rx2(i,j))>(I3,f8.2,I3),
!     +           e12.4,I3,<max(1,ncp_rx3(i,j))>(I3,f8.2,I3))')
!     +  rkaq(i,j),ncp_rx1(i,j),(ia1(i,j,n),s_rx1(i,j,n),
!     +  icprx1(i,j,n),n=1,ncp_rx1(i,j)),       
!     +  rkaq(i,j),ncp_rx2(i,j),(ia2(i,j,n),s_rx2(i,j,n),
!     +  icprx2(i,j,n),n=1,ncp_rx2(i,j)),       
!     +  rkaq(i,j),ncp_rx3(i,j),(ia3(i,j,n),s_rx3(i,j,n),
!     +  icprx3(i,j,n),n=1,ncp_rx3(i,j))
!      enddo       
!      if(i_mod(i)==2)
!     +  write(28,'(5e16.8)') (coef_rx(i,n),n=1,5)  
!      enddo
!
      read(27,*,err=300,end=310)
      read(27,*,err=300,end=310)
      read(27,*,err=300,end=310)
      do i=1,naqx
          ii=i+npri
      read(27,'(a)') line
      read(line(38:40),'(I3)') nsp
      write(form(1:42),'(a23,I2,a17)') '(x,a20,f6.2,2x,f8.4,I3,', &
                                 nsp,'(f6.2,I3),5e14.6)' 
!      read(27,'(x,a20,f6.2,2x,f8.4,I3,<SSpe(i)%ncps>(f6.2,I3),5e14.6)',
      read(line,form,err=300,end=310) &
           SSpe(i)%naaqx,ASpe(ii)%z,ASpe(ii)%a0,SSpe(i)%ncps, &
           (SSpe(i)%stqs(n),SSpe(i)%icps(n),n=1,SSpe(i)%ncps), &
           (SSpe(i)%akcoes(n),n=1,5)  
      enddo
!
      write(28,*)
      write(28,*) 'Secondary Species:'
      write(28,*)
      do i=1,naqx
          ii=i+npri
      write(form(1:42),'(a23,I2,a17)') '(x,a20,f6.2,2x,f8.4,I3,', &
            SSpe(i)%ncps,'(f6.2,I3),5e14.6)' 
!      write(28,'(x,a20,f6.2,2x,f8.4,I3,<SSpe(i)%ncps>(f6.2,I3),5e14.6)') 
      write(28,form) &
            SSpe(i)%naaqx,ASpe(ii)%z,ASpe(ii)%a0,SSpe(i)%ncps,&
           (SSpe(i)%stqs(n),SSpe(i)%icps(n),n=1,SSpe(i)%ncps),&
           (SSpe(i)%akcoes(n),n=1,5)
      enddo
!
      if(nmin.ne.0) then
      read(27,*,err=400,end=410)
      read(27,*,err=400,end=410)
      read(27,*,err=400,end=410)
      do i=1,nmin
          ii=i+npri+naqx
      read(27,'(a)') line
      read(line(46:48),'(I3)') nsp
      write(form(1:46),'(a27,I2,a17)') '(x,a20,3I2,f8.2,2x,f8.4,I3,',&
                                nsp,'(f6.2,I3),5e14.6)' 
!      read(27,'(x,a20,3I2,f8.2,2x,f8.4,I3,<MGen(i)%ncpm>(f6.2,I3),
!     +                            5e14.6)
      read(line,form,err=400,end=410) &
            MGen(i)%namin,MGen(i)%ikin,MGen(i)%idispre,MGen(i)%iss,&
            MGen(i)%dmolwm,MGen(i)%vmin,MGen(i)%ncpm,&
           (MGen(i)%stqm(n),MGen(i)%icpm(n),n=1,MGen(i)%ncpm),&
           (MGen(i)%akcoem(n),n=1,5)
      enddo
!
      write(28,*)
      write(28,*) 'Mineral Thermodynamics:'
      write(28,*)
      do i=1,nmin
          ii=i+npri+naqx
      write(form(1:46),'(a27,I2,a17)') '(x,a20,3I2,f8.2,2x,f8.4,I3,',&
                                   MGen(i)%ncpm,'(f6.2,I3),5e14.6)' 
!      write(28,'(x,a20,3I2,f8.2,2x,f8.4,I3,<MGen(i)%ncpm>(f6.2,I3),
!     +                            5e14.6)
      write(28,form)&
            MGen(i)%namin,MGen(i)%ikin,MGen(i)%idispre,MGen(i)%iss,&
            MGen(i)%dmolwm,MGen(i)%vmin,MGen(i)%ncpm,&
           (MGen(i)%stqm(n),MGen(i)%icpm(n),n=1,MGen(i)%ncpm),&
           (MGen(i)%akcoem(n),n=1,5)
      enddo
      endif
!
      if(nmkin.ne.0) then
      read(27,*,err=400,end=410)
      read(27,*,err=400,end=410)
      read(27,*,err=400,end=410)
      do i=1,nmkin
       ii=i+nmequ

      read(27,'(a)') line
      read(line(141:143),'(I3)') nsp
      read(line(max(1,nsp)*19+144:max(1,nsp)*19+146),'(I3)') nmp

      write(form(1:72),'(a40,I2,a15,I2,a13)') &
          "(x,a20,I2,e12.4,I2,2f6.2,f7.2,7e12.4,I3,",&
                         max(1,nsp),"(e12.4,f7.2,I3,",& 
                         max(1,nmp),"(I3,f7.2)))')" 
!       read(27,'(x,a20,I2,e12.4,I2,2f6.2,f7.2,7e12.4,I3,
!     + <max(1,MKin(i)%ndis)>(e12.4,f7.2,I3,
!     + <max(1,MKin(i)%nspds(j))>(I3,f7.2)))') 


       read(line,form,err=400,end=420) &
       MGen(ii)%namin,MGen(i)%kineq,&
       MKin(i)%rkf,MKin(i)%idep,MKin(i)%ck1,MKin(i)%ck2,&
       MKin(i)%ea,MKin(i)%acfdiss,MKin(i)%bcfdiss,MKin(i)%ccfdiss,&
       MKin(i)%aH1,MKin(i)%aH2,MKin(i)%aHexp,MKin(i)%aOHexp,&
       MKin(i)%ndis,(MKin(i)%rkds(j),MKin(i)%eads(j),MKin(i)%nspds(j),&
      (MKin(i)%ids(j,isp),MKin(i)%expds(j,isp),isp=1,MKin(i)%nspds(j)),&
       j=1,MKin(i)%ndis)

      read(27,'(a)') line
      read(line(141:143),'(I3)') nsp
      read(line(max(1,nsp)*19+144:max(1,nsp)*19+146),'(I3)') nmp
      write(form(1:67),'(a35,I2,a15,I2,a13)') &
          "(23x,e12.4,I2,2f6.2,f7.2,7e12.4,I3,",&
                         max(1,nsp),"(e12.4,f7.2,I3,", &
                         max(1,nmp),"(I3,f7.2)))')" 
!       read(27,'(23x,e12.4,I2,2f6.2,f7.2,7e12.4,I3,
!     + <max(1,MKin(i)%npre)>(e12.4,f7.2,I3,
!     + <max(1,MKin(i)%nsppr(j))>(I3,f7.2)))')
       read(line,form,err=400,end=420) &
            MKin(i)%rkprec,MKin(i)%ideprec,MKin(i)%ck1prec,MKin(i)%ck2prec,&
            MKin(i)%eaprec,MKin(i)%acfprec,MKin(i)%bcfprec,MKin(i)%ccfprec,&
            MKin(i)%aH1p,MKin(i)%aH2p,MKin(i)%aOHexpp,MKin(i)%aOHexpp,     &      
            MKin(i)%npre,(MKin(i)%rkdsp(j),MKin(i)%eadsp(j),MKin(i)%nsppr(j),&
            (MKin(i)%idsp(j,isp),MKin(i)%expdsp(j,isp),                    &
            isp=1,MKin(i)%nsppr(j)),j=1,MKin(i)%npre)                      

       read(27,'(23x,4e12.4,i5)') MKin(i)%ssqk0, &
            MKin(i)%sstk1,MKin(i)%sstk2,MKin(i)%rnucl,MKin(i)%nplaw
      enddo
!
      write(28,*)
      write(28,*) 'Mineral Kinetics:'
      write(28,*)
      do i=1,nmkin
       ii=i+nmequ

!     I changed something here. In last line, MKin(i)%ndis was MKin(i)%nspds(j)
      write(form(1:72),'(a40,I2,a15,I2,a13)') &
          "(x,a20,I2,e12.4,I2,2f6.2,f7.2,7e12.4,I3,", &
               max(1,MKin(i)%ndis),"(e12.4,f7.2,I3,", &
               max(1,MKin(i)%ndis),"(I3,f7.2)))')" 

!       write(28,'(x,a20,I2,e12.4,I2,2f6.2,f7.2,7e12.4,I3,
!     + <max(1,MKin(i)%ndis)>(e12.4,f7.2,I3,
!     + <max(1,MKin(i)%nspds(j))>(I3,f7.2)))') 
       write(28,form)&
       MGen(ii)%namin,MGen(i)%kineq,&
       MKin(i)%rkf,MKin(i)%idep,MKin(i)%ck1,MKin(i)%ck2,&
       MKin(i)%ea,MKin(i)%acfdiss,MKin(i)%bcfdiss,MKin(i)%ccfdiss,&
       MKin(i)%aH1,MKin(i)%aH2,MKin(i)%aHexp,MKin(i)%aOHexp,&
       MKin(i)%ndis,(MKin(i)%rkds(j),MKin(i)%eads(j),MKin(i)%nspds(j),&
      (MKin(i)%ids(j,isp),MKin(i)%expds(j,isp),isp=1,MKin(i)%nspds(j)),&
       j=1,MKin(i)%ndis)

!     I changed: MKin(i)%nsppr(j) to ,MKin(i)%npre
      write(form(1:67),'(a35,I2,a15,I2,a13)') &
          "(23x,e12.4,I2,2f6.2,f7.2,7e12.4,I3,",&
              max(1,MKin(i)%npre),"(e12.4,f7.2,I3,", &
              max(1,MKin(i)%npre),"(I3,f7.2)))')" 

!       write(28,'(23x,e12.4,I2,2f6.2,f7.2,7e12.4,I3,
!     + <max(1,MKin(i)%npre)>(e12.4,f7.2,I3,
!     + <max(1,MKin(i)%nsppr(j))>(I3,f7.2)))')
       write(28,form)&
       MKin(i)%rkprec,MKin(i)%ideprec,MKin(i)%ck1prec,MKin(i)%ck2prec,&
       MKin(i)%eaprec,MKin(i)%acfprec,MKin(i)%bcfprec,MKin(i)%ccfprec,&
       MKin(i)%aH1p,MKin(i)%aH2p,MKin(i)%aOHexpp,MKin(i)%aOHexpp,      &      
       MKin(i)%npre,(MKin(i)%rkdsp(j),MKin(i)%eadsp(j),MKin(i)%nsppr(j),&
       (MKin(i)%idsp(j,isp),MKin(i)%expdsp(j,isp),&
       isp=1,MKin(i)%nsppr(j)),j=1,MKin(i)%npre)                      

       write(28,'(23x,4e12.4,i5)') MKin(i)%ssqk0,&
       MKin(i)%sstk1,MKin(i)%sstk2, MKin(i)%rnucl,MKin(i)%nplaw
      enddo
      endif
!
      if(nss.ne.0) then
      read(27,*)
      read(27,*)
      read(27,*)
      do i=1,nss
       read(27,'(a)',err=400,end=430) line
       read(line(1:3),'(I3)',err=400,end=430) nsp
       write(form(1:8),'(a4,I2,a3)') "(i3,",nsp,"I4)"
       read(line,form,err=400,end=430) &
            ncpss(i),(icpss(i,j),j=1,ncpss(i))
      enddo
!
      write(28,*)
      write(28,'(a16,I7)') 'Solid Solutions:', nss
      write(28,*)
      do i=1,nss
       write(form(1:8),'(a4,I2,a3)') "(i3,",ncpss(i),"I4)"
       write(28,form) ncpss(i),(icpss(i,j),j=1,ncpss(i))
      enddo
      endif
!
      if(ngas.ne.0) then
      read(27,*,err=500,end=510)
      read(27,*,err=500,end=510)
      read(27,*,err=500,end=510)
      do i=1,ngas

      read(27,'(a)') line
      read(line(38:40),'(I3)') nsp
      write(form(1:42),'(a23,I2,a17)') "(x,a20,f6.2,2x,f8.4,I3,", &
                                nsp,"(f6.2,I3),5e14.6)" 

!      read(27,'(x,a20,f6.2,2x,f8.4,I3,<GSpe(i)%ncpg>(f6.2,I3),5e14.6)'
!     +         ,err=500,end=510) 
      read(line,form,err=500,end=510) &
            GSpe(i)%nagas,GSpe(i)%dmwgas,GSpe(i)%diamol,GSpe(i)%ncpg,&
           (GSpe(i)%stqg(n),GSpe(i)%icpg(n),n=1,GSpe(i)%ncpg),&
                       (GSpe(i)%akcoeg(n),n=1,5)  
      enddo
!
      write(28,*)
      write(28,*) 'Gases:'
      write(28,*)

      do i=1,ngas
      write(form(1:42),'(a23,I2,a17)') "(x,a20,f6.2,2x,f8.4,I3,", &
                                GSpe(i)%ncpg,"(f6.2,I3),5e14.6)" 
!      write(28,'(x,a20,f6.2,2x,f8.4,I3,<GSpe(i)%ncpg>(f6.2,I3),5e14.6)') 
      write(28,form) &
            GSpe(i)%nagas,GSpe(i)%dmwgas,GSpe(i)%diamol,GSpe(i)%ncpg,&
           (GSpe(i)%stqg(n),GSpe(i)%icpg(n),n=1,GSpe(i)%ncpg),&
                       (GSpe(i)%akcoeg(n),n=1,5)   
      enddo
      endif

      close(27)
      close(28)

!
      ntot = npri + nmin + ngas    
!
! ----------
!.....Aqueous species name index
! ----------
!  
      do i=1,npri
         ASpe(i)%naaqt = PSpe(i)%napri
      end do
!
      do j=1,naqx
         ASpe(npri+j)%naaqt = SSpe(j)%naaqx
      end do
!
      do i=1,naqt
         ASpe(i)%zz2 = ASpe(i)%z*ASpe(i)%z
      end do                                  
!
! ----------
!.....Derivative of activity coefficients (gam) is not used for DH model
!.....is used for pitzer model. Set zero here, when Pitzer model is used 
!.....this arrays will be overwritten
! ----------
!  
      do j=1,npri
         do i=1,npri                              
            PSpe(j)%dgamp(i) = 0.0d0
            ASpe(j)%dgamt(i) = 0.0d0
         end do
      end do
!
      do j=1,naqx
         do i=1,npri 
            SSpe(j)%dgams(i)      = 0.0d0
            ASpe(npri+j)%dgamt(i) = 0.0d0
         enddo
      end do

      return

 100  write(*,*)
      write(*,*) " Error found in file 'chemsystem.dat' while reading"
      write(*,*) "   chemical dimenssion and index varibles, stopped!"
       stop
 110  write(*,*)
      write(*,*) " File 'chemsystem.dat' ended, while reading"
      write(*,*) "   chemical dimenssion and index varibles, stopped!"
       stop
 200  write(*,*)
      write(*,'(a90)') " Error found in file 'chemsystem.dat', &
       while reading primary species, stopped!"
       stop
 210  write(*,*)
      write(*,'(a90)') " File 'chemsystem.dat' ended, &
      while reading primary species, stopped!"
       stop
 300  write(*,*)
      write(*,'(a90)') " Error found in file 'chemsystem.dat', &
      while reading secondary species, stopped!"
       stop
 310  write(*,*)
      write(*,'(a90)') " File 'chemsystem.dat' ended, &
      while reading secondary species, stopped!"
       stop
 400  write(*,*)
      write(*,'(a90)') " Error found in file 'chemsystem.dat', &
      while reading minerlas, stopped!"
       stop
 410  write(*,*)
      write(*,'(a90)') " File 'chemsystem.dat' ended, &
      while reading mineral thermodynamics, stopped!"
       stop
 420  write(*,*)
      write(*,'(a90)') " File 'chemsystem.dat' ended, &
      while reading mineral kinetics, stopped!"
       stop
 430  write(*,*)
      write(*,'(a90)') " File 'chemsystem.dat' ended, &
      while reading solid solution, stopped!"
       stop
 500  write(*,*)
      write(*,'(a90)') " Error found in file 'chemsystem.dat', &
      while reading gases, stopped!"
       stop
 510  write(*,*)
      write(*,'(a90)') " File 'chemsystem.dat' ended, &
      while reading gases, stopped!"
       stop

      end subroutine Read_Chem 
