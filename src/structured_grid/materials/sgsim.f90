!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                      %
! Copyright (C) 2003, Statios Software and Services Incorporated.  All %
! rights reserved.                                                     %
!                                                                      %
! This program has been modified from the one distributed in 1996 (see %
! below).  This version is also distributed in the hope that it will   %
! be useful, but WITHOUT ANY WARRANTY. Compiled programs based on this %
! code may be redistributed without restriction; however, this code is %
! for one developer only. Each developer or user of this source code   %
! must purchase a separate copy from Statios.                          %
!                                                                      %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                      %
! Copyright (C) 1996, The Board of Trustees of the Leland Stanford     %
! Junior University.  All rights reserved.                             %
!                                                                      %
! The programs in GSLIB are distributed in the hope that they will be  %
! useful, but WITHOUT ANY WARRANTY.  No author or distributor accepts  %
! responsibility to anyone for the consequences of using them or for   %
! whether they serve any particular purpose or work at all, unless he  %
! says so in writing.  Everyone is granted permission to copy, modify  %
! and redistribute the programs in GSLIB, but only under the condition %
! that this notice and the above copyright notice remain intact.       %
!                                                                      %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine readparm(paramfl,nx_in,xsiz_in,xmn_in,rseed)
! -----------------------------------------------------------------------
!
!                  Initialization and Read Parameters
!                  **********************************
!
! The input parameters and data are read in from their files. Some quick
! error checking is performed and the statistics of all the variables
! being considered are written to standard output.
!
!
!
! -----------------------------------------------------------------------

      use       geostat
      include  'sgsim.inc'

      integer   rseed
      integer   MAXSB
      real*8    var(50)
      real*8    p,acorni,cp,oldcp,w
      character transfl*512,smthfl*512,tmpfl*512
      character datafl*512,outfl*512
      character dbgfl*512,lvmfl*512,str*512
!      character paramfl*paramsz
      character paramfl*512
      logical   testfl,trans
      
      integer nx_in(3)
      real*8  xsiz_in(3), xmn_in(3)

      integer do_cond
      integer ierr,iend,istart,nt,icolwt,icolvr,j,nvari,idum,i
      integer MAXX,MAXY,MAXZ,MXYZ,MAXTMP,MAXDAT,MAXKR2,MAXSAM,MAXXYZ
      integer iswt,isvr,ismooth,isecvr
      integer iwt,ixl,iyl,izl,ivrl
      real*8  powint,vrr,vrg,twt,av,ss
      real*8  aa1,aa2,xx,radius1,radius2
      real*8  tmin,tmax,sill

      call bl_abort()
!
! Input/Output units used:
!
      lin  = 1
      lout = 2
      ldbg = 3
      llvm = 4

      open(lin,file=paramfl,status='OLD')
!
! Find Start of Parameters:
!
 1    read(lin,'(a4)',end=98) str(1:4)
      if(str(1:4).ne.'STAR') go to 1
!
! Read Input Parameters:
!
!      read(lin,*,err=98) do_cond
!      write(*,*) ' do_cond = ',do_cond
      do_cond = 1

      read(lin,'(a512)',err=98) datafl
      call chknam(datafl,512)
      write(*,*) ' data file = ',datafl(1:40)

      read(lin,*,err=98) ixl,iyl,izl,ivrl,iwt,isecvr
      write(*,*) ' input columns = ',ixl,iyl,izl,ivrl,iwt,isecvr

      read(lin,*,err=98) tmin,tmax
      write(*,*) ' trimming limits = ',tmin,tmax

      read(lin,*,err=98) itrans
      write(*,*) ' transformation flag = ',itrans

      read(lin,'(a512)',err=98) transfl
      call chknam(transfl,512)
      write(*,*) ' transformation file = ',transfl(1:40)

      read(lin,*,err=98) ismooth
      write(*,*) ' consider smoothed distribution (1=yes) = ',ismooth

      read(lin,'(a512)',err=98) smthfl
      call chknam(smthfl,512)
      write(*,*) ' file with smoothed distribution = ',smthfl(1:40)

      read(lin,*,err=98) isvr,iswt
      write(*,*) ' columns = ',isvr,iswt

      read(lin,*,err=98) zmin,zmax
      write(*,*) ' data limits (tails) = ',zmin,zmax

      read(lin,*,err=98) ltail,ltpar
      write(*,*) ' lower tail = ',ltail,ltpar

      read(lin,*,err=98) utail,utpar
      write(*,*) ' upper tail = ',utail,utpar

      read(lin,*,err=98) idbg
      write(*,*) ' debugging level = ',idbg

      read(lin,'(a512)',err=98) dbgfl
      call chknam(dbgfl,512)
      write(*,*) ' debugging file = ',dbgfl(1:40)
      open(ldbg,file=dbgfl,status='UNKNOWN')

      read(lin,'(a512)',err=98) outfl
      call chknam(outfl,512)
      write(*,*) ' output file ',outfl(1:40)

      read(lin,*,err=98) nsim
      write(*,*) ' number of realizations = ',nsim

      read(lin,*,err=98) nx,xmn,xsiz
      write(*,*) ' X grid specification = ',nx,xmn,xsiz

      read(lin,*,err=98) ny,ymn,ysiz
      write(*,*) ' Y grid specification = ',ny,ymn,ysiz

      read(lin,*,err=98) nz,zmn,zsiz
      write(*,*) ' Z grid specification = ',nz,zmn,zsiz
 
      nx   = nx_in(1)
      ny   = nx_in(2)
      nz   = nx_in(3)
      xmn  = xmn_in(1)
      ymn  = xmn_in(2)
      zmn  = xmn_in(3)
      xsiz = xsiz_in(1)
      ysiz = xsiz_in(2)
      zsiz = xsiz_in(3)

      nxy  = nx*ny
      nxyz = nx*ny*nz

      read(lin,*,err=98) ixv(1)
      rseed = ixv(1)
      write(*,*) ' random number seed = ',ixv(1)

      call bl_pd_myproc(myprocid)
      call blutilinitrand(rseed + myprocid)

!      do i=1,1000
!         p = acorni(idum)
!      end do

      read(lin,*,err=98) ndmin,ndmax
      write(*,*) ' min and max data = ',ndmin,ndmax

      read(lin,*,err=98) nodmax
      write(*,*) ' maximum previous nodes = ',nodmax

      read(lin,*,err=98) sstrat
      write(*,*) ' two-part search flag = ',sstrat
      if(sstrat.eq.1) ndmax = 0

      read(lin,*,err=98) mults,nmult
      write(*,*) ' multiple grid search flag = ',mults,nmult

      read(lin,*,err=98) noct
      write(*,*) ' number of octants = ',noct

      read(lin,*,err=98) radius,radius1,radius2
      write(*,*) ' search radii = ',radius,radius1,radius2
      if(radius.lt.EPSLON) stop 'radius must be greater than zero'
      radsqd = radius  * radius
      sanis1 = radius1 / radius
      sanis2 = radius2 / radius

      read(lin,*,err=98) sang1,sang2,sang3
      write(*,*) ' search anisotropy angles = ',sang1,sang2,sang3

      read(lin,*,err=98) mxctx,mxcty,mxctz
      write(*,*) ' size of covariance lookup = ',mxctx,mxcty,mxctz
      
      read(lin,*,err=98) ktype
      write(*,*) ' kriging type = ',ktype
      
      trans = .true.
      if(ktype.lt.0) then
            trans = .false.
            ktype = abs(ktype)
      end if

      colocorr = 0.d0
      if(ktype.eq.4) then
            backspace lin
            read(lin,*,err=98) i,colocorr
            varred = 1.d0
            backspace lin
            read(lin,*,err=9990) i,xx,varred
 9990       continue
            write(*,*) ' correlation coefficient = ',colocorr
            write(*,*) ' secondary variable varred = ',varred
      end if

      read(lin,'(a512)',err=98) lvmfl
      call chknam(lvmfl,512)
      write(*,*) ' secondary model file = ',lvmfl(1:40)

      read(lin,*,err=98) icollvm
      write(*,*) ' column in secondary model file = ',icollvm

      read(lin,*,err=98) nst(1),c0(1)
      sill = c0(1)
      write(*,*) ' nst, c0 = ',nst(1),c0(1)

      if(nst(1).le.0) then
            write(*,9997) nst(1)
 9997       format(' nst must be at least 1, it has been set to ',i4,/, &
                   ' The c or a values can be set to zero')
            stop
      endif

      do i=1,nst(1)
            read(lin,*,err=98) it(i),cc(i),ang1(i),ang2(i),ang3(i)
            read(lin,*,err=98) aa(i),aa1,aa2
            anis1(i) = aa1 / max(aa(i),EPSLON)
            anis2(i) = aa2 / max(aa(i),EPSLON)
            sill     = sill + cc(i)
            if(it(i).eq.4) then
                  write(*,*) ' A power model is NOT allowed '
                  write(*,*) ' Choose a different model and re start '
                  stop
            endif
            write(*,*) ' it,cc,ang[1,2,3]; ',it(i),cc(i),&
                         ang1(i),ang2(i),ang3(i)
            write(*,*) ' a1 a2 a3: ',aa(i),aa1,aa2
      end do
      write(*,*)
      close(lin)
!
! Find the needed parameters:
!
      MAXCTX = mxctx
      MAXCTY = mxcty
      MAXCTZ = mxctz
      MAXCXY = MAXCTX * MAXCTY
      MAXXYZ = MAXCTX * MAXCTY * MAXCTZ
      MAXX   = nx
      MAXY   = ny
      MAXZ   = nz
      MXYZ   = MAXX * MAXY * MAXZ
      if(MXYZ.lt.100) MXYZ = 100
      MAXNOD = nodmax
      MAXSAM = ndmax
      MAXKR1 = MAXNOD + MAXSAM + 1
      MAXKR2 = MAXKR1 * MAXKR1
      MAXSBX = 1
      if(nx.gt.1)then
            MAXSBX = int(nx/2)
            if(MAXSBX.gt.50)MAXSBX=50
      end if

      MAXSBY = 1
      if(ny.gt.1)then
            MAXSBY = int(ny/2)
            if(MAXSBY.gt.50)MAXSBY=50
      end if

      MAXSBZ = 1
      if(nz.gt.1)then
            MAXSBZ = int(nz/2)
            if(MAXSBZ.gt.50)MAXSBZ=50
      end if
      MAXSB = MAXSBX*MAXSBY*MAXSBZ

!
! Find MAXDAT:
!
      MAXDAT = 100
      inquire(file=datafl,exist=testfl)
      if(testfl)then
            open(lin,file=datafl,status='UNKNOWN')
            read(lin,*,err=98)
            read(lin,*,err=99) nvari
            do i=1,nvari
                  read(lin,*)
            end do
            MAXDAT = 0
 33         read(lin,*,end=66,err=98)(var(j),j=1,nvari)
            MAXDAT = MAXDAT + 1
            go to 33
 66         continue
            rewind(lin)
            close(lin)
      end if

      MAXTMP = 1
      inquire(file=smthfl,exist=testfl)
      if(testfl)then
            open(lin,file=smthfl,status='UNKNOWN')
            read(lin,*,err=98)
            read(lin,*,err=97) nvari
            do i=1,nvari
                  read(lin,*)
            end do
            MAXTMP = 0
 22         read(lin,*,end=55,err=97)(var(j),j=1,nvari)
            MAXTMP = MAXTMP + 1
            go to 22
 55         continue
            rewind(lin)
            close(lin)
      end if
      if(MAXTMP.gt.MAXDAT)MAXDAT = MAXTMP

!
! Allocate the needed memory:
!
      allocate(x(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 1: Allocation failed',&
                       ' due to insufficient memory.'
                  stop
            end if

      allocate(y(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 2: Allocation failed',&
                       ' due to insufficient memory.'
                  stop
            end if

      allocate(z(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 3: Allocation failed',&
                       ' due to insufficient memory.'
                  stop
            end if

      allocate(vr(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 4: Allocation failed',&
                       ' due to insufficient memory.'
                  stop
            end if

      allocate(wt(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 5: Allocation failed',&
                       ' due to insufficient memory.'
                  stop
            end if

      allocate(vrtr(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 6: Allocation failed',&
                       ' due to insufficient memory.'
                  stop
            end if

      allocate(vrgtr(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 7: Allocation failed',&
                       ' due to insufficient memory.'
                  stop
            end if

      allocate(close(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 8: Allocation failed',&
                       ' due to insufficient memory.'
                  stop
            end if

      allocate(sec(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 9: Allocation failed',&
                       ' due to insufficient memory.'
                  stop
            end if

      allocate(covtab(MAXCTX,MAXCTY,MAXCTZ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 14: Allocation failed',&
                       ' due to insufficient memory.'
                  stop
            end if

      allocate(cnodex(MAXNOD),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 15: Allocation failed',&
                       ' due to insufficient memory.'
                        stop
            end if

      allocate(cnodey(MAXNOD),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 16: Allocation failed',&
                       ' due to insufficient memory.'
                  stop
            end if

      allocate(cnodez(MAXNOD),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 17: Allocation failed',&
                       ' due to insufficient memory.'
                  stop
            end if

      allocate(cnodev(MAXNOD),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 18: Allocation failed',&
                       ' due to insufficient memory.'
                        stop
            end if

      allocate(r(MAXKR1),stat = test)

            if(test.ne.0)then
                  write(*,*)'ERROR 21: Allocation failed',&
                       ' due to insufficient memory.'
                  stop
            end if

      allocate(rr(MAXKR1),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 22: Allocation failed',&
                       ' due to insufficient memory.'
                  stop
            end if

      allocate(s(MAXKR1),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 23: Allocation failed',&
                       ' due to insufficient memory.'
                  stop
            end if

      allocate(a(MAXKR2),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 24: Allocation failed',&
                       ' due to insufficient memory.'
                  stop
            end if

      allocate(nisb(MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 25: Allocation failed',&
                       ' due to insufficient memory.'
                  stop
            end if

      allocate(icnode(MAXNOD),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 26: Allocation failed',&
                       ' due to insufficient memory.'
                  stop
            end if

      allocate(ixnode(MAXXYZ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 27: Allocation failed',&
                       ' due to insufficient memory.'
                  stop
            end if

      allocate(iynode(MAXXYZ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 28: Allocation failed',&
                       ' due to insufficient memory.'
                  stop
            end if

      allocate(iznode(MAXXYZ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 29: Allocation failed',&
                       ' due to insufficient memory.'
                  stop
            end if

      allocate(ixsbtosr(8*MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 30: Allocation failed',&
                       ' due to insufficient memory.'
                  stop
            end if

      allocate(iysbtosr(8*MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 31: Allocation failed',&
                       ' due to insufficient memory.'
                  stop
            end if

      allocate(izsbtosr(8*MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 32: Allocation failed',&
                       ' due to insufficient memory.'
                  stop
            end if

!
! Warn the user if the sill is different than 1.0:
!
      if(sill.gt.(1.d0+EPSLON).or.sill.lt.(1.d0-EPSLON)) then
            write(*,*) 'WARNING the sill of your variogram is not 1.0!'
            write(*,*) '        the sill = ',sill
            write(*,*)
      end if
!
! Perform some quick error checking:
!
      testfl = .false.
      if(nx.gt.MAXX.or.ny.gt.MAXY.or.nz.gt.MAXZ) then
            write(*,*) 'ERROR: available grid size: ',MAXX,MAXY,MAXZ
            write(*,*) '       you have asked for : ',nx,ny,nz
            testfl = .true.
      end if
      if(ltail.ne.1.and.ltail.ne.2) then
            write(*,*) 'ERROR invalid lower tail option ',ltail
            write(*,*) '      only allow 1 or 2 - see manual '
            testfl = .true.
      endif
      if(utail.ne.1.and.utail.ne.2.and.utail.ne.4) then
            write(*,*) 'ERROR invalid upper tail option ',ltail
            write(*,*) '      only allow 1,2 or 4 - see manual '
            testfl = .true.
      endif
      if(utail.eq.4.and.utpar.lt.1.d0) then
            write(*,*) 'ERROR invalid power for hyperbolic tail',utpar
            write(*,*) '      must be greater than 1.0!'
            testfl = .true.
      endif
      if(ltail.eq.2.and.ltpar.lt.0.d0) then
            write(*,*) 'ERROR invalid power for power model',ltpar
            write(*,*) '      must be greater than 0.0!'
            testfl = .true.
      endif
      if(utail.eq.2.and.utpar.lt.0.d0) then
            write(*,*) 'ERROR invalid power for power model',utpar
            write(*,*) '      must be greater than 0.0!'
            testfl = .true.
      endif
      if(testfl) stop
!
! Check to make sure the data file exists:
!
      nd = 0
      av = 0.d0
      ss = 0.d0
      inquire(file=datafl,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'WARNING data file ',datafl,' does not exist!'
            write(*,*) '   - Hope your intention was to create an ',&
                            'unconditional simulation'
            write(*,*) '   - Resetting ndmin, ndmax, and itrans  to 0 '
            write(*,*) '   - Resetting sstrat to 1 '
            ndmin  = 0
            ndmax  = 0
            sstrat = 1
      end if
!
! Establish the reference histogram for the simulation (provided that
! we have data, and we are transforming the data):
!
      if(itrans.eq.1) then
            write(*,*) 'Setting up transformation table'
!
! Decide which file to use for establishing the transformation table:
!
            if(ismooth.eq.1) then
                  tmpfl  = smthfl
                  icolvr = isvr
                  icolwt = iswt
            else
                  tmpfl  = datafl
                  icolvr = ivrl
                  icolwt = iwt
            end if
            inquire(file=tmpfl,exist=testfl)
            if(.not.testfl) then
                  write(*,*) 'ERROR: ',tmpfl,' does not exist'
                  write(*,*) '       this file is needed! '
                  stop
            endif
!
! Open up the file with reference distribution:
!
            open(lin,file=tmpfl,status='UNKNOWN')
            read(lin,'(a40)',err=98) str(1:40)
            read(lin,*,err=99) nvari
            do i=1,nvari
                  read(lin,*,err=98)
            end do
!
! Now, read in the actual data:
!
            nt     = 0
            ntr    = 0
            twt    = 0.d0
 3          read(lin,*,end=4,err=99) (var(j),j=1,nvari)
!
! Trim this data?
!
            if(var(icolvr).lt.tmin.or.var(icolvr).ge.tmax) then
                  nt = nt + 1
                  go to 3
            endif
            ntr = ntr + 1
!
! Exceeded available storage?
!
            if(icolvr.gt.nvari.or.icolwt.gt.nvari) then
                  write(*,*) ' ERROR: too few columns in ref data '
                  stop
            endif
!
! Keep this data: Assign the data value and coordinate location:
!
            vrtr(ntr) = var(icolvr)
            if(icolwt.le.0) then
                  vrgtr(ntr) = 1.d0
            else
                  vrgtr(ntr) = var(icolwt)
            endif
            if(vrgtr(ntr).le.0.d0) then
                  ntr = ntr - 1
                  nt  = nt  + 1
                  go to 3
            end if
            twt = twt + vrgtr(ntr)
!
! Go back for another datum:
!
            go to 3
 4          close(lin)
            if(ntr.le.1) then
                  write(*,*) 'ERROR: too few data for transformation'
                  stop
            endif
!
! Write transformation table:
!
            open(lout,file=transfl,status='UNKNOWN')
!
! Sort data by value:
!
            istart = 1
            iend   = ntr
!            call sortem(istart,iend,vrtr,1,vrgtr,c,d,e,f,g,h)
            call dsortem(istart,iend,vrtr,1,vrgtr)
!
! Compute the cumulative probabilities and write transformation table
!
            twt   = max(twt,EPSLON)
            oldcp = 0.d0
            cp    = 0.d0
            do j=istart,iend
                  cp =  cp + dble(vrgtr(j)/twt)
                  w  = (cp + oldcp)*0.5
!                  call gauinv(w,vrg,ierr)
!                  if(ierr.eq.1) vrg = UNEST
                  call blinvnormdist(vrg)
                  write(lout,201) vrtr(j),vrg
 201              format(f12.5,1x,f12.5)
                  oldcp =  cp
!
! Now, reset the weight to the normal scores value:
!
                  vrgtr(j) = vrg
            end do
            close(lout)
      end if
!
! Now, read the data if the file exists:
!
      inquire(file=datafl,exist=testfl)
      if(testfl) then
            write(*,*) 'SGSIM: Reading input data'
            open(lin,file=datafl,status='OLD')
            read(lin,*,err=99)
            read(lin,*,err=99) nvari
            do i=1,nvari
                  read(lin,*,err=99)
            end do
            if(ixl.gt.nvari.or.iyl.gt.nvari.or.izl.gt.nvari.or.&
               ivrl.gt.nvari.or.isecvr.gt.nvari.or.iwt.gt.nvari) then
                  write(*,*) 'ERROR: you have asked for a column number'
                  write(*,*) '       greater than available in file'
                  stop
            end if
!
! Read all the data until the end of the file:
!
            twt = 0.d0
            nd  = 0
            nt  = 0
 5          read(lin,*,end=6,err=99) (var(j),j=1,nvari)
            if(var(ivrl).lt.tmin.or.var(ivrl).ge.tmax) then
                  nt = nt + 1
                  go to 5
            end if
            nd = nd + 1
!
! Acceptable data, assign the value, X, Y, Z coordinates, and weight:
!
            vr(nd) = var(ivrl)
            if(ixl.le.0) then
                  x(nd) = xmn
            else
                  x(nd) = var(ixl)
            endif
            if(iyl.le.0) then
                  y(nd) = ymn
            else
                  y(nd) = var(iyl)
            endif
            if(izl.le.0) then
                  z(nd) = zmn
            else
                  z(nd) = var(izl)
            endif
            if(iwt.le.0) then
                  wt(nd) = 1.d0
            else
                  wt(nd) = var(iwt)
            endif
            if(isecvr.le.0) then
                  sec(nd) = UNEST
            else
                  sec(nd) = var(isecvr)
            endif
!
! Normal scores transform?
!
            if(itrans.eq.1) then
                  vrr = vr(nd)
                  call locate(vrtr,ntr,1,ntr,vrr,j)
                  j   = min(max(1,j),(ntr-1))
                  vrg = powint(vrtr(j),vrtr(j+1),vrgtr(j),vrgtr(j+1),&
                               vrr,1.d0)
                  if(vrg.lt.vrgtr(1)  ) vrg = vrgtr(1)
                  if(vrg.gt.vrgtr(ntr)) vrg = vrgtr(ntr)
                  vr(nd) = vrg
            end if
            twt = twt + wt(nd)
            av  = av  + var(ivrl)*wt(nd)
            ss  = ss  + var(ivrl)*var(ivrl)*wt(nd)
            go to 5
 6          close(lin)
!
! Compute the averages and variances as an error check for the user:
!
            av = av / max(twt,EPSLON)
            ss =(ss / max(twt,EPSLON)) - av * av
            write(ldbg,111) nd,nt,av,ss
            write(*,   111) nd,nt,av,ss
 111  format(/,' Data for SGSIM: Number of acceptable data  = ',i8,/,&
               '                 Number trimmed             = ',i8,/,&
               '                 Weighted Average           = ',f12.4,/,&
               '                 Weighted Variance          = ',f12.4,/)
      endif
!
! Read secondary attribute model if necessary:
!
!!$      if(ktype.ge.2) then
!!$            write(*,*) 'Reading secondary attribute file'
!!$            inquire(file=lvmfl,exist=testfl)
!!$            if(.not.testfl) then
!!$                  write(*,104) lvmfl
!!$ 104              format('WARNING secondary attribute file ',a40,
!!$     +             ' does not exist!')
!!$                  stop
!!$            end if
!!$            open(llvm,file=lvmfl,status='OLD')
!!$            read(llvm,*,err=97)
!!$            read(llvm,*,err=97) nvaril
!!$            do i=1,nvaril
!!$                  read(llvm,*,err=97)
!!$            end do
!!$            index = 0
!!$             
!!$            av = 0.d0
!!$            ss = 0.d0
!!$            ns = 0
!!$            do iz=1,nz
!!$                  do iy=1,ny
!!$                        do ix=1,nx
!!$                           index = index + 1
!!$                           read(llvm,*,err=97) (var(j),j=1,nvaril)
!!$                           vrr = var(icollvm)
!!$                           lvm(index) = vrr
!!$                           sim(index) = dble(index)
!
! Do we to transform the secondary variable for a local mean?
!
!!$                           if(trans.and.ktype.eq.2.and.itrans.eq.1) then
!!$                                 if(vrr.le.tmin.or.vrr.ge.tmax) then
!!$                                       lvm(index) = -1.0e21
!!$                                 else   
!!$                                 call locate(vrtr,ntr,1,ntr,vrr,j)
!!$                                 j   =min(max(1,j),(ntr-1))
!!$                                 vrg =powint(vrtr(j),vrtr(j+1),vrgtr(j),
!!$     +                                       vrgtr(j+1),vrr,1.0)
!!$                                 if(vrg.lt.vrgtr(1)  ) vrg = vrgtr(1)
!!$                                 if(vrg.gt.vrgtr(ntr)) vrg = vrgtr(nd)
!!$                                 lvm(index) = vrg
!!$                                 end if
!!$                           end if
!!$                           if(vrr.ge.tmin.or.vrr.le.tmax) then
!!$                                 av = av + vrr
!!$                                 ss = ss + vrr*vrr
!!$                                 ns = ns + 1
!!$                           end if   
!!$                        end do
!!$                  end do
!!$            end do
!!$            ns = max(ns,1)
!!$            av = av / dble(ns)
!!$            ss =(ss / dble(ns)) - av * av
!!$            write(ldbg,112) ns,av,ss
!!$            write(*,   112) ns,av,ss
!!$ 112  format(/,' Secondary Data: Number of data             = ',i8,/,
!!$     +         '                 Equal Weighted Average     = ',f12.4,/,
!!$     +         '                 Equal Weighted Variance    = ',f12.4,/)
!
! Do we need to work with data residuals? (Locally Varying Mean)
!
!!$            if(ktype.eq.2) then
!!$                  do i=1,nd
!!$                        call getindx(nx,xmn,xsiz,x(i),ix,testind)
!!$                        call getindx(ny,ymn,ysiz,y(i),iy,testind)
!!$                        call getindx(nz,zmn,zsiz,z(i),iz,testind)
!!$                        index = ix + (iy-1)*nx + (iz-1)*nxy
!!$                        sec(i) = lvm(index)
!
! Calculation of residual moved to krige subroutine: vr(i)=vr(i)-sec(i)
!
!!$                  end do
!!$            end if
!
! Do we need to get an external drift attribute for the data?
!
!!$            if(ktype.eq.3) then
!!$                  do i=1,nd
!!$                        if(sec(i).eq.UNEST) then
!!$                              call getindx(nx,xmn,xsiz,x(i),ix,testind)
!!$                              call getindx(ny,ymn,ysiz,y(i),iy,testind)
!!$                              call getindx(nz,zmn,zsiz,z(i),iz,testind)
!!$                              index = ix + (iy-1)*nx + (iz-1)*nxy
!!$                              sec(i) = lvm(index)
!!$                        end if
!!$                  end do
!!$            end if
!
! Transform the secondary attribute to normal scores?
!
!!$            if(trans.and.ktype.eq.4) then
!!$                  write(ldbg,113) varred
!!$ 113              format(/,' Transforming Secondary Data with',
!!$     +                     ' variance reduction of ',f12.4,/)
!!$                  write(*,*) 'Transforming secondary variable'
!!$                  write(*,*)
!!$                  call sortem(1,nxyz,lvm,1,sim,c,d,e,f,g,h)
!!$                  oldcp = 0.d0
!!$                  cp    = 0.d0
!!$                  do i=1,nxyz
!!$                        if(lvm(i).gt.tmin.and.lvm(i).le.tmax) then
!!$                              cp =  cp + dble(1.0/dble(ns))
!!$                              w  = (cp + oldcp)/2.0
!!$                              call gauinv(w,lvm(i),ierr)
!!$                              lvm(i) = lvm(i) * varred
!!$                              oldcp  =  cp
!!$                        end if      
!!$                  end do
!!$                  call sortem(1,nxyz,sim,1,lvm,c,d,e,f,g,h)
!!$            end if
!!$      end if
!
! Open the output file:
!

         open(lout,file=outfl,status='UNKNOWN')
         write(lout,210)
 210     format('SGSIM Realizations')
         write(lout,211) 1,nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz,nsim
 211     format(i2,3(1x,i4),3(1x,g14.8),3(1x,g12.6),i4) 
         write(lout,212)
 212     format('value')
      return
!
! Error in an Input File Somewhere:
!
 97   stop 'ERROR in secondary data file!'
 98   stop 'ERROR in parameter file!'
 99   stop 'ERROR in data file!'

      end


      subroutine readparm2(paramfl,cdata_sz,cdata_idx,c_idx_siz,&
                           real_sz,real_idx,r_idx_siz,int_sz,int_idx,&
                           i_idx_siz,cond_option,rseed)
!-----------------------------------------------------------------------
!
!                  Initialization and Read Parameters
!                  **********************************
!
! The input parameters and data are read in from their files. Some quick
! error checking is performed and the statistics of all the variables
! being considered are written to standard output.
!
!
!
!-----------------------------------------------------------------------

      use geostat2
      implicit none
      include  'sgsim2.inc'

      integer   rseed
      integer   MAXSB,cond_option
      real*8    p,acorni
      character paramfl*512
      logical   testfl,trans

      integer cdata_sz,int_sz, real_sz
      integer c_idx_siz,r_idx_siz,i_idx_siz
      integer cdata_idx(c_idx_siz),int_idx(i_idx_siz),real_idx(r_idx_siz)

      integer lin,lout
      integer idum,i
      integer MAXDAT,MAXXYZ
      integer nxt,nyt,nzt
      real*8  aa1,aa2,xx,radius1,radius2
      real*8  xmnt,ymnt,zmnt,xsizt,ysizt,zsizt,sill
      character outfl*64,dbgfl*64,str*64
      integer myprocid

!
! Input/Output units used:
!
      lin  = 1
      lout = 2
      ldbg = 3
      llvm = 4

      open(lin,file=paramfl,status='OLD')
!
! Find Start of Parameters:
!
 1    read(lin,'(a4)',end=98) str(1:4)
      if(str(1:4).ne.'STAR') go to 1
!
! Read Input Parameters:
!
!      read(lin,*,err=98) do_cond
      do_cond = 1
      read(lin,'(a64)',err=98) datafl
      call chknam(datafl,64)

      read(lin,*,err=98) ixl,iyl,izl,ivrl,iwt,isecvr
      read(lin,*,err=98) tmin,tmax
      read(lin,*,err=98) itrans
      read(lin,'(a64)',err=98) transfl
      call chknam(transfl,64)

      read(lin,*,err=98) ismooth
      read(lin,'(a64)',err=98) smthfl
      call chknam(smthfl,64)

      read(lin,*,err=98) isvr,iswt
      read(lin,*,err=98) zmin,zmax
      read(lin,*,err=98) ltail,ltpar
      read(lin,*,err=98) utail,utpar
      read(lin,*,err=98) idbg
      read(lin,'(a64)',err=98) dbgfl
      call chknam(dbgfl,64)
      open(ldbg,file=dbgfl,status='UNKNOWN')

      read(lin,'(a64)',err=98) outfl
      call chknam(outfl,64)

      read(lin,*,err=98) nsim
      read(lin,*,err=98) nxt,xmnt,xsizt
      read(lin,*,err=98) nyt,ymnt,ysizt
      read(lin,*,err=98) nzt,zmnt,zsizt
      read(lin,*,err=98) ixv(1)
      rseed = ixv(1)
      call bl_pd_myproc(myprocid)
      call blutilinitrand(rseed + myprocid)
!      do i=1,1000
!         p = acorni(idum)
!      end do

      read(lin,*,err=98) ndmin,ndmax
      read(lin,*,err=98) nodmax
      read(lin,*,err=98) sstrat
      if(sstrat.eq.1) ndmax = 0
      read(lin,*,err=98) mults,nmult
      read(lin,*,err=98) noct
      read(lin,*,err=98) radius,radius1,radius2
      if(radius.lt.EPSLON) stop 'radius must be greater than zero'
      radsqd = radius  * radius
      sanis1 = radius1 / radius
      sanis2 = radius2 / radius

      read(lin,*,err=98) sang1,sang2,sang3
      read(lin,*,err=98) mxctx,mxcty,mxctz
      read(lin,*,err=98) ktype

      if (ktype.ge.2) stop 'ktype >= 2 not yet supported'
      
      trans = .true.
      if(ktype.lt.0) then
         trans = .false.
         ktype = abs(ktype)
      end if
      colocorr = 0.d0
      if(ktype.eq.4) then
         backspace lin
         read(lin,*,err=98) i,colocorr
         varred = 1.d0
         backspace lin
         read(lin,*,err=9990) i,xx,varred
 9990    continue
      end if
      read(lin,'(a64)',err=98) lvmfl
      call chknam(lvmfl,64)
      read(lin,*,err=98) icollvm
      read(lin,*,err=98) nst(1),c0(1)
      sill = c0(1)

      if(nst(1).le.0) then
            write(*,9997) nst(1)
 9997       format(' nst must be at least 1, it has been set to ',i4,/,&
                   ' The c or a values can be set to zero')
            stop
      endif

      do i=1,nst(1)
            read(lin,*,err=98) it(i),cc(i),ang1(i),ang2(i),ang3(i)
            read(lin,*,err=98) aa(i),aa1,aa2
            anis1(i) = aa1 / max(aa(i),EPSLON)
            anis2(i) = aa2 / max(aa(i),EPSLON)
            sill     = sill + cc(i)
            if(it(i).eq.4) then
               write(*,*) ' A power model is NOT allowed '
               write(*,*) ' Choose a different model and re start '
               stop
            endif
      end do

      if (idbg.gt.1 .and. myprocid.eq.0) then
!         write(*,*) ' do_cond = ',do_cond
         write(*,*) ' data file = ',datafl(1:40)
         write(*,*) ' input columns = ',ixl,iyl,izl,ivrl,iwt,isecvr
         write(*,*) ' trimming limits = ',tmin,tmax
         write(*,*) ' transformation flag = ',itrans
         write(*,*) ' transformation file = ',transfl(1:40)
         write(*,*) ' consider smoothed distribution (1=yes) = ',ismooth
         write(*,*) ' file with smoothed distribution = ',smthfl(1:40)
         write(*,*) ' columns = ',isvr,iswt
         write(*,*) ' data limits (tails) = ',zmin,zmax
         write(*,*) ' lower tail = ',ltail,ltpar
         write(*,*) ' upper tail = ',utail,utpar
         write(*,*) ' debugging level = ',idbg
         write(*,*) ' debugging file = ',dbgfl(1:40)
         write(*,*) ' output file ',outfl(1:40)
         write(*,*) ' number of realizations = ',nsim
         write(*,*) ' X grid specification = ',nxt,xmnt,xsizt
         write(*,*) ' Y grid specification = ',nyt,ymnt,ysizt
         write(*,*) ' Z grid specification = ',nzt,zmnt,zsizt
         write(*,*) ' random number seed = ',ixv(1)
         write(*,*) ' min and max data = ',ndmin,ndmax
         write(*,*) ' maximum previous nodes = ',nodmax
         write(*,*) ' two-part search flag = ',sstrat
         write(*,*) ' multiple grid search flag = ',mults,nmult
         write(*,*) ' number of octants = ',noct
         write(*,*) ' search radii = ',radius,radius1,radius2
         write(*,*) ' search anisotropy angles = ',sang1,sang2,sang3
         write(*,*) ' size of covariance lookup = ',mxctx,mxcty,mxctz
         write(*,*) ' kriging type = ',ktype
         if (ktype .eq. 4) then
            write(*,*) ' correlation coefficient = ',colocorr
            write(*,*) ' secondary variable varred = ',varred
         end if
         write(*,*) ' secondary model file = ',lvmfl(1:40)
         write(*,*) ' column in secondary model file = ',icollvm
         write(*,*) ' nst, c0 = ',nst(1),c0(1)
         do i=1,nst(1)
            write(*,*) ' it,cc,ang[1,2,3]; ',it(i),cc(i),&
                 ang1(i),ang2(i),ang3(i)
            write(*,*) ' a1 a2 a3: ',aa(i),aa1,aa2
         end do
!
! Warn the user if the sill is different than 1.0:
!
         if(sill.gt.(1.d0+EPSLON).or.sill.lt.(1.d0-EPSLON)) then
            write(*,*) 'WARNING the sill of your variogram is not 1.0!'
            write(*,*) '        the sill = ',sill
            write(*,*)
         end if
      end if

      close(lin)

!
! Perform some quick error checking:
!
      testfl = .false.
      if(ltail.ne.1.and.ltail.ne.2) then
            write(*,*) 'ERROR invalid lower tail option ',ltail
            write(*,*) '      only allow 1 or 2 - see manual '
            testfl = .true.
      endif
      if(utail.ne.1.and.utail.ne.2.and.utail.ne.4) then
            write(*,*) 'ERROR invalid upper tail option ',ltail
            write(*,*) '      only allow 1,2 or 4 - see manual '
            testfl = .true.
      endif
      if(utail.eq.4.and.utpar.lt.1.d0) then
            write(*,*) 'ERROR invalid power for hyperbolic tail',utpar
            write(*,*) '      must be greater than 1.0!'
            testfl = .true.
      endif
      if(ltail.eq.2.and.ltpar.lt.0.d0) then
            write(*,*) 'ERROR invalid power for power model',ltpar
            write(*,*) '      must be greater than 0.0!'
            testfl = .true.
      endif
      if(utail.eq.2.and.utpar.lt.0.d0) then
            write(*,*) 'ERROR invalid power for power model',utpar
            write(*,*) '      must be greater than 0.0!'
            testfl = .true.
      endif
      if(testfl) stop

      if (do_cond .eq. 0) then
         ndmin  = 0
         ndmax  = 0
         sstrat = 1
      endif

!
! Find the needed parameters:
!
      MAXCTX = mxctx
      MAXCTY = mxcty
      MAXCTZ = mxctz
      MAXCXY = MAXCTX * MAXCTY
      MAXXYZ = MAXCTX * MAXCTY * MAXCTZ
      MAXNOD = nodmax
      MAXKR1 = MAXNOD + ndmax + 1
      MAXSBX = 50
      MAXSBY = 50
      MAXSBZ = 50
      MAXSB  = MAXSBX*MAXSBY*MAXSBZ

      allocate(covtab(MAXCTX,MAXCTY,MAXCTZ),stat = test)
      if(test.ne.0)then
         write(*,*) MAXCTX, MAXCTY,MAXCTZ
         write(*,*)'ERROR 14: Allocation failed',&
                   ' due to insufficient memory.'
         stop
      end if
!
! Determine length of scratch space for nisb,ixsbtosr,ysbtosr,izsbtosr,
! ixnode,iynode,iznode
!
      real_sz = 13 
      if (r_idx_siz .lt. real_sz) then
         print *,'ERROR: r_idx_siz too small, increase to at least',real_sz
      endif
      do i = 1,real_sz
         real_idx(i) = i
      end do
      !TODO: if ktype >= 2, need space for lvm, which will be as big as the biggest box

      int_sz = 6 + MAXSB + 3*8*MAXSB + 3*MAXXYZ + 2*3
      if (i_idx_siz .lt. 15) then
         print *,'ERROR: i_idx_siz too small, increase to at least',15
      endif
      int_idx(1)  = 1
      int_idx(2)  = 2
      int_idx(3)  = 3
      int_idx(4)  = 4
      int_idx(5)  = 5
      int_idx(6)  = 6
      int_idx(7)  = 7
      int_idx(8)  = int_idx(7)  + MAXSB
      int_idx(9)  = int_idx(8)  + 8*MAXSB
      int_idx(10) = int_idx(9)  + 8*MAXSB
      int_idx(11) = int_idx(10) + 8*MAXSB
      int_idx(12) = int_idx(11) + MAXXYZ
      int_idx(13) = int_idx(12) + MAXXYZ
      int_idx(14) = int_idx(13) + MAXXYZ  ! scratch_i(int_idx(14):int_idx(14)+5) = dlo1,dlo2,dlo3,dhi1,dhi2,dhi3
      int_idx(15) = int_idx(14) + 6

      cdata_sz = 0
      cdata_idx = 0
      cond_option = do_cond
      if (do_cond .eq. 1) then

!        Find MAXDAT
         call get_maxdat(MAXDAT,datafl,smthfl)

!        Allocate the needed memory
         cdata_sz = MAXDAT*9 + MAXXYZ
         if (c_idx_siz .lt. 10) then
            print *,'ERROR: c_idx_siz too small, increase to at least',10
         endif
         cdata_idx(1)  = 1
         cdata_idx(2)  = cdata_idx(1) + MAXDAT
         cdata_idx(3)  = cdata_idx(2) + MAXDAT
         cdata_idx(4)  = cdata_idx(3) + MAXDAT
         cdata_idx(5)  = cdata_idx(4) + MAXDAT
         cdata_idx(6)  = cdata_idx(5) + MAXDAT
         cdata_idx(7)  = cdata_idx(6) + MAXDAT
         cdata_idx(8)  = cdata_idx(7) + MAXDAT
         cdata_idx(9)  = cdata_idx(8) + MAXDAT
         cdata_idx(10) = cdata_idx(9) + MAXDAT

         nd  = MAXDAT
         ntr = MAXDAT

      end if
      return
!
! Error in an Input File Somewhere:
!
 97   stop 'ERROR in secondary data file!'
 98   stop 'ERROR in parameter file!'
 99   stop 'ERROR in data file!'

      end

      subroutine read_cond_data(scratch_c,c_sz,c_idx,c_idx_siz)
!-----------------------------------------------------
! Read conditioning data from GSLIB-style files
!-----------------------------------------------------
      use geostat2
      implicit none
      include 'sgsim2.inc'
       
      integer   c_idx_siz
      integer   c_sz,c_idx(c_idx_siz)
      character tmpfl*64
      
      real*8, target ::  scratch_c(c_sz)
      integer :: ndt
      
      x     => scratch_c(c_idx(1):c_idx(2)-1)
      y     => scratch_c(c_idx(2):c_idx(3)-1)
      z     => scratch_c(c_idx(3):c_idx(4)-1)
      vr    => scratch_c(c_idx(4):c_idx(5)-1)
      wt    => scratch_c(c_idx(5):c_idx(6)-1)
      vrtr  => scratch_c(c_idx(6):c_idx(7)-1)
      vrgtr => scratch_c(c_idx(7):c_idx(8)-1)
      close => scratch_c(c_idx(8):c_idx(9)-1)
      sec   => scratch_c(c_idx(9):c_idx(10)-1)

      if (do_cond .eq. 1) then      
!        Read the data if the file exists
         ndt = nd
         call read_datafl(x,y,z,vr,wt,sec,ndt)
!
! Build transformation table and perform normal score transform
!
         if (itrans .eq. 1) then

            if(ismooth.eq.1) then
               tmpfl  = smthfl
               icolvr = isvr
               icolwt = iswt
            else
               tmpfl  = datafl
               icolvr = ivrl
               icolwt = iwt
            end if
            call setup_trans_table(vrtr,vrgtr,nd,ntr,transfl,&
                 tmpfl,icolvr,icolwt,tmin,tmax)
            call norm_score_trans(vr,vrtr,vrgtr,nd,ntr)

         end if
         
      end if

      end


      subroutine get_maxdat(ndat,datafl,smthfl)
      implicit none
      integer   ndat
      character datafl*64, smthfl*64

      integer i,j,lin,nvari,ntmp
      logical testfl

      real*8, allocatable :: var(:)

      lin = 1

      ndat = 100
      inquire(file=datafl,exist=testfl)
      if(testfl)then
         open(lin,file=datafl,status='UNKNOWN')
         read(lin,*,err=98)
         read(lin,*,err=99) nvari
         do i=1,nvari
            read(lin,*)
         end do

         ndat = 0
         allocate(var(nvari))
 33      read(lin,*,end=66,err=98)(var(j),j=1,nvari)
         ndat = ndat + 1
         go to 33
 66      continue
         deallocate(var)
         rewind(lin)
         close(lin)
      end if

      ntmp = 1
      inquire(file=smthfl,exist=testfl)
      if(testfl)then
         open(lin,file=smthfl,status='UNKNOWN')
         read(lin,*,err=98)
         read(lin,*,err=99) nvari
         do i=1,nvari
            read(lin,*)
         end do
         ntmp = 0
         allocate(var(nvari))
 22      read(lin,*,end=55,err=99)(var(j),j=1,nvari)
         ntmp = ntmp + 1
         go to 22
 55      continue
         deallocate(var)
         rewind(lin)
         close(lin)
      end if
      if(ntmp.gt.ndat) ndat = ntmp
      return
 98   stop 'ERROR in parameter file!'
 99   stop 'ERROR in data file!'
      
      end 

      subroutine read_datafl(x,y,z,vr,wt,sec,nvr)

      implicit none
      include 'sgsim2.inc'

      integer nvr
      real*8  x(nvr),y(nvr),z(nvr),vr(nvr),wt(nvr),sec(nvr)

      logical testfl      
      integer i,j,nvari,nt,lin,myprocid
      real*8  twt,av,ss
      real*8,allocatable :: var(:)

      lin = 1

      inquire(file=datafl,exist=testfl)
      if(.not.testfl)then
         write(*,*) "datafl is wrong: ",datafl
         stop
      end if
      call bl_pd_myproc(myprocid)
      if (myprocid.eq.0) then
         write(*,*) 'Reading input data'
      endif
      open(lin,file=datafl,status='OLD')
      read(lin,*,err=99)
      read(lin,*,err=99) nvari
      do i=1,nvari
         read(lin,*,err=99)
      end do
      if(ixl.gt.nvari.or.iyl.gt.nvari.or.izl.gt.nvari.or.&
         ivrl.gt.nvari.or.isecvr.gt.nvari.or.iwt.gt.nvari) then
         write(*,*) 'ERROR: you have asked for a column number'
         write(*,*) '       greater than available in file'
         stop
      end if

      av  = 0.d0
      ss  = 0.d0
      twt = 0.d0
      nd  = 0
      nt  = 0
      allocate(var(nvari))
      do i = 1,nvr
         read(lin,*,err=99) (var(j),j=1,nvari)
         if(var(ivrl).lt.tmin.or.var(ivrl).ge.tmax) then
            nt = nt + 1
         else
            nd = nd + 1
            vr(nd) = var(ivrl)
            if(ixl.le.0) then
               x(nd) = 0.0
            else
               x(nd) = var(ixl)
            endif
            if(iyl.le.0) then
               y(nd) = 0.0
            else
               y(nd) = var(iyl)
            endif
            if(izl.le.0) then
               z(nd) = 0.0
            else
               z(nd) = var(izl)
            endif
            if(iwt.le.0) then
               wt(nd) = 1.d0
            else
               wt(nd) = var(iwt)
            endif
            if(isecvr.le.0) then
               sec(nd) = UNEST
            else
               sec(nd) = var(isecvr)
            endif
            twt = twt + wt(nd)
            av  = av  + var(ivrl)*wt(nd)
            ss  = ss  + var(ivrl)*var(ivrl)*wt(nd)
         end if
      end do

      deallocate(var)
!
! Compute the averages and variances as an error check for the user:
!
      if (myprocid.eq.0) then
         av = av / max(twt,EPSLON)
         ss =(ss / max(twt,EPSLON)) - av * av
         write(*,   111) nd,nt,av,ss
111      format(/,' Data for SGSIM: Number of acceptable data  = ',i8,/,&
              '                 Number trimmed             = ',i8,/,&
              '                 Weighted Average           = ',f12.4,/,&
              '                 Weighted Variance          = ',f12.4,/)
      endif
      close(lin)
      return
 99   stop 'ERROR in data file!'

      end


      subroutine setup_trans_table(vrtr,vrgtr,nvr,ntr,transfl,&
                                   tmpfl,icolvr,icolwt,tmin,tmax)

      implicit none
      character tmpfl*64,transfl*64,str*64
      integer   icolvr,icolwt,nvr,ntr
      real*8    vrtr(nvr),vrgtr(nvr)
      real*8    tmin,tmax

      integer   i,j,lin,lout,nt,nvari,ierr
      real*8    twt,oldcp,cp,w,vrg
      logical   testfl

      real*8,allocatable :: var(:)
      real*8,parameter   :: UNEST=-1.0d20
      real*8,parameter   :: EPSLON=1.0d-20

      lin  = 1
      lout = 2
      
      write(*,*) 'Setting up transformation table'

      inquire(file=tmpfl,exist=testfl)
      if(.not.testfl) then
         write(*,*) 'ERROR: ',tmpfl,' does not exist'
         write(*,*) '       this file is needed! '
         stop
      end if
!
! Open up the file with reference distribution:
!
      open(lin,file=tmpfl,status='UNKNOWN')
      read(lin,'(a40)',err=98) str(1:40)
      read(lin,*,err=99) nvari
      do i=1,nvari
         read(lin,*,err=98)
      end do
      if(icolvr.gt.nvari.or.icolwt.gt.nvari) then
         write(*,*) ' ERROR: too few columns in ref data '
         stop
      endif
!
! Now, read in the actual data:
!
      nt     = 0
      ntr    = 0
      twt    = 0.d0
      allocate(var(nvari))
      do i = 1,nvr
         read(lin,*,err=99) (var(j),j=1,nvari)
         if(var(icolvr).lt.tmin.or.var(icolvr).ge.tmax) then
            nt = nt + 1
         else
            ntr = ntr + 1   
!
! Assign the data value and coordinate location:
!
            vrtr(ntr) = var(icolvr)
            if(icolwt.le.0) then
               vrgtr(ntr) = 1.d0
            else
               vrgtr(ntr) = var(icolwt)
            endif
            if(vrgtr(ntr).le.0.d0) then
               ntr = ntr - 1
               nt  = nt  + 1
            else
               twt = twt + vrgtr(ntr)
            end if
         end if
      end do

      deallocate(var)
      close(lin)
      if(ntr.le.1) then
         write(*,*) 'ERROR: too few data for transformation'
         stop
      end if
!
! Write transformation table:
!
      open(lout,file=transfl,status='UNKNOWN')
!
! Sort data by value:
!
      call dsortem(1,ntr,vrtr,1,vrgtr)
!
! Compute the cumulative probabilities and write transformation table
!
      twt   = max(twt,EPSLON)
      oldcp = 0.d0
      cp    = 0.d0
      do j=1,ntr
         cp =  cp + dble(vrgtr(j)/twt)
         w  = (cp + oldcp)*0.5
!         call gauinv(w,vrg,ierr)
!         if(ierr.eq.1) vrg = UNEST
         call blinvnormdist(vrg)
         write(lout,201) vrtr(j),vrg
201      format(f12.5,1x,f12.5)
         oldcp =  cp
!
! Now, reset the weight to the normal scores value:
!
         vrgtr(j) = vrg
      end do
      close(lout)
      return
 98   stop 'ERROR in parameter file!'
 99   stop 'ERROR in data file!'
      
      end

      subroutine norm_score_trans(vr,vrtr,vrgtr,nvr,ntr)
!-----------------------------------------------------------------------
! Normal scores transform
!-----------------------------------------------------------------------
      implicit none
      integer nvr,ntr
      real*8  vr(nvr),vrtr(ntr),vrgtr(ntr)      

      integer i,idx
      real*8  powint,vrr,vrg

      do i = 1,nvr
         vrr = vr(i)
         call locate(vrtr,ntr,1,ntr,vrr,idx)
         idx = min(max(1,idx),(ntr-1))
         vrg = powint(vrtr(idx),vrtr(idx+1),vrgtr(idx),&
                      vrgtr(idx+1),vrr,1.d0)
         if(vrg.lt.vrgtr(1)) vrg = vrgtr(1)
         if(vrg.gt.vrgtr(ntr)) vrg = vrgtr(ntr)
         vr(i) = vrg
      end do

      end 

      subroutine read_sec()
!-----------------------------------------------------------------------
! Read secondary attribute model if necessary:
!-----------------------------------------------------------------------

      use geostat2
      implicit none
      include "sgsim2.inc"

      logical testfl,testindx,testindy,testindz,trans
      integer i,j,index,ns
      integer iz,iy,ix,ierr
      real*8  av,ss,vrr,vrg,powint,oldcp,cp,w
      real*8  tmp(nxyz)

      real*8, allocatable :: var(:)

      trans = .true.

      write(*,*) 'Reading secondary attribute file'
      inquire(file=lvmfl,exist=testfl)
      if(.not.testfl) then
         write(*,104) lvmfl
 104     format('WARNING secondary attribute file ',a40,&
                ' does not exist!')
         stop
      end if
      open(llvm,file=lvmfl,status='OLD')
      read(llvm,*,err=97)
      read(llvm,*,err=97) nvaril
      do i=1,nvaril
         read(llvm,*,err=97)
      end do
         
      allocate(var(nvaril))
      index = 0
      av = 0.d0
      ss = 0.d0
      ns = 0
      do iz=1,nz
         do iy=1,ny
            do ix=1,nx
               index = index + 1
               read(llvm,*,err=97) (var(j),j=1,nvaril)
               vrr = var(icollvm)
               lvm(index) = vrr
               tmp(index) = dble(index)
!
! Do we to transform the secondary variable for a local mean?
!
               if(trans.and.ktype.eq.2.and.itrans.eq.1) then
                  if(vrr.le.tmin.or.vrr.ge.tmax) then
                     lvm(index) = -1.0e21
                  else   
                     call locate(vrtr,ntr,1,ntr,vrr,j)
                     j=min(max(1,j),(ntr-1))
                     vrg =powint(vrtr(j),vrtr(j+1),vrgtr(j),&
                                 vrgtr(j+1),vrr,1.0d0)
                     if(vrg.lt.vrgtr(1)  ) vrg = vrgtr(1)
                     if(vrg.gt.vrgtr(ntr)) vrg = vrgtr(nd)
                     lvm(index) = vrg
                  end if
               end if
               if(vrr.ge.tmin.or.vrr.le.tmax) then
                  av = av + vrr
                  ss = ss + vrr*vrr
                  ns = ns + 1
               end if
            end do
         end do
      end do
      ns = max(ns,1)
      av = av / dble(ns)
      ss =(ss / dble(ns)) - av * av
      write(*,   112) ns,av,ss
 112  format(/,' Secondary Data: Number of data             = ',i8,/,&
               '                 Equal Weighted Average     = ',f12.4,/,&
               '                 Equal Weighted Variance    = ',f12.4,/)
!
! Do we need to work with data residuals? (Locally Varying Mean)
!
      if(ktype.eq.2) then
         do i=1,nd
            !call getindx(nx,xmn,xsiz,x(i),ix,testind)
            !call getindx(ny,ymn,ysiz,y(i),iy,testind)
            !call getindx(nz,zmn,zsiz,z(i),iz,testind)
            !index = ix + (iy-1)*nx + (iz-1)*nxy
            call getindxmod(dlo(1),dhi(1),xmn,xsiz,x(i),ix,testindx)
            call getindxmod(dlo(2),dhi(2),ymn,ysiz,y(i),iy,testindy)
            call getindxmod(dlo(3),dhi(3),zmn,zsiz,z(i),iz,testindz)
            if (testindx .and. testindy .and. testindz) then
               call get1Didx(ix,iy,iz,dlo,dhi,index)
               sec(i) = lvm(index)
            endif
!
! Calculation of residual moved to krige subroutine: vr(i)=vr(i)-sec(i)
!
         end do
      end if
!
! Do we need to get an external drift attribute for the data?
!
      if(ktype.eq.3) then
         do i=1,nd
            if(sec(i).eq.UNEST) then
               !call getindx(nx,xmn,xsiz,x(i),ix,testind)
               !call getindx(ny,ymn,ysiz,y(i),iy,testind)
               !call getindx(nz,zmn,zsiz,z(i),iz,testind)
               !index = ix + (iy-1)*nx + (iz-1)*nxy
               call getindxmod(dlo(1),dhi(1),xmn,xsiz,x(i),ix,testindx)
               call getindxmod(dlo(2),dhi(2),ymn,ysiz,y(i),iy,testindy)
               call getindxmod(dlo(3),dhi(3),zmn,zsiz,z(i),iz,testindz)
               if (testindx .and. testindy .and. testindz) then
                  call get1Didx(ix,iy,iz,dlo,dhi,index)
                  sec(i) = lvm(index)
               endif
            end if
         end do
      end if
!
! Transform the secondary attribute to normal scores?
!
      if(trans.and.ktype.eq.4) then
         write(ldbg,113) varred
 113      format(/,' Transforming Secondary Data with',&
                   ' variance reduction of ',f12.4,/)
         write(*,*) 'Transforming secondary variable'
         write(*,*)
         call dsortem(1,nxyz,lvm,1,tmp)
         oldcp = 0.d0
         cp    = 0.d0
         do i=1,nxyz
            if(lvm(i).gt.tmin.and.lvm(i).le.tmax) then
               cp =  cp + dble(1.0/dble(ns))
               w  = (cp + oldcp)/2.0
!               call gauinv(w,lvm(i),ierr)
               call blinvnormdist(lvm(i))
               lvm(i) = lvm(i) * varred
               oldcp  =  cp
            end if
         end do
         call dsortem(1,nxyz,tmp,1,lvm)
      end if

!
! Read in the secondary data distribution for this realization:
! There was a condition where isim > 1.  Not sure why.
!
      if(ktype.eq.4) then
         write(*,*)
         write(*,*) ' Reading next secondary model'
         index = 0
         do iz=1,nz
            do iy=1,ny
               do ix=1,nx
                  index = index + 1
                  read(llvm,*,end=977)(var(j),j=1,nvaril)
                  lvm(index) = var(icollvm)
                  tmp(index) = dble(index)
               end do
            end do
         end do
         write(*,*) ' Building CDF from  secondary model'
         call dsortem(1,nxyz,lvm,1,tmp)
         oldcp = 0.d0
         cp    = 0.d0
         do i=1,nxyz
            cp =  cp + 1.d0/dble(nxyz)
            w  = (cp + oldcp)/2.0
!            call gauinv(w,lvm(i),ierr)
            call blinvnormdist(lvm(i))
            lvm(i) = lvm(i) * varred
            oldcp  =  cp
         end do
         write(*,*) ' Restoring order of secondary model'
         call dsortem(1,nxyz,tmp,1,lvm)
 977     continue
         write(*,*)
      end if
      
      deallocate(var)
      return
 97   stop 'ERROR in secondary data file!'

      end


      subroutine sgsim(sim)
!-----------------------------------------------------------------------
!
!           Conditional Simulation of a 3-D Rectangular Grid
!           ************************************************
!
! This subroutine generates 3-D realizations of a Gaussian process with
! a given autocovariance model, and conditional to input Gaussian data.
! The conditional simulation is achieved by sequential simulation of all
! the nodes visited by a random path.
!
!
!
! PROGRAM NOTES:
!
!  1. The three dimensional anisotropy parameters, i.e., of the search
!     ellipse and variogram ranges are described in section 2.3 of the
!     manual.   The variogram parameters are described in the same place
!
!  2. The original data and previously simulated grid nodes can be
!     searched separately.  There can be a different maximum number of
!     each and a minimum number of original data can be specified
!     to restrict simulation beyond the limits of the data.  The
!     closeness of previously simulated grid nodes is measured according
!     to the variogram structural distance.
!
!
!
! INPUT VARIABLES:
!
!   nd               Number of data (no missing values)
!   x,y,z(nd)        coordinates of the data
!   vr(nd)           gaussian data (normal scores)
!
!   nx,ny,nz         Number of blocks in X,Y, and Z
!   xmn,ymn,zmn      Coordinate at the center of the first Block
!   xsiz,ysiz,zsiz   spacing of the grid nodes (block size)
!
!   nsim             number of simulations
!   ktype            =1, ordinary kriging; =0, simple kriging
!   sim              the current realization
!   idbg             integer debugging level (0=none,2=normal,4=serious)
!   ldbg             unit number for the debugging output
!   lout             unit number for the output
!
!   radius           Maximum search radius
!   sang1            Azimuth angle of the principal search direction
!   sang2            Dip angle of the principal search direction
!   sang3            Third rotation angle of the search ellipse
!   sanis1           Anisotropy for the dip angle
!   sanis2           Anisotropy for the plunge angle
!   ndmin            Minimum number of data required before sim
!   ndmax            Maximum number of samples for simulation
!   noct             Maximum number per octant if an octant search is
!                      desired (if <= 0, then no octant search)
!
!   nodmax           Maximum number of previously simulated grid nodes
!                      to consider in the simulation.  The structural
!                      variogram distance is used to identify close ones
!
!   c0               Nugget constant (isotropic).
!   cc(nst)          Multiplicative factor of each nested structure.
!   aa(nst)          Parameter "a" of each nested structure.
!   it(nst)          Type of nested structures (1=sph,2=exp,3=gau,4=pow)
!   ang1(nst)        Azimuth angle for the principal direction
!   ang2(nst)        Dip angle for the principal direction
!   ang3(nst)        Third rotation angle to rotate the two minor
!                      directions around the principal direction
!   anis1(nst)       Anisotropy (radius in minor direction at 90
!                      degrees from "ang1" divided by the principal
!                      radius in direction "ang1")
!   anis2(nst)       Anisotropy (radius in minor direction at 90 degrees
!                      vertical from "ang1" divided by the principal
!                      radius in direction "ang1")
!
!
! OUTPUT VARIABLES:  Simulated Values are written to "lout"
!
!
!
! EXTERNAL REFERENCES:
!
!   super            Sets up the super block search of original data
!   search           Search for nearby data values
!   ctable           Builds a covariance table and "spiral" search
!   srchnd           Search for nearby simulated grid nodes
!   sqdist           computes anisotropi! squared distance
!   sortem           sorts multiple arrays in ascending order (separate)
!   cova3            Calculates the covariance given a variogram model
!   krige            Sets up and solves either the SK or OK system
!   ksol             Linear system solver using Gaussian elimination
!
!
!
! Concepts taken from F. Alabert and E. Isaaks
!
!-----------------------------------------------------------------------

      use       geostat
      include  'sgsim.inc'

      real*8  sim(nxyz)

      integer ix,iy,iz,id2,id,ind,i,j,nsbtosr,nsec,is,ierr,index,ne,isim
      integer lktype,infoct,idum,in,irepo,jx,jy,jz,nnx,nny,nnz,imult
      integer nxsup,nysup,nzsup
      real*8  ss,av,xx,yy,zz
      real*8  var(10),sec2,sec3
      real*8  p,acorni,cp,oldcp,w,test2,xp,cstdev,cmean,gmean
      real*8  backtr,simval,TINY
      real*8  xmnsup,ymnsup,zmnsup,xsizsup,ysizsup,zsizsup
      logical testind
      allocate(order(nxyz))
      allocate(lvm(nxyz))

!
! Set up the rotation/anisotropy matrices that are needed for the
! variogram and search.
!
      write(*,*) 'Setting up rotation matrices for variogram and search'
      do is=1,nst(1)
            call setrot(ang1(is),ang2(is),ang3(is),anis1(is),anis2(is),&
                        is,MAXROT,rotmat)
      end do
      isrot = MAXNST + 1
      call setrot(sang1,sang2,sang3,sanis1,sanis2,isrot,MAXROT,rotmat)
!
! Set up the super block search:
!
      if(sstrat.eq.0) then
            write(*,*) 'Setting up super block search strategy'
            nsec = 1
            call setsupr(nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,nd,x,y,z,&
                         vr,wt,nsec,sec,sec2,sec3,MAXSBX,MAXSBY,MAXSBZ,&
                         nisb,nxsup,xmnsup,xsizsup,nysup,ymnsup,ysizsup,&
                         nzsup,zmnsup,zsizsup)
            call picksup(nxsup,xsizsup,nysup,ysizsup,nzsup,zsizsup,&
                         isrot,MAXROT,rotmat,radsqd,nsbtosr,ixsbtosr,&
                         iysbtosr,izsbtosr)
      end if
!
! Set up the covariance table and the spiral search:
!
      call ctable()
!
! MAIN LOOP OVER ALL THE SIMULATIONS:
!
      do isim=1,nsim
!
! Read in the secondary data distribution for this realization:
!
            if(isim.gt.1.and.ktype.eq.4) then
                  write(*,*)
                  write(*,*) ' Reading next secondary model'
                  index = 0
                  do iz=1,nz
                        do iy=1,ny
                              do ix=1,nx
                                 index = index + 1
                                 read(llvm,*,end=977)(var(j),j=1,nvaril)
                                 lvm(index) = var(icollvm)
                                 sim(index) = dble(index)
                              end do
                        end do
                  end do
                  write(*,*) ' Building CDF from  secondary model'
!                  call sortem(1,nxyz,lvm,1,sim,c,d,e,f,g,h)
                  call sortem(1,nxyz,lvm,1,sim)
                  oldcp = 0.d0
                  cp    = 0.d0
                  do i=1,nxyz
                        cp =  cp + dble(1.d0/dble(nxyz))
                        w  = (cp + oldcp)/2.0
!                        call gauinv(w,lvm(i),ierr)
                        call blinvnormdist(lvm(i))
                        lvm(i) = lvm(i) * varred
                        oldcp  =  cp
                  end do
                  write(*,*) ' Restoring order of secondary model'
!                  call sortem(1,nxyz,sim,1,lvm,c,d,e,f,g,h)
                  call sortem(1,nxyz,sim,1,lvm)
 977              continue
                  write(*,*)
            end if
!
! Work out a random path for this realization:
!
            do ind=1,nxyz
!                  sim(ind)   = acorni(idum)
                  call blutilrand(sim(ind))

                  order(ind) = ind
            end do

!
! The multiple grid search works with multiples of 4 (yes, that is
! somewhat arbitrary):
!
            if(mults.eq.1) then
               do imult=1,nmult
                  nnz = max(1,nz/(imult*4))
                  nny = max(1,ny/(imult*4))
                  nnx = max(1,nx/(imult*4))
                  jz  = 1
                  jy  = 1
                  jx  = 1
                  do iz=1,nnz
                     if(nnz.gt.1) jz = iz*imult*4
                     do iy=1,nny
                        if(nny.gt.1) jy = iy*imult*4
                        do ix=1,nnx
                           if(nnx.gt.1) jx = ix*imult*4
                           index = jx + (jy-1)*nx + (jz-1)*nxy
                           sim(index) = sim(index) - imult
                        end do
                     end do
                  end do
               end do
            end if
!            call sortem(1,nxyz,sim,1,order,c,d,e,f,g,h)
            call sortem(1,nxyz,sim,1,order)
!
! Initialize the simulation:
!
            do ind=1,nxyz
               sim(ind) = UNEST
            end do
            write(*,*)
            write(*,*) 'Working on realization number ',isim
!
! Assign the data to the closest grid node:
!
            TINY = 0.0001
            do id=1,nd
               call getindx(nx,xmn,xsiz,x(id),ix,testind)
               call getindx(ny,ymn,ysiz,y(id),iy,testind)
               call getindx(nz,zmn,zsiz,z(id),iz,testind)
               ind = ix + (iy-1)*nx + (iz-1)*nxy
               xx  = xmn + dble(ix-1)*xsiz
               yy  = ymn + dble(iy-1)*ysiz
               zz  = zmn + dble(iz-1)*zsiz
               test = abs(xx-x(id)) + abs(yy-y(id)) + abs(zz-z(id))
!
! Assign this data to the node (unless there is a closer data):
!
               if(sstrat.eq.1) then
                  if(sim(ind).ge.0.d0) then
                     id2 = int(sim(ind)+0.5)
                     test2 = abs(xx-x(id2)) + abs(yy-y(id2)) &
                             + abs(zz-z(id2))
                      if(test.le.test2) sim(ind) = dble(id)
                      write(ldbg,102) id,id2
                   else
                      sim(ind) = dble(id)
                   end if
                end if
!
! Assign a flag so that this node does not get simulated:
!
                  if(sstrat.eq.0.and.test.le.TINY) sim(ind)=10.0*UNEST
            end do
 102        format(' WARNING data values ',2i5,' are both assigned to ',&
                 /,'         the same node - taking the closest')
!
! Now, enter data values into the simulated grid:
!
            do ind=1,nxyz
               id = int(sim(ind)+0.5)
               if(id.gt.0) sim(ind) = vr(id)
            end do
            irepo = max(1,min((nxyz/10),10000))
!
! MAIN LOOP OVER ALL THE NODES:
!
            do in=1,nxyz
                  if((int(in/irepo)*irepo).eq.in) write(*,103) in
 103              format('   currently on node ',i9)
!
! Figure out the location of this point and make sure it has
! not been assigned a value already:
!
                  index = order(in)
                  if(sim(index).gt.(UNEST+EPSLON).or.&
                     sim(index).lt.(UNEST*2.d0)) go to 5
                  iz = ((index-1)/nxy) + 1
                  iy = ((index-(iz-1)*nxy-1)/nx) + 1
                  ix = (index) - (iz-1)*nxy - (iy-1)*nx
                  xx = xmn + dble(ix-1)*xsiz
                  yy = ymn + dble(iy-1)*ysiz
                  zz = zmn + dble(iz-1)*zsiz
!
! Now, we'll simulate the point ix,iy,iz.  First, get the close data
! and make sure that there are enough to actually simulate a value,
! we'll only keep the closest "ndmax" data, and look for previously
! simulated grid nodes:
!
                  if(sstrat.eq.0) then
                        call srchsupr(xx,yy,zz,radsqd,isrot,MAXROT,&
                               rotmat,nsbtosr,ixsbtosr,iysbtosr,&
                               izsbtosr,noct,nd,x,y,z,wt,nisb,nxsup,&
                               xmnsup,xsizsup,nysup,ymnsup,ysizsup,&
                               nzsup,zmnsup,zsizsup,nclose,close,&
                               infoct)
                        if(nclose.lt.ndmin) go to 5
                        if(nclose.gt.ndmax) nclose = ndmax
                  endif
                  call srchnd(sim,ix,iy,iz)
!
! Calculate the conditional mean and standard deviation.  This will be
! done with kriging if there are data, otherwise, the global mean and
! standard deviation will be used:
!
                  if(ktype.eq.2) then
                        gmean = lvm(index)
                  else
                        gmean = 0.d0
                  end if
                  if((nclose+ncnode).lt.1) then
                        cmean  = gmean
                        cstdev = 1.d0
                  else
!
! Perform the kriging.  Note that if there are fewer than four data
! then simple kriging is prefered so that the variance of the
! realization does not become artificially inflated:
!
                        lktype = ktype
                        if(ktype.eq.1.and.(nclose+ncnode).lt.4)lktype=0
                        call krige(ix,iy,iz,xx,yy,zz,lktype,gmean,&
                                   cmean,cstdev)
                  endif
!
! Draw a random number and assign a value to this node:
!
!                  p = acorni(idum)
!                  call gauinv(p,xp,ierr)
                  call blinvnormdist(xp)
                  sim(index) = xp * cstdev + cmean
!                  write(*,*), index, sim(index),xp, cstdev, cmean
                  if(idbg.ge.3) write(ldbg,141) p,sim(index)
 141              format(' random number ',f6.4,' realization ',f7.4)
!
! Quick check for far out results:
!
                  if(abs(cmean).gt.5.0.or.abs(cstdev).gt.5.0.or.&
                     abs(sim(index)).gt.6.0) then
                  write(ldbg,104) ix,iy,iz,cmean,cstdev,sim(index)
  104             format('WARNING: grid node location: ',3i5,/,&
                        '         conditional mean:   ',f12.5,/,&
                        '         conditional stdev:  ',f12.5,/,&
                        '         simulated value:    ',f12.5)
                  endif
!
! END MAIN LOOP OVER NODES:
!
 5                continue
            end do
!
! Do we need to reassign the data to the grid nodes?
!
            if(sstrat.eq.0) then
                  do id=1,nd
                        call getindx(nx,xmn,xsiz,x(id),ix,testind)
                        call getindx(ny,ymn,ysiz,y(id),iy,testind)
                        call getindx(nz,zmn,zsiz,z(id),iz,testind)
                        xx  = xmn + dble(ix-1)*xsiz
                        yy  = ymn + dble(iy-1)*ysiz
                        zz  = zmn + dble(iz-1)*zsiz
                        ind = ix + (iy-1)*nx + (iz-1)*nxy
                        test=abs(xx-x(id))+abs(yy-y(id))+abs(zz-z(id))
                        if(test.le.TINY) sim(ind) = vr(id)
                  end do
            end if
!
! Back transform each value and write results:
!
            ne = 0
            av = 0.d0
            ss = 0.d0
            do ind=1,nxyz
                  simval = sim(ind)
                  if(simval.gt.-9.0.and.simval.lt.9.0) then
                        ne = ne + 1
                        av = av + simval
                        ss = ss + simval*simval
                  end if
                  if(itrans.eq.1.and.simval.gt.(UNEST+EPSLON)) then
                        simval = backtr(simval,ntr,vrtr,vrgtr,zmin,&
                                        zmax,ltail,ltpar,utail,utpar)
                        if(simval.lt.zmin) simval = zmin
                        if(simval.gt.zmax) simval = zmax
                  end if
                  sim(ind) = simval
                  write(lout,'(g14.8)') simval
            end do
            av = av / max(dble(ne),1.d0)
            ss =(ss / max(dble(ne),1.d0)) - av * av
            write(ldbg,111) isim,ne,av,ss
!            write(*,   111) isim,ne,av,ss
 111        format(/,' Realization ',i3,': number   = ',i8,/,&
                     '                  mean     = ',f12.4,&
                     ' (close to 0.0?)',/,&
                     '                  variance = ',f12.4,&
                     ' (close to gammabar(V,V)? approx. 1.0)',/)
!
! END MAIN LOOP OVER SIMULATIONS:
!
      end do
!
! Return to the main program:
!
      return
      end

      subroutine getindxmod(ilo,ihi,min,siz,loc,index,inflag)
! -----------------------------------------------------------------------
!      Gets the coordinate index location of a point within a grid cell
!      -MSD: Modified orig to augment arglist with allowable idx range
!      ***********************************************************
!      ilo     minimum allowable index
!      ihi     maximum allowable index
!      min     location of left side of ilo
!      siz     size of the cells
!      loc     location of the point being considered
!      index   output index within [ilo,ihi]
!      inflag  true if the location is actually in the grid (false otherwise
!              e.g., if the location is outside then index will be set to
!              nearest boundary
! -----------------------------------------------------------------------
      implicit none
      integer   ilo, ihi,index
      real*8    min,siz,loc
      logical   inflag

      ! Compute the index of "loc":
      index = int( (loc-min)/siz ) + ilo

      ! Check to see if in or out:
      if(index.lt.ilo) then
            index  = ilo
            inflag = .false.
      else if(index.gt.ihi) then
            index  = ihi
            inflag = .false.
      else
            inflag = .true.
      end if
      end

      subroutine get3Didx(ix,iy,iz,dlo,dhi,idx1D)
      implicit none
      integer idx1D,ix,iy,iz,dlo(3),dhi(3),nx,ny
      nx = (dhi(1)-dlo(1)+1)
      ny = (dhi(2)-dlo(2)+1)
      iz = int((idx1D-1)/(nx*ny)) + dlo(3)
      iy = int((idx1D-(iz-dlo(3))*nx*ny-1)/nx) + dlo(2)
      ix = idx1D - (iz-dlo(3))*nx*ny - (iy-dlo(2))*nx + dlo(1) - 1
      end

      subroutine get1Didx(ix,iy,iz,dlo,dhi,idx1D)
      implicit none
      integer idx1D,ix,iy,iz,dlo(3),dhi(3),nx,ny
      nx = (dhi(1)-dlo(1)+1)
      ny = (dhi(2)-dlo(2)+1)
      idx1D = ix-dlo(1) + (iy-dlo(2))*nx + (iz-dlo(3))*nx*ny + 1
      end

      double precision function getCCloc(idx,mn,lo,sz)
      implicit none
      integer idx,lo
      double precision mn, sz
      getCCloc = mn + dble(idx-lo+0.5d0)*sz
      end

      subroutine sgsim_setup(sim,sim_sz,scratch_c,c_sz,c_idx,&
           c_idx_siz,scratch_r,real_sz,real_idx,r_idx_siz,&
           scratch_i,int_sz,int_idx,i_idx_siz)
!-----------------------------------------------------------------------
!
!    SGSIM setup: Only for a single instance. 
!
!-----------------------------------------------------------------------

      use       geostat2
      implicit none
      include  'sgsim2.inc'

      integer, intent(in) ::  sim_sz,c_sz,real_sz,int_sz
      integer   c_idx_siz, r_idx_siz, i_idx_siz
      integer   c_idx(c_idx_siz),real_idx(r_idx_siz),int_idx(i_idx_siz)
      real*8    sim(sim_sz)

      integer, target :: scratch_i(int_sz)
      real*8 , target :: scratch_c(c_sz)
      real*8 , target :: scratch_r(real_sz)
      
      logical   testindx, testindy, testindz
      integer   ix,iy,iz,id2,id,ind,nsbtosr,nsec,is
      real*8    xx,yy,zz
      real*8    sec2,sec3
      real*8    test2,TINY
      double precision getCCloc

      nx       => scratch_i(int_idx(1))
      ny       => scratch_i(int_idx(2))
      nz       => scratch_i(int_idx(3))
      nxsup    => scratch_i(int_idx(4))
      nysup    => scratch_i(int_idx(5))
      nzsup    => scratch_i(int_idx(6))
      nisb     => scratch_i(int_idx(7):int_idx(8)-1)
      ixsbtosr => scratch_i(int_idx(8):int_idx(9)-1)
      iysbtosr => scratch_i(int_idx(9):int_idx(10)-1)
      izsbtosr => scratch_i(int_idx(10):int_idx(11)-1)
      ixnode   => scratch_i(int_idx(11):int_idx(12)-1)
      iynode   => scratch_i(int_idx(12):int_idx(13)-1)
      iznode   => scratch_i(int_idx(13):int_idx(14)-1)
      dlo      => scratch_i(int_idx(14):int_idx(14)+2)
      dhi      => scratch_i(int_idx(14)+3:int_idx(14)+5)

      xmn      => scratch_r(real_idx(1))
      ymn      => scratch_r(real_idx(2))
      zmn      => scratch_r(real_idx(3))
      xsiz     => scratch_r(real_idx(4))
      ysiz     => scratch_r(real_idx(5))
      zsiz     => scratch_r(real_idx(6))
      xmnsup   => scratch_r(real_idx(7))
      ymnsup   => scratch_r(real_idx(8))
      zmnsup   => scratch_r(real_idx(9))
      xsizsup  => scratch_r(real_idx(10))
      ysizsup  => scratch_r(real_idx(11))
      zsizsup  => scratch_r(real_idx(12))
      if (ktype .ge. 2) then
         lvm      => scratch_r(real_idx(13):real_idx(13)+nx*ny*nz-1)
      endif

      if (do_cond .gt. 0) then
         x     => scratch_c(c_idx(1):c_idx(2)-1)
         y     => scratch_c(c_idx(2):c_idx(3)-1)
         z     => scratch_c(c_idx(3):c_idx(4)-1)
         vr    => scratch_c(c_idx(4):c_idx(5)-1)
         wt    => scratch_c(c_idx(5):c_idx(6)-1)
         vrtr  => scratch_c(c_idx(6):c_idx(7)-1)
         vrgtr => scratch_c(c_idx(7):c_idx(8)-1)
         close => scratch_c(c_idx(8):c_idx(9)-1)
         sec   => scratch_c(c_idx(9):c_idx(10)-1)
      end if
! 
! Set up the rotation/anisotropy matrices that are needed for the
! variogram and search.
!
      do is=1,nst(1)
         call setrot(ang1(is),ang2(is),ang3(is),anis1(is),anis2(is),&
                     is,MAXROT,rotmat)
      end do
      isrot = MAXNST + 1
      call setrot(sang1,sang2,sang3,sanis1,sanis2,isrot,MAXROT,rotmat)
!
! Set up the super block search:
!
      if(sstrat.eq.0) then
         write(*,*) 'Setting up super block search strategy'
         nsec = 1
         call setsupr(nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,nd,x,y,z,&
                     vr,wt,nsec,sec,sec2,sec3,MAXSBX,MAXSBY,MAXSBZ,&
                     nisb,nxsup,xmnsup,xsizsup,nysup,ymnsup,ysizsup,&
                     nzsup,zmnsup,zsizsup)
         call picksup(nxsup,xsizsup,nysup,ysizsup,nzsup,zsizsup,&
                     isrot,MAXROT,rotmat,radsqd,nsbtosr,ixsbtosr,&
                     iysbtosr,izsbtosr)
      end if
!
! Set up the covariance table and the spiral search:
!      
!      write(*,*) 'Constructing covariance matrix'
      call ctable2()
!
! Assign the data to the closest grid node:
!
      sim = UNEST
      TINY = 0.0001
      do id=1,nd

         call getindxmod(dlo(1),dhi(1),xmn,xsiz,x(id),ix,testindx)
         call getindxmod(dlo(2),dhi(2),ymn,ysiz,y(id),iy,testindy)
         call getindxmod(dlo(3),dhi(3),zmn,zsiz,z(id),iz,testindz)

         if (testindx .and. testindy .and. testindz) then

            xx = getCCloc(ix,xmn,dlo(1),xsiz)
            yy = getCCloc(iy,ymn,dlo(2),ysiz)
            zz = getCCloc(iz,zmn,dlo(3),zsiz)

            test = abs(xx-x(id)) + abs(yy-y(id)) + abs(zz-z(id))
!
! Assign this data to the node (unless there is a closer data):
!
            call get1Didx(ix,iy,iz,dlo,dhi,ind)
            if(sstrat.eq.1) then
               if (ind .gt. sim_sz) then
                  print *,'dlo,dhi',dlo,dhi
                  print *,'ix,iy,iz,ind:',ix,iy,iz,ind
                  print *,'bust in sgsim_setup',ind,sim_sz
                  stop
               endif
               if(sim(ind).ge.0.d0) then
                  id2 = int(sim(ind)+0.5)
                  test2 = abs(xx-x(id2)) + abs(yy-y(id2))&
                       + abs(zz-z(id2))
                  if(test.le.test2) sim(ind) = dble(id)
                  write(ldbg,102) id,id2
               else
                  sim(ind) = dble(id)
               end if
            end if
!
! If not assign to node, but coincidently lands on node,
! assign a flag so that this node does not get simulated:
!
            if(sstrat.eq.0.and.test.le.TINY) sim(ind)=10.0*UNEST

         else
            ! print *,'conditioning data outside domain'
            ! print *,'   (x,y,z) = (',x(id),y(id),z(id),')'
            ! print *,'   domain_lo = (',xmn,ymn,zmn,')'
            ! print *,'   domain_hi = (',&
            !      xmn+dble(dhi(1)-dlo(1)+1)*xsiz,&
            !      ymn+dble(dhi(2)-dlo(2)+1)*ysiz,&
            !      zmn+dble(dhi(3)-dlo(3)+1)*zsiz,')'
            ! print *,testindx,testindy,testindz
         endif
      end do

 102  format(' WARNING data values ',2i5,' are both assigned to ',&
             /,'         the same node - taking the closest')

!
! Now, enter data values into the simulated grid:
!
      do ind=1,sim_sz
         if(sim(ind).gt.0) then
            id = int(sim(ind)+0.5)
            if(id.gt.0) sim(ind) = vr(id)
         endif
      end do
      end

      subroutine sgsim_single_iter(index,sim,sim_sz,&
           scratch_c,c_sz,c_idx,c_idx_siz,&
           scratch_r,real_sz,real_idx,r_idx_siz,&
           scratch_i,int_sz,int_idx,i_idx_siz)
!-----------------------------------------------------------------------
!
!    SGSIM iteration: Only for a single iteration. 
!
!-----------------------------------------------------------------------      

      use       geostat2
      implicit none
      include  'sgsim2.inc'

      integer index
      integer sim_sz,c_sz,real_sz,int_sz
      integer c_idx_siz, r_idx_siz, i_idx_siz
      integer c_idx(c_idx_siz),real_idx(r_idx_siz),int_idx(i_idx_siz)
      real*8  sim(sim_sz)

      integer, target :: scratch_i(int_sz) 
      real*8 , target :: scratch_c(c_sz)
      real*8 , target :: scratch_r(real_sz)

      integer ix,iy,iz,ierr,lktype,infoct,nsbtosr,idum
      integer icnode(MAXNOD)
      real*8  xx,yy,zz
      real*8  p,acorni,xp,cstdev,cmean,gmean
      real*8  cnodex(MAXNOD),cnodey(MAXNOD),cnodez(MAXNOD)
      real*8  cnodev(MAXNOD)
      double precision getCCloc

      nx       => scratch_i(int_idx(1))
      ny       => scratch_i(int_idx(2))
      nz       => scratch_i(int_idx(3))
      nxsup    => scratch_i(int_idx(4))
      nysup    => scratch_i(int_idx(5))
      nzsup    => scratch_i(int_idx(6))
      nisb     => scratch_i(int_idx(7):int_idx(8)-1)
      ixsbtosr => scratch_i(int_idx(8):int_idx(9)-1)
      iysbtosr => scratch_i(int_idx(9):int_idx(10)-1)
      izsbtosr => scratch_i(int_idx(10):int_idx(11)-1)
      ixnode   => scratch_i(int_idx(11):int_idx(12)-1)
      iynode   => scratch_i(int_idx(12):int_idx(13)-1)
      iznode   => scratch_i(int_idx(13):int_idx(14)-1)
      dlo      => scratch_i(int_idx(14):int_idx(14)+2)
      dhi      => scratch_i(int_idx(14)+3:int_idx(14)+5)

      xmn      => scratch_r(real_idx(1))
      ymn      => scratch_r(real_idx(2))
      zmn      => scratch_r(real_idx(3))
      xsiz     => scratch_r(real_idx(4))
      ysiz     => scratch_r(real_idx(5))
      zsiz     => scratch_r(real_idx(6))
      xmnsup   => scratch_r(real_idx(7))
      ymnsup   => scratch_r(real_idx(8))
      zmnsup   => scratch_r(real_idx(9))
      xsizsup  => scratch_r(real_idx(10))
      ysizsup  => scratch_r(real_idx(11))
      zsizsup  => scratch_r(real_idx(12))
      if (ktype .ge. 2) then
         lvm      => scratch_r(real_idx(13):real_idx(13)+nx*ny*nz-1)
      endif

      x        => scratch_c(c_idx(1):c_idx(2)-1)
      y        => scratch_c(c_idx(2):c_idx(3)-1)
      z        => scratch_c(c_idx(3):c_idx(4)-1)
      vr       => scratch_c(c_idx(4):c_idx(5)-1)
      wt       => scratch_c(c_idx(5):c_idx(6)-1)
      vrtr     => scratch_c(c_idx(6):c_idx(7)-1)
      vrgtr    => scratch_c(c_idx(7):c_idx(8)-1)
      close    => scratch_c(c_idx(8):c_idx(9)-1)
      sec      => scratch_c(c_idx(9):c_idx(10)-1)

!
! Figure out the location of this point and make sure it has
! not been assigned a value already:
!      
      if(sim(index).gt.(UNEST+EPSLON).or.&
         sim(index).lt.(UNEST*2.0)) return

      call get3Didx(ix,iy,iz,dlo,dhi,index)
      xx = getCCloc(ix,xmn,dlo(1),xsiz)
      yy = getCCloc(iy,ymn,dlo(2),ysiz)
      zz = getCCloc(iz,zmn,dlo(3),zsiz)
!
! Now, we'll simulate the point ix,iy,iz.  First, get the close data
! and make sure that there are enough to actually simulate a value,
! we'll only keep the closest "ndmax" data, and look for previously
! simulated grid nodes:
!
      if(sstrat.eq.0) then
         call srchsupr(xx,yy,zz,radsqd,isrot,MAXROT,&
                      rotmat,nsbtosr,ixsbtosr,iysbtosr,&
                      izsbtosr,noct,nd,x,y,z,wt,nisb,nxsup,&
                      xmnsup,xsizsup,nysup,ymnsup,ysizsup,&
                      nzsup,zmnsup,zsizsup,nclose,close,&
                      infoct)
         if(nclose.lt.ndmin) return
         if(nclose.gt.ndmax) nclose = ndmax
      endif
      call srchnd2(icnode,cnodex,cnodey,cnodez,cnodev,sim,sim_sz,ix,iy,iz)
!
! Calculate the conditional mean and standard deviation.  This will be
! done with kriging if there are data, otherwise, the global mean and
! standard deviation will be used:
!
      if(ktype.eq.2) then
         gmean = lvm(index)
      else
         gmean = 0.d0
      end if

      if((nclose+ncnode).lt.1) then
         cmean  = gmean
         cstdev = 1.d0
      else
!
! Perform the kriging.  Note that if there are fewer than four data
! then simple kriging is prefered so that the variance of the
! realization does not become artificially inflated:
!

         lktype = ktype
         if(ktype.eq.1.and.(nclose+ncnode).lt.4)lktype=0

         call krige2(icnode,cnodex,cnodey,cnodez,cnodev,&
                   ix,iy,iz,xx,yy,zz,lktype,gmean,&
                   cmean,cstdev)
     endif
!
! Draw a random number and assign a value to this node:
!
!      p = acorni(idum)
!      call gauinv(p,xp,ierr)
      call blinvnormdist(xp)
      sim(index) = xp * cstdev + cmean
!      write(*,*) index,sim(index),xp,cstdev,cmean
      if(idbg.ge.3) write(ldbg,141) p,sim(index)
 141  format(' random number ',f6.4,' realization ',f7.4)
!
! Quick check for far out results:
!
      if(abs(cmean).gt.5.0.or.abs(cstdev).gt.5.0.or.&
        abs(sim(index)).gt.6.0) then
         write(ldbg,104) ix,iy,iz,cmean,cstdev,sim(index)
  104    format('WARNING: grid node location: ',3i5,/,&
               '         conditional mean:   ',f12.5,/,&
               '         conditional stdev:  ',f12.5,/,&
               '         simulated value:    ',f12.5)
      endif

      end
 
      subroutine sgsim_post(sim,sim_sz, scratch_c,c_sz,c_idx,c_idx_siz,&
           scratch_r,real_sz,real_idx,r_idx_siz,scratch_i,int_sz,&
           int_idx,i_idx_siz)
!-----------------------------------------------------------------------
!
!    SGSIM postprocessing: for a single instance. 
!
!-----------------------------------------------------------------------


      use       geostat2
      implicit none
      include  'sgsim2.inc'

      integer   sim_sz,c_sz,real_sz,int_sz
      integer   c_idx_siz, r_idx_siz, i_idx_siz
      integer   c_idx(c_idx_siz),real_idx(r_idx_siz),int_idx(i_idx_siz)
      real*8    sim(sim_sz)

      integer,  target :: scratch_i(int_sz) 
      real*8,   target :: scratch_c(c_sz)
      real*8,   target :: scratch_r(real_sz)

      logical   testindx, testindy, testindz
      integer   ne,ind,id,ix,iy,iz
      real*8    ss,av,xx,yy,zz
      real*8    backtr,simval,TINY
      double precision getCCloc

      nx       => scratch_i(int_idx(1))
      ny       => scratch_i(int_idx(2))
      nz       => scratch_i(int_idx(3))
      dlo      => scratch_i(int_idx(14):int_idx(14)+2)
      dhi      => scratch_i(int_idx(14)+3:int_idx(14)+5)

      xmn      => scratch_r(real_idx(1))
      ymn      => scratch_r(real_idx(2))
      zmn      => scratch_r(real_idx(3))
      xsiz     => scratch_r(real_idx(4))
      ysiz     => scratch_r(real_idx(5))
      zsiz     => scratch_r(real_idx(6))
      if (ktype .ge. 2) then
         lvm      => scratch_r(real_idx(13):real_idx(13)+nx*ny*nz-1)
      endif

      x        => scratch_c(c_idx(1):c_idx(2)-1)
      y        => scratch_c(c_idx(2):c_idx(3)-1)
      z        => scratch_c(c_idx(3):c_idx(4)-1)
      vr       => scratch_c(c_idx(4):c_idx(5)-1)
      wt       => scratch_c(c_idx(5):c_idx(6)-1)
      vrtr     => scratch_c(c_idx(6):c_idx(7)-1)
      vrgtr    => scratch_c(c_idx(7):c_idx(8)-1)

!
! Do we need to reassign the data to the grid nodes?
!
      TINY = 0.0001
      if(sstrat.eq.0) then
         do id=1,nd
            call getindxmod(dlo(1),dhi(1),xmn,xsiz,x(id),ix,testindx)
            call getindxmod(dlo(2),dhi(2),ymn,ysiz,y(id),iy,testindy)
            call getindxmod(dlo(3),dhi(3),zmn,zsiz,z(id),iz,testindz)

            if (testindx .and. testindy .and. testindz) then
               call get1Didx(ix,iy,iz,dlo,dhi,ind)
               xx = getCCloc(ix,xmn,dlo(1),xsiz)
               yy = getCCloc(iy,ymn,dlo(2),ysiz)
               zz = getCCloc(iz,zmn,dlo(3),zsiz)

               test=abs(xx-x(id))+abs(yy-y(id))+abs(zz-z(id))
               if(test.le.TINY) sim(ind) = vr(id)
            endif
         end do
      end if
!
! Back transform each value and write results:
!
      ne = 0
      av = 0.d0
      ss = 0.d0
      do ind=1,sim_sz
         simval = sim(ind)
         if(simval.gt.-9.0.and.simval.lt.9.0) then
            ne = ne + 1
            av = av + simval
            ss = ss + simval*simval
         end if
         if(itrans.eq.1.and.simval.gt.(UNEST+EPSLON)) then
            simval = backtr(simval,ntr,vrtr,vrgtr,zmin,&
                            zmax,ltail,ltpar,utail,utpar)
            if(simval.lt.zmin) simval = zmin
            if(simval.gt.zmax) simval = zmax
         end if
         sim(ind) = simval
      end do
      av = av / max(dble(ne),1.d0)
      ss =(ss / max(dble(ne),1.d0)) - av * av
      write(ldbg,111) ne,av,ss
!      write(*,   111) ne,av,ss
 111  format(/,' Realization : number   = ',i8,/,&
              '               mean     = ',f12.4,&
              ' (close to 0.0?)',/,&
              '               variance = ',f12.4,&
              ' (close to gammabar(V,V)? approx. 1.0)',/)

      end

      subroutine sgsim_deallocate()

      use geostat2
      
      deallocate(covtab)

      end
        

      subroutine ctable()
!-----------------------------------------------------------------------
!
!               Establish the Covariance Look up Table
!               **************************************
!
! The idea is to establish a 3-D network that contains the covariance
! value for a range of grid node offsets that should be at as large
! as twice the search radius in each direction.  The reason it has to
! be twice as large as the search radius is because we want to use it
! to compute the data covariance matrix as well as the data-point
! covariance matrix.
!
! Secondly, we want to establish a search for nearby nodes that
! in order of closeness as defined by the variogram.
!
!
!
! INPUT VARIABLES:
!
!   xsiz,ysiz,zsiz  Definition of the grid being considered
!   MAXCTX,Y,Z      Number of blocks in covariance table
!
!   covariance table parameters
!
!
!
! OUTPUT VARIABLES:  covtab()         Covariance table
!
! EXTERNAL REFERENCES:
!
!   sqdist          Computes 3-D anisotropic squared distance
!   sortem          Sorts multiple arrays in ascending order
!   cova3           Computes the covariance according to a 3-D model
!
!
!
!-----------------------------------------------------------------------


      use       geostat
      implicit none
      include  'sgsim.inc'

      integer ix,iy,iz,loc
      integer il,i,j,k,ic,jc,kc
      integer order_tmp(MAXCTX*MAXCTY*MAXCTZ)
      real*8  hsqd,sqdist,TINY
      real*8  xx,yy,zz
      real*8  tmp(MAXCTX*MAXCTY*MAXCTZ)

      parameter(TINY=1.0d-10)

!
! Size of the look-up table:
!
      nctx = min(((MAXCTX-1)/2),(nx-1))
      ncty = min(((MAXCTY-1)/2),(ny-1))
      nctz = min(((MAXCTZ-1)/2),(nz-1))
!
! Debugging output:
!
      write(ldbg,*)
      write(ldbg,*) 'Covariance Look up table and search for previously'
      write(ldbg,*) 'simulated grid nodes.  The maximum range in each '
      write(ldbg,*) 'coordinate direction for covariance look up is:'
      write(ldbg,*) '          X direction: ',nctx*xsiz
      write(ldbg,*) '          Y direction: ',ncty*ysiz
      write(ldbg,*) '          Z direction: ',nctz*zsiz
      write(ldbg,*) 'Node Values are not searched beyond this distance!'
      write(ldbg,*)
!
! NOTE: If dynamically allocating memory, and if there is no shortage
!       it would a good idea to go at least as far as the radius and
!       twice that far if you wanted to be sure that all covariances
!       in the left hand covariance matrix are within the table look-up.
!
! Initialize the covariance subroutine and cbb at the same time:
!
      call cova3m(0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,1,nst,MAXNST,c0,it,cc,aa,&
                1,MAXROT,rotmat,cmax,cbb)
!
! Now, set up the table and keep track of the node offsets that are
! within the search radius:
!
      nlooku = 0
      do i=-nctx,nctx
      xx = i * xsiz
      ic = nctx + 1 + i
      do j=-ncty,ncty
      yy = j * ysiz
      jc = ncty + 1 + j
      do k=-nctz,nctz
      zz = k * zsiz
      kc = nctz + 1 + k
            call cova3m(0.d0,0.d0,0.d0,xx,yy,zz,1,nst,MAXNST,c0,it,cc,aa,&
                      1,MAXROT,rotmat,cmax,covtab(ic,jc,kc))
            hsqd = sqdist(0.d0,0.d0,0.d0,xx,yy,zz,isrot,MAXROT,rotmat)
            
            if(hsqd.le.radsqd) then
                  nlooku         = nlooku + 1
!
! We want to search by closest variogram distance (and use the
! anisotropic Euclidean distance to break ties:
!
                  tmp(nlooku)   = - covtab(ic,jc,kc)-TINY*hsqd
                  order_tmp(nlooku)=((kc-1)*MAXCXY+(jc-1)*MAXCTX+ic)
            endif
      end do
      end do
      end do
!
! Finished setting up the look-up table, now order the nodes such
! that the closest ones, according to variogram distance, are searched
! first. Note: the "loc" array is used because I didn't want to make
! special allowance for 2 byte integers in the sorting subroutine:
!
!      call sortem(1,nlooku,tmp,1,order_tmp,c,d,e,f,g,h)
      call sortem(1,nlooku,tmp,1,order_tmp)
      do il=1,nlooku
            loc = order_tmp(il)
            iz  = int((loc-1)/MAXCXY) + 1
            iy  = int((loc-(iz-1)*MAXCXY-1)/MAXCTX) + 1
            ix  = loc-(iz-1)*MAXCXY - (iy-1)*MAXCTX
            iznode(il) = iz
            iynode(il) = iy
            ixnode(il) = ix
      end do
      if(nodmax.gt.MAXNOD) then
            write(ldbg,*)
            write(ldbg,*) 'The maximum number of close nodes = ',nodmax
            write(ldbg,*) 'this was reset from your specification due '
            write(ldbg,*) 'to storage limitations.'
            nodmax = MAXNOD
      endif
!
! Debugging output if requested:
!
      if(idbg.lt.2) return
      write(ldbg,*)
      write(ldbg,*) 'There are ',nlooku,' nearby nodes that will be '
      write(ldbg,*) 'checked until enough close data are found.'
      write(ldbg,*)
      if(idbg.lt.14) return
      do i=1,nlooku
            xx = dble(ixnode(i) - nctx - 1) * xsiz
            yy = dble(iynode(i) - ncty - 1) * ysiz
            zz = dble(iznode(i) - nctz - 1) * zsiz
            write(ldbg,100) i,xx,yy,zz
      end do
 100  format('Point ',i3,' at ',3f12.4)
!
! All finished:
!
      return
      end



      subroutine srchnd(sim,ix,iy,iz)
!-----------------------------------------------------------------------
!
!               Search for nearby Simulated Grid nodes
!               **************************************
!
! The idea is to spiral away from the node being simulated and note all
! the nearby nodes that have been simulated.
!
!
!
! INPUT VARIABLES:
!
!   ix,iy,iz        index of the point currently being simulated
!   sim             the realization so far
!   nodmax          the maximum number of nodes that we want
!   nlooku          the number of nodes in the look up table
!   i[x,y,z]node    the relative indices of those nodes.
!   [x,y,z]mn       the origin of the global grid netwrok
!   [x,y,z]siz      the spacing of the grid nodes.
!
!
!
! OUTPUT VARIABLES:
!
!   ncnode          the number of close nodes
!   icnode()        the number in the look up table
!   cnode[x,y,z]()  the location of the nodes
!   cnodev()        the values at the nodes
!
!
!
!-----------------------------------------------------------------------


      use geostat
      implicit none
      include  'sgsim.inc'
      integer   ninoct(8)

      integer ix,iy,iz
      real*8  sim(nxyz)

      integer i,j,k,il,ind,idx,idy,idz,iq
!
! Consider all the nearby nodes until enough have been found:
!
      ncnode = 0
      if(noct.gt.0) then
         do i=1,8
            ninoct(i) = 0
         end do
      end if
      do 2 il=2,nlooku
            if(ncnode.eq.nodmax) return
            i = ix + (int(ixnode(il))-nctx-1)
            j = iy + (int(iynode(il))-ncty-1)
            k = iz + (int(iznode(il))-nctz-1)
            if(i.lt. 1.or.j.lt. 1.or.k.lt. 1) go to 2
            if(i.gt.nx.or.j.gt.ny.or.k.gt.nz) go to 2
            ind = i + (j-1)*nx + (k-1)*nxy
            if(sim(ind).gt.UNEST) then
!
! Check the number of data already taken from this octant:
!
                  if(noct.gt.0) then
                        idx = ix - i
                        idy = iy - j
                        idz = iz - k
                        if(idz.gt.0) then
                              iq = 4
                              if(idx.le.0 .and. idy.gt.0) iq = 1
                              if(idx.gt.0 .and. idy.ge.0) iq = 2
                              if(idx.lt.0 .and. idy.le.0) iq = 3
                        else
                              iq = 8
                              if(idx.le.0 .and. idy.gt.0) iq = 5
                              if(idx.gt.0 .and. idy.ge.0) iq = 6
                              if(idx.lt.0 .and. idy.le.0) iq = 7
                        end if
                        ninoct(iq) = ninoct(iq) + 1
                        if(ninoct(iq).gt.noct) go to 2
                  end if
                  ncnode = ncnode + 1
                  icnode(ncnode) = il
                  cnodex(ncnode) = xmn + dble(i-1)*xsiz
                  cnodey(ncnode) = ymn + dble(j-1)*ysiz
                  cnodez(ncnode) = zmn + dble(k-1)*zsiz
                  cnodev(ncnode) = sim(ind)
            endif
 2    continue
!
! Return to !alling program:
!
      return
      end



      subroutine krige(ix,iy,iz,xx,yy,zz,lktype,gmean,cmean,cstdev)
!-----------------------------------------------------------------------
!
!            Builds and Solves the SK or OK Kriging System
!            *********************************************
!
! INPUT VARIABLES:
!
!   ix,iy,iz        index of the point currently being simulated
!   xx,yy,zz        location of the point currently being simulated
!
!
!
! OUTPUT VARIABLES:
!
!   cmean           kriged estimate
!   cstdev          kriged standard deviation
!
!
!
! EXTERNAL REFERENCES: ksol   Gaussian elimination system solution
!
!
!
!-----------------------------------------------------------------------

      use      geostat
      implicit none
      include 'sgsim.inc'

      logical first
      
      integer ix,iy,iz
      integer lktype
      real*8  gmean, cmean, cstdev
      real*8  xx,yy,zz

      integer na,neq,ind,in,i,j,ii,jj,kk,ising,ie,is
      integer ix1,ix2,iy1,iy2,iz1,iz2,index
      real*8  x1,y1,z1,x2,y2,z2
      real*8  edmax,edmin,sfmin,sfmax
      real*8  cov,sumwts
      real*8  vra(MAXKR1),vrea(MAXKR1)

!
! Size of the kriging system:
!
      first = .false.
      na    = nclose + ncnode
 33   continue      
      if(lktype.eq.0) neq = na
      if(lktype.eq.1) neq = na + 1
      if(lktype.eq.2) neq = na
      if(lktype.eq.3) neq = na + 2
      if(lktype.eq.4) neq = na + 1
      if(lktype.ge.3) then
            ind = ix + (iy-1)*nx + (iz-1)*nxy
            if(lvm(ind).le.-6.d0.or.lvm(ind).ge.6.d0) then
                  lktype = 0
                  go to 33
            end if      
      end if
!
! Set up kriging matri!es:
!
      in=0
      do j=1,na
!
! Sort out the actual location of point "j"
!
            if(j.le.nclose) then
                  index  = int(close(j))
                  x1     = x(index)
                  y1     = y(index)
                  z1     = z(index)
                  vra(j) = vr(index)
                  vrea(j)= sec(index)
                  if(lktype.eq.2) vra(j) = vra(j) - vrea(j)
            else
!
! It is a previously simulated node (keep index for table look-up):
!
                  index  = j-nclose
                  x1     = cnodex(index)
                  y1     = cnodey(index)
                  z1     = cnodez(index)
                  vra(j) = cnodev(index)
                  ind    = icnode(index)
                  ix1    = ix + (int(ixnode(ind))-nctx-1)
                  iy1    = iy + (int(iynode(ind))-ncty-1)
                  iz1    = iz + (int(iznode(ind))-nctz-1)
                  index  = ix1 + (iy1-1)*nx + (iz1-1)*nxy
                  if (ktype.eq.2) then
                     vrea(j)= lvm(index)
                     vra(j) = vra(j) - vrea(j)
                  endif
            endif
            do i=1,j
!
! Sort out the actual lo cation of point "i"
!
                  if(i.le.nclose) then
                        index  = int(close(i))
                        x2     = x(index)
                        y2     = y(index)
                        z2     = z(index)
                  else
!
! It is a previously simulated node (keep index for table look-up):
!
                        index  = i-nclose
                        x2     = cnodex(index)
                        y2     = cnodey(index)
                        z2     = cnodez(index)
                        ind    = icnode(index)
                        ix2    = ix + (int(ixnode(ind))-nctx-1)
                        iy2    = iy + (int(iynode(ind))-ncty-1)
                        iz2    = iz + (int(iznode(ind))-nctz-1)
                  endif
!
! Now, get the covariance value:
!
                  in = in + 1
!
! Decide whether or not to use the covariance look-up table:
!
                  if(j.le.nclose.or.i.le.nclose) then
                        call cova3m(x1,y1,z1,x2,y2,z2,1,nst,MAXNST,c0,it,&
                                  cc,aa,1,MAXROT,rotmat,cmax,cov)
                        a(in) = cov
                  else
!
! Try to use the covariance look-up (if the distance is in range):
!
                        ii = nctx + 1 + (ix1 - ix2)
                        jj = ncty + 1 + (iy1 - iy2)
                        kk = nctz + 1 + (iz1 - iz2)
                        if(ii.lt.1.or.ii.gt.MAXCTX.or.&
                          jj.lt.1.or.jj.gt.MAXCTY.or.&
                          kk.lt.1.or.kk.gt.MAXCTZ) then
                              call cova3m(x1,y1,z1,x2,y2,z2,1,nst,MAXNST,&
                                  c0,it,cc,aa,1,MAXROT,rotmat,cmax,cov)
                        else
                              cov = covtab(ii,jj,kk)
                        endif
                        a(in) = cov
                  endif
            end do
!
! Get the RHS value (possibly with covariance look-up table):
!
            if(j.le.nclose) then
                  call cova3m(xx,yy,zz,x1,y1,z1,1,nst,MAXNST,c0,it,cc,aa,&
                            1,MAXROT,rotmat,cmax,cov)
                  r(j) = cov
            else
!
! Try to use the covariance look-up (if the distance is in range):
!
                  ii = nctx + 1 + (ix - ix1)
                  jj = ncty + 1 + (iy - iy1)
                  kk = nctz + 1 + (iz - iz1)
                  if(ii.lt.1.or.ii.gt.MAXCTX.or.&
                    jj.lt.1.or.jj.gt.MAXCTY.or.&
                    kk.lt.1.or.kk.gt.MAXCTZ) then
                        call cova3m(xx,yy,zz,x1,y1,z1,1,nst,MAXNST,c0,it,&
                                  cc,aa,1,MAXROT,rotmat,cmax,cov)
                  else
                        cov = covtab(ii,jj,kk)
                  endif
                  r(j) = cov
            endif
            rr(j) = r(j)
      end do
!
! Addition of OK constraint:
!
      if(lktype.eq.1.or.lktype.eq.3) then
            do i=1,na
                  in    = in + 1
                  a(in) = 1.d0
            end do
            in       = in + 1
            a(in)    = 0.d0
            r(na+1)  = 1.d0
            rr(na+1) = 1.d0
      endif
!
! Addition of the External Drift Constraint:
!
      if(lktype.eq.3) then
            edmin =  9.99999d5
            edmax = -9.99999d5
            do i=1,na
                  in    = in + 1
                  a(in) = vrea(i)
                  if(a(in).lt.edmin) edmin = a(in)
                  if(a(in).gt.edmax) edmax = a(in)
            end do
            in       = in + 1
            a(in)    = 0.d0
            in       = in + 1
            a(in)    = 0.d0
            ind      = ix + (iy-1)*nx + (iz-1)*nxy
            r(na+2)  = lvm(ind)
            rr(na+2) = r(na+2)
            if((edmax-edmin).lt.EPSLON) neq = neq - 1
      endif
!
! Addition of Collocated Cosimulation Constraint:
!
      if(lktype.eq.4) then
            sfmin =  1.0d21
            sfmax = -1.0d21
            do i=1,na
                  in    = in + 1
                  a(in) = dble(colocorr)*r(i)
                  if(a(in).lt.sfmin) sfmin = a(in)
                  if(a(in).gt.sfmax) sfmax = a(in)
            end do
            in    = in + 1
            a(in) = 1.d0
            ii    = na + 1
            r(ii) = dble(colocorr)
            rr(ii)= r(ii)
!           if((sfmax-sfmin).lt.EPSLON) neq = neq - 1
      end if
!
! Write out the kriging Matrix if Seriously Debugging:
!
      if(idbg.ge.3) then
            write(ldbg,100) ix,iy,iz
            is = 1
            do i=1,neq
                  ie = is + i - 1
                  write(ldbg,101) i,r(i),(a(j),j=is,ie)
                  is = is + i
            end do
 100        format(/,'Kriging Matrices for Node: ',3i4,' RHS first')
 101        format('    r(',i2,') =',f7.4,'  a= ',99f7.4)
      endif
!
! Solve the Kriging System:
!
      if(neq.eq.1.and.lktype.ne.3) then
            s(1)  = r(1) / a(1)
            ising = 0
      else
            call ksol(1,neq,1,a,r,s,ising)
      endif
!
! Write a warning if the matrix is singular:
!
      if(ising.ne.0) then
            if(idbg.ge.1) then
                  write(ldbg,*) 'WARNING SGSIM: singular matrix'
                  write(ldbg,*) '               for node',ix,iy,iz
            endif
            cmean  = gmean
            cstdev = 1.d0
            return
      endif
!
! Compute the estimate and kriging variance.  Recall that kriging type
!     0 = Simple Kriging:
!     1 = Ordinary Kriging:
!     2 = Locally Varying Mean:
!     3 = External Drift:
!     4 = Collocated Cosimulation:
!
      cmean  = 0.d0
      cstdev = cbb
      sumwts = 0.d0
      do i=1,na
            cmean  = cmean  + dble(s(i))*vra(i)
            cstdev = cstdev - dble(s(i)*rr(i))
            sumwts = sumwts + dble(s(i))
      end do

      if(lktype.eq.1) cstdev = cstdev - dble(s(na+1))

      if(lktype.eq.2) cmean  = cmean + gmean

      if(lktype.eq.4) then
            ind    = ix + (iy-1)*nx + (iz-1)*nxy
            cmean  = cmean  + dble(s(na+1))*lvm(ind)
            cstdev = cstdev - dble(s(na+1) *rr(na+1))
      end if
!
! Error message if negative variance:
!
      if(cstdev.lt.0.d0) then
            write(ldbg,*) 'ERROR: Negative Variance: ',cstdev
            cstdev = 0.d0
      endif
      cstdev = dsqrt(max(cstdev,0.d0))
!
! Write out the kriging Weights if Seriously Debugging:
!
      if(idbg.ge.3) then
            do i=1,na
                  write(ldbg,140) i,vra(i),s(i)
            end do
 140        format(' Data ',i4,' value ',f8.4,' weight ',f8.4)
            if(lktype.eq.4) write(ldbg,141) lvm(ind),s(na+1)
 141        format(' Sec Data  value ',f8.4,' weight ',f8.4)
            write(ldbg,142) gmean,cmean,cstdev
 142        format(' Global mean ',f8.4,' conditional ',f8.4,&
                   ' std dev ',f8.4)
      end if


!
! Finished Here:
!
      return
      end


      subroutine ctable2()
!-----------------------------------------------------------------------
!
!               Establish the covariance Look up Table
!               **************************************
!
! The idea is to establish a 3-D network that contains the covariance
! value for a range of grid node offsets that should be at as large
! as twice the search radius in each direction.  The reason it has to
! be twice as large as the search radius is because we want to use it
! to compute the data covariance matrix as well as the data-point
! covariance matrix.
!
! Secondly, we want to establish a search for nearby nodes that
! in order of closeness as defined by the variogram.
!
!
!
! INPUT VARIABLES:
!
!   xsiz,ysiz,zsiz  Definition of the grid being considered
!   MAXCTX,Y,Z      Number of blocks in covariance table
!
!   covariance table parameters
!
!
!
! OUTPUT VARIABLES:  covtab()         Covariance table
!
! EXTERNAL REFERENCES:
!
!   sqdist          Computes 3-D anisotropic squared distance
!   sortem          Sorts multiple arrays in ascending order
!   cova3           Computes the covariance according to a 3-D model
!
!
!
!-----------------------------------------------------------------------


      use       geostat2
      implicit none
      include  'sgsim2.inc'


      integer ix,iy,iz,loc
      integer il,i,j,k,ic,jc,kc
      integer order_tmp(MAXCTX*MAXCTY*MAXCTZ)
      real*8  hsqd,sqdist,TINY
      real*8  xx,yy,zz
      real*8  tmp(MAXCTX*MAXCTY*MAXCTZ)

      parameter(TINY=1.0d-10)

!
! Size of the look-up table:
!
      nctx = min(((MAXCTX-1)/2),(nx-1))
      ncty = min(((MAXCTY-1)/2),(ny-1))
      nctz = min(((MAXCTZ-1)/2),(nz-1))

!
! Debugging output:
!
      write(ldbg,*)
      write(ldbg,*) 'Covariance Look up table and search for previously'
      write(ldbg,*) 'simulated grid nodes.  The maximum range in each '
      write(ldbg,*) 'coordinate direction for covariance look up is:'
      write(ldbg,*) '          X direction: ',nctx*xsiz
      write(ldbg,*) '          Y direction: ',ncty*ysiz
      write(ldbg,*) '          Z direction: ',nctz*zsiz
      write(ldbg,*) 'Node Values are not searched beyond this distance!'
      write(ldbg,*)
!
! NOTE: If dynamically allocating memory, and if there is no shortage
!       it would a good idea to go at least as far as the radius and
!       twice that far if you wanted to be sure that all covariances
!       in the left hand covariance matrix are within the table look-up.
!
! Initialize the covariance subroutine and cbb at the same time:
!
      call cova3m(0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,1,nst,MAXNST,c0,it,cc,aa,&
                 1,MAXROT,rotmat,cmax,cbb)
!
! Now, set up the table and keep track of the node offsets that are
! within the search radius:
!
      nlooku = 0
      do i=-nctx,nctx
      xx = i * xsiz
      ic = nctx + 1 + i
      do j=-ncty,ncty
      yy = j * ysiz
      jc = ncty + 1 + j
      do k=-nctz,nctz
      zz = k * zsiz
      kc = nctz + 1 + k
            call cova3m(0.d0,0.d0,0.d0,xx,yy,zz,1,nst,MAXNST,c0,it,cc,aa,&
                      1,MAXROT,rotmat,cmax,covtab(ic,jc,kc))
            hsqd = sqdist(0.d0,0.d0,0.d0,xx,yy,zz,isrot,MAXROT,rotmat)
            
            if(hsqd.le.radsqd) then
                  nlooku         = nlooku + 1
!
! We want to search by closest variogram distance (and use the
! anisotropic Euclidean distance to break ties:
!
                  tmp(nlooku)   = - covtab(ic,jc,kc)-TINY*hsqd
                  order_tmp(nlooku)=((kc-1)*MAXCXY+(jc-1)*MAXCTX+ic)
            endif
      end do
      end do
      end do
!
! Finished setting up the look-up table, now order the nodes such
! that the closest ones, according to variogram distance, are searched
! first. Note: the "loc" array is used because I didn't want to make
! special allowance for 2 byte integers in the sorting subroutine:
!
      call sortem(1,nlooku,tmp,1,order_tmp)
      do il=1,nlooku
            loc = order_tmp(il)
            iz  = int((loc-1)/MAXCXY) + 1
            iy  = int((loc-(iz-1)*MAXCXY-1)/MAXCTX) + 1
            ix  = loc-(iz-1)*MAXCXY - (iy-1)*MAXCTX
            iznode(il) = iz
            iynode(il) = iy
            ixnode(il) = ix
      end do
      if(nodmax.gt.MAXNOD) then
            write(ldbg,*)
            write(ldbg,*) 'The maximum number of close nodes = ',nodmax
            write(ldbg,*) 'this was reset from your specification due '
            write(ldbg,*) 'to storage limitations.'
            nodmax = MAXNOD
      endif
!
! Debugging output if requested:
!
      if(idbg.lt.2) return
      write(ldbg,*)
      write(ldbg,*) 'There are ',nlooku,' nearby nodes that will be '
      write(ldbg,*) 'checked until enough close data are found.'
      write(ldbg,*)
      if(idbg.lt.14) return
      do i=1,nlooku
            xx = dble(ixnode(i) - nctx - 1) * xsiz
            yy = dble(iynode(i) - ncty - 1) * ysiz
            zz = dble(iznode(i) - nctz - 1) * zsiz
            write(ldbg,100) i,xx,yy,zz
      end do
 100  format('Point ',i3,' at ',3f12.4)
!
! All finished:
!
      return
      end

      subroutine srchnd2(icnode,cnodex,cnodey,cnodez,cnodev,&
                        sim,sim_sz,ix,iy,iz)
!-----------------------------------------------------------------------
!
!               Search for nearby Simulated Grid nodes
!               **************************************
!
! The idea is to spiral away from the node being simulated and note all
! the nearby nodes that have been simulated.
!
!
!
! INPUT VARIABLES:
!
!   ix,iy,iz        index of the point currently being simulated
!   sim             the realization so far
!   nodmax          the maximum number of nodes that we want
!   nlooku          the number of nodes in the look up table
!   i[x,y,z]node    the relative indices of those nodes.
!   [x,y,z]mn       the origin of the global grid netwrok
!   [x,y,z]siz      the spacing of the grid nodes.
!
!
!
! OUTPUT VARIABLES:
!
!   ncnode          the number of close nodes
!   icnode()        the number in the look up table
!   cnode[x,y,z]()  the location of the nodes
!   cnodev()        the values at the nodes
!
!
!
!-----------------------------------------------------------------------

      use geostat2
      implicit none
      include 'sgsim2.inc'

      integer sim_sz,ix,iy,iz
      integer icnode(MAXNOD)
      real*8  cnodex(MAXNOD),cnodey(MAXNOD),cnodez(MAXNOD)
      real*8  cnodev(MAXNOD)
      real*8  sim(sim_sz)

      integer i,j,k,il,ind,idx,idy,idz,iq
      integer ninoct(8)
      
!
! Consider all the nearby nodes until enough have been found:
!
      ncnode = 0
      if(noct.gt.0) then
         do il=1,8
            ninoct(il) = 0
         end do
      end if
      do 2 il=2,nlooku
            if(ncnode.eq.nodmax) return
            i = ix + (ixnode(il)-nctx-1)
            j = iy + (iynode(il)-ncty-1)
            k = iz + (iznode(il)-nctz-1)

            if (i.lt.dlo(1).or.j.lt.dlo(2).or.k.lt.dlo(3)) goto 2
            if (i.gt.dhi(1).or.j.gt.dhi(2).or.k.gt.dhi(3)) goto 2

            call get1Didx(i,j,k,dlo,dhi,ind)
            if (ind.gt.sim_sz) then
               print *,'bust in srchnd2',ind,sim_sz
            endif


            if(sim(ind).gt.UNEST) then
!
! Check the number of data already taken from this octant:
!
                  if(noct.gt.0) then
                        idx = ix - i
                        idy = iy - j
                        idz = iz - k
                        if(idz.gt.0) then
                           iq = 4
                           if(idx.le.0 .and. idy.gt.0) iq = 1
                           if(idx.gt.0 .and. idy.ge.0) iq = 2
                           if(idx.lt.0 .and. idy.le.0) iq = 3
                        else
                           iq = 8
                           if(idx.le.0 .and. idy.gt.0) iq = 5
                           if(idx.gt.0 .and. idy.ge.0) iq = 6
                           if(idx.lt.0 .and. idy.le.0) iq = 7
                        end if
                        ninoct(iq) = ninoct(iq) + 1
                        if(ninoct(iq).gt.noct) go to 2
                  end if
                  ncnode = ncnode + 1
                  icnode(ncnode) = il
                  cnodex(ncnode) = xmn + dble(i-1)*xsiz
                  cnodey(ncnode) = ymn + dble(j-1)*ysiz
                  cnodez(ncnode) = zmn + dble(k-1)*zsiz
                  cnodev(ncnode) = sim(ind)
            endif
 2    continue
!
! Return to calling program:
!
      return
      end



      subroutine krige2(icnode,cnodex,cnodey,cnodez,cnodev,&
                       ix,iy,iz,xx,yy,zz,&
                       lktype,gmean,cmean,cstdev)
!-----------------------------------------------------------------------
!
!            Builds and Solves the SK or OK Kriging System
!            *********************************************
!
! INPUT VARIABLES:
!
!   ix,iy,iz        index of the point currently being simulated
!   xx,yy,zz        location of the point currently being simulated
!
!
!
! OUTPUT VARIABLES:
!
!   cmean           kriged estimate
!   cstdev          kriged standard deviation
!
!
!
! EXTERNAL REFERENCES: ksol   Gaussian elimination system solution
!
!
!
!-----------------------------------------------------------------------

      use      geostat2
      implicit none
      include 'sgsim2.inc'

      integer ix,iy,iz
      integer lktype
      integer icnode(MAXNOD)
      real*8  xx,yy,zz
      real*8  gmean, cmean, cstdev
      real*8  cnodex(MAXNOD),cnodey(MAXNOD),cnodez(MAXNOD)
      real*8  cnodev(MAXNOD)     
 
      integer na,neq,ind,in,i,j,ii,jj,kk,ising,ie,is
      integer ix1,ix2,iy1,iy2,iz1,iz2,index
      real*8  x1,y1,z1,x2,y2,z2
      real*8  edmax,edmin,sfmin,sfmax
      real*8  cov,sumwts
      real*8  vra(MAXKR1),vrea(MAXKR1)
      real*8  r(MAXKR1),rr(MAXKR1),s(MAXKR1),a(MAXKR1*MAXKR1)

      logical first

!
! Size of the kriging system:
!
      first = .false.
      na    = nclose + ncnode
 33   continue      
      if(lktype.eq.0) neq = na
      if(lktype.eq.1) neq = na + 1
      if(lktype.eq.2) neq = na
      if(lktype.eq.3) neq = na + 2
      if(lktype.eq.4) neq = na + 1
      if(lktype.eq.5) neq = na + 1
      if(lktype.ge.3) then
         call get1Didx(ix,iy,iz,dlo,dhi,ind)
         !ind = ix + (iy-1)*nx + (iz-1)*nxy
         if(lvm(ind).le.-6.0.or.lvm(ind).ge.6.0) then
            lktype = 0
            go to 33
         end if      
      end if

!
! Set up kriging matrices:
!
      in=0
      do j=1,na
!
! Sort out the actual location of point "j"
!
         if(j.le.nclose) then
            index  = int(close(j))
            x1     = x(index)
            y1     = y(index)
            z1     = z(index)
            vra(j) = vr(index)
            vrea(j)= sec(index)
            if(lktype.eq.2) vra(j) = vra(j) - vrea(j)
         else
!
! It is a previously simulated node (keep index for table look-up):
!
            index  = j-nclose
            x1     = cnodex(index)
            y1     = cnodey(index)
            z1     = cnodez(index)
            vra(j) = cnodev(index)
            ind    = icnode(index)
            ix1    = ix + (int(ixnode(ind))-nctx-1)
            iy1    = iy + (int(iynode(ind))-ncty-1)
            iz1    = iz + (int(iznode(ind))-nctz-1)
            call get1Didx(ix1,iy1,iz1,dlo,dhi,index)
            !index  = ix1 + (iy1-1)*nx + (iz1-1)*nxy
            if (lktype.eq.2) then
               vrea(j)= lvm(index)
               vra(j) = vra(j) - vrea(j)
            endif
         endif
         do i=1,j
!
! Sort out the actual location of point "i"
!
            if(i.le.nclose) then
               index  = int(close(i))
               x2     = x(index)
               y2     = y(index)
               z2     = z(index)
            else
!
! It is a previously simulated node (keep index for table look-up):
!
               index  = i-nclose
               x2     = cnodex(index)
               y2     = cnodey(index)
               z2     = cnodez(index)
               ind    = icnode(index)
               ix2    = ix + (int(ixnode(ind))-nctx-1)
               iy2    = iy + (int(iynode(ind))-ncty-1)
               iz2    = iz + (int(iznode(ind))-nctz-1)
            endif
!
! Now, get the covariance value:
!
             in = in + 1
!
! Decide whether or not to use the covariance look-up table:
!


             if(j.le.nclose.or.i.le.nclose) then
                call cova3m(x1,y1,z1,x2,y2,z2,1,nst,MAXNST,c0,it,&
                          cc,aa,1,MAXROT,rotmat,cmax,cov)
                a(in) = dble(cov)
             else
!
! Try to use the covariance look-up (if the distance is in range):
!
                ii = nctx + 1 + (ix1 - ix2)
                jj = ncty + 1 + (iy1 - iy2)
                kk = nctz + 1 + (iz1 - iz2)
                if(ii.lt.1.or.ii.gt.MAXCTX.or.&
                  jj.lt.1.or.jj.gt.MAXCTY.or.&
                  kk.lt.1.or.kk.gt.MAXCTZ) then
                   call cova3m(x1,y1,z1,x2,y2,z2,1,nst,MAXNST,&
                             c0,it,cc,aa,1,MAXROT,rotmat,cmax,cov)
                else
                   cov = covtab(ii,jj,kk)
                endif
                a(in) = dble(cov)
             endif
          end do
!
! Get the RHS value (possibly with covariance look-up table):
!
          if(j.le.nclose) then
             call cova3m(xx,yy,zz,x1,y1,z1,1,nst,MAXNST,c0,it,cc,aa,&
                       1,MAXROT,rotmat,cmax,cov)
             r(j) = dble(cov)
          else
!
! Try to use the covariance look-up (if the distance is in range):
!
             ii = nctx + 1 + (ix - ix1)
             jj = ncty + 1 + (iy - iy1)
             kk = nctz + 1 + (iz - iz1)

             if(ii.lt.1.or.ii.gt.MAXCTX.or.&
               jj.lt.1.or.jj.gt.MAXCTY.or.&
               kk.lt.1.or.kk.gt.MAXCTZ) then
                call cova3m(xx,yy,zz,x1,y1,z1,1,nst,MAXNST,c0,it,&
                          cc,aa,1,MAXROT,rotmat,cmax,cov)
             else
                cov = covtab(ii,jj,kk)
             endif
             r(j) = cov
          endif
          rr(j) = r(j)
       end do
!
! Addition of OK constraint:
!
      if(lktype.eq.1.or.lktype.eq.3) then
         do i=1,na
            in    = in + 1
            a(in) = 1.d0
         end do
         in       = in + 1
         a(in)    = 0.d0
         r(na+1)  = 1.d0
         rr(na+1) = 1.d0
      endif
!
! Addition of mean constraint:
!
      if (lktype.eq.5) then
        
         do i = 1,na
            in = in + 1
            a(in) = - vra(i)/dble(na+1)

         end do
         in = in + 1
         a(in) = -2*dble(na)
         r(na+1) = -gmean
         do i = 1,na
            r(na+1) = r(na+1) + vra(i)/dble(na+1) 
         end do
!         r(na+1) = 10.d0/na
         rr(na+1) = r(na+1)
!          write(*,*) "na = ", na, r(na+1) 
      end if
!
! Addition of the External Drift Constraint:
!
      if(lktype.eq.3) then
         edmin =  999999.
         edmax = -999999.
         do i=1,na
            in    = in + 1
            a(in) = vrea(i)
            if(a(in).lt.edmin) edmin = a(in)
            if(a(in).gt.edmax) edmax = a(in)
         end do
         in       = in + 1
         a(in)    = 0.d0
         in       = in + 1
         a(in)    = 0.d0
         call get1Didx(ix,iy,iz,dlo,dhi,ind)
         !ind      = ix + (iy-1)*nx + (iz-1)*nxy
         r(na+2)  = lvm(ind)
         rr(na+2) = r(na+2)
         if((edmax-edmin).lt.EPSLON) neq = neq - 1
      endif
!
! Addition of Collocated Cosimulation Constraint:
!
      if(lktype.eq.4) then
         sfmin =  1.0d21
         sfmax = -1.0d21
         do i=1,na
            in    = in + 1
            a(in) = dble(colocorr)*r(i)
            if(a(in).lt.sfmin) sfmin = a(in)
            if(a(in).gt.sfmax) sfmax = a(in)
         end do
         in    = in + 1
         a(in) = 1.d0
         ii    = na + 1
         r(ii) = dble(colocorr)
         rr(ii)= r(ii)
!        if((sfmax-sfmin).lt.EPSLON) neq = neq - 1
      end if
!
! Write out the kriging Matrix if Seriously Debugging:
!
      if(idbg.ge.3) then
         write(ldbg,100) ix,iy,iz
         is = 1
         do i=1,neq
            ie = is + i - 1
            write(ldbg,101) i,r(i),(a(j),j=is,ie)
            is = is + i
         end do
 100     format(/,'Kriging Matrices for Node: ',3i4,' RHS first')
 101     format('    r(',i2,') =',f7.4,'  a= ',99f7.4)
      endif
!
! Solve the Kriging System:
!
      if(neq.eq.1.and.lktype.ne.3) then
         s(1)  = r(1) / a(1)
         ising = 0
      else
         call ksol(1,neq,1,a,r,s,ising)
      endif
!
! Write a warning if the matrix is singular:
!
      if(ising.ne.0) then
         if(idbg.ge.1) then
            write(ldbg,*) 'WARNING SGSIM: singular matrix'
            write(ldbg,*) '               for node',ix,iy,iz
         endif
         cmean  = gmean
         cstdev = 1.d0
         return
      endif
!
! Compute the estimate and kriging variance.  Recall that kriging type
!     0 = Simple Kriging:
!     1 = Ordinary Kriging:
!     2 = Locally Varying Mean:
!     3 = External Drift:
!     4 = Collocated Cosimulation:
!
      cmean  = 0.d0
      cstdev = cbb
      sumwts = 0.d0
      do i=1,na
         cmean  = cmean  + dble(s(i))*vra(i)
         cstdev = cstdev - dble(s(i)*rr(i))
         sumwts = sumwts + dble(s(i))
      end do

      if(lktype.eq.1) cstdev = cstdev - dble(s(na+1))

      if(lktype.eq.2) cmean  = cmean + gmean

      if(lktype.eq.4) then
            call get1Didx(ix,iy,iz,dlo,dhi,ind)
            !ind    = ix + (iy-1)*nx + (iz-1)*nxy
            cmean  = cmean  + dble(s(na+1))*lvm(ind)
            cstdev = cstdev - dble(s(na+1) *rr(na+1))
      end if

!      write(*,*) cmean, cstdev
!
! Error message if negative variance:
!
      if(cstdev.lt.0.d0) then
         write(ldbg,*) 'ERROR: Negative Variance: ',cstdev
         cstdev = 0.d0
      endif
      cstdev = dsqrt(max(cstdev,0.d0))

!
! Write out the kriging Weights if Seriously Debugging:
!
      if(idbg.ge.3) then
         do i=1,na
            write(ldbg,140) i,vra(i),s(i)
         end do
 140     format(' Data ',i4,' value ',f8.4,' weight ',f8.4)
         if(lktype.eq.4) write(ldbg,141) lvm(ind),s(na+1)
 141     format(' Sec Data  value ',f8.4,' weight ',f8.4)
         write(ldbg,142) gmean,cmean,cstdev
 142     format(' Global mean ',f8.4,' conditional ',f8.4,&
                ' std dev ',f8.4)
      end if


!
! Finished Here:
!
      return
      end


      subroutine cova3m(x1,y1,z1,x2,y2,z2,ivarg,nst,MAXNST,c0,it,cc,aa,&
           irot,MAXROT,rotmat,cmax,cova)
! -----------------------------------------------------------------------

!                     Covariance Between Two Points
!                     *****************************

!  This subroutine calculated the covariance associated with a variogram
!  model specified by a nugget effect and nested varigoram structures.
!  The anisotropy definition can be different for each nested structure.



!  INPUT VARIABLES:

!    x1,y1,z1         coordinates of first point
!    x2,y2,z2         coordinates of second point
!    nst(ivarg)       number of nested structures (maximum of 4)
!    ivarg            variogram number (set to 1 unless doing cokriging
!                        or indicator kriging)
!    MAXNST           size of variogram parameter arrays
!    c0(ivarg)        isotropic nugget constant
!    it(i)            type of each nested structure:
!                       1. spherical model of range a;
!                       2. exponential model of parameter a;
!                            i.e. practical range is 3a
!                       3. gaussian model of parameter a;
!                            i.e. practical range is a*sqrt(3)
!                       4. power model of power a (a must be gt. 0  and
!                            lt. 2).  if linear model, a=1,c=slope.
!                       5. hole effect model
!    cc(i)            multiplicative factor of each nested structure.
!                       (sill-c0) for spherical, exponential,and gaussian
!                       slope for linear model.
!    aa(i)            parameter "a" of each nested structure.
!    irot             index of the rotation matrix for the first nested 
!                     structure (the second nested structure will use
!                     irot+1, the third irot+2, and so on)
!    MAXROT           size of rotation matrix arrays
!    rotmat           rotation matrices


!  OUTPUT VARIABLES:

!    cmax             maximum covariance
!    cova             covariance between (x1,y1,z1) and (x2,y2,z2)



!  EXTERNAL REFERENCES: sqdist    computes anisotropic squared distance
!                       rotmat    computes rotation matrix for distance
! -----------------------------------------------------------------------

      implicit none

      double precision PI, PMX, EPSLON
      parameter(PI=3.14159265,PMX=999.,EPSLON=1.e-20)

      double precision x1,y1,z1,x2,y2,z2,cmax,cova
      integer ivarg,MAXNST,irot,MAXROT
      integer   nst(*),it(*)
      double precision c0(*),cc(*),aa(*)
      double precision rotmat(MAXROT,3,3),hsqd,sqdist,h,hr
      integer ir,is,ist,istart

 ! Calculate the maximum covariance value (used for zero distances and
 ! for power model covariance):

      istart = 1 + (ivarg-1)*MAXNST
      cmax   = c0(ivarg)
      do is=1,nst(ivarg)
            ist = istart + is - 1
            if(it(ist).eq.4) then
                  cmax = cmax + PMX
            else
                  cmax = cmax + cc(ist)
            endif
      end do

 ! Check for "zero" distance, return with cmax if so:

      hsqd = sqdist(x1,y1,z1,x2,y2,z2,irot,MAXROT,rotmat)
      if(real(hsqd).lt.EPSLON) then
            cova = cmax
            return
      endif

 ! Loop over all the structures:

      cova = 0.d0
      do is=1,nst(ivarg)
            ist = istart + is - 1

 ! Compute the appropriate distance:

            if(ist.ne.1) then
                  ir = min((irot+is-1),MAXROT)
                  hsqd=sqdist(x1,y1,z1,x2,y2,z2,ir,MAXROT,rotmat)
            end if
            h = real(dsqrt(hsqd))

 ! Spherical Variogram Model?

            if(it(ist).eq.1) then
                  hr = h/aa(ist)
                  if(hr.lt.1.) cova=cova+cc(ist)*(1.-hr*(1.5-.5*hr*hr))

 ! Exponential Variogram Model?

            else if(it(ist).eq.2) then
                  cova = cova + cc(ist)*exp(-3.0*h/aa(ist))

 ! Gaussian Variogram Model?

            else if(it(ist).eq.3) then
                  cova = cova + cc(ist)*exp(-3.*(h/aa(ist))*(h/aa(ist)))

 ! Power Variogram Model?

            else if(it(ist).eq.4) then
                  cova = cova + cmax - cc(ist)*(h**aa(ist))

 ! Hole Effect Model?

            else if(it(ist).eq.5) then
                 ! d = 10.0 * aa(ist)
                 ! cova = cova + cc(ist)*exp(-3.0*h/d)*cos(h/aa(ist)*PI)
                  cova = cova + cc(ist)*cos(h/aa(ist)*PI)
            endif
      end do

 ! Finished:

      return
      end

