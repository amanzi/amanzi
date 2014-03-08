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
!-----------------------------------------------------------------------
!
! The following Parameters control static dimensioning for SGSIM:
!
!
!-----------------------------------------------------------------------
!
! Parameters:
!
      implicit integer (i-n)
      implicit double precision (a-h,o-z)

      integer,parameter :: MAXNST = 4
      integer,parameter :: MAXROT = MAXNST + 1
      real*8,parameter  :: UNEST=-1.0d20
      real*8,parameter  :: EPSLON=1.0d-20
!
! ACORN parameters:
!
      integer,parameter :: KORDEI = 12
      integer,parameter :: MAXOP1 = KORDEI + 1
      integer,parameter :: MAXINT = 2**30
!
! Variable Declaration:
!
      real*8    tmin,tmax,zmin,zmax,varred
      real*8    radius,radsqd,sang1,sang2,sang3,sanis1,sanis2
      real*8    ltpar,utpar,c0(1),cc(MAXNST),aa(MAXNST)
      real*8    ang1(MAXNST),ang2(MAXNST),ang3(MAXNST)
      real*8    anis1(MAXNST),anis2(MAXNST)
      real*8    rotmat(MAXROT,3,3)
      real*8    cbb

      integer   ixv(MAXOP1),do_cond
      integer   nsim,nxy,nxyz,nd,itrans,ntr,idbg,ismooth
      integer   llvm,ldbg,icollvm,nvaril,ktype,colocorr
      integer   noct,nclose,ndmin,ndmax,nmult,mxctx,mxcty,mxctz
      integer   cmax,isrot,nctx,ncty,nctz,nlooku,ncnode,nodmax
      integer   ltail,utail,sstrat,mults,nst(1),it(MAXNST),test
      integer   MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ
      integer   MAXSBX,MAXSBY,MAXSBZ,MAXKR1
      integer   ixl,iyl,izl,ivrl,isecvr,iwt,icolvr,icolwt,isvr,iswt
      character datafl*64, lvmfl*64, smthfl*64, transfl*64

!
! Common Blocks:
!
      common /piaco/   ixv

      common /pgenerl1/ nsim,nxy,nxyz,do_cond,ismooth
      common /pgenerl2/ nd,itrans,ntr,idbg,ldbg,icollvm
      common /pgenerl3/ nvaril,ktype,colocorr,test,ltail,utail
      common /pgenerr/  ltpar,utpar,zmin,zmax,varred

      common /psearch1/ radius,radsqd,sang1,sang2,sang3,sanis1,sanis2
      common /psearch2/ noct,nclose,ndmin,ndmax,sstrat,mults,nmult
      common /psearch3/ mxctx,mxcty,mxctz

      common /pcova3d1/ nst,it,cmax,c0,cc,aa,ang1,ang2,ang3
      common /pcova3d2/ anis1,anis2,isrot

      common /pclooku/ nctx,ncty,nctz,nlooku,ncnode,nodmax

      common /pkrigev/ rotmat,cbb

      common /pMAXMAX1/ MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ
      common /pMAXMAX2/ MAXSBX,MAXSBY,MAXSBZ,MAXKR1
  
      common /pidxcol1/ ixl,iyl,izl,ivrl,isecvr,iwt,icolvr,icolwt
      common /pidxcol2/ isvr,iswt

      common /pflname/ lvmfl, datafl, smthfl, transfl
