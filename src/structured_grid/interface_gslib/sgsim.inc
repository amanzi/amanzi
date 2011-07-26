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
!      implicit integer (i-n)
!      implicit double precision (a-h,o-z)

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
      real*8    xsiz,ysiz,zsiz,xmn,ymn,zmn,zmin,zmax,varred
      real*8    radius,radsqd,sang1,sang2,sang3,sanis1,sanis2
      real*8    ltpar,utpar,c0(1),cc(MAXNST),aa(MAXNST)
      real*8    ang1(MAXNST),ang2(MAXNST),ang3(MAXNST)
      real*8    anis1(MAXNST),anis2(MAXNST)
      real*8    rotmat(MAXROT,3,3)
      real*8    cbb

      integer   ixv(MAXOP1)
      integer   nsim,nx,ny,nz,nxy,nxyz,nd,itrans,ntr,idbg
      integer   lin,lout,ldbg,llvm,icollvm,nvaril,ktype,colocorr
      integer   noct,nclose,ndmin,ndmax,nmult,mxctx,mxcty,mxctz
      integer   cmax,isrot,nctx,ncty,nctz,nlooku,ncnode,nodmax
      integer   ltail,utail,sstrat,mults,nst(1),it(MAXNST),test
      integer   MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ
      integer   MAXSBX,MAXSBY,MAXSBZ
      integer   MAXKR1
!
! Common Blocks:
!
      common /iaco/   ixv

      common /generl1/ nsim,nx,ny,nz,nxy,nxyz,xsiz,ysiz,zsiz,xmn,ymn,zmn
      common /generl2/ nd,itrans,ntr,idbg,lin,lout,ldbg,llvm,icollvm
      common /generl3/ ktype,colocorr,test,ltail,utail 
      common /generr/  ltpar,utpar,zmin,zmax,varred,nvaril

      common /search1/ radius,radsqd,sang1,sang2,sang3,sanis1,sanis2
      common /search2/ noct,nclose,ndmin,ndmax,sstrat,mults,nmult,isrot
      common /search3/ mxctx,mxcty,mxctz

      common /cova3d/ nst,it,cmax,c0,cc,aa,ang1,ang2,ang3,anis1,anis2

      common /clooku/ nctx,ncty,nctz,nlooku,ncnode,nodmax

      common /krigev/ rotmat,cbb

      common /MAXMAX1/ MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ
      common /MAXMAX2/ MAXSBX,MAXSBY,MAXSBZ,MAXKR1
