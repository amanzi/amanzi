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
!
! Module to declare dynamic arrays in multiple subroutines:
!
module geostat

  real*8,allocatable      :: x(:),y(:),z(:),vr(:),wt(:),lvm(:), &
                             vrtr(:),vrgtr(:),close(:),sec(:), & 
                             covtab(:,:,:), &
                             cnodex(:),cnodey(:),cnodez(:),cnodev(:)
  real*8,   allocatable :: r(:),rr(:),s(:),a(:)
  integer,  allocatable :: nisb(:),icnode(:),order(:)
  integer*2,allocatable :: ixnode(:),iynode(:),iznode(:), &
            ixsbtosr(:),iysbtosr(:),izsbtosr(:)

      end module


!
! Module for parallel implementation: 
!
      module geostat2

      integer, pointer :: nx,ny,nz,nxsup,nysup,nzsup
      integer, pointer :: nisb(:)
      integer, pointer :: ixsbtosr(:), iysbtosr(:), izsbtosr(:)
      integer, pointer :: ixnode(:), iynode(:), iznode(:)
      integer, pointer :: order(:)

      real*8,  pointer :: xmn,xsiz,xmnsup,xsizsup
      real*8,  pointer :: ymn,ysiz,ymnsup,ysizsup
      real*8,  pointer :: zmn,zsiz,zmnsup,zsizsup
      real*8,  pointer :: lvm(:)


      real*8,  pointer :: x(:),y(:),z(:),vr(:),wt(:), &
                          vrtr(:),vrgtr(:),close(:),sec(:)
      
      real*8,  allocatable :: covtab(:,:,:)

      end module

