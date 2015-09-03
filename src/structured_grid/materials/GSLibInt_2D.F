#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "GSLibInt_F.H" 
#include "ArrayLim.H"

#define SDIM 2


c :: ----------------------------------------------------------
c :: Read GSLIB parameters
c :: ----------------------------------------------------------
      subroutine FORT_INIT_GSLIB(paramfl,paramsz,nx_in,xsiz_in,xmn_in,rseed)
      implicit none
      integer   paramsz, rseed
      character paramfl(paramsz)

      integer nx_in(3)
      REAL_T  xsiz_in(3),xmn_in(3)

      call readparm(paramfl,nx_in,xsiz_in,xmn_in,rseed)

      end

c :: ----------------------------------------------------------
c :: Read GSLIB parameters version 2
c :: ----------------------------------------------------------
      subroutine FORT_INIT_GSLIB2(paramfl,paramsz,c_sz,c_idx,
     &     c_idx_siz,real_sz,real_idx,r_idx_siz,int_sz,int_idx,
     &     i_idx_siz,cond_option,rseed)
      implicit none
      integer   paramsz,cond_option,rseed
      integer   c_sz,real_sz,int_sz
      integer   c_idx_siz,r_idx_siz,i_idx_siz
      integer   c_idx(c_idx_siz),real_idx(r_idx_siz),int_idx(i_idx_siz)
      character paramfl(paramsz)

      call readparm2(paramfl,c_sz,c_idx,c_idx_siz,real_sz,real_idx,
     &     r_idx_siz,int_sz,int_idx,i_idx_siz,cond_option,rseed)

      end

c :: ----------------------------------------------------------
c :: Conditioned the data using internal values. One big fab.
c :: ----------------------------------------------------------
      subroutine FORT_INTERNAL_DATA(scratch_c,c_sz,c_idx,c_idx_siz)
      implicit none
#include "sgsim2.inc"
      integer c_idx_siz
      integer c_sz,c_idx(c_idx_siz)
      REAL_T  scratch_c(c_sz)

      if (do_cond .eq. 1) then
         call read_cond_data(scratch_c,c_sz,c_idx)
      end if
      end

c :: ----------------------------------------------------------
c :: Conditioned the data using secondary values. One big fab.
c :: ----------------------------------------------------------
      subroutine FORT_SEC_DATA(dat,DIMS(dat),
     &                         scratch_r,real_sz,real_idx,r_idx_siz,
     &                         ratios,condg,sig,dx)
      implicit none
#include "sgsim2.inc"

      integer DIMDEC(dat)
      integer r_idx_siz
      integer real_sz, real_idx(r_idx_siz)
      integer ratios(SDIM)
      REAL_T  scratch_r(real_sz)
      REAL_T  dat(DIMV(dat)),dx(SDIM),condg,sig

      if (ktype .ge. 2) then
         if (do_cond .eq. 1) then
            call read_sec()
         end if
      end if

      end


c :: ----------------------------------------------------------
c :: Conditioned the data using secondary values. One big fab.
c :: ----------------------------------------------------------
      subroutine FORT_AVE_DATA(dat,DIMS(dat),sdat,DIMS(sdat),
     &     ratios,condg,sig,dx)
      implicit none
      integer DIMDEC(dat),DIMDEC(sdat),ratios(SDIM)
      REAL_T  dat(DIMV(dat)),sdat(DIMV(sdat)),dx(SDIM),condg,sig

      integer i,j,ioff,joff,ic,jc,idx,nnx,nny
      integer l1,l2,h1,h2
      real*8  lvmtmp
      
      l1 = 0
      l2 = 0
      h1 = (ARG_H1(dat)-2)*ratios(1)-1
      h2 = (ARG_H2(dat)-2)*ratios(2)-1
      
         nnx = h1 - l1 + 1
         nny = h2 - l2 + 1
         do jc = ARG_L2(dat),ARG_H2(dat)
            do ic = ARG_L1(dat),ARG_H1(dat)
               lvmtmp = -log(dat(ic,jc)/condg)/sig
               do joff = 0,ratios(2)-1
               j = (jc)*ratios(2) + joff
                  do ioff = 0,ratios(1)-1
                     i = (ic)*ratios(1) + ioff
                     if (i .ge. l1 .and. j .ge. l2 .and.
     &                   i .le. h1 .and. j .le. h2) then
                        idx = (j-l2)*nnx + (i-l1+1)
                        lvmtmp = lvmtmp + sdat(i,j)/64.d0
                     end if
                  end do
               end do
               do joff = 0,ratios(2)-1
               j = (jc)*ratios(2) + joff
                  do ioff = 0,ratios(1)-1
                     i = (ic)*ratios(1) + ioff
                     if (i .ge. l1 .and. j .ge. l2 .and.
     &                   i .le. h1 .and. j .le. h2) then
                        sdat(i,j) = sdat(i,j) - lvmtmp
                        dat(ic,jc) = condg*exp(sig*lvmtmp)
                      end if
                  end do
               end do
            end do 
         end do
         

      end

c :: ----------------------------------------------------------
c :: Sequential full SGSIM code
c :: ----------------------------------------------------------
      subroutine FORT_SGSIM(dat,DIMS(dat),dx)

      implicit none
#include "sgsim2.inc"

      integer DIMDEC(dat)
      REAL_T  dat(DIMV(dat)),dx(SDIM)
      integer nx, ny, nz
      REAL_T  xsiz,ysiz,zsiz,xmn,ymn,zmn

      nx = ARG_H1(dat)-ARG_L1(dat)+1
      ny = ARG_H2(dat)-ARG_L2(dat)+1
      nz = 1
      nxy  = nx*ny
      nxyz = nx*ny*nz

      xsiz = dx(1)
      ysiz = dx(2)
      zsiz = dx(1)

      xmn = (ARG_L1(dat)+0.5d0)*xsiz
      ymn = (ARG_L2(dat)+0.5d0)*ysiz
      zmn = 0.5d0*zsiz

      call sgsim(dat)

      end

c :: ----------------------------------------------------------
c :: Read GSLIB parameters
c :: ----------------------------------------------------------
      subroutine FORT_SGSIM_SETUP(dat,DIMS(dat),
     &     order,nvalid,
     &     scratch_c,c_sz,c_idx,c_idx_siz,
     &     scratch_r,real_sz,real_idx,r_idx_siz,
     &     scratch_i,int_sz,int_idx,i_idx_siz,
     &     lo,hi,dx,plo,domnlo,domnhi,iseed)

      implicit none
#include "sgsim2.inc"

      integer DIMDEC(dat)
      integer nvalid
      integer c_sz,real_sz,int_sz
      integer c_idx_siz,r_idx_siz,i_idx_siz
      integer c_idx(c_idx_siz),real_idx(r_idx_siz),int_idx(i_idx_siz)
      integer iseed
      integer lo(SDIM),hi(SDIM),domnlo(SDIM),domnhi(SDIM)
      integer scratch_i(int_sz)
      integer order(nvalid)
      REAL_T  dat(DIMV(dat))
      REAL_T  scratch_r(real_sz)
      REAL_T  scratch_c(c_sz)
      REAL_T  dx(SDIM),plo(SDIM)
      REAL_T  random_val(nvalid)

      integer i,j,idx,idx2,imult,index
      integer nnx,nny,ix,iy,jx,jy
      real*8  acorni

      scratch_i(int_idx(1)) = ARG_H1(dat)-ARG_L1(dat)+1
      scratch_i(int_idx(2)) = 1
      scratch_i(int_idx(3)) = ARG_H2(dat)-ARG_L2(dat)+1
      scratch_i(int_idx(14)+0) = ARG_L1(dat)
      scratch_i(int_idx(14)+1) = 0
      scratch_i(int_idx(14)+2) = ARG_L2(dat)
      scratch_i(int_idx(14)+3) = ARG_H1(dat)
      scratch_i(int_idx(14)+4) = 0
      scratch_i(int_idx(14)+5) = ARG_H2(dat)

      scratch_r(int_idx(1)) = plo(1)
      scratch_r(int_idx(2)) = 0.d0
      scratch_r(int_idx(3)) = plo(2)
      scratch_r(int_idx(4)) = dx(1)
      scratch_r(int_idx(5)) = 1.d0
      scratch_r(int_idx(6)) = dx(2)

      if (ixv(1) .eq. 0) ixv(1) = iseed

      nxyz = 1
      do i = 1,3
         nxyz = nxyz * scratch_i(int_idx(i))
      enddo

      ! Work out a random path over cells to fill (valid), use indexing 
      ! based on all cells, including grow region.  This ensures that
      ! the grow cells have an id that we can use when passing known
      ! values to the krige algorithms
      idx  = 0
      idx2 = 0
      do j = ARG_L2(dat),ARG_H2(dat)
         do i = ARG_L1(dat),ARG_H1(dat)
            idx = idx + 1
            if (i .ge. lo(1) .and. i .le. hi(1) .and.
     &          j .ge. lo(2) .and. j .le. hi(2)) then
               idx2 = idx2 + 1
               call blutilrand(random_val(idx2))
               order(idx2) = idx
            end if
         end do
      end do

      if(mults.eq.1) then
         do imult=1,nmult
            nny = max(1,(hi(2)-lo(2)+1)/(imult*4))
            nnx = max(1,(hi(1)-lo(1)+1)/(imult*4))
            jy  = 1
            jx  = 1

            do iy=1,nny
               if(nny.gt.1) jy = iy*imult*4
               do ix=1,nnx
                  if(nnx.gt.1) jx = ix*imult*4
                  index = jx + (jy-1)*(hi(1)-lo(1)+1) 
                  random_val(index) = random_val(index)-imult
               end do
            end do

         end do
      end if

      call sortem(1,idx2,random_val,1,order)
      call sgsim_setup(dat,nxyz,
     &     scratch_c,c_sz,c_idx,c_idx_siz,
     &     scratch_r,real_sz,real_idx,r_idx_siz,
     &     scratch_i,int_sz,int_idx,i_idx_siz)
      end


c :: ----------------------------------------------------------
c :: Read GSLIB parameters
c :: ----------------------------------------------------------
      subroutine FORT_SGSIM_ITER(dat,DIMS(dat),  
     &     scratch_c,c_sz,c_idx,c_idx_siz,
     &     scratch_r,real_sz,real_idx,r_idx_siz,
     &     scratch_i,int_sz,int_idx,i_idx_siz,
     &     idx)
      implicit none
      integer DIMDEC(dat)
      integer c_sz,real_sz,int_sz,idx,nxyz
      integer c_idx_siz,r_idx_siz,i_idx_siz
      integer c_idx(c_idx_siz),real_idx(r_idx_siz),int_idx(i_idx_siz)
      integer scratch_i(int_sz)
      REAL_T  dat(DIMV(dat))
      REAL_T  scratch_r(real_sz)
      REAL_T  scratch_c(c_sz)

      nxyz = (ARG_H1(dat)-ARG_L1(dat)+1)*
     &       (ARG_H2(dat)-ARG_L2(dat)+1)

      call sgsim_single_iter(idx,dat,nxyz,
     &     scratch_c,c_sz,c_idx,c_idx_siz,
     &     scratch_r,real_sz,real_idx,r_idx_siz,
     &     scratch_i,int_sz,int_idx,i_idx_siz)

      end

c :: ----------------------------------------------------------
c :: Read GSLIB parameters
c :: ----------------------------------------------------------
      subroutine FORT_SGSIM_POST(dat,DIMS(dat),
     &     scratch_c,c_sz,c_idx,c_idx_siz,
     &     scratch_r,real_sz,real_idx,r_idx_siz,
     &     scratch_i,int_sz,int_idx,i_idx_siz)
      implicit none
      integer DIMDEC(dat)
      integer c_sz,real_sz,int_sz
      integer c_idx_siz,r_idx_siz,i_idx_siz
      integer c_idx(c_idx_siz),real_idx(r_idx_siz),int_idx(i_idx_siz)
      integer scratch_i(int_sz)
      REAL_T  dat(DIMV(dat))
      REAL_T  scratch_r(real_sz)
      REAL_T  scratch_c(c_sz)
      
      integer nxyz

      nxyz = (ARG_H1(dat)-ARG_L1(dat)+1)*
     &       (ARG_H2(dat)-ARG_L2(dat)+1)

      call sgsim_post(dat,nxyz,
     &     scratch_c,c_sz,c_idx,c_idx_siz,
     &     scratch_r,real_sz,real_idx,r_idx_siz,
     &     scratch_i,int_sz,int_idx,i_idx_siz)

      end

c :: ----------------------------------------------------------
c :: Turn normal field into lognormal field
c :: ----------------------------------------------------------
      subroutine FORT_SGSIM_DEALLOC()

      call sgsim_deallocate()

      end

c :: ----------------------------------------------------------
c :: Turn normal field into lognormal field
c :: ----------------------------------------------------------
      subroutine FORT_LGNORM(dat,DIMS(dat),condg,lo,hi)
      implicit none
      integer i,j,lo(SDIM),hi(SDIM)
      integer DIMDEC(dat)
      REAL_T  condg, lc
      REAL_T  dat(DIMV(dat))

      lc = log10(condg)
      do i=lo(1),hi(1)
         do j=lo(2),hi(2)
            dat(i,j) = 10**(dat(i,j)+lc)
         end do
      end do

      end
c ::: -------------------------------------------------
c ::: This routine create non-uniform layered structure
c ::: -------------------------------------------------
      subroutine FORT_INTERP_LAYER(yl,DIMS(yl),
     &     xs,ys,ns,nl,dx,problo,probhi)

      integer DIMDEC(yl)
      integer ns,nl
      REAL_T  yl(DIMV(yl))
      REAL_T  xs(ns),ys(ns)
      REAL_T  dx(SDIM)
      REAL_T  problo(SDIM),probhi(SDIM)

      integer i,n
      REAL_T  xc
      
      if (ns .gt. 1) then         

         do i=ARG_L1(yl),ARG_H1(yl)
            xc = problo(1) + (dble(i)+half)*dx(1)
            do n = 1,ns-1
               if (xc .gt. xs(n) .and. xc .le. xs(n+1)) then
                  yl(i,ARG_L2(yl)+nl) = ys(n) + 
     &                 (ys(n+1)-ys(n))/(xs(n+1)-xs(n))*(xc-xs(n))
               end if
            end do
         end do

      else
         do i=ARG_L1(yl),ARG_H1(yl)
               yl(i,ARG_L2(yl)+nl) = ys(1)
         end do
      end if
         
      end

c ::: -------------------------------------------------
c ::: This routine create non-uniform layered structure
c ::: -------------------------------------------------
      subroutine FORT_ASSIGN_LAYER(kappa,DIMS(kappa),
     &     yl,DIMS(yl),kpval,nkpval,dkp,dx,problo,probhi)

      integer DIMDEC(kappa)
      integer DIMDEC(yl)
      integer nkpval
      REAL_T  dkp
      REAL_T  kappa(DIMV(kappa))
      REAL_T  yl(DIMV(yl))
      REAL_T  kpval(nkpval)
      REAL_T  dx(SDIM)
      REAL_T  problo(SDIM),probhi(SDIM)

      integer i,j,n
      REAL_T  yc

      do j=ARG_L2(kappa),ARG_H2(kappa)
         yc = problo(2) + (dble(j)+half)*dx(2)
         do i=ARG_L1(kappa),ARG_H1(kappa)
            if (yc .gt. yl(i,ARG_L2(kappa)+nkpval)) then
               kappa(i,j) = 1.d200
            else if (yc .le. yl(i,ARG_L2(kappa))) then
               kappa(i,j) = 1.d200
            else
               do n = 0,nkpval-1
                  if (yc .gt. yl(i,ARG_L2(kappa)+n) .and.
     &                yc .le. yl(i,ARG_L2(kappa)+n+1)) then
                     kappa(i,j) = kpval(n+1)
                  end if
               end do
            end if
         end do
      end do

      end

c ::: --------------------------------------------------------------------------
c ::: This routine create random permeability according to some layer structures
c ::: --------------------------------------------------------------------------
      subroutine FORT_PHIRAND(kappa,DIMS(kappa),kpval,nkpval,dkp,
     &     domlo,domhi,iuc)

      integer DIMDEC(kappa)
      integer nkpval,iuc
      integer domlo(SDIM),domhi(SDIM)
      REAL_T  dkp
      REAL_T  kappa(DIMV(kappa))
      REAL_T  kpval(nkpval)

      integer i,j,n
      REAL_T  phirnd
      REAL_T  y,ydiv(nkpval+1)

      call blutilinitrand(iuc)

      ydiv(nkpval+1) = dble(domhi(2))
      do n = 1,nkpval
         ydiv(n) = dble(domlo(2)) + 
     &        dble(domhi(2)-domlo(2))*dble(n-1)/dble(nkpval)
      end do 
      
      do i=ARG_L1(kappa),ARG_H1(kappa)
         do j=ARG_L2(kappa),ARG_H2(kappa)
            y = dble(j) + half
            call blutilrand(phirnd)
            do n = 1,nkpval
               if (y .ge. ydiv(n) .and. y .lt. ydiv(n+1)) then
                  kappa(i,j) = kpval(n)*(1+2*dkp*(0.5-phirnd))
               end if
            end do

            if (y .lt. ydiv(1)) then
               kappa(i,j) = kpval(1)*(1+2*dkp*(0.5-phirnd))
            else if (y .ge. ydiv(nkpval+1)) then
               kappa(i,j) = kpval(nkpval)*(1+2*dkp*(0.5-phirnd))
            end if
         end do
      end do

      end
