#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "BC_TYPES_PM.H"
#include "GSLIB_F.H" 
#include "ArrayLim.H"

#define SDIM 3


c :: ----------------------------------------------------------
c :: Read GSLIB parameters
c :: ----------------------------------------------------------
      subroutine FORT_INIT_GSLIB(paramfl,paramsz,nx_in,xsiz_in,xmn_in)

      integer   paramsz
      character paramfl(paramsz)

      integer nx_in(3)
      REAL_T  xsiz_in(3),xmn_in(3)

      call readparm(paramfl,paramsz,nx_in,xsiz_in,xmn_in)

      end

c :: ----------------------------------------------------------
c :: Read GSLIB parameters version 2
c :: ----------------------------------------------------------
      subroutine FORT_INIT_GSLIB2(paramfl,paramsz,c_sz,c_idx, 
     &     real_sz,real_idx,int_sz,int_idx,cond_option)

      integer   paramsz,cond_option
      integer   c_sz,real_sz,int_sz
      integer   c_idx(10),real_idx(20),int_idx(20)
      character paramfl(paramsz)

      call readparm2(paramfl,c_sz,c_idx,
     &     real_sz,real_idx,int_sz,int_idx,cond_option)
      
      end

      
c :: ----------------------------------------------------------
c :: Conditioned the data using internal values. One big fab.  WRONG!!!
c :: ----------------------------------------------------------
      subroutine FORT_INTERNAL_DATA(dat,DIMS(dat),
     &     scratch_c,c_sz,c_idx,condg,sig,dx,domlo,domhi)

#include "sgsim2.inc"

      integer DIMDEC(dat)
      integer c_sz,c_idx(10)
      integer domlo(SDIM),domhi(SDIM)
      REAL_T  dat(DIMV(dat))
      REAL_T  scratch_c(c_sz)
      REAL_T  dx(SDIM),condg,sig

      integer i,j,k,idx,nnx,nny,nnz,nnxyz

      if (do_cond .eq. 1) then
         call read_cond_data(scratch_c,c_sz,c_idx)

      else if (do_cond .eq. 2) then

         nnx  = domhi(1)-domlo(1)+1
         nny  = domhi(2)-domlo(2)+1
         nnz  = domhi(3)-domlo(3)+1
         nnxyz= nnx*nny*nnz
         do k = ARG_L2(dat),ARG_H2(dat)
            do j = ARG_L2(dat),ARG_H2(dat)
               do i = ARG_L1(dat),ARG_H1(dat)
                  idx =  (k-domlo(3))*nnx*nny + (j-domlo(2))*nnx + (i-domlo(1))
                  scratch_c(c_idx(1)+idx) = (dble(i)+0.5d0)*dx(1)
                  scratch_c(c_idx(2)+idx) = (dble(j)+0.5d0)*dx(2)
                  scratch_c(c_idx(3)+idx) = (dble(k)+0.5d0)*dx(3)
                  scratch_c(c_idx(4)+idx) = log(dat(i,j,k)/condg)/sig
                  scratch_c(c_idx(5)+idx) = 1
                  scratch_c(c_idx(6)+idx) = log(dat(i,j,k)/condg)/sig
                  scratch_c(c_idx(7)+idx) = 1
                  scratch_c(c_idx(8)+idx) = log(dat(i,j,k)/condg)/sig
                  scratch_c(c_idx(9)+idx) = 1
               end do
            end do 
         end do

         nd  = nnxyz
         ntr = nnxyz

      end if
      
      end

c :: ----------------------------------------------------------
c :: Conditioned the data using secondary values. One big fab.
c :: ----------------------------------------------------------
      subroutine FORT_SEC_DATA(dat,DIMS(dat),
     &                         scratch_r,real_sz,real_idx,
     &                         ratios,condg,sig,dx)

#include "sgsim2.inc"

      integer DIMDEC(dat)
      integer real_sz, real_idx(20)
      integer ratios(SDIM)
      REAL_T  scratch_r(real_sz)
      REAL_T  dat(DIMV(dat)),dx(SDIM),condg,sig

      integer i,j,k,ioff,joff,koff,ic,jc,kc,idx,nnx,nny,nnz
      integer l1,l2,l3,h1,h2,h3,idxlvm
      
      if (ktype .ge. 2) then
         if (do_cond .eq. 1) then
            call read_sec()

         else 

            idxlvm = real_idx(13)
         
            l1 = -3
            l2 = -3
            l3 = -3 
            h1 = (ARG_H1(dat)-2)*ratios(1)+2
            h2 = (ARG_H2(dat)-2)*ratios(2)+2
            h3 = (ARG_H3(dat)-2)*ratios(3)+2            
            nnx = h1 - l1 + 1
            nny = h2 - l2 + 1
            nnz = h3 - l3 + 1

            do kc = ARG_L3(dat),ARG_H3(dat)
            do koff = 0,ratios(3)-1
               k = (kc)*ratios(3) + koff            
               do jc = ARG_L2(dat),ARG_H2(dat)
               do joff = 0,ratios(2)-1
                  j = (jc)*ratios(2) + joff
                  do ic = ARG_L1(dat),ARG_H1(dat)
                  do ioff = 0,ratios(1)-1
                     i = (ic)*ratios(1) + ioff
                     if (i .ge. l1 .and. j .ge. l2 .and. k .ge. l3 .and.
     &                   i .le. h1 .and. j .le. h2 .and. k .le. h3) then
                        idx = idxlvm + (k-l3)*nnx*nny + (j-l2)*nnx + (i-l1)
                        scratch_r(idx) = log(dat(ic,jc,kc)/condg)/sig
                     end if
                  end do
                  end do
               end do 
               end do
            end do
            end do

         end if

      end if

      end


c :: ----------------------------------------------------------
c :: Conditioned the data using secondary values. One big fab.
c :: ----------------------------------------------------------
      subroutine FORT_AVE_DATA(dat,DIMS(dat),sdat,DIMS(sdat),
     &     ratios,condg,sig,dx)

      integer DIMDEC(dat),DIMDEC(sdat),ratios(SDIM)
      REAL_T  dat(DIMV(dat)),sdat(DIMV(sdat)),dx(SDIM),condg,sig

      integer i,j,k,ioff,joff,koff,ic,jc,kc,idx,nnx,nny,nnz
      integer l1,l2,l3,h1,h2,h3
      real*8  lvmtmp
      
      l1 = 0
      l2 = 0
      l3 = 0
      h1 = (ARG_H1(dat)-2)*ratios(1)-1
      h2 = (ARG_H2(dat)-2)*ratios(2)-1
      h3 = (ARG_H3(dat)-2)*ratios(3)-1

      nnx = h1 - l1 + 1
      nny = h2 - l2 + 1
      nnz = h3 - l3 + 1

      do kc = ARG_L3(dat),ARG_H3(dat)
         do jc = ARG_L2(dat),ARG_H2(dat)
            do ic = ARG_L1(dat),ARG_H1(dat)
               lvmtmp = -log(dat(ic,jc,kc)/condg)/sig
               do koff = 0,ratios(3)-1
                  k = kc*ratios(3) + koff
                  do joff = 0,ratios(2)-1
                     j = jc*ratios(2) + joff
                     do ioff = 0,ratios(1)-1
                        i = ic*ratios(1) + ioff
                        if (i .ge. l1 .and. j .ge. l2 .and. k .ge. l3 .and.
     &                      i .le. h1 .and. j .le. h2 .and. k .le. h3) then
                           idx = (k-l3)*nnx*nny + (j-l2)*nnx + (i-l1+1)
                           lvmtmp = lvmtmp + sdat(i,j,k)/64.d0
                        end if
                     end do
                  end do
               end do
               do koff = 0,ratios(3)-1
                  k = kc*ratios(3) + koff
                  do joff = 0,ratios(2)-1
                     j = (jc)*ratios(2) + joff
                     do ioff = 0,ratios(1)-1
                        i = (ic)*ratios(1) + ioff
                        if (i .ge. l1 .and. j .ge. l2 .and. k .ge. l3 .and.
     &                      i .le. h1 .and. j .le. h2 .and. k .le. h3) then
                           sdat(i,j,k) = sdat(i,j,k) - lvmtmp
                           dat(ic,jc,kc) = condg*exp(sig*lvmtmp)
                        end if
                     end do
                  end do
               end do
            end do 
         end do
      end do

      end

c :: ----------------------------------------------------------
c :: Sequential full SGSIM code
c :: ----------------------------------------------------------
      subroutine FORT_SGSIM(dat,DIMS(dat),dx)

      include 'sgsim.inc'

      integer DIMDEC(dat)
      REAL_T  dat(DIMV(dat)),dx(SDIM)

      nx = ARG_H1(dat)-ARG_L1(dat)+1
      ny = ARG_H2(dat)-ARG_L2(dat)+1
      nz = ARG_H3(dat)-ARG_L3(dat)+1
      nxy  = nx*ny
      nxyz = nx*ny*nz

      xsiz = dx(1)
      ysiz = dx(2)
      zsiz = dx(3)

      xmn = (ARG_L1(dat)+0.5d0)*xsiz
      ymn = (ARG_L2(dat)+0.5d0)*ysiz
      zmn = (ARG_L3(dat)+0.5d0)*zsiz

      call sgsim(dat)

      end

c :: ----------------------------------------------------------
c :: Read GSLIB parameters
c :: ----------------------------------------------------------
      subroutine FORT_SGSIM_SETUP(dat,DIMS(dat),
     &     order,nvalid,
     &     scratch_c,c_sz,c_idx,
     &     scratch_r,real_sz,real_idx,
     &     scratch_i,int_sz,int_idx,
     &     lo,hi,dx,iseed)

#include "sgsim2.inc"

      integer DIMDEC(dat)
      integer nvalid
      integer c_sz,real_sz,int_sz
      integer c_idx(10),real_idx(20),int_idx(20)
      integer iseed
      integer lo(SDIM),hi(SDIM)
      integer scratch_i(int_sz)
      integer order(nvalid)
      REAL_T  dat(DIMV(dat))
      REAL_T  scratch_r(real_sz)
      REAL_T  scratch_c(c_sz)
      REAL_T  dx(SDIM)
      REAL_T  random_val(nvalid)

      integer i,j,k,idx,idx2,imult,index
      integer nnx,nny,nnz,ix,iy,iz,jx,jy,jz
      real    c,d,e,f,g,h
      real*8  acorni

      scratch_i(int_idx(1)) = ARG_H1(dat)-ARG_L1(dat)+1
      scratch_i(int_idx(2)) = ARG_H2(dat)-ARG_L2(dat)+1
      scratch_i(int_idx(3)) = ARG_H3(dat)-ARG_L3(dat)+1

      scratch_r(int_idx(1)) = (ARG_L1(dat)+0.5d0)*dx(1)
      scratch_r(int_idx(2)) = (ARG_L2(dat)+0.5d0)*dx(2)
      scratch_r(int_idx(3)) = (ARG_L3(dat)+0.5d0)*dx(3)
      scratch_r(int_idx(4)) = dx(1)
      scratch_r(int_idx(5)) = dx(2)
      scratch_r(int_idx(6)) = dx(3)

      if (ixv(1) .eq. 0) ixv(1) = iseed
      nxyz = (ARG_H1(dat)-ARG_L1(dat)+1)*(ARG_H2(dat)-ARG_L2(dat)+1)*
     &       (ARG_H3(dat)-ARG_L3(dat)+1)  
c     work out a random path within search region
      idx  = 0
      idx2 = 0
      do k = ARG_L3(dat),ARG_H3(dat)
         do j = ARG_L2(dat),ARG_H2(dat)
            do i = ARG_L1(dat),ARG_H1(dat)
               idx = idx + 1
               if (i .ge. lo(1) .and. i .le. hi(1) .and. 
     &             j .ge. lo(2) .and. j .le. hi(2) .and.
     &             k .ge. lo(3) .and. k .le. hi(3)) then
                  idx2 = idx2 + 1
                  random_val(idx2) = acorni(iseed)
                  order(idx2) = idx
               end if
            end do
         end do
      end do

      if(mults.eq.1) then
         do imult=1,nmult
            nnx = max(1,(hi(1)-lo(1)+1)/(imult*4))            
            nny = max(1,(hi(2)-lo(2)+1)/(imult*4))
            nnz = max(1,(hi(3)-lo(3)+1)/(imult*4))
            jx  = 1
            jy  = 1
            jz  = 1

            do iz = 1,nnz
               if (nnz .gt. 1) jz = iz*imult*4
               do iy=1,nny
                  if(nny.gt.1) jy = iy*imult*4
                  do ix=1,nnx
                     if(nnx.gt.1) jx = ix*imult*4
                     index = jx + (jy-1)*(hi(1)-lo(1)+1) + 
     &                       (jz-1)*(hi(1)-lo(1)+1)*(hi(2)-lo(2)+1) 
                     random_val(index) = random_val(index)-imult
                  end do
               end do
            end do
         end do
      end if

      call sortem(1,idx2,random_val,1,order)
      call sgsim_setup(dat,nxyz,
     &     scratch_c,c_sz,c_idx,
     &     scratch_r,real_sz,real_idx,
     &     scratch_i,int_sz,int_idx)

      end


c :: ----------------------------------------------------------
c :: Read GSLIB parameters
c :: ----------------------------------------------------------
      subroutine FORT_SGSIM_ITER(dat,DIMS(dat),  
     &     scratch_c,c_sz,c_idx,
     &     scratch_r,real_sz,real_idx,
     &     scratch_i,int_sz,int_idx,
     &     idx)

      integer DIMDEC(dat)
      integer c_sz,real_sz,int_sz,idx,nxyz
      integer c_idx(10),real_idx(20),int_idx(20)
      integer scratch_i(int_sz)
      REAL_T  dat(DIMV(dat))
      REAL_T  scratch_r(real_sz)
      REAL_T  scratch_c(c_sz)

      nxyz = (ARG_H1(dat)-ARG_L1(dat)+1)*
     &       (ARG_H2(dat)-ARG_L2(dat)+1)*
     &       (ARG_H3(dat)-ARG_L3(dat)+1)

      call sgsim_single_iter(idx,dat,nxyz,
     &     scratch_c,c_sz,c_idx,
     &     scratch_r,real_sz,real_idx,
     &     scratch_i,int_sz,int_idx)

      end
c :: ----------------------------------------------------------
c :: Read GSLIB parameters
c :: ----------------------------------------------------------
      subroutine FORT_SGSIM_POST(dat,DIMS(dat),
     &     scratch_c,c_sz,c_idx,
     &     scratch_r,real_sz,real_idx,
     &     scratch_i,int_sz,int_idx)

      integer DIMDEC(dat)
      integer c_sz,real_sz,int_sz
      integer c_idx(10),real_idx(20),int_idx(20)
      integer scratch_i(int_sz)
      REAL_T  dat(DIMV(dat))
      REAL_T  scratch_r(real_sz)
      REAL_T  scratch_c(c_sz)
      
      integer nxyz

      nxyz = (ARG_H1(dat)-ARG_L1(dat)+1)*
     &       (ARG_H2(dat)-ARG_L2(dat)+1)*
     &       (ARG_H3(dat)-ARG_L3(dat)+1)

      call sgsim_post(dat,nxyz,
     &     scratch_c,c_sz,c_idx,
     &     scratch_r,real_sz,real_idx,
     &     scratch_i,int_sz,int_idx)

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
      subroutine FORT_LGNORM(dat,DIMS(dat),condg,sig)

      integer i,j,k
      integer DIMDEC(dat)
      REAL_T  condg,sig
      REAL_T  dat(DIMV(dat))

      do k = ARG_L3(dat),ARG_H3(dat)      
         do j=ARG_L2(dat),ARG_H2(dat)
            do i=ARG_L1(dat),ARG_H1(dat)
               dat(i,j,k)=condg*exp(sig*dat(i,j,k))
            end do
         end do
      end do

      end


c ::: -------------------------------------------------
c ::: This routine create non-uniform layered structure
c ::: -------------------------------------------------
      subroutine FORT_INTERP_LAYER(zl,DIMS(zl),
     &     xs,ys,zs,ns,nl,dx,problo,probhi)

      integer DIMDEC(zl)
      integer ns,nl
      REAL_T  zl(DIMV(zl))
      REAL_T  xs(ns),ys(ns),zs(ns)
      REAL_T  dx(SDIM)
      REAL_T  problo(SDIM),probhi(SDIM)

      integer i,j,n
      integer md,ncp,ndp,nip,nipx,nipy
      REAL_T  xc,yc,zc,x0,y0
      
      integer*4, allocatable :: iwk(:)
      REAL_T,  allocatable :: wk(:)
      REAL_T,  allocatable :: x(:),y(:),z(:)

c     shifting the reference point
      x0 = 500.d0
      y0 = 1000.d0

      if (ns .gt. 1) then
         md   = 1
         ncp  = 4
         ndp  = ns
         nipx = ARG_H1(zl)-ARG_L1(zl) + 1
         nipy = ARG_H2(zl)-ARG_L2(zl) + 1
         nip  = nipx*nipy

         allocate(iwk(max(31,27+ncp)*ndp+nip))
         allocate(wk(8*ndp))
         allocate(x(nip))
         allocate(y(nip))
         allocate(z(nip))
         
         n = 0
         do i=ARG_L1(zl),ARG_H1(zl)
            xc = problo(1) + x0 + (dble(i)+half)*dx(1)
            do j=ARG_L2(zl),ARG_H2(zl)
               yc = problo(2)+ y0 + (dble(j)+half)*dx(2)
               n  = n+1
               x(n) = xc
               y(n) = yc
            end do
         end do

         call idbvip(md,ncp,ndp,xs,ys,zs,nip,x,y,z,iwk,wk)
           
         n = 0
         do i=ARG_L1(zl),ARG_H1(zl)
            do j=ARG_L2(zl),ARG_H2(zl)
               n = n+1
               zl(i,j,ARG_L3(zl)+nl) = z(n)
               if (nl .eq. 1 .and. z(n) .lt. one) then
                  print *, i,j,x(n), y(n), z(n)
               end if
            end do
         end do
                  
         deallocate(iwk)
         deallocate(wk)
         deallocate(x)
         deallocate(y)
         deallocate(z)

      else
         do j=ARG_L2(zl),ARG_H2(zl)
            do i=ARG_L1(zl),ARG_H1(zl)
               zl(i,j,ARG_L3(zl)+nl) = zs(1)
            end do
         end do
      end if
         
      end

c ::: -------------------------------------------------
c ::: This routine create non-uniform layered structure
c ::: -------------------------------------------------
      subroutine FORT_ASSIGN_LAYER(kappa,DIMS(kappa),
     &     zl,DIMS(zl),kpval,nkpval,dkp,dx,problo,probhi)

      integer DIMDEC(kappa)
      integer DIMDEC(zl)
      integer nkpval
      REAL_T  dkp
      REAL_T  kappa(DIMV(kappa))
      REAL_T  zl(DIMV(zl))
      REAL_T  kpval(nkpval)
      REAL_T  dx(SDIM)
      REAL_T  problo(SDIM),probhi(SDIM)

      integer i,j,k,n
      REAL_T  zc, z0

c     shifting the reference point
      z0 = 0.0d0

      do k = ARG_L3(kappa),ARG_H3(kappa)
         zc = problo(3) + z0 + (dble(k)+half)*dx(3)
         do i=ARG_L1(kappa),ARG_H1(kappa)
            do j=ARG_L2(kappa),ARG_H2(kappa)
               if (zc .gt. zl(i,j,ARG_L3(kappa)+nkpval)) then
                  kappa(i,j,k) = 1.d200
               else if (zc .le. zl(i,j,ARG_L3(kappa))) then
                  kappa(i,j,k) = 1.d200
               else
                  do n = 0,nkpval-1
                     if (zc .gt. zl(i,j,ARG_L3(kappa)+n) .and.
     &                   zc .le. zl(i,j,ARG_L3(kappa)+n+1)) then
                        kappa(i,j,k) = kpval(n+1)
                     end if
                  end do
               end if
            end do
         end do
      end do

      end
