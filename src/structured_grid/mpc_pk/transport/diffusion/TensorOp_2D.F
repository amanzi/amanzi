#include <REAL.H>

#include "TensorOp_F.H"
#include <LO_BCTYPES.H>
#include "ArrayLim.H"

      subroutine FORT_APPLY(
     $     u, DIMS(u),
     $     alpha, beta,
     $     a, DIMS(a),
     $     muX, DIMS(muX),
     $     mu1X, DIMS(mu1X),
     $     muY, DIMS(muY),
     $     mu1Y, DIMS(mu1Y),
     $     out, DIMS(out),
     $     maskn,DIMS(maskn),
     $     maske,DIMS(maske),
     $     maskw,DIMS(maskw),
     $     masks,DIMS(masks),
     $     trandern,DIMS(trandern),
     $     trandere,DIMS(trandere),
     $     tranderw,DIMS(tranderw),
     $     tranders,DIMS(tranders),
     $     lo,hi,h
     $     )
      implicit none
      REAL_T alpha, beta
      integer DIMDEC(u)
      integer DIMDEC(a)
      integer DIMDEC(muX)
      integer DIMDEC(mu1X)
      integer DIMDEC(muY)
      integer DIMDEC(mu1Y)
      integer DIMDEC(out)
      integer DIMDEC(maskn)
      integer DIMDEC(maske)
      integer DIMDEC(maskw)
      integer DIMDEC(masks)
      integer DIMDEC(trandern)
      integer DIMDEC(trandere)
      integer DIMDEC(tranderw)
      integer DIMDEC(tranders)
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      REAL_T h(BL_SPACEDIM)
      REAL_T u(DIMV(u))
      REAL_T a(DIMV(a))
      REAL_T muX(DIMV(muX))
      REAL_T mu1X(DIMV(mu1X))
      REAL_T muY(DIMV(muY))
      REAL_T mu1Y(DIMV(mu1Y))
      REAL_T out(DIMV(out))
      integer  maskn(DIMV(maskn))
      integer  maske(DIMV(maske))
      integer  maskw(DIMV(maskw))
      integer  masks(DIMV(masks))

      REAL_T trandern(DIMV(trandern))
      REAL_T trandere(DIMV(trandere))
      REAL_T tranderw(DIMV(tranderw))
      REAL_T tranders(DIMV(tranders))

      integer i,j
      REAL_T dudye,dudyw,dudyn,dudys
      REAL_T dudxe,dudxw,dudxn,dudxs
      REAL_T xfluxw,xfluxe,yfluxs,yfluxn
      REAL_T hx, hxx, hx4
      REAL_T hy, hyy, hy4

      hx = 1.d0/h(1)
      hxx = hx*hx
      hx4 = 0.25*hx
      hy = 1.d0/h(2)
      hyy = hy*hy
      hy4 = 0.25*hy

c ::: cases:
c     ::: 1) interior.  All normal derivative expressions
c     ::: 2) edge but not on corner.  check mask value to see if use outside
c     :::    value or use interpolations
c     ::: 3) corners.  Must check mask values to see if use outside values
c     :::    values or use interpolations ON TWO SIDES.  But must also
c     :::    check to see if need to do one-sided derivs for tangential
c     :::    derivatives.

c ::: case 1
      do j=lo(2)+1,hi(2)-1
         do i=lo(1)+1,hi(1)-1
            dudxw = (U(i  ,j  ) - U(i-1,j  ))*hx
            dudxe = (U(i+1,j  ) - U(i  ,j  ))*hx
            dudys = (U(i  ,j  ) - U(i  ,j-1))*hy
            dudyn = (U(i  ,j+1) - U(i  ,j  ))*hy

            dudyw = (U(i-1,j+1) + U(i  ,j+1) - U(i-1,j-1) - U(i  ,j-1))*hy4
            dudye = (U(i+1,j+1) + U(i  ,j+1) - U(i+1,j-1) - U(i  ,j-1))*hy4
            dudxs = (U(i+1,j  ) + U(i+1,j-1) - U(i-1,j  ) - U(i-1,j-1))*hx4
            dudxn = (U(i+1,j  ) + U(i+1,j+1) - U(i-1,j  ) - U(i-1,j+1))*hx4

            xfluxw = - dudxw*muX(i  ,j) - dudyw*mu1X(i  ,j)
            xfluxe = - dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
            yfluxs = - dudys*muY(i  ,j) - dudxs*mu1Y(i,j  )
            yfluxn = - dudyn*muY(i,j+1) - dudxn*mu1Y(i,j+1)

            out(i,j) = alpha*a(i,j)*U(i,j) + beta*(hx*(xfluxe-xfluxw)+hy*(yfluxn-yfluxs))
         enddo
      enddo

c     ::: ::: case 2a: north side
      j = hi(2)
      do i=lo(1)+1,hi(1)-1
         dudxw = (U(i  ,j  ) - U(i-1,j  ))*hx
         dudxe = (U(i+1,j  ) - U(i  ,j  ))*hx
         dudys = (U(i  ,j  ) - U(i  ,j-1))*hy
         dudyn = (U(i  ,j+1) - U(i  ,j  ))*hy
         
         dudxs = (U(i+1,j  ) + U(i+1,j-1) - U(i-1,j  ) - U(i-1,j-1))*hx4
         dudxn = trandern(i,j+1)
         
         if( maskn(i-1,j+1).gt.0 .or. maskn(i,j+1).gt.0 ) then
            dudyw = (U(i-1,j-2) - 4*U(i-1,j-1) + 3*U(i-1,j  )
     &           +   U(i  ,j-2) - 4*U(i  ,j-1) + 3*U(i  ,j  ))*hy4
         else
            dudyw = (U(i-1,j+1) + U(i  ,j+1) - U(i-1,j-1) - U(i  ,j-1))*hy4
         endif
         if( maskn(i+1,j+1).gt.0 .or. maskn(i,j+1).gt.0 ) then
            dudye = (U(i+1,j-2) - 4*U(i+1,j-1) + 3*U(i+1,j  )
     &           +   U(i  ,j-2) - 4*U(i  ,j-1) + 3*U(i  ,j  ))*hy4
         else
            dudye = (U(i+1,j+1) + U(i  ,j+1) - U(i+1,j-1) - U(i  ,j-1))*hy4
         endif            

         xfluxw = - dudxw*muX(i  ,j) - dudyw*mu1X(i  ,j)
         xfluxe = - dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
         yfluxs = - dudys*muY(i  ,j) - dudxs*mu1Y(i,j  )
         yfluxn = - dudyn*muY(i,j+1) - dudxn*mu1Y(i,j+1)
         
         out(i,j) = alpha*a(i,j)*U(i,j) + beta*(hx*(xfluxe-xfluxw)+hy*(yfluxn-yfluxs))
      enddo
      
c     ::: case 2b south side
      j=lo(2)
      do i=lo(1)+1,hi(1)-1
         dudxw = (U(i  ,j  ) - U(i-1,j  ))*hx
         dudxe = (U(i+1,j  ) - U(i  ,j  ))*hx
         dudys = (U(i  ,j  ) - U(i  ,j-1))*hy
         dudyn = (U(i  ,j+1) - U(i  ,j  ))*hy

         dudxs = tranders(i,j-1)
         dudxn = (U(i+1,j  ) + U(i+1,j+1) - U(i-1,j  ) - U(i-1,j+1))*hx4

         if( masks(i-1,j-1).gt.0 .or. masks(i,j-1).gt.0 ) then
            dudyw = (-U(i-1,j+2) + 4*U(i-1,j+1) - 3*U(i-1,j  )
     &           -    U(i  ,j+2) + 4*U(i  ,j+1) - 3*U(i  ,j  ))*hy4
         else
            dudyw = (U(i-1,j+1) + U(i  ,j+1) - U(i-1,j-1) - U(i  ,j-1))*hy4
         endif
         if( masks(i+1,j-1).gt.0 .or. masks(i,j-1).gt.0 ) then
            dudye = (-U(i+1,j+2) + 4*U(i+1,j+1) - 3*U(i+1,j  )
     &           -    U(i  ,j+2) + 4*U(i  ,j+1) - 3*U(i  ,j  ))*hy4
         else
            dudye = (U(i+1,j+1) + U(i  ,j+1) - U(i+1,j-1) - U(i  ,j-1))*hy4
         endif            

         xfluxw = - dudxw*muX(i  ,j) - dudyw*mu1X(i  ,j)
         xfluxe = - dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
         yfluxs = - dudys*muY(i  ,j) - dudxs*mu1Y(i,j  )
         yfluxn = - dudyn*muY(i,j+1) - dudxn*mu1Y(i,j+1)

         out(i,j) = alpha*a(i,j)*U(i,j) + beta*(hx*(xfluxe-xfluxw)+hy*(yfluxn-yfluxs))
      enddo

c     ::: case 2c: west side
      i=lo(1)
      do j=lo(2)+1,hi(2)-1

         dudxw = (U(i  ,j  ) - U(i-1,j  ))*hx
         dudxe = (U(i+1,j  ) - U(i  ,j  ))*hx
         dudys = (U(i  ,j  ) - U(i  ,j-1))*hy
         dudyn = (U(i  ,j+1) - U(i  ,j  ))*hy

         dudyw = tranderw(i-1,j)
         dudye = (U(i  ,j+1) + U(i+1,j+1) - U(i  ,j-1) - U(i+1,j-1))*hy4

         if( maskw(i-1,j-1).gt.0 .or. maskw(i-1,j).gt.0 ) then
            dudxs = (-U(i+2,j-1) + 4*U(i+1,j-1) - 3*U(i  ,j-1)
     &           -    U(i+2,j  ) + 4*U(i+1,j  ) - 3*U(i  ,j  ))*hx4
         else
            dudxs = (U(i+1,j-1) + U(i+1,j  ) - U(i-1,j-1) - U(i-1,j  ))*hx4
         endif
         if( maskw(i-1,j  ).gt.0 .or. maskw(i-1,j+1).gt.0 ) then
            dudxn = (-U(i+2,j+1) + 4*U(i+1,j+1) - 3*U(i  ,j+1)
     &           -    U(i+2,j  ) + 4*U(i+1,j  ) - 3*U(i  ,j  ))*hx4
         else
            dudxn = (U(i+1,j+1) + U(i+1,j  ) - U(i-1,j+1) - U(i-1,j  ))*hx4
         endif            

         xfluxw = - dudxw*muX(i  ,j) - dudyw*mu1X(i  ,j)
         xfluxe = - dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
         yfluxs = - dudys*muY(i  ,j) - dudxs*mu1Y(i,j  )
         yfluxn = - dudyn*muY(i,j+1) - dudxn*mu1Y(i,j+1)

         out(i,j) = alpha*a(i,j)*U(i,j) + beta*(hx*(xfluxe-xfluxw)+hy*(yfluxn-yfluxs))
      enddo

c     ::: case 2d: east side
      i=hi(1)
      do j=lo(2)+1,hi(2)-1
         dudxw = (U(i  ,j  ) - U(i-1,j  ))*hx
         dudxe = (U(i+1,j  ) - U(i  ,j  ))*hx
         dudys = (U(i  ,j  ) - U(i  ,j-1))*hy
         dudyn = (U(i  ,j+1) - U(i  ,j  ))*hy

         dudyw = (U(i  ,j+1) + U(i-1,j+1) - U(i  ,j-1) - U(i-1,j-1))*hy4
         dudye = trandere(i+1,j)

         if( maske(i+1,j-1).gt.0 .or. maske(i+1,j).gt.0 ) then
            dudxs = (U(i-2,j-1) - 4*U(i-1,j-1) + 3*U(i  ,j-1)
     &           +   U(i-2,j  ) - 4*U(i-1,j  ) + 3*U(i  ,j  ))*hx4
         else
            dudxs = (U(i+1,j-1) + U(i+1,j  ) - U(i-1,j-1) - U(i-1,j  ))*hx4
         endif
         if( maske(i+1,j  ).gt.0 .or. maske(i+1,j+1).gt.0 ) then
            dudxn = (U(i-2,j+1) - 4*U(i-1,j+1) + 3*U(i  ,j+1)
     &           +   U(i-2,j  ) - 4*U(i-1,j  ) + 3*U(i  ,j  ))*hx4
         else
            dudxn = (U(i+1,j+1) + U(i+1,j  ) - U(i-1,j+1) - U(i-1,j  ))*hx4
         endif            

         xfluxw = - dudxw*muX(i  ,j) - dudyw*mu1X(i  ,j)
         xfluxe = - dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
         yfluxs = - dudys*muY(i  ,j) - dudxs*mu1Y(i,j  )
         yfluxn = - dudyn*muY(i,j+1) - dudxn*mu1Y(i,j+1)

         out(i,j) = alpha*a(i,j)*U(i,j) + beta*(hx*(xfluxe-xfluxw)+hy*(yfluxn-yfluxs))
      enddo

c     ::: case 3a: north-east corner
      j=hi(2)
      i=hi(1)
      dudxw = (U(i  ,j  ) - U(i-1,j  ))*hx
      dudxe = (U(i+1,j  ) - U(i  ,j  ))*hx
      dudys = (U(i  ,j  ) - U(i  ,j-1))*hy
      dudyn = (U(i  ,j+1) - U(i  ,j  ))*hy
      
      dudxn = trandern(i,j+1)
      dudye = trandere(i+1,j)

      if( maske(i+1,j-1).gt.0 .or. maske(i+1,j).gt.0 ) then
         dudxs = (U(i-2,j  ) - 4*U(i-1,j  ) + 3*U(i  ,j  )
     &        +   U(i-2,j-1) - 4*U(i-1,j-1) + 3*U(i  ,j-1))*hx4
      else
         dudxs = (U(i+1,j  ) + U(i+1,j-1) - U(i-1,j  ) - U(i-1,j-1))*hx4
      endif
      if( maskn(i-1,j+1).gt.0 .or. maskn(i,j+1).gt.0 ) then
         dudyw = (U(i-1,j-2) - 4*U(i-1,j-1) + 3*U(i-1,j  )
     &        +   U(i  ,j-2) - 4*U(i  ,j-1) + 3*U(i  ,j  ))*hy4
      else
         dudyw = (U(i-1,j+1) + U(i  ,j+1) - U(i-1,j-1) - U(i  ,j-1))*hy4
      endif
      
      xfluxw = - dudxw*muX(i  ,j) - dudyw*mu1X(i  ,j)
      xfluxe = - dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
      yfluxs = - dudys*muY(i  ,j) - dudxs*mu1Y(i,j  )
      yfluxn = - dudyn*muY(i,j+1) - dudxn*mu1Y(i,j+1)
      
      out(i,j) = alpha*a(i,j)*U(i,j) + beta*(hx*(xfluxe-xfluxw)+hy*(yfluxn-yfluxs))

c     ::: case 3b: south-east corner, plus i, minus j, plus i, plus j
      j=lo(2)
      i=hi(1)
      dudxw = (U(i  ,j  ) - U(i-1,j  ))*hx
      dudxe = (U(i+1,j  ) - U(i  ,j  ))*hx
      dudys = (U(i  ,j  ) - U(i  ,j-1))*hy
      dudyn = (U(i  ,j+1) - U(i  ,j  ))*hy
      
      dudxs = tranders(i,j-1)
      dudye = trandere(i+1,j)

      if( maske(i+1,j+1).gt.0 .or. maske(i+1,j).gt.0 ) then
         dudxn = (U(i-2,j  ) - 4*U(i-1,j  ) + 3*U(i  ,j  )
     &        +   U(i-2,j+1) - 4*U(i-1,j+1) + 3*U(i  ,j+1))*hx4
      else
         dudxn = (U(i+1,j  ) + U(i+1,j+1) - U(i-1,j  ) - U(i-1,j+1))*hx4
      endif
      if( masks(i-1,j-1).gt.0 .or. masks(i,j-1).gt.0 ) then
         dudyw = (-U(i-1,j+2) + 4*U(i-1,j+1) - 3*U(i-1,j  )
     &        -    U(i  ,j+2) + 4*U(i  ,j+1) - 3*U(i  ,j  ))*hy4
      else
         dudyw = (U(i-1,j+1) + U(i  ,j+1) - U(i-1,j-1) - U(i  ,j-1))*hy4
      endif
      
      xfluxw = - dudxw*muX(i  ,j) - dudyw*mu1X(i  ,j)
      xfluxe = - dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
      yfluxs = - dudys*muY(i  ,j) - dudxs*mu1Y(i,j  )
      yfluxn = - dudyn*muY(i,j+1) - dudxn*mu1Y(i,j+1)
      
      out(i,j) = alpha*a(i,j)*U(i,j) + beta*(hx*(xfluxe-xfluxw)+hy*(yfluxn-yfluxs))

c     ::: case 3c: south-west corner, minus i, minus j
      j=lo(2)
      i=lo(1)
      dudxw = (U(i  ,j  ) - U(i-1,j  ))*hx
      dudxe = (U(i+1,j  ) - U(i  ,j  ))*hx
      dudys = (U(i  ,j  ) - U(i  ,j-1))*hy
      dudyn = (U(i  ,j+1) - U(i  ,j  ))*hy
      
      dudxs = tranders(i,j-1)
      dudyw = tranderw(i-1,j)

      if( maskw(i-1,j+1).gt.0 .or. maskw(i-1,j).gt.0 ) then
         dudxn = (-U(i+2,j  ) + 4*U(i+1,j  ) - 3*U(i  ,j  )
     &        -    U(i+2,j+1) + 4*U(i+1,j+1) - 3*U(i  ,j+1))*hx4
      else
         dudxn = (U(i+1,j  ) + U(i+1,j+1) - U(i-1,j  ) - U(i-1,j+1))*hx4
      endif
      if( masks(i  ,j-1).gt.0 .or. masks(i+1,j-1).gt.0 ) then
         dudye = (-U(i-1,j+2) + 4*U(i-1,j+1) - 3*U(i-1,j  )
     &        -    U(i  ,j+2) + 4*U(i  ,j+1) - 3*U(i  ,j  ))*hy4
      else
         dudye = (U(i-1,j+1) + U(i  ,j+1) - U(i-1,j-1) - U(i  ,j-1))*hy4
      endif
      
      xfluxw = - dudxw*muX(i  ,j) - dudyw*mu1X(i  ,j)
      xfluxe = - dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
      yfluxs = - dudys*muY(i  ,j) - dudxs*mu1Y(i,j  )
      yfluxn = - dudyn*muY(i,j+1) - dudxn*mu1Y(i,j+1)
      
      out(i,j) = alpha*a(i,j)*U(i,j) + beta*(hx*(xfluxe-xfluxw)+hy*(yfluxn-yfluxs))


c     ::: case 3d: north-west corner, minus i, plus j
      j=hi(2)
      i=lo(1)
      dudxw = (U(i  ,j  ) - U(i-1,j  ))*hx
      dudxe = (U(i+1,j  ) - U(i  ,j  ))*hx
      dudys = (U(i  ,j  ) - U(i  ,j-1))*hy
      dudyn = (U(i  ,j+1) - U(i  ,j  ))*hy
      
      dudxn = trandern(i,j+1)
      dudyw = tranderw(i-1,j)

      if( maskw(i-1,j-1).gt.0 .or. maskw(i-1,j).gt.0 ) then
         dudxs = (-U(i+2,j  ) + 4*U(i+1,j  ) - 3*U(i  ,j  )
     &        -    U(i+2,j-1) + 4*U(i+1,j-1) - 3*U(i  ,j-1))*hx4
      else
         dudxs = (U(i+1,j  ) + U(i+1,j-1) - U(i-1,j  ) - U(i-1,j-1))*hx4
      endif
      if( maskn(i  ,j+1).gt.0 .or. maskn(i+1,j+1).gt.0 ) then
         dudye = ( U(i+1,j-2) - 4*U(i+1,j-1) + 3*U(i+1,j  )
     &        +    U(i  ,j-2) - 4*U(i  ,j-1) + 3*U(i  ,j  ))*hy4
      else
         dudye = (U(i+1,j+1) + U(i  ,j+1) - U(i+1,j-1) - U(i  ,j-1))*hy4
      endif
      
      xfluxw = - dudxw*muX(i  ,j) - dudyw*mu1X(i  ,j)
      xfluxe = - dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
      yfluxs = - dudys*muY(i  ,j) - dudxs*mu1Y(i,j  )
      yfluxn = - dudyn*muY(i,j+1) - dudxn*mu1Y(i,j+1)
      
      out(i,j) = alpha*a(i,j)*U(i,j) + beta*(hx*(xfluxe-xfluxw)+hy*(yfluxn-yfluxs))
      end


      subroutine FORT_FLUX(
     $     u, DIMS(u),
     $     muX, DIMS(muX),
     $     mu1X, DIMS(mu1X),
     $     muY, DIMS(muY),
     $     mu1Y, DIMS(mu1Y),
     $     xflux,DIMS(xflux),
     $     yflux,DIMS(yflux),
     $     maskn,DIMS(maskn),
     $     maske,DIMS(maske),
     $     maskw,DIMS(maskw),
     $     masks,DIMS(masks),
     $     trandern,DIMS(trandern),
     $     trandere,DIMS(trandere),
     $     tranderw,DIMS(tranderw),
     $     tranders,DIMS(tranders),
     $     lo,hi,h
     $     )
      implicit none
      integer DIMDEC(u)
      integer DIMDEC(muX)
      integer DIMDEC(mu1X)
      integer DIMDEC(muY)
      integer DIMDEC(mu1Y)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(maskn)
      integer DIMDEC(maske)
      integer DIMDEC(maskw)
      integer DIMDEC(masks)
      integer DIMDEC(trandern)
      integer DIMDEC(trandere)
      integer DIMDEC(tranderw)
      integer DIMDEC(tranders)
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      REAL_T h(BL_SPACEDIM)
      REAL_T u(DIMV(u))
      REAL_T muX(DIMV(muX))
      REAL_T mu1X(DIMV(mu1X))
      REAL_T muY(DIMV(muY))
      REAL_T mu1Y(DIMV(mu1Y))
      REAL_T xflux(DIMV(xflux))
      REAL_T yflux(DIMV(yflux))
      integer  maskn(DIMV(maskn))
      integer  maske(DIMV(maske))
      integer  maskw(DIMV(maskw))
      integer  masks(DIMV(masks))

      REAL_T trandern(DIMV(trandern))
      REAL_T trandere(DIMV(trandere))
      REAL_T tranderw(DIMV(tranderw))
      REAL_T tranders(DIMV(tranders))

      integer i,j
      REAL_T hx,hy,hx4,hy4
      REAL_T dudxe,dudxw,dudxn,dudxs
      REAL_T dudye,dudyw,dudyn,dudys

      hx  = 1.d0/h(1)
      hx4 = hx*0.25d0
      hy  = 1.d0/h(2)
      hy4 = hy*0.25d0

c ::: cases:
c     ::: 1) interior.  All normal derivative expressions
c     ::: 2) edge but not on corner.  check mask value to see if use outside
c     :::    value or use interpolations
c     ::: 3) corners.  Must check mask values to see if use outside values
c     :::    values or use interpolations ON TWO SIDES.  But must also
c     :::    check to see if need to do one-sided derivs for tangential
c     :::    derivatives.

c ::: case 1
      do j=lo(2)+1,hi(2)-1
         do i=lo(1)+1,hi(1)-1
            dudxw = (U(i  ,j  ) - U(i-1,j  ))*hx
            dudys = (U(i  ,j  ) - U(i  ,j-1))*hy
            dudyw = (U(i-1,j+1) + U(i  ,j+1) - U(i-1,j-1) - U(i  ,j-1))*hy4
            dudxs = (U(i+1,j  ) + U(i+1,j-1) - U(i-1,j  ) - U(i-1,j-1))*hx4
            xflux(i,j) = -dudxw*muX(i,j) - dudyw*mu1X(i,j)
            yflux(i,j) = -dudys*muY(i,j) - dudxs*mu1Y(i,j)
         enddo
      enddo

c ::: case 2a: north side, do fluxes at (i-1/2,j), (i,j-1/2) and (i,j+1/2)
      j = hi(2)
      do i=lo(1)+1,hi(1)-1
         dudxw = (U(i  ,j  ) - U(i-1,j  ))*hx
         dudys = (U(i  ,j  ) - U(i  ,j-1))*hy
         dudyn = (U(i  ,j+1) - U(i  ,j  ))*hy
         dudxn = trandern(i,j+1)
         if( maskn(i-1,j+1).gt.0 .or. maskn(i,j+1).gt.0 ) then
            dudyw = (U(i-1,j-2) - 4*U(i-1,j-1) + 3*U(i-1,j  )
     &           +   U(i  ,j-2) - 4*U(i  ,j-1) + 3*U(i  ,j  ))*hy4
         else
            dudyw = (U(i-1,j+1) + U(i  ,j+1) - U(i-1,j-1) - U(i  ,j-1))*hy4
         endif
         xflux(i,j) = -dudxw*muX(i,j) - dudyw*mu1X(i,j)
         yflux(i,j) = -dudxs*mu1Y(i,j) - dudys*muY(i,j)
         yflux(i,j+1) = -dudxn*mu1Y(i,j+1) - dudyn*muY(i,j+1)
      enddo


c ::: case 2b: south side, do fluxes at (i-1/2,j), (i,j-1/2)
      j=lo(2)
      do i=lo(1)+1,hi(1)-1
         dudxw = (U(i  ,j  ) - U(i-1,j  ))*hx
         dudys = (U(i  ,j  ) - U(i  ,j-1))*hy
         dudxs = tranders(i,j-1)
         if( masks(i-1,j-1).gt.0 .or. masks(i,j-1).gt.0 ) then
            dudyw = (-U(i-1,j+2) + 4*U(i-1,j+1) - 3*U(i-1,j  )
     &           -    U(i  ,j+2) + 4*U(i  ,j+1) - 3*U(i  ,j  ))*hy4
         else
            dudyw = (U(i-1,j+1) + U(i  ,j+1) - U(i-1,j-1) - U(i  ,j-1))*hy4
         endif
         xflux(i,j) = -dudxw*muX(i,j) - dudyw*mu1X(i,j)
         yflux(i,j) = -dudxs*mu1Y(i,j) - dudys*muY(i,j)
      enddo


c ::: case 2c: west side, do fluxes at (i-1/2,j), (i,j-1/2)
      i=lo(1)
      do j=lo(2)+1,hi(2)-1
         dudxw = (U(i  ,j  ) - U(i-1,j  ))*hx
         dudys = (U(i  ,j  ) - U(i  ,j-1))*hx
         dudyw = tranderw(i-1,j)
         if( maskw(i-1,j-1).gt.0 .or. maskw(i-1,j).gt.0 ) then
            dudxs = (-U(i+2,j  ) + 4*U(i+1,j  ) - 3*U(i  ,j  )
     &           -    U(i+2,j-1) + 4*U(i+1,j-1) - 3*U(i  ,j-1))*hx4
         else
            dudxs = (U(i  ,j  ) + U(i  ,j-1) - U(i-1,j  ) - U(i-1,j-1))*hx4
         endif
         xflux(i,j) = -dudxw*muX(i,j) - dudyw*mu1X(i,j)
         yflux(i,j) = -dudxs*mu1Y(i,j) - dudys*muY(i,j)
      enddo

c ::: case 2d: east side, do fluxes at (i-1/2,j), (i,j-1/2) and (i+1/2,j)
      i=hi(1)
      do j=lo(2)+1,hi(2)-1
         dudxw = (U(i  ,j  ) - U(i-1,j  ))*hx
         dudxe = (U(i+1,j  ) - U(i  ,j  ))*hx
         dudys = (U(i  ,j  ) - U(i  ,j-1))*hy
         dudyw = (U(i-1,j+1) + U(i  ,j+1) - U(i-1,j-1) - U(i  ,j-1))*hy4
         dudye = trandere(i+1,j)
         if( maske(i+1,j-1).gt.0 .or. maske(i+1,j).gt.0 ) then
            dudxs = (U(i-2,j  ) - 4*U(i-1,j  ) + 3*U(i  ,j  )
     &           +   U(i-2,j-1) - 4*U(i-1,j-1) + 3*U(i  ,j-1))*hx4
         else
            dudxs = (U(i+1,j  ) + U(i+1,j-1) - U(i-1,j  ) - U(i-1,j-1))*hx4
         endif
         xflux(i,j) = -dudxw*muX(i,j) - dudyw*mu1X(i,j)
         xflux(i+1,j) = -dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
         yflux(i,j) = -dudxs*mu1Y(i,j) - dudys*muY(i,j)
      enddo

c ::: case 3a: north-east corner, do fluxes at (i-1/2,j), (i,j-1/2), (i+1/2,j) and (i+1/2,j+1/2)
      j=hi(2)
      i=hi(1)
      dudxe = (U(i+1,j  ) - U(i  ,j  ))*hx
      dudyn = (U(i  ,j+1) - U(i  ,j  ))*hy
      dudys = (U(i  ,j  ) - U(i  ,j-1))*hy
      dudye = trandere(i+1,j)
      dudxn = trandern(i,j+1)
      if( maske(i+1,j-1).gt.0 .or. maske(i+1,j).gt.0 ) then
         dudxs = (U(i-2,j  ) - 4*U(i-1,j  ) + 3*U(i  ,j  )
     &        +   U(i-2,j-1) - 4*U(i-1,j-1) + 3*U(i  ,j-1))*hx4
      else
         dudxs = (U(i+1,j  ) + U(i+1,j-1) - U(i-1,j  ) - U(i-1,j-1))*hx4
      endif
      if( maskn(i-1,j+1).gt.0 .or. maskn(i,j+1).gt.0 ) then
         dudyw = (U(i-1,j-2) - 4*U(i-1,j-1) + 3*U(i-1,j  )
     &        +   U(i  ,j-2) - 4*U(i  ,j-1) + 3*U(i  ,j  ))*hy4
      else
         dudyw = (U(i-1,j+1) + U(i  ,j+1) - U(i-1,j-1) - U(i  ,j-1))*hy4
      endif
      xflux(i,j) = -dudxw*muX(i,j) - dudyw*mu1X(i,j)
      xflux(i+1,j) = -dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
      yflux(i,j) = -dudxs*mu1Y(i,j) - dudys*muY(i,j)
      yflux(i,j+1) = -dudxn*mu1Y(i,j+1) - dudyn*muY(i,j+1)


c ::: case 3b: south-east corner, do fluxes at (i-1/2,j), (i,j-1/2), (i+1/2,j)
      j=lo(2)
      i=hi(1)
      dudxe = (U(i+1,j  ) - U(i  ,j  ))*hx
      dudys = (U(i  ,j  ) - U(i  ,j-1))*hy
      dudye = trandere(i+1,j)
      dudxs = tranders(i,j-1)
      if( masks(i-1,j-1).gt.0 .or. masks(i,j-1).gt.0 ) then
         dudyw = (-U(i-1,j+2) + 4*U(i-1,j+1) - 3*U(i-1,j  )
     &        -    U(i  ,j+2) + 4*U(i  ,j+1) - 3*U(i  ,j  ))*hy4
      else
         dudyw = (U(i-1,j+1) + U(i  ,j+1) - U(i-1,j-1) - U(i  ,j-1))*hy4
      endif
      xflux(i,j) = -dudxw*muX(i,j) - dudyw*mu1X(i,j)
      xflux(i+1,j) = -dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
      yflux(i,j) = -dudxs*mu1Y(i,j) - dudys*muY(i,j)

c ::: case 3c: south-west corner, , do fluxes at (i-1/2,j), (i,j-1/2), (i+1/2,j) and (i+1/2,j+1/2)
      j=lo(2)
      i=lo(1)
      dudxw = (U(i  ,j  ) - U(i-1,j  ))*hx
      dudys = (U(i  ,j  ) - U(i  ,j-1))*hy
      dudyw = tranderw(i-1,j)
      dudxs = tranders(i,j-1)
      xflux(i,j) = -dudxw*muX(i,j) - dudyw*mu1X(i,j)
      yflux(i,j) = -dudxs*mu1Y(i,j) - dudys*muY(i,j)


c ::: case 3d: north-west corner, , do fluxes at (i-1/2,j), (i,j-1/2), (i+1/2,j) and (i+1/2,j+1/2)
      j=hi(2)
      i=lo(1)
      dudxw = (U(i  ,j  ) - U(i-1,j  ))*hx
      dudyn = (U(i  ,j+1) - U(i  ,j  ))*hy
      dudys = (U(i  ,j  ) - U(i  ,j-1))*hy
      dudyw = tranderw(i-1,j)
      dudxn = trandern(i,j+1)
      if( maskw(i-1,j-1).gt.0 .or. maskw(i-1,j).gt.0 ) then
         dudxs = (-U(i+2,j  ) + 4*U(i+1,j  ) - 3*U(i  ,j  )
     &        -    U(i+2,j-1) + 4*U(i+1,j-1) - 3*U(i  ,j-1))*hx4
      else
         dudxs = (U(i+1,j  ) + U(i+1,j-1) - U(i-1,j  ) - U(i-1,j-1))*hx4
      endif
      if( maskn(i-1,j+1).gt.0 .or. maskn(i,j+1).gt.0 ) then
         dudye = (U(i+1,j-2) - 4*U(i+1,j-1) + 3*U(i+1,j  )
     &        +   U(i  ,j-2) - 4*U(i  ,j-1) + 3*U(i  ,j  ))*hy4
      else
         dudye = (U(i  ,j+1) + U(i+1,j+1) - U(i  ,j-1) - U(i+1,j-1))*hy4
      endif
      xflux(i,j) = -dudxw*muX(i,j) - dudyw*mu1X(i,j)
      yflux(i,j) = -dudxs*mu1Y(i,j) - dudys*muY(i,j)
      yflux(i,j+1) = -dudxn*mu1Y(i,j+1) - dudyn*muY(i,j+1)

      end

      subroutine FORT_GSRB (
     $     u, DIMS(u),
     $     rhs, DIMS(rhs),
     $     alpha, beta,
     $     a, DIMS(a),
     $     muX, DIMS(muX),
     $     mu1X, DIMS(mu1X),
     $     muY, DIMS(muY),
     $     mu1Y, DIMS(mu1Y),
     $     maskn,DIMS(maskn),
     $     fn, DIMS(fn),
     $     maske,DIMS(maske),
     $     fe, DIMS(fe),
     $     maskw,DIMS(maskw),
     $     fw, DIMS(fw),
     $     masks,DIMS(masks),
     $     fs, DIMS(fs),
     $     trandern,DIMS(trandern),
     $     trandere,DIMS(trandere),
     $     tranderw,DIMS(tranderw),
     $     tranders,DIMS(tranders),
     $     lo,hi,h,phaseflag
     $     )

      implicit none
      REAL_T alpha, beta
      integer DIMDEC(u)
      integer DIMDEC(rhs)
      integer DIMDEC(a)
      integer DIMDEC(muX)
      integer DIMDEC(mu1X)
      integer DIMDEC(muY)
      integer DIMDEC(mu1Y)
      integer DIMDEC(maskn)
      integer DIMDEC(maske)
      integer DIMDEC(maskw)
      integer DIMDEC(masks)
      integer DIMDEC(fn)
      integer DIMDEC(fe)
      integer DIMDEC(fw)
      integer DIMDEC(fs)
      integer DIMDEC(trandern)
      integer DIMDEC(trandere)
      integer DIMDEC(tranderw)
      integer DIMDEC(tranders)
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      REAL_T h(BL_SPACEDIM)
      integer phaseflag

      REAL_T u(DIMV(u))
      REAL_T rhs(DIMV(rhs))
      REAL_T a(DIMV(a))
      REAL_T muX(DIMV(muX))
      REAL_T mu1X(DIMV(mu1X))
      REAL_T muY(DIMV(muY))
      REAL_T mu1Y(DIMV(mu1Y))
      integer maskn(DIMV(maskn))
      integer maske(DIMV(maske))
      integer maskw(DIMV(maskw))
      integer masks(DIMV(masks))
      REAL_T fn(DIMV(fn))
      REAL_T fe(DIMV(fe))
      REAL_T fw(DIMV(fw))
      REAL_T fs(DIMV(fs))
      REAL_T trandern(DIMV(trandern))
      REAL_T trandere(DIMV(trandere))
      REAL_T tranderw(DIMV(tranderw))
      REAL_T tranders(DIMV(tranders))

      integer i,j,modx,mody
      integer istart,jstart
      REAL_T dudye,dudyw,dudyn,dudys
      REAL_T dudxe,dudxw,dudxn,dudxs
      REAL_T xfluxw,xfluxe,yfluxs,yfluxn
      REAL_T diag, oper
      REAL_T hx, hxx, hx4
      REAL_T hy, hyy, hy4
      REAL_T dw, de, ds, dn
      REAL_T fac


      fac = 1.d0
      hx = 1.d0/h(1)
      hxx = hx*hx
      hx4 = 0.25*hx
      hy = 1.d0/h(2)
      hyy = hy*hy
      hy4 = 0.25*hy

c ::: convert phase flag into modx and mody
      if(    phaseflag.eq.0 ) then
         modx = 0
         mody = 0
      elseif(phaseflag.eq.1) then
         modx = 1
         mody = 0
      elseif(phaseflag.eq.2) then
         modx = 0
         mody = 1
      elseif(phaseflag.eq.3) then
         modx = 1
         mody = 1
      else
         write(6,*)'FORT_GSRB: bad phaseflag', phaseflag
         stop
      endif

c ::: cases:
c     ::: 1) interior.  All normal derivative expressions
c     ::: 2) edge but not on corner.  check mask value to see if use outside
c     :::    value or use interpolations
c     ::: 3) corners.  Must check mask values to see if use outside values
c     :::    values or use interpolations ON TWO SIDES.  But must also
c     :::    check to see if need to do one-sided derivs for tangential
c     :::    derivatives.

c ::: case 1
      istart = lo(1)+1
      if( mod(istart,2) .ne. modx ) istart = istart+1
      jstart = lo(2)+1
      if( mod(jstart,2) .ne. mody ) jstart = jstart+1
      do j = jstart,hi(2)-1,2
         do i = istart,hi(1)-1,2
            dudxw = (U(i  ,j  ) - U(i-1,j  ))*hx
            dudxe = (U(i+1,j  ) - U(i  ,j  ))*hx
            dudys = (U(i  ,j  ) - U(i  ,j-1))*hy
            dudyn = (U(i  ,j+1) - U(i  ,j  ))*hy

            dudyw = (U(i-1,j+1) + U(i  ,j+1) - U(i-1,j-1) - U(i  ,j-1))*hy4
            dudye = (U(i+1,j+1) + U(i  ,j+1) - U(i+1,j-1) - U(i  ,j-1))*hy4
            dudxs = (U(i+1,j  ) + U(i+1,j-1) - U(i-1,j  ) - U(i-1,j-1))*hx4
            dudxn = (U(i+1,j  ) + U(i+1,j+1) - U(i-1,j  ) - U(i-1,j+1))*hx4

            xfluxw = - dudxw*muX(i  ,j) - dudyw*mu1X(i  ,j)
            xfluxe = - dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
            yfluxs = - dudys*muY(i  ,j) - dudxs*mu1Y(i,j  )
            yfluxn = - dudyn*muY(i,j+1) - dudxn*mu1Y(i,j+1)

            diag = alpha*a(i,j) + 
     $           beta*( ( muX(i,j) + muX(i+1,j) )*hxx
     $           +      ( muY(i,j) + muY(i,j+1) )*hyy )

            oper = alpha*a(i,j)*U(i,j) + beta*(hx*(xfluxe-xfluxw)+hy*(yfluxn-yfluxs))
            U(i,j) = U(i,j) + fac*(rhs(i,j) - oper)/diag
         enddo
      enddo

c ::: case 2a, north side
      j = hi(2)
      if( mod(j,2) .eq. mody ) then
         do i = istart,hi(1)-1,2
            dudxw = (U(i  ,j  ) - U(i-1,j  ))*hx
            dudxe = (U(i+1,j  ) - U(i  ,j  ))*hx
            dudys = (U(i  ,j  ) - U(i  ,j-1))*hy
            dudyn = (U(i  ,j+1) - U(i  ,j  ))*hy

            dudxs = (U(i+1,j  ) + U(i+1,j-1) - U(i-1,j  ) - U(i-1,j-1))*hx4
            dudxn = trandern(i,j+1)

            if( maskn(i-1,j+1).gt.0 .or. maskn(i,j+1).gt.0 ) then
               dudyw = (U(i-1,j-2) - 4*U(i-1,j-1) + 3*U(i-1,j  )
     &              +   U(i  ,j-2) - 4*U(i  ,j-1) + 3*U(i  ,j  ))*hy4
               dw = 3*hy4
            else
               dudyw = (U(i-1,j+1) + U(i  ,j+1) - U(i-1,j-1) - U(i  ,j-1))*hy4
               dw = 0.d0
            endif
            if( maskn(i+1,j+1).gt.0 .or. maskn(i,j+1).gt.0 ) then
               dudye = (U(i+1,j-2) - 4*U(i+1,j-1) + 3*U(i+1,j  )
     &              +   U(i  ,j-2) - 4*U(i  ,j-1) + 3*U(i  ,j  ))*hy4
               de = 3*hy4
            else
               dudye = (U(i+1,j+1) + U(i  ,j+1) - U(i+1,j-1) - U(i  ,j-1))*hy4
               de = 0.d0
            endif            
	    if( maskn(i,j+1).gt. 0) then
	       dn = - fn(i,j)
            else
               dn = 0.d0
            endif

            xfluxw = - dudxw*muX(i  ,j) - dudyw*mu1X(i  ,j)
            xfluxe = - dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
            yfluxs = - dudys*muY(i  ,j) - dudxs*mu1Y(i,j  )
            yfluxn = - dudyn*muY(i,j+1) - dudxn*mu1Y(i,j+1)

            oper = alpha*a(i,j)*U(i,j) + beta*(hx*(xfluxe-xfluxw)+hy*(yfluxn-yfluxs))
            diag = alpha*a(i,j) + 
     $           beta*( muX(i  ,j)*hxx - mu1X(i  ,j)*hx*dw
     $           +      muX(i+1,j)*hxx + mu1X(i+1,j)*hx*de
     $           +      muY(i,j  )*hyy
     $           +      muY(i,j+1)*hyy*(1.d0 + dn) )
            U(i,j) = U(i,j) + fac*(rhs(i,j) - oper)/diag
         enddo
      endif

c ::: case 2b, south side
      j=lo(2)
      if( mod(j,2) .eq. mody ) then
         do i=istart,hi(1)-1,2
            dudxw = (U(i  ,j  ) - U(i-1,j  ))*hx
            dudxe = (U(i+1,j  ) - U(i  ,j  ))*hx
            dudys = (U(i  ,j  ) - U(i  ,j-1))*hy
            dudyn = (U(i  ,j+1) - U(i  ,j  ))*hy

            dudxs = tranders(i,j-1)
            dudxn = (U(i+1,j  ) + U(i+1,j+1) - U(i-1,j  ) - U(i-1,j+1))*hx4

            if( masks(i-1,j-1).gt.0 .or. masks(i,j-1).gt.0 ) then
               dudyw = (-U(i-1,j+2) + 4*U(i-1,j+1) - 3*U(i-1,j  )
     &              -    U(i  ,j+2) + 4*U(i  ,j+1) - 3*U(i  ,j  ))*hy4
               dw = - 3*hy4
            else
               dudyw = (U(i-1,j+1) + U(i  ,j+1) - U(i-1,j-1) - U(i  ,j-1))*hy4
               dw = 0.d0
            endif
            if( masks(i+1,j-1).gt.0 .or. masks(i,j-1).gt.0 ) then
               dudye = (-U(i+1,j+2) + 4*U(i+1,j+1) - 3*U(i+1,j  )
     &              -    U(i  ,j+2) + 4*U(i  ,j+1) - 3*U(i  ,j  ))*hy4
               de = - 3*hy4
            else
               dudye = (U(i+1,j+1) + U(i  ,j+1) - U(i+1,j-1) - U(i  ,j-1))*hy4
               de = 0.d0
            endif            
	    if( masks(i,j-1).gt. 0) then
	       ds = - fs(i,j)
            else
               ds = 0.d0
            endif

            xfluxw = - dudxw*muX(i  ,j) - dudyw*mu1X(i  ,j)
            xfluxe = - dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
            yfluxs = - dudys*muY(i  ,j) - dudxs*mu1Y(i,j  )
            yfluxn = - dudyn*muY(i,j+1) - dudxn*mu1Y(i,j+1)

            oper = alpha*a(i,j)*U(i,j) + beta*(hx*(xfluxe-xfluxw)+hy*(yfluxn-yfluxs))
            diag = alpha*a(i,j) + 
     $           beta*( muX(i  ,j)*hxx - mu1X(i  ,j)*hx*dw
     $           +      muX(i+1,j)*hxx + mu1X(i+1,j)*hx*de
     $           +      muY(i,j  )*hyy*(1.d0 + ds)
     $           +      muY(i,j+1)*hyy )
            U(i,j) = U(i,j) + fac*(rhs(i,j) - oper)/diag
         enddo
      endif


c ::: case 2c: west side
      i=lo(1)
      if( mod(i,2) .eq. modx ) then
         do j=jstart,hi(2)-1,2
            dudxw = (U(i  ,j  ) - U(i-1,j  ))*hx
            dudxe = (U(i+1,j  ) - U(i  ,j  ))*hx
            dudys = (U(i  ,j  ) - U(i  ,j-1))*hy
            dudyn = (U(i  ,j+1) - U(i  ,j  ))*hy

            dudyw = tranderw(i-1,j)
            dudye = (U(i  ,j+1) + U(i+1,j+1) - U(i  ,j-1) - U(i+1,j-1))*hy4

            if( maskw(i-1,j-1).gt.0 .or. maskw(i-1,j).gt.0 ) then
               dudxs = (-U(i+2,j-1) + 4*U(i+1,j-1) - 3*U(i  ,j-1)
     &              -    U(i+2,j  ) + 4*U(i+1,j  ) - 3*U(i  ,j  ))*hx4
               ds = - 3*hx4
            else
               dudxs = (U(i+1,j-1) + U(i+1,j  ) - U(i-1,j-1) - U(i-1,j  ))*hx4
               ds = 0.d0
            endif
            if( maskw(i-1,j  ).gt.0 .or. maskw(i-1,j+1).gt.0 ) then
               dudxn = (-U(i+2,j+1) + 4*U(i+1,j+1) - 3*U(i  ,j+1)
     &              -    U(i+2,j  ) + 4*U(i+1,j  ) - 3*U(i  ,j  ))*hx4
               dn = - 3*hx4
            else
               dudxn = (U(i+1,j+1) + U(i+1,j  ) - U(i-1,j+1) - U(i-1,j  ))*hx4
               dn = 0.d0
            endif            
	    if( maskw(i-1,j  ).gt. 0) then
	       dw = - fw(i,j)
            else
               dw = 0.d0
            endif

            xfluxw = - dudxw*muX(i  ,j) - dudyw*mu1X(i  ,j)
            xfluxe = - dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
            yfluxs = - dudys*muY(i  ,j) - dudxs*mu1Y(i,j  )
            yfluxn = - dudyn*muY(i,j+1) - dudxn*mu1Y(i,j+1)

            oper = alpha*a(i,j)*U(i,j) + beta*(hx*(xfluxe-xfluxw)+hy*(yfluxn-yfluxs))
            diag = alpha*a(i,j) + 
     $           beta*( muX(i  ,j)*hxx*(1.d0 + dw)
     $           +      muX(i+1,j)*hxx
     $           +      muY(i,j  )*hyy - mu1Y(i  ,j)*hy*ds
     $           +      muY(i,j+1)*hyy + mu1Y(i,j+1)*hy*dn)
            U(i,j) = U(i,j) + fac*(rhs(i,j) - oper)/diag
         enddo
      endif

c ::: case 2d: east side
      i=hi(1)
      if( mod(i,2) .eq. modx ) then
         do j=jstart,hi(2)-1,2
            dudxw = (U(i  ,j  ) - U(i-1,j  ))*hx
            dudxe = (U(i+1,j  ) - U(i  ,j  ))*hx
            dudys = (U(i  ,j  ) - U(i  ,j-1))*hy
            dudyn = (U(i  ,j+1) - U(i  ,j  ))*hy

            dudyw = (U(i  ,j+1) + U(i-1,j+1) - U(i  ,j-1) - U(i-1,j-1))*hy4
            dudye = trandere(i+1,j)

            if( maske(i+1,j-1).gt.0 .or. maske(i+1,j).gt.0 ) then
               dudxs = (U(i-2,j-1) - 4*U(i-1,j-1) + 3*U(i  ,j-1)
     &              +   U(i-2,j  ) - 4*U(i-1,j  ) + 3*U(i  ,j  ))*hx4
               ds = 3*hx4
            else
               dudxs = (U(i+1,j-1) + U(i+1,j  ) - U(i-1,j-1) - U(i-1,j  ))*hx4
               ds = 0.d0
            endif
            if( maske(i+1,j  ).gt.0 .or. maske(i+1,j+1).gt.0 ) then
               dudxn = (U(i-2,j+1) - 4*U(i-1,j+1) + 3*U(i  ,j+1)
     &              +   U(i-2,j  ) - 4*U(i-1,j  ) + 3*U(i  ,j  ))*hx4
               dn = 3*hx4
            else
               dudxn = (U(i+1,j+1) + U(i+1,j  ) - U(i-1,j+1) - U(i-1,j  ))*hx4
               dn = 0.d0
            endif            
	    if( maske(i+1,j  ).gt. 0) then
	       de = - fe(i,j)
            else
               de = 0.d0
            endif

            xfluxw = - dudxw*muX(i  ,j) - dudyw*mu1X(i  ,j)
            xfluxe = - dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
            yfluxs = - dudys*muY(i  ,j) - dudxs*mu1Y(i,j  )
            yfluxn = - dudyn*muY(i,j+1) - dudxn*mu1Y(i,j+1)

            oper = alpha*a(i,j)*U(i,j) + beta*(hx*(xfluxe-xfluxw)+hy*(yfluxn-yfluxs))
            diag = alpha*a(i,j) + 
     $           beta*( muX(i  ,j)*hxx
     $           +      muX(i+1,j)*hxx*(1.d0 + de)
     $           +      muY(i,j  )*hyy - mu1Y(i  ,j)*hy*ds
     $           +      muY(i,j+1)*hyy + mu1Y(i,j+1)*hy*dn)
            U(i,j) = U(i,j) + fac*(rhs(i,j) - oper)/diag
         enddo
      endif

c ::: case 3a: north-east corner
      j=hi(2)
      i=hi(1)
      if( mod(i,2).eq.modx .and. mod(j,2).eq.mody) then
         dudxw = (U(i  ,j  ) - U(i-1,j  ))*hx
         dudxe = (U(i+1,j  ) - U(i  ,j  ))*hx
         dudys = (U(i  ,j  ) - U(i  ,j-1))*hy
         dudyn = (U(i  ,j+1) - U(i  ,j  ))*hy
         
         dudxn = trandern(i,j+1)
         dudye = trandere(i+1,j)

         if( maske(i+1,j-1).gt.0 .or. maske(i+1,j).gt.0 ) then
            dudxs = (U(i-2,j  ) - 4*U(i-1,j  ) + 3*U(i  ,j  )
     &           +   U(i-2,j-1) - 4*U(i-1,j-1) + 3*U(i  ,j-1))*hx4
            ds = 3*hx4
         else
            dudxs = (U(i+1,j  ) + U(i+1,j-1) - U(i-1,j  ) - U(i-1,j-1))*hx4
            ds = 0.d0
         endif
         if( maskn(i-1,j+1).gt.0 .or. maskn(i,j+1).gt.0 ) then
            dudyw = (U(i-1,j-2) - 4*U(i-1,j-1) + 3*U(i-1,j  )
     &           +   U(i  ,j-2) - 4*U(i  ,j-1) + 3*U(i  ,j  ))*hy4
            dw = 3*hy4
         else
            dudyw = (U(i-1,j+1) + U(i  ,j+1) - U(i-1,j-1) - U(i  ,j-1))*hy4
            dw = 0.d0
         endif
         if( maske(i+1,j  ).gt. 0) then
            de = - fe(i,j)
         else
            de = 0.d0
         endif
         if( maskn(i,j+1).gt. 0) then
            dn = - fn(i,j)
         else
            dn = 0.d0
         endif
         
         xfluxw = - dudxw*muX(i  ,j) - dudyw*mu1X(i  ,j)
         xfluxe = - dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
         yfluxs = - dudys*muY(i  ,j) - dudxs*mu1Y(i,j  )
         yfluxn = - dudyn*muY(i,j+1) - dudxn*mu1Y(i,j+1)
         
         oper = alpha*a(i,j)*U(i,j) + beta*(hx*(xfluxe-xfluxw)+hy*(yfluxn-yfluxs))
         diag = alpha*a(i,j) + 
     $        beta*( muX(i  ,j)*hxx - mu1X(i,j  )*hx*dw
     $        +      muX(i+1,j)*hxx*(1.d0 + de)
     $        +      muY(i,j  )*hyy - mu1Y(i  ,j)*hy*ds
     $        +      muY(i,j+1)*hyy*(1.d0 + dn))
         U(i,j) = U(i,j) + fac*(rhs(i,j) - oper)/diag
      endif


c ::: case 3b: south-east corner
      j=lo(2)
      i=hi(1)      
      if( mod(i,2).eq.modx .and. mod(j,2).eq.mody) then
         dudxw = (U(i  ,j  ) - U(i-1,j  ))*hx
         dudxe = (U(i+1,j  ) - U(i  ,j  ))*hx
         dudys = (U(i  ,j  ) - U(i  ,j-1))*hy
         dudyn = (U(i  ,j+1) - U(i  ,j  ))*hy
         
         dudxs = tranders(i,j-1)
         dudye = trandere(i+1,j)

         if( maske(i+1,j+1).gt.0 .or. maske(i+1,j).gt.0 ) then
            dudxn = (U(i-2,j  ) - 4*U(i-1,j  ) + 3*U(i  ,j  )
     &           +   U(i-2,j+1) - 4*U(i-1,j+1) + 3*U(i  ,j+1))*hx4
            dn = 3*hx4
         else
            dudxn = (U(i+1,j  ) + U(i+1,j+1) - U(i-1,j  ) - U(i-1,j+1))*hx4
            dn = 0.d0
         endif
         if( masks(i-1,j-1).gt.0 .or. masks(i,j-1).gt.0 ) then
            dudyw = (-U(i-1,j+2) + 4*U(i-1,j+1) - 3*U(i-1,j  )
     &           -    U(i  ,j+2) + 4*U(i  ,j+1) - 3*U(i  ,j  ))*hy4
            dw = - 3*hy4
         else
            dudyw = (U(i-1,j+1) + U(i  ,j+1) - U(i-1,j-1) - U(i  ,j-1))*hy4
            dw = 0.d0
         endif
         if( maske(i+1,j  ).gt. 0) then
            de = - fe(i,j)
         else
            de = 0.d0
         endif
         if( masks(i,j-1).gt. 0) then
            ds = - fs(i,j)
         else
            ds = 0.d0
         endif
         
         xfluxw = - dudxw*muX(i  ,j) - dudyw*mu1X(i  ,j)
         xfluxe = - dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
         yfluxs = - dudys*muY(i  ,j) - dudxs*mu1Y(i,j  )
         yfluxn = - dudyn*muY(i,j+1) - dudxn*mu1Y(i,j+1)
         
         oper = alpha*a(i,j)*U(i,j) + beta*(hx*(xfluxe-xfluxw)+hy*(yfluxn-yfluxs))
         diag = alpha*a(i,j) + 
     $        beta*( muX(i  ,j)*hxx - mu1X(i,j  )*hx*dw
     $        +      muX(i+1,j)*hxx*(1.d0 + de)
     $        +      muY(i,j  )*hyy*(1.d0 + ds)
     $        +      muY(i,j+1)*hyy + mu1Y(i,j+1)*hy*dn)
         U(i,j) = U(i,j) + fac*(rhs(i,j) - oper)/diag
      endif

c ::: case 3c: south-west corner
      j=lo(2)
      i=lo(1)
      if( mod(i,2).eq.modx .and. mod(j,2).eq.mody) then
         dudxw = (U(i  ,j  ) - U(i-1,j  ))*hx
         dudxe = (U(i+1,j  ) - U(i  ,j  ))*hx
         dudys = (U(i  ,j  ) - U(i  ,j-1))*hy
         dudyn = (U(i  ,j+1) - U(i  ,j  ))*hy
         
         dudxs = tranders(i,j-1)
         dudyw = tranderw(i-1,j)

         if( maskw(i-1,j+1).gt.0 .or. maskw(i-1,j).gt.0 ) then
            dudxn = (-U(i+2,j  ) + 4*U(i+1,j  ) - 3*U(i  ,j  )
     &           -    U(i+2,j+1) + 4*U(i+1,j+1) - 3*U(i  ,j+1))*hx4
            dn = - 3*hx4
         else
            dudxn = (U(i+1,j  ) + U(i+1,j+1) - U(i-1,j  ) - U(i-1,j+1))*hx4
            dn = 0.d0
         endif
         if( masks(i  ,j-1).gt.0 .or. masks(i+1,j-1).gt.0 ) then
            dudye = (-U(i-1,j+2) + 4*U(i-1,j+1) - 3*U(i-1,j  )
     &           -    U(i  ,j+2) + 4*U(i  ,j+1) - 3*U(i  ,j  ))*hy4
            de = - 3*hy4
         else
            dudye = (U(i-1,j+1) + U(i  ,j+1) - U(i-1,j-1) - U(i  ,j-1))*hy4
            de = 0.d0
         endif
         if( maskw(i-1,j  ).gt. 0) then
            dw = - fw(i,j)
         else
            dw = 0.d0
         endif
         if( masks(i,j-1).gt. 0) then
            ds = - fs(i,j)
         else
            ds = 0.d0
         endif
         
         xfluxw = - dudxw*muX(i  ,j) - dudyw*mu1X(i  ,j)
         xfluxe = - dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
         yfluxs = - dudys*muY(i  ,j) - dudxs*mu1Y(i,j  )
         yfluxn = - dudyn*muY(i,j+1) - dudxn*mu1Y(i,j+1)
         
         oper = alpha*a(i,j)*U(i,j) + beta*(hx*(xfluxe-xfluxw)+hy*(yfluxn-yfluxs))
         diag = alpha*a(i,j) + 
     $        beta*( muX(i  ,j)*hxx*(1.d0 + dw)
     $        +      muX(i+1,j)*hxx + mu1X(i+1,j)*hx*de
     $        +      muY(i,j  )*hyy*(1.d0 + ds)
     $        +      muY(i,j+1)*hyy + mu1Y(i,j+1)*hy*dn)
         U(i,j) = U(i,j) + fac*(rhs(i,j) - oper)/diag
      endif

c ::: case 3d: north-west corner
      j=hi(2)
      i=lo(1)
      if( mod(i,2).eq.modx .and. mod(j,2).eq.mody) then
         dudxw = (U(i  ,j  ) - U(i-1,j  ))*hx
         dudxe = (U(i+1,j  ) - U(i  ,j  ))*hx
         dudys = (U(i  ,j  ) - U(i  ,j-1))*hy
         dudyn = (U(i  ,j+1) - U(i  ,j  ))*hy
         
         dudxn = trandern(i,j+1)
         dudyw = tranderw(i-1,j)

         if( maskw(i-1,j-1).gt.0 .or. maskw(i-1,j).gt.0 ) then
            dudxs = (-U(i+2,j  ) + 4*U(i+1,j  ) - 3*U(i  ,j  )
     &           -    U(i+2,j-1) + 4*U(i+1,j-1) - 3*U(i  ,j-1))*hx4
            ds = - 3*hx4
         else
            dudxs = (U(i+1,j  ) + U(i+1,j-1) - U(i-1,j  ) - U(i-1,j-1))*hx4
            ds = 0.d0
         endif
         if( maskn(i  ,j+1).gt.0 .or. maskn(i+1,j+1).gt.0 ) then
            dudye = ( U(i+1,j-2) - 4*U(i+1,j-1) + 3*U(i+1,j  )
     &           +    U(i  ,j-2) - 4*U(i  ,j-1) + 3*U(i  ,j  ))*hy4
            de = 3*hy4
         else
            dudye = (U(i+1,j+1) + U(i  ,j+1) - U(i+1,j-1) - U(i  ,j-1))*hy4
            de = 0.d0
         endif
         if( maskw(i-1,j  ).gt. 0) then
            dw = - fw(i,j)
         else
            dw = 0.d0
         endif
         if( maskn(i,j+1).gt. 0) then
            dn = - fn(i,j)
         else
            dn = 0.d0
         endif
         
         xfluxw = - dudxw*muX(i  ,j) - dudyw*mu1X(i  ,j)
         xfluxe = - dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
         yfluxs = - dudys*muY(i  ,j) - dudxs*mu1Y(i,j  )
         yfluxn = - dudyn*muY(i,j+1) - dudxn*mu1Y(i,j+1)
         
         oper = alpha*a(i,j)*U(i,j) + beta*(hx*(xfluxe-xfluxw)+hy*(yfluxn-yfluxs))
         diag = alpha*a(i,j) + 
     $        beta*( muX(i  ,j)*hxx*(1.d0 + dw)
     $        +      muX(i+1,j)*hxx + mu1X(i+1,j)*hx*de
     $        +      muY(i,j  )*hyy - mu1Y(i,j  )*hy*ds
     $        +      muY(i,j+1)*hyy*(1.d0 + dn))        
         U(i,j) = U(i,j) + fac*(rhs(i,j) - oper)/diag
      endif
      end
