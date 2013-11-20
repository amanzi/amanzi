#include <REAL.H>

#include "TensorOp_F.H"
#include <LO_BCTYPES.H>
#include "ArrayLim.H"

      subroutine FORT_TOAPPLY(
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
      if( masks(i+1,j-1).gt.0 .or. masks(i,j-1).gt.0 ) then
         dudye = (-U(i+1,j+2) + 4*U(i+1,j+1) - 3*U(i+1,j  )
     &        -    U(i  ,j+2) + 4*U(i  ,j+1) - 3*U(i  ,j  ))*hy4
      else
         dudye = (U(i+1,j+1) + U(i  ,j+1) - U(i+1,j-1) - U(i  ,j-1))*hy4
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


      subroutine FORT_TOFLUX(
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
      if( maskn(i+1,j+1).gt.0 .or. maskn(i,j+1).gt.0 ) then
         dudye = (U(i+1,j-2) - 4*U(i+1,j-1) + 3*U(i+1,j  )
     &        +   U(i  ,j-2) - 4*U(i  ,j-1) + 3*U(i  ,j  ))*hy4
      else
         dudye = (U(i  ,j+1) + U(i+1,j+1) - U(i  ,j-1) - U(i+1,j-1))*hy4
      endif
      xflux(i,j) = -dudxw*muX(i,j) - dudyw*mu1X(i,j)
      yflux(i,j) = -dudxs*mu1Y(i,j) - dudys*muY(i,j)
      yflux(i,j+1) = -dudxn*mu1Y(i,j+1) - dudyn*muY(i,j+1)

      end

      subroutine FORT_TOGSRB (
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
     $     ftn, DIMS(ftn),
     $     maske,DIMS(maske),
     $     fe, DIMS(fe),
     $     fte, DIMS(fte),
     $     maskw,DIMS(maskw),
     $     fw, DIMS(fw),
     $     ftw, DIMS(ftw),
     $     masks,DIMS(masks),
     $     fs, DIMS(fs),
     $     fts, DIMS(fts),
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
      integer DIMDEC(ftn)
      integer DIMDEC(fte)
      integer DIMDEC(ftw)
      integer DIMDEC(fts)
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
      REAL_T ftn(DIMV(ftn))
      REAL_T fte(DIMV(fte))
      REAL_T ftw(DIMV(ftw))
      REAL_T fts(DIMV(fts))
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
      REAL_T dtw, dte, dts, dtn


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
         write(6,*)'FORT_TOGSRB: bad phaseflag', phaseflag
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
            U(i,j) = U(i,j) + (rhs(i,j) - oper)/diag
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
	       dn  = - fn(i,j)
	       dtn =   ftn(i,j)
            else
               dn  = 0.d0
               dtn = 0.d0
            endif

            xfluxw = - dudxw*muX(i  ,j) - dudyw*mu1X(i  ,j)
            xfluxe = - dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
            yfluxs = - dudys*muY(i  ,j) - dudxs*mu1Y(i,j  )
            yfluxn = - dudyn*muY(i,j+1) - dudxn*mu1Y(i,j+1)

            oper = alpha*a(i,j)*U(i,j) + beta*(hx*(xfluxe-xfluxw)+hy*(yfluxn-yfluxs))
            diag = alpha*a(i,j) + 
     $           beta*( muX(i  ,j)*hxx             + mu1X(i  ,j)*hx*dw
     $           +      muX(i+1,j)*hxx             - mu1X(i+1,j)*hx*de
     $           +      muY(i,j  )*hyy
     $           +      muY(i,j+1)*hyy*(1.d0 + dn) - mu1Y(i,j+1)*hy*dtn )
            U(i,j) = U(i,j) + (rhs(i,j) - oper)/diag
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
	       ds  = - fs(i,j)
	       dts =   fts(i,j)
            else
               ds  = 0.d0
               dts = 0.d0
            endif

            xfluxw = - dudxw*muX(i  ,j) - dudyw*mu1X(i  ,j)
            xfluxe = - dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
            yfluxs = - dudys*muY(i  ,j) - dudxs*mu1Y(i,j  )
            yfluxn = - dudyn*muY(i,j+1) - dudxn*mu1Y(i,j+1)

            oper = alpha*a(i,j)*U(i,j) + beta*(hx*(xfluxe-xfluxw)+hy*(yfluxn-yfluxs))
            diag = alpha*a(i,j) + 
     $           beta*( muX(i  ,j)*hxx             + mu1X(i  ,j)*hx*dw
     $           +      muX(i+1,j)*hxx             - mu1X(i+1,j)*hx*de
     $           +      muY(i,j  )*hyy*(1.d0 + ds) + mu1Y(i  ,j)*hy*dts
     $           +      muY(i,j+1)*hyy )
            U(i,j) = U(i,j) + (rhs(i,j) - oper)/diag
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
	       dw  = - fw(i,j)
	       dtw =   ftw(i,j)
            else
               dw  = 0.d0
               dtw = 0.d0
            endif

            xfluxw = - dudxw*muX(i  ,j) - dudyw*mu1X(i  ,j)
            xfluxe = - dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
            yfluxs = - dudys*muY(i  ,j) - dudxs*mu1Y(i,j  )
            yfluxn = - dudyn*muY(i,j+1) - dudxn*mu1Y(i,j+1)

            oper = alpha*a(i,j)*U(i,j) + beta*(hx*(xfluxe-xfluxw)+hy*(yfluxn-yfluxs))
            diag = alpha*a(i,j) + 
     $           beta*( muX(i  ,j)*hxx*(1.d0 + dw) + mu1X(i  ,j)*hx*dtw
     $           +      muX(i+1,j)*hxx
     $           +      muY(i,j  )*hyy             + mu1Y(i  ,j)*hy*ds
     $           +      muY(i,j+1)*hyy             - mu1Y(i,j+1)*hy*dn )
            U(i,j) = U(i,j) + (rhs(i,j) - oper)/diag
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
	       de  = - fe(i,j)
	       dte =   fte(i,j)
            else
               de  = 0.d0
               dte = 0.d0
            endif

            xfluxw = - dudxw*muX(i  ,j) - dudyw*mu1X(i  ,j)
            xfluxe = - dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
            yfluxs = - dudys*muY(i  ,j) - dudxs*mu1Y(i,j  )
            yfluxn = - dudyn*muY(i,j+1) - dudxn*mu1Y(i,j+1)

            oper = alpha*a(i,j)*U(i,j) + beta*(hx*(xfluxe-xfluxw)+hy*(yfluxn-yfluxs))
            diag = alpha*a(i,j) + 
     $           beta*( muX(i  ,j)*hxx
     $           +      muX(i+1,j)*hxx*(1.d0 + de) - mu1X(i+1,j)*hx*dte
     $           +      muY(i,j  )*hyy             + mu1Y(i  ,j)*hy*ds
     $           +      muY(i,j+1)*hyy             - mu1Y(i,j+1)*hy*dn)
            U(i,j) = U(i,j) + (rhs(i,j) - oper)/diag
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
            de  = - fe(i,j)
            dte =   fte(i,j)
         else
            de  = 0.d0
            dte = 0.d0
         endif
         if( maskn(i,j+1).gt. 0) then
            dn  = - fn(i,j)
            dtn =   ftn(i,j)
         else
            dn  = 0.d0
            dtn = 0.d0
         endif
         
         xfluxw = - dudxw*muX(i  ,j) - dudyw*mu1X(i  ,j)
         xfluxe = - dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
         yfluxs = - dudys*muY(i  ,j) - dudxs*mu1Y(i,j  )
         yfluxn = - dudyn*muY(i,j+1) - dudxn*mu1Y(i,j+1)
         
         oper = alpha*a(i,j)*U(i,j) + beta*(hx*(xfluxe-xfluxw)+hy*(yfluxn-yfluxs))
         diag = alpha*a(i,j) + 
     $        beta*( muX(i  ,j)*hxx             + mu1X(i  ,j)*hx*dw
     $        +      muX(i+1,j)*hxx*(1.d0 + de) - mu1X(i+1,j)*hx*dte
     $        +      muY(i,j  )*hyy             + mu1Y(i  ,j)*hy*ds
     $        +      muY(i,j+1)*hyy*(1.d0 + dn) - mu1Y(i,j+1)*hy*dtn )
         U(i,j) = U(i,j) + (rhs(i,j) - oper)/diag
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
            de  = - fe(i,j)
            dte =   fte(i,j)
         else
            de  = 0.d0
            dte = 0.d0
         endif
         if( masks(i,j-1).gt. 0) then
            ds  = - fs(i,j)
            dts =   fts(i,j)
         else
            ds  = 0.d0
            dts = 0.d0
         endif
         
         xfluxw = - dudxw*muX(i  ,j) - dudyw*mu1X(i  ,j)
         xfluxe = - dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
         yfluxs = - dudys*muY(i  ,j) - dudxs*mu1Y(i,j  )
         yfluxn = - dudyn*muY(i,j+1) - dudxn*mu1Y(i,j+1)
         
         oper = alpha*a(i,j)*U(i,j) + beta*(hx*(xfluxe-xfluxw)+hy*(yfluxn-yfluxs))
         diag = alpha*a(i,j) + 
     $        beta*( muX(i  ,j)*hxx             + mu1X(i  ,j)*hx*dw
     $        +      muX(i+1,j)*hxx*(1.d0 + de) - mu1X(i+1,j)*hx*dte
     $        +      muY(i,j  )*hyy*(1.d0 + ds) + mu1Y(i,j  )*hy*dts
     $        +      muY(i,j+1)*hyy             - mu1Y(i,j+1)*hy*dn )
         U(i,j) = U(i,j) + (rhs(i,j) - oper)/diag
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
            dudye = (-U(i+1,j+2) + 4*U(i+1,j+1) - 3*U(i+1,j  )
     &           -    U(i  ,j+2) + 4*U(i  ,j+1) - 3*U(i  ,j  ))*hy4
            de = - 3*hy4
         else
            dudye = (U(i+1,j+1) + U(i  ,j+1) - U(i+1,j-1) - U(i  ,j-1))*hy4
            de = 0.d0
         endif
         if( maskw(i-1,j  ).gt. 0) then
            dw  = - fw(i,j)
            dtw =   ftw(i,j)
         else
            dw  = 0.d0
            dtw = 0.d0
         endif
         if( masks(i,j-1).gt. 0) then
            ds  = - fs(i,j)
            dts =   fts(i,j)
         else
            ds  = 0.d0
            dts = 0.d0
         endif
         
         xfluxw = - dudxw*muX(i  ,j) - dudyw*mu1X(i  ,j)
         xfluxe = - dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
         yfluxs = - dudys*muY(i  ,j) - dudxs*mu1Y(i,j  )
         yfluxn = - dudyn*muY(i,j+1) - dudxn*mu1Y(i,j+1)
         
         oper = alpha*a(i,j)*U(i,j) + beta*(hx*(xfluxe-xfluxw)+hy*(yfluxn-yfluxs))
         diag = alpha*a(i,j) + 
     $        beta*( muX(i  ,j)*hxx*(1.d0 + dw) + mu1X(i  ,j)*hx*dtw
     $        +      muX(i+1,j)*hxx             - mu1X(i+1,j)*hx*de
     $        +      muY(i,j  )*hyy*(1.d0 + ds) + mu1Y(i,j  )*hy*dts
     $        +      muY(i,j+1)*hyy             - mu1Y(i,j+1)*hy*dn)
         U(i,j) = U(i,j) + (rhs(i,j) - oper)/diag
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
            dw  = - fw(i,j)
            dtw =   ftw(i,j)
         else
            dw  = 0.d0
            dtw = 0.d0
         endif
         if( maskn(i,j+1).gt. 0) then
            dn  = - fn(i,j)
            dtn =   ftn(i,j)
         else
            dn  = 0.d0
            dtn = 0.d0
         endif
         
         xfluxw = - dudxw*muX(i  ,j) - dudyw*mu1X(i  ,j)
         xfluxe = - dudxe*muX(i+1,j) - dudye*mu1X(i+1,j)
         yfluxs = - dudys*muY(i  ,j) - dudxs*mu1Y(i,j  )
         yfluxn = - dudyn*muY(i,j+1) - dudxn*mu1Y(i,j+1)
         
         oper = alpha*a(i,j)*U(i,j) + beta*(hx*(xfluxe-xfluxw)+hy*(yfluxn-yfluxs))
         diag = alpha*a(i,j) + 
     $        beta*( muX(i  ,j)*hxx*(1.d0 + dw) + mu1X(i  ,j)*hx*dtw
     $        +      muX(i+1,j)*hxx             - mu1X(i+1,j)*hx*de
     $        +      muY(i,j  )*hyy             + mu1Y(i,j  )*hy*ds
     $        +      muY(i,j+1)*hyy*(1.d0 + dn) - mu1Y(i,j+1)*hy*dtn )
         U(i,j) = U(i,j) + (rhs(i,j) - oper)/diag
      endif
      end

      subroutine FORT_TOAPPLYBC (
     $     flagden, flagbc, maxorder,
     $     u,   DIMS(u),
     $     cdir, bct, bcl,
     $     bcval, DIMS(bcval),
     $     mask,  DIMS(mask),
     $     maskp, DIMS(maskp),
     $     maskm, DIMS(maskm),
     $     den,   DIMS(den),
     $     dent,  DIMS(dent),
     $     exttd, DIMS(exttd),
     $     tander,DIMS(tander),
     $     lo, hi, nc,
     $     h
     $     )
c
c     If the boundary is of Neumann type, set the ghost cell value to
c     that of the outermost point in the valid data (2nd order accurate)
c     and then fill the "den" array with the value "1"
c     
c     
c     If flagbc==1:
c     
c     If the boundary is of Dirichlet type, construct a polynomial
c     interpolation through the boundary location and internal points
c     (at locations x(-1:len-2) that generates the ghost cell value (at
c     location xInt).  Then fill the ghost cell with the interpolated value.
c     If flagden==1, load the "den" array with the interpolation
c     coefficient corresponding to outermost point in the valid region
c     ( the coef(0) corresponding to the location x(0) )
c      
c     Note: 
c     The bc type = LO_REFLECT_ODD is a special type of boundary condition.

c ::: the old FORT_APPLYBC put values out in ghost cells.  The new
c ::: FORT_APPLYBC makes a distinction between normal and tangential derivs.
c ::: the normal derivatives are still evaluated with ghost cells ( and
c ::: the den array for adjusting the diagonal element).  Tangential
c ::: derivatives are NOT computed from the ghost cells (note: we are
c ::: discussing here the tangential derivates which are centered on the
c ::: external cells walls of the rectangular domain.  Internal tangential
c ::: derivatives whose stencil extends outside the rectangular domain have
c ::: still ANOTHER calculational trick, one sided derivatives, which is 
c ::: implemented in the apply operator).  For these tangential derivatives,
c ::: if some of the cells in the stencil are uncovered, the derivative is
c ::: computed half from cells which are internal, plus an externally supplied
c ::: tangential derivative centered at location bcl.
c ??? exttd is the externally supplied tangential derivative which is
c ??? centered at bcl from wall.  Since the location of exttd is
c ??? indeterminate (determined by bcl), we arbitrarily index it in the
c ??? cell just outside the wall.
c ??? tander is the calculated tangential derivative which is centered at
c ??? at the wall.  It is indexed as if it were a cell centered quantity
c ??? just outside the wall.  This is anomalous since it is really edge
c ??? centered on the outer wall.
c      
c ::: other notes since previous developers didn't bother to document
c ::: cdir is mnemonic for coordinate direction, i.e. which side this is
c ::: cdir==0 -> left wall
c ::: cdir==2 -> right
c ::: cdir==1 -> bottom
c ::: cdir==3 -> top

c ::: notes about mask
c ::: since we now require other masks, there are 3.
c ::: 1) mask -- the main one covering the cells off the side of valid domain.
c ::: 2) maskp -- the plus end of transverse direction
c ::: 3) maskn -- the minus end of transverse direction

      integer maxorder
      integer nc, cdir, flagden, flagbc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(u)
      REAL_T u(DIMV(u),nc)
      integer DIMDEC(den)
      REAL_T den(DIMV(den),nc)
      integer DIMDEC(dent)
      REAL_T dent(DIMV(dent),nc)
      integer DIMDEC(exttd)
      REAL_T exttd(DIMV(exttd),nc)
      integer DIMDEC(tander)
      REAL_T tander(DIMV(tander),nc)
      integer DIMDEC(bcval)
      REAL_T bcval(DIMV(bcval),nc)
      integer DIMDEC(mask)
      integer mask(DIMV(mask))
      integer DIMDEC(maskp)
      integer maskp(DIMV(maskp))
      integer DIMDEC(maskm)
      integer maskm(DIMV(maskm))
      integer bct(nc)
      REAL_T bcl
      REAL_T h(BL_SPACEDIM)
c
      REAL_T hx,hy,lambda
      REAL_T innder,outder
      REAL_T innloc,outloc,edgloc,innden
      integer i
      integer j
      integer n
      logical is_dirichlet
      logical is_neumann
c
      integer lenx
      integer leny
      integer m
c
      integer Lmaxorder
      integer maxmaxorder
      parameter(maxmaxorder=4)
      REAL_T x(-1:maxmaxorder-2)
      REAL_T coef(-1:maxmaxorder-2)
      REAL_T xInt
c
      logical False, True
      parameter( False=.false.)
      parameter( True=.true. )
c
      is_dirichlet(i) = ( i .eq. LO_DIRICHLET   )
      is_neumann(i)   = ( i .eq. LO_NEUMANN )
      hx = h(1)
      hy = h(2)
c
      if ( maxorder .eq. -1 ) then
         Lmaxorder = maxmaxorder
      else
         Lmaxorder = MIN(maxorder,maxmaxorder)
      endif
      lenx = MIN(hi(1)-lo(1), Lmaxorder-2)
      leny = MIN(hi(2)-lo(2), Lmaxorder-2)
c
c     TODO:
c     In order for this to work with growing multigrid, must
c     sort xa[] because it is possible for the xb value to lay
c     within this range.
c     
c     The Left face of the grid
c
      if(cdir .eq. 0) then
         do n = 1, nc
            if (is_neumann(bct(n))) then
               do j = lo(2), hi(2)
                  u(lo(1)-1,j,n) = merge(
     $                 u(lo(1),j,n),
     $                 u(lo(1)-1,j,n),
     $                 mask(lo(1)-1,j) .gt. 0)
               enddo
               if ( flagden .eq. 1) then
                  do j = lo(2), hi(2)
                     den(lo(1),j,n) = 1.0
                  enddo
               endif
c     :::::: now do the tangential derivative part
c     :::::: Null
c     :::::: interior part of left side
               i = lo(1)
               do j=lo(2)+1,hi(2)-1
                  innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
                  tander(i-1,j,n) = innder
                  if ( flagden .eq. 1) then
                     dent(i,j,n) = 0.d0
                  endif
               enddo
c     :::::: now end points
c     :::::: lower leftside
               j = lo(2)
               if( maskm(i,-1 + j).gt.0 ) then
                  innder = (-3*U(i,j,n) + 4*U(i,1 + j,n) - U(i,2 + j,n))/(2.*hy)
                  innden = - 3.d0/(2.*hy)
               else
                  innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
                  innden = 0.d0
               endif
               tander(i-1,j,n) = innder
               if ( flagden .eq. 1) then
                  dent(i,j,n) = innden
               endif
c     :::::: upper left side
                j = hi(2)
                if( maskp(i,1 + j).gt.0 ) then
                   innder = (U(i,-2 + j,n) - 4*U(i,-1 + j,n) + 3*U(i,j,n))/(2.*hy)
                   innden = 3.d0/(2.*hy)
                else
                   innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
                   innden = 0.d0
                endif
                tander(i-1,j,n) = innder
                if ( flagden .eq. 1) then
                   dent(i,j,n) = innden
                endif
            else if (is_dirichlet(bct(n))) then
               do m=0,lenx
                  x(m) = m + 0.5
               enddo
               x(-1) = - bcl/h(1)
               xInt = - 0.5
               call polyInterpCoeff(xInt, x, lenx+2, coef)
               if ( flagbc .eq. 1 ) then
                  do j = lo(2), hi(2)
                     u(lo(1)-1, j, n) = merge(
     $                    bcval(lo(1)-1,j,n)*coef(-1),
     $                    u(lo(1)-1, j, n),
     $                    mask(lo(1)-1,j) .gt. 0)
                  enddo
               else
                  do j = lo(2), hi(2)
                     u(lo(1)-1, j, n) = merge(
     $                    0.0d0,
     $                    u(lo(1)-1, j, n),
     $                    mask(lo(1)-1,j) .gt. 0)
                  enddo
               endif
               do m = 0, lenx
                  do j = lo(2), hi(2)
                     u(lo(1)-1,j,n) = merge(
     $                    u(lo(1)-1,j,n)
     $                    + u(lo(1)+m, j, n)*coef(m),
     $                    u(lo(1)-1,j,n),
     $                    mask(lo(1)-1,j) .gt. 0)
                  enddo
               enddo
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     den(lo(1),j,n) = merge(coef(0), 0.0d0,
     $                    mask(lo(1)-1,j) .gt. 0)
                  enddo
               endif
c     :::::: now do the tangential derivative part
c     :::::: Null
c     :::::: interior part of left side
               edgloc = 0.
               i = lo(1)
               do j=lo(2)+1,hi(2)-1
                  if( mask(-1 + i,-1 + j).eq.0.and.mask(-1 + i,1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-U(-1 + i,-1 + j,n) + U(-1 + i,1 + j,n))/(2.*hy)
                  elseif( mask(-1 + i,j).eq.0.and.mask(-1 + i,1 + j).eq.0.and.mask(-1 + i,2 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-3*U(-1 + i,j,n) + 4*U(-1 + i,1 + j,n) - U(-1 + i,2 + j,n))/(2.*hy)
                  elseif( mask(-1 + i,-2 + j).eq.0.and.mask(-1 + i,-1 + j).eq.0.and.mask(-1 + i,j).eq.0 ) then
                     outloc = -0.5
                     outder = (U(-1 + i,-2 + j,n) - 4*U(-1 + i,-1 + j,n) + 3*U(-1 + i,j,n))/(2.*hy)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1) then
                        outder = exttd(i-1,j,n)
                     else
                        outder = 0.
                     endif
                  endif
                  innloc = 0.5
                  innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
            
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  tander(i-1,j,n) = lambda*innder+(1-lambda)*outder
                  if (flagden .eq. 1 ) then
                     dent(i,j,n) = 0.d0
                  endif
               enddo
c     :::::: now end points
c     :::::: lower leftside
               j = lo(2)
               if( mask(-1 + i,-1 + j).eq.0.and.mask(-1 + i,1 + j).eq.0 ) then
                   outloc = -0.5
                   outder = (-U(-1 + i,-1 + j,n) + U(-1 + i,1 + j,n))/(2.*hy)
               elseif( mask(-1 + i,j).eq.0.and.mask(-1 + i,1 + j).eq.0.and.mask(-1 + i,2 + j).eq.0 ) then
                   outloc = -0.5
                   outder = (-3*U(-1 + i,j,n) + 4*U(-1 + i,1 + j,n) - U(-1 + i,2 + j,n))/(2.*hy)
               else
                   outloc = x(-1)
                   if( flagbc .eq. 1) then
                      outder = exttd(i-1,j,n)
                   else
                      outder = 0.
                   endif
               endif
               if( maskm(i,-1 + j).gt.0 ) then
                  innder = (-3*U(i,j,n) + 4*U(i,1 + j,n) - U(i,2 + j,n))/(2.*hy)
                  innden = -3.d0/(2.*hy)
               else
                  innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
                  innden = 0.d0
               endif
               innloc = 0.5
               lambda = (edgloc-outloc)/(innloc-outloc)
               tander(i-1,j,n) = lambda*innder+(1-lambda)*outder
               if (flagden .eq. 1) then
                  dent(i,j,n) = lambda*innden
               endif
c     :::::: upper left side
                j = hi(2)
                if( mask(-1 + i,-1 + j).eq.0.and.mask(-1 + i,1 + j).eq.0 ) then
                   outloc = -0.5
                   outder = (-U(-1 + i,-1 + j,n) + U(-1 + i,1 + j,n))/(2.*hy)
                elseif( mask(-1 + i,-2 + j).eq.0.and.mask(-1 + i,-1 + j).eq.0.and.mask(-1 + i,j).eq.0 ) then
                   outloc = -0.5
                   outder = (U(-1 + i,-2 + j,n) - 4*U(-1 + i,-1 + j,n) + 3*U(-1 + i,j,n))/(2.*hy)
                else
                   outloc = x(-1)
                   if( flagbc .eq. 1) then
                      outder = exttd(i-1,j,n)
                   else
                      outder = 0.
                   endif
                endif
                if( maskp(i,1 + j).gt.0 ) then
                   innder = (U(i,-2 + j,n) - 4*U(i,-1 + j,n) + 3*U(i,j,n))/(2.*hy)
                   innden = 3.d0/(2.*hy)
                else
                   innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
                   innden = 0.d0
                endif
                innloc = 0.5
                lambda = (edgloc-outloc)/(innloc-outloc)
                tander(i-1,j,n) = lambda*innder+(1-lambda)*outder
                if (flagden .eq. 1) then
                   dent(i,j,n) = lambda*innden
                endif
               
            else if ( bct(n) .eq. LO_REFLECT_ODD ) then
               
               do j = lo(2), hi(2)
                  u(lo(1)-1, j, n) = merge(
     $                 -u(lo(1),j,n),
     $                 u(lo(1)-1, j, n),
     $                 mask(lo(1)-1,j) .gt. 0)
                  tander(lo(1)-1,j,n) = 0.0
               enddo
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     den(lo(1),j,n) = merge(-1.0d0, 0.0d0,
     $                    mask(lo(1)-1,j) .gt. 0)
                     dent(lo(1),j,n) = 0.d0
                  enddo
               endif
               
            else 
               print *,'UNKNOWN BC ON LEFT FACE IN APPLYBC'
               stop
            endif
         enddo
      endif
c     
c     The Right face of the grid
c
      if(cdir .eq. 2) then
         do n = 1, nc
            if(is_neumann(bct(n))) then
               do j = lo(2), hi(2)
                  u(hi(1)+1,j,n) = merge(
     $                 u(hi(1), j, n),
     $                 u(hi(1)+1, j, n),
     $                 mask(hi(1)+1,j) .gt. 0)
               enddo
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     den(hi(1),j,n) = 1.0
                  enddo
               endif
c     :::::: now do the tangential derivative part
c     :::::: Null
c     :::::: interior part of right side
               i = hi(1)
               do j=lo(2)+1,hi(2)-1
                  innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
                  tander(i+1,j,n) = innder
                  if (flagden .eq. 1) then
                     dent(i,j,n) = 0.d0
                  endif
               enddo
c :::::: now end points
c :::::: lower right side
               j = lo(2)
               if( maskm(i,-1 + j).gt.0 ) then
                  innder = (-3*U(i,j,n) + 4*U(i,1 + j,n) - U(i,2 + j,n))/(2.*hy)
                  innden = -3.d0/(2.*hy)
               else
                  innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
                  innden = 0.d0
               endif
               tander(i+1,j,n) = innder
               if (flagden .eq. 1) then
                  dent(i,j,n) = innden
               endif
c     :::::: upper right side
               j = hi(2)
               if( maskp(i,1 + j).gt.0 ) then
                  innder = (U(i,-2 + j,n) - 4*U(i,-1 + j,n) + 3*U(i,j,n))/(2.*hy)
                  innden = 3.d0/(2.*hy)
               else
                  innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
                  innden = 0.d0
               endif
               tander(i+1,j,n) = innder
               if (flagden .eq. 1) then
                  dent(i,j,n) = innden
               endif
            else if (is_dirichlet(bct(n))) then
               do m=0,lenx
                  x(m) = m + 0.5
               enddo
               x(-1) = - bcl/h(1)
               xInt = - 0.5
               call polyInterpCoeff(xInt, x, lenx+2, coef)
               if ( flagbc .eq. 1 ) then
                  do j = lo(2), hi(2)
                     u(hi(1)+1,j,n) = merge(
     $                    bcval(hi(1)+1,j,n)*coef(-1),
     $                    u(hi(1)+1,j,n),
     $                    mask(hi(1)+1,j) .gt. 0)
                  enddo
               else
                  do j = lo(2), hi(2)
                     u(hi(1)+1,j,n) = merge(
     $                    0.0d0,
     $                    u(hi(1)+1,j,n),
     $                    mask(hi(1)+1,j) .gt. 0)
                  enddo
               endif
               do m = 0, lenx
                  do j = lo(2), hi(2)
                     u(hi(1)+1,j,n) = merge(
     $                    u(hi(1)+1,j,n)
     $                    + u(hi(1)-m,j,n)*coef(m),
     $                    u(hi(1)+1,j,n),
     $                    mask(hi(1)+1,j) .gt. 0)
                  enddo
               enddo
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     den(hi(1),j,n)   = merge(coef(0), 0.0d0,
     $                    mask(hi(1)+1,j) .gt. 0)
                  enddo
               endif
c     :::::: now do the tangential derivative part
c     :::::: Null
c     :::::: interior part of right side
               edgloc = 0.
               i = hi(1)
               do j=lo(2)+1,hi(2)-1
                  if( mask(1 + i,-1 + j).eq.0.and.mask(1 + i,1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-U(1 + i,-1 + j,n) + U(1 + i,1 + j,n))/(2.*hy)
                  elseif( mask(1 + i,j).eq.0.and.mask(1 + i,1 + j).eq.0.and.mask(1 + i,2 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-3*U(1 + i,j,n) + 4*U(1 + i,1 + j,n) - U(1 + i,2 + j,n))/(2.*hy)
                  elseif( mask(1 + i,-2 + j).eq.0.and.mask(1 + i,-1 + j).eq.0.and.mask(1 + i,j).eq.0 ) then
                     outloc = -0.5
                     outder = (U(1 + i,-2 + j,n) - 4*U(1 + i,-1 + j,n) + 3*U(1 + i,j,n))/(2.*hy)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1) then
                        outder = exttd(i+1,j,n)
                     else
                        outder = 0.
                     endif
                  endif
                  innloc = 0.5
                  innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)

                  lambda = (edgloc-outloc)/(innloc-outloc)
                  tander(i+1,j,n) = lambda*innder+(1-lambda)*outder
                  if (flagden .eq. 1) then
                     dent(i,j,n) = 0.d0
                  endif
               enddo
c :::::: now end points
c :::::: lower right side
               j = lo(2)
               if( mask(1 + i,-1 + j).eq.0.and.mask(1 + i,1 + j).eq.0 ) then
                  outloc = -0.5
                  outder = (-U(1 + i,-1 + j,n) + U(1 + i,1 + j,n))/(2.*hy)
               elseif( mask(1 + i,j).eq.0.and.mask(1 + i,1 + j).eq.0.and.mask(1 + i,2 + j).eq.0 ) then
                  outloc = -0.5
                  outder = (-3*U(1 + i,j,n) + 4*U(1 + i,1 + j,n) - U(1 + i,2 + j,n))/(2.*hy)
               else
                  outloc = x(-1)
                  if( flagbc .eq. 1) then
                     outder = exttd(i+1,j,n)
                  else
                     outder = 0.
                  endif
               endif
               if( maskm(i,-1 + j).gt.0 ) then
                  innder = (-3*U(i,j,n) + 4*U(i,1 + j,n) - U(i,2 + j,n))/(2.*hy)
                  innden = -3.d0/(2.*hy)
               else
                  innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
                  innden = 0.d0
               endif
               innloc = 0.5
               lambda = (edgloc-outloc)/(innloc-outloc)
               tander(i+1,j,n) = lambda*innder+(1-lambda)*outder
               if (flagden .eq. 1) then
                  dent(i,j,n) = lambda*innden
               endif
c     :::::: upper right side
               j = hi(2)
               if( mask(1 + i,-1 + j).eq.0.and.mask(1 + i,1 + j).eq.0 ) then
                  outloc = -0.5
                  outder = (-U(1 + i,-1 + j,n) + U(1 + i,1 + j,n))/(2.*hy)
               elseif( mask(1 + i,-2 + j).eq.0.and.mask(1 + i,-1 + j).eq.0.and.mask(1 + i,j).eq.0 ) then
                  outloc = -0.5
                  outder = (U(1 + i,-2 + j,n) - 4*U(1 + i,-1 + j,n) + 3*U(1 + i,j,n))/(2.*hy)
               else
                  outloc = x(-1)
                  if( flagbc .eq. 1) then
                     outder = exttd(i+1,j,n)
                  else
                     outder = 0.
                  endif
               endif
               if( maskp(i,1 + j).gt.0 ) then
                  innder = (U(i,-2 + j,n) - 4*U(i,-1 + j,n) + 3*U(i,j,n))/(2.*hy)
                  innden = 3.d0/(2.*hy)
               else
                  innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
                  innden = 0.d0
               endif
               innloc = 0.5
               lambda = (edgloc-outloc)/(innloc-outloc)
               tander(i+1,j,n) = lambda*innder+(1-lambda)*outder
               if (flagden .eq. 1) then
                  dent(i,j,n) = lambda*innden
               endif
               
            else if ( bct(n) .eq. LO_REFLECT_ODD ) then
               
               do j = lo(2), hi(2)
                  u(hi(1)+1, j, n) = merge(
     $                 -u(hi(1),j,n),
     $                 u(hi(1)+1, j, n),
     $                 mask(hi(1)+1,j) .gt. 0)
c                 ::: tangential der is trivial
                  tander(hi(1)+1,j,n) = 0.0
               enddo
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     den(hi(1),j,n) = merge(-1.0d0, 0.0d0,
     $                    mask(hi(1)+1,j) .gt. 0)
                     dent(hi(1),j,n) = 0.d0
                  enddo
               endif
               
            else
               print *,'UNKNOWN BC ON RIGHT FACE IN APPLYBC'
               stop
            endif
         enddo
      endif
c
c     The Bottom of the Grid
c
      if(cdir .eq. 1) then
         do n = 1, nc
            if(is_neumann(bct(n))) then
               do i = lo(1),hi(1)
                  u(i,lo(2)-1,n) = merge(
     $                 u(i,lo(2),n),
     $                 u(i,lo(2)-1,n),
     $                 mask(i,lo(2)-1) .gt. 0)
               enddo
               if ( flagden .eq. 1 ) then
                  do i = lo(1),hi(1)
                     den(i,lo(2),n)   = 1.0
                  enddo
               endif
c     :::::: now do the tangential derivative part
c     :::::: interior part of bottom side
               j = lo(2)
               do i=lo(1)+1,hi(1)-1
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
                  tander(i,j-1,n) = innder
                  if (flagden .eq. 1) then
                     dent(i,j,n) = 0.d0
                  endif
               enddo
c :::::: now end points
c :::::: left bottom side
               i = lo(1)
               if( maskm(-1 + i,j).gt.0 ) then
                  innder = (-3*U(i,j,n) + 4*U(1 + i,j,n) - U(2 + i,j,n))/(2.*hx)
                  innden = -3.d0/(2.*hx)
               else
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
                  innden = 0.d0
               endif
               tander(i,j-1,n) = innder
               if (flagden .eq. 1) then
                  dent(i,j,n) = innden
               endif
c :::::: right bottom side
               i = hi(1)
               if( maskp(1 + i,j).gt.0 ) then
                  innder = (U(-2 + i,j,n) - 4*U(-1 + i,j,n) + 3*U(i,j,n))/(2.*hx)
                  innden = 3.d0/(2.*hx)
               else
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
                  innden = 0.d0
               endif
               tander(i,j-1,n) = innder
               if (flagden .eq. 1) then
                  dent(i,j,n) = innden
               endif
            else if (is_dirichlet(bct(n))) then
               do m=0,leny
                  x(m) = m + 0.5
               enddo
               x(-1) = - bcl/h(2)
               xInt = - 0.5
               call polyInterpCoeff(xInt, x, leny+2, coef)
               if ( flagbc .eq. 1 ) then
                  do i = lo(1), hi(1)
                     u(i,lo(2)-1,n) = merge(
     $                    bcval(i,lo(2)-1,n)*coef(-1),
     $                    u(i,lo(2)-1,n),
     $                    mask(i,lo(2)-1) .gt. 0)
                  enddo
               else
                  do i = lo(1), hi(1)
                     u(i,lo(2)-1,n) = merge(
     $                    0.0d0,
     $                    u(i,lo(2)-1,n),
     $                    mask(i,lo(2)-1) .gt. 0)
                  enddo
               endif
               do m = 0, leny
                  do i = lo(1), hi(1)
                     u(i, lo(2)-1, n) = merge(
     $                    u(i, lo(2)-1,n)
     $                    + u(i, lo(2)+m,n)*coef(m),
     $                    u(i, lo(2)-1, n),
     $                    mask(i, lo(2)-1) .gt. 0)
                  enddo
               enddo
               if ( flagden .eq. 1 ) then
                  do i = lo(1), hi(1)
                     den(i, lo(2),n)   = merge(coef(0), 0.0d0,
     $                    mask(i, lo(2)-1) .gt. 0)
                  enddo
               endif
c     :::::: now do the tangential derivative part
c     :::::: interior part of bottom side
               edgloc = 0.
               j = lo(2)
               do i=lo(1)+1,hi(1)-1
                  if(  mask(-1 + i,-1 + j).eq.0.and.mask(1 + i,-1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-U(-1 + i,-1 + j,n) + U(1 + i,-1 + j,n))/(2.*hx)
                  else if(  mask(i,-1 + j).eq.0.and.mask(1 + i,-1 + j).eq.0.and.mask(2 + i,-1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-3*U(i,-1 + j,n) + 4*U(1 + i,-1 + j,n) - U(2 + i,-1 + j,n))/(2.*hx)
                  else if(  mask(-2 + i,-1 + j).eq.0.and.mask(-1 + i,-1 + j).eq.0.and.mask(i,-1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (U(-2 + i,-1 + j,n) - 4*U(-1 + i,-1 + j,n) + 3*U(i,-1 + j,n))/(2.*hx)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1) then
                        outder = exttd(i,j-1,n)
                     else
                        outder = 0.
                     endif
                  endif
                  innloc = 0.5
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)

                  lambda = (edgloc-outloc)/(innloc-outloc)
                  tander(i,j-1,n) = lambda*innder+(1-lambda)*outder
                  if (flagden .eq. 1) then
                     dent(i,j,n) = 0.d0
                  endif
               enddo
c :::::: now end points
c :::::: left bottom side
               i = lo(1)
               if(  mask(-1 + i,-1 + j).eq.0.and.mask(1 + i,-1 + j).eq.0 ) then
                  outloc = -0.5
                  outder = (-U(-1 + i,-1 + j,n) + U(1 + i,-1 + j,n))/(2.*hx)
               else if( mask(i,-1 + j).eq.0.and.mask(1 + i,-1 + j).eq.0.and.mask(2 + i,-1 + j).eq.0 ) then
                  outloc = -0.5
                  outder = (-3*U(i,-1 + j,n) + 4*U(1 + i,-1 + j,n) - U(2 + i,-1 + j,n))/(2.*hx)
               else
                  outloc = x(-1)
                  if( flagbc .eq. 1) then
                     outder = exttd(i,j-1,n)
                  else
                     outder = 0.
                  endif
               endif
               if( maskm(-1 + i,j).gt.0 ) then
                  innder = (-3*U(i,j,n) + 4*U(1 + i,j,n) - U(2 + i,j,n))/(2.*hx)
                  innden = -3.d0/(2.*hx)
               else
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
                  innder = 0.d0
               endif
               innloc = 0.5
               lambda = (edgloc-outloc)/(innloc-outloc)
               tander(i,j-1,n) = lambda*innder+(1-lambda)*outder
               if (flagden .eq. 1) then
                  dent(i,j,n) = lambda*innden
               endif
c :::::: right bottom side
               i = hi(1)
               if( mask(-1 + i,-1 + j).eq.0.and.mask(1 + i,-1 + j).eq.0 ) then
                  outloc = -0.5
                  outder = (-U(-1 + i,-1 + j,n) + U(1 + i,-1 + j,n))/(2.*hx)
               else if( mask(-2 + i,-1 + j).eq.0.and.mask(-1 + i,-1 + j).eq.0.and.mask(i,-1 + j).eq.0 ) then
                  outloc = -0.5
                  outder = (U(-2 + i,-1 + j,n) - 4*U(-1 + i,-1 + j,n) + 3*U(i,-1 + j,n))/(2.*hx)
               else
                  outloc = x(-1)
                  if( flagbc .eq. 1) then
                     outder = exttd(i,j-1,n)
                  else
                     outder = 0.
                  endif
               endif
               if( maskp(1 + i,j).gt.0 ) then
                  innder = (U(-2 + i,j,n) - 4*U(-1 + i,j,n) + 3*U(i,j,n))/(2.*hx)
                  innden = 3.d0/(2.*hx)
               else
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
                  innden = 0.d0
               endif
               innloc = 0.5
               lambda = (edgloc-outloc)/(innloc-outloc)
               tander(i,j-1,n) = lambda*innder+(1-lambda)*outder
               if (flagden .eq. 1) then
                  dent(i,j,n) = lambda*innden
               endif
               
            else if ( bct(n) .eq. LO_REFLECT_ODD ) then
               
               do i = lo(1), hi(1)
                  u(i,lo(2)-1,n) = merge(
     $                 -u(i,lo(2),n),
     $                 u(i,lo(2)-1,n),
     $                 mask(i,lo(2)-1) .gt. 0)
                  tander(i,lo(2)-1,n) = 0.0
               enddo
               if ( flagden .eq. 1 ) then
                  do i = lo(1), hi(1)
                     den(i,lo(2),n) = merge(-1.0d0, 0.0d0,
     $                    mask(i,lo(2)-1) .gt. 0)
                     dent(i,lo(2),n) = 0.d0
                  enddo
               endif
               
            else
               print *,'UNKNOWN BC ON BOTTOM FACE IN APPLYBC'
               stop
            endif
        enddo
      endif
c     
c     The top of the grid
c
      if (cdir .eq. 3) then
         do n = 1, nc
            if(is_neumann(bct(n))) then
               do i = lo(1), hi(1)
                  u(i,hi(2)+1,n) = merge(
     $                 u(i,hi(2),n),
     $                 u(i,hi(2)+1,n),
     $                 mask(i,hi(2)+1) .gt. 0)
               enddo
               if ( flagden .eq. 1 ) then
                  do i = lo(1), hi(1)
                     den(i,hi(2),n)   = 1.0
                  enddo
               endif
c     :::::: now do the tangential derivative part
c     :::::: interior part of top side
               j = hi(2)
               do i=lo(1)+1,hi(1)-1
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
                  tander(i,j+1,n) = innder
                  if (flagden .eq. 1) then
                     dent(i,j,n) = 0.d0
                  endif
               enddo
c     :::::: now end points
c     :::::: left top side
               i = lo(1)
               if( maskm(-1 + i,j).gt.0 ) then
                  innder = (-3*U(i,j,n) + 4*U(1 + i,j,n) - U(2 + i,j,n))/(2.*hx)
                  innden = -3/(2.*hx)
               else
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
                  innden = 0.d0
               endif
               tander(i,j+1,n) = innder
               if (flagden .eq. 1) then
                  dent(i,j,n) = innden
               endif
c     :::::: right top side
               i = hi(1)
               if( maskp(1 + i,j).gt.0 ) then
                  innder = (U(-2 + i,j,n) - 4*U(-1 + i,j,n) + 3*U(i,j,n))/(2.*hx)
                  innden = 3.d0/(2.*hx)
               else
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
                  innden = 0.d0
               endif
               tander(i,j+1,n) = innder
               if (flagden .eq. 1) then
                  dent(i,j,n) = innden
               endif
            else if (is_dirichlet(bct(n))) then
               if ( bct(n) .eq. LO_REFLECT_ODD ) leny = 0
               do m=0,leny
                  x(m) = m + 0.5
               enddo
               x(-1) = - bcl/h(2)
               xInt = - 0.5
               call polyInterpCoeff(xInt, x, leny+2, coef)
               if ( flagbc .eq. 1 ) then
                  do i = lo(1), hi(1)
                     u(i,hi(2)+1,n) = merge(
     $                    bcval(i,hi(2)+1,n)*coef(-1),
     $                    u(i,hi(2)+1,n),
     $                    mask(i,hi(2)+1) .gt. 0)
                  enddo
               else
                  do i = lo(1), hi(1)
                     u(i,hi(2)+1,n) = merge(
     $                    0.0d0,
     $                    u(i,hi(2)+1,n),
     $                    mask(i,hi(2)+1) .gt. 0)
                  enddo
               endif
               do m = 0, leny
                  do i = lo(1), hi(1)
                     u(i, hi(2)+1,n) = merge(
     $                    u(i,hi(2)+1,n)
     $                    + u(i, hi(2)-m,n)*coef(m),
     $                    u(i,hi(2)+1,n),
     $                    mask(i,hi(2)+1) .gt. 0)
                  enddo
               enddo
               if ( flagden .eq. 1 ) then
                  do i = lo(1), hi(1)
                     den(i,hi(2),n)   = merge(coef(0), 0.0d0,
     $                    mask(i,hi(2)+1) .gt. 0)
                  enddo
               endif
c     :::::: now do the tangential derivative part
c     :::::: interior part of top side
               edgloc = 0.
               j = hi(2)
               do i=lo(1)+1,hi(1)-1
                  if(  mask(-1 + i,1 + j).eq.0.and.mask(1 + i,1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-U(-1 + i,1 + j,n) + U(1 + i,1 + j,n))/(2.*hx)
                  else if(  mask(i,1 + j).eq.0.and.mask(1 + i,1 + j).eq.0.and.mask(2 + i,1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-3*U(i,1 + j,n) + 4*U(1 + i,1 + j,n) - U(2 + i,1 + j,n))/(2.*hx)
                  else if( mask(-2 + i,1 + j).eq.0.and.mask(-1 + i,1 + j).eq.0.and.mask(i,1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (U(-2 + i,1 + j,n) - 4*U(-1 + i,1 + j,n) + 3*U(i,1 + j,n))/(2.*hx)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1) then
                        outder = exttd(i,j+1,n)
                     else
                        outder = 0.
                     endif
                  endif
                  innloc = 0.5
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
                  
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  tander(i,j+1,n) = lambda*innder+(1-lambda)*outder
                  if (flagden .eq. 1) then
                     dent(i,j,n) = 0.d0
                  endif
               enddo
c     :::::: now end points
c     :::::: left top side
               i = lo(1)
               if( mask(-1 + i,1 + j).eq.0.and.mask(1 + i,1 + j).eq.0 ) then
                  outloc = -0.5
                  outder = (-U(-1 + i,1 + j,n) + U(1 + i,1 + j,n))/(2.*hx)
               else if( mask(i,1 + j).eq.0.and.mask(1 + i,1 + j).eq.0.and.mask(2 + i,1 + j).eq.0 ) then
                  outloc = -0.5
                  outder = (-3*U(i,1 + j,n) + 4*U(1 + i,1 + j,n) - U(2 + i,1 + j,n))/(2.*hx)
               else
                  outloc = x(-1)
                  if( flagbc .eq. 1) then
                     outder = exttd(i,j+1,n)
                  else
                     outder = 0.
                  endif
               endif
               if( maskm(-1 + i,j).gt.0 ) then
                  innder = (-3*U(i,j,n) + 4*U(1 + i,j,n) - U(2 + i,j,n))/(2.*hx)
                  innden = -3/(2.*hx)
               else
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
                  innden = 0.d0
               endif
               innloc = 0.5
               lambda = (edgloc-outloc)/(innloc-outloc)
               tander(i,j+1,n) = lambda*innder+(1-lambda)*outder
               if (flagden .eq. 1) then
                  dent(i,j,n) = lambda*innden
               endif
c     :::::: right top side
               i = hi(1)
               if( mask(-1 + i,1 + j).eq.0.and.mask(1 + i,1 + j).eq.0 ) then
                  outloc = -0.5
                  outder = (-U(-1 + i,1 + j,n) + U(1 + i,1 + j,n))/(2.*hx)
               else if(  mask(-2 + i,1 + j).eq.0.and.mask(-1 + i,1 + j).eq.0.and.mask(i,1 + j).eq.0 ) then
                  outloc = -0.5
                  outder = (U(-2 + i,1 + j,n) - 4*U(-1 + i,1 + j,n) + 3*U(i,1 + j,n))/(2.*hx)
               else
                  outloc = x(-1)
                  if( flagbc .eq. 1) then
                     outder = exttd(i,j+1,n)
                  else
                     outder = 0.
                  endif
               endif
               if( maskp(1 + i,j).gt.0 ) then
                  innder = (U(-2 + i,j,n) - 4*U(-1 + i,j,n) + 3*U(i,j,n))/(2.*hx)
                  innden = 3.d0/(2.*hx)
               else
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
                  innden = 0.d0
               endif
               innloc = 0.5
               lambda = (edgloc-outloc)/(innloc-outloc)
               tander(i,j+1,n) = lambda*innder+(1-lambda)*outder
               if (flagden .eq. 1) then
                  dent(i,j,n) = lambda*innden
               endif
               
            else if ( bct(n) .eq. LO_REFLECT_ODD ) then
               
               do i = lo(1), hi(1)
                  u(i,hi(2)+1,n) = merge(
     $                 -u(i,hi(2),n),
     $                 u(i,hi(2)+1,n),
     $                 mask(i,hi(2)+1) .gt. 0)
                  tander(i,hi(2)+1,n) = 0.0
               enddo
               if ( flagden .eq. 1 ) then
                  do i = lo(1), hi(1)
                     den(i,hi(2),n) = merge(-1.0d0, 0.0d0,
     $                    mask(i,hi(2)+1) .gt. 0)
                     dent(i,hi(2),n) = 0.d0
                  enddo
               endif
               
            else
               print *,'UNKNOWN BC ON TOP FACE IN APPLYBC'
               stop
            endif
         enddo
      endif
c
      end
