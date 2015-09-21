subroutine adv_update(Cold,Cnew,ng_C,Sold,Snew,ng_S,aofs,ng_a,phi,ng_phi,dt,num_comp,lo,hi)

  use bl_types
  implicit none

  integer        , intent(in   ) :: ng_C,ng_S,ng_a,ng_phi,num_comp
  integer        , intent(in   ) :: lo(2), hi(2)
  real(kind=dp_t), intent(inout) :: Cold(lo(1)-ng_C:hi(1)+ng_C,lo(2)-ng_C:hi(2)+ng_C,0:num_comp-1)
  real(kind=dp_t), intent(inout) :: Cnew(lo(1)-ng_C:hi(1)+ng_C,lo(2)-ng_C:hi(2)+ng_C,0:num_comp-1)
  real(kind=dp_t), intent(inout) :: Sold(lo(1)-ng_S:hi(1)+ng_S,lo(2)-ng_S:hi(2)+ng_S)
  real(kind=dp_t), intent(inout) :: Snew(lo(1)-ng_S:hi(1)+ng_S,lo(2)-ng_S:hi(2)+ng_S)
  real(kind=dp_t), intent(inout) :: aofs(lo(1)-ng_a:hi(1)+ng_a,lo(2)-ng_a:hi(2)+ng_a,0:num_comp-1)
  real(kind=dp_t), intent(inout) :: phi(lo(1)-ng_phi:hi(1)+ng_phi,lo(2)-ng_phi:hi(2)+ng_phi)
  real(kind=dp_t), intent(inout) :: dt

  integer :: i,j,comp

  do comp=0,num_comp-1
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)           
           Cnew(i,j,comp) = (Cold(i,j,comp)*Sold(i,j) - dt*aofs(i,j,comp)/phi(i,j))/Snew(i,j)
        end do
     end do
  end do ! loop over components

end subroutine adv_update


subroutine adv_fdiv(aofs,ng_a,edgex,edgey,ng_se,umac,vmac,ng_u,ax,ay,vol,num_comp,lo,hi)

  use bl_types
  implicit none

  integer        , intent(in   ) :: ng_a,ng_se,ng_u,num_comp
  integer        , intent(in   ) :: lo(2), hi(2)
  real(kind=dp_t), intent(inout) :: aofs(lo(1)-ng_a:hi(1)+ng_a,lo(2)-ng_a:hi(2)+ng_a,0:num_comp-1)
  real(kind=dp_t), intent(inout) :: edgex(lo(1)-ng_se:hi(1)+ng_se+1,lo(2)-ng_se:hi(2)+ng_se,0:num_comp-1)
  real(kind=dp_t), intent(inout) :: edgey(lo(1)-ng_se:hi(1)+ng_se,lo(2)-ng_se:hi(2)+ng_se+1,0:num_comp-1)
  real(kind=dp_t), intent(in   ) :: umac(lo(1)-ng_u:hi(1)+ng_u+1,lo(2)-ng_u:hi(2)+ng_u)
  real(kind=dp_t), intent(in   ) :: vmac(lo(1)-ng_u:hi(1)+ng_u,lo(2)-ng_u:hi(2)+ng_u+1)
  real(kind=dp_t), intent(in   ) :: ax(lo(1):hi(1)+1,lo(2):hi(2))
  real(kind=dp_t), intent(in   ) :: ay(lo(1):hi(1),lo(2):hi(2)+1)
  real(kind=dp_t), intent(in   ) :: vol(lo(1):hi(1),lo(2):hi(2))

  integer :: i,j,comp

  do comp=0,num_comp-1

     do j = lo(2),hi(2)
        do i = lo(1),hi(1)+1
           edgex(i,j,comp) = edgex(i,j,comp) * umac(i,j) * ax(i,j)
        end do
     end do

     do j = lo(2),hi(2)+1
        do i = lo(1),hi(1)
           edgey(i,j,comp) = edgey(i,j,comp) * vmac(i,j) * ay(i,j)
        end do
     end do

     do j = lo(2),hi(2)
        do i = lo(1),hi(1)
           aofs(i,j,comp) = (edgex(i+1,j,comp) - edgex(i,j,comp) &
                +            edgey(i,j+1,comp) - edgey(i,j,comp) ) / vol(i,j)
        end do
     end do

  end do ! loop over components

end subroutine adv_fdiv


subroutine bds_est_eigen(s,ng_s,umac,vmac,ng_u,phi,ng_phi,lo,hi,eigenval)

  use bl_types
  implicit none

  integer        , intent(in   ) :: ng_s,ng_u,ng_phi
  integer        , intent(in   ) :: lo(2), hi(2)
  real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_s:hi(1)+ng_s,lo(2)-ng_s:hi(2)+ng_s)
  real(kind=dp_t), intent(in   ) :: umac(lo(1)-ng_u:hi(1)+ng_u+1,lo(2)-ng_u:hi(2)+ng_u)
  real(kind=dp_t), intent(in   ) :: vmac(lo(1)-ng_u:hi(1)+ng_u,lo(2)-ng_u:hi(2)+ng_u+1)
  real(kind=dp_t), intent(in   ) :: phi(lo(1)-ng_phi:hi(1)+ng_phi,lo(2)-ng_phi:hi(2)+ng_phi)
  real(kind=dp_t), intent(inout) :: eigenval(2)

  integer         :: i,j,comp
  real(kind=dp_t) :: esat, ephi, eigtmp

  eigenval = 0.d0
  do j = lo(2),hi(2)
     do i = lo(1),hi(1)+1
        esat = 0.5d0*(   s(i,j) +   s(i-1,j) )
        ephi = 0.5d0*( phi(i,j) + phi(i-1,j) )
        eigtmp = umac(i,j)/(esat * ephi)
        eigenval(1) = max(eigenval(1),dabs(eigtmp))
     end do
  end do
  do j = lo(2),hi(2)+1
     do i = lo(1),hi(1)
        esat = 0.5d0*(   s(i,j) +   s(i,j-1) )
        ephi = 0.5d0*( phi(i,j) + phi(i,j-1) )
        eigtmp = vmac(i,j)/(esat * ephi)
        eigenval(2) = max(eigenval(2),dabs(eigtmp))
     end do
  end do
  
end subroutine bds_est_eigen



subroutine bds_set_ngrowhyp(bds_ngrowhyp)
  implicit none
  integer, intent(out) :: bds_ngrowhyp
  bds_ngrowhyp = 3
end subroutine bds_set_ngrowhyp


subroutine bds_set_ngrowforce(bds_ngrowforce)
  implicit none
  integer, intent(out) :: bds_ngrowforce
  bds_ngrowforce = 1
end subroutine bds_set_ngrowforce


subroutine bds_set_nwork(bds_nwork)
  implicit none
  integer, intent(out) :: bds_nwork
  bds_nwork = 3;
end subroutine bds_set_nwork


subroutine bds_edge_states(s,ng_s,ci,ng_ci,sedgex,sedgey,ng_se,umac,vmac,ng_u,force,ng_f,&
     slope,ng_c,nw,dx,dt,num_comp,is_cons,lo,hi,bc)

  use bl_types
  implicit none

  real(kind=dp_t), intent(in   ) :: dx(2),dt
  integer        , intent(in   ) :: num_comp,is_cons
  integer        , intent(in   ) :: lo(2), hi(2)
  real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_s: hi(1)+ng_s, lo(2)-ng_s: hi(2)+ng_s ,0:num_comp-1)
  real(kind=dp_t), intent(in   ) :: ci(lo(1)-ng_ci: hi(1)+ng_ci, lo(2)-ng_ci: hi(2)+ng_ci)
  real(kind=dp_t), intent(inout) :: sedgex(lo(1)-ng_se:hi(1)+ng_se+1,lo(2)-ng_se:hi(2)+ng_se,0:num_comp-1)
  real(kind=dp_t), intent(inout) :: sedgey(lo(1)-ng_se:hi(1)+ng_se,lo(2)-ng_se:hi(2)+ng_se+1,0:num_comp-1)
  real(kind=dp_t), intent(inout) :: umac(lo(1)-ng_u:hi(1)+ng_u+1,lo(2)-ng_u:hi(2)+ng_u)
  real(kind=dp_t), intent(inout) :: vmac(lo(1)-ng_u:hi(1)+ng_u,lo(2)-ng_u:hi(2)+ng_u+1)
  real(kind=dp_t), intent(in   ) :: force(lo(1)-ng_f:hi(1)+ng_f,lo(2)-ng_f: hi(2)+ng_f,0:num_comp-1)
  real(kind=dp_t), intent(inout) :: slope(lo(1)-ng_c: hi(1)+ng_c, lo(2)-ng_c: hi(2)+ng_c ,1:nw)
  integer        , intent(in   ) :: bc(2,2) ! bc(dir,lohi)

  integer :: dm,ng_s,ng_c,ng_u,ng_se,ng_f,nw,ng_ci
  integer :: i,comp,ioff,joff,j
  logical :: is_conservative

  interface
     subroutine bdsslope_2d(lo,hi,s,ng_s,slope,ng_c,nw,dx,bc)
       use bl_types
       implicit none
       integer        , intent(in   ) :: lo(2),hi(2),ng_s,ng_c,nw
       real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_s: hi(1)+ng_s, lo(2)-ng_s: hi(2)+ng_s)
       real(kind=dp_t), intent(inout) :: slope(lo(1)-ng_c: hi(1)+ng_c, lo(2)-ng_c: hi(2)+ng_c ,1:nw)
       real(kind=dp_t), intent(in   ) :: dx(2)
       integer        , intent(in   ) :: bc(2,2) ! bc(dir,lohi)
     end subroutine bdsslope_2d

     subroutine bdsconc_2d(lo,hi,s,ng_s,ci,ng_ci,slope,ng_c,umac,vmac,ng_u,force,ng_f, &
          sedgex,sedgey,ng_se,nw,dx,dt,bc,is_conservative)
       use bl_types
       implicit none
       integer        , intent(in   ) :: lo(2),hi(2),ng_s,ng_c,ng_u,ng_f,ng_se,nw,ng_ci
       real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_s: hi(1)+ng_s, lo(2)-ng_s: hi(2)+ng_s)
       real(kind=dp_t), intent(in   ) :: ci(lo(1)-ng_ci: hi(1)+ng_ci, lo(2)-ng_ci: hi(2)+ng_ci)
       real(kind=dp_t), intent(in   ) :: slope(lo(1)-ng_c: hi(1)+ng_c, lo(2)-ng_c: hi(2)+ng_c ,1:nw)
       real(kind=dp_t), intent(in   ) :: umac(lo(1)-ng_u:hi(1)+ng_u+1,lo(2)-ng_u:hi(2)+ng_u)
       real(kind=dp_t), intent(in   ) :: vmac(lo(1)-ng_u:hi(1)+ng_u,lo(2)-ng_u:hi(2)+ng_u+1)
       real(kind=dp_t), intent(in   ) :: force(lo(1)-ng_f:hi(1)+ng_f,lo(2)-ng_f: hi(2)+ng_f)
       real(kind=dp_t), intent(inout) :: sedgex(lo(1)-ng_se:hi(1)+ng_se+1,lo(2)-ng_se:hi(2)+ng_se)
       real(kind=dp_t), intent(inout) :: sedgey(lo(1)-ng_se:hi(1)+ng_se,lo(2)-ng_se:hi(2)+ng_se+1)
       real(kind=dp_t), intent(in   ) :: dx(2),dt
       integer        , intent(in   ) :: bc(2,2) ! bc(dir,lohi)
       logical        , intent(in   ) :: is_conservative
     end subroutine bdsconc_2d
  end interface

  is_conservative = is_cons > 0

  do j = lo(2),hi(2)
     do i = lo(1),hi(1)+1
        if (umac(i,j) .gt. 0) then
           ioff = -1
        else
           ioff = 0
        endif
        umac(i,j) = umac(i,j) * ci(i+ioff,j)
     end do
  end do

  do j = lo(2),hi(2)+1
     do i = lo(1),hi(1)
        if (vmac(i,j) .gt. 0) then
           joff = -1
        else
           joff = 0
        endif
        vmac(i,j) = vmac(i,j) * ci(i,j+joff)
     end do
  end do

  do comp=0,num_comp-1

     call bdsslope_2d(lo, hi, s(:,:,comp), ng_s, slope, ng_c, nw, dx, bc) 

     call bdsconc_2d(lo, hi, s(:,:,comp), ng_s, ci, ng_ci,&
          slope(:,:,:), ng_c, &
          umac(:,:), vmac(:,:), ng_u, &
          force(:,:,comp), ng_f, &
          sedgex(:,:,comp), sedgey(:,:,comp), ng_se, nw, &
          dx, dt, bc, is_conservative)

  end do ! loop over components

  do j = lo(2),hi(2)
     do i = lo(1),hi(1)+1
        if (umac(i,j) .gt. 0) then
           ioff = -1
        else
           ioff = 0
        endif
        umac(i,j) = umac(i,j) / ci(i+ioff,j)
     end do
  end do

  do j = lo(2),hi(2)+1
     do i = lo(1),hi(1)
        if (vmac(i,j) .gt. 0) then
           joff = -1
        else
           joff = 0
        endif
        vmac(i,j) = vmac(i,j) / ci(i,j+joff)
     end do
  end do

end subroutine bds_edge_states


subroutine bdsslope_2d(lo,hi,s,ng_s,slope,ng_c,nw,dx,bc)

  use bl_types
  implicit none

  integer        , intent(in   ) :: lo(2),hi(2),ng_s,ng_c,nw
  real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_s: hi(1)+ng_s, lo(2)-ng_s: hi(2)+ng_s)
  real(kind=dp_t), intent(inout) :: slope(lo(1)-ng_c: hi(1)+ng_c, lo(2)-ng_c: hi(2)+ng_c ,1:nw)
  real(kind=dp_t), intent(in   ) :: dx(2)
  integer        , intent(in   ) :: bc(2,2) ! bc(dir,lohi)

  ! local variables
  real(kind=dp_t), allocatable :: sint(:,:)

  real(kind=dp_t) :: diff(4)
  real(kind=dp_t) :: smin(4)
  real(kind=dp_t) :: smax(4)
  real(kind=dp_t) :: sc(4)

  real(kind=dp_t) :: hx,hy,eps
  real(kind=dp_t) :: sumloc,redfac,redmax,div,kdp,sumdif,sgndif
  integer         :: i,j,ll,mm

  real(kind=dp_t) :: t1,t2,teps

  if (nw < 3 .or. ng_s < 1) then
     call bl_abort('Not enough work space for bdsslope_2d')
  end if

  ! nodal with one ghost cell
  allocate(sint(lo(1)-ng_c:hi(1)+ng_c+1,lo(2)-ng_c:hi(2)+ng_c+1))

  hx = dx(1)
  hy = dx(2)

  eps = 1.d-10

  ! bicubic interpolation to corner points
  ! (i,j) refers to lower corner of cell
  do j = lo(2)-ng_c,hi(2)+ng_c+1
     do i = lo(1)-ng_c,hi(1)+ng_c+1
        sint(i,j) = (s(i-2,j-2) + s(i-2,j+1) + s(i+1,j-2) + s(i+1,j+1) &
             - 7.d0*(s(i-2,j-1) + s(i-2,j  ) + s(i-1,j-2) + s(i  ,j-2) + & 
             s(i-1,j+1) + s(i  ,j+1) + s(i+1,j-1) + s(i+1,j  )) &
             + 49.d0*(s(i-1,j-1) + s(i  ,j-1) + s(i-1,j  ) + s(i  ,j  )) ) / 144.d0
     enddo
  enddo

  ! NOTE: By convention, bc(dir,lohi) values mean the following:
  !
  !   -1    Reflect odd
  !    0    Periodic
  !    1    Reflect even
  !    2    First-order extrap
  !    3    Dirichlet
  !    4    High-order extrap
  !
  if (bc(1,1) .eq. 3) then
     do j=lo(2)-ng_c,hi(2)+ng_c
        do i=lo(1)-ng_c,lo(1)
           sint(i,j) = (s(lo(1)-1,j)+s(lo(1)-1,j-1)) * 0.5d0
        end do
        i = lo(1) + 1
        sint(i,j) = (s(i,j)+s(i-1,j)+s(i-1,j-1)+s(i,j-1)) * 0.25d0
     end do
  else if (bc(1,1) .eq. -1) then
     i = lo(1)
     do j=lo(2)-ng_c,hi(2)+ng_c
        sint(i,j) = 0.d0
     end do
     do j=lo(2),hi(2)
        do i=1,ng_c
           sint(lo(1)-i,j) = sint(lo(1)+i,j)
        end do
     end do
  end if

  if (bc(1,2) .eq. 3) then
     do j=lo(2)-ng_c,hi(2)+ng_c
        do i=hi(1)+2,hi(1)+ng_c+1
           sint(i,j) = (s(hi(1)+1,j) + s(hi(1)+1,j-1)) * 0.5d0
        end do
        i = hi(1)
        sint(i,j) = (s(i,j)+s(i+1,j)+s(i+1,j+1)+s(i,j+1)) * 0.25d0
     end do
  else if (bc(1,2) .eq. -1) then
     i = hi(1)+1
     do j=lo(2)-ng_c,hi(2)+ng_c
        sint(i,j) = 0.d0
     end do
     do j=lo(2),hi(2)
        do i=1,ng_c
           sint(hi(1)+1+i,j) = sint(hi(1)+1-i,j)
        end do
     end do
  end if

  if (bc(2,1) .eq. 3) then
     do i=lo(1)-ng_c,hi(1)+ng_c
        do j=lo(2)-ng_c,lo(2)
           sint(i,j) = (s(i,lo(2)-1) + s(i-1,lo(2)-1)) * 0.5d0
        end do
        j = lo(2) + 1
        sint(i,j) = (s(i,j)+s(i-1,j)+s(i-1,j-1)+s(i,j-1)) * 0.25d0
     end do
  else if (bc(2,1) .eq. -1) then
     j = lo(2)
     do i=lo(1)-ng_c,hi(1)+ng_c
        sint(i,j) = 0.d0
     end do
     do i=lo(1),hi(1)
        do j=1,ng_c
           sint(i,lo(2)-j) = sint(i,lo(2)+j)
        end do
     end do
  end if

  if (bc(2,2) .eq. 3) then
     do i=lo(1)-ng_c,hi(1)+ng_c
        do j=hi(2)+2,hi(2)+ng_c+1
           sint(i,j) = (s(i,hi(2)+1) + s(i-1,hi(2)+1)) * 0.5d0
        end do
        j = hi(2)
        sint(i,j) = (s(i,j)+s(i+1,j)+s(i+1,j+1)+s(i,j+1)) * 0.25d0
     end do
  else if (bc(2,2) .eq. -1) then
     j = hi(2)+1
     do i=lo(1)-ng_c,hi(1)+ng_c
        sint(i,j) = 0.d0
     end do
     do i=lo(1),hi(1)
        do j=1,ng_c
           sint(i,hi(2)+1+j) = sint(i,hi(2)+1-j)
        end do
     end do
  end if

  do j = lo(2)-ng_c,hi(2)+ng_c
     do i = lo(1)-ng_c,hi(1)+ng_c

        ! compute initial estimates of slopes from unlimited corner points

        ! sx
        slope(i,j,1) = 0.5d0*(sint(i+1,j+1) + sint(i+1,j  ) - &
             sint(i  ,j+1) - sint(i  ,j  ) ) / hx

        ! sy
        slope(i,j,2) = 0.5d0*(sint(i+1,j+1) - sint(i+1,j  ) + &
             sint(i  ,j+1) - sint(i  ,j  ) ) / hy

        ! sxy
        slope(i,j,3) = ( sint(i+1,j+1) - sint(i+1,j  ) &
             -sint(i  ,j+1) + sint(i  ,j  ) ) / (hx*hy)

        ! ++ / sint(i+1,j+1)
        sc(4) = s(i,j) + 0.5d0*(hx*slope(i,j,1) + hy*slope(i,j,2))  &
             + 0.25d0*hx*hy*slope(i,j,3)

        ! +- / sint(i+1,j  )
        sc(3) = s(i,j) + 0.5d0*(hx*slope(i,j,1) - hy*slope(i,j,2))  &
             - 0.25d0*hx*hy*slope(i,j,3)

        ! -+ / sint(i  ,j+1)
        sc(2) = s(i,j) - 0.5d0*(hx*slope(i,j,1) - hy*slope(i,j,2)) &
             - 0.25d0*hx*hy*slope(i,j,3)

        ! -- / sint(i  ,j  )
        sc(1) = s(i,j) - 0.5d0*(hx*slope(i,j,1) + hy*slope(i,j,2)) &
             + 0.25d0*hx*hy*slope(i,j,3)

        ! enforce max/min bounds
        smin(4) = min(s(i,j), s(i+1,j), s(i,j+1), s(i+1,j+1))
        smax(4) = max(s(i,j), s(i+1,j), s(i,j+1), s(i+1,j+1))

        smin(3) = min(s(i,j), s(i+1,j), s(i,j-1), s(i+1,j-1))
        smax(3) = max(s(i,j), s(i+1,j), s(i,j-1), s(i+1,j-1))

        smin(2) = min(s(i,j), s(i-1,j), s(i,j+1), s(i-1,j+1))
        smax(2) = max(s(i,j), s(i-1,j), s(i,j+1), s(i-1,j+1))

        smin(1) = min(s(i,j), s(i-1,j), s(i,j-1), s(i-1,j-1))
        smax(1) = max(s(i,j), s(i-1,j), s(i,j-1), s(i-1,j-1))

        do mm=1,4
           sc(mm) = max(min(sc(mm), smax(mm)), smin(mm))
        enddo

        ! iterative loop
        do ll = 1,3
           sumloc = 0.25d0*(sc(4) + sc(3) + sc(2) + sc(1))
           sumdif = (sumloc - s(i,j))*4.d0
           sgndif = sign(1.d0,sumdif)

           do mm=1,4
              diff(mm) = (sc(mm) - s(i,j))*sgndif
           enddo

           kdp = 0

           do mm=1,4
              if (diff(mm) .gt. eps) then
                 kdp = kdp+1
              end if
           end do

           do mm = 1,4 
              if (kdp.lt.1) then 
                 div = 1.d0
              else
                 div = dble(kdp)
              end if

              if (diff(mm).gt.eps) then
                 redfac = sumdif*sgndif/div
                 kdp = kdp-1
              else
                 redfac = 0.d0
              end if

              if (sgndif .gt. 0.d0) then
                 redmax = sc(mm) - smin(mm)
              else
                 redmax = smax(mm) - sc(mm)
              end if

              redfac = min(redfac,redmax)
              sumdif = sumdif - redfac*sgndif
              sc(mm) = sc(mm) - redfac*sgndif
           enddo
        enddo

        ! ! check that results are in bounds
        ! teps = 1.e-18
        ! smin(4) = min(s(i,j), s(i+1,j), s(i,j+1), s(i+1,j+1))
        ! smax(4) = max(s(i,j), s(i+1,j), s(i,j+1), s(i+1,j+1))

        ! smin(3) = min(s(i,j), s(i+1,j), s(i,j-1), s(i+1,j-1))
        ! smax(3) = max(s(i,j), s(i+1,j), s(i,j-1), s(i+1,j-1))

        ! smin(2) = min(s(i,j), s(i-1,j), s(i,j+1), s(i-1,j+1))
        ! smax(2) = max(s(i,j), s(i-1,j), s(i,j+1), s(i-1,j+1))

        ! smin(1) = min(s(i,j), s(i-1,j), s(i,j-1), s(i-1,j-1))
        ! smax(1) = max(s(i,j), s(i-1,j), s(i,j-1), s(i-1,j-1))
        ! do ll=1,4
        !    if ( sc(ll)-smax(ll) .gt. teps  .or.  smin(ll)-sc(ll) .gt. teps ) then
        !       print *,'sc OOB',i,j,ll,sc(ll),smin(ll),smax(ll)
        !       stop
        !    endif
        ! enddo

        ! final slopes

        ! sx
        slope(i,j,1) = 0.5d0*( sc(4) + sc(3) &
             -sc(1) - sc(2))/hx

        ! sy
        slope(i,j,2) = 0.5d0*( sc(4) + sc(2) &
             -sc(1) - sc(3))/hy

        ! sxy
        slope(i,j,3) = ( sc(1) + sc(4) &
             -sc(2) - sc(3) ) / (hx*hy)

     enddo
  enddo

  deallocate(sint)

end subroutine bdsslope_2d

subroutine bdsconc_2d(lo,hi,s,ng_s,ci,ng_ci,slope,ng_c,umac,vmac,ng_u,force,ng_f, &
     sedgex,sedgey,ng_se,nw,dx,dt,bc,is_conservative)

  use bl_types
  implicit none

  integer        , intent(in   ) :: lo(2),hi(2),ng_s,ng_c,ng_u,ng_f,ng_se,nw,ng_ci
  real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_s: hi(1)+ng_s, lo(2)-ng_s: hi(2)+ng_s)
  real(kind=dp_t), intent(in   ) :: ci(lo(1)-ng_ci: hi(1)+ng_ci, lo(2)-ng_ci: hi(2)+ng_ci)
  real(kind=dp_t), intent(in   ) :: slope(lo(1)-ng_c: hi(1)+ng_c, lo(2)-ng_c: hi(2)+ng_c ,1:nw)
  real(kind=dp_t), intent(in   ) :: umac(lo(1)-ng_u:hi(1)+ng_u+1,lo(2)-ng_u:hi(2)+ng_u)
  real(kind=dp_t), intent(in   ) :: vmac(lo(1)-ng_u:hi(1)+ng_u,lo(2)-ng_u:hi(2)+ng_u+1)
  real(kind=dp_t), intent(in   ) :: force(lo(1)-ng_f:hi(1)+ng_f,lo(2)-ng_f: hi(2)+ng_f)
  real(kind=dp_t), intent(inout) :: sedgex(lo(1)-ng_se:hi(1)+ng_se+1,lo(2)-ng_se:hi(2)+ng_se)
  real(kind=dp_t), intent(inout) :: sedgey(lo(1)-ng_se:hi(1)+ng_se,lo(2)-ng_se:hi(2)+ng_se+1)
  real(kind=dp_t), intent(in   ) :: dx(2),dt
  integer        , intent(in   ) :: bc(2,2) ! bc(dir,lohi)
  logical        , intent(in   ) :: is_conservative

  ! local variables
  integer i,j,ioff,joff,ll

  real(kind=dp_t), allocatable ::   ux(:,:)
  real(kind=dp_t), allocatable ::   vy(:,:)
  real(kind=dp_t), allocatable :: divu(:,:)

  real(kind=dp_t) :: isign,jsign,hx,hy
  real(kind=dp_t) :: del(2),p1(2),p2(2),p3(2)
  real(kind=dp_t) :: val1,val2,val3
  real(kind=dp_t) :: u,v,gamma
  real(kind=dp_t) :: dt2,dt3,half

  interface
     subroutine eval_2d(s,slope,del,val)
       use bl_types
       implicit none
       real(kind=dp_t), intent(in   ) :: s
       real(kind=dp_t), intent(in   ) :: slope(:)
       real(kind=dp_t), intent(in   ) :: del(:)
       real(kind=dp_t), intent(  out) :: val
     end subroutine eval_2d
  end interface

  allocate(  ux(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))
  allocate(  vy(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))
  allocate(divu(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))

  hx = dx(1)
  hy = dx(2)

  dt2 = dt/2.d0
  dt3 = dt/3.d0

  half = 0.5d0

  ! compute cell-centered ux and vy
  do j=lo(2)-1,hi(2)+1
     do i=lo(1)-1,hi(1)+1
        ux(i,j) = (umac(i+1,j) - umac(i,j)) / hx
        vy(i,j) = (vmac(i,j+1) - vmac(i,j)) / hy
        divu(i,j) = ux(i,j) + vy(i,j)
     end do
  end do

  ! compute sedgex on x-faces
  do j=lo(2),hi(2)
     do i=lo(1),hi(1)+1

        ! compute sedgex without transverse corrections
        if (umac(i,j) .gt. 0) then
           isign = 1.d0
           ioff = -1
        else
           isign = -1.d0
           ioff = 0
        endif

        ! centroid of rectangular volume
        del(1) = isign*0.5d0*hx - 0.5d0*umac(i,j)*dt
        del(2) = 0.d0

        call eval_2d(s(i+ioff,j),slope(i+ioff,j,:),del,sedgex(i,j))

        ! source term
        if (is_conservative) then
           sedgex(i,j) = sedgex(i,j)*(1.d0 - dt2*ux(i+ioff,j)) + dt2*force(i+ioff,j)
        else
           sedgex(i,j) = sedgex(i,j)*(1.d0 + dt2*vy(i+ioff,j)) + dt2*force(i+ioff,j)
        end if

        ! compute \Gamma^{y+}
        if (vmac(i+ioff,j+1) .gt. 0) then
           jsign = 1.d0
           joff = 0
        else
           jsign = -1.d0
           joff = 1
        endif

        u = 0.d0
        if (umac(i,j)*umac(i,j+joff) .gt. 0) then
           u = umac(i,j+joff)
        endif

        p1(1) = isign*0.5d0*hx
        p1(2) = jsign*0.5d0*hy

        p2(1) = isign*0.5d0*hx - umac(i,j)*dt
        p2(2) = jsign*0.5d0*hy

        p3(1) = isign*0.5d0*hx - u*dt
        p3(2) = jsign*0.5d0*hy - vmac(i+ioff,j+1)*dt

        do ll=1,2
           del(ll) = (p2(ll)+p3(ll))/2.d0
        end do
        call eval_2d(s(i+ioff,j+joff),slope(i+ioff,j+joff,:),del,val1)

        do ll=1,2
           del(ll) = (p1(ll)+p3(ll))/2.d0
        end do
        call eval_2d(s(i+ioff,j+joff),slope(i+ioff,j+joff,:),del,val2)

        do ll=1,2
           del(ll) = (p1(ll)+p2(ll))/2.d0
        end do
        call eval_2d(s(i+ioff,j+joff),slope(i+ioff,j+joff,:),del,val3)

        ! average these centroid values to get the average value
        gamma = (val1+val2+val3)/3.d0

        ! source term
        if (is_conservative) then
           gamma = gamma*(1.d0 - dt3*divu(i+ioff,j+joff))
        end if

        ! correct sedgex with \Gamma^{y+}
        gamma = gamma * vmac(i+ioff,j+1)
        sedgex(i,j) = sedgex(i,j) - dt*gamma/(2.d0*hy)

        ! compute \Gamma^{y-}
        if (vmac(i+ioff,j) .gt. 0) then
           jsign = 1.d0
           joff = -1
        else
           jsign = -1.d0
           joff = 0
        endif

        u = 0.d0
        if (umac(i,j)*umac(i,j+joff) .gt. 0) then
           u = umac(i,j+joff)
        endif

        p1(1) = isign*0.5d0*hx
        p1(2) = jsign*0.5d0*hy

        p2(1) = isign*0.5d0*hx - umac(i,j)*dt
        p2(2) = jsign*0.5d0*hy

        p3(1) = isign*0.5d0*hx - u*dt
        p3(2) = jsign*0.5d0*hy - vmac(i+ioff,j)*dt

        do ll=1,2
           del(ll) = (p2(ll)+p3(ll))/2.d0
        end do
        call eval_2d(s(i+ioff,j+joff),slope(i+ioff,j+joff,:),del,val1)

        do ll=1,2
           del(ll) = (p1(ll)+p3(ll))/2.d0
        end do
        call eval_2d(s(i+ioff,j+joff),slope(i+ioff,j+joff,:),del,val2)

        do ll=1,2
           del(ll) = (p1(ll)+p2(ll))/2.d0
        end do
        call eval_2d(s(i+ioff,j+joff),slope(i+ioff,j+joff,:),del,val3)

        ! average these centroid values to get the average value
        gamma = (val1+val2+val3)/3.d0

        ! source term
        if (is_conservative) then
           gamma = gamma*(1.d0 - dt3*divu(i+ioff,j+joff))
        end if

        ! correct sedgex with \Gamma^{y-}
        gamma = gamma * vmac(i+ioff,j)
        sedgex(i,j) = sedgex(i,j) + dt*gamma/(2.d0*hy)

     end do
  end do

  ! compute sedgey on y-faces    
  do j=lo(2),hi(2)+1
     do i=lo(1),hi(1)

        ! compute sedgey without transverse corrections

        ! centroid of rectangular volume
        if (vmac(i,j) .gt. 0) then
           jsign = 1.d0
           joff = -1
        else
           jsign = -1.d0
           joff = 0
        endif

        del(1) = 0.d0
        del(2) = jsign*0.5d0*hy - 0.5d0*vmac(i,j)*dt
        call eval_2d(s(i,j+joff),slope(i,j+joff,:),del,sedgey(i,j))

        ! source term
        if (is_conservative) then
           sedgey(i,j) = sedgey(i,j)*(1.d0 - dt2*vy(i,j+joff)) + dt2*force(i,j+joff)
        else
           sedgey(i,j) = sedgey(i,j)*(1.d0 + dt2*ux(i,j+joff)) + dt2*force(i,j+joff)
        end if

        ! compute \Gamma^{x+} without corner corrections
        if (umac(i+1,j+joff) .gt. 0) then
           isign = 1.d0
           ioff = 0
        else
           isign = -1.d0
           ioff = 1
        endif

        v = 0.d0
        if (vmac(i,j)*vmac(i+ioff,j) .gt. 0) then
           v = vmac(i+ioff,j)
        endif

        p1(1) = isign*0.5d0*hx
        p1(2) = jsign*0.5d0*hy

        p2(1) = isign*0.5d0*hx
        p2(2) = jsign*0.5d0*hy - vmac(i,j)*dt

        p3(1) = isign*0.5d0*hx - umac(i+1,j+joff)*dt
        p3(2) = jsign*0.5d0*hy - v*dt

        do ll=1,2
           del(ll) = (p2(ll)+p3(ll))/2.d0
        end do
        call eval_2d(s(i+ioff,j+joff),slope(i+ioff,j+joff,:),del,val1)

        do ll=1,2
           del(ll) = (p1(ll)+p3(ll))/2.d0
        end do
        call eval_2d(s(i+ioff,j+joff),slope(i+ioff,j+joff,:),del,val2)

        do ll=1,2
           del(ll) = (p1(ll)+p2(ll))/2.d0
        end do
        call eval_2d(s(i+ioff,j+joff),slope(i+ioff,j+joff,:),del,val3)

        ! average these centroid values to get the average value
        gamma = (val1+val2+val3)/3.d0

        ! source term
        if (is_conservative) then
           gamma = gamma*(1.d0 - dt3*divu(i+ioff,j+joff))
        end if

        ! correct sedgey with \Gamma^{x+}
        gamma = gamma * umac(i+1,j+joff)
        sedgey(i,j) = sedgey(i,j) - dt*gamma/(2.d0*hx)

        ! compute \Gamma^{x-}
        if (umac(i,j+joff) .gt. 0) then
           isign = 1.d0
           ioff = -1
        else
           isign = -1.d0
           ioff = 0
        endif

        v = 0.d0
        if (vmac(i,j)*vmac(i+ioff,j) .gt. 0) then
           v = vmac(i+ioff,j)
        endif

        p1(1) = isign*0.5d0*hx
        p1(2) = jsign*0.5d0*hy

        p2(1) = isign*0.5d0*hx
        p2(2) = jsign*0.5d0*hy - vmac(i,j)*dt

        p3(1) = isign*0.5d0*hx - umac(i,j+joff)*dt
        p3(2) = jsign*0.5d0*hy - v*dt

        do ll=1,2
           del(ll) = (p2(ll)+p3(ll))/2.d0
        end do
        call eval_2d(s(i+ioff,j+joff),slope(i+ioff,j+joff,:),del,val1)

        do ll=1,2
           del(ll) = (p1(ll)+p3(ll))/2.d0
        end do
        call eval_2d(s(i+ioff,j+joff),slope(i+ioff,j+joff,:),del,val2)

        do ll=1,2
           del(ll) = (p1(ll)+p2(ll))/2.d0
        end do
        call eval_2d(s(i+ioff,j+joff),slope(i+ioff,j+joff,:),del,val3)

        ! average these centroid values to get the average value
        gamma = (val1+val2+val3)/3.d0

        ! source term
        if (is_conservative) then
           gamma = gamma*(1.d0 - dt3*divu(i+ioff,j+joff))
        end if

        ! correct sedgey with \Gamma^{x-}
        gamma = gamma * umac(i,j+joff)
        sedgey(i,j) = sedgey(i,j) + dt*gamma/(2.d0*hx)

     end do
  end do

  deallocate(ux,vy,divu)

  ! NOTE: By convention, bc(dir,lohi) values mean the following:
  !
  !   -1    Reflect odd
  !    0    Periodic
  !    1    Reflect even
  !    2    First-order extrap
  !    3    Dirichlet
  !    4    High-order extrap
  !
  i = lo(1)
  if (bc(1,1) .eq. 3) then
     do j=lo(2),hi(2)
        sedgex(i,j) = s(i-1,j)
     end do
  else if (bc(1,1) .eq. -1) then
     do j=lo(2),hi(2)
        sedgex(i,j) = 0.d0
     end do
  end if

  i = hi(1)+1
  if (bc(1,2) .eq. 3) then
     do j=lo(2),hi(2)
        sedgex(i,j) = s(i,j)
     end do
  else if (bc(1,2) .eq. -1) then
     do j=lo(2),hi(2)
        sedgex(i,j) = 0.d0
     end do
  end if

  j = lo(2)
  if (bc(2,1) .eq. 3) then
     do i=lo(1),hi(1)
        sedgey(i,j) = s(i,j-1)
     end do
  else if (bc(2,1) .eq. -1) then
     do i=lo(1),hi(1)
        sedgey(i,j) = 0.d0
     end do
  end if

  j = hi(2)+1
  if (bc(2,2) .eq. 3) then
     do i=lo(1),hi(1)
        sedgey(i,j) = s(i,j)
     end do
  else if (bc(2,2) .eq. -1) then
     do i=lo(1),hi(1)
        sedgey(i,j) = 0.d0
     end do
  end if

end subroutine bdsconc_2d

subroutine eval_2d(s,slope,del,val)

  use bl_types
  implicit none

  real(kind=dp_t), intent(in   ) :: s
  real(kind=dp_t), intent(in   ) :: slope(:)
  real(kind=dp_t), intent(in   ) :: del(:)
  real(kind=dp_t), intent(  out) :: val

  val = s + del(1)*slope(1) + del(2)*slope(2) + del(1)*del(2)*slope(3)

end subroutine eval_2d
