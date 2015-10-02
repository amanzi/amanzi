subroutine adv_update(Cold,Cnew,ng_C,Sold,Snew,ng_S,aofs,ng_a,phi,ng_phi,dt,num_comp,lo,hi)

  use bl_types
  implicit none

  integer        , intent(in   ) :: ng_C,ng_S,ng_a,ng_phi,num_comp
  integer        , intent(in   ) :: lo(3), hi(3)
  real(kind=dp_t), intent(inout) :: Cold(lo(1)-ng_C:hi(1)+ng_C,lo(2)-ng_C:hi(2)+ng_C,lo(3)-ng_C:hi(3)+ng_C,0:num_comp-1)
  real(kind=dp_t), intent(inout) :: Cnew(lo(1)-ng_C:hi(1)+ng_C,lo(2)-ng_C:hi(2)+ng_C,lo(3)-ng_C:hi(3)+ng_C,0:num_comp-1)
  real(kind=dp_t), intent(inout) :: Sold(lo(1)-ng_S:hi(1)+ng_S,lo(2)-ng_S:hi(2)+ng_S,lo(3)-ng_S:hi(3)+ng_S)
  real(kind=dp_t), intent(inout) :: Snew(lo(1)-ng_S:hi(1)+ng_S,lo(2)-ng_S:hi(2)+ng_S,lo(3)-ng_S:hi(3)+ng_S)
  real(kind=dp_t), intent(inout) :: aofs(lo(1)-ng_a:hi(1)+ng_a,lo(2)-ng_a:hi(2)+ng_a,lo(3)-ng_a:hi(3)+ng_a,0:num_comp-1)
  real(kind=dp_t), intent(inout) :: phi(lo(1)-ng_phi:hi(1)+ng_phi,lo(2)-ng_phi:hi(2)+ng_phi,lo(3)-ng_phi:hi(3)+ng_phi)
  real(kind=dp_t), intent(inout) :: dt

  integer :: i,j,k,comp

  do comp=0,num_comp-1
     do k = lo(3),hi(3)
        do j = lo(2),hi(2)
           do i = lo(1),hi(1)           
              Cnew(i,j,k,comp) = (Cold(i,j,k,comp)*Sold(i,j,k) - dt*aofs(i,j,k,comp)/phi(i,j,k))/Snew(i,j,k)
           end do
        end do
     end do
  end do ! loop over components

end subroutine adv_update


subroutine adv_fdiv(aofs,ng_a,edgex,edgey,edgez,ng_se,umac,vmac,wmac,ng_u,ax,ay,az,vol,num_comp,lo,hi)

  use bl_types
  implicit none

  integer        , intent(in   ) :: ng_a,ng_se,ng_u,num_comp
  integer        , intent(in   ) :: lo(3), hi(3)
  real(kind=dp_t), intent(inout) :: aofs(lo(1)-ng_a:hi(1)+ng_a,lo(2)-ng_a:hi(2)+ng_a,lo(3)-ng_a:hi(3)+ng_a,0:num_comp-1)
  real(kind=dp_t), intent(inout) :: edgex(lo(1)-ng_se:hi(1)+ng_se+1,lo(2)-ng_se:hi(2)+ng_se,lo(3)-ng_se:hi(3)+ng_se,0:num_comp-1)
  real(kind=dp_t), intent(inout) :: edgey(lo(1)-ng_se:hi(1)+ng_se,lo(2)-ng_se:hi(2)+ng_se+1,lo(3)-ng_se:hi(3)+ng_se,0:num_comp-1)
  real(kind=dp_t), intent(inout) :: edgez(lo(1)-ng_se:hi(1)+ng_se,lo(2)-ng_se:hi(2)+ng_se,lo(3)-ng_se:hi(3)+ng_se+1,0:num_comp-1)
  real(kind=dp_t), intent(in   ) :: umac(lo(1)-ng_u:hi(1)+ng_u+1,lo(2)-ng_u:hi(2)+ng_u,lo(3)-ng_u:hi(3)+ng_u)
  real(kind=dp_t), intent(in   ) :: vmac(lo(1)-ng_u:hi(1)+ng_u,lo(2)-ng_u:hi(2)+ng_u+1,lo(3)-ng_u:hi(3)+ng_u)
  real(kind=dp_t), intent(in   ) :: wmac(lo(1)-ng_u:hi(1)+ng_u,lo(2)-ng_u:hi(2)+ng_u,lo(3)-ng_u:hi(3)+ng_u+1)
  real(kind=dp_t), intent(in   ) :: ax(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3))
  real(kind=dp_t), intent(in   ) :: ay(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3))
  real(kind=dp_t), intent(in   ) :: az(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1)
  real(kind=dp_t), intent(in   ) :: vol(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

  integer :: i,j,k,comp

  do comp=0,num_comp-1

     do k = lo(3),hi(3)
        do j = lo(2),hi(2)
           do i = lo(1),hi(1)+1
              edgex(i,j,k,comp) = edgex(i,j,k,comp) * umac(i,j,k) * ax(i,j,k)
           end do
        end do
     end do

     do k = lo(3),hi(3)
        do j = lo(2),hi(2)+1
           do i = lo(1),hi(1)
              edgey(i,j,k,comp) = edgey(i,j,k,comp) * vmac(i,j,k) * ay(i,j,k)
           end do
        end do
     end do

     do k = lo(3),hi(3)+1
        do j = lo(2),hi(2)
           do i = lo(1),hi(1)
              edgez(i,j,k,comp) = edgez(i,j,k,comp) * wmac(i,j,k) * az(i,j,k)
           end do
        end do
     end do

     do k = lo(3),hi(3)
        do j = lo(2),hi(2)
           do i = lo(1),hi(1)
              aofs(i,j,k,comp) = (edgex(i+1,j,k,comp) - edgex(i,j,k,comp) &
                   +              edgey(i,j+1,k,comp) - edgey(i,j,k,comp) &
                   +              edgez(i,j,k+1,comp) - edgez(i,j,k,comp) ) / vol(i,j,k)
           end do
        end do
     end do

  end do ! loop over components

end subroutine adv_fdiv


subroutine bds_est_eigen(s,ng_s,umac,vmac,wmac,ng_u,phi,ng_phi,lo,hi,eigenval)

  use bl_types
  implicit none

  integer        , intent(in   ) :: ng_s,ng_u,ng_phi
  integer        , intent(in   ) :: lo(3), hi(3)
  real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_s:hi(1)+ng_s,lo(2)-ng_s:hi(2)+ng_s,lo(3)-ng_s:hi(3)+ng_s)
  real(kind=dp_t), intent(in   ) :: umac(lo(1)-ng_u:hi(1)+ng_u+1,lo(2)-ng_u:hi(2)+ng_u,lo(3)-ng_u:hi(3)+ng_u)
  real(kind=dp_t), intent(in   ) :: vmac(lo(1)-ng_u:hi(1)+ng_u,lo(2)-ng_u:hi(2)+ng_u+1,lo(3)-ng_u:hi(3)+ng_u)
  real(kind=dp_t), intent(in   ) :: wmac(lo(1)-ng_u:hi(1)+ng_u,lo(2)-ng_u:hi(2)+ng_u,lo(3)-ng_u:hi(3)+ng_u+1)
  real(kind=dp_t), intent(in   ) :: phi(lo(1)-ng_phi:hi(1)+ng_phi,lo(2)-ng_phi:hi(2)+ng_phi,lo(3)-ng_phi:hi(3)+ng_phi)
  real(kind=dp_t), intent(inout) :: eigenval(3)

  integer         :: i,j,k,comp
  real(kind=dp_t) :: esat, ephi, eigtmp

  eigenval = 0.d0
  do k = lo(3),hi(3)
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)+1
           esat = 0.5d0*(   s(i,j,k) +   s(i-1,j,k) )
           ephi = 0.5d0*( phi(i,j,k) + phi(i-1,j,k) )
           eigtmp = umac(i,j,k)/(esat * ephi)
           eigenval(1) = max(eigenval(1),dabs(eigtmp))
        end do
     end do
  end do
  do k = lo(3),hi(3)
     do j = lo(2),hi(2)+1
        do i = lo(1),hi(1)
           esat = 0.5d0*(   s(i,j,k) +   s(i,j-1,k) )
           ephi = 0.5d0*( phi(i,j,k) + phi(i,j-1,k) )
           eigtmp = vmac(i,j,k)/(esat * ephi)
           eigenval(2) = max(eigenval(2),dabs(eigtmp))
        end do
     end do
  end do
  do k = lo(3),hi(3)+1
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)
           esat = 0.5d0*(   s(i,j,k) +   s(i,j,k-1) )
           ephi = 0.5d0*( phi(i,j,k) + phi(i,j,k-1) )
           eigtmp = wmac(i,j,k)/(esat * ephi)
           eigenval(3) = max(eigenval(3),dabs(eigtmp))
        end do
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
  bds_nwork = 7; ! Derives x, y, z, xy, xz, yz, xyz
end subroutine bds_set_nwork


subroutine bds_edge_states(s,ng_s,ci,ng_ci,sedgex,sedgey,sedgez,ng_se,umac,vmac,wmac,ng_u,force,ng_f,&
     slope,ng_c,nw,dx,dt,num_comp,is_cons,lo,hi,bc)

  use bl_types
  implicit none

  real(kind=dp_t), intent(in   ) :: dx(3),dt
  integer        , intent(in   ) :: num_comp,is_cons
  integer        , intent(in   ) :: lo(3), hi(3)
  real(kind=dp_t), intent(in   ) :: s(lo(1)-ng_s: hi(1)+ng_s, lo(2)-ng_s: hi(2)+ng_s , lo(3)-ng_s: hi(3)+ng_s ,0:num_comp-1)
  real(kind=dp_t), intent(in   ) :: ci(lo(1)-ng_ci: hi(1)+ng_ci, lo(2)-ng_ci: hi(2)+ng_ci, lo(3)-ng_ci: hi(3)+ng_ci)
  real(kind=dp_t), intent(inout) :: sedgex(lo(1)-ng_se:hi(1)+ng_se+1,lo(2)-ng_se:hi(2)+ng_se,lo(3)-ng_se:hi(3)+ng_se,0:num_comp-1)
  real(kind=dp_t), intent(inout) :: sedgey(lo(1)-ng_se:hi(1)+ng_se,lo(2)-ng_se:hi(2)+ng_se+1,lo(3)-ng_se:hi(3)+ng_se,0:num_comp-1)
  real(kind=dp_t), intent(inout) :: sedgez(lo(1)-ng_se:hi(1)+ng_se,lo(2)-ng_se:hi(2)+ng_se,lo(3)-ng_se:hi(3)+ng_se+1,0:num_comp-1)
  real(kind=dp_t), intent(inout) :: umac(lo(1)-ng_u:hi(1)+ng_u+1,lo(2)-ng_u:hi(2)+ng_u,lo(3)-ng_u:hi(3)+ng_u)
  real(kind=dp_t), intent(inout) :: vmac(lo(1)-ng_u:hi(1)+ng_u,lo(2)-ng_u:hi(2)+ng_u+1,lo(3)-ng_u:hi(3)+ng_u)
  real(kind=dp_t), intent(inout) :: wmac(lo(1)-ng_u:hi(1)+ng_u,lo(2)-ng_u:hi(2)+ng_u,lo(3)-ng_u:hi(3)+ng_u+1)
  real(kind=dp_t), intent(in   ) :: force(lo(1)-ng_f:hi(1)+ng_f,lo(2)-ng_f:hi(2)+ng_f,lo(3)-ng_f:hi(3)+ng_f,0:num_comp-1)
  real(kind=dp_t), intent(inout) :: slope(lo(1)-ng_c:hi(1)+ng_c,lo(2)-ng_c:hi(2)+ng_c,lo(3)-ng_c:hi(3)+ng_c,1:nw)
  integer        , intent(in   ) :: bc(3,2) ! bc(dir,lohi)

  integer :: dm,ng_s,ng_c,ng_u,ng_se,ng_f,nw,ng_ci
  integer :: comp,i,ioff,j,joff,k,koff
  logical :: is_conservative

  interface
     subroutine bdsslope_3d(lo,hi,s,ng_s,slope,ng_c,nw,dx,bc)
       use bl_types
       implicit none
       integer        , intent(in   ) :: lo(3),hi(3),ng_s,ng_c,nw
       real(kind=dp_t), intent(in   ) ::     s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
       real(kind=dp_t), intent(inout) :: slope(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:,1:)
       real(kind=dp_t), intent(in   ) :: dx(:)
       integer        , intent(in   ) :: bc(:,:) ! bc(dir,lohi)
     end subroutine bdsslope_3d

     subroutine bdsconc_3d(lo,hi,s,ng_s,ci,ng_ci,slope,ng_c,umac,vmac,wmac,ng_u,force,ng_f, &
          sedgex,sedgey,sedgez,ng_se,nw,dx,dt,bc,is_conservative)       
       use bl_types
       implicit none
       integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_c,ng_u,ng_f,ng_se,nw,ng_ci
       real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :)
       real(kind=dp_t), intent(in   ) ::     ci(lo(1)-ng_ci:,lo(2)-ng_ci:,lo(3)-ng_ci:)
       real(kind=dp_t), intent(in   ) ::  slope(lo(1)-ng_c :,lo(2)-ng_c :,lo(3)-ng_c :,1:)
       real(kind=dp_t), intent(in   ) ::   umac(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :)
       real(kind=dp_t), intent(in   ) ::   vmac(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :)
       real(kind=dp_t), intent(in   ) ::   wmac(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :)
       real(kind=dp_t), intent(in   ) ::  force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :)
       real(kind=dp_t), intent(inout) :: sedgex(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:)
       real(kind=dp_t), intent(inout) :: sedgey(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:)
       real(kind=dp_t), intent(inout) :: sedgez(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:)
       real(kind=dp_t), intent(in   ) :: dx(3),dt
       integer        , intent(in   ) :: bc(3,2) ! bc(dir,lohi)
       logical        , intent(in   ) :: is_conservative
     end subroutine bdsconc_3d
  end interface

  is_conservative = is_cons > 0

  do k = lo(3),hi(3)
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)+1
           if (umac(i,j,k) .gt. 0) then
              ioff = -1
           else
              ioff = 0
           endif
           umac(i,j,k) = umac(i,j,k) * ci(i+ioff,j,k)
        end do
     end do
  end do

  do k = lo(3),hi(3)
     do j = lo(2),hi(2)+1
        do i = lo(1),hi(1)
           if (vmac(i,j,k) .gt. 0) then
              joff = -1
           else
              joff = 0
           endif
           vmac(i,j,k) = vmac(i,j,k) * ci(i,j+joff,k)
        end do
     end do
  end do

  do k = lo(3),hi(3)+1
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)
           if (wmac(i,j,k) .gt. 0) then
              koff = -1
           else
              koff = 0
           endif
           wmac(i,j,k) = wmac(i,j,k) * ci(i,j,k+koff)
        end do
     end do
  end do

  do comp=0,num_comp-1

     call bdsslope_3d(lo, hi, s(:,:,:,comp), ng_s, slope, ng_c, nw, dx, bc) 

     call bdsconc_3d(lo, hi, s(:,:,:,comp), ng_s, ci, ng_ci,&
          slope(:,:,:,:), ng_c, &
          umac(:,:,:), vmac(:,:,:), wmac(:,:,:), ng_u, &
          force(:,:,:,comp), ng_f, &
          sedgex(:,:,:,comp), sedgey(:,:,:,comp), sedgez(:,:,:,comp), ng_se, nw, &
          dx, dt, bc, is_conservative)

  end do ! loop over components

  do k = lo(3),hi(3)
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)+1
           if (umac(i,j,k) .gt. 0) then
              ioff = -1
           else
              ioff = 0
           endif
           umac(i,j,k) = umac(i,j,k) / ci(i+ioff,j,k)
        end do
     end do
  end do

  do k = lo(3),hi(3)
     do j = lo(2),hi(2)+1
        do i = lo(1),hi(1)
           if (vmac(i,j,k) .gt. 0) then
              joff = -1
           else
              joff = 0
           endif
           vmac(i,j,k) = vmac(i,j,k) / ci(i,j+joff,k)
        end do
     end do
  end do

  do k = lo(3),hi(3)+1
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)
           if (wmac(i,j,k) .gt. 0) then
              koff = -1
           else
              koff = 0
           endif
           wmac(i,j,k) = wmac(i,j,k) / ci(i,j,k+koff)
        end do
     end do
  end do

end subroutine bds_edge_states


subroutine bdsslope_3d(lo,hi,s,ng_s,slope,ng_c,nw,dx,bc)

  use bl_types
  implicit none

  integer        , intent(in   ) :: lo(3),hi(3),ng_s,ng_c,nw
  real(kind=dp_t), intent(in   ) ::     s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
  real(kind=dp_t), intent(inout) :: slope(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:,1:)
  real(kind=dp_t), intent(in   ) :: dx(:)
  integer        , intent(in   ) :: bc(:,:) ! bc(dir,lohi)

  ! local variables
  real(kind=dp_t), allocatable :: sint(:,:,:)

  real(kind=dp_t) :: diff(8)
  real(kind=dp_t) :: smin(8)
  real(kind=dp_t) :: smax(8)
  real(kind=dp_t) :: sc(8)

  real(kind=dp_t) :: c1,c2,c3,c4
  real(kind=dp_t) :: hx,hy,hz,eps
  real(kind=dp_t) :: sumloc,redfac,redmax,div,kdp,sumdif,sgndif
  integer         :: i,j,k,ll,mm

  interface
     subroutine eval_3d(s,slope,del,val)
       use bl_types
       implicit none
       real(kind=dp_t), intent(in   ) :: s
       real(kind=dp_t), intent(in   ) :: slope(:)
       real(kind=dp_t), intent(in   ) :: del(:)
       real(kind=dp_t), intent(  out) :: val
     end subroutine eval_3d
  end interface

  ! nodal with one ghost cell
  allocate(sint(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+2,lo(3)-1:hi(3)+2))

  hx = dx(1)
  hy = dx(2)
  hz = dx(3)

  eps = 1.d-10

  c1 = (343.d0/1728.d0)
  c2 = (49.d0 /1728.d0)
  c3 = (7.d0  /1728.d0)
  c4 = (1.d0  /1728.d0)

  ! tricubic interpolation to corner points
  ! (i,j,k) refers to lower corner of cell
  do k = lo(3)-1,hi(3)+2
     do j = lo(2)-1,hi(2)+2
        do i = lo(1)-1,hi(1)+2
           sint(i,j,k) = c1*( s(i  ,j  ,k  ) + s(i-1,j  ,k  ) + s(i  ,j-1,k  ) &
                +s(i  ,j  ,k-1) + s(i-1,j-1,k  ) + s(i-1,j  ,k-1) &
                +s(i  ,j-1,k-1) + s(i-1,j-1,k-1) ) &
                -c2*( s(i-1,j  ,k+1) + s(i  ,j  ,k+1) + s(i-1,j-1,k+1) &
                +s(i  ,j-1,k+1) + s(i-1,j+1,k  ) + s(i  ,j+1,k  ) &
                +s(i-2,j  ,k  ) + s(i+1,j  ,k  ) + s(i-2,j-1,k  ) &
                +s(i+1,j-1,k  ) + s(i-1,j-2,k  ) + s(i  ,j-2,k  ) &
                +s(i-1,j+1,k-1) + s(i  ,j+1,k-1) + s(i-2,j  ,k-1) &
                +s(i+1,j  ,k-1) + s(i-2,j-1,k-1) + s(i+1,j-1,k-1) &
                +s(i-1,j-2,k-1) + s(i  ,j-2,k-1) + s(i-1,j  ,k-2) &
                +s(i  ,j  ,k-2) + s(i-1,j-1,k-2) + s(i  ,j-1,k-2) ) &
                +c3*( s(i-1,j+1,k+1) + s(i  ,j+1,k+1) + s(i-2,j  ,k+1) &
                +s(i+1,j  ,k+1) + s(i-2,j-1,k+1) + s(i+1,j-1,k+1) &
                +s(i-1,j-2,k+1) + s(i  ,j-2,k+1) + s(i-2,j+1,k  ) &
                +s(i+1,j+1,k  ) + s(i-2,j-2,k  ) + s(i+1,j-2,k  ) &
                +s(i-2,j+1,k-1) + s(i+1,j+1,k-1) + s(i-2,j-2,k-1) &
                +s(i+1,j-2,k-1) + s(i-1,j+1,k-2) + s(i  ,j+1,k-2) &
                +s(i-2,j  ,k-2) + s(i+1,j  ,k-2) + s(i-2,j-1,k-2) &
                +s(i+1,j-1,k-2) + s(i-1,j-2,k-2) + s(i  ,j-2,k-2) ) &
                -c4*( s(i-2,j+1,k+1) + s(i+1,j+1,k+1) + s(i-2,j-2,k+1) &
                +s(i+1,j-2,k+1) + s(i-2,j+1,k-2) + s(i+1,j+1,k-2) &
                +s(i-2,j-2,k-2) + s(i+1,j-2,k-2) )
        enddo
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
     call bl_abort('dirichlet bc not yet implmented for 3D bds')
  else if (bc(1,1) .eq. -1) then
     call bl_abort('reflect odd bc not yet implmented for 3D bds')
  end if

  if (bc(1,2) .eq. 3) then
     call bl_abort('dirichlet bc not yet implmented for 3D bds')
  else if (bc(1,2) .eq. -1) then
     call bl_abort('reflect odd bc not yet implmented for 3D bds')
  end if

  if (bc(2,1) .eq. 3) then
     call bl_abort('dirichlet bc not yet implmented for 3D bds')
  else if (bc(2,1) .eq. -1) then
     call bl_abort('reflect odd bc not yet implmented for 3D bds')
  end if

  if (bc(2,2) .eq. 3) then
     call bl_abort('dirichlet bc not yet implmented for 3D bds')
  else if (bc(2,2) .eq. -1) then
     call bl_abort('reflect odd bc not yet implmented for 3D bds')
  end if

  if (bc(3,1) .eq. 3) then
     call bl_abort('dirichlet bc not yet implmented for 3D bds')
  else if (bc(3,1) .eq. -1) then
     call bl_abort('reflect odd bc not yet implmented for 3D bds')
  end if

  if (bc(3,2) .eq. 3) then
     call bl_abort('dirichlet bc not yet implmented for 3D bds')
  else if (bc(3,2) .eq. -1) then
     call bl_abort('reflect odd bc not yet implmented for 3D bds')
  end if


  do k = lo(3)-1,hi(3)+1
     do j = lo(2)-1,hi(2)+1
        do i = lo(1)-1,hi(1)+1 

           ! compute initial estimates of slopes from unlimited corner points

           ! sx
           slope(i,j,k,1) = 0.25d0*( ( sint(i+1,j  ,k  ) + sint(i+1,j+1,k  ) &
                +sint(i+1,j  ,k+1) + sint(i+1,j+1,k+1)) &
                -( sint(i  ,j  ,k  ) + sint(i  ,j+1,k  ) &
                +sint(i  ,j  ,k+1) + sint(i  ,j+1,k+1)) ) / hx

           ! sy
           slope(i,j,k,2) = 0.25d0*( ( sint(i  ,j+1,k  ) + sint(i+1,j+1,k  ) &
                +sint(i  ,j+1,k+1) + sint(i+1,j+1,k+1)) &
                -( sint(i  ,j  ,k  ) + sint(i+1,j  ,k  ) &
                +sint(i  ,j  ,k+1) + sint(i+1,j  ,k+1)) ) / hy

           ! sz
           slope(i,j,k,3) = 0.25d0*( ( sint(i  ,j  ,k+1) + sint(i+1,j  ,k+1) &
                +sint(i  ,j+1,k+1) + sint(i+1,j+1,k+1)) &
                -( sint(i  ,j  ,k  ) + sint(i+1,j  ,k  ) &
                +sint(i  ,j+1,k  ) + sint(i+1,j+1,k  )) ) / hz

           ! sxy
           slope(i,j,k,4) = 0.5d0*( ( sint(i  ,j  ,k  ) + sint(i  ,j  ,k+1) &
                +sint(i+1,j+1,k  ) + sint(i+1,j+1,k+1)) &
                -( sint(i+1,j  ,k  ) + sint(i+1,j  ,k+1) &
                +sint(i  ,j+1,k  ) + sint(i  ,j+1,k+1)) ) / (hx*hy)

           ! sxz
           slope(i,j,k,5) = 0.5d0*( ( sint(i  ,j  ,k  ) + sint(i  ,j+1,k  ) &
                +sint(i+1,j  ,k+1) + sint(i+1,j+1,k+1)) &
                -( sint(i+1,j  ,k  ) + sint(i+1,j+1,k  ) &
                +sint(i  ,j  ,k+1) + sint(i  ,j+1,k+1)) ) / (hx*hz)

           ! syz
           slope(i,j,k,6) = 0.5d0*( ( sint(i  ,j  ,k  ) + sint(i+1,j  ,k  ) &
                +sint(i  ,j+1,k+1) + sint(i+1,j+1,k+1)) &
                -( sint(i  ,j  ,k+1) + sint(i+1,j  ,k+1) &
                +sint(i  ,j+1,k  ) + sint(i+1,j+1,k  )) ) / (hy*hz)

           ! sxyz
           slope(i,j,k,7) = (-sint(i  ,j  ,k  ) + sint(i+1,j  ,k  ) + sint(i  ,j+1,k  ) &
                +sint(i  ,j  ,k+1) - sint(i+1,j+1,k  ) - sint(i+1,j  ,k+1) &
                -sint(i  ,j+1,k+1) + sint(i+1,j+1,k+1) ) / (hx*hy*hz)

           ! +++ / sint(i+1,j+1,k+1)
           sc(8) = s(i,j,k) &
                +0.5d0*( hx*slope(i,j,k,1)+hy*slope(i,j,k,2)+hz*slope(i,j,k,3)) &
                +0.25d0*( hx*hy*slope(i,j,k,4)+hx*hz*slope(i,j,k,5)+hy*hz*slope(i,j,k,6)) &
                +0.125d0*hx*hy*hz*slope(i,j,k,7)

           ! ++- / sint(i+1,j+1,k  )
           sc(7) = s(i,j,k) &
                +0.5d0*( hx*slope(i,j,k,1)+hy*slope(i,j,k,2)-hz*slope(i,j,k,3)) &
                +0.25d0*( hx*hy*slope(i,j,k,4)-hx*hz*slope(i,j,k,5)-hy*hz*slope(i,j,k,6)) &
                -0.125d0*hx*hy*hz*slope(i,j,k,7)

           ! +-+ / sint(i+1,j  ,k+1)
           sc(6) = s(i,j,k) &
                +0.5d0*( hx*slope(i,j,k,1)-hy*slope(i,j,k,2)+hz*slope(i,j,k,3)) &
                +0.25d0*(-hx*hy*slope(i,j,k,4)+hx*hz*slope(i,j,k,5)-hy*hz*slope(i,j,k,6)) &
                -0.125d0*hx*hy*hz*slope(i,j,k,7)

           ! +-- / sint(i+1,j  ,k  )
           sc(5) = s(i,j,k) &
                +0.5d0*( hx*slope(i,j,k,1)-hy*slope(i,j,k,2)-hz*slope(i,j,k,3)) &
                +0.25d0*(-hx*hy*slope(i,j,k,4)-hx*hz*slope(i,j,k,5)+hy*hz*slope(i,j,k,6)) &
                +0.125d0*hx*hy*hz*slope(i,j,k,7)

           ! -++ / sint(i  ,j+1,k+1)
           sc(4) = s(i,j,k) &
                +0.5d0*(-hx*slope(i,j,k,1)+hy*slope(i,j,k,2)+hz*slope(i,j,k,3)) &
                +0.25d0*(-hx*hy*slope(i,j,k,4)-hx*hz*slope(i,j,k,5)+hy*hz*slope(i,j,k,6)) &
                -0.125d0*hx*hy*hz*slope(i,j,k,7)

           ! -+- / sint(i  ,j+1,k  )
           sc(3) = s(i,j,k) &
                +0.5d0*(-hx*slope(i,j,k,1)+hy*slope(i,j,k,2)-hz*slope(i,j,k,3)) &
                +0.25d0*(-hx*hy*slope(i,j,k,4)+hx*hz*slope(i,j,k,5)-hy*hz*slope(i,j,k,6)) &
                +0.125d0*hx*hy*hz*slope(i,j,k,7)

           ! --+ / sint(i  ,j  ,k+1)
           sc(2) = s(i,j,k) &
                +0.5d0*(-hx*slope(i,j,k,1)-hy*slope(i,j,k,2)+hz*slope(i,j,k,3)) &
                +0.25d0*( hx*hy*slope(i,j,k,4)-hx*hz*slope(i,j,k,5)-hy*hz*slope(i,j,k,6)) &
                +0.125d0*hx*hy*hz*slope(i,j,k,7)

           ! ---/ sint(i  ,j  ,k  )
           sc(1) = s(i,j,k) &
                +0.5d0*(-hx*slope(i,j,k,1)-hy*slope(i,j,k,2)-hz*slope(i,j,k,3)) &
                +0.25d0*( hx*hy*slope(i,j,k,4)+hx*hz*slope(i,j,k,5)+hy*hz*slope(i,j,k,6)) &
                -0.125d0*hx*hy*hz*slope(i,j,k,7)

           ! enforce max/min bounds
           smin(8) = min(s(i,j,k),s(i+1,j,k),s(i,j+1,k),s(i,j,k+1), &
                s(i+1,j+1,k),s(i+1,j,k+1),s(i,j+1,k+1),s(i+1,j+1,k+1))
           smax(8) = max(s(i,j,k),s(i+1,j,k),s(i,j+1,k),s(i,j,k+1), &
                s(i+1,j+1,k),s(i+1,j,k+1),s(i,j+1,k+1),s(i+1,j+1,k+1))

           smin(7) = min(s(i,j,k-1),s(i+1,j,k-1),s(i,j+1,k-1),s(i,j,k), &
                s(i+1,j+1,k-1),s(i+1,j,k),s(i,j+1,k),s(i+1,j+1,k))
           smax(7) = max(s(i,j,k-1),s(i+1,j,k-1),s(i,j+1,k-1),s(i,j,k), &
                s(i+1,j+1,k-1),s(i+1,j,k),s(i,j+1,k),s(i+1,j+1,k))

           smin(6) = min(s(i,j-1,k),s(i+1,j-1,k),s(i,j,k),s(i,j-1,k+1), &
                s(i+1,j,k),s(i+1,j-1,k+1),s(i,j,k+1),s(i+1,j,k+1))
           smax(6) = max(s(i,j-1,k),s(i+1,j-1,k),s(i,j,k),s(i,j-1,k+1), &
                s(i+1,j,k),s(i+1,j-1,k+1),s(i,j,k+1),s(i+1,j,k+1))

           smin(5) = min(s(i,j-1,k-1),s(i+1,j-1,k-1),s(i,j,k-1),s(i,j-1,k), &
                s(i+1,j,k-1),s(i+1,j-1,k),s(i,j,k),s(i+1,j,k))
           smax(5) = max(s(i,j-1,k-1),s(i+1,j-1,k-1),s(i,j,k-1),s(i,j-1,k), &
                s(i+1,j,k-1),s(i+1,j-1,k),s(i,j,k),s(i+1,j,k))

           smin(4) = min(s(i-1,j,k),s(i,j,k),s(i-1,j+1,k),s(i-1,j,k+1), &
                s(i,j+1,k),s(i,j,k+1),s(i-1,j+1,k+1),s(i,j+1,k+1))
           smax(4) = max(s(i-1,j,k),s(i,j,k),s(i-1,j+1,k),s(i-1,j,k+1), &
                s(i,j+1,k),s(i,j,k+1),s(i-1,j+1,k+1),s(i,j+1,k+1))

           smin(3) = min(s(i-1,j,k-1),s(i,j,k-1),s(i-1,j+1,k-1),s(i-1,j,k), &
                s(i,j+1,k-1),s(i,j,k),s(i-1,j+1,k),s(i,j+1,k))
           smax(3) = max(s(i-1,j,k-1),s(i,j,k-1),s(i-1,j+1,k-1),s(i-1,j,k), &
                s(i,j+1,k-1),s(i,j,k),s(i-1,j+1,k),s(i,j+1,k))

           smin(2) = min(s(i-1,j-1,k),s(i,j-1,k),s(i-1,j,k),s(i-1,j-1,k+1), &
                s(i,j,k),s(i,j-1,k+1),s(i-1,j,k+1),s(i,j,k+1))
           smax(2) = max(s(i-1,j-1,k),s(i,j-1,k),s(i-1,j,k),s(i-1,j-1,k+1), &
                s(i,j,k),s(i,j-1,k+1),s(i-1,j,k+1),s(i,j,k+1))

           smin(1) = min(s(i-1,j-1,k-1),s(i,j-1,k-1),s(i-1,j,k-1),s(i-1,j-1,k), &
                s(i,j,k-1),s(i,j-1,k),s(i-1,j,k),s(i,j,k))
           smax(1) = max(s(i-1,j-1,k-1),s(i,j-1,k-1),s(i-1,j,k-1),s(i-1,j-1,k), &
                s(i,j,k-1),s(i,j-1,k),s(i-1,j,k),s(i,j,k))

           do mm=1,8
              sc(mm) = max(min(sc(mm), smax(mm)), smin(mm))
           enddo

           ! iterative loop
           do ll = 1,3 
              sumloc = 0.125d0*(sc(1)+sc(2)+sc(3)+sc(4)+sc(5)+sc(6)+sc(7)+sc(8))
              sumdif = (sumloc - s(i,j,k))*8.d0
              sgndif = sign(1.d0,sumdif)

              do mm=1,8
                 diff(mm) = (sc(mm) - s(i,j,k))*sgndif
              enddo

              kdp = 0

              do mm=1,8
                 if (diff(mm) .gt. eps) then
                    kdp = kdp+1
                 end if
              end do

              do mm = 1,8
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

           ! final slopes

           ! sx
           slope(i,j,k,1) = 0.25d0*( ( sc(5) + sc(7) &
                +sc(6) + sc(8)) &
                -( sc(1) + sc(3) &
                +sc(2) + sc(4)) ) / hx

           ! sy
           slope(i,j,k,2) = 0.25d0*( ( sc(3) + sc(7) &
                +sc(4) + sc(8)) &
                -( sc(1) + sc(5) &
                +sc(2) + sc(6)) ) / hy

           ! sz
           slope(i,j,k,3) = 0.25d0*( ( sc(2) + sc(6) &
                +sc(4) + sc(8)) &
                -( sc(1) + sc(5) &
                +sc(3) + sc(7)) ) / hz

           ! sxy
           slope(i,j,k,4) = 0.5d0*( ( sc(1) + sc(2) &
                +sc(7) + sc(8)) &
                -( sc(5) + sc(6) &
                +sc(3) + sc(4)) ) / (hx*hy)

           ! sxz
           slope(i,j,k,5) = 0.5d0*( ( sc(1) + sc(3) &
                +sc(6) + sc(8)) &
                -( sc(5) + sc(7) &
                +sc(2) + sc(4)) ) / (hx*hz)

           ! syz
           slope(i,j,k,6) = 0.5d0*( ( sc(1) + sc(5) &
                +sc(4) + sc(8)) &
                -( sc(2) + sc(6) &
                +sc(3) + sc(7)) ) / (hy*hz)

           ! sxyz
           slope(i,j,k,7) = (-sc(1) + sc(5) + sc(3) &
                +sc(2) - sc(7) - sc(6) &
                -sc(4) + sc(8) ) / (hx*hy*hz)

        enddo
     enddo
  enddo


end subroutine bdsslope_3d

subroutine bdsconc_3d(lo,hi,s,ng_s,ci,ng_ci,slope,ng_c,umac,vmac,wmac,ng_u,force,ng_f, &
     sedgex,sedgey,sedgez,ng_se,nw,dx,dt,bc,is_conservative)

  use bl_types
  implicit none

  integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_c,ng_u,ng_f,ng_se,nw,ng_ci
  real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :)
  real(kind=dp_t), intent(in   ) ::     ci(lo(1)-ng_ci:,lo(2)-ng_ci:,lo(3)-ng_ci:)
  real(kind=dp_t), intent(in   ) ::  slope(lo(1)-ng_c :,lo(2)-ng_c :,lo(3)-ng_c :,1:)
  real(kind=dp_t), intent(in   ) ::   umac(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :)
  real(kind=dp_t), intent(in   ) ::   vmac(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :)
  real(kind=dp_t), intent(in   ) ::   wmac(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :)
  real(kind=dp_t), intent(in   ) ::  force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :)
  real(kind=dp_t), intent(inout) :: sedgex(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:)
  real(kind=dp_t), intent(inout) :: sedgey(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:)
  real(kind=dp_t), intent(inout) :: sedgez(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:)
  real(kind=dp_t), intent(in   ) :: dx(:),dt
  integer        , intent(in   ) :: bc(3,2) ! bc(dir,lohi)
  logical        , intent(in   ) :: is_conservative

  ! local variables
  integer i,j,k,ioff,joff,koff,ll

  real(kind=dp_t), allocatable ::   ux(:,:,:)
  real(kind=dp_t), allocatable ::   vy(:,:,:)
  real(kind=dp_t), allocatable ::   wz(:,:,:)
  real(kind=dp_t), allocatable :: divu(:,:,:)

  real(kind=dp_t) :: isign,jsign,ksign,hx,hy,hz
  real(kind=dp_t) :: del(3),p1(3),p2(3),p3(3),p4(3)
  real(kind=dp_t) :: val1,val2,val3,val4,val5
  real(kind=dp_t) :: u,v,w,uu,vv,ww,gamma,gamma2
  real(kind=dp_t) :: dt2,dt3,dt4,half,sixth

  interface
     subroutine eval_3d(s,slope,del,val)
       use bl_types
       implicit none
       real(kind=dp_t), intent(in   ) :: s
       real(kind=dp_t), intent(in   ) :: slope(:)
       real(kind=dp_t), intent(in   ) :: del(:)
       real(kind=dp_t), intent(  out) :: val
     end subroutine eval_3d
  end interface

  allocate(  ux(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
  allocate(  vy(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
  allocate(  wz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
  allocate(divu(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))

  hx = dx(1)
  hy = dx(2)
  hz = dx(3)

  dt2 = dt/2.d0
  dt3 = dt/3.d0
  dt4 = dt/4.d0

  half = 0.5d0
  sixth = 1.d0/6.d0

  ! compute cell-centered ux, vy, and wz
  do k=lo(3)-1,hi(3)+1
     do j=lo(2)-1,hi(2)+1
        do i=lo(1)-1,hi(1)+1
           ux(i,j,k) = (umac(i+1,j,k) - umac(i,j,k)) / hx
           vy(i,j,k) = (vmac(i,j+1,k) - vmac(i,j,k)) / hy
           wz(i,j,k) = (wmac(i,j,k+1) - wmac(i,j,k)) / hz
           divu(i,j,k) = ux(i,j,k) + vy(i,j,k) + wz(i,j,k)
        end do
     end do
  end do

  ! compute sedgex on x-faces
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)+1

           ! compute sedgex without transverse corrections
           if (umac(i,j,k) .gt. 0) then
              isign = 1.d0
              ioff = -1
           else
              isign = -1.d0
              ioff = 0
           endif

           ! centroid of rectangular volume
           del(1) = isign*0.5d0*hx - 0.5d0*umac(i,j,k)*dt
           del(2) = 0.d0
           del(3) = 0.d0
           call eval_3d(s(i+ioff,j,k),slope(i+ioff,j,k,:),del,sedgex(i,j,k))

           ! source term
           if (is_conservative) then
              sedgex(i,j,k) = sedgex(i,j,k)* &
                   (1.d0 - dt2*ux(i+ioff,j,k)) + dt2*force(i+ioff,j,k)
           else
              sedgex(i,j,k) = sedgex(i,j,k)* &
                   (1.d0 + dt2*(vy(i+ioff,j,k)+wz(i+ioff,j,k))) + dt2*force(i+ioff,j,k)
           end if

           ! compute \Gamma^{y+} without corner corrections
           if (vmac(i+ioff,j+1,k) .gt. 0) then
              jsign = 1.d0
              joff = 0
           else
              jsign = -1.d0
              joff = 1
           endif

           u = 0.d0
           if (umac(i,j,k)*umac(i,j+joff,k) .gt. 0) then
              u = umac(i,j+joff,k)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = 0.d0

           p2(1) = isign*0.5d0*hx - umac(i,j,k)*dt
           p2(2) = jsign*0.5d0*hy
           p2(3) = 0.d0

           p3(1) = isign*0.5d0*hx - u*dt
           p3(2) = jsign*0.5d0*hy - vmac(i+ioff,j+1,k)*dt
           p3(3) = 0.d0

           do ll=1,3
              del(ll) = (p2(ll)+p3(ll))/2.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val1)

           do ll=1,3
              del(ll) = (p1(ll)+p3(ll))/2.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val2)

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll))/2.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val3)

           ! average these centroid values to get the average value
           gamma = (val1+val2+val3)/3.d0

           ! source term
           if (is_conservative) then
              gamma = gamma*(1.d0 - dt3*(ux(i+ioff,j+joff,k)+vy(i+ioff,j+joff,k)))
           else
              gamma = gamma*(1.d0 + dt3*wz(i+ioff,j+joff,k))
           end if

           ! correct \Gamma^{y+} with \Gamma^{y+,z+}
           if (wmac(i+ioff,j+joff,k+1) .gt. 0) then
              ksign = 1.d0
              koff = 0
           else
              ksign = -1.d0
              koff = 1
           endif

           uu = 0.d0
           if (umac(i,j,k)*umac(i,j+joff,k+koff) .gt. 0) then
              uu = umac(i,j+joff,k+koff)
           endif

           vv = 0.d0
           if (vmac(i+ioff,j+1,k)*vmac(i+ioff,j+1,k+koff) .gt. 0) then
              vv = vmac(i+ioff,j+1,k+koff)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx - umac(i,j,k)*dt
           p2(2) = jsign*0.5d0*hy
           p2(3) = ksign*0.5d0*hz

           p3(1) = isign*0.5d0*hx - umac(i,j,k)*dt
           p3(2) = jsign*0.5d0*hy - vmac(i+ioff,j+1,k)*dt
           p3(3) = ksign*0.5d0*hz

           p4(1) = isign*0.5d0*hx - uu*dt
           p4(2) = jsign*0.5d0*hy - vv*dt
           p4(3) = ksign*0.5d0*hz - wmac(i+ioff,j+joff,k+1)*dt

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

           do ll=1,3
              del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

           do ll=1,3
              del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

           gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

           ! source term
           if (is_conservative) then
              gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
           end if

           gamma2 = gamma2 * wmac(i+ioff,j+joff,k+1)

           gamma = gamma - dt*gamma2/(3.d0*hz)

           ! correct \Gamma^{y+} with \Gamma^{y+,z-}
           if (wmac(i+ioff,j+joff,k) .gt. 0) then
              ksign = 1.d0
              koff = -1
           else
              ksign = -1.d0
              koff = 0
           endif

           uu = 0.d0
           if (umac(i,j,k)*umac(i,j+joff,k+koff) .gt. 0) then
              uu = umac(i,j+joff,k+koff)
           endif

           vv = 0.d0
           if (vmac(i+ioff,j+1,k)*vmac(i+ioff,j+1,k+koff) .gt. 0) then
              vv = vmac(i+ioff,j+1,k+koff)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx - umac(i,j,k)*dt
           p2(2) = jsign*0.5d0*hy
           p2(3) = ksign*0.5d0*hz

           p3(1) = isign*0.5d0*hx - umac(i,j,k)*dt
           p3(2) = jsign*0.5d0*hy - vmac(i+ioff,j+1,k)*dt
           p3(3) = ksign*0.5d0*hz

           p4(1) = isign*0.5d0*hx - uu*dt
           p4(2) = jsign*0.5d0*hy - vv*dt
           p4(3) = ksign*0.5d0*hz - wmac(i+ioff,j+joff,k)*dt

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

           do ll=1,3
              del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

           do ll=1,3
              del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

           gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

           ! source term
           if (is_conservative) then
              gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
           end if

           gamma2 = gamma2 * wmac(i+ioff,j+joff,k)

           gamma = gamma + dt*gamma2/(3.d0*hz)

           ! correct sedgex with \Gamma^{y+}
           gamma = gamma * vmac(i+ioff,j+1,k)
           sedgex(i,j,k) = sedgex(i,j,k) - dt*gamma/(2.d0*hy)

           ! compute \Gamma^{y-} without corner corrections
           if (vmac(i+ioff,j,k) .gt. 0) then
              jsign = 1.d0
              joff = -1
           else
              jsign = -1.d0
              joff = 0
           endif

           u = 0.d0
           if (umac(i,j,k)*umac(i,j+joff,k) .gt. 0) then
              u = umac(i,j+joff,k)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = 0.d0

           p2(1) = isign*0.5d0*hx - umac(i,j,k)*dt
           p2(2) = jsign*0.5d0*hy
           p2(3) = 0.d0

           p3(1) = isign*0.5d0*hx - u*dt
           p3(2) = jsign*0.5d0*hy - vmac(i+ioff,j,k)*dt
           p3(3) = 0.d0

           do ll=1,3
              del(ll) = (p2(ll)+p3(ll))/2.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val1)

           do ll=1,3
              del(ll) = (p1(ll)+p3(ll))/2.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val2)

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll))/2.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val3)

           ! average these centroid values to get the average value
           gamma = (val1+val2+val3)/3.d0

           ! source term
           if (is_conservative) then
              gamma = gamma*(1.d0 - dt3*(ux(i+ioff,j+joff,k)+vy(i+ioff,j+joff,k)))
           else
              gamma = gamma*(1.d0 + dt3*wz(i+ioff,j+joff,k))
           end if

           ! correct \Gamma^{y-} with \Gamma^{y-,z+}
           if (wmac(i+ioff,j+joff,k+1) .gt. 0) then
              ksign = 1.d0
              koff = 0
           else
              ksign = -1.d0
              koff = 1
           endif

           uu = 0.d0
           if (umac(i,j,k)*umac(i,j+joff,k+koff) .gt. 0) then
              uu = umac(i,j+joff,k+koff)
           endif

           vv = 0.d0
           if (vmac(i+ioff,j,k)*vmac(i+ioff,j,k+koff) .gt. 0) then
              vv = vmac(i+ioff,j,k+koff)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx - umac(i,j,k)*dt
           p2(2) = jsign*0.5d0*hy
           p2(3) = ksign*0.5d0*hz

           p3(1) = isign*0.5d0*hx - umac(i,j,k)*dt
           p3(2) = jsign*0.5d0*hy - vmac(i+ioff,j,k)*dt
           p3(3) = ksign*0.5d0*hz

           p4(1) = isign*0.5d0*hx - uu*dt
           p4(2) = jsign*0.5d0*hy - vv*dt
           p4(3) = ksign*0.5d0*hz - wmac(i+ioff,j+joff,k+1)*dt

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

           do ll=1,3
              del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

           do ll=1,3
              del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

           gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

           ! source term
           if (is_conservative) then
              gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
           end if

           gamma2 = gamma2 * wmac(i+ioff,j+joff,k+1)

           gamma = gamma - dt*gamma2/(3.d0*hz)

           ! correct \Gamma^{y-} with \Gamma^{y-,z-}
           if (wmac(i+ioff,j+joff,k) .gt. 0) then
              ksign = 1.d0
              koff = -1
           else
              ksign = -1.d0
              koff = 0
           endif

           uu = 0.d0
           if (umac(i,j,k)*umac(i,j+joff,k+koff) .gt. 0) then
              uu = umac(i,j+joff,k+koff)
           endif

           vv = 0.d0
           if (vmac(i+ioff,j,k)*vmac(i+ioff,j,k+koff) .gt. 0) then
              vv = vmac(i+ioff,j,k+koff)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx - umac(i,j,k)*dt
           p2(2) = jsign*0.5d0*hy
           p2(3) = ksign*0.5d0*hz

           p3(1) = isign*0.5d0*hx - umac(i,j,k)*dt
           p3(2) = jsign*0.5d0*hy - vmac(i+ioff,j,k)*dt
           p3(3) = ksign*0.5d0*hz

           p4(1) = isign*0.5d0*hx - uu*dt
           p4(2) = jsign*0.5d0*hy - vv*dt
           p4(3) = ksign*0.5d0*hz - wmac(i+ioff,j+joff,k)*dt

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

           do ll=1,3
              del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

           do ll=1,3
              del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

           gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

           ! source term
           if (is_conservative) then
              gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
           end if

           gamma2 = gamma2 * wmac(i+ioff,j+joff,k)

           gamma = gamma + dt*gamma2/(3.d0*hz)

           ! correct sedgex with \Gamma^{y-}
           gamma = gamma * vmac(i+ioff,j,k)
           sedgex(i,j,k) = sedgex(i,j,k) + dt*gamma/(2.d0*hy)

           ! compute \Gamma^{z+} without corner corrections
           if (wmac(i+ioff,j,k+1) .gt. 0) then
              ksign = 1.d0
              koff = 0
           else
              ksign = -1.d0
              koff = 1
           endif

           u = 0.d0
           if (umac(i,j,k)*umac(i,j,k+koff) .gt. 0) then
              u = umac(i,j,k+koff)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = 0.d0
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx - umac(i,j,k)*dt
           p2(2) = 0.d0
           p2(3) = ksign*0.5d0*hz

           p3(1) = isign*0.5d0*hx - u*dt
           p3(2) = 0.d0
           p3(3) = ksign*0.5d0*hz - wmac(i+ioff,j,k+1)*dt

           do ll=1,3
              del(ll) = (p2(ll)+p3(ll))/2.d0
           end do
           call eval_3d(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = (p1(ll)+p3(ll))/2.d0
           end do
           call eval_3d(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll))/2.d0
           end do
           call eval_3d(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val3)

           ! average these centroid values to get the average value
           gamma = (val1+val2+val3)/3.d0

           ! source term
           if (is_conservative) then
              gamma = gamma*(1.d0 - dt3*(ux(i+ioff,j,k+koff)+wz(i+ioff,j,k+koff)))
           else
              gamma = gamma*(1.d0 + dt3*vy(i+ioff,j,k+koff))

           end if

           ! correct \Gamma^{z+} with \Gamma^{z+,y+}
           if (vmac(i+ioff,j+1,k+koff) .gt. 0) then
              jsign = 1.d0
              joff = 0
           else
              jsign = -1.d0
              joff = 1
           endif

           uu = 0.d0
           if (umac(i,j,k)*umac(i,j+joff,k+koff) .gt. 0) then
              uu = umac(i,j+joff,k+koff)
           endif

           ww = 0.d0
           if (wmac(i+ioff,j,k+1)*wmac(i+ioff,j+joff,k+1) .gt. 0) then
              ww = wmac(i+ioff,j+joff,k+1)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx - umac(i,j,k)*dt
           p2(2) = jsign*0.5d0*hy
           p2(3) = ksign*0.5d0*hz

           p3(1) = isign*0.5d0*hx - umac(i,j,k)*dt
           p3(2) = jsign*0.5d0*hy
           p3(3) = ksign*0.5d0*hz - wmac(i+ioff,j,k+1)*dt

           p4(1) = isign*0.5d0*hx - uu*dt
           p4(2) = jsign*0.5d0*hy - vmac(i+ioff,j+1,k+koff)*dt
           p4(3) = ksign*0.5d0*hz - ww*dt

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

           do ll=1,3
              del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

           do ll=1,3
              del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

           gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

           ! source term
           if (is_conservative) then
              gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
           end if

           gamma2 = gamma2 * vmac(i+ioff,j+1,k+koff)

           gamma = gamma - dt*gamma2/(3.d0*hy)

           ! correct \Gamma^{z+} with \Gamma^{z+,y-}
           if (vmac(i+ioff,j,k+koff) .gt. 0) then
              jsign = 1.d0
              joff = -1
           else
              jsign = -1.d0
              joff = 0
           endif

           uu = 0.d0
           if (umac(i,j,k)*umac(i,j+joff,k+koff) .gt. 0) then
              uu = umac(i,j+joff,k+koff)
           endif

           ww = 0.d0
           if (wmac(i+ioff,j,k+1)*wmac(i+ioff,j+joff,k+1) .gt. 0) then
              ww = wmac(i+ioff,j+joff,k+1)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx - umac(i,j,k)*dt
           p2(2) = jsign*0.5d0*hy
           p2(3) = ksign*0.5d0*hz

           p3(1) = isign*0.5d0*hx - umac(i,j,k)*dt
           p3(2) = jsign*0.5d0*hy
           p3(3) = ksign*0.5d0*hz - wmac(i+ioff,j,k+1)*dt

           p4(1) = isign*0.5d0*hx - uu*dt
           p4(2) = jsign*0.5d0*hy - vmac(i+ioff,j,k+koff)*dt
           p4(3) = ksign*0.5d0*hz - ww*dt

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

           do ll=1,3
              del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

           do ll=1,3
              del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

           gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

           ! source term
           if (is_conservative) then
              gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
           end if

           gamma2 = gamma2 * vmac(i+ioff,j,k+koff)

           gamma = gamma + dt*gamma2/(3.d0*hy)

           ! correct sedgex with \Gamma^{z+}
           gamma = gamma * wmac(i+ioff,j,k+1)
           sedgex(i,j,k) = sedgex(i,j,k) - dt*gamma/(2.d0*hz)

           ! compute \Gamma^{z-} without corner corrections
           if (wmac(i+ioff,j,k) .gt. 0) then
              ksign = 1.d0
              koff = -1
           else
              ksign = -1.d0
              koff = 0
           endif

           u = 0.d0
           if (umac(i,j,k)*umac(i,j,k+koff) .gt. 0) then
              u = umac(i,j,k+koff)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = 0.d0
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx - umac(i,j,k)*dt
           p2(2) = 0.d0
           p2(3) = ksign*0.5d0*hz

           p3(1) = isign*0.5d0*hx - u*dt
           p3(2) = 0.d0
           p3(3) = ksign*0.5d0*hz - wmac(i+ioff,j,k)*dt

           do ll=1,3
              del(ll) = (p2(ll)+p3(ll))/2.d0
           end do
           call eval_3d(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = (p1(ll)+p3(ll))/2.d0
           end do
           call eval_3d(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll))/2.d0
           end do
           call eval_3d(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val3)

           ! average these centroid values to get the average value
           gamma = (val1+val2+val3)/3.d0

           ! source term
           if (is_conservative) then
              gamma = gamma*(1.d0 - dt3*(ux(i+ioff,j,k+koff)+wz(i+ioff,j,k+koff)))
           else
              gamma = gamma*(1.d0 + dt3*vy(i+ioff,j,k+koff))
           end if

           ! correct \Gamma^{z-} with \Gamma^{z-,y+}
           if (vmac(i+ioff,j+1,k+koff) .gt. 0) then
              jsign = 1.d0
              joff = 0
           else
              jsign = -1.d0
              joff = 1
           endif

           uu = 0.d0
           if (umac(i,j,k)*umac(i,j+joff,k+koff) .gt. 0) then
              uu = umac(i,j+joff,k+koff)
           endif

           ww = 0.d0
           if (wmac(i+ioff,j,k)*wmac(i+ioff,j+joff,k) .gt. 0) then
              ww = wmac(i+ioff,j+joff,k)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx - umac(i,j,k)*dt
           p2(2) = jsign*0.5d0*hy
           p2(3) = ksign*0.5d0*hz

           p3(1) = isign*0.5d0*hx - umac(i,j,k)*dt
           p3(2) = jsign*0.5d0*hy
           p3(3) = ksign*0.5d0*hz - wmac(i+ioff,j,k)*dt

           p4(1) = isign*0.5d0*hx - uu*dt
           p4(2) = jsign*0.5d0*hy - vmac(i+ioff,j+1,k+koff)*dt
           p4(3) = ksign*0.5d0*hz - ww*dt

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

           do ll=1,3
              del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

           do ll=1,3
              del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

           gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

           ! source term
           if (is_conservative) then
              gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
           end if

           gamma2 = gamma2 * vmac(i+ioff,j+1,k+koff)

           gamma = gamma - dt*gamma2/(3.d0*hy)

           ! correct \Gamma^{z-} with \Gamma^{z-,y-}
           if (vmac(i+ioff,j,k+koff) .gt. 0) then
              jsign = 1.d0
              joff = -1
           else
              jsign = -1.d0
              joff = 0
           endif

           uu = 0.d0
           if (umac(i,j,k)*umac(i,j+joff,k+koff) .gt. 0) then
              uu = umac(i,j+joff,k+koff)
           endif

           ww = 0.d0
           if (wmac(i+ioff,j,k)*wmac(i+ioff,j+joff,k) .gt. 0) then
              ww = wmac(i+ioff,j+joff,k)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx - umac(i,j,k)*dt
           p2(2) = jsign*0.5d0*hy
           p2(3) = ksign*0.5d0*hz

           p3(1) = isign*0.5d0*hx - umac(i,j,k)*dt
           p3(2) = jsign*0.5d0*hy
           p3(3) = ksign*0.5d0*hz - wmac(i+ioff,j,k)*dt

           p4(1) = isign*0.5d0*hx - uu*dt
           p4(2) = jsign*0.5d0*hy - vmac(i+ioff,j,k+koff)*dt
           p4(3) = ksign*0.5d0*hz - ww*dt

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

           do ll=1,3
              del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

           do ll=1,3
              del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

           gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

           ! source term
           if (is_conservative) then
              gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
           end if

           gamma2 = gamma2 * vmac(i+ioff,j,k+koff)

           gamma = gamma + dt*gamma2/(3.d0*hy)

           ! correct sedgex with \Gamma^{z-}
           gamma = gamma * wmac(i+ioff,j,k)
           sedgex(i,j,k) = sedgex(i,j,k) + dt*gamma/(2.d0*hz)

        enddo
     enddo
  enddo

  ! compute sedgey on y-faces    
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)+1
        do i=lo(1),hi(1)

           ! compute sedgey without transverse corrections
           ! centroid of rectangular volume
           if (vmac(i,j,k) .gt. 0) then
              jsign = 1.d0
              joff = -1
           else
              jsign = -1.d0
              joff = 0
           endif

           del(1) = 0.d0
           del(2) = jsign*0.5d0*hy - 0.5d0*vmac(i,j,k)*dt
           del(3) = 0.d0
           call eval_3d(s(i,j+joff,k),slope(i,j+joff,k,:),del,sedgey(i,j,k))

           ! source term
           if (is_conservative) then
              sedgey(i,j,k) = sedgey(i,j,k)* &
                   (1.d0 - dt2*vy(i,j+joff,k)) + dt2*force(i,j+joff,k)
           else
              sedgey(i,j,k) = sedgey(i,j,k)* &
                   (1.d0 + dt2*(ux(i,j+joff,k)+wz(i,j+joff,k))) + dt2*force(i,j+joff,k)
           end if

           ! compute \Gamma^{x+} without corner corrections
           if (umac(i+1,j+joff,k) .gt. 0) then
              isign = 1.d0
              ioff = 0
           else
              isign = -1.d0
              ioff = 1
           endif

           v = 0.d0
           if (vmac(i,j,k)*vmac(i+ioff,j,k) .gt. 0) then
              v = vmac(i+ioff,j,k)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = 0.d0

           p2(1) = isign*0.5d0*hx
           p2(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
           p2(3) = 0.d0

           p3(1) = isign*0.5d0*hx - umac(i+1,j+joff,k)*dt
           p3(2) = jsign*0.5d0*hy - v*dt
           p3(3) = 0.d0

           do ll=1,3
              del(ll) = (p2(ll)+p3(ll))/2.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val1)

           do ll=1,3
              del(ll) = (p1(ll)+p3(ll))/2.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val2)

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll))/2.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val3)

           ! average these centroid values to get the average value
           gamma = (val1+val2+val3)/3.d0

           ! source term
           if (is_conservative) then
              gamma = gamma*(1.d0 - dt3*(vy(i+ioff,j+joff,k)+ux(i+ioff,j+joff,k)))
           else
              gamma = gamma*(1.d0 + dt3*wz(i+ioff,j+joff,k))
           end if

           ! correct \Gamma^{x+} with \Gamma^{x+,z+}
           if (wmac(i+ioff,j+joff,k+1) .gt. 0) then
              ksign = 1.d0
              koff = 0
           else
              ksign = -1.d0
              koff = 1
           endif

           vv = 0.d0
           if (vmac(i,j,k)*vmac(i+ioff,j,k+koff) .gt. 0) then
              vv = vmac(i+ioff,j,k+koff)
           endif

           uu = 0.d0
           if (umac(i+1,j+joff,k)*umac(i+1,j+joff,k+koff) .gt. 0) then
              uu = umac(i+1,j+joff,k+koff)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx
           p2(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
           p2(3) = ksign*0.5d0*hz

           p3(1) = isign*0.5d0*hx - umac(i+1,j+joff,k)*dt
           p3(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
           p3(3) = ksign*0.5d0*hz

           p4(1) = isign*0.5d0*hx - uu*dt
           p4(2) = jsign*0.5d0*hy - vv*dt
           p4(3) = ksign*0.5d0*hz - wmac(i+ioff,j+joff,k+1)*dt

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

           do ll=1,3
              del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

           do ll=1,3
              del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

           gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

           ! source term
           if (is_conservative) then
              gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
           end if

           gamma2 = gamma2 * wmac(i+ioff,j+joff,k+1)

           gamma = gamma - dt*gamma2/(3.d0*hz)

           ! correct \Gamma^{x+} with \Gamma^{x+,z-}
           if (wmac(i+ioff,j+joff,k) .gt. 0) then
              ksign = 1.d0
              koff = -1
           else
              ksign = -1.d0
              koff = 0
           endif

           vv = 0.d0
           if (vmac(i,j,k)*vmac(i+ioff,j,k+koff) .gt. 0) then
              vv = vmac(i+ioff,j,k+koff)
           endif

           uu = 0.d0
           if (umac(i+1,j+joff,k)*umac(i+1,j+joff,k+koff) .gt. 0) then
              uu = umac(i+1,j+joff,k+koff)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx
           p2(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
           p2(3) = ksign*0.5d0*hz

           p3(1) = isign*0.5d0*hx - umac(i+1,j+joff,k)*dt
           p3(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
           p3(3) = ksign*0.5d0*hz

           p4(1) = isign*0.5d0*hx - uu*dt
           p4(2) = jsign*0.5d0*hy - vv*dt
           p4(3) = ksign*0.5d0*hz - wmac(i+ioff,j+joff,k)*dt

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

           do ll=1,3
              del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

           do ll=1,3
              del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

           gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

           ! source term
           if (is_conservative) then
              gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
           end if

           gamma2 = gamma2 * wmac(i+ioff,j+joff,k)

           gamma = gamma + dt*gamma2/(3.d0*hz)

           ! correct sedgey with \Gamma^{x+}
           gamma = gamma * umac(i+1,j+joff,k)
           sedgey(i,j,k) = sedgey(i,j,k) - dt*gamma/(2.d0*hx)

           ! compute \Gamma^{x-} without corner corrections
           if (umac(i,j+joff,k) .gt. 0) then
              isign = 1.d0
              ioff = -1
           else
              isign = -1.d0
              ioff = 0
           endif

           v = 0.d0
           if (vmac(i,j,k)*vmac(i+ioff,j,k) .gt. 0) then
              v = vmac(i+ioff,j,k)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = 0.d0

           p2(1) = isign*0.5d0*hx
           p2(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
           p2(3) = 0.d0

           p3(1) = isign*0.5d0*hx - umac(i,j+joff,k)*dt
           p3(2) = jsign*0.5d0*hy - v*dt
           p3(3) = 0.d0

           do ll=1,3
              del(ll) = (p2(ll)+p3(ll))/2.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val1)

           do ll=1,3
              del(ll) = (p1(ll)+p3(ll))/2.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val2)

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll))/2.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val3)

           ! average these centroid values to get the average value
           gamma = (val1+val2+val3)/3.d0

           ! source term
           if (is_conservative) then
              gamma = gamma*(1.d0 - dt3*(vy(i+ioff,j+joff,k)+ux(i+ioff,j+joff,k)))
           else
              gamma = gamma*(1.d0 + dt3*wz(i+ioff,j+joff,k))
           end if

           ! correct \Gamma^{x-} with \Gamma^{x-,z+}
           if (wmac(i+ioff,j+joff,k+1) .gt. 0) then
              ksign = 1.d0
              koff = 0
           else
              ksign = -1.d0
              koff = 1
           endif

           vv = 0.d0
           if (vmac(i,j,k)*vmac(i+ioff,j,k+koff) .gt. 0) then
              vv = vmac(i+ioff,j,k+koff)
           endif

           uu = 0.d0
           if (umac(i,j+joff,k)*umac(i,j+joff,k+koff) .gt. 0) then
              uu = umac(i,j+joff,k+koff)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx
           p2(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
           p2(3) = ksign*0.5d0*hz

           p3(1) = isign*0.5d0*hx - umac(i,j+joff,k)*dt
           p3(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
           p3(3) = ksign*0.5d0*hz

           p4(1) = isign*0.5d0*hx - uu*dt
           p4(2) = jsign*0.5d0*hy - vv*dt
           p4(3) = ksign*0.5d0*hz - wmac(i+ioff,j+joff,k+1)*dt

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

           do ll=1,3
              del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

           do ll=1,3
              del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

           gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

           ! source term
           if (is_conservative) then
              gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
           end if

           gamma2 = gamma2 * wmac(i+ioff,j+joff,k+1)

           gamma = gamma - dt*gamma2/(3.d0*hz)

           ! correct \Gamma^{x-} with \Gamma^{x-,z-}
           if (wmac(i+ioff,j+joff,k) .gt. 0) then
              ksign = 1.d0
              koff = -1
           else
              ksign = -1.d0
              koff = 0
           endif

           vv = 0.d0
           if (vmac(i,j,k)*vmac(i+ioff,j,k+koff) .gt. 0) then
              vv = vmac(i+ioff,j,k+koff)
           endif

           uu = 0.d0
           if (umac(i,j+joff,k)*umac(i,j+joff,k+koff) .gt. 0) then
              uu = umac(i,j+joff,k+koff)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx
           p2(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
           p2(3) = ksign*0.5d0*hz

           p3(1) = isign*0.5d0*hx - umac(i,j+joff,k)*dt
           p3(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
           p3(3) = ksign*0.5d0*hz

           p4(1) = isign*0.5d0*hx - uu*dt
           p4(2) = jsign*0.5d0*hy - vv*dt
           p4(3) = ksign*0.5d0*hz - wmac(i+ioff,j+joff,k)*dt

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

           do ll=1,3
              del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

           do ll=1,3
              del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

           gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

           ! source term
           if (is_conservative) then
              gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
           end if

           gamma2 = gamma2 * wmac(i+ioff,j+joff,k)

           gamma = gamma + dt*gamma2/(3.d0*hz)

           ! correct sedgey with \Gamma^{x-}
           gamma = gamma * umac(i,j+joff,k)
           sedgey(i,j,k) = sedgey(i,j,k) + dt*gamma/(2.d0*hx)

           ! compute \Gamma^{z+} without corner corrections
           if (wmac(i,j+joff,k+1) .gt. 0) then
              ksign = 1.d0
              koff = 0
           else
              ksign = -1.d0
              koff = 1
           endif

           v = 0.d0
           if (vmac(i,j,k)*vmac(i,j,k+koff) .gt. 0) then
              v = vmac(i,j,k+koff)
           endif

           p1(1) = 0.d0
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = 0.d0
           p2(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
           p2(3) = ksign*0.5d0*hz

           p3(1) = 0.d0
           p3(2) = jsign*0.5d0*hy - v*dt
           p3(3) = ksign*0.5d0*hz - wmac(i,j+joff,k+1)*dt

           do ll=1,3
              del(ll) = (p2(ll)+p3(ll))/2.d0
           end do
           call eval_3d(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = (p1(ll)+p3(ll))/2.d0
           end do
           call eval_3d(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll))/2.d0
           end do
           call eval_3d(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val3)

           ! average these centroid values to get the average value
           gamma = (val1+val2+val3)/3.d0

           ! source term
           if (is_conservative) then
              gamma = gamma*(1.d0 - dt3*(vy(i,j+joff,k+koff)+wz(i,j+joff,k+koff)))
           else
              gamma = gamma*(1.d0 + dt3*ux(i,j+joff,k+koff))
           end if

           ! correct \Gamma^{z+} with \Gamma^{z+,x+}
           if (umac(i+1,j+joff,k+koff) .gt. 0) then
              isign = 1.d0
              ioff = 0
           else
              isign = -1.d0
              ioff = 1
           endif

           vv = 0.d0
           if (vmac(i,j,k)*vmac(i+ioff,j,k+koff) .gt. 0) then
              vv = vmac(i+ioff,j,k+koff)
           endif

           ww = 0.d0
           if (wmac(i,j+joff,k+1)*wmac(i+ioff,j+joff,k+1) .gt. 0) then
              ww = wmac(i+ioff,j+joff,k+1)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx
           p2(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
           p2(3) = ksign*0.5d0*hz

           p3(1) = isign*0.5d0*hx
           p3(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
           p3(3) = ksign*0.5d0*hz - wmac(i,j+joff,k+1)*dt

           p4(1) = isign*0.5d0*hx - umac(i+1,j+joff,k+koff)*dt
           p4(2) = jsign*0.5d0*hy - vv*dt
           p4(3) = ksign*0.5d0*hz - ww*dt

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

           do ll=1,3
              del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

           do ll=1,3
              del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

           gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

           ! source term
           if (is_conservative) then
              gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
           end if

           gamma2 = gamma2 * umac(i+1,j+joff,k+koff)

           gamma = gamma - dt*gamma2/(3.d0*hx)

           ! correct \Gamma^{z+} with \Gamma^{z+,x-}
           if (umac(i,j+joff,k+koff) .gt. 0) then
              isign = 1.d0
              ioff = -1
           else
              isign = -1.d0
              ioff = 0
           endif

           vv = 0.d0
           if (vmac(i,j,k)*vmac(i+ioff,j,k+koff) .gt. 0) then
              vv = vmac(i+ioff,j,k+koff)
           endif

           ww = 0.d0
           if (wmac(i,j+joff,k+1)*wmac(i+ioff,j+joff,k+1) .gt. 0) then
              ww = wmac(i+ioff,j+joff,k+1)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx
           p2(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
           p2(3) = ksign*0.5d0*hz

           p3(1) = isign*0.5d0*hx
           p3(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
           p3(3) = ksign*0.5d0*hz - wmac(i,j+joff,k+1)*dt

           p4(1) = isign*0.5d0*hx - umac(i,j+joff,k+koff)*dt
           p4(2) = jsign*0.5d0*hy - vv*dt
           p4(3) = ksign*0.5d0*hz - ww*dt

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

           do ll=1,3
              del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

           do ll=1,3
              del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

           gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

           ! source term
           if (is_conservative) then
              gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
           end if

           gamma2 = gamma2 * umac(i,j+joff,k+koff)

           gamma = gamma + dt*gamma2/(3.d0*hx)

           ! correct sedgey with \Gamma^{z+}
           gamma = gamma * wmac(i,j+joff,k+1)
           sedgey(i,j,k) = sedgey(i,j,k) - dt*gamma/(2.d0*hz)

           ! compute \Gamma^{z-} without corner corrections
           if (wmac(i,j+joff,k) .gt. 0) then
              ksign = 1.d0
              koff = -1
           else
              ksign = -1.d0
              koff = 0
           endif

           v = 0.d0
           if (vmac(i,j,k)*vmac(i,j,k+koff) .gt. 0) then
              v = vmac(i,j,k+koff)
           endif

           p1(1) = 0.d0
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = 0.d0
           p2(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
           p2(3) = ksign*0.5d0*hz

           p3(1) = 0.d0
           p3(2) = jsign*0.5d0*hy - v*dt
           p3(3) = ksign*0.5d0*hz - wmac(i,j+joff,k)*dt

           do ll=1,3
              del(ll) = (p2(ll)+p3(ll))/2.d0
           end do
           call eval_3d(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = (p1(ll)+p3(ll))/2.d0
           end do
           call eval_3d(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll))/2.d0
           end do
           call eval_3d(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val3)

           ! average these centroid values to get the average value
           gamma = (val1+val2+val3)/3.d0

           ! source term
           if (is_conservative) then
              gamma = gamma*(1.d0 - dt3*(vy(i,j+joff,k+koff)+wz(i,j+joff,k+koff)))
           else
              gamma = gamma*(1.d0 + dt3*ux(i,j+joff,k+koff))
           end if

           ! correct \Gamma^{z-} with \Gamma^{z-,x+}
           if (umac(i+1,j+joff,k+koff) .gt. 0) then
              isign = 1.d0
              ioff = 0
           else
              isign = -1.d0
              ioff = 1
           endif

           vv = 0.d0
           if (vmac(i,j,k)*vmac(i+ioff,j,k+koff) .gt. 0) then
              vv = vmac(i+ioff,j,k+koff)
           endif

           ww = 0.d0
           if (wmac(i,j+joff,k)*wmac(i+ioff,j+joff,k) .gt. 0) then
              ww = wmac(i+ioff,j+joff,k)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx
           p2(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
           p2(3) = ksign*0.5d0*hz

           p3(1) = isign*0.5d0*hx
           p3(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
           p3(3) = ksign*0.5d0*hz - wmac(i,j+joff,k)*dt

           p4(1) = isign*0.5d0*hx - umac(i+1,j+joff,k+koff)*dt
           p4(2) = jsign*0.5d0*hy - vv*dt
           p4(3) = ksign*0.5d0*hz - ww*dt

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

           do ll=1,3
              del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

           do ll=1,3
              del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

           gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

           ! source term
           if (is_conservative) then
              gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
           end if

           gamma2 = gamma2 * umac(i+1,j+joff,k+koff)

           gamma = gamma - dt*gamma2/(3.d0*hx)

           ! correct \Gamma^{z-} with \Gamma^{z-,x-}
           if (umac(i,j+joff,k+koff) .gt. 0) then
              isign = 1.d0
              ioff = -1
           else
              isign = -1.d0
              ioff = 0
           endif

           vv = 0.d0
           if (vmac(i,j,k)*vmac(i+ioff,j,k+koff) .gt. 0) then
              vv = vmac(i+ioff,j,k+koff)
           endif

           ww = 0.d0
           if (wmac(i,j+joff,k)*wmac(i+ioff,j+joff,k) .gt. 0) then
              ww = wmac(i+ioff,j+joff,k)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx
           p2(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
           p2(3) = ksign*0.5d0*hz

           p3(1) = isign*0.5d0*hx
           p3(2) = jsign*0.5d0*hy - vmac(i,j,k)*dt
           p3(3) = ksign*0.5d0*hz - wmac(i,j+joff,k)*dt

           p4(1) = isign*0.5d0*hx - umac(i,j+joff,k+koff)*dt
           p4(2) = jsign*0.5d0*hy - vv*dt
           p4(3) = ksign*0.5d0*hz - ww*dt

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

           do ll=1,3
              del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

           do ll=1,3
              del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

           gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

           ! source term
           if (is_conservative) then
              gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
           end if

           gamma2 = gamma2 * umac(i,j+joff,k+koff)

           gamma = gamma + dt*gamma2/(3.d0*hx)

           ! correct sedgey with \Gamma^{z-}
           gamma = gamma * wmac(i,j+joff,k)
           sedgey(i,j,k) = sedgey(i,j,k) + dt*gamma/(2.d0*hz)

        enddo
     enddo
  enddo

  ! compute sedgez on z-faces
  do k=lo(3),hi(3)+1
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)

           ! compute sedgez without transverse corrections
           ! centroid of rectangular volume
           if (wmac(i,j,k) .gt. 0) then
              ksign = 1.d0
              koff = -1
           else
              ksign = -1.d0
              koff = 0
           endif

           del(1) = 0.d0
           del(2) = 0.d0
           del(3) = ksign*0.5d0*hz - 0.5d0*wmac(i,j,k)*dt
           call eval_3d(s(i,j,k+koff),slope(i,j,k+koff,:),del,sedgez(i,j,k))

           ! source term
           if (is_conservative) then
              sedgez(i,j,k) = sedgez(i,j,k)* &
                   (1.d0 - dt2*wz(i,j,k+koff)) + dt2*force(i,j,k+koff)
           else
              sedgez(i,j,k) = sedgez(i,j,k)* &
                   (1.d0 + dt2*(ux(i,j,k+koff)+vy(i,j,k+koff))) + dt2*force(i,j,k+koff)
           end if

           ! compute \Gamma^{x+} without corner corrections
           if (umac(i+1,j,k+koff) .gt. 0) then
              isign = 1.d0
              ioff = 0
           else
              isign = -1.d0
              ioff = 1
           endif

           w = 0.d0
           if (wmac(i,j,k)*wmac(i+ioff,j,k) .gt. 0) then
              w = wmac(i+ioff,j,k)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = 0.d0
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx
           p2(2) = 0.d0
           p2(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

           p3(1) = isign*0.5d0*hx - umac(i+1,j,k+koff)*dt
           p3(2) = 0.d0
           p3(3) = ksign*0.5d0*hz - w*dt

           do ll=1,3
              del(ll) = (p2(ll)+p3(ll))/2.d0
           end do
           call eval_3d(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = (p1(ll)+p3(ll))/2.d0
           end do
           call eval_3d(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll))/2.d0
           end do
           call eval_3d(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val3)

           ! average these centroid values to get the average value
           gamma = (val1+val2+val3)/3.d0

           ! source term
           if (is_conservative) then
              gamma = gamma*(1.d0 - dt3*(wz(i+ioff,j,k+koff)+ux(i+ioff,j,k+koff)))
           else
              gamma = gamma*(1.d0 + dt3*vy(i+ioff,j,k+koff))
           end if

           ! correct \Gamma^{x+} with \Gamma^{x+,y+}
           if (vmac(i+ioff,j+1,k+koff) .gt. 0) then
              jsign = 1.d0
              joff = 0
           else
              jsign = -1.d0
              joff = 1
           endif

           ww = 0.d0
           if (wmac(i,j,k)*wmac(i+ioff,j+joff,k) .gt. 0) then
              ww = wmac(i+ioff,j+joff,k)
           endif

           uu = 0.d0
           if (umac(i+1,j,k+koff)*umac(i+1,j+joff,k+koff) .gt. 0) then
              uu = umac(i+1,j+joff,k+koff)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx
           p2(2) = jsign*0.5d0*hy
           p2(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

           p3(1) = isign*0.5d0*hx - umac(i+1,j+joff,k)*dt
           p3(2) = jsign*0.5d0*hy
           p3(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

           p4(1) = isign*0.5d0*hx - uu*dt
           p4(2) = jsign*0.5d0*hy - vmac(i+ioff,j+1,k+koff)*dt
           p4(3) = ksign*0.5d0*hz - ww*dt

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

           do ll=1,3
              del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

           do ll=1,3
              del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

           gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

           ! source term
           if (is_conservative) then
              gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
           end if

           gamma2 = gamma2 * vmac(i+ioff,j+1,k+koff)

           gamma = gamma - dt*gamma2/(3.d0*hy)

           ! correct \Gamma^{x+} with \Gamma^{x+,y-}
           if (vmac(i+ioff,j,k+koff) .gt. 0) then
              jsign = 1.d0
              joff = -1
           else
              jsign = -1.d0
              joff = 0
           endif

           ww = 0.d0
           if (wmac(i,j,k)*wmac(i+ioff,j+joff,k) .gt. 0) then
              ww = wmac(i+ioff,j+joff,k)
           endif

           uu = 0.d0
           if (umac(i+1,j,k+koff)*umac(i+1,j+joff,k+koff) .gt. 0) then
              uu = umac(i+1,j+joff,k+koff)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx
           p2(2) = jsign*0.5d0*hy
           p2(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

           p3(1) = isign*0.5d0*hx - umac(i+1,j+joff,k)*dt
           p3(2) = jsign*0.5d0*hy
           p3(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

           p4(1) = isign*0.5d0*hx - uu*dt
           p4(2) = jsign*0.5d0*hy - vmac(i+ioff,j,k+koff)*dt
           p4(3) = ksign*0.5d0*hz - ww*dt

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

           do ll=1,3
              del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

           do ll=1,3
              del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

           gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

           ! source term
           if (is_conservative) then
              gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
           end if

           gamma2 = gamma2 * vmac(i+ioff,j,k+koff)

           gamma = gamma + dt*gamma2/(3.d0*hy)

           ! correct sedgez with \Gamma^{x+}
           gamma = gamma * umac(i+1,j,k+koff)
           sedgez(i,j,k) = sedgez(i,j,k) - dt*gamma/(2.d0*hx)

           ! compute \Gamma^{x-} without corner corrections
           if (umac(i,j,k+koff) .gt. 0) then
              isign = 1.d0
              ioff = -1
           else
              isign = -1.d0
              ioff = 0
           endif

           w = 0.d0
           if (wmac(i,j,k)*wmac(i+ioff,j,k) .gt. 0) then
              w = wmac(i+ioff,j,k)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = 0.d0
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx
           p2(2) = 0.d0
           p2(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

           p3(1) = isign*0.5d0*hx - umac(i,j,k+koff)*dt
           p3(2) = 0.d0
           p3(3) = ksign*0.5d0*hz - w*dt

           do ll=1,3
              del(ll) = (p2(ll)+p3(ll))/2.d0
           end do
           call eval_3d(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = (p1(ll)+p3(ll))/2.d0
           end do
           call eval_3d(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll))/2.d0
           end do
           call eval_3d(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val3)

           ! average these centroid values to get the average value
           gamma = (val1+val2+val3)/3.d0

           ! source term
           if (is_conservative) then
              gamma = gamma*(1.d0 - dt3*(wz(i+ioff,j,k+koff)+ux(i+ioff,j,k+koff)))
           else
              gamma = gamma*(1.d0 + dt3*vy(i+ioff,j,k+koff))
           end if

           ! correct \Gamma^{x-} with \Gamma^{x-,y+}
           if (vmac(i+ioff,j+1,k+koff) .gt. 0) then
              jsign = 1.d0
              joff = 0
           else
              jsign = -1.d0
              joff = 1
           endif

           ww = 0.d0
           if (wmac(i,j,k)*wmac(i+ioff,j+joff,k) .gt. 0) then
              ww = wmac(i+ioff,j+joff,k)
           endif

           uu = 0.d0
           if (umac(i,j,k+koff)*umac(i,j+joff,k+koff) .gt. 0) then
              uu = umac(i,j+joff,k+koff)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx
           p2(2) = jsign*0.5d0*hy
           p2(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

           p3(1) = isign*0.5d0*hx - umac(i,j+joff,k)*dt
           p3(2) = jsign*0.5d0*hy
           p3(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

           p4(1) = isign*0.5d0*hx - uu*dt
           p4(2) = jsign*0.5d0*hy - vmac(i+ioff,j+1,k+koff)*dt
           p4(3) = ksign*0.5d0*hz - ww*dt

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

           do ll=1,3
              del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

           do ll=1,3
              del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

           gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

           ! source term
           if (is_conservative) then
              gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
           end if

           gamma2 = gamma2 * vmac(i+ioff,j+1,k+koff)

           gamma = gamma - dt*gamma2/(3.d0*hy)

           ! correct \Gamma^{x-} with \Gamma^{x-,y-}
           if (vmac(i+ioff,j,k+koff) .gt. 0) then
              jsign = 1.d0
              joff = -1
           else
              jsign = -1.d0
              joff = 0
           endif

           ww = 0.d0
           if (wmac(i,j,k)*wmac(i+ioff,j+joff,k) .gt. 0) then
              ww = wmac(i+ioff,j+joff,k)
           endif

           uu = 0.d0
           if (umac(i,j,k+koff)*umac(i,j+joff,k+koff) .gt. 0) then
              uu = umac(i,j+joff,k+koff)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx
           p2(2) = jsign*0.5d0*hy
           p2(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

           p3(1) = isign*0.5d0*hx - umac(i,j+joff,k)*dt
           p3(2) = jsign*0.5d0*hy
           p3(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

           p4(1) = isign*0.5d0*hx - uu*dt
           p4(2) = jsign*0.5d0*hy - vmac(i+ioff,j,k+koff)*dt
           p4(3) = ksign*0.5d0*hz - ww*dt

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

           do ll=1,3
              del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

           do ll=1,3
              del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

           gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

           ! source term
           if (is_conservative) then
              gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
           end if

           gamma2 = gamma2 * vmac(i+ioff,j,k+koff)

           gamma = gamma + dt*gamma2/(3.d0*hy)

           ! correct sedgez with \Gamma^{x-}
           gamma = gamma * umac(i,j,k+koff)
           sedgez(i,j,k) = sedgez(i,j,k) + dt*gamma/(2.d0*hx)

           ! compute \Gamma^{y+} without corner corrections
           if (vmac(i,j+1,k+koff) .gt. 0) then
              jsign = 1.d0
              joff = 0
           else
              jsign = -1.d0
              joff = 1
           endif

           w = 0.d0
           if (wmac(i,j,k)*wmac(i,j+joff,k) .gt. 0) then
              w = wmac(i,j+joff,k)
           endif

           p1(1) = 0.d0
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = 0.d0
           p2(2) = jsign*0.5d0*hy
           p2(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

           p3(1) = 0.d0
           p3(2) = jsign*0.5d0*hy - vmac(i,j+1,k+koff)*dt
           p3(3) = ksign*0.5d0*hz - w*dt

           do ll=1,3
              del(ll) = (p2(ll)+p3(ll))/2.d0
           end do
           call eval_3d(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = (p1(ll)+p3(ll))/2.d0
           end do
           call eval_3d(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll))/2.d0
           end do
           call eval_3d(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val3)

           ! average these centroid values to get the average value
           gamma = (val1+val2+val3)/3.d0

           ! source term
           if (is_conservative) then
              gamma = gamma*(1.d0 - dt3*(wz(i,j+joff,k+koff)+vy(i,j+joff,k+koff)))
           else
              gamma = gamma*(1.d0 + dt3*ux(i,j+joff,k+koff))
           end if

           ! correct \Gamma^{y+} with \Gamma^{y+,x+}
           if (umac(i+1,j+joff,k+koff) .gt. 0) then
              isign = 1.d0
              ioff = 0
           else
              isign = -1.d0
              ioff = 1
           endif

           ww = 0.d0
           if (wmac(i,j,k)*wmac(i+ioff,j+joff,k) .gt. 0) then
              ww = wmac(i+ioff,j+joff,k)
           endif

           vv = 0.d0
           if (vmac(i,j+1,k+koff)*vmac(i+ioff,j+1,k+koff) .gt. 0) then
              vv = vmac(i+ioff,j+1,k+koff)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx
           p2(2) = jsign*0.5d0*hy
           p2(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

           p3(1) = isign*0.5d0*hx
           p3(2) = jsign*0.5d0*hy - vmac(i+ioff,j+1,k)*dt
           p3(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

           p4(1) = isign*0.5d0*hx - umac(i+1,j+joff,k+koff)*dt
           p4(2) = jsign*0.5d0*hy - vv*dt
           p4(3) = ksign*0.5d0*hz - ww*dt

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

           do ll=1,3
              del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

           do ll=1,3
              del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

           gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

           ! source term
           if (is_conservative) then
              gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
           end if

           gamma2 = gamma2 * umac(i+1,j+joff,k+koff)

           gamma = gamma - dt*gamma2/(3.d0*hx)

           ! correct \Gamma^{y+} with \Gamma^{y+,x-}
           if (umac(i,j+joff,k+koff) .gt. 0) then
              isign = 1.d0
              ioff = -1
           else
              isign = -1.d0
              ioff = 0
           endif

           ww = 0.d0
           if (wmac(i,j,k)*wmac(i+ioff,j+joff,k) .gt. 0) then
              ww = wmac(i+ioff,j+joff,k)
           endif

           vv = 0.d0
           if (vmac(i,j+1,k+koff)*vmac(i+ioff,j+1,k+koff) .gt. 0) then
              vv = vmac(i+ioff,j+1,k+koff)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx
           p2(2) = jsign*0.5d0*hy
           p2(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

           p3(1) = isign*0.5d0*hx
           p3(2) = jsign*0.5d0*hy - vmac(i+ioff,j+1,k)*dt
           p3(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

           p4(1) = isign*0.5d0*hx - umac(i,j+joff,k+koff)*dt
           p4(2) = jsign*0.5d0*hy - vv*dt
           p4(3) = ksign*0.5d0*hz - ww*dt

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

           do ll=1,3
              del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

           do ll=1,3
              del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

           gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

           ! source term
           if (is_conservative) then
              gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
           end if

           gamma2 = gamma2 * umac(i,j+joff,k+koff)

           gamma = gamma + dt*gamma2/(3.d0*hx)

           ! correct sedgez with \Gamma^{y+}
           gamma = gamma * vmac(i,j+1,k+koff)
           sedgez(i,j,k) = sedgez(i,j,k) - dt*gamma/(2.d0*hy)

           ! compute \Gamma^{y-} without corner corrections
           if (vmac(i,j,k+koff) .gt. 0) then
              jsign = 1.d0
              joff = -1
           else
              jsign = -1.d0
              joff = 0
           endif

           w = 0.d0
           if (wmac(i,j,k)*wmac(i,j+joff,k) .gt. 0) then
              w = wmac(i,j+joff,k)
           endif

           p1(1) = 0.d0
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = 0.d0
           p2(2) = jsign*0.5d0*hy
           p2(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

           p3(1) = 0.d0
           p3(2) = jsign*0.5d0*hy - vmac(i,j,k+koff)*dt
           p3(3) = ksign*0.5d0*hz - w*dt

           do ll=1,3
              del(ll) = (p2(ll)+p3(ll))/2.d0
           end do
           call eval_3d(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = (p1(ll)+p3(ll))/2.d0
           end do
           call eval_3d(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll))/2.d0
           end do
           call eval_3d(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val3)

           ! average these centroid values to get the average value
           gamma = (val1+val2+val3)/3.d0

           ! source term
           if (is_conservative) then
              gamma = gamma*(1.d0 - dt3*(wz(i,j+joff,k+koff)+vy(i,j+joff,k+koff)))
           else
              gamma = gamma*(1.d0 + dt3*ux(i,j+joff,k+koff))
           end if

           ! correct \Gamma^{y-} with \Gamma^{y-,x+}
           if (umac(i+1,j+joff,k+koff) .gt. 0) then
              isign = 1.d0
              ioff = 0
           else
              isign = -1.d0
              ioff = 1
           endif

           ww = 0.d0
           if (wmac(i,j,k)*wmac(i+ioff,j+joff,k) .gt. 0) then
              ww = wmac(i+ioff,j+joff,k)
           endif

           vv = 0.d0
           if (vmac(i,j,k+koff)*vmac(i+ioff,j,k+koff) .gt. 0) then
              vv = vmac(i+ioff,j,k+koff)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx
           p2(2) = jsign*0.5d0*hy
           p2(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

           p3(1) = isign*0.5d0*hx
           p3(2) = jsign*0.5d0*hy - vmac(i+ioff,j,k)*dt
           p3(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

           p4(1) = isign*0.5d0*hx - umac(i+1,j+joff,k+koff)*dt
           p4(2) = jsign*0.5d0*hy - vv*dt
           p4(3) = ksign*0.5d0*hz - ww*dt

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

           do ll=1,3
              del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

           do ll=1,3
              del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

           gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

           ! source term
           if (is_conservative) then
              gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
           end if

           gamma2 = gamma2 * umac(i+1,j+joff,k+koff)

           gamma = gamma - dt*gamma2/(3.d0*hx)

           ! correct \Gamma^{y-} with \Gamma^{y-,x-}
           if (umac(i,j+joff,k+koff) .gt. 0) then
              isign = 1.d0
              ioff = -1
           else
              isign = -1.d0
              ioff = 0
           endif

           ww = 0.d0
           if (wmac(i,j,k)*wmac(i+ioff,j+joff,k) .gt. 0) then
              ww = wmac(i+ioff,j+joff,k)
           endif

           vv = 0.d0
           if (vmac(i,j,k+koff)*vmac(i+ioff,j,k+koff) .gt. 0) then
              vv = vmac(i+ioff,j,k+koff)
           endif

           p1(1) = isign*0.5d0*hx
           p1(2) = jsign*0.5d0*hy
           p1(3) = ksign*0.5d0*hz

           p2(1) = isign*0.5d0*hx
           p2(2) = jsign*0.5d0*hy
           p2(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

           p3(1) = isign*0.5d0*hx
           p3(2) = jsign*0.5d0*hy - vmac(i+ioff,j,k)*dt
           p3(3) = ksign*0.5d0*hz - wmac(i,j,k)*dt

           p4(1) = isign*0.5d0*hx - umac(i,j+joff,k+koff)*dt
           p4(2) = jsign*0.5d0*hy - vv*dt
           p4(3) = ksign*0.5d0*hz - ww*dt

           do ll=1,3
              del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

           do ll=1,3
              del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

           do ll=1,3
              del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

           do ll=1,3
              del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

           do ll=1,3
              del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
           end do
           call eval_3d(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

           gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

           ! source term
           if (is_conservative) then
              gamma2 = gamma2*(1.d0 - dt4*divu(i+ioff,j+joff,k+koff))
           end if

           gamma2 = gamma2 * umac(i,j+joff,k+koff)

           gamma = gamma + dt*gamma2/(3.d0*hx)

           ! correct sedgez with \Gamma^{y-}
           gamma = gamma * vmac(i,j,k+koff)
           sedgez(i,j,k) = sedgez(i,j,k) + dt*gamma/(2.d0*hy)

        enddo
     enddo
  enddo

  deallocate(ux,vy,wz,divu)

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
     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           sedgex(i,j,k) = s(i-1,j,k)
        end do
     end do
  else if (bc(1,1) .eq. -1) then
     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           sedgex(i,j,k) = 0.d0
        end do
     end do
  end if

  i = hi(1)+1
  if (bc(1,2) .eq. 3) then
     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           sedgex(i,j,k) = s(i,j,k)
        end do
     end do
  else if (bc(1,2) .eq. -1) then
     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           sedgex(i,j,k) = 0.d0
        end do
     end do
  end if

  j = lo(2)
  if (bc(2,1) .eq. 3) then
     do k=lo(3),hi(3)
        do i=lo(1),hi(1)
           sedgey(i,j,k) = s(i,j-1,k)
        end do
     end do
  else if (bc(2,1) .eq. -1) then
     do k=lo(3),hi(3)
        do i=lo(1),hi(1)
           sedgey(i,j,k) = 0.d0
        end do
     end do
  end if

  j = hi(2)+1
  if (bc(2,2) .eq. 3) then
     do k=lo(3),hi(3)
        do i=lo(1),hi(1)
           sedgey(i,j,k) = s(i,j,k)
        end do
     end do
  else if (bc(2,2) .eq. -1) then
     do k=lo(3),hi(3)
        do i=lo(1),hi(1)
           sedgey(i,j,k) = 0.d0
        end do
     end do
  end if

  k = lo(3)
  if (bc(3,1) .eq. 3) then
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           sedgez(i,j,k) = s(i,j,k-1)
        end do
     end do
  else if (bc(3,1) .eq. -1) then
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           sedgez(i,j,k) = 0.d0
        end do
     end do
  end if

  k = hi(3)+1
  if (bc(3,2) .eq. 3) then
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           sedgez(i,j,k) = s(i,j,k)
        end do
     end do
  else if (bc(3,2) .eq. -1) then
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           sedgez(i,j,k) = 0.d0
        end do
     end do
  end if

end subroutine bdsconc_3d

subroutine eval_3d(s,slope,del,val)

  use bl_types
  implicit none

  real(kind=dp_t), intent(in   ) :: s
  real(kind=dp_t), intent(in   ) :: slope(:)
  real(kind=dp_t), intent(in   ) :: del(:)
  real(kind=dp_t), intent(  out) :: val

  val = s + del(1)*slope(1) + del(2)*slope(2) + del(3)*slope(3) &
       + del(1)*del(2)*slope(4) + del(1)*del(3)*slope(5) + del(2)*del(3)*slope(6) &
       + del(1)*del(2)*del(3)*slope(7)

end subroutine eval_3d
