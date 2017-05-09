! Read the boundary conditions and modify the RHS vector accordingly. Stiffness
! matrix will be modified in the later stage.
! DEVELOPER
!   Hom Nath Gharti
!   hngharti_AT_gmail_DOT_com
!
! HISTORY
!   Oct 24,2006, HNG: original package created
!   Apr 25,2017, HNG: modified package
module boundary_condition

contains

subroutine apply_boundary_condition
use global
implicit none
integer :: i,j,inode

! define u(0)
inode=1 ! first node
node(inode)%uid=1
node(inode)%u=bc_u0

! define q(L)
inode=nnode ! last node
global_fid(inode)=1
global_fvec(inode)=global_fvec(inode)+bc_qL

! modify rhs for defined primary variable bcs
! following is the general algorithm
! we can modify the dof associated with only the first node!
ngdof=nndof*nnode
do i=1,nnode
  if(node(i)%uid==1)then
    ngdof=ngdof-1
    do j=1,nnode
      if(j/=i)then
        global_fvec(j)=global_fvec(j)-global_kmat(j,i)*node(i)%u
      endif
    enddo
  endif
enddo

return
end subroutine apply_boundary_condition
!------------------------------------------------------------------------------
end module boundary_condition
!==============================================================================
