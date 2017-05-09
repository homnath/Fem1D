! Create mesh
! DEVELOPER
!   Hom Nath Gharti
!   hngharti_AT_gmail_DOT_com
!
! HISTORY
!   Oct 24,2006, HNG: original package created
!   Apr 25,2017, HNG: modified package
module mesh
contains
subroutine create_mesh
use global
implicit none
integer :: i,j
real(kind=kreal) :: dx

! number of elemental degrees of freedoms
nedof=nndof*nenode

! number of nodes
nnode=nelmt*(nenode-1)+1

! length of the small part
dx=length/(nnode-1)

! nodal co-ordinates
print*,'----------------------------------------'
print*,'Coordinates'
print*,'----------------------------------------'
do i=1,nnode
  node(i)%x=(i-1)*dx
  print*,node(i)%x
end do
print*,'----------------------------------------'

! connectivity
print*,'Connectivity'
print*,'----------------------------------------'

do i=1,nelmt
  do j=1,nenode
    elmt(i)%node(j)=(i-1)*(nenode-1)+j
  enddo
  print*,(elmt(i)%node(j),j=1,nenode)
enddo
print*,'----------------------------------------'

return
end subroutine create_mesh
end module mesh
!===============================================
