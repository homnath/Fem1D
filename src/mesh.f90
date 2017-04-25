! Create mesh
! FEM analysis for 1-D second order ordinary differential equation
! Date    Programmer      Status
!========== =================== ================
!2006/10/24   Hom Nath Gharti   General Lagrange elements
!-----------------------------------------------
! This program solves the 1-D second order differential equation:
! a(d2u/dx2)+b(du/dx)+cu=f, u(0)=0, a(du/dx) at (x=L) = q at (x=L)
! where, both u and f are function of x only

subroutine meshgen
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
end subroutine meshgen
!===============================================
