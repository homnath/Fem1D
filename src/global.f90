!This module defines global variables/parameters
!DEVELOPER
!  Hom Nath Gharti
!  hngharti_AT_gmail_DOT_com
!
!HISTORY
!  Apr 25,2017, HNG: modified package
!  Oct 24,2006, HNG: original package created
module global
use set_precision
implicit none
save
real(kind=kreal),parameter :: zero=0.0_kreal

! differential equation
! a(d2u/dx2) + b(du/dx) + cu = f
real(kind=kreal) :: de_a,de_b,de_c,de_f
! boundary conditions
real(kind=kreal) :: bc_u0,bc_qL ! u(0), q(L)

! plot exact result
integer :: pexact

integer,parameter :: maxnode=5000 ! Maximum number of nodes
integer,parameter :: enode=10 ! Maximum nodes in an element
integer,parameter :: eldf=20 ! Local degrees of freedom in an element
integer :: nnode,nelmt ! Number of nodes and number of elements
integer :: nenode ! Element type
integer,parameter :: nndof=1 ! Number of nodal degrees of freedoms
integer :: nedof ! Number of elemental degrees of freedoms
integer :: ngdof ! Maximum degrees of freedom
real(kind=kreal),parameter :: eps=1.0E-12_kreal
real(kind=kreal),parameter :: PI=3.141592653589793238462643383279502884197_kreal
real(kind=kreal) :: length ! total length of the domain

! Nodal specifications
type :: nodenum
  real(kind=kreal) :: x ! Coordinates of the node
  real(kind=kreal) :: u ! Displacement of the node
  integer :: uid ! 1=Known,0=Unknown
end type nodenum
type (nodenum),dimension(maxnode) :: node

! Elemental specifications
type :: element
  integer, dimension(enode) :: node ! Connectivity
end type element
type (element),dimension(maxnode) :: elmt

! Global matrices
real(kind=kreal),dimension(maxnode,maxnode) :: global_kmat
real(kind=kreal),dimension(maxnode) :: global_fvec
integer,dimension(maxnode) :: global_fid ! 1=Known, 0=Unknown

integer :: stdout=6
end module global
!===============================================
