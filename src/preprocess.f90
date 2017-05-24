! Formation of local and global matrices and vectors
! DEVELOPER
!   Hom Nath Gharti
!   hngharti_AT_gmail_DOT_com
!
! HISTORY
!   Oct 24,2006, HNG: original package created
!   Apr 25,2017, HNG: modified package
module preprocess

contains

subroutine form_matrix_vector
use global
use integration
implicit none
integer :: i,j,k,kk,elm,ngpoint,row,col,gdfr,gdfc,gnoder,gnodec
real(kind=kreal) :: jac,jaci
real(kind=kreal),dimension(10) :: phi,dphi_dxi,xi,w
real(kind=kreal),allocatable :: local_fvec(:,:),local_kmat(:,:,:)

allocate(local_fvec(nedof,nelmt),local_kmat(nedof,nedof,nelmt))

! integration points
ngpoint=nenode+1 ! Number of gauss points, maximum=6
call gauss(ngpoint,xi,w)

! Local matrices
do elm=1,nelmt
  do i=1,nedof
    local_fvec(i,elm)=0.0_kreal
    do j=1,nedof
      local_kmat(i,j,elm)=0.0_kreal
      do k=1,ngpoint
        call lagrange(nenode,xi(k),phi,dphi_dxi)
        call jacobian(elm,nenode,dphi_dxi,jac)
        jaci=1/jac ! For scalar jac
        local_kmat(i,j,elm)=local_kmat(i,j,elm)+(de_a*dphi_dxi(i)*dphi_dxi(j)*jaci)*w(k) &
        +(-de_b*phi(i)*dphi_dxi(j))*w(k)+(-de_c*phi(i)*phi(j)*jac)*w(k)
        if(j==i)then
          local_fvec(i,elm)=local_fvec(i,elm)+(de_f*phi(i)*jac)*w(k)
        endif
      enddo
    enddo
  enddo
enddo
!stop
do i=1,nelmt
  write(*,*)'----------------------------------------'
  write(*,*)'Local stiffness matrix for element:',i
  write(*,*)'----------------------------------------'
  do j=1,nedof
    write(*,'(3(f15.10,1x))')local_kmat(:,j,i)
  enddo
enddo

do i=1,nelmt
  write(*,*)'----------------------------------------'
  write(*,*)'Local force vector for element:',i
  write(*,*)'----------------------------------------'
  write(*,'(f10.4)')local_fvec(:,i)
enddo

! Global matrices
do elm=1,nelmt
  row=0
  do i=1,nedof ! Nodes of an element
    gnoder=elmt(elm)%node(i)
    gdfr=1*(gnoder-1)

    do k=1,nndof ! Degrees of freedom at a node
      row=row+1
      col=0
      do j=1,nedof
        gnodec=elmt(elm)%node(j)
        gdfc=1*(gnodec-1)
        do kk=1,nndof
          col=col+1
          global_kmat(gdfr+k,gdfc+kk)=global_kmat(gdfr+k,gdfc+kk)+local_kmat(row,col,elm)
        enddo
      enddo
      global_fvec(gdfr+k)=global_fvec(gdfr+k)+local_fvec(row,elm)
    enddo
  enddo
enddo

deallocate(local_fvec,local_kmat)

write(*,*)'----------------------------------------'
write(*,*)'Global stiffness matrix'
write(*,*)'----------------------------------------'
do i=1,nnode
  write(*,'(11(f7.4,1x))')global_kmat(i,1:nnode)
enddo

write(*,*)'----------------------------------------'
write(*,*)'Global force vector'
write(*,*)'----------------------------------------'
write(*,'(f7.4)')global_fvec(1:nnode)

return

end subroutine form_matrix_vector
!------------------------------------------------------------------------------
end module preprocess
!==============================================================================
