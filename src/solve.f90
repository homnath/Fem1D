! Solution of matrix equation A x = b
! ngdof : total global degrees of freedoms (size of the vector)
! amat  : stiffness matrix of size ngdof X ngdof
! bvec  : load vector of size ngdof which becomes a solution vector at the end
! DEVELOPER
!   Hom Nath Gharti
!   hngharti_AT_gmail_DOT_com
!
! HISTORY
!   Oct 24,2006, HNG: original package created
!   Apr 25,2017, HNG: modified package
subroutine solve
use global
implicit none
real(kind=kreal), dimension(ngdof,ngdof) :: amat
real(kind=kreal), dimension(ngdof) :: bvec

integer :: i,j,k,kk,row,col,peak,num
real(kind=kreal) :: temp, factor
real(kind=kreal), dimension(ngdof) :: temp1

! form the coefficient matrix and load vector
row=0
col=0
do i=1,nnode
  if(node(i)%uid/=1)then
    row=row+1
    col=0
    do j=1,nnode
      if(node(j)%uid/=1)then
        col=col+1
        amat(row,col)=global_kmat(i,j)
        bvec(row)=global_fvec(i)
      endif
    enddo
  endif
enddo
print*,'----------------------------------------'
print*,'Global stiffness matrix'
print*,'----------------------------------------'
do i=1,row
  print 30,amat(i,1:col)
30  format(1x,20(f7.4))
enddo

print*,'----------------------------------------'
print*,'Global force vector'
print*,'----------------------------------------'
print 40,bvec(1:row)
40  format(1x,f7.4)

print*,'----------------------------------------'

! solver starts here
! forward elimination
do k=1,ngdof
  ! find the pivotal row
  peak=k
  do kk=k+1,ngdof
    if(abs(amat(kk,k))>abs(amat(peak,k)))then
      peak=kk
    endif
  enddo

  ! check for singularity
  if(abs(amat(peak,k))<eps)then
    write(*,*)'Equation is singular'
    stop
  endif

  ! swap
  if(peak/=k)then
    temp1=amat(peak,1:ngdof)
    amat(peak,1:ngdof)=amat(k,1:ngdof)
    amat(k,1:ngdof)=temp1
    temp=bvec(peak)
    bvec(peak)=bvec(k)
    bvec(k)=temp
  endif

  ! forward elimination
  do i=1,ngdof
    if(i/=k)then
      factor=amat(i,k)/amat(k,k)
      amat(i,1:ngdof)=amat(i,1:ngdof)-factor*amat(k,1:ngdof)
      bvec(i)=bvec(i)-factor*bvec(k)
    endif
  enddo
enddo

! substitution
do i=1,ngdof
  bvec(i)=bvec(i)/amat(i,i)
enddo
! solver ends here

! map to nodal freedoms
num=0
do i=1,nnode
  if(node(i)%uid/=1)then
    num=num+1
    node(i)%u=bvec(num)
  endif
enddo

return
end subroutine solve
!===============================================
