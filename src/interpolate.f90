! Lagrange interpolation functions
! DEVELOPER
!   Hom Nath Gharti
!   hngharti_AT_gmail_DOT_com
!
! HISTORY
!   Oct 24,2006, HNG: original package created
!   Apr 25,2017, HNG: modified package
subroutine lagrange(nenode,xi,phi,dphi_dxi)
use set_precision
implicit none
integer,intent(in) :: nenode
integer :: i,j,k,ldf ! ldf=Local degrees of freedom
real(kind=kreal),intent(in) :: xi
real(kind=kreal),dimension(nenode),intent(out) :: phi,dphi_dxi
real(kind=kreal),dimension(nenode) :: xii,term,dterm,sum_term
real(kind=kreal) :: dx

! Compute natural coordnates
dx=2.0_kreal/(nenode-1)! Length =2.0 as xi is taken -1 to +1
do i=1,nenode
  ! Coordinates when origin is in the left
  xii(i)=(i-1)*dx
enddo

! Origin is tranformed to mid point
xii=xii-1.0_kreal

ldf=nenode ! For lagrange straight element
do i=1,ldf
  k=0
  phi(i)=1.0_kreal
  do j=1,ldf
    if(j/=i)then
      k=k+1
      term(k)=(xi-xii(j))/(xii(i)-xii(j))
      dterm(k)=1.0_kreal/(xii(i)-xii(j)) ! Derivative of the term wrt xi

      phi(i)=phi(i)*(xi-xii(j))/(xii(i)-xii(j))
    endif
  enddo

  sum_term=1.0_kreal
  do j=1,ldf-1
    do k=1,ldf-1
      if(k==j)then
        sum_term(j)=sum_term(j)*dterm(k)
      else
        sum_term(j)=sum_term(j)*term(k)
      endif
    enddo
  enddo
  dphi_dxi(i)=0.0_kreal
  do j=1,nenode-1
    dphi_dxi(i)=dphi_dxi(i)+sum_term(j)
  enddo
enddo

return
end subroutine lagrange
!===============================================

! gauss quadrature rule
! formation of local and global matrices
subroutine gauss(ngpoint,xi,w)
use set_precision
implicit none
integer,intent(in) :: ngpoint ! Number of Gauss points
real(kind=kreal),dimension(ngpoint) :: xi,w

! gauss points and weights
if(ngpoint<1)then
  print*,ngpoint,' Integration points impossible!'
  stop
elseif(ngpoint==1)then
  xi(1)=0.0_kreal
  w(1)=2.0_kreal
elseif(ngpoint==2)then
  xi(1)=-0.5773502692_kreal
  w(1)=1.0_kreal
  xi(1)=-xi(1)
  w(2)=w(1)
elseif(ngpoint==3)then
  xi(1)=-0.7745966692_kreal
  w(1)=0.5555555555_kreal
  xi(2)=0.0_kreal
  w(2)=0.8888888889_kreal
  xi(3)=-xi(1)
  w(3)=w(1)
elseif(ngpoint==4)then
  xi(1)=-0.8611363116_kreal
  w(1)=0.3478548451_kreal
  xi(2)=-0.3399810435_kreal
  w(2)=0.6521451548_kreal
  xi(3)=-xi(2)
  w(3)=w(2)
  xi(4)=-xi(1)
  w(4)=w(1)
elseif(ngpoint==5)then
  xi(1)=-0.9061798459_kreal
  w(1)=0.2369268850_kreal
  xi(2)=-0.5384693101_kreal
  w(2)=0.4786286705_kreal
  xi(3)=0.0_kreal
  w(3)=0.5688888889_kreal
  xi(4)=-xi(2)
  w(4)=w(2)
  xi(5)=-xi(1)
  w(5)=w(1)
elseif(ngpoint==6)then
  xi(1)=-0.9324695142_kreal
  w(1)=0.1713244924_kreal
  xi(2)=-0.6612093865_kreal
  w(2)=0.3607615730_kreal
  xi(3)=-0.2386191861_kreal
  w(3)=0.4679139346_kreal
  xi(4)=-xi(3)
  w(4)=w(3)
  xi(5)=-xi(2)
  w(3)=w(2)
  xi(6)=-xi(1)
  w(6)=w(1)
else
  print*,'Modify the gauss integration for more than six points!'
  stop
endif
end subroutine gauss
!===============================================

subroutine jacobian(elm,nenod,dphi_dxi,jac)
use global
implicit none
integer,intent(in) :: elm,nenod
integer :: i
real(kind=kreal),intent(out) :: jac
real(kind=kreal),dimension(nenod) :: dphi_dxi

jac=0.0_kreal
do i=1,nenod
  jac=jac+node(elmt(elm)%node(i))%x*dphi_dxi(i)
enddo
return
end subroutine jacobian
!===============================================
