! Plot the result
! DEVELOPER
!   Hom Nath Gharti
!   hngharti_AT_gmail_DOT_com
!
! HISTORY
!   Oct 24,2006, HNG: original package created
!   Apr 25,2017, HNG: modified package
subroutine plot
use global
implicit none
integer :: i
real(kind=kreal) :: u,exactu,rer,ernorm,erave

open(unit=20,file='./output/plot.dat',status='replace',action='write')
if(pexact.gt.0)then
  write(20,*)'# x-coordinates, Computed values, Analytical values'
else
  write(20,*)'# x-coordinates, Computed values'
endif

open(unit=10,file='./output/plotTP.dat',status='replace',action='write')
write(10,*)'VARIABLES = "X", "Computed"'
write(10,10)nnode
10  format(1x,'ZONE I=',i3,', J=1, K=1,F=POINT')

erave=zero; ernorm=zero; rer=zero
do i=1,nnode
  write(10,'(1x,2(f12.6,1x))')node(i)%x,node(i)%u
  
  if(pexact.gt.0)then
    u=exactu(node(i)%x)
    write(20,'(1x,3(f12.6,1x))')node(i)%x,node(i)%u,u
    if(u.ne.zero)then
      rer=abs(node(i)%u-exactu(node(i)%x))/exactu(node(i)%x)
    endif
    erave=erave+rer
    ernorm=ernorm+(rer*rer)
  else
    write(20,'(1x,2(f12.6,1x))')node(i)%x,node(i)%u
  endif
enddo
close(20)

if(pexact.gt.0)then
  write(*,30)erave/nnode,sqrt(ernorm)
  30  format(1x,'Average error = ',f10.6,' Error norm = ',f10.6)

  write(10,*)'VARIABLES = "X", "Exact"'
  write(10,40)nnode
  40  format(1x,'ZONE I=',i3,', J=1, K=1,F=POINT')

  do i=1,nnode
    write(10,'(1x,2(f18.6,1x))')node(i)%x,exactu(node(i)%x)
  enddo
endif
close(10)
return
end subroutine plot
!===============================================

! function to compute exact displacement u
real(kind=kreal) function exactu(x)
use set_precision
implicit none
real(kind=kreal),intent(in) :: x

! this expression is only valid for
! a d2u/dx2 + b du/dx + c u = f, u(0) = 0, a du/dx at (x = L) = q
! where, a=1, b=0, c=-1, f=10, u0=0, q=10, L=10
exactu=10.0_8*exp(-x)*(-exp(10.0_8)-exp(20.0_8)+exp(x)- &
 exp(2.0*x)+exp(2*(5+x))+exp(20+x))/(1.0_8+exp(20.0_8))
return
end function exactu
!===============================================
