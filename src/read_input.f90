! Read input information
! DEVELOPER
!   Hom Nath Gharti
!   hngharti_AT_gmail_DOT_com
!
! HISTORY
!   Oct 24,2006, HNG: original package created
!   Apr 25,2017, HNG: modified package
subroutine read_input(inp_fname)
use global
implicit none
character(len=250),intent(in) :: inp_fname

integer :: ios,slen
character(len=250) :: line

! open file to read
open(unit=11,file=trim(inp_fname),status='old', action='read',iostat=ios)
if (ios /= 0)then
  write(*,*)'ERROR: input file "'//trim(inp_fname)//'" cannot be opened!'
  stop
endif

call read_nextline
slen=len_trim(line)
read(line(3:slen),*)de_a

call read_nextline
slen=len_trim(line)
read(line(3:slen),*)de_b

call read_nextline
slen=len_trim(line)
read(line(3:slen),*)de_c

call read_nextline
slen=len_trim(line)
read(line(3:slen),*)de_f

call read_nextline
slen=len_trim(line)
read(line(8:slen),*)length

call read_nextline
slen=len_trim(line)
read(line(20:slen),*)nelmt

call read_nextline
slen=len_trim(line)
read(line(29:slen),*)nenode

call read_nextline
slen=len_trim(line)
read(line(4:slen),*)bc_u0

call read_nextline
slen=len_trim(line)
read(line(4:slen),*)bc_qL

call read_nextline
slen=len_trim(line)
read(line(31:slen),*)pexact
close(11)

write(stdout,*)'Differential equation parameters'
write(stdout,*)'a:',de_a
write(stdout,*)'b:',de_b
write(stdout,*)'c:',de_c
write(stdout,*)'f:',de_f
write(stdout,*)'Boundary conditions'
write(stdout,*)'u(0):',bc_u0
write(stdout,*)'q(L):',bc_qL
write(stdout,*)'Model and mesh parameters'
write(stdout,*)'Length:',length
write(stdout,*)'Number of elements:',nelmt
write(stdout,*)'Number of nodes per element:',nenode
write(stdout,*)'----------------------------------------'
return
contains
subroutine read_nextline
implicit none
integer,parameter :: nmax_comment=20
integer :: i
do i=1,nmax_comment
  read(11,'(a)',iostat=ios)line
  if(ios.ne.0)then
    write(*,*)'ERROR: reading input file!'
    stop
  endif
  if(line(1:1).eq.'#' .or. len(trim(line)).lt.1)then
    cycle 
  else
    return
  endif
enddo
end subroutine read_nextline
!===============================================
end subroutine read_input
!===============================================
