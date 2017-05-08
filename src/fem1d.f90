!FEM analysis of 1-D second order ordinary differential equationi using
!general Lagrange elements
!
!PURPOSE
!  This program solves the 1-D second order differential equation:
!  a d2u/dx2 + b du/dx + c u = f, u(0) = 0, a du/dx at (x = L) = q
!  where, both u and f are functions of x only.
!
!DEVELOPER
!  Hom Nath Gharti
!  formerly at Institute of Engineering, Tribhuvan University, Nepal
!  formerly at NORSAR, Norway
!  Department of Geosciences, Princeton University, USA
!  hngharti_AT_gmail_DOT_com
!
!HISTORY
!  Oct 24,2006, HNG: original package creates
!  Apr 25,2017, HNG: modified package

!lines below are helpful to compile a single file if Makefile is not used
!include './global.f90' ! Global variables
!include './mesh.f90' ! Mesh information
!include './interpolate.f90'
!include './form.f90' ! Formation of matrices
!include './bc.f90' ! Boundary conditions
!include './solve.f90' ! Solution of matrices
!include './plot.f90' ! Plotting the result

program fem1d
use global
implicit none
character(len=250) :: arg1,inp_fname,prog

call get_command_argument(0, prog)
!----input and initialisation----
if (command_argument_count() <= 0) then
  write(*,*)'ERROR: no input file!'
  stop
endif

call get_command_argument(1, arg1)
if(trim(arg1)==('--help'))then
  write(stdout,'(a)')'Usage: '//trim(prog)//' [Options] [input_file]'
  write(stdout,'(a)')'Options:'
  write(stdout,'(a)')'    --help        : Display this information.'
  write(stdout,'(a)')'    --version     : Display version information.'
  stop
elseif(trim(arg1)==('--version'))then
  write(stdout,'(a)')'Fem1D 1.0.0'
  write(stdout,'(a)')'This is free software; see the source for copying '
  write(stdout,'(a)')'conditions.  There is NO warranty; not even for '
  write(stdout,'(a)')'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.'
  stop
endif

! get input file name
call get_command_argument(1, inp_fname)

!-----------------
! Preprocess
!-----------------
print*,'Reading input....'
call read_input(inp_fname)
print*,'Creating mesh....'
call mesh

!-----------------
! Process
!-----------------
print*,'Forming matrices.....'
call form
print*,'Applying boundary conditions.......'
call apply_bc
print*,'Solving..........'
call solve

!-----------------
! Post process
!-----------------
print*,'Plotting..................'
call plot

print*,'Finished successfully'
print*,'----------------------------------------'
end program fem1d
!===============================================
