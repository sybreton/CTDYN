module cio

  !---------------------------------------------------------
  !
  ! Common variables are listed in this module.
  !
  !---------------------------------------------------------

  integer :: np ! mesh points
  integer :: nt ! matrix dimension 
  integer :: na !     na = 1,3,5,... radial order in an and an
  integer :: nb !     nb = 2,4,6,... radial order in bn and bn
  
  parameter(na=49)     ! harmonics
  parameter(np=50)     ! mesh points 
  
  parameter(nb=na+1)           ! harmonics b number
  parameter(nt=nb*np)          ! global matrix dimension
  
  integer :: lwork
  parameter(lwork=540000)      
  
  integer :: nueg
  parameter(nueg=8)
  
  integer :: n_theta
  parameter(n_theta=301) 
  
  integer :: n_time
  parameter(n_time=300)
  
  real :: pi
  parameter(pi=3.14159265359)

  ! Additional variables that are used by several 
  ! subroutines. These were included in "common" 
  ! blocks in the old version of the code. 
  real :: mmm
  real :: vtu, rt, imag, co, c_u, beta, ffree, betb
  real :: etep, etet, xbt, xbo
  real :: reg(10),ieg(10)
  real :: eep
  integer :: ii, it
  real :: xa1, xa2, xa3, xb, xda1, xda2
  real :: edr, xe1, xde1, hd
  real :: x_in, bct, c3, sr, rotp, gd, aqu, flg
  real :: dd1, rc1, rc2, oco
  real :: s0, s2, s4, s6, a2p, a4p, xm
  real :: gam, zeta_r
  character(len=128) :: dir
  character*8 :: ans1, ans2, ans3, ans4
  character*2 :: jobvr, jobvl


end module cio
