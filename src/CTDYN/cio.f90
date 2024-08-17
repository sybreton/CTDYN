module cio

  !---------------------------------------------------------
  !
  ! Kind parameters and common variables are listed 
  ! in this module.
  !
  !---------------------------------------------------------

   implicit none
   public

   !> Single precision real numbers, 6 digits, range 10^-307 to 10^307-1; 32 bits
   integer, parameter :: sp = selected_real_kind(6, 37)
   !> Double precision real numbers, 15 digits, range 10-307 to 10^307-1; 64 bits
   integer, parameter :: dp = selected_real_kind(15, 307)
   !> Quadruple precision real numbers, 33 digits, range 10^-4931 to 10^4931-1; 128 bits
   integer, parameter :: qp = selected_real_kind(33, 4931)

   !> Char length for integers, range -2^7 to 2^7-1; 8 bits
   integer, parameter :: i1 = selected_int_kind(2)
   !> Short length for integers, range -2^15 to 2^15-1; 16 bits
   integer, parameter :: i2 = selected_int_kind(4)
   !> Length of default integers, range -2^31 to 2^31-1; 32 bits
   integer, parameter :: i4 = selected_int_kind(9)
   !> Long length for integers, range -2^63 to 2^63-1; 64 bits
   integer, parameter :: i8 = selected_int_kind(18)

  integer, parameter :: np=50 ! mesh points
  integer, parameter :: na=49 !     na = 1,3,5,... radial order in an and an
  integer, parameter :: nb=na+1 !     nb = 2,4,6,... radial order in bn and bn
  integer, parameter :: nt=nb*np ! matrix dimension 
  
  integer, parameter :: lwork=540000
  integer, parameter :: nueg=8
  integer, parameter :: n_theta=301
  integer, parameter :: n_time=300
  real(dp), parameter :: pi=3.14159265359_dp

  ! Additional variables that are used by several 
  ! subroutines. These were included in "common" 
  ! blocks in the old version of the code. 
  integer :: ii, it
  real(dp) :: mmm
  real(dp) :: rt, imag, co, c_u, beta, ffree, betb
  real(dp) :: etep, etet, xbt, xbo
  real(dp) :: reg(10),ieg(10)
  real(dp) :: eep
  real(dp) :: xa1, xa2, xa3, xb, xda1, xda2
  real(dp) :: edr, xe1, xde1, hd
  real(dp) :: x_in, bct, c3, sr, rotp, gd, aqu, flg
  real(dp) :: dd1, rc1, rc2, oco
  real(dp) :: s0, s2, s4, s6, a2p, a4p, xm
  real(dp) :: gam, zeta_r
  character(len=128) :: dir
  character*8 :: ans1, ans2, ans3, ans4
  character*2 :: jobvr, jobvl

end module cio
