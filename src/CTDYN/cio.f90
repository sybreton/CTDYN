module cio

  !---------------------------------------------------------
  !
  ! Common variables are listed in this module.
  !
  !---------------------------------------------------------

  use kind_parameters

  implicit none
  public

  integer, parameter :: np=50 ! mesh points
  integer, parameter :: na=49 !     na = 1,3,5,... radial order in an and an
  integer, parameter :: nb=na+1 !     nb = 2,4,6,... radial order in bn and bn
  integer, parameter :: nt=nb*np ! matrix dimension 
  
  integer, parameter :: lwork=540000
  integer, parameter :: nueg=8
  integer, parameter :: n_theta=301
  integer, parameter :: n_time=300
  real(dp), parameter :: pi=3.14159265359_dp

  !-------------------------------------------
  !-------------------------------------------
  ! Namelist variables and default values
  !-------------------------------------------
  !-------------------------------------------
  !> global
  !-------------------------------------------
  real(dp) :: sr, rotp
 
  !-------------------------------------------
  !> profiles
  !-------------------------------------------
  character(len=8) :: regime
  real(dp) :: s0, s2, s4, s6, a2p, a4p

  !-------------------------------------------
  !> brent 
  !-------------------------------------------
  real(dp) :: al_i, al_f, accu

  !-------------------------------------------
  !> outputs
  !-------------------------------------------
  character(len=128) :: dir
  logical :: write_vectors

  !-------------------------------------------
  !> boundaries
  !-------------------------------------------
  real(dp) :: x_in

  !-------------------------------------------
  !> fields
  !-------------------------------------------
  character(len=8) :: degree
  real(dp) :: mmm

  !-------------------------------------------
  !> physics
  !-------------------------------------------
  real(dp) :: hd, ffree, xm

  !-------------------------------------------
  !> other 
  !-------------------------------------------
  integer :: nso
  character(len=2) :: jobvr, jobvl
  real(dp) :: zeta_r
  real(dp) :: xbt, xbo
  real(dp) :: xa1, xa2, xa3, xb, xda1, xda2
  real(dp) :: bct, c3, aqu, flg, gd
  real(dp) :: edr, xe1, xde1
  real(dp) :: dd1, rc1, rc2, oco
  real(dp) :: beta_f, beta_i, beta_s 
  real(dp) :: cm_f, cm_i, rm_f, rm_i

  ! Additional variables that are used by several 
  ! subroutines. These were included in "common" 
  ! blocks in the old version of the code. 
  integer :: ii, it
  real(dp) :: etep, etet, eep
  real(dp) :: rt, imag, co, c_u, beta, betb
  real(dp) :: reg(10),ieg(10)
  real(dp) :: gam

contains

  subroutine initialise_namelist_values 
    !-------------------------------------------
    ! Initialise namelist to default values
    !-------------------------------------------

    !-------------------------------------------
    !> global
    !-------------------------------------------
    sr = 1. 
    rotp = 1.
   
    !-------------------------------------------
    !> profiles
    !-------------------------------------------
    regime  = 'h2'
    s0      =  0.85
    s2      =  0.08
    s4      =  0.
    s6      =  0.
    a2p     =  1.
    a4p     =  0.
  
    !-------------------------------------------
    !> brent 
    !-------------------------------------------
    al_i =  0. 
    al_f =  10.      
    accu =  1d-3 
  
    !-------------------------------------------
    !> outputs
    !-------------------------------------------
    dir  = 'tests/test_default'
    write_vectors = .true.
  
    !-------------------------------------------
    !> boundaries
    !-------------------------------------------
    x_in = 0.58
  
    !-------------------------------------------
    !> fields
    !-------------------------------------------
    degree = 'd'   
    mmm = 0
  
    !-------------------------------------------
    !> physics
    !-------------------------------------------
    hd = 1
    ffree = 0
    xm = -0.45
  
    !-------------------------------------------
    !> other 
    !-------------------------------------------
    xa1     =  0.64   
    xa2     =  0.72   
    xa3     =  2.0000000   
    xb      =  0.65
    xda1    =  0.025
    xda2    =  0.025
    rm_i    =  400.0       
    rm_f    =  400.0     
    cm_i    =  40000.0       
    cm_f    =  40000.0        
    nso     =  0    
    edr     =  0.1     
    xe1     =  0.70     
    xde1    =  0.025   
    c3      =  0
    bct     =  1.
    gd      =  1.2 
    aqu     =  1
    flg     =  0.
    dd1     =  0.05
    rc1     =  0.70
    rc2     =  0.2
    oco     =  0.9
    beta_i  =  0
    beta_f  =  0
    beta_s  =  0.1
    zeta_r  = 1.0
    xbt     =  0.75
    xbo     =  1.5

  end subroutine 
  
end module cio


