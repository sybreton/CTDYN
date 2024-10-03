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
  real(dp) :: edr, xe1, xde1
  real(dp) :: xa1, xa2, xda1, xda2
  real(dp) :: xb, gd
  real(dp) :: dd1, xc1, c2_h2, oco

  !-------------------------------------------
  !> brent 
  !-------------------------------------------
  real(dp) :: al_i, al_f, accu

  !-------------------------------------------
  !> outputs
  !-------------------------------------------
  character(len=128) :: dir
  integer :: nj
  logical :: write_vectors
  real(dp) :: zeta_r
  real(dp) :: xbt, xbo

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
  real(dp) :: hd, ffree, xm, c3, bct
  real(dp) :: aqu 
  real(dp) :: beta_i

  !-------------------------------------------
  !> controls 
  !-------------------------------------------
  integer :: nso
  real(dp) :: flg
  real(dp) :: cm_f, cm_i, rm_f, rm_i

  ! Additional variables that are used by several 
  ! subroutines. These were included in "common" 
  ! blocks in the old version of the code. 
  integer :: ii, it
  real(dp) :: etep, etet, eep
  real(dp) :: rt, imag, co, c_u, beta, betb
  real(dp) :: reg(10),ieg(10)
  real(dp) :: gam
  character(len=1) :: jobvr

  namelist /global/ sr, rotp

  namelist /profiles/ regime, s0, s2, s4, s6, &
                      a2p, a4p, xa1, xa2, xda1, &
                      xda2, xb, gd, edr, xe1, xde1, &
                      dd1, xc1, c2_h2, oco

  namelist /brent/ al_i, al_f, accu

  namelist /outputs/ dir, nj, write_vectors, zeta_r, &
                     xbt, xbo

  namelist /fields/ degree, mmm

  namelist /physics/ hd, ffree, xm, aqu, c3, bct, &
                     beta_i

  namelist /boundaries/ x_in

  namelist /controls/ rm_i, rm_f, cm_i, cm_f, nso, flg

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
    s0      =  0.85d0
    s2      =  0.08d0
    s4      =  0.d0
    s6      =  0.d0
    a2p     =  1.d0
    a4p     =  0.d0
    xa1     =  0.64d0   
    xa2     =  0.72d0   
    xda1    =  0.025d0
    xda2    =  0.025d0
    xb      =  0.65d0
    gd      =  1.2d0 
    edr     =  0.1d0     
    xe1     =  0.70d0     
    xde1    =  0.025d0   
    dd1     =  0.05d0
    xc1     =  0.70d0
    c2_h2   =  0.2d0
    oco     =  0.9d0
  
    !-------------------------------------------
    !> brent 
    !-------------------------------------------
    al_i =  0.d0 
    al_f =  10.d0      
    accu =  1d-3 
  
    !-------------------------------------------
    !> outputs
    !-------------------------------------------
    dir  = 'tests/test_default'
    write_vectors = .true.
    zeta_r  = 1.3d0
    xbt     = 0.75d0
    xbo     = 1.5d0
    nj      = 8
  
    !-------------------------------------------
    !> boundaries
    !-------------------------------------------
    x_in = 0.58d0
  
    !-------------------------------------------
    !> fields
    !-------------------------------------------
    degree = 'd'   
    mmm    = 0
  
    !-------------------------------------------
    !> physics
    !-------------------------------------------
    hd      = 1
    ffree   = 0
    xm      = -0.45d0
    aqu     =  1
    c3      =  0
    bct     =  1.d0
    beta_i  =  0
  
    !-------------------------------------------
    !> controls 
    !-------------------------------------------
    rm_i    =  0.d0      
    rm_f    =  400.d0     
    cm_i    =  1000.d0      
    cm_f    =  40000.d0        
    nso     =  0    
    flg     =  0.d0 

  end subroutine 

  subroutine read_namelist (inlist)
    !> Read namelists variables provided in input file.
    !> Namelist variables are declared in cio.f90 and
    !> are therefore accessible to any module that import
    !> it.
    character(len=128) :: inlist

    integer :: fu

    call initialise_namelist_values
    call getarg(1, inlist)
    ! Open and close to read in any order
    open(newunit=fu, file=inlist, status="old", action="read")
    read(unit=fu, nml=global)
    close(fu)
    open(newunit=fu, file=inlist, status="old", action="read")
    read(unit=fu, nml=profiles)
    close(fu)
    open(newunit=fu, file=inlist, status="old", action="read")
    read(unit=fu, nml=brent)
    close(fu)
    open(newunit=fu, file=inlist, status="old", action="read")
    read(unit=fu, nml=boundaries)
    close(fu)
    open(newunit=fu, file=inlist, status="old", action="read")
    read(unit=fu, nml=fields)
    close(fu)
    open(newunit=fu, file=inlist, status="old", action="read")
    read(unit=fu, nml=physics)
    close(fu)
    open(newunit=fu, file=inlist, status="old", action="read")
    read(unit=fu, nml=outputs)
    close(fu)
    open(newunit=fu, file=inlist, status="old", action="read")
    read(unit=fu, nml=controls)
    close(fu)

  end subroutine
  
  subroutine write_namelist (inlist)
    !> Write full set of namelists variables.
    !> Namelist variables are declared in cio.f90 and
    !> are therefore accessible to any module that import
    !> it.
    character(len=128) :: inlist

    integer :: fu

    ! Open and close to read in any order
    open(newunit=fu, file=inlist, status="replace", & 
         action="write")
    write(unit=fu, nml=global, delim="quote")
    write(unit=fu, nml=profiles, delim="quote")
    write(unit=fu, nml=brent, delim="quote")
    write(unit=fu, nml=boundaries, delim="quote")
    write(unit=fu, nml=fields, delim="quote")
    write(unit=fu, nml=physics, delim="quote")
    write(unit=fu, nml=outputs, delim="quote")
    write(unit=fu, nml=controls, delim="quote")
    close(fu)

  end subroutine

end module cio


