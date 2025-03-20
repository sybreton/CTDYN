module cio

  !---------------------------------------------------------
  !
  ! Common variables are listed in this module.
  !
  !---------------------------------------------------------

  use kind_parameters

  implicit none
  public

  integer(i4), parameter :: lwork   = 540000
  integer(i4), parameter :: nueg    = 8
  integer(i4), parameter :: n_time  = 300

  ! Timer variables
  integer (i8) :: count_rate, count_max

  ! Some constants
  real(dp), parameter :: pi   = 3.14159265359d0 
  real(dp), parameter :: sr0  = 6.955e10  ! solar radius (cm)
  real(dp), parameter :: nueq = 460.7e-9 ! solar rotation frequency (s^-1)

  ! Angular parameter for outputs
  real(dp), parameter :: theta1 = 0.d0
  real(dp), parameter :: theta2 = 3.1415926d0
  real(dp) :: dt

  ! Arrays that will need to be allocated
  integer(i4), allocatable :: indeg(:)     ! index the eigenvalues      
  real(dp), allocatable :: inde(:)         ! index the eigenvalues      
  real(dp), allocatable :: vr(:,:)
  complex(dp), allocatable :: cvr(:,:)
  complex(dp), allocatable :: ww(:)
  real(dp), allocatable :: wr(:), wi(:)  ! eigenvalues

  complex(dp), allocatable :: ala(:,:)              ! alpha and gamma
  complex(dp), allocatable :: gaa(:,:)
  complex(dp), allocatable :: bea(:,:)              ! array n-dep beta 
  complex(dp), allocatable :: alb(:,:)              ! alpha and gamma
  complex(dp), allocatable :: gab(:,:)
  complex(dp), allocatable :: beb(:,:)               ! array n-dep beta 

  complex(dp), allocatable :: fanp2(:,:)
  complex(dp), allocatable :: fanm2(:,:)
  complex(dp), allocatable :: fanp2p(:,:)
  complex(dp), allocatable :: fanm2p(:,:)
  complex(dp), allocatable :: fbnp2(:,:)
  complex(dp), allocatable :: fbnm2(:,:)
  complex(dp), allocatable :: fbnp2p(:,:)
  complex(dp), allocatable :: fbnm2p(:,:)

  complex(dp), allocatable :: fanp4(:,:)
  complex(dp), allocatable :: fanm4(:,:)
  complex(dp), allocatable :: fanp4p(:,:)
  complex(dp), allocatable :: fanm4p(:,:)
  complex(dp), allocatable :: fbnp4(:,:)
  complex(dp), allocatable :: fbnm4(:,:)
  complex(dp), allocatable :: fbnp4p(:,:)
  complex(dp), allocatable :: fbnm4p(:,:)

  complex(dp), allocatable :: fanp6(:,:)
  complex(dp), allocatable :: fanm6(:,:)
  complex(dp), allocatable :: fanp6p(:,:)
  complex(dp), allocatable :: fanm6p(:,:)
  complex(dp), allocatable :: fbnp6(:,:)
  complex(dp), allocatable :: fbnm6(:,:)
  complex(dp), allocatable :: fbnp6p(:,:)
  complex(dp), allocatable :: fbnm6p(:,:)

  complex(dp), allocatable :: omep1(:,:)
  complex(dp), allocatable :: omem1(:,:)
  complex(dp), allocatable :: omep1p(:,:)
  complex(dp), allocatable :: omem1p(:,:)
  complex(dp), allocatable :: omep3(:,:)
  complex(dp), allocatable :: omem3(:,:)
  complex(dp), allocatable :: omep3p(:,:)
  complex(dp), allocatable :: omem3p(:,:)
  complex(dp), allocatable :: omep5(:,:)
  complex(dp), allocatable :: omem5(:,:)
  complex(dp), allocatable :: omep5p(:,:)
  complex(dp), allocatable :: omem5p(:,:)

  complex(dp), allocatable :: alqm1(:,:)
  complex(dp), allocatable :: alqp1(:,:)
  complex(dp), allocatable :: gaqm1(:,:)
  complex(dp), allocatable :: gaqp1(:,:)
  complex(dp), allocatable :: beqm1(:,:)
  complex(dp), allocatable :: beqp1(:,:)

  complex(dp), allocatable :: alqm3(:,:)
  complex(dp), allocatable :: alqp3(:,:)
  complex(dp), allocatable :: gaqm3(:,:)
  complex(dp), allocatable :: gaqp3(:,:)
  complex(dp), allocatable :: beqm3(:,:)
  complex(dp), allocatable :: beqp3(:,:)
  

  !-------------------------------------------
  !-------------------------------------------
  ! Namelist variables and default values
  !-------------------------------------------
  !-------------------------------------------
  !> global
  !-------------------------------------------
  real(dp) :: sr, rotp

  !-------------------------------------------
  !> grid
  !-------------------------------------------
  integer(i4) :: np    ! mesh points
  integer(i4) :: na    ! na = 1,3,5,... radial order in an
  integer(i4) :: n_theta
  ! nb and nt are not actually in the namelist but depend on 
  ! np and na
  integer(i4) :: nb ! nb = 2,4,6,... radial order in bn
  integer(i4) :: nt ! matrix dimension 
 
  !-------------------------------------------
  !> profiles
  !-------------------------------------------
  character(len=8) :: regime
  real(dp) :: s0, s2
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
  integer(i4) :: nj
  logical :: write_vectors
  real(dp) :: zeta_r
  real(dp) :: xbt, xbo
  logical :: show_timer

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
  integer(i4) :: nso
  real(dp) :: flg
  real(dp) :: cm_f, cm_i, rm_f, rm_i

  ! Additional variables that are used by several 
  ! subroutines. These were included in "common" 
  ! blocks in the old version of the code. 
  integer(i4)  :: ii, it
  real(dp) :: etep, etet, eep
  real(dp) :: co, c_u, beta, betb
  real(dp) :: reg(10),ieg(10)
  real(dp) :: gam
  character(len=1) :: jobvr

  namelist /global/ sr, rotp

  namelist /grid/ np, n_theta, na 

  namelist /profiles/ regime, s0, s2, xa1, xa2, xda1, &
                      xda2, xb, gd, edr, xe1, xde1, &
                      dd1, xc1, c2_h2, oco

  namelist /brent/ al_i, al_f, accu

  namelist /outputs/ dir, nj, write_vectors, zeta_r, &
                     xbt, xbo, show_timer

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
    sr = 1.d0 
    rotp = 1.d0

    !-------------------------------------------
    !> grid
    !-------------------------------------------
    np = 50 
    na = 49
    n_theta = 301
   
    !-------------------------------------------
    !> profiles
    !-------------------------------------------
    regime  = 'h2'
    xa1     =  0.64d0   
    xa2     =  0.72d0   
    xda1    =  0.025d0
    xda2    =  0.025d0
    s0      =  0.85d0
    s2      =  0.08d0
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
    show_timer = .false.
  
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
    !> it. Array depending on grid parameters are also
    !> allocated at this step.
    character(len=128) :: inlist

    integer(i4) :: fu

    call initialise_namelist_values
    call getarg(1, inlist)
    ! Open and close to read in any order
    open(newunit=fu, file=inlist, status="old", action="read")
    read(unit=fu, nml=global)
    close(fu)
    open(newunit=fu, file=inlist, status="old", action="read")
    read(unit=fu, nml=grid)
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

    ! We need to initialise the final grid parameters
    nb = na+1  ! nb = 2,4,6,... radial order in bn and bn
    nt = nb*np ! matrix dimension 
    dt = (theta2 -  theta1) / n_theta
    ! Arrays are then allocated
    allocate (indeg(nt))
    allocate (inde(nt))
    allocate (vr(nt, nt))
    allocate (cvr(nt, nt))
    allocate (ww(nt))
    allocate (wr(nt))
    allocate (wi(nt))

    allocate(ala(np,nb)) ! alpha and gamma
    allocate(gaa(np,nb))
    allocate(bea(np,nb)) ! array n-dep beta 
    allocate(alb(np,nb)) ! alpha and gamma
    allocate(gab(np,nb))
    allocate(beb(np,nb)) ! array n-dep beta 

    allocate(fanp2(np, nb))
    allocate(fanm2(np, nb))
    allocate(fanp2p(np,nb))
    allocate(fanm2p(np,nb))
    allocate(fbnp2(np, nb))
    allocate(fbnm2(np, nb))
    allocate(fbnp2p(np,nb))
    allocate(fbnm2p(np,nb))

    allocate(fanp4(np, nb))
    allocate(fanm4(np, nb))
    allocate(fanp4p(np,nb))
    allocate(fanm4p(np,nb))
    allocate(fbnp4(np, nb))
    allocate(fbnm4(np, nb))
    allocate(fbnp4p(np,nb))
    allocate(fbnm4p(np,nb))

    allocate(fanp6(np, nb))
    allocate(fanm6(np, nb))
    allocate(fanp6p(np,nb))
    allocate(fanm6p(np,nb))
    allocate(fbnp6(np, nb))
    allocate(fbnm6(np, nb))
    allocate(fbnp6p(np,nb))
    allocate(fbnm6p(np,nb))

    allocate(omep1(np, nb))
    allocate(omem1(np, nb))
    allocate(omep1p(np,nb))
    allocate(omem1p(np,nb))
    allocate(omep3(np, nb))
    allocate(omem3(np, nb))
    allocate(omep3p(np,nb))
    allocate(omem3p(np,nb))
    allocate(omep5(np, nb))
    allocate(omem5(np, nb))
    allocate(omep5p(np,nb))
    allocate(omem5p(np,nb))

    allocate(alqm1(np,nb))
    allocate(alqp1(np,nb))
    allocate(gaqm1(np,nb))
    allocate(gaqp1(np,nb))
    allocate(beqm1(np,nb))
    allocate(beqp1(np,nb))

    allocate(alqm3(np,nb))
    allocate(alqp3(np,nb))
    allocate(gaqm3(np,nb))
    allocate(gaqp3(np,nb))
    allocate(beqm3(np,nb))
    allocate(beqp3(np,nb))

    if (show_timer) then
      call system_clock (count_rate=count_rate)
      call system_clock (count_max=count_max)
    endif

  end subroutine
  
  subroutine write_namelist (inlist)
    !> Write full set of namelists variables.
    !> Namelist variables are declared in cio.f90 and
    !> are therefore accessible to any module that import
    !> it.
    character(len=128) :: inlist

    integer(i4) :: fu

    ! Open and close to read in any order
    open(newunit=fu, file=inlist, status="replace", & 
         action="write")
    write(unit=fu, nml=global, delim="quote")
    write(unit=fu, nml=grid, delim="quote")
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


