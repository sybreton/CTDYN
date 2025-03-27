!     Linear dynamo code CTDYN 
!                     double precision/change intrinsic es: erf -> derf 
!
!   The induction equation is written dt B = rot ( v x B "+" alpha B) -rot eta rot B
!
!   Therefore the turbulent elf is e = + alpha B
!   
!   This alpha is "minus" the Raedler an 307 89 1986 paper (r86)
!   
!   The B field is written as  
!
!   B \propto exp ( \nu \tau)
!   
!   and in the R86 paper instead we have 
!   
!   B \propto exp (im \phi + (\lambda - i \omega) \tau )      
!
!   therefore im(nu) =  - \omega
!

program main
 
  use kind_parameters
  use cio
  use dyna
  use zbrent

  implicit none   

  real(dp) :: critical ! critical coefficient for alpha effect
  real(dp) :: eta      ! turbulent diffusivity coefficient
  real(dp) :: period   ! cycle period in year
  real(dp) :: rate, imag ! real and imaginary part of the selected eigenvalue

  real(dp) :: ca(10,10,4)

  character(len=128) :: inlist, inlist_full
  character(len=5)   :: char_m         ! convert m in char
  character(len=128) :: fileout        ! file crtal_am
  character(len=5)   :: qq             ! _am

  integer(i4) :: mm, i, iome
  integer(i4) :: fu 

  ! Timer variables
  integer (i8) :: count_0, count_1
  
  ! Read input namelist
  call read_namelist (inlist) 

  if (show_timer) then
    call system_clock (count_0)
  endif

  ! Write the full namelist for reference
  inlist_full = trim(trim(dir)//'/inlist_full')
  call write_namelist (inlist_full) 

  ! Convert mm in char and create a string _am or _sm
  mm=int(mmm)
  write(char_m,'(i2)')  mm
  if (mod(mm,2) .eq. 0) then 
     if (degree.eq. 'd') then
       qq=trim('_a'//adjustl(char_m))
     else
       qq=trim('_s'//adjustl(char_m))
     endif
  else if (mod(mm,2) .eq. 1) then
     if(degree .eq. 'd' ) then
       qq=trim('_s'//adjustl(char_m))
     else
       qq=trim('_a'//adjustl(char_m))
     endif
  endif
  
  ! Main ff-beta parameter
  beta = beta_i
  eep = -99

  ! Setting global variables for the loop
  jobvr = 'n'
  fileout = trim(trim(dir)//"/critical"//qq)//".dat"

  ! Entering the main loop calling the bisection
  ! for different set of rotation/circulation regimes
  open(newunit=fu, status='unknown', file=fileout)
  write(fu, '(a)') "# C_alpha, C_omega, C_meridional, omega_cycle, period_cycle, eta, beta, Etor, Epol"  
  write(fu, '(a)') "# (adim), (adim), (adim), (adim), (yr), (cm^2/s), (-), (-), (-)"  
  c_u = rm_i+rm_f*co**xm
  call zbr(al_i, al_f, accu, critical) 
  if (write_vectors) then
    jobvr='v'
    call dynamo (critical, rate, imag, eta, period)
  endif
  write (fu,'(9es12.4)') critical, co, c_u, abs(imag), abs(period), eta, beta, etep, etet
  close(fu)

  if (show_timer) then
    call system_clock (count_1)
    write (*,'(a,1x,f0.3,1x,a)') "Total elapsed time::", real(count_1-count_0)/count_rate, "seconds."
  endif

end program main

