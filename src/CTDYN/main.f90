!     Linear dynamo code CTDYN 
!                     double precision/change intrinsic es: erf -> derf 
!
!     one cell   s0 = -1/2    s2 =  3/2     s4=0    s6=0  
!     two cells  s0 =  3/8    s2 = -15/4    s4=35/8 
!
!
!   The induction equation is written dt B = rot ( v x B "+" alpha B) -rot eta rot B
!
!   Therefore the turbulent elf is e = + alpha B
!   
!   This alpha is "minus" the Raedler an 307 89 1986 paper (r86)
!   
!   The B field is written as  
!
!   B \propto exp (  \nu \tau)
!   
!   and in the R86 paper instead we have 
!   
!   B \propto exp (im \phi + (\lambda - i \omega) \tau )      
!
!   therefore im(nu) =  - \omega
!
!  as a test example you have for  n = 14 m = 30
!
!    a0 comega = 20      -> calpha = 4.51   \nu = 0
!    a1 comega = 20      -> calpha = 5.01   \nu = 7.43 = - \omega r86
!    

program main
 
  use kind_parameters
  use cio
  use dyna
  use zbrent

  implicit none   

  real(dp) :: astep, critical, rate 
  real(dp) :: ca(10,10,4)

  character(len=128) :: inlist, inlist_full
  character(len=5)   :: char_m         ! convert m in char
  character(len=128) :: fileout        ! file crtal_am
  character(len=5)   :: qq             ! _am

  integer :: mm, i, iome
  integer :: fu 


  
  ! Read input namelist
  call read_namelist (inlist) 

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

  ! Setting global variables
  ii = 0
  jobvr = 'n'
  fileout = trim(trim(dir)//'/critical'//qq)

  ! Entering the main loop calling the bisection
  ! for different set of rotation/circulation
  ! regimes
  open(newunit=fu, status='unknown', file=fileout)
  write(fu, '(a,i4,a,i4,a)') '#  n, C_alpha, C_omega, r, beta, Etor, Epol, z = (', np, ',', na+1, ')'  
  do iome=0, nso
    co = cm_i + iome*(cm_f-cm_i)/(nso+1) 
    c_u = rm_i+rm_f*co**xm
    ii = ii+1    
    call zbr(al_i, al_f, accu, critical) 
    write (fu,'(i4,16e12.4)') ii, critical, co, c_u, beta, etep, etet, zeta_r, (reg(i), abs(ieg(i)),i=1,1), eep
    if (write_vectors) then
      jobvr='v'
      call dynamo (critical, rate)
      jobvr='n'
    endif
  enddo
  close(fu)

end program main

