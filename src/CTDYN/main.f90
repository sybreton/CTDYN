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
 
  use cio
  use dyna
  use zbrent

  implicit none   

  real(dp) :: astep, beta_f, beta_i, &
          & beta_s, cm_f, cm_i, &
          & critical, rate, rm_f, rm_i 
  
  real(dp) :: ca(10,10,4)
  real(dp) :: accu
  real(dp) :: al_i, al_f

  character*30 :: inp
  character*52 :: mach, ddt
  character*43 :: version, ver
  external version
  
  character(len=5)   :: char_m  ! converti m in char
  character(len=128) :: char_fm ! file input_am
  character(len=128) :: char_cm ! file crtal_am
  character(len=128) :: char_cf ! file crtaf_am
  character(len=128) :: fom     ! crea un format dinamico 
  character(len=5)   :: qq      ! _am
  integer :: mm
  
  integer :: i
  integer :: iome, nso
  integer :: ialo, nsa
  integer :: ires, nsr
  
  include 'cint.f90'
  
  call getarg(1, inp)
  open(31, file=inp, status='old')
  read(31, namp)
  close(31)
  !convert mm in char and create a string _am or _sm
  mm=int(mmm)
  write(char_m,'(i2)')  mm
  if (mod(mm,2) .eq. 0) then 
     if (ans2.eq. 'd') then
       qq=trim('_a'//adjustl(char_m))
     else
       qq=trim('_s'//adjustl(char_m))
     endif
  else if (mod(mm,2) .eq. 1) then
     if(ans2 .eq. 'd' ) then
       qq=trim('_s'//adjustl(char_m))
     else
       qq=trim('_a'//adjustl(char_m))
     endif
  endif
  
  char_fm=trim(trim(dir)//'/input'//qq)
  char_cm=trim(trim(dir)//'/crtal'//qq)
  char_cf=trim(trim(dir)//'/crtaf'//qq)
  
  open(34,status='unknown',file=char_fm)
  write(34,namp)
  close(34) 
  open(35,status='unknown',file=char_cm)
  
  !main ff-beta parameter
  beta=beta_i
  eep=-99
  ii=0
  
  jobvr='n'
  ans3='cr'
  
  write(35,'(a,2i4)')  '#n c_alpha,c_omega,r,beta,etor,epol,z=', np, na+1 
  
  do iome=0, nso
    if ( cm_i.ne.0 ) then
      astep=(abs(cm_f/cm_i))**(1./nso) 
      co = astep**iome*cm_i
    else
      co = cm_i + iome*(cm_f-cm_i)/(nso+1) 
    endif
    c_u=rm_i+rm_f*co**xm
    ii=ii+1    
    call zbr(al_i, al_f, accu, critical) 
    write (35,'(i4,16e12.4)') ii, critical, co, c_u, beta, etep, etet, zeta_r, (reg(i), abs(ieg(i)),i=1,1), eep
    if (ans4 .eq. 'v') then
      jobvr='v'
      call dynamo(critical,rate)
      jobvr='n'
    endif
  enddo
  
  close(35)

end

