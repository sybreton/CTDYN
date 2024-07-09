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

program advection
 
  use cio
  implicit none   

  real :: imag
  common/part/vtu,rt,imag,co,c_u,beta,ffree,betb,etep,etet,xbt,xbo
  
  real :: ca(10,10,4)
  real :: reg(10),ieg(10)
  real :: eep
  common/eira/reg,ieg,eep
  integer :: ii,it
  common/ipar/ii,it
  
  real :: xa1,xa2,xa3,xb,xda1,xda2
  common/apar/xa1,xa2,xa3,xb,xda1,xda2
  
  real :: edr,xe1,xde1
  common/epar/edr,xe1,xde1,hd
  
  real :: x_in, bct, c3, mmm,sr,rotp,gd, aqu, flg
  common/ppar/x_in,bct,c3,mmm,sr,rotp, gd,aqu, flg
  
  real :: dd1,rc1,rc2,oco
  common/dpar/dd1,rc1,rc2,oco
  
  real :: s0,s2,s4,s6, a2p, a4p,xm
  common/psi/s0,s2,s4,s6,a2p,a4p,xm
  
  character(len=128) :: dir
  character*8 :: ans1,ans2,ans3,ans4
  common/var3/ans1,ans2,ans3,ans4, dir
  
  character*2 :: jobvr,jobvl
  common/lap/jobvr,jobvl
  
  character*30 :: inp
  character*52 :: mach, ddt
  character*43 :: version, ver
  
  character(len=5)   :: char_m  ! converti m in char
  character(len=128) :: char_fm ! file input_am
  character(len=128) :: char_cm ! file crtal_am
  character(len=128) :: char_cf ! file crtaf_am
  character(len=128) :: fom     ! crea un format dinamico 
  character(len=5)   :: qq      ! _am
  external version
  integer :: mm
  
  integer :: iome,nso
  integer :: ialo,nsa
  integer :: ires,nsr
  
  common/parker/gam,zeta_r,ratio
  
  include 'cint.f90'
  
  call getarg(1,inp)
  open(31,file=inp,status='old')
  read(31,namp)
  close(31)
  !converti mm in char e crea stringa _am o _sm
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
    call zbr(al_i,al_f,accu,critical) 
    write (35,'(i4,16e12.4)') ii, critical,co,c_u,beta,etep,etet,zeta_r,(reg(i),abs(ieg(i)),i=1,1), eep
    if (ans4 .eq. 'v') then
      jobvr='v'
      call dynamo(critical,rate)
      jobvr='n'
    endif
  enddo
  
  close(35)

end

