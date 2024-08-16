!
!   Radler73 Zur DynamoTheorie / Notes 
!
!   the general structure of the equation is:
!
!   dPhi/dt = Gn + Dn + Nabla Phi  3.17a
!   dPsi/dt = Fn + Cn + Nabla Psi  3.17b
!   
!   Phi - a-coefficients 
!   Psi - b-coefficients
!
!   Gn : alpha-omega  m =1   
!   Dn : alpha B term: written as: alpha * cos * ( 1 - c3 * cos^2)  if c3=0,  standard , if c3=1 we have cos * sin^2
!   c3=1 only correct for m=0 !!!!!!!! 
!   note that alpha = -alpha1 of Rad73
!
!   General structure of the type:
!
!  a Phi'' + b Phi' - c Phi ....
!
!  a ( Phi_i+1 -2 Phi_i +Phi_i-1)/hh**2 + b (Phi_i+1 -Phi_i-1)/(2 hh) - c Phi ....
!
!  (a+h2 b)  Phi_i+1  - (2 a + hh**2 c) Phi_i +(a-h2 b) Phi_i-1 
!
!  alpha_i Phi_i-1  - Beta_i Phi_i +gamma_i Phi_i+1
!
!  where h2 = hh/2
!
!
!-------------------------------------------------------------
subroutine dynamo(turb,rate)
!-------------------------------------------------------------

  use cio
  use cvect
  use ccov
  use b2func
  use b3func
  use b4func
  use b5func
  use c1s2func
  use c3s2func
  use cdfunc
  use bessel
  use write_outputs

  implicit none
  
  real :: imag, turb, rate, a2, alp, alp_p, &
          & am, apu, apup, ar, arp, bapu, bapup, &
          & beta, betb, btp, c_u, ca, co, ct, &
          & deltami, deltapi, e1, e2, e3, etep, &
          & etet, facp, ffc, ffree, gam, h2, hh, &
          & om0, om0p, om2, om2p, om4, om4p, psip, &
          & psipp, ratio, rho, rt, sigmami, sigmapi, &
          & um, vtu, x, x1, x2, xbo, xbt, zeta_r, &
          & bt, nabp1, nabp3, psi
  common/part/vtu,rt,imag,co,c_u,beta,ffree,betb,etep,etet,xbt,xbo
  real :: eep
  real :: reg(10),ieg(10)
  common/eira/reg,ieg,eep
  
  !      definitions:
  !      turb= c_alpha*c_omega
  !      c_alpha = alpha_zero*sr/eta
  !      c_omega = omega_zero*sr**2/eta
  !      common param
  
  integer :: i, j, k, i2, i3, nabm1, &
             & k1, ki, kk, nabm3, napb3, &
             & naf, nafm2, nafp2, nai, nas, &
             & nas_count, nbam1, nbam3, nbap1, &
             & nbap3, nbf, nbfm2, nbfp2, &
             & nbi, nc                
  integer :: info
  integer :: ii, it ! loop counter
  common/ipar/ii,it
  
  real :: xa1,xa2,xa3,xb,xda1,xda2
  common/apar/xa1,xa2,xa3,xb,xda1,xda2
  
  real :: edr,xe1,xde1,hd
  common/epar/edr,xe1,xde1,hd
  
  real :: x_in, bct, c3,mmm,sr,rotp, gd, aqu, flg
  common/ppar/x_in,bct,c3,mmm,sr,rotp,gd, aqu, flg
  
  real :: s0,s2,s4,s6, a2p, a4p,xm
  common/psi/s0,s2,s4,s6,a2p,a4p,xm 
  
  character(len=128) :: dir
  character*8 :: ans1,ans2,ans3,ans4
  common/var3/ans1,ans2,ans3,ans4, dir
  
  character*2 :: jobvr,jobvl
  common/lap/jobvr,jobvl
  
  external version
  
  !------------- temporary variables
  
  real :: omegaeq
  real :: nueq, eta, sr0
  parameter(sr0=6.955e10,nueq=460.7e-9)
  
  real :: egr(np,nb), egi(np,nb)
  real :: agr(np,nb/2), agi(np,nb/2)
  real :: bgr(np,nb/2), bgi(np,nb/2)
  
  character*1 :: jobpr    
  character*1 :: balanc    
  character*1 :: sense    
  real :: rwork(2*nt)
  complex*16 :: work(lwork)  
  character*43 :: ver
  integer :: qa, qb
  complex*16 :: axx(nt,nt)
  real :: axr(nt,nt)
  complex*16 :: drb, dra, sigmeno, sigpiu
  common/parker/gam,zeta_r,ratio
  common/vec/cvr
  
  real :: c1(np)                  ! coefficient dyn.  
  real :: c2(np)
  
  !------------------------------------------------------
  jobpr = 'n'    ! y = write all pij matrix       !!!!
  jobvl = 'n'    ! v = calculate r and l eig.
  !----------------  boundaries ------------------------------
  
  x1 = x_in           !inner boundary
  x2 = 1.0            !outer boundary
  hh =(x2-x1)/(np+1)  !stepsize: radial accuracy parm.
  h2 = hh/2.e0
  a2 = c3       ! cos^3 term in alpha  
  
  ! rotp = periodo di rotazione in unita' di quello solare
  
  omegaeq = nueq*2.e0*pi/rotp
  
  if(co.ne.0)then            
    eta = (sr*sr0)**2*omegaeq/co
  else
    eta=1
  endif
  
  !-------------------------------------------------------------
  
  um = c_u*eta/(sr*sr0)
  
  !-------------------------------------------------------------
  
  ca = turb             
  
  !-------------------------------------------------------------
  !
  !              Write iterations on screen
  !
  !-------------------------------------------------------------
  
  ! In this case the coefficients are not computed....
  if(c3.ne.0.and.mmm.ne.0)then
    print*, 'not available'
    stop
  endif 
  
  if(ans2.eq.'d'.and.mmm.eq.0) write(*,*) '********** a0-mode *********'
  if(ans2.eq.'d'.and.mmm.ne.0) write(*,*) '********** s1-mode *********'
  if(ans2.eq.'q'.and.mmm.eq.0) write(*,*) '********** s0-mode *********'
  if(ans2.eq.'q'.and.mmm.ne.0) write(*,*) '********** a1-mode *********'
  
  write(*,'(1x,a,i4,1x,a,i4,1x,a,i4)') ' it=', it,' na=', nb, ' np=', np
  write(*,'(a,e12.5,a,e12.5)' ) ' c_alpha=', turb, ' c_omega=', co 
  write(*,'(a,e12.5,a,e12.5)' ) ' r_flow =', c_u,  ' f-f=    ', beta 
  
  
  !
  !-------------------------------------------------------------
  !     Write alpha and eta the first iteration/test run
  !-------------------------------------------------------------
  !
  
  if(ii.eq.1)then     
    open(22,status='unknown',file=trim(trim(dir)//'/alpha.dat'))
    do i=1, np
      x = x1+i*hh
      call rot(x, om0, om0p, om2, om2p, om4, om4p) 
      call alt(x, alp, alp_p)     
      call eta_turb(x, e1, e2, e3)
      call stream(x,ar,arp,bt,btp, psi)
      write(22,'(15f14.6)') x,alp,alp_p,e1,e2,e3,psi,ar,arp,bt,btp,om0,om0p,om2,om2p
    enddo
    close(22)
  endif
  
  !----------------- Load arrays in m and n
  
  nc                =  2 ! number of latitudinal cells /2 set to 2 by def.
  if (a4p.ne.0) nc  =  4 ! if a4p not zero there are two cells-> nc=4
  if (s6.ne.0)  nc  =  6       
  !
  ! These are the conversion factors for p1n for utheta
  !
  !      s4 = (35.e0/8.e0) * a4p
  !      s2 = (3.e0/2.e0) * a2p -(15.e0/4.e0)*a4p
  !
  !      ss2 = (9.e0/(nc*(nc+1.e0)))*a2p-(75.e0/(nc*(nc+1.e0)))*a4p
  !      s0= 15.e0/(2*(nc*(nc+1.e0)))*a4p-(3.e0/(nc*(nc+1.e0)))*a2p
  !
   
  am = abs(mmm)
  
  do j=1, nb
    x=0.0
    do i=1, np
      x=x1+i*hh
      
      call rot(x, om0, om0p, om2, om2p, om4, om4p)
      call alt(x, alp, alp_p)
      call eta_turb(x, e1, e2, e3)
      call stream(x,ar,arp, bt,btp, psi)
      
      rho=x**gd 
      rho=0
      psi=0
      psip=0
      psipp=0
      
      c1(i) = ca*alp*hh**2
      c2(i) = c_u*ar*hh**2
      
      betb=beta
      
      !  alpha beta gamma of the notes for a-vectors
      !  in the matrix we write  
      !  
      !  alpha_i b_i-1  - beta_i b_i  + gamma_i b_i+1    
      !
      !  -> alpha_i , beta_i, gamma_i  .... = axx
      ! 
      !  beta_i = 2+h^2 n(n+1)/x^2  + h^2 c_u f ...
      !
      !
      dra = h2*((-2*e1/x - hd * e2/2) &
           & +ca*alp*mmm*(0,1)/j/(j+1.) &  ! check sign
           & -c_u*((j*(j+1)-3*am**2)/j/(j+1.)/(2*j-1.)/(2*j+3.)) &
           & *(2*(j*(j+1)-3)*ar)) 
      
      bea(i,j) =2*e1+hh**2*((e1*j*(j+1)/x**2 -hd*e2/x  ) & 
           &  +ca*alp*mmm*(0,1)/j/(j+1.)/x &  ! check sign
           & +(co)*(+om0*mmm*(0,1) + (+1)*om2*mmm*(0,1)*(j*(j+1)-3*mmm**2)/(2*j-1)/(2*j-3)) & 
           & -c_u*((j*(j+1)-3*am**2)/j/(j+1.)/(2*j-1.)/(2*j+3.)) &
           & *(2*(j*(j+1)-3)*ar/x-j*(j+1)*bt/x)) 
      
      gaa(i,j) =  e1-dra
      ala(i,j) =  e1+dra
      
      ! set outer boundary conditions for a-vectors
      
      !        bea(np,j) = bea(np,j) - 4*gaa(np,j)/(2*hh*(j+1)+3)
      !        ala(np,j) = ala(np,j) -  gaa(np,j)/(2*hh*(j+1)+3)     
      
      !  alpha beta gamma of the notes for b-vectors
      
      beb(i,j) = e1*2.e+00 + hh**2*((e1*j*(j+1.e+00)/x**2-(e2/x)-hd*e2/x/2 -hd*e3/2  )  &
           & +ca*(alp/x+alp_p)*mmm*(0,1)/j/(j+1.e+00) & !check sign
           & + (co)*(+om0*mmm*(0,1) + (+1)*om2*(0,1)*mmm*(j*(j+1.e+00)-3.e+00*mmm**2)/(2*j-1.e+00)/(2*j+3.e+00)) &
           & -c_u*((j*(j+1)-3.*am**2)/j/(j+1)/(2*j-1.)/(2*j+3.)) &
           & *(2*(j*(j+1)-3)*(arp+ar/x)+bt*j*(j+1)/x)) 
      
      drb =  h2*((-2.*e1/x-e2 -hd*e2/2.) &
           &  +ca*(alp)*mmm*(0,1)/j/(j+1.) &  !check sign
      ! check sign c_u
           & -c_u*((j*(j+1)-3.*am**2)/j/(j+1)/(2*j-1.)/(2*j+3.))*(2*(j*(j+1)-3)*ar)) 
      
      gab(i,j) = e1 -  drb
      alb(i,j) = e1 +  drb
      
      !--------------------  a-2 a+2 ------------------------------------------
      
      fanm2(i,j) = -c_u*hh**2*((j-2.)*(j-am-1.)*(j-am)/j/(2*j-3.)/(2*j-1)) &
           & * (3*ar/x-(j-1.)*bt/x)
        
      fanm2p(i,j) = -c_u*h2*((j-2.)*(j-am-1.)*(j-am)/j/(2*j-3.)/(2*j-1))*3*ar  
      
      fanp2(i,j) = -c_u*hh**2*((j+3.)*(j+am+1)*(j+am+2)/(j+1)/(2*j+3.)/(2*j+5)) &
           & * (3*ar/x+(j+2)*bt/x)
      
      fanp2p(i,j) = -c_u*h2*((j+3.)*(j+am+1)*(j+am+2)/(j+1)/(2*j+3.)/(2*j+5.))*3*ar
      
      !---------------------- a-4 a+4
      
      fanm4(i,j)  = 0
      fanp4(i,j)  = 0
      fanm4p(i,j) = 0
      fanp4p(i,j) = 0
      fanm6(i,j)  = 0
      fanp6(i,j)  = 0
      fanm6p(i,j) = 0
      fanp6p(i,j) = 0
      
      !--------------------------------- b-2 b+2---------------------------------
      
      fbnm2(i,j) = -c_u*hh**2* &
           & ((j-2.)*(j-am-1.)*(j-am)/j/(2.*j-3.)/(2.*j-1)) &
           & *(3.*(arp+ar/x)-j*bt/x)
      
      fbnm2p(i,j) = -c_u*h2* &
           & ((j-2.)*(j-am-1.)*(j-am)/j/(2.*j-3.)/(2.*j-1))*3.*ar
      
      fbnp2(i,j) = -c_u*hh**2* &
           & ((j+3.)*(j+am+1.)*(j+am+2)/(j+1.)/(2.*j+3.)/(2.*j+5.)) &
           & *(3.*(arp+ar/x)+(j+1.)*bt/x)
      
      fbnp2p(i,j) = -c_u*h2* &
           & ((j+3.)*(j+am+1.)*(j+am+2)/(j+1.)/(2*j+3.)/(2*j+5.))*3.*ar
      
      fbnm4(i,j)  = 0
      fbnp4(i,j)  = 0
      fbnm4p(i,j) = 0
      fbnp4p(i,j) = 0
      fbnm6(i,j)  = 0
      fbnp6(i,j)  = 0
      fbnm6p(i,j) = 0
      fbnp6p(i,j) = 0
      
      ! --------------- omega arrays ---------------------------------------
      
      if(ans1.eq.'h4')then
        omep1(i,j) = co*hh**2* (b0m1(j+1)*om0p+b2m1(j+1)*om2p &
           &+b4m1(j+1)*om4p+2*c1s2m1(j+1)*om2/x+4*c3s2m1(j+1)*om4/x)
      
        omem1(i,j) = co*hh**2*(b0p1(j-1)*om0p+b2p1(j-1)*om2p &
           &+b4p1(j-1)*om4p+2*c1s2p1(j-1)*om2/x+4*c3s2p1(j-1)*om4/x)
      
        omep3(i,j)= co*hh**2* &
           &       (b2m3(j+3)*om2p+2*om2*c1s2m3(j+3)/x+4*om4*c3s2m3(j+3)/x)
      
        omem3(i,j)= co*hh**2* &
           &       (b2p3(j-3)*om2p+2*om2*c1s2p3(j-3)/x+4*om4*c3s2p3(j-3)/x)
      
        omep1p(i,j)= co*h2*(2*om2*c1s2m1(j+1)+4*om4*c3s2m1(j+1))
        omem1p(i,j)= co*h2*(2*om2*c1s2p1(j-1)+4*om4*c3s2p1(j-1))
        omep3p(i,j)= co*h2*(2*om2*c1s2m3(j+3)+4*om4*c3s2m3(j+3))
        omem3p(i,j)= co*h2*(2*om2*c1s2p3(j-3)+4*om4*c3s2p3(j-3))
      
        omep5(i,j)= co*hh**2*(4*om4*c3s2m5(j+5)/x)
        omem5(i,j)= co*hh**2*(4*om4*c3s2p5(j-5)/x)
      
        omep5p(i,j)= co*h2*(4*om4*c3s2m5(j+5))
        omem5p(i,j)= co*h2*(4*om4*c3s2p5(j-5))
  
      else
        omem1(i,j) = co*hh**2*( &
             &  om0p*(j-1)*(j-abs(mmm))/(2*j-1) &
             &  +(j-abs(mmm))*(6*j*((j-1)*(j+1)-5*mmm**2)*om2/x - (j-1)*((j-3)*j*(j+1)+3*(j+6)*mmm**2)*om2p) &
             & /2/(j+1)/(2*j-3)/(2*j-1)/(2*j+3) )
        
        omep1(i,j) = co*hh**2*( &
             &  -om0p*(j+2)*(j+abs(mmm)+1)/(2*j+3) &
             &  +(abs(j+mmm)+1)*(6*(j+1)*(j*(j+2)-5*mmm**2)*om2/x+(j+2)*(j*(j+1)*(j+4)+3*(j-5)*mmm**2)*om2p) &
             &  /2/j/(2*j-1)/(2*j+3)/(2*j+5))
        
        omem3(i,j) = co*hh**2*( &
             & -3*(j-3)*(j-abs(mmm)-2)*(j-abs(mmm)-1)*(j-abs(mmm))*(om2/x-(j-2)*om2p/2) &
             &  /j/(2*j-5)/(2*j-3)/(2*j-1))      
         
          
        omep3(i,j) = co*hh**2*( &
             & -3*(j+4)*(j+abs(mmm)+1)*(j+abs(mmm)+2)*(j+abs(mmm)+3)*(om2/x+(j+3)*om2p/2) &
             & /(j+1)/(2*j+3)/(2*j+5)/(2*j+7))
             
        omem1p(i,j) = co*h2*((j-abs(mmm))*(6*j*((j-1)*(j+1)-5*mmm**2)*om2) &
             & /2/(j+1)/(2*j-3)/(2*j-1)/(2*j+3) )
           
        omep1p(i,j) = co*h2*((abs(j+mmm)+1)*(6*(j+1)*(j*(j+2)-5*mmm**2)*om2) &
             &  /2/j/(2*j-1)/(2*j+3)/(2*j+5))
        
        omem3p(i,j) = co*h2*( &
             & -3*(j-3)*(j-abs(mmm)-2)*(j-abs(mmm)-1)*(j-abs(mmm))*(om2) &
             &  /j/(2*j-5)/(2*j-3)/(2*j-1))
        
        omep3p(i,j) = co*h2*( &
             & -3*(j+4)*(j+abs(mmm)+1)*(j+abs(mmm)+2)*(j+abs(mmm)+3)*(om2) &
             &  /(j+1)/(2*j+3)/(2*j+5)/(2*j+7))
      endif
      
      !     check the relative sign: in the old version "+" because they are imput with "-" in the matrix!!!!
      !     check omega hh**2    
      !
      !     4.13a pnm2   
      fbnm2(i,j) = fbnm2(i,j) + 3*(0,1)*mmm*(j-am-1)*(j-am)*hh**2*om2/2./(2*j-3.)/(2*j-1.)
      !     4.13a pnp2
      fbnp2(i,j) = fbnp2(i,j) + 3*(0,1)*mmm*(j+am+1)*(j+am+2)*hh**2*om2/2./(2*j+3.)/(2*j+5.)
      fanm2(i,j) = fanm2(i,j) + 3*(0,1)*mmm*(j-2)*(j-1)*(j-am-1)*(j-am)*hh**2*om2/2./j/(j+1)/(2.*j-3.)/(2*j-1.)
      fanp2(i,j) = fanp2(i,j) + 3*(0,1)*mmm*(j+2)*(j+3)*(j+am+1)*(j+am+2)*hh**2*om2/2./j/(j+1)/(2.*j+3.)/(2*j+5.)
      
      ! ------------------- alpha^2 -----------------------------------
      
      
      ! include the full alpha^2 term for aqu=1
      
      ct = aqu*ca   ! org
      
      alqm1(i,j)=-ct*((alp-h2*(2*alp/x+alp_p))*((j-1.)*(j-abs(mmm))/(2*j-1.)/j &
           & -c3*(j**2-3.)*(j-1.)*3/(2*j-3.)/(2*j+3.)/(2*j-1.)))
      
      alqp1(i,j)=-ct*((alp-h2*(2*alp/x+alp_p))*((j+2.)*(j+abs(mmm)+1)/(2*j+3.)/(j+1.) &
           & -c3*(j**2+2*j-2)*(j+2)*3/(2*j+5.)/(2*j+3.)/(2*j-1.)))
      
      gaqm1(i,j)=-ct*((alp+h2*(2*alp/x+alp_p))*((j-1.)*(j-abs(mmm))/(2*j-1.)/j &
           & -c3*(j**2-3.)*(j-1.)*3/(2*j-3.)/(2*j+3.)/(2*j-1.))) 
      
      gaqp1(i,j)=-ct*((alp+h2*(2*alp/x+alp_p))*((j+2.)*(j+abs(mmm)+1)/(2*j+3.)/(j+1) &
           & -c3*(j**2+2*j-2)*(j+2)*3/(2*j+5.)/(2*j+3.)/(2*j-1.)))
      
      beqm1(i,j)=-ct*(-(2*alp+hh**2*(alp*j*(j-1)/x**2-alp_p/x))*((j-1.)*(j-abs(mmm))/(2*j-1.)/j &
           & +(-c3)* (j**2-3)*(j-1.)*3/(2*j-3.)/(2*j+3.)/(2*j-1.)) &
           & -hh**2*(alp/x/x)*((j-1.)*(j-abs(mmm))/(2*j-1.) &
           & -3*c3*((j**2+j-3)*j*(j-1.)/(2*j-3.)/(2*j+3.)/(2*j-1)))) 
      
      beqp1(i,j)=-ct*(-(2*alp+hh**2*(alp*(j+1)*(j+2.)/x**2-alp_p/x)) &
           & *((j+2.)*(j+abs(mmm)+1)/(2*j+3.)/(j+1.) &
           & +(-c3)*(j**2+2*j-2)*(j+2)*3/(2*j+5.)/(2*j+3.)/(2*j-1.)) &
           & -hh**2*(alp/x/x)*(-(j+2.)*(j+abs(mmm)+1)/(2*j+3.) &
           & -3*c3*(-(j**2+j-3)*(j+1.)*(j+2)/(2*j+5.)/(2*j-1.)/(2*j+3.))))   
      
      alqm3(i,j) = &
           & -ct*(-c3*(alp-h2*(2*alp/x+alp_p)) &
           & *((j-1)*(j-2)*(j-3)/(2*j-5.)/(2*j-3.)/(2*j-1.)))
       
      alqp3(i,j) = &
           & -ct*(-c3*(alp-h2*(2*alp/x+alp_p)) &
           & *((j+2)*(j+3)*(j+4)/(2*j+3.)/(2*j+5.)/(2*j+7.)))
      
      gaqm3(i,j) = &
           & -ct*(-c3*(alp+h2*(2*alp/x+alp_p)) &
           & *((j-1)*(j-2)*(j-3)/(2*j-5.)/(2*j-3.)/(2*j-1.)))
      
      gaqp3(i,j) = &
           & -ct*(-c3*(alp+h2*(2*alp/x+alp_p)) &
           & *((j+2)*(j+3)*(j+4)/(2*j+3.)/(2*j+5.)/(2*j+7.)))
      
      beqm3(i,j) = &
           & -ct*(-(2*alp+hh**2*(j*(j-1)/x**2-alp_p/x))*( &
           & -c3)*((j-1)*(j-2)*(j-3)/(2*j-5.)/(2*j-3.)/(2*j-1.))) &
           & -ct*(hh**2*(alp/x/x)*3*c3*(j-2)*(j-1) &
           & *(j-3.)*(j-2.)/(2*j-3.)/(2*j-1.)/(2*j-5.))
      
      beqp3(i,j) = &
           & -ct*(-(2*alp+hh**2*((j+1)*(j+2.)/x**2-alp_p/x))*( &
           & -c3)*((j+2)*(j+3)*(j+4)/(2*j+3.)/(2*j+5.)/(2*j+7.))) &
           & -ct*(hh**2*(alp/x/x)*3* &
           & c3*(-(j+3)*(j+2)*(j+3.)*(j+4.)/(2*j+3.)/(2*j+5.)/(2*j+7)))
      
      !
      ! define sigmeno complex !!!! ..... otherwise it will be wrong.....!!!!
      !
      
      apu = ar
      apup =  arp
      bapu =  bt
      bapup  = btp
      
      !  check sign with the m-independent terms 
      sigmeno  = -c_u*(0,1)*mmm*(j-am)/j/(j+1.)/(2*j-1.)
      sigpiu   = -c_u*(0,1)*mmm*(j+am+1)/j/(j+1.)/(2*j+3.)
      
      sigmami  =  12*apu/x+6*apup-2*j*bapu/x  ! ok
      
      ! coefficient of the non-derivative term changed of sign - the overall - goes in the definition of "beta" below -
      deltami  =  6*apu*j*(j+1)/x/x - 6*apup/x + 2*j*bapu/x/x-(j-1)* &
           & j*(6*apu/x/x+bapup/x-bapu/x/x)   !ok 
      
      sigmapi  = 12*apu/x+6*apup+2*(j+1)*bapu/x  !ok
      
      deltapi  = 6*apu*j*(j+1)/x/x - 6*apup/x - 2*(j+1)*bapu/x/x-(j+1)* &
           & (j+2)*(6*apu/x/x+bapup/x-bapu/x/x)  ! ok
      
      alqm1(i,j) = alqm1(i,j) + sigmeno*(6*apu-h2*sigmami)
      beqm1(i,j) = beqm1(i,j) - sigmeno*(2*6*apu+hh**2*deltami)
      gaqm1(i,j) = gaqm1(i,j) + sigmeno*(6*apu+h2*sigmami)
      
      alqp1(i,j) = alqp1(i,j) + sigpiu*(6*apu-h2*sigmapi)
      beqp1(i,j) = beqp1(i,j) - sigpiu*(2*6*apu+hh**2*deltapi)
      gaqp1(i,j) = gaqp1(i,j) + sigpiu*(6*apu+h2*sigmapi)
    enddo 
    
    ! outer boundary conditions
    ! you now have two helmholtz equations with different parameters
    ! 
    
    call bess(beta, j, ffc)
    bea(np,j) = bea(np,j) - 4*gaa(np,j)/(2*hh*(j+1+ffc )+3)
    ala(np,j) = ala(np,j) -  gaa(np,j)/(2*hh*(j+1+ffc )+3)     
    
    ! set outer boundary conditions for B-vectors
    if(beta.ne.0)then 
      call bess(betb, j, ffc)
    if(betb.gt.10) ffc=betb
      beb(np,j) = beb(np,j) - 4*gab(np,j)/(2*hh*(j+1+ffc)+3)
      alb(np,j) = alb(np,j) -  gab(np,j)/(2*hh*(j+1+ffc)+3)     
    endif
  enddo
  
  !     initialize the entries in the matrix 
  axx = 0 
  
  ! -------------------------------------------
  
  if (ans2 .eq. 'd') then       ! vector is a1 b2 a3 b4 ...
  
    naf    = na
    nai    = 0
    nbf    = nb
    nbi    = 1
    qa     = 1
    qb     = 2
    
    nabp1  = na
    nabm1  = 3
    nabp3  = na-2
    nabm3  = 5
    nafm2  = 3
    nafp2  = na-2
    
    nbam1 =  2 
    nbap1 =  nb-2
    nbam3 =  4
    nbap3 =  nb-4
    nbfm2 =  4
    nbfp2 =  nb-2
  
  else if ( ans2.eq.'q') then ! vector is b1 a2 b3 a4 ...
  
    naf    =  nb
    nai    =  1 
    nbf    =  na          
    nbi    =  0
    qa     =  2
    qb     =  1 
    
    nabp1  =  na-1
    nabm1  =  2
    nabp3  =  na-3
    nabm3  =  4
    nafm2  =  4
    nafp2  =  na-1
    
    nbam1 =  3
    nbap1 =  nb-1
    nbam3 =  5
    nbap3 =  nb-3
    nbfm2 =  3
    nbfp2 =  na-2
  
  else
    write(*,*) 'not allowed', ans2
    stop
  endif
  
  !---- load matrix entries 
  ! runs col 2*nb  ! inner 
  
  do j = 1, 3*nb, nb
    nas = 0
    nas_count = 0
    do k1 = nai, naf, 2   
      nas_count = nas_count+1
      nas = nas_count-1
      nas=2*nas+qa
      if(j/nb.eq.0)axx(k1+1,j+k1) = -bea(1,nas)
      if(j/nb.eq.1)axx(k1+1,j+k1) = +gaa(1,nas)
  
      !bn+1
      if (nas.le.nabp1) then
        if(j/nb.eq.0)axx(1+k1,j+k1+1)= c1(1)*(nas+2.)*(nas+abs(mmm)+1)/(2*nas+3.)/(nas+1.) &
          & -a2*c1(1)*3*(nas+2)*(nas**2+2*nas-2.) &
          & /(2*nas+5.)/(2*nas+3.)/(2*nas-1.) &
          & +6*c2(1)*(0,1)*mmm*(nas+abs(mmm)+1.0)/nas/(nas+1)/(2*nas+3)  
      endif 
      !bn-1
      if (nas.ge.nabm1) then
        if(j/nb.eq.0)axx(1+k1,j+k1-1)= c1(1)*(nas-1.)*(nas-abs(mmm))/(2*nas-1.)/nas &
          & -a2*c1(1)*3*(nas**2-3.)*(nas-1.) &
          & /(2*nas-3.)/(2*nas+3.)/(2*nas-1.) &
          &  +6*c2(1)*(0,1)*mmm*(nas-abs(mmm))/nas/(nas+1)/(2*nas-1.0)
      endif
      !bn+3
      if (nas.le.nabp3) then
        if(j/nb.eq.0)axx(1+k1,j+k1+3)= &
          & -a2*c1(1)*(nas+2.) &
          & *(nas+3.)*(nas+4.)/(2*nas+3.)/(2*nas+5.)/(2*nas+7.)
      endif
      !bn-3
      if (nas.ge.nabm3) then
       if(j/nb.eq.0)axx(1+k1,j+k1-3)= & 
         & -a2*c1(1)*(nas-1.)*(nas-2.)*(nas-3.) &
         & /(2*nas-1.)/(2*nas-5.)/(2*nas-3.)  
      endif
  
      !fanm2
      if (nas.ge.nafm2) then
        if(j/nb.eq.0)axx(1+k1,j+k1-2)=-fanm2(1,nas)
      endif
      !fanp2          
      if (nas.le.nafp2) then
        if(j/nb.eq.0)axx(1+k1,j+k1+2)=-fanp2(1,nas)
      endif
      !fanm2p 
      if (nas.ge.nafm2) then
        if(j/nb.eq.1)axx(1+k1,j+k1-2)= -fanm2p(1,nas)
      endif
      !fanp2p 
      if (nas.le.nafp2) then
        if(j/nb.eq.1)axx(1+k1,j+k1+2)= -fanp2p(1,nas)
      endif
      !fanm4
      if (nas.ge.nafm2+2) then
        if(j/nb.eq.0)axx(1+k1,j+k1-4)=-fanm4(1,nas)
      endif
      !fanp4          
      if (nas.le.nafp2-2) then
        if(j/nb.eq.0)axx(1+k1,j+k1+4)=-fanp4(1,nas)
      endif
      !fanm4p 
      if (nas.ge.nafm2+2) then
        if(j/nb.eq.1)axx(1+k1,j+k1-4)= -fanm4p(1,nas)
      endif
      !fanp4p 
      if (nas.le.nafp2-2) then
        if(j/nb.eq.1)axx(1+k1,j+k1+4)= -fanp4p(1,nas)
      endif
      !fanm6
      if (nas.ge.nafm2+4) then
        if(j/nb.eq.0)axx(1+k1,j+k1-6)=-fanm6(1,nas)
      endif
      !fanp6          
      if (nas.le.nafp2-4) then
        if(j/nb.eq.0)axx(1+k1,j+k1+6)=-fanp6(1,nas)
      endif
      !fanm6p 
      if (nas.ge.nafm2+4)  then
        if(j/nb.eq.1)axx(1+k1,j+k1-6)= -fanm6p(1,nas)
      endif
      !fanp6p 
      if (nas.le.nafp2-4) then
        if(j/nb.eq.1)axx(1+k1,j+k1+6)= -fanp6p(1,nas)
      endif
    enddo
  enddo
  
  !     runs col 2*nb  inner bn
  do j = 1, 3*nb, nb   
    nas = 0
    nas_count = 0
    do k1 = nbi, nbf, 2 
      nas_count = nas_count+1
      nas = nas_count-1
      nas=2*nas+qb
      if(j/nb.eq.0)axx(k1+1,j+k1) =-beb(1,nas) &
       & -bct*4*x1*alb(1,nas)/(2*hh-3*x1)           ! bct=1 perf con         
      if(j/nb.eq.1)axx(k1+1,j+k1) = gab(1,nas) &    ! bct =0 effect. perf.
       & +bct*x1*alb(1,nas)/(2*hh-3*x1)      
  
      !an-1     !alpha-omega + alpha^2
      if (nas.ge.nbam1) then
        if(j/nb.eq.0)axx(1+k1,j+k1-1) = beqm1(1,nas) + omem1(1, nas)
        if(j/nb.eq.1)axx(1+k1,j+k1-1) = gaqm1(1,nas) + omem1p(1,nas)
      endif
      !an+1     !alpha-omega + alpha^2
      if (nas.le.nbap1) then
        if(j/nb.eq.0)axx(1+k1,j+k1+1) = beqp1(1,nas) + omep1(1, nas)
        if(j/nb.eq.1)axx(1+k1,j+k1+1) = gaqp1(1,nas) + omep1p(1,nas)
      endif
      !an-3     !omega
      if (nas.ge.nbam3) then
        if(j/nb.eq.0)axx(1+k1,j+k1-3)= beqm3(1,nas) + omem3(1, nas)
        if(j/nb.eq.1)axx(1+k1,j+k1-3)= gaqm3(1,nas) + omem3p(1,nas)
      endif
      !an+3     !omega
      if (nas.le.nbap3) then
        if(j/nb.eq.0)axx(1+k1,j+k1+3)= beqp3(1,nas) + omep3(1, nas)
        if(j/nb.eq.1)axx(1+k1,j+k1+3)= gaqp3(1,nas) + omep3p(1,nas)
      endif
      !an-5
      if (nas.ge.nbam3+2) then
        if(j/nb.eq.0)axx(1+k1,j+k1-5)= omem5(1, nas)
        if(j/nb.eq.1)axx(1+k1,j+k1-5)= omem5p(1,nas)
      endif
      !an+3     !omega
      if (nas.le.nbap3-2) then
        if(j/nb.eq.0)axx(1+k1,j+k1+5)= omep5(1, nas)
        if(j/nb.eq.1)axx(1+k1,j+k1+5)= omep5p(1,nas)
      endif
  
      !fbnm4
      if (nas.ge.nbfm2+2) then
        if(j/nb.eq.0)axx(1+k1,j+k1-4)=-fbnm4(1,nas) &
          & -bct*4*x1*fbnm4p(1,nas)/(2*hh-3.*x1)
      endif  
      !fbnp4
      if (nas.le.nbfp2-2) then
        if(j/nb.eq.0)axx(1+k1,j+k1+4)=-fbnp4(1,nas) &
          & -bct*4*x1*fbnp4p(1,nas)/(2*hh-3.*x1)
      endif  
      !fbnm4p
      if (nas.ge.nbfm2+2) then
        if(j/nb.eq.1)axx(1+k1,j+k1-4)=-fbnm4p(1,nas) + &
          & bct*x1*fbnm4p(1,nas)/(2*hh-3.*x1)
      endif  
      !fbnp4p
      if (nas.le.nbfp2-2) then
        if(j/nb.eq.1)axx(1+k1,j+k1+4)=-fbnp4p(1,nas) + &
          & bct*x1*fbnp4p(1,nas)/(2*hh-3.*x1)
      endif  
  
      !fbnm6
      if (nas.ge.nbfm2+4) then
        if(j/nb.eq.0)axx(1+k1,j+k1-6)=-fbnm6(1,nas) &  
          & -bct*4*x1*fbnm6p(1,nas)/(2*hh-3.*x1)
      endif  
      !fbnp6
      if (nas.le.nbfp2-4) then
        if(j/nb.eq.0)axx(1+k1,j+k1+6)=-fbnp6(1,nas) &  
          & -bct*4*x1*fbnp6p(1,nas)/(2*hh-3.*x1)
      endif  
      !fbnm6p
      if (nas.ge.nbfm2+4) then
        if(j/nb.eq.1)axx(1+k1,j+k1-6)=-fbnm6p(1,nas) + &
          & bct*x1*fbnm6p(1,nas)/(2*hh-3.*x1)
      endif  
      !fbnp6p
      if (nas.le.nbfp2-4) then
        if(j/nb.eq.1)axx(1+k1,j+k1+6)=-fbnp6p(1,nas) + &
          & bct*x1*fbnp6p(1,nas)/(2*hh-3.*x1)
      endif  
      !fbnm2
      if (nas.ge.nbfm2) then
        if(j/nb.eq.0)axx(1+k1,j+k1-2)=-fbnm2(1,nas) & 
          & -bct*4*x1*fbnm2p(1,nas)/(2*hh-3.*x1)
      endif  
      !fbnp2
      if (nas.le.nbfp2) then
        if(j/nb.eq.0)axx(1+k1,j+k1+2)=-fbnp2(1,nas) &  
          & -bct*4*x1*fbnp2p(1,nas)/(2*hh-3.*x1)
      endif  
      !fbnm2p
      if (nas.ge.nbfm2) then
        if(j/nb.eq.1)axx(1+k1,j+k1-2)=-fbnm2p(1,nas) + &
          & bct*x1*fbnm2p(1,nas)/(2*hh-3.*x1)
      endif  
      !fbnp2p
      if (nas.le.nbfp2) then
        if(j/nb.eq.1)axx(1+k1,j+k1+2)=-fbnp2p(1,nas) + &
         & bct*x1*fbnp2p(1,nas)/(2*hh-3.*x1)
      endif  
    enddo
  enddo
  
  i = nt-nb+1
  k =i/nb      
  do j = 1, 3*nb, nb   ! runs col 2*nb  !outer an
    nas=0
    nas_count=0
    do k1 = nai, naf, 2   
      nas_count = nas_count+1
      nas = nas_count-1
      nas=2*nas+qa
      if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1)=  ala(k+1,nas)
      if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1)= -bea(k+1,nas)
  
      !bn+1
      if (nas.le.nabp1) then
        if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1+1)= c1(k+1)*(nas+2.)*(nas+abs(mmm)+1)/(2*nas+3.)/(nas+1.) &   
          & -a2*c1(k+1)*3*(nas+2)*(nas**2+2*nas-2.) &
          & /(2*nas+5.)/(2*nas+3.)/(2*nas-1.) &
          & +6*c2(k+1)*(0,1)*mmm*(nas+abs(mmm)+1.0)/nas/(nas+1)/(2*nas+3)
      endif
      !bn-1
      if (nas.ge.nabm1) then
        if(j/nb.eq.1)axx(i+k1,j+(k-1)*nb+k1-1)= c1(k+1)*(nas-1.)*(nas-abs(mmm))/(2*nas-1.)/nas  &
          & -a2*c1(k+1)*3*(nas**2-3.)*(nas-1.)/(2*nas-3.)/(2*nas+3.)/(2*nas-1.) &
          &  +6*c2(k+1)*(0,1)*mmm*(nas-abs(mmm))/nas/(nas+1)/(2*nas-1.0)
      endif
      !bn+3
      if (nas.le.nabp3) then
        if(j/nb.eq.1)axx(i+k1,j+(k-1)*nb+k1+3)= &
          & -a2*c1(k+1)*(nas+2)*(nas+3)*(nas+4)/(2*nas+3.)/(2*nas+5.)/(2*nas+7.)
      endif
      !bn-3
      if (nas.ge.nabm3) then
        if(j/nb.eq.1)axx(i+k1,j+(k-1)*nb+k1-3)= &
          & -a2*c1(k+1)*(nas-1.)*(nas-2.)*(nas-3.)/(2*nas-1.)/(2*nas-5.)/(2*nas-3.)  
      endif
  
      !fanm2p 
      if (nas.ge.nafm2) then
        call bess(beta, nas-2, ffc)
        if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1-2)= -(-fanm2p(k+1,nas)-fanm2p(k+1,nas)/(2*hh*(nas-1.+ffc)+3))
        if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1-2)= -( fanm2(k+1,nas) +fanm2p(k+1,nas)*4/(2*hh*(nas-1.+ffc)+3))
      endif
      !fanp2p 
      if (nas.le.nafp2) then
        call bess(beta, nas+2, ffc)
        if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1+2)= -(-fanp2p(k+1,nas)-fanp2p(k+1,nas)/(2*hh*(nas+3.+ffc)+3))
        if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1+2)= -( fanp2(k+1,nas) +fanp2p(k+1,nas)*4/(2*hh*(nas+3.+ffc)+3))
      endif
  
      !fanm4p 
      if (nas.ge.nafm2+2) then
        call bess(beta, nas-2, ffc)
        if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1-4)=-(-fanm4p(k+1,nas)-fanm4p(k+1,nas)/(2*hh*(nas-1.+ffc)+3))
        if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1-4)=-( fanm4(k+1,nas) +fanm4p(k+1,nas)*4/(2*hh*(nas-1.+ffc)+3))
      endif
      !fanp4p 
      if (nas.le.nafp2-2) then
        call bess(beta, nas+2, ffc)
        if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1+4)=-(-fanp4p(k+1,nas)-fanp4p(k+1,nas)/(2*hh*(nas+3.+ffc)+3))
        if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1+4)=-( fanp4(k+1,nas) +fanp4p(k+1,nas)*4/(2*hh*(nas+3.+ffc)+3))
      endif
  
      !fanm6p 
      if (nas.ge.nafm2+4) then
        call bess(beta, nas-2, ffc)
        if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1-6)= -(-fanm6p(k+1,nas)-fanm6p(k+1,nas)/(2*hh*(nas-1.+ffc)+3))
        if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1-6)= -( fanm6(k+1,nas) +fanm6p(k+1,nas)*4/(2*hh*(nas-1.+ffc)+3))
      endif
  
      !fanp6p 
      if (nas.le.nafp2-4) then
        call bess(beta, nas+2, ffc)
        if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1+6)= -(-fanp6p(k+1,nas)-fanp6p(k+1,nas)/(2*hh*(nas+3.+ffc)+3))
        if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1+6)= -( fanp6(k+1,nas) +fanp6p(k+1,nas)*4/(2*hh*(nas+3.+ffc)+3))
      endif
    enddo
  enddo
  
  !runs col 2*nb ! outer bn
  do j = 1, 3*nb, nb   
    nas=0
    nas_count=0
    do k1 = nbi, nbf, 2   
      nas_count = nas_count+1
      nas = nas_count-1
      nas = 2*nas+qb   
      if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1)=+alb(k+1,nas)
      if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1)=-beb(k+1,nas)
  
      !an-1     !alpha-omega + alpha^2
      if (nas.ge.nbam1) then
        call bess(beta,nas-1,ffc)
        if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1-1)= &
          & alqm1(k+1,nas)-gaqm1(k+1,nas)/(2*hh*(nas+ffc)+3)-omem1p(k+1,nas)*(1+1/(2*hh*(nas+ffc)+3))
        if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1-1)=beqm1(k+1,nas) &
          & +4*gaqm1(k+1,nas)/(2*hh*(nas+ffc)+3)+omem1(k+1,nas)+4*omem1p(k+1,nas)/(2*hh*(nas+ffc)+3)
      endif
  
      !an+1     !alpha-omega + alpha^2
      if (nas.le.nbap1) then
        call bess(beta,nas+1,ffc)
        if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1+1)= &
          & alqp1(k+1,nas) - gaqp1(k+1,nas)/(2*hh*(nas+2+ffc)+3)-omep1p(k+1,nas)*(1.+1./(2*hh*(nas+2+ffc)+3) )
        if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1+1)= &
          & beqp1(k+1,nas) + 4*gaqp1(k+1,nas)/(2*hh*(nas+2+ffc)+3)+omep1(k+1,nas)+4*omep1p(k+1,nas)/(2*hh*(nas+2+ffc)+3)
      endif
  
      !an-3     !alpha-omega + alpha^2
      if (nas.ge.nbam3) then
        call bess(beta,nas-3,ffc)
        if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1-3)= &
          & -omem3p(k+1,nas)*(1.+1./(2*hh*(nas-2+ffc)+3))+alqm3(k+1,nas)*(1.+1./(2*hh*(nas-2+ffc)+3))             
        if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1-3)= &
          & +omem3(k+1, nas)+4*omem3p(k+1,nas)/(2*hh*(nas-2+ffc)+3)+beqm3(k+1, nas)-4*alqm3(k+1,nas)/(2*hh*(nas-2+ffc)+3)    
      endif
  
      !an+3     !alpha-omega + alpha^2
      if (nas.le.nbap3) then
        call bess(beta,nas+3,ffc)
        if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1+3)= &
          & -omep3p(k+1,nas)*(1.+1./(2*hh*(nas+4+ffc)+3))+alqp3(k+1,nas)*(1.+1./(2*hh*(nas+4+ffc)+3))             
        if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1+3)= &
          & +omep3(k+1, nas)+4*omep3p(k+1,nas)/(2*hh*(nas+4+ffc)+3) &
          & +beqp3(k+1, nas)-4*alqp3(k+1,nas)/(2*hh*(nas+4+ffc)+3)     
      endif
  
      !an-5     !alpha-omega 
      if (nas.ge.nbam3+2) then
        call bess(beta,nas-3,ffc)
        if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1-5)=-omem5p(k+1,nas)*(1.+1./(2*hh*(nas-2+ffc)+3) )
        if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1-5)=+omem5(k+1,nas)+4*omem5p(k+1,nas)/(2*hh*(nas-2+ffc)+3)
      endif
  
      !an+5     !alpha-omega 
      if (nas.le.nbap3-2) then
        call bess(beta,nas+3,ffc)
        if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1+5)=-omep5p(k+1,nas)*(1.+1./(2*hh*(nas+4+ffc)+3) )
        if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1+5)=+omep5(k+1,nas)+4*omep5p(k+1,nas)/(2*hh*(nas+4+ffc)+3)
      endif
  
      !fbnm2
      if (nas.ge.nbfm2) then
        if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1-2)=-fbnm2(k+1,nas) 
      endif  
      !fbnp2
      if (nas.le.nbfp2) then
        if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1+2)=-fbnp2(k+1,nas) !here
      endif  
      !fbnm2p
      if (nas.ge.nbfm2) then
        if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1-2)=fbnm2p(k+1,nas) 
      endif  
      !fbnp2p
      if (nas.le.nbfp2) then
        if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1+2)=fbnp2p(k+1,nas) 
      endif  
  
      !fbnm4
      if (nas.ge.nbfm2+2) then
        if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1-4)=-fbnm4(k+1,nas) 
      endif  
      !fbnp4
      if (nas.le.nbfp2-2) then
        if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1+4)=-fbnp4(k+1,nas)
      endif  
      !fbnm4p
      if (nas.ge.nbfm2+2) then
        if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1-4)=fbnm4p(k+1,nas) 
      endif  
      !fbnp4p
      if (nas.le.nbfp2-2) then
        if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1+4)=fbnp4p(k+1,nas) 
      endif  
  
      !fbnm6
      if (nas.ge.nbfm2+4) then
        if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1-6)=-fbnm6(k+1,nas) 
      endif  
      !fbnp6
      if (nas.le.nbfp2-4) then
        if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1+6)=-fbnp6(k+1,nas)
      endif  
      !fbnm6p
      if (nas.ge.nbfm2+4) then
        if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1-6)= fbnm6p(k+1,nas) 
      endif  
      !fbnp6p
      if (nas.le.nbfp2-4) then
        if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1+6)= fbnp6p(k+1,nas) 
      endif  
    enddo
  enddo
  
  !-------     write the inner blocks of the matrix------- an
      
  do i = nb+1, nt-nb, nb   ! start at each nb row 
    k = i/nb                    ! count the mesh-point
    do j = 1, 3*nb, nb       ! runs col 
      nas_count = 0
      nas = 0
      do k1 = nai, naf, 2  
        nas_count = nas_count+1
        nas = nas_count-1
        nas=2*nas+qa
  
        if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1)=  ala(k+1,nas)
        if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1)= -bea(k+1,nas)
        if (j/nb.eq.2) axx(i+k1,j+(k-1)*nb+k1)=  gaa(k+1,nas)
  
        !bn+1
        if (nas.le.nabp1) then
          if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1+1)= &
            & c1(k+1)*(nas+2.)*(nas+abs(mmm)+1)/(2*nas+3.)/(nas+1.) &   
            &-a2*c1(k+1)*3*(nas+2)*(nas**2+2*nas-2.) &
            & /(2*nas+5.)/(2*nas+3.)/(2*nas-1.) &
            & +6*c2(k+1)*(0,1)*mmm*(nas+abs(mmm)+1.0)/nas/(nas+1)/(2*nas+3) ! 4.11b ! check sign
        endif
  
        !bn-1
        if (nas.ge.nabm1) then
          if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1-1)= c1(k+1)*(nas-1.)*(nas-abs(mmm))/(2*nas-1.)/nas &    
            &-a2*c1(k+1)*3*(nas**2-3.)*(nas-1.)/(2*nas-3.)/(2*nas+3)/(2*nas-1.) &  !check
            &  +6*c2(k+1)*(0,1)*mmm*(nas-abs(mmm))/nas/(nas+1)/(2*nas-1.0)
        endif
  
        !bn+3
        if (nas.le.nabp3) then
          if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1+3)= &
            & -a2*c1(k+1)*(nas+2.)*(nas+3.)*(nas+4.) &
            & /(2*nas+3.)/(2*nas+5.)/(2*nas+7.)
        endif
  
        !bn-3
        if (nas.ge.nabm3) then
          if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1-3)= &
            & -a2*c1(k+1)*(nas-1.)*(nas-2.)*(nas-3.)/(2*nas-1.) &
            & /(2*nas-5.)/(2*nas-3.)  
        endif
  
        !fanm2
        if (nas.ge.nafm2) then
          if(j/nb.eq.1)axx(i+k1,j+(k-1)*nb+k1-2)=-fanm2(k+1,nas)
        endif
        !fanp2          
        if (nas.le.nafp2) then
          if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1+2)=-fanp2(k+1,nas)
        endif
        !fanm2p 
        if (nas.ge.nafm2) then
          if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1-2)= fanm2p(k+1,nas)
          if (j/nb.eq.2) axx(i+k1,j+(k-1)*nb+k1-2)=-fanm2p(k+1,nas)
        endif
        !fanp2p 
        if (nas.le.nafp2) then
          if (j/(nb).eq.0) axx(i+k1,j+(k-1)*nb+k1+2)=  fanp2p(k+1,nas)
          if (j/(nb).eq.2) axx(i+k1,j+(k-1)*nb+k1+2)= -fanp2p(k+1,nas)
        endif
  
        !fanm4
        if (nas.ge.nafm2+2) then
          if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1-4)=-fanm4(k+1,nas)
        endif
        !fanp4          
        if (nas.le.nafp2-2) then
          if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1+4)=-fanp4(k+1,nas)
        endif
        !fanm4p 
        if (nas.ge.nafm2+2) then
          if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1-4)=+fanm4p(k+1,nas)
          if (j/nb.eq.2) axx(i+k1,j+(k-1)*nb+k1-4)=-fanm4p(k+1,nas)
        endif
        !fanp4p 
        if (nas.le.nafp2-2) then
          if (j/(nb).eq.0) axx(i+k1,j+(k-1)*nb+k1+4)=+fanp4p(k+1,nas)
          if (j/(nb).eq.2) axx(i+k1,j+(k-1)*nb+k1+4)=-fanp4p(k+1,nas)
        endif
  
        !fanm6
        if (nas.ge.nafm2+4) then
          if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1-6)=-fanm6(k+1,nas)
        endif
        !fanp6          
        if (nas.le.nafp2-4) then
          if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1+6)=-fanp6(k+1,nas)
        endif
        !fanm6 
        if (nas.ge.nafm2+4) then
          if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1-6)=+fanm6p(k+1,nas)
          if (j/nb.eq.2) axx(i+k1,j+(k-1)*nb+k1-6)=-fanm6p(k+1,nas)
        endif
        !fanp6p 
        if (nas.le.nafp2-4) then
          if (j/(nb).eq.0) axx(i+k1,j+(k-1)*nb+k1+6)=+fanp6p(k+1,nas)
          if (j/(nb).eq.2) axx(i+k1,j+(k-1)*nb+k1+6)=-fanp6p(k+1,nas)
        endif
  
      enddo
    enddo
  enddo
  
  !-------     write the inner blocks of the matrix------- bn
  
  do i = nb+1, nt-nb, nb   ! start at each 2*nb row 
    k = i/(nb)                   ! count the mesh-point
    do j = 1, 3*nb, nb       ! runs col 2*nb*i
    nas_count = 0
    nas =0
      do k1 = nbi, nbf, 2   
        nas_count = nas_count+1
        nas = nas_count-1
        nas=2*nas+qb
  
        if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1) =  alb(k+1,nas)
        if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1) = -beb(k+1,nas)
        if (j/nb.eq.2) axx(i+k1,j+(k-1)*nb+k1) =  gab(k+1,nas)
        !an-1     !alpha-omega + alpha^2
        if (nas.ge.nbam1) then
          if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1-1)= alqm1(k+1,nas) - omem1p(k+1,nas)
          if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1-1)= beqm1(k+1,nas) + omem1(k+1, nas)
          if (j/nb.eq.2) axx(i+k1,j+(k-1)*nb+k1-1)= gaqm1(k+1,nas) + omem1p(k+1,nas)
        endif
        !an+1     !alpha-omega + alpha^2
        if( nas.le.nbap1) then
          if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1+1)= alqp1(k+1,nas) - omep1p(k+1,nas)     
          if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1+1)= beqp1(k+1,nas) + omep1(k+1, nas)
          if (j/nb.eq.2) axx(i+k1,j+(k-1)*nb+k1+1)= gaqp1(k+1,nas) + omep1p(k+1,nas)
        endif
        !an-3     !alpha-omega + alpha^2
        if (nas.ge.nbam3) then
          if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1-3)= alqm3(k+1,nas) - omem3p(k+1,nas)
          if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1-3)= beqm3(k+1,nas) + omem3(k+1, nas)
          if (j/nb.eq.2) axx(i+k1,j+(k-1)*nb+k1-3)= gaqm3(k+1,nas) + omem3p(k+1,nas)
        endif
        !an+3     !alpha-omega + alpha^2
        if (nas.le.nbap3) then
          if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1+3)= alqp3(k+1,nas) - omep3p(k+1,nas)     
          if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1+3)= beqp3(k+1,nas) + omep3(k+1, nas)
          if (j/nb.eq.2) axx(i+k1,j+(k-1)*nb+k1+3)= gaqp3(k+1,nas) + omep3p(k+1,nas)
        endif
        !an-5     !alpha-omega + alpha^2
        if (nas.ge.nbam3+2) then
          if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1-5)= - omem5p(k+1,nas)
          if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1-5)= + omem5(k+1, nas)
          if (j/nb.eq.2) axx(i+k1,j+(k-1)*nb+k1-5)= + omem5p(k+1,nas)
        endif
        !an+5     !alpha-omega + alpha^2
        if (nas.le.nbap3-2) then
          if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1+5)= - omep5p(k+1,nas)     
          if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1+5)= + omep5(k+1, nas)
          if (j/nb.eq.2) axx(i+k1,j+(k-1)*nb+k1+5)= + omep5p(k+1,nas)
        endif
        !fbnm2
        if (nas.ge.nbfm2) then
          if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1-2)=-fbnm2(k+1,nas)
        endif  
        !fbnp2
        if (nas.le.nbfp2) then
          if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1+2)=-fbnp2(k+1,nas) !here
        endif  
        !fbnm2p
        if (nas.ge.nbfm2) then
          if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1-2)= fbnm2p(k+1,nas)
          if (j/nb.eq.2) axx(i+k1,j+(k-1)*nb+k1-2)=-fbnm2p(k+1,nas)
        endif  
        !fbnp2p
        if (nas.le.nbfp2) then
          if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1+2)= fbnp2p(k+1,nas)  
          if (j/nb.eq.2) axx(i+k1,j+(k-1)*nb+k1+2)=-fbnp2p(k+1,nas) 
        endif  
        !fbnm4
        if (nas.ge.nbfm2+2) then
          if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1-4)=-fbnm4(k+1,nas)
        endif  
        !fbnp4
        if (nas.le.nbfp2-2) then
          if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1+4)=-fbnp4(k+1,nas)
        endif  
        !fbnm4p
        if (nas.ge.nbfm2+2) then
          if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1-4)= fbnm4p(k+1,nas)
          if (j/nb.eq.2) axx(i+k1,j+(k-1)*nb+k1-4)=-fbnm4p(k+1,nas)
        endif  
        !fbnp4p
        if (nas.le.nbfp2-2) then
          if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1+4)= fbnp4p(k+1,nas) 
          if (j/nb.eq.2) axx(i+k1,j+(k-1)*nb+k1+4)=-fbnp4p(k+1,nas) 
        endif  
        !fbnm6
        if (nas.ge.nbfm2+4) then
          if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1-6)=-fbnm6(k+1,nas)
        endif  
        !fbnp6
        if (nas.le.nbfp2-4) then
          if (j/nb.eq.1) axx(i+k1,j+(k-1)*nb+k1+6)=-fbnp6(k+1,nas)
        endif  
        !fbnm6p
        if (nas.ge.nbfm2+4) then
          if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1-6)= fbnm6p(k+1,nas)
          if (j/nb.eq.2) axx(i+k1,j+(k-1)*nb+k1-6)=-fbnm6p(k+1,nas)
        endif  
        !fbnp6p
        if (nas.le.nbfp2-4) then
          if (j/nb.eq.0) axx(i+k1,j+(k-1)*nb+k1+6)= fbnp6p(k+1,nas) 
          if (j/nb.eq.2) axx(i+k1,j+(k-1)*nb+k1+6)=-fbnp6p(k+1,nas) 
        endif  
      enddo
    enddo
  enddo
  
  !----- lapack solvers
  
  if(mmm.ne.0) flg=1
  if(flg.eq.1) then
     call zgeev(jobvl, jobvr, nt, axx, nt, ww, cvr, nt, cvr, nt, &
                &  work, lwork, rwork, info )
     ! rescale the eigenvalues 
     do i=1, nt
       wr(i) = real(ww(i))/hh**2    
       wi(i) = aimag(ww(i))/hh**2
       inde(i) = i
     enddo
  else
    forall (i=1:nt,j=1:nt) axr(i,j) = real(axx(i,j))
    axr = real(axx)
  
    call dgeev(jobvl, jobvr, nt, axr, nt, wr, wi, vr, nt, & 
               &  vr, nt, work, lwork, info )
  
    do i=1, nt
      wr(i) = wr(i)/hh**2    
      wi(i) = wi(i)/hh**2
      inde(i) = i
    enddo
  endif
  
  call sort2(nt,wr,inde)     
  
  indeg=int(inde)
  
  rate = wr(nt) 
  imag = wi(int(indeg(nt)))
  
  do i2 =1,nb
    k = 1          
    k = 0          
    k = -1         
    do i3 = i2,nt,nb
      k = k + 1   
      egr(k+1,i2) = vr(i3,indeg(nt))
      egi(k+1,i2) = vr(i3,indeg(nt-1)) 
    enddo
  enddo
  
  ki=0            
  do i=1,na,2 
    ki=ki+1 
    do j = 1,np
      agr(j,ki) = egr(j,i)
      agi(j,ki) = egi(j,i)
    enddo
  enddo
  
  ki=0
  do i=2,nb,2 
   ki=ki+1
   do j = 1,np
     bgr(j,ki) = egr(j,i)
     bgi(j,ki) = egi(j,i)
    enddo
  enddo
  
  eep=(maxval(abs(bgr)) + maxval(abs(bgi)))/(maxval(abs(agr)) + maxval(abs(agi)))
  
  !------------------------------------------------------------------
  rt = rate  
  vtu = turb 
  !------------------------------------------------------------------
  
  !------------ write the period in years ---------------------------
  
  facp= 2.e0*pi*(sr*sr0)**2.e0/eta/3.1536e7
  do kk = 0,5
    write(*,'(i4, 1p,3e13.5)') kk, wr(nt-kk), wi(int(indeg(nt-kk))), facp/wi(int(indeg(nt-kk)))
  enddo
  
  write(*,*)
  
  do kk = 0,9
    reg(kk+1) =  wr(nt-kk)
    ieg(kk+1) =  wi(int(indeg(nt-kk)))
  enddo
  
  if (jobvr.eq.'v') then
    call writefield
  endif
  
  return
end

