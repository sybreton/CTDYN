c
c   Radler73 Zur DynamoTheorie / Notes 
c
c   the general structure of the equation is:
c
c   dPhi/dt = Gn + Dn + Nabla Phi  3.17a
c   dPsi/dt = Fn + Cn + Nabla Psi  3.17b
c   
c   Phi - a-coefficients 
c   Psi - b-coefficients
c
c   Gn : alpha-omega  m =1   
c   Dn : alpha B term: written as: alpha * cos * ( 1 - c3 * cos^2)  if c3=0,  standard , if c3=1 we have cos * sin^2
c   c3=1 only correct for m=0 !!!!!!!! 
c   note that alpha = -alpha1 of Rad73
c
c   General structure of the type:
c
c  a Phi'' + b Phi' - c Phi ....
c
c  a ( Phi_i+1 -2 Phi_i +Phi_i-1)/hh**2 + b (Phi_i+1 -Phi_i-1)/(2 hh) - c Phi ....
c
c  (a+h2 b)  Phi_i+1  - (2 a + hh**2 c) Phi_i +(a-h2 b) Phi_i-1 
c
c  alpha_i Phi_i-1  - Beta_i Phi_i +gamma_i Phi_i+1
c
c  where h2 = hh/2
c
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DYNAMO(turb,rate)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT REAL(A-H,O-Z)
 
      REAL imag
      COMMON/part/vtu,rt,imag,co,c_u,beta,ffree,betb,etep,etet,xbt,xbo
      real eep
      REAL REG(10),IEG(10)
      COMMON/eira/REG,IEG,eep

C      DEFINITIONS:
C      turb= C_alpha*C_omega
C      C_alpha = alpha_zero*SR/eta
C      C_omega = omega_zero*SR**2/eta
C      COMMON param

      INTEGER II,IT                 ! LOOP COUNTER
      COMMON/ipar/II,IT

      REAL xa1,xa2,xa3,xb,xda1,xda2
      COMMON/apar/xa1,xa2,xa3,xb,xda1,xda2

      REAL edr,xe1,xde1,hd
      COMMON/epar/edr,xe1,xde1,hd

      REAL x_in, bct, c3,mmm,SR,rotp, gd, aqu, flg
      COMMON/ppar/x_in,bct,c3,mmm,SR,rotp,gd, aqu, flg
      
      REAL s0,s2,s4,s6, A2P, A4P,xm
      COMMON/PSI/s0,s2,s4,s6,A2P,A4P,xm 

      character(len=128)    ::  dir
      CHARACTER*8 ANS1,ANS2,ANS3,ANS4
      COMMON/var3/ANS1,ANS2,ANS3,ANS4, dir

      CHARACTER*2 JOBVR,JOBVL
      COMMON/lap/JOBVR,JOBVL

      EXTERNAL version

C------------- temporary variables

      REAL OMEGAEQ
      REAL  NUEQ, ETA, SR0
      PARAMETER(SR0=6.955e10,NUEQ=460.7e-9)
c      PARAMETER(NUEQ=460.7e-9)
C
C     1./freq / seconds in a day = eq. rotation period for the Sun
C     1./(1.d-9*460) / 86400. =  25.161031       
C --- array section
C     C_\Omega = R^2*omg/eta = 2247 for eta = 1.d12
C
      include 'cio'
      include 'cvect'
      include 'ccov'

      REAL EGR(NP,NB), EGI(NP,NB)
      REAL AGR(NP,NB/2), AGI(NP,NB/2)
      REAL BGR(NP,NB/2), BGI(NP,NB/2)

      CHARACTER*1 JOBPR    
      CHARACTER*1 BALANC    
      CHARACTER*1 SENSE    
      REAL RWORK(2*NT)
      COMPLEX*16 WORK(LWORK)  
ccc      COMPLEX WORK(LWORK)  
      
      CHARACTER*43 ver
      INTEGER qa, qb

      COMPLEX*16 AXX(NT,NT)
      REAL   AXR(NT,NT)

      COMPLEX*16 drb, dra, sigmeno, sigpiu

      COMMON/parker/gam,zeta_r,ratio

      COMMON/vec/CVR

ccc      COMPLEX AXX(NT,NT)

      REAL c1(NP)                  ! coefficient dyn.  
      REAL c2(NP)

!------------------------------------------------------
      JOBPR = 'N'    ! Y = write all pij matrix       !!!!
      JOBVL = 'N'    ! V = calculate R and L eig.
!----------------  boundaries ------------------------------

c      t1 = secnds(0.0)
c      write(*,*) sr
c      stop
 

      x1 = x_in           !inner boundary
      x2 = 1.0            !outer boundary
      hh =(x2-x1)/(NP+1)  !stepsize: radial accuracy parm.
      h2 = hh/2.e0
      a2 = c3       ! cos^3 term in alpha  

C
C
C Check!!!! verifica la definizione di c3
C           verifica se i coeefficienti per c3 non zero sono validi anche per m=1, sembra di no
C
C rotp = periodo di rotazione in unita' di quello solare
C

      OMEGAEQ = NUEQ*2.e0*PI/rotp

      if(co.ne.0)then            
      ETA = (SR*SR0)**2*OMEGAEQ/co
      else
      ETA=1
      endif

c      write(*,*) ETA
c      write(*,*) SR, SR0, OMEGAEQ, co
c      stop


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      um = C_U*eta/(SR*SR0)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      ca = turb             

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C              WRITE ITERATIONS on screeen
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c in this case the coefficients are not computed....
      if(c3.ne.0.and.mmm.ne.0)then
        print*, 'NOT AVAILABLE'
        stop
      endif 

      if(ans2.eq.'D'.and.mmm.eq.0) write(*,*) '********** A0-MODE *********'
      if(ans2.eq.'D'.and.mmm.ne.0) write(*,*) '********** S1-MODE *********'
      if(ans2.eq.'Q'.and.mmm.eq.0) write(*,*) '********** S0-MODE *********'
      if(ans2.eq.'Q'.and.mmm.ne.0) write(*,*) '********** A1-MODE *********'

      write(*,'(1x,a,I4,1x,a,I4,1x,a,I4)') ' it=', it,' NA=', NB, ' NP=',NP

c      write(*,'(a,1p,e12.5,a,1p,e12.5,a,1p,e12.5,a,1p,
c     &   e12.5,a,1p,e12.5,a,1p,e12.5)' )
c     & ' c_alpha=',turb, '  c_omega=', co, '  reynolds=',c_u, ' f-f=', beta
c     & ,'  eta=',eta,'  um=',um 

      write(*,'(a,e12.5,a,e12.5)' ) ' c_alpha=',turb, ' c_omega=', co 
      write(*,'(a,e12.5,a,e12.5)' ) ' R_flow =',c_u,  ' f-f=    ', beta 

c      write(*,'(a,1p,e12.5,a,1p,e12.5,a,1p,e12.5,a,1p,
c     &   e12.5,a,1p,e12.5,a,1p,e12.5)' )
c     & ' c_alpha=',turb, '  c_omega=', co, '  reynolds=',c_u, ' f-f=', beta
c     & ,'  eta=',eta,'  um=',um 


C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     WRITE ALPHA AND ETA the first iteration/test run
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C


      IF(II.eq.1)THEN     
      open(22,status='unknown',file=trim(trim(dir)//'/alpha.dat'))
      DO 199 I=1, NP
          x = x1+I*HH
          call rot(x, om0, om0p, om2, om2p, om4, om4p) 
          call alt(x, alp, alp_p)     
          call eta_turb(x, e1, e2, e3)
          call stream(x,ar,arp,bt,btp, psi)

          write(22,'(15f14.6)') x,alp,alp_p,e1,e2,e3,psi,ar,arp,bt,btp,om0,om0p,om2,om2p

 199      continue
       close(22)
          ENDIF

C----------------- LOAD ARRAYS IN M and N

       nc              =  2 ! number of latitudinal cells /2 set to 2 by def.
       if(A4P.ne.0)nc  =  4 ! if A4P not zero there are two cells-> nc=4
       if(s6.ne.0) nc  =  6       
c
C THESE ARE THE CONVERSION FACTORS for P1N for Utheta
c
c      s4 = (35.e0/8.e0) * A4P
c      s2 = (3.e0/2.e0) * A2P -(15.e0/4.e0)*A4P
c
c      ss2 = (9.e0/(nc*(nc+1.e0)))*A2P-(75.e0/(nc*(nc+1.e0)))*A4P
c      s0= 15.e0/(2*(nc*(nc+1.e0)))*A4P-(3.e0/(nc*(nc+1.e0)))*A2P
c
 
      am = abs(mmm)

      do 211, J=1, NB
      x=0.0
      do 210, I=1, NP
      x=x1+i*hh

      call rot(x, om0, om0p, om2, om2p, om4, om4p)
      call alt(x, alp, alp_p)
      call eta_turb(x, e1, e2, e3)
      call stream(x,ar,arp, bt,btp, psi)


c bench test newtest

      rho=x**gd 
      rho=0
      psi=0
      psip=0
      psipp=0

      c1(i) = ca*alp*hh**2

c      c2(i) = -c_u*(3*psi/rho/x/x)*hh**2
c check sign 

      c2(i) =  c_u*ar*hh**2

c      betb=10*beta

       betb=beta

cc      betb=2000

c  alpha beta gamma of the notes for a-vectors
c  in the matrix we write  
c  
c  alpha_i B_i-1  - beta_i B_i  + gamma_i B_i+1    
c
c  -> alpha_i , beta_i, gamma_i  .... = AXX
c 
c  beta_i = 2+h^2 n(n+1)/x^2  + h^2 c_u f ...
c
c
       dra = h2*((-2*e1/x - hd * e2/2)
     & +ca*alp*mmm*(0,1)/J/(J+1.)   ! check sign
     & -c_u*((J*(J+1)-3*am**2)/J/(J+1.)/(2*J-1.)/(2*J+3.))
     & *(2*(J*(J+1)-3)*ar)) 

         bea(i,j) =2*e1+hh**2*((e1*J*(J+1)/x**2 -hd*e2/x  ) 
     &  +ca*alp*mmm*(0,1)/J/(J+1.)/x   ! check sign
     & +(co)*(+om0*mmm*(0,1) + (+1)*om2*mmm*(0,1)*(J*(J+1)-3*mmm**2)/(2*J-1)/(2*J-3)) 
     & -c_u*((J*(J+1)-3*am**2)/J/(J+1.)/(2*J-1.)/(2*J+3.))
     & *(2*(J*(J+1)-3)*ar/x-J*(J+1)*bt/x)) 

        gaa(i,j) =  e1-dra
        ala(i,j) =  e1+dra

c set outer boundary conditions for a-vectors

c        bea(NP,J) = bea(NP,J) - 4*gaa(NP,J)/(2*hh*(J+1)+3)
c        ala(NP,J) = ala(NP,J) -  gaa(NP,J)/(2*hh*(J+1)+3)     

c  alpha beta gamma of the notes for b-vectors

      beb(i,j) = e1*2.e+00 + hh**2*((e1*J*(J+1.e+00)/x**2-(e2/x)-hd*e2/x/2 -hd*e3/2  )  
     & +ca*(alp/x+alp_p)*mmm*(0,1)/J/(J+1.e+00) !check sign
     & + (co)*(+om0*mmm*(0,1) + (+1)*om2*(0,1)*mmm*(J*(J+1.e+00)-3.e+00*mmm**2)/(2*J-1.e+00)/(2*J+3.e+00)) 
     & -c_u*((J*(J+1)-3.*am**2)/J/(J+1)/(2*J-1.)/(2*J+3.))
     & *(2*(J*(J+1)-3)*(arp+ar/x)+bt*J*(J+1)/x)) 

      drb =  h2*((-2.*e1/x-e2 -hd*e2/2.)
     &  +ca*(alp)*mmm*(0,1)/J/(J+1.)   !check sign
c check sign c_u
     & -c_u*((J*(J+1)-3.*am**2)/J/(J+1)/(2*J-1.)/(2*J+3.))*(2*(J*(J+1)-3)*ar)) 

      gab(i,j) = e1 -  drb
      alb(i,j) = e1 +  drb

C--------------------  a-2 a+2 ------------------------------------------

      fanm2(i,j) = -C_U*hh**2*((J-2.)*(J-am-1.)*(J-am)/J/(2*J-3.)/(2*J-1))
     & * (3*ar/x-(J-1.)*bt/x)
  
      fanm2p(i,j) = -c_u*h2*((J-2.)*(J-am-1.)*(J-am)/J/(2*J-3.)/(2*J-1))*3*ar  

      fanp2(i,j) = -C_U*hh**2*((J+3.)*(J+am+1)*(J+am+2)/(J+1)/(2*J+3.)/(2*J+5))  
     & * (3*ar/x+(J+2)*bt/x)

      fanp2p(i,j) = -C_U*h2*((J+3.)*(J+am+1)*(J+am+2)/(J+1)/(2*J+3.)/(2*J+5.))*3*ar

c---------------------- a-4 a+4

      fanm4(i,j)  = 0
      fanp4(i,j)  = 0
      fanm4p(i,j) = 0
      fanp4p(i,j) = 0
      fanm6(i,j)  = 0
      fanp6(i,j)  = 0
      fanm6p(i,j) = 0
      fanp6p(i,j) = 0

c--------------------------------- b-2 b+2---------------------------------

      fbnm2(i,j) = -c_u*hh**2*
     & ((J-2.)*(J-am-1.)*(J-am)/J/(2.*J-3.)/(2.*J-1))
     & *(3.*(arp+ar/x)-J*bt/x)

      fbnm2p(i,j) = -c_u*h2*
     & ((J-2.)*(J-am-1.)*(J-am)/J/(2.*J-3.)/(2.*J-1))*3.*ar

      fbnp2(i,j) = -c_u*hh**2*
     & ((J+3.)*(J+am+1.)*(J+am+2)/(J+1.)/(2.*J+3.)/(2.*J+5.))
     & *(3.*(arp+ar/x)+(J+1.)*bt/x)

      fbnp2p(i,j) = -c_u*h2*
     & ((J+3.)*(J+am+1.)*(J+am+2)/(J+1.)/(2*J+3.)/(2*J+5.))*3.*ar

      fbnm4(i,j)  = 0
      fbnp4(i,j)  = 0
      fbnm4p(i,j) = 0
      fbnp4p(i,j) = 0
      fbnm6(i,j)  = 0
      fbnp6(i,j)  = 0
      fbnm6p(i,j) = 0
      fbnp6p(i,j) = 0

C --------------- omega arrays ---------------------------------------

      if(ans1.eq.'H4')then

      OMEP1(i,j) = co*hh**2* (B0M1(j+1)*om0p+B2M1(J+1)*om2p 
     &+B4M1(J+1)*om4p+2*C1S2M1(J+1)*om2/x+4*C3S2M1(J+1)*om4/x)

      OMEM1(i,j) = co*hh**2*(B0P1(j-1)*om0p+B2P1(J-1)*om2p 
     &+B4P1(J-1)*om4p+2*C1S2P1(J-1)*om2/x+4*C3S2P1(J-1)*om4/x)

      OMEP3(i,j)= co*hh**2*
     &       (B2M3(j+3)*om2p+2*om2*C1S2M3(j+3)/x+4*om4*C3S2M3(j+3)/x)

      OMEM3(i,j)= co*hh**2*
     &       (B2P3(j-3)*om2p+2*om2*C1S2P3(j-3)/x+4*om4*C3S2P3(j-3)/x)

      OMEP1P(i,j)= co*h2*(2*om2*C1S2M1(J+1)+4*om4*C3S2M1(J+1))
      OMEM1P(i,j)= co*h2*(2*om2*C1S2P1(J-1)+4*om4*C3S2P1(J-1))
      OMEP3P(i,j)= co*h2*(2*om2*C1S2M3(J+3)+4*om4*C3S2M3(J+3))
      OMEM3P(i,j)= co*h2*(2*om2*C1S2P3(J-3)+4*om4*C3S2P3(J-3))

      OMEP5(i,j)= co*hh**2*(4*om4*C3S2M5(j+5)/x)
      OMEM5(i,j)= co*hh**2*(4*om4*C3S2P5(j-5)/x)

      OMEP5P(i,j)= co*h2*(4*om4*C3S2M5(J+5))
      OMEM5P(i,j)= co*h2*(4*om4*C3S2P5(J-5))

      else

      OMEM1(i,j) = co*hh**2*(
     &  om0p*(J-1)*(J-abs(mmm))/(2*J-1)
     &  +(J-abs(mmm))*(6*J*((J-1)*(j+1)-5*mmm**2)*om2/x - (j-1)*((J-3)*J*(J+1)+3*(J+6)*mmm**2)*om2p)
     & /2/(J+1)/(2*J-3)/(2*J-1)/(2*J+3) )

      OMEP1(i,j) = co*hh**2*(
     &  -om0p*(J+2)*(J+abs(mmm)+1)/(2*J+3)
     &  +(abs(J+mmm)+1)*(6*(J+1)*(J*(J+2)-5*mmm**2)*om2/x+(J+2)*(J*(J+1)*(J+4)+3*(J-5)*mmm**2)*om2p) 
     &  /2/J/(2*J-1)/(2*J+3)/(2*J+5))

      OMEM3(I,J) = co*hh**2*(
     & -3*(J-3)*(J-abs(mmm)-2)*(J-abs(mmm)-1)*(J-abs(mmm))*(om2/x-(J-2)*om2p/2)
     &  /J/(2*J-5)/(2*J-3)/(2*J-1))      
 
  
      OMEP3(I,J) = co*hh**2*(
     & -3*(J+4)*(J+abs(mmm)+1)*(J+abs(mmm)+2)*(J+abs(mmm)+3)*(om2/x+(J+3)*om2p/2)
     & /(J+1)/(2*J+3)/(2*J+5)/(2*J+7))
     
      OMEM1P(I,J) = co*h2*((J-abs(mmm))*(6*J*((J-1)*(j+1)-5*mmm**2)*om2)
     & /2/(J+1)/(2*J-3)/(2*J-1)/(2*J+3) )
   

      OMEP1P(I,J) = co*h2*((abs(J+mmm)+1)*(6*(J+1)*(J*(J+2)-5*mmm**2)*om2)
     &  /2/J/(2*J-1)/(2*J+3)/(2*J+5))

      OMEM3P(I,J) = co*h2*(
     & -3*(J-3)*(J-abs(mmm)-2)*(J-abs(mmm)-1)*(J-abs(mmm))*(om2)
     &  /J/(2*J-5)/(2*J-3)/(2*J-1))

       OMEP3P(I,J) = co*h2*(
     & -3*(J+4)*(J+abs(mmm)+1)*(J+abs(mmm)+2)*(J+abs(mmm)+3)*(om2)
     &  /(J+1)/(2*J+3)/(2*J+5)/(2*J+7))

      endif

!     CHECK THE RELATIVE SIGN: in the OLD version "+" because they are imput with "-" in the matrix!!!!
C     CHECK OMEGA hh**2    
C
C     4.13a pnm2   
      fbnm2(I,J) = fbnm2(I,J) + 3*(0,1)*mmm*(J-am-1)*(J-am)*hh**2*om2/2./(2*J-3.)/(2*J-1.)
C     4.13a pnp2
      fbnp2(I,J) = fbnp2(I,J) + 3*(0,1)*mmm*(J+am+1)*(J+am+2)*hh**2*om2/2./(2*J+3.)/(2*J+5.)

      fanm2(I,J) = fanm2(I,J) + 3*(0,1)*mmm*(J-2)*(J-1)*(J-am-1)*(J-am)*hh**2*om2/2./J/(J+1)/(2.*J-3.)/(2*J-1.)

      fanp2(I,J) = fanp2(I,J) + 3*(0,1)*mmm*(J+2)*(J+3)*(J+am+1)*(J+am+2)*hh**2*om2/2./J/(J+1)/(2.*J+3.)/(2*J+5.)

C ------------------- alpha^2 -----------------------------------


C include the full alpha^2 term for aqu=1

      ct = aqu*ca   ! ORG
c     

      alqm1(i,j)=-ct*((alp-h2*(2*alp/x+alp_p))*((j-1.)*(J-abs(mmm))/(2*J-1.)/J  
     & -c3*(j**2-3.)*(j-1.)*3/(2*j-3.)/(2*j+3.)/(2*j-1.)))

c    
      alqp1(i,j)=-ct*((alp-h2*(2*alp/x+alp_p))*((j+2.)*(J+abs(mmm)+1)/(2*J+3.)/(J+1.)
     & -c3*(j**2+2*j-2)*(j+2)*3/(2*j+5.)/(2*j+3.)/(2*j-1.)))
c
      gaqm1(i,j)=-ct*((alp+h2*(2*alp/x+alp_p))*((j-1.)*(J-abs(mmm))/(2*j-1.)/J
     & -c3*(j**2-3.)*(j-1.)*3/(2*j-3.)/(2*j+3.)/(2*j-1.))) 
c

      gaqp1(i,j)=-ct*((alp+h2*(2*alp/x+alp_p))*((j+2.)*(J+abs(mmm)+1)/(2*j+3.)/(J+1)
     & -c3*(j**2+2*j-2)*(j+2)*3/(2*j+5.)/(2*j+3.)/(2*j-1.)))
c
      beqm1(i,j)=
     & -ct*(-(2*alp+hh**2*(alp*j*(j-1)/x**2-alp_p/x))*((j-1.)*(J-abs(mmm))/(2*j-1.)/J
     & +(-c3)* (j**2-3)*(j-1.)*3/(2*j-3.)/(2*j+3.)/(2*j-1.))
     & -hh**2*(alp/x/x)*((j-1.)*(J-abs(mmm))/(2*j-1.)
     & -3*c3*((j**2+j-3)*J*(j-1.)/(2*j-3.)/(2*j+3.)/(2*j-1))))
c
      beqp1(i,j)=-ct*(-(2*alp+hh**2*(alp*(j+1)*(j+2.)/x**2-alp_p/x))
     & *((j+2.)*(j+abs(mmm)+1)/(2*j+3.)/(J+1.)
     & +(-c3)*(j**2+2*j-2)*(j+2)*3/(2*j+5.)/(2*j+3.)/(2*j-1.))
     & -hh**2*(alp/x/x)*(-(j+2.)*(J+abs(mmm)+1)/(2*j+3.)
     & -3*c3*(-(j**2+j-3)*(j+1.)*(j+2)/(2*j+5.)/(2*j-1.)/(2*j+3.))))   


      alqm3(i,j) = 
     & -ct*(-c3*(alp-h2*(2*alp/x+alp_p))
     & *((j-1)*(j-2)*(j-3)/(2*j-5.)/(2*j-3.)/(2*j-1.)))
 
      alqp3(i,j) = 
     & -ct*(-c3*(alp-h2*(2*alp/x+alp_p))
     & *((j+2)*(j+3)*(j+4)/(2*j+3.)/(2*j+5.)/(2*j+7.)))

      gaqm3(i,j) = 
     & -ct*(-c3*(alp+h2*(2*alp/x+alp_p))
     & *((j-1)*(j-2)*(j-3)/(2*j-5.)/(2*j-3.)/(2*j-1.)))

      gaqp3(i,j) = 
     & -ct*(-c3*(alp+h2*(2*alp/x+alp_p))
     & *((j+2)*(j+3)*(j+4)/(2*j+3.)/(2*j+5.)/(2*j+7.)))

      beqm3(i,j) = 
     & -ct*(-(2*alp+hh**2*(j*(j-1)/x**2-alp_p/x))*(
     & -c3)*((j-1)*(j-2)*(j-3)/(2*j-5.)/(2*j-3.)/(2*j-1.)))
     & -ct*(hh**2*(alp/x/x)*3*c3*(j-2)*(j-1)
     & *(j-3.)*(j-2.)/(2*j-3.)/(2*j-1.)/(2*j-5.))

      beqp3(i,j) = 
     & -ct*(-(2*alp+hh**2*((j+1)*(j+2.)/x**2-alp_p/x))*(
     & -c3)*((j+2)*(j+3)*(j+4)/(2*j+3.)/(2*j+5.)/(2*j+7.)))
     & -ct*(hh**2*(alp/x/x)*3*
     & c3*(-(j+3)*(j+2)*(j+3.)*(j+4.)/(2*j+3.)/(2*j+5.)/(2*j+7)))


C
C DEFINE SIGMENO COMPLEX !!!! ..... otherwise it will be wrong.....!!!!
C

c       apu = -3*psi/x/x/rho
       apu = ar
c       apup = -3*psip/x/x/rho +(2+gd)*psi/x/x 
       apup =  arp
c       bapu = -3*psip/x/rho 
       bapu =  bt
c       bapup = -3*psipp/x/x/rho + (2+gd)*3*psip/x/x/x/rho
       bapup  = btp

c  check sign with the m-independent terms 
       sigmeno  = -c_u*(0,1)*mmm*(J-am)/J/(J+1.)/(2*J-1.)
       sigpiu   = -c_u*(0,1)*mmm*(J+am+1)/J/(J+1.)/(2*J+3.)

       sigmami  =  12*apu/x+6*apup-2*J*bapu/x  ! ok

c coefficient of the non-derivative term CHANGED of sign - the overall - goes in the definition of "beta" below -
       deltami  =  6*apu*j*(J+1)/x/x - 6*apup/x + 2*J*bapu/x/x-(J-1)*
     & J*(6*apu/x/x+bapup/x-bapu/x/x)   !ok 

       sigmapi  = 12*apu/x+6*apup+2*(J+1)*bapu/x  !ok

corg       deltapi  = 6*apu*J*(J+1)/x/x - 6*apup/x + 2*J*bapu/x/x+(J+1)*
corg     & (J+2)*(6*apu/x/x+bapup/x-bapu/x/x)

       deltapi  = 6*apu*J*(J+1)/x/x - 6*apup/x - 2*(J+1)*bapu/x/x-(J+1)*
     & (J+2)*(6*apu/x/x+bapup/x-bapu/x/x)  ! ok

       alqm1(I,J) = alqm1(I,J) + sigmeno*(6*apu-h2*sigmami)
       beqm1(I,J) = beqm1(I,J) - sigmeno*(2*6*apu+hh**2*deltami)
       gaqm1(I,J) = gaqm1(I,J) + sigmeno*(6*apu+h2*sigmami)
       
       alqp1(I,J) = alqp1(I,J) + sigpiu*(6*apu-h2*sigmapi)
       beqp1(I,J) = beqp1(I,J) - sigpiu*(2*6*apu+hh**2*deltapi)
       gaqp1(I,J) = gaqp1(I,J) + sigpiu*(6*apu+h2*sigmapi)
 
 210  continue

c outer boundary conditions
c you now have two helmholtz equations with different parameters
c 

           call bess(beta, J, ffc)
           bea(NP,J) = bea(NP,J) - 4*gaa(NP,J)/(2*hh*(J+1+ffc )+3)
           ala(NP,J) = ala(NP,J) -  gaa(NP,J)/(2*hh*(J+1+ffc )+3)     

c set outer boundary conditions for b-vectors
 
           if(beta.ne.0)then 

           call bess(betb, J, ffc)

          if(betb.gt.10) ffc=betb

          beb(NP,J) = beb(NP,J) - 4*gab(NP,J)/(2*hh*(J+1+ffc)+3)
          alb(NP,J) = alb(NP,J) -  gab(NP,J)/(2*hh*(J+1+ffc)+3)     

           endif
 211  continue

c     initialize the entries in the matrix 
    
          AXX = 0 

C -------------------------------------------

       if(ans2.eq.'D')then       ! VECTOR IS a1 b2 a3 b4 ...

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

       else if(ans2.eq.'Q')then ! VECTOR IS b1 a2 b3 a4 ...

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

C---- LOAD MATRIX ENTRIES 
C runs col 2*nb  ! INNER 

      do 62 j = 1, 3*nb, nb
      nas = 0
      nas_count = 0
      do 52 k1 = nai, naf, 2   
      nas_count = nas_count+1
      nas = nas_count-1	
      nas=2*nas+qa
      if(j/nb.eq.0)AXX(k1+1,j+k1) = -bea(1,nas)
      if(j/nb.eq.1)AXX(k1+1,j+k1) = +gaa(1,nas)

      !BN+1
      if(nas.le.nabp1)then
      if(j/nb.eq.0)AXX(1+k1,j+k1+1)= c1(1)*(nas+2.)*(nas+abs(mmm)+1)/(2*nas+3.)/(nas+1.) 
     & -a2*c1(1)*3*(nas+2)*(nas**2+2*nas-2.)
     & /(2*nas+5.)/(2*nas+3.)/(2*nas-1.)
     & +6*c2(1)*(0,1)*mmm*(nas+abs(mmm)+1.0)/nas/(nas+1)/(2*nas+3) 
      endif 
      !BN-1
      if(nas.ge.nabm1)then
      if(j/nb.eq.0)AXX(1+k1,j+k1-1)= c1(1)*(nas-1.)*(nas-abs(mmm))/(2*nas-1.)/nas
     & -a2*c1(1)*3*(nas**2-3.)*(nas-1.)
     & /(2*nas-3.)/(2*nas+3.)/(2*nas-1.)
     &  +6*c2(1)*(0,1)*mmm*(nas-abs(mmm))/nas/(nas+1)/(2*nas-1.0)
      endif
      !BN+3
      if(nas.le.nabp3)then
      if(j/nb.eq.0)AXX(1+k1,j+k1+3)= 
     & -a2*c1(1)*(nas+2.)
     & *(nas+3.)*(nas+4.)/(2*nas+3.)/(2*nas+5.)/(2*nas+7.)
      endif
      !BN-3
      if(nas.ge.nabm3)then
      if(j/nb.eq.0)AXX(1+k1,j+k1-3)= 
     & -a2*c1(1)*(nas-1.)*(nas-2.)*(nas-3.)
     & /(2*nas-1.)/(2*nas-5.)/(2*nas-3.)  
      endif

      !FANM2
      if(nas.ge.nafm2)then
      if(j/nb.eq.0)AXX(1+k1,j+k1-2)=-fanm2(1,nas)
      endif
      !FANP2          
      if(nas.le.nafp2)then
      if(j/nb.eq.0)AXX(1+k1,j+k1+2)=-fanp2(1,nas)
      endif
      !FANM2P 
      if(nas.ge.nafm2)then
      if(j/nb.eq.1)AXX(1+k1,j+k1-2)= -fanm2p(1,nas)
      endif
      !FANP2P 
      if(nas.le.nafp2)then
      if(j/nb.eq.1)AXX(1+k1,j+k1+2)= -fanp2p(1,nas)
      endif

      !FANM4
      if(nas.ge.nafm2+2)then
      if(j/nb.eq.0)AXX(1+k1,j+k1-4)=-fanm4(1,nas)
      endif
      !FANP4          
      if(nas.le.nafp2-2)then
      if(j/nb.eq.0)AXX(1+k1,j+k1+4)=-fanp4(1,nas)
      endif
      !FANM4P 
      if(nas.ge.nafm2+2)then
      if(j/nb.eq.1)AXX(1+k1,j+k1-4)= -fanm4p(1,nas)
      endif
      !FANP4P 
      if(nas.le.nafp2-2)then
      if(j/nb.eq.1)AXX(1+k1,j+k1+4)= -fanp4p(1,nas)
      endif

      !FANM6
      if(nas.ge.nafm2+4)then
      if(j/nb.eq.0)AXX(1+k1,j+k1-6)=-fanm6(1,nas)
      endif
      !FANP6          
      if(nas.le.nafp2-4)then
      if(j/nb.eq.0)AXX(1+k1,j+k1+6)=-fanp6(1,nas)
      endif
      !FANM6P 
      if(nas.ge.nafm2+4)then
      if(j/nb.eq.1)AXX(1+k1,j+k1-6)= -fanm6p(1,nas)
      endif
      !FANP6P 
      if(nas.le.nafp2-4)then
      if(j/nb.eq.1)AXX(1+k1,j+k1+6)= -fanp6p(1,nas)
      endif


 52   CONTINUE     
 62   CONTINUE

C     runs col 2*nb  INNER BN
      do 66 j = 1, 3*nb, nb   
      nas = 0
      nas_count = 0
      do 56 k1 = nbi, nbf, 2 
      nas_count = nas_count+1
      nas = nas_count-1	
      nas=2*nas+qb
      if(j/nb.eq.0)AXX(k1+1,j+k1) =-beb(1,nas) 
     & -bct*4*x1*alb(1,nas)/(2*hh-3*x1)          ! bct=1 perf con         
      if(j/nb.eq.1)AXX(k1+1,j+k1) = gab(1,nas)     ! bct =0 effect. perf.
     & +bct*x1*alb(1,nas)/(2*hh-3*x1)      

      !AN-1     !ALPHA-OMEGA + ALPHA^2
      if(nas.ge.nbam1)then
      if(j/nb.eq.0)AXX(1+k1,j+k1-1) = beqm1(1,nas) + omem1(1, nas)
      if(j/nb.eq.1)AXX(1+k1,j+k1-1) = gaqm1(1,nas) + omem1p(1,nas)
      endif
      !AN+1     !ALPHA-OMEGA + ALPHA^2
      if(nas.le.nbap1)then
      if(j/nb.eq.0)AXX(1+k1,j+k1+1) = beqp1(1,nas) + omep1(1, nas)
      if(j/nb.eq.1)AXX(1+k1,j+k1+1) = gaqp1(1,nas) + omep1p(1,nas)
      endif
      !AN-3     !OMEGA
      if(nas.ge.nbam3)then
      if(j/nb.eq.0)AXX(1+k1,j+k1-3)= beqm3(1,nas) + omem3(1, nas)
      if(j/nb.eq.1)AXX(1+k1,j+k1-3)= gaqm3(1,nas) + omem3p(1,nas)
      endif
      !AN+3     !OMEGA
      if(nas.le.nbap3)then
      if(j/nb.eq.0)AXX(1+k1,j+k1+3)= beqp3(1,nas) + omep3(1, nas)
      if(j/nb.eq.1)AXX(1+k1,j+k1+3)= gaqp3(1,nas) + omep3p(1,nas)
      endif
      !AN-5
      if(nas.ge.nbam3+2)then
      if(j/nb.eq.0)AXX(1+k1,j+k1-5)= omem5(1, nas)
      if(j/nb.eq.1)AXX(1+k1,j+k1-5)= omem5p(1,nas)
      endif
      !AN+3     !OMEGA
      if(nas.le.nbap3-2)then
      if(j/nb.eq.0)AXX(1+k1,j+k1+5)= omep5(1, nas)
      if(j/nb.eq.1)AXX(1+k1,j+k1+5)= omep5p(1,nas)
      endif

      !FBNM4
      if(nas.ge.nbfm2+2)then
      if(j/nb.eq.0)AXX(1+k1,j+k1-4)=-fbnm4(1,nas)  
     & -bct*4*x1*fbnm4p(1,nas)/(2*hh-3.*x1)
      endif  
      !FBNP4
      if(nas.le.nbfp2-2)then
      if(j/nb.eq.0)AXX(1+k1,j+k1+4)=-fbnp4(1,nas)  
     & -bct*4*x1*fbnp4p(1,nas)/(2*hh-3.*x1)
      endif  
      !FBNM4P
      if(nas.ge.nbfm2+2)then
      if(j/nb.eq.1)AXX(1+k1,j+k1-4)=-fbnm4p(1,nas) + 
     & bct*x1*fbnm4p(1,nas)/(2*hh-3.*x1)
      endif  
      !FBNP4P
      if(nas.le.nbfp2-2)then
      if(j/nb.eq.1)AXX(1+k1,j+k1+4)=-fbnp4p(1,nas) + 
     & bct*x1*fbnp4p(1,nas)/(2*hh-3.*x1)
      endif  

      !FBNM6
      if(nas.ge.nbfm2+4)then
      if(j/nb.eq.0)AXX(1+k1,j+k1-6)=-fbnm6(1,nas)  
     & -bct*4*x1*fbnm6p(1,nas)/(2*hh-3.*x1)
      endif  
      !FBNP6
      if(nas.le.nbfp2-4)then
      if(j/nb.eq.0)AXX(1+k1,j+k1+6)=-fbnp6(1,nas)  
     & -bct*4*x1*fbnp6p(1,nas)/(2*hh-3.*x1)
      endif  
      !FBNM6P
      if(nas.ge.nbfm2+4)then
      if(j/nb.eq.1)AXX(1+k1,j+k1-6)=-fbnm6p(1,nas) + 
     & bct*x1*fbnm6p(1,nas)/(2*hh-3.*x1)
      endif  
      !FBNP6P
      if(nas.le.nbfp2-4)then
      if(j/nb.eq.1)AXX(1+k1,j+k1+6)=-fbnp6p(1,nas) + 
     & bct*x1*fbnp6p(1,nas)/(2*hh-3.*x1)
      endif  
      !FBNM2
      if(nas.ge.nbfm2)then
      if(j/nb.eq.0)AXX(1+k1,j+k1-2)=-fbnm2(1,nas)  
     & -bct*4*x1*fbnm2p(1,nas)/(2*hh-3.*x1)
      endif  
      !FBNP2
      if(nas.le.nbfp2)then
      if(j/nb.eq.0)AXX(1+k1,j+k1+2)=-fbnp2(1,nas)  
     & -bct*4*x1*fbnp2p(1,nas)/(2*hh-3.*x1)
      endif  
      !FBNM2P
      if(nas.ge.nbfm2)then
      if(j/nb.eq.1)AXX(1+k1,j+k1-2)=-fbnm2p(1,nas) + 
     & bct*x1*fbnm2p(1,nas)/(2*hh-3.*x1)
      endif  
      !FBNP2P
      if(nas.le.nbfp2)then
      if(j/nb.eq.1)AXX(1+k1,j+k1+2)=-fbnp2p(1,nas) + 
     & bct*x1*fbnp2p(1,nas)/(2*hh-3.*x1)
      endif  

 56   CONTINUE
 66   CONTINUE

      i = nt-nb+1 	
      k =i/nb      
      do 672 j = 1, 3*nb, nb   ! runs col 2*nb  !OUTER AN
      nas=0
      nas_count=0
      do 572 k1 = nai, naf, 2   
      nas_count = nas_count+1
      nas = nas_count-1	
      nas=2*nas+qa
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1)=  ala(k+1,nas)
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1)= -bea(k+1,nas)

      !BN+1
      if(nas.le.nabp1)then
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1+1)= c1(k+1)*(nas+2.)*(nas+abs(mmm)+1)/(2*nas+3.)/(nas+1.)    
     & -a2*c1(k+1)*3*(nas+2)*(nas**2+2*nas-2.)
     & /(2*nas+5.)/(2*nas+3.)/(2*nas-1.)
     & +6*c2(k+1)*(0,1)*mmm*(nas+abs(mmm)+1.0)/nas/(nas+1)/(2*nas+3)
      endif
      !BN-1
      if(nas.ge.nabm1)then
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1-1)= c1(k+1)*(nas-1.)*(nas-abs(mmm))/(2*nas-1.)/nas  
     & -a2*c1(k+1)*3*(nas**2-3.)*(nas-1.)/(2*nas-3.)/(2*nas+3.)/(2*nas-1.)
     &  +6*c2(k+1)*(0,1)*mmm*(nas-abs(mmm))/nas/(nas+1)/(2*nas-1.0)
      endif
      !BN+3
      if(nas.le.nabp3)then
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1+3)= 
     & -a2*c1(k+1)*(nas+2)*(nas+3)*(nas+4)/(2*nas+3.)/(2*nas+5.)/(2*nas+7.)
      endif
      !BN-3
      if(nas.ge.nabm3)then
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1-3)= 
     & -a2*c1(k+1)*(nas-1.)*(nas-2.)*(nas-3.)/(2*nas-1.)/(2*nas-5.)/(2*nas-3.)  
      endif

      !FANM2P 
      if(nas.ge.nafm2)then
      call bess(beta, nas-2, ffc)
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1-2)= -(-fanm2p(k+1,nas)-fanm2p(k+1,nas)/(2*hh*(nas-1.+ffc)+3))
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1-2)= -( fanm2(k+1,nas) +fanm2p(k+1,nas)*4/(2*hh*(nas-1.+ffc)+3))
      endif
      !FANP2P 
      if(nas.le.nafp2)then
      call bess(beta, nas+2, ffc)
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1+2)= -(-fanp2p(k+1,nas)-fanp2p(k+1,nas)/(2*hh*(nas+3.+ffc)+3))
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1+2)= -( fanp2(k+1,nas) +fanp2p(k+1,nas)*4/(2*hh*(nas+3.+ffc)+3))
      endif
CCC             CHECK!!!!!!!! BENE !!!!!!!!!!!!!!!!!

      !FANM4P 
      if(nas.ge.nafm2+2)then
      call bess(beta, nas-2, ffc)
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1-4)=-(-fanm4p(k+1,nas)-fanm4p(k+1,nas)/(2*hh*(nas-1.+ffc)+3))
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1-4)=-( fanm4(k+1,nas) +fanm4p(k+1,nas)*4/(2*hh*(nas-1.+ffc)+3))
      endif
      !FANP4P 
      if(nas.le.nafp2-2)then
      call bess(beta, nas+2, ffc)
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1+4)=-(-fanp4p(k+1,nas)-fanp4p(k+1,nas)/(2*hh*(nas+3.+ffc)+3))
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1+4)=-( fanp4(k+1,nas) +fanp4p(k+1,nas)*4/(2*hh*(nas+3.+ffc)+3))
      endif

      !FANM6P 
      if(nas.ge.nafm2+4)then
      call bess(beta, nas-2, ffc)
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1-6)= -(-fanm6p(k+1,nas)-fanm6p(k+1,nas)/(2*hh*(nas-1.+ffc)+3))
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1-6)= -( fanm6(k+1,nas) +fanm6p(k+1,nas)*4/(2*hh*(nas-1.+ffc)+3))
      endif
      !FANP6P 
      if(nas.le.nafp2-4)then
      call bess(beta, nas+2, ffc)
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1+6)= -(-fanp6p(k+1,nas)-fanp6p(k+1,nas)/(2*hh*(nas+3.+ffc)+3))
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1+6)= -( fanp6(k+1,nas) +fanp6p(k+1,nas)*4/(2*hh*(nas+3.+ffc)+3))
      endif

 572   CONTINUE
 672   CONTINUE

C     !runs col 2*nb ! OUTER BN
      do 676 j = 1, 3*nb, nb   
      nas=0
      nas_count=0
      do 576 k1 = nbi, nbf, 2   
      nas_count = nas_count+1
      nas = nas_count-1	
      nas = 2*nas+qb
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1)=+alb(k+1,nas)
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1)=-beb(k+1,nas)

      !AN-1     !ALPHA-OMEGA + ALPHA^2
      if(nas.ge.nbam1)then
      call bess(beta,nas-1,ffc)
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1-1)=
     & alqm1(k+1,nas)-gaqm1(k+1,nas)/(2*hh*(nas+ffc)+3)-omem1p(k+1,nas)*(1+1/(2*hh*(nas+ffc)+3))
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1-1)=beqm1(k+1,nas)
     & +4*gaqm1(k+1,nas)/(2*hh*(nas+ffc)+3)+omem1(k+1,nas)+4*omem1p(k+1,nas)/(2*hh*(nas+ffc)+3)
      endif

      !AN+1     !ALPHA-OMEGA + ALPHA^2
      if(nas.le.nbap1)then
      call bess(beta,nas+1,ffc)
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1+1)=
     & alqp1(k+1,nas) - gaqp1(k+1,nas)/(2*hh*(nas+2+ffc)+3)-omep1p(k+1,nas)*(1.+1./(2*hh*(nas+2+ffc)+3) )
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1+1)=
     & beqp1(k+1,nas) + 4*gaqp1(k+1,nas)/(2*hh*(nas+2+ffc)+3)+omep1(k+1,nas)+4*omep1p(k+1,nas)/(2*hh*(nas+2+ffc)+3)
      endif

      !AN-3     !ALPHA-OMEGA + ALPHA^2
      if(nas.ge.nbam3)then
      call bess(beta,nas-3,ffc)
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1-3)= 
     & -omem3p(k+1,nas)*(1.+1./(2*hh*(nas-2+ffc)+3))+alqm3(k+1,nas)*(1.+1./(2*hh*(nas-2+ffc)+3))             
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1-3)= 
     & +omem3(k+1, nas)+4*omem3p(k+1,nas)/(2*hh*(nas-2+ffc)+3)+beqm3(k+1, nas)-4*alqm3(k+1,nas)/(2*hh*(nas-2+ffc)+3)    
      endif

      !AN+3     !ALPHA-OMEGA + ALPHA^2
      if(nas.le.nbap3)then
      call bess(beta,nas+3,ffc)
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1+3)= 
     & -omep3p(k+1,nas)*(1.+1./(2*hh*(nas+4+ffc)+3))+alqp3(k+1,nas)*(1.+1./(2*hh*(nas+4+ffc)+3))             
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1+3)= 
     & +omep3(k+1, nas)+4*omep3p(k+1,nas)/(2*hh*(nas+4+ffc)+3)
     & +beqp3(k+1, nas)-4*alqp3(k+1,nas)/(2*hh*(nas+4+ffc)+3)     
      endif

      !AN-5     !ALPHA-OMEGA 
      if(nas.ge.nbam3+2)then
      call bess(beta,nas-3,ffc)
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1-5)=-omem5p(k+1,nas)*(1.+1./(2*hh*(nas-2+ffc)+3) )
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1-5)=+omem5(k+1,nas)+4*omem5p(k+1,nas)/(2*hh*(nas-2+ffc)+3)
      endif

      !AN+5     !ALPHA-OMEGA 
      if(nas.le.nbap3-2)then
      call bess(beta,nas+3,ffc)
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1+5)=-omep5p(k+1,nas)*(1.+1./(2*hh*(nas+4+ffc)+3) )
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1+5)=+omep5(k+1,nas)+4*omep5p(k+1,nas)/(2*hh*(nas+4+ffc)+3)
      endif


      !FBNM2
      if(nas.ge.nbfm2)then
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1-2)=-fbnm2(k+1,nas) 
      endif  
      !FBNP2
      if(nas.le.nbfp2)then
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1+2)=-fbnp2(k+1,nas) !HERE
      endif  
      !FBNM2P
      if(nas.ge.nbfm2)then
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1-2)=fbnm2p(k+1,nas) 
      endif  
      !FBNP2P
      if(nas.le.nbfp2)then
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1+2)=fbnp2p(k+1,nas) 
      endif  

      !FBNM4
      if(nas.ge.nbfm2+2)then
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1-4)=-fbnm4(k+1,nas) 
      endif  
      !FBNP4
      if(nas.le.nbfp2-2)then
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1+4)=-fbnp4(k+1,nas)
      endif  
      !FBNM4P
      if(nas.ge.nbfm2+2)then
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1-4)=fbnm4p(k+1,nas) 
      endif  
      !FBNP4P
      if(nas.le.nbfp2-2)then
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1+4)=fbnp4p(k+1,nas) 
      endif  

      !FBNM6
      if(nas.ge.nbfm2+4)then
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1-6)=-fbnm6(k+1,nas) 
      endif  
      !FBNP6
      if(nas.le.nbfp2-4)then
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1+6)=-fbnp6(k+1,nas)
      endif  
      !FBNM6P
      if(nas.ge.nbfm2+4)then
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1-6)= fbnm6p(k+1,nas) 
      endif  
      !FBNP6P
      if(nas.le.nbfp2-4)then
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1+6)= fbnp6p(k+1,nas) 
      endif  


 576   CONTINUE
 676   CONTINUE       

C-------     write the inner blocks of the matrix------- AN
    
      do 43 i = nb+1, NT-nb, nb   ! start at each nb row 
      k = i/nb                    ! count the mesh-point
      do 42 j = 1, 3*nb, nb       ! runs col 
      nas_count = 0
      nas = 0
      do 39 k1 = nai, naf, 2  
      nas_count = nas_count+1
      nas = nas_count-1	
      nas=2*nas+qa

      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1)=  ala(k+1,nas)
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1)= -bea(k+1,nas)
      if(j/nb.eq.2)AXX(i+k1,j+(k-1)*nb+k1)=  gaa(k+1,nas)

      !BN+1
      if(nas.le.nabp1)then
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1+1)= 
c     4.15b
     & c1(k+1)*(nas+2.)*(nas+abs(mmm)+1)/(2*nas+3.)/(nas+1.)    
     &-a2*c1(k+1)*3*(nas+2)*(nas**2+2*nas-2.)
     & /(2*nas+5.)/(2*nas+3.)/(2*nas-1.)
     & +6*c2(k+1)*(0,1)*mmm*(nas+abs(mmm)+1.0)/nas/(nas+1)/(2*nas+3) ! 4.11B ! check sign

      endif
      !BN-1
      if(nas.ge.nabm1)then
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1-1)= c1(k+1)*(nas-1.)*(nas-abs(mmm))/(2*nas-1.)/nas    
     &-a2*c1(k+1)*3*(nas**2-3.)*(nas-1.)/(2*nas-3.)/(2*nas+3)/(2*nas-1.)  !CHECK
     &  +6*c2(k+1)*(0,1)*mmm*(nas-abs(mmm))/nas/(nas+1)/(2*nas-1.0)

      endif
      !BN+3
      if(nas.le.nabp3)then
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1+3)= 
     & -a2*c1(k+1)*(nas+2.)*(nas+3.)*(nas+4.)
     & /(2*nas+3.)/(2*nas+5.)/(2*nas+7.)
      endif
      !BN-3
      if(nas.ge.nabm3)then
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1-3)= 
     & -a2*c1(k+1)*(nas-1.)*(nas-2.)*(nas-3.)/(2*nas-1.)
     & /(2*nas-5.)/(2*nas-3.)  
      endif

      !FANM2
      if(nas.ge.nafm2)then
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1-2)=-fanm2(k+1,nas)
      endif
      !FANP2          
      if(nas.le.nafp2)then
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1+2)=-fanp2(k+1,nas)
      endif
      !FANM2P 
      if(nas.ge.nafm2)then
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1-2)= fanm2p(k+1,nas)
      if(j/nb.eq.2)AXX(i+k1,j+(k-1)*nb+k1-2)=-fanm2p(k+1,nas)
      endif
      !FANP2P 
      if(nas.le.nafp2)then
      if(j/(nb).eq.0)AXX(i+k1,j+(k-1)*nb+k1+2)=  fanp2p(k+1,nas)
      if(j/(nb).eq.2)AXX(i+k1,j+(k-1)*nb+k1+2)= -fanp2p(k+1,nas)
      endif

      !FANM4
      if(nas.ge.nafm2+2)then
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1-4)=-fanm4(k+1,nas)
      endif
      !FANP4          
      if(nas.le.nafp2-2)then
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1+4)=-fanp4(k+1,nas)
      endif
      !FANM4P 
      if(nas.ge.nafm2+2)then
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1-4)=+fanm4p(k+1,nas)
      if(j/nb.eq.2)AXX(i+k1,j+(k-1)*nb+k1-4)=-fanm4p(k+1,nas)
      endif
      !FANP4P 
      if(nas.le.nafp2-2)then
      if(j/(nb).eq.0)AXX(i+k1,j+(k-1)*nb+k1+4)=+fanp4p(k+1,nas)
      if(j/(nb).eq.2)AXX(i+k1,j+(k-1)*nb+k1+4)=-fanp4p(k+1,nas)
      endif

      !FANM6
      if(nas.ge.nafm2+4)then
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1-6)=-fanm6(k+1,nas)
      endif
      !FANP6          
      if(nas.le.nafp2-4)then
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1+6)=-fanp6(k+1,nas)
      endif
      !FANM6 
      if(nas.ge.nafm2+4)then
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1-6)=+fanm6p(k+1,nas)
      if(j/nb.eq.2)AXX(i+k1,j+(k-1)*nb+k1-6)=-fanm6p(k+1,nas)
      endif
      !FANP6P 
      if(nas.le.nafp2-4)then
      if(j/(nb).eq.0)AXX(i+k1,j+(k-1)*nb+k1+6)=+fanp6p(k+1,nas)
      if(j/(nb).eq.2)AXX(i+k1,j+(k-1)*nb+k1+6)=-fanp6p(k+1,nas)
      endif


 39   CONTINUE
 42   CONTINUE
 43   CONTINUE



C-------     write the inner blocks of the matrix------- BN

      do 243 i = nb+1, NT-nb, nb   ! start at each 2*nb row 
      k = i/(nb)                   ! count the mesh-point
      do 242 j = 1, 3*nb, nb       ! runs col 2*nb*i
      nas_count = 0
      nas =0
      do 239 k1 = nbi, nbf, 2   
      nas_count = nas_count+1
      nas = nas_count-1	
      nas=2*nas+qb

check segno di beta nel boundary!!!!!!!!!

      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1) =  alb(k+1,nas)
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1) = -beb(k+1,nas)
      if(j/nb.eq.2)AXX(i+k1,j+(k-1)*nb+k1) =  gab(k+1,nas)

      !AN-1     !ALPHA-OMEGA + ALPHA^2
      if(nas.ge.nbam1)then
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1-1)= alqm1(k+1,nas) - omem1p(k+1,nas)
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1-1)= beqm1(k+1,nas) + omem1(k+1, nas)
      if(j/nb.eq.2)AXX(i+k1,j+(k-1)*nb+k1-1)= gaqm1(k+1,nas) + omem1p(k+1,nas)
      endif
      !AN+1     !ALPHA-OMEGA + ALPHA^2
      if(nas.le.nbap1)then
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1+1)= alqp1(k+1,nas) - omep1p(k+1,nas)     
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1+1)= beqp1(k+1,nas) + omep1(k+1, nas)
      if(j/nb.eq.2)AXX(i+k1,j+(k-1)*nb+k1+1)= gaqp1(k+1,nas) + omep1p(k+1,nas)
      endif
      !AN-3     !ALPHA-OMEGA + ALPHA^2
      if(nas.ge.nbam3)then
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1-3)= alqm3(k+1,nas) - omem3p(k+1,nas)
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1-3)= beqm3(k+1,nas) + omem3(k+1, nas)
      if(j/nb.eq.2)AXX(i+k1,j+(k-1)*nb+k1-3)= gaqm3(k+1,nas) + omem3p(k+1,nas)
      endif
      !AN+3     !ALPHA-OMEGA + ALPHA^2
      if(nas.le.nbap3)then
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1+3)= alqp3(k+1,nas) - omep3p(k+1,nas)     
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1+3)= beqp3(k+1,nas) + omep3(k+1, nas)
      if(j/nb.eq.2)AXX(i+k1,j+(k-1)*nb+k1+3)= gaqp3(k+1,nas) + omep3p(k+1,nas)
      endif
      !AN-5     !ALPHA-OMEGA + ALPHA^2
      if(nas.ge.nbam3+2)then
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1-5)= - omem5p(k+1,nas)
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1-5)= + omem5(k+1, nas)
      if(j/nb.eq.2)AXX(i+k1,j+(k-1)*nb+k1-5)= + omem5p(k+1,nas)
      endif
      !AN+5     !ALPHA-OMEGA + ALPHA^2
      if(nas.le.nbap3-2)then
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1+5)= - omep5p(k+1,nas)     
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1+5)= + omep5(k+1, nas)
      if(j/nb.eq.2)AXX(i+k1,j+(k-1)*nb+k1+5)= + omep5p(k+1,nas)
      endif

      !FBNM2
      if(nas.ge.nbfm2)then
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1-2)=-fbnm2(k+1,nas)
      endif  
      !FBNP2
      if(nas.le.nbfp2)then
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1+2)=-fbnp2(k+1,nas) !HERE
      endif  
      !FBNM2P
      if(nas.ge.nbfm2)then
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1-2)= fbnm2p(k+1,nas)
      if(j/nb.eq.2)AXX(i+k1,j+(k-1)*nb+k1-2)=-fbnm2p(k+1,nas)
      endif  
      !FBNP2P
      if(nas.le.nbfp2)then
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1+2)= fbnp2p(k+1,nas)  
      if(j/nb.eq.2)AXX(i+k1,j+(k-1)*nb+k1+2)=-fbnp2p(k+1,nas) 
      endif  

      !FBNM4
      if(nas.ge.nbfm2+2)then
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1-4)=-fbnm4(k+1,nas)
      endif  
      !FBNP4
      if(nas.le.nbfp2-2)then
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1+4)=-fbnp4(k+1,nas)
      endif  
      !FBNM4P
      if(nas.ge.nbfm2+2)then
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1-4)= fbnm4p(k+1,nas)
      if(j/nb.eq.2)AXX(i+k1,j+(k-1)*nb+k1-4)=-fbnm4p(k+1,nas)
      endif  
      !FBNP4P
      if(nas.le.nbfp2-2)then
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1+4)= fbnp4p(k+1,nas) 
      if(j/nb.eq.2)AXX(i+k1,j+(k-1)*nb+k1+4)=-fbnp4p(k+1,nas) 
      endif  

      !FBNM6
      if(nas.ge.nbfm2+4)then
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1-6)=-fbnm6(k+1,nas)
      endif  
      !FBNP6
      if(nas.le.nbfp2-4)then
      if(j/nb.eq.1)AXX(i+k1,j+(k-1)*nb+k1+6)=-fbnp6(k+1,nas)
      endif  
      !FBNM6P
      if(nas.ge.nbfm2+4)then
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1-6)= fbnm6p(k+1,nas)
      if(j/nb.eq.2)AXX(i+k1,j+(k-1)*nb+k1-6)=-fbnm6p(k+1,nas)
      endif  
      !FBNP6P
      if(nas.le.nbfp2-4)then
      if(j/nb.eq.0)AXX(i+k1,j+(k-1)*nb+k1+6)= fbnp6p(k+1,nas) 
      if(j/nb.eq.2)AXX(i+k1,j+(k-1)*nb+k1+6)=-fbnp6p(k+1,nas) 
      endif  

 239   CONTINUE
 242   CONTINUE
 243   CONTINUE

c----- LAPACK solvers
c       stop  

        if(mmm.ne.0) flg=1

        if(flg.eq.1) then
       CALL ZGEEV(JOBVL, JOBVR, NT, AXX, NT, WW, CVR, NT, CVR, NT,
     &  WORK, LWORK, RWORK, INFO )

c------- rescale the eigenvalues 

      do 121 i=1, NT
      wr(i) = REAL(WW(i))/hh**2    
      wi(i) = AIMAG(WW(i))/hh**2
      INDE(i) = i
 121  continue

       else

       forall (i=1:nt,j=1:nt) axr(i,j) = real(axx(i,j))
       axR = real(axx)

C CHECK CALL VR or VC, it should be OK, but maybe better use another
C name?

       CALL DGEEV(JOBVL, JOBVR, NT, AXR, NT, WR, WI, VR, NT, 
     &  VR, NT, WORK, LWORK, INFO )
C
      do 122 i=1, NT
      wr(i) = Wr(i)/hh**2    
      wi(i) = Wi(i)/hh**2
      INDE(i) = i
 122   continue     

       endif

 
c-------  order in increasing order for the real part
C------   vedi se esiste una nuova subr in F90 

      call sort2(NT,wr,INDE)     

      INDEG=INT(INDE)

c      wr(NT) = wr(NT)-10
    
      rate = wr(NT) 

c      rate = wr(NT-1)
c      if(abs(wr(NT-1)) .eq. abs(wr(NT))) rate=wr(NT-2)  
      
      imag = wi(INT(INDEG(NT)))
 
c       print*, wi(INT(INDEG(NT))), wi(INT(INDEG(NT-1)))
c       print*, wr(nt), wr(nt-1)

        do i2 =1,NB
        k = 1          
        k = 0          
        k = -1         
        do i3 = i2,NT,NB
        k = k + 1   

c        write(13,'(I4, 3e15.5)') i2, x1+(k+1)*hh, VR(i3,INDEG(NT)), VR(i3,INDEG(NT-1)) 

        ! check sort2 - 
c        do ii2=1,NT 
c        write(44,*) ii2, INDEG(ii2), WR(ii2)
c        enddo
c        stop
c        write(*,*) nt, indeg(nt)
c        stop
c        print*, INFO, WORK(1), max(1,LWORK)

        EGR(k+1,i2) = VR(i3,INDEG(NT))
        EGI(k+1,i2) = VR(i3,INDEG(NT-1)) 
        enddo
        enddo

         ki=0            
         do i=1,na,2 
         ki=ki+1 
         do j = 1,NP
         AGR(j,ki) = EGR(j,i)
         AGI(j,ki) = EGI(j,i)
         enddo
         enddo

         ki=0
         do i=2,nb,2 
         ki=ki+1
         do j = 1,NP
         BGR(j,ki) = EGR(j,i)
         BGI(j,ki) = EGI(j,i)
         enddo
         enddo

c max_Btor / max_Bpol

         eep=(maxval(abs(bgr)) + maxval(abs(bgi)))/(maxval(abs(agr)) + maxval(abs(agi)))

c         write(*,*) maxval(abs(bgr))
c         write(*,*) maxval(abs(agr))

         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       rt = rate  
       vtu = turb 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C------------ write the period in years ---------------------------

      facp= 2.e0*PI*(SR*SR0)**2.e0/ETA/3.1536e7
c      facp= SR**2.e0/eta/3.1536e7
c        write(*,*) (REAL(WW(J))/hh**2,AIMAG(WW(J))/hh**2,J=1,NT) 
c      print*, facp 
c      do kk = 0,9
      do kk = 0,5
      WRITE(*,'(I4, 1P,3e13.5)') kk, WR(NT-KK), WI(INT(INDEG(NT-KK))),   facp/WI(INT(INDEG(NT-KK)))
      enddo

      write(*,*)

      do kk = 0,9
        REG(kk+1) =  wr(nt-kk)
        IEG(kk+1) =  wi(INT(INDEG(nt-kk)))
      enddo
 
      if(JOBVR.eq.'V')then
c      write(*,*) 'time elapsed before writing ', secnds(t1)
      call writefield
      endif


      return
      END

