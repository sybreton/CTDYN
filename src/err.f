
      SUBROUTINE ROT(x, om0, om0p, om2, om2p, om4, om4p)
      IMPLICIT REAL(A-H,O-Z)

      include 'cio'

      REAL x,x1,x2,x3,xd,xd2,xd3
      REAL om0, om2, om0p, om2p, om4, om4p, a1, a2, xmax

      REAL omc,omeq,hzn,c2,c4,omegac,oms,xbi,omegaeq
     
c---------- common main

      REAL xa1,xa2,xa3,xb,xda1,xda2
      COMMON/apar/xa1,xa2,xa3,xb,xda1,xda2

      REAL x_in, bct, c3,mmm,SR,rotp,gd,aqu,flg
      COMMON/ppar/x_in,bct,c3,mmm,SR,rotp,gd,aqu,flg

      REAL dd1,rc1,rc2,oco
      COMMON/dpar/dd1,rc1,rc2,oco

      character*52 dir
      CHARACTER*8 ANS1,ANS2,ANS3,ANS4
      COMMON/var3/ANS1,ANS2,ANS3,ANS4,dir

      REAL s0,s2,s4,s6,A2P,A4P,xm
      COMMON/PSI/s0,s2,s4,s6,A2P,A4P,xm

c--------- definition of omega see Dikpati & Carbonneau Apj99.
               
      facn      =  1 
      omc	=  432.8e0*facn
      omeq   	=  460.7e0*facn
c   from namelist
c      rc1	=  0.7e0
c      dd1	=  2.5e-2
c      dd1	=  0.025 
c      hzn	=  1.0e-9
      hzn	=  1.0

C see Dikpati & Charbonneau

      c2	=  +62.69e0*hzn*2.0*pi*facn
      c4	=  +67.13e0*hzn*2.0*pi

c      c4	=  +63.13e0*hzn*2.0*pi*facn

c from namelist
      omegac	=  2.0e0*pi*omc*hzn
      omegaeq	=  2.0e0*pi*omeq*hzn

         om2    = 0.0e0
         om2p   = 0.0e0
         om4  =   0.0e0
         om4p =   0.0e0

       if(ANS1.eq.'R86') then
C
C eq (9) pag 91, Rad86 AN paper, appendix tables
C-----------  if ANS1 = R86  then we test Radler '86  AN paper --------
C    dw = xda1
C    xw = xda1
C
         dw = xda2
         xw = xa2 

         xi = (x-xw)/dw

         if(x.le.xw-dw)then
         om0 = -1
         om0p = 0
         else if(x.gt.xw-dw.and.x.le.xw+dw)then
         om0 = -(1-3*xi/2+xi**3/2)/2.
         om0p =-(-3+3*xi**2 )/4./dw
         else 
         om0 =0
         om0p =0
         endif

       else if(ANS1.eq.'NS')then
c PNS model
         om0  = 1.+x**2
         om0p = 2*x

c for cylindrical rotation law
         om0  = 1+2*x**2/3
         om0p = 4*x/3 
         om2  = -2*x**2/3 
         om2p = -4*x/3 
         om4  = 0.0e+00
         om4p = 0.0e+00

       else if(ANS1.eq.'DJ')then
c Dudley and James
         om0  = x*sin(x*pi)
         om0p = sin(pi*x) + cos(pi*x)*pi*x

         else if(ANS1.eq.'vega') then
         om0 = x
         om0p = 1.

         else if(ANS1.eq.'Rp72'.OR.ANS1.eq.'Rm72'.or.ans1.eq.'ekeri') then
C Roberts 72 ALPHA OMEGA model 1 eq.4.3

         om0 = (x-x_in)
         om0p = 1.

c         om0 = (x-x_in)**2
c         om0p =2*(x-x_in) 
c      
c         else if(ANS1.eq.'ST76') then
C Stix 76
c        om2  = -2.*x**2/3.
c        om2p = -4.*x/3.
c        om0 = 2.*x**2/3.
c        om0p = 4.*x/3.

         else if(ANS1.eq.'corona') then
 
c        d1 = 0.075
c        r1 = rc1

c        om0 = (1-derf((x-r1)/d1))/2.
c        om0p =-exp(-((x-r1)/d1)**2)/sqrt(PI)/d1

         xx1=0.4
         xx2=0.9
         xxd1=0.1
         xxd2=0.02


       om0 = 1

       om0p= (-2*(1 + Erf((x - xx1)/xxd1)))/(exp((x - xx2)**2/xxd2**2)*Sqrt(Pi)*xxd2) + 
     -  (2*(1 - Erf((x - xx2)/xxd2)))/(exp((x - xx1)**2/xxd1**2)*Sqrt(Pi)*xxd1)

      else if(ANS1.eq.'R72'.or.ANS1.eq.'R72M')then
C Roberts 72 ALPHA OMEGA model 2 eq.4.5  and model 5.6
         om0 = -19683.*(1-x**2)**5/40960.
         om0p = -19683.0*5*(1-x**2)**4*(-2*x)/40960.

         else if(ANS1.eq.'R72B') then
C Roberts 72 ALPHA OMEGA model eq.4.5
         om0 = -3.*sqrt(3.)*(1-x**2)**2/8.
         om0p = 3.*sqrt(3.)*(1-x**2)*4*x/8.

         else if(ANS1.eq.'BR') then
c
C  Roberts 72 ALPHA OMEGA model 5.2  (Bragisnki)
c         om0 = (3.*sqrt(3.)/8)*(-x**2*(1-x**2)**2)
c         om0p = (3.*sqrt(3.)/8)*(-2*x*(1-x**2)**2 -x**2*(1-x**2)*(-2*x)*2)
c
         om0 = -(3.*sqrt(3.)/8)*((1-x**2)**2)
         om0p = (3.*sqrt(3.)/8)*(+2*x*(1-x**2)*2)

         else if(ANS1.eq.'SK69'.or.ANS1.eq.'STM') then
c model 1 steenbeck-krause SK69 paper, ALPHA OMEGA

        d1 = 0.075
        r1 = rc1
        om0 = (1-derf((x-r1)/d1))/2.
        om0p =-exp(-((x-r1)/d1)**2)/sqrt(PI)/d1

         else if(ANS1.eq.'ST74') then
c  Stix 74 paper, ALPHA OMEGA

        d1 = 0.075
        r1 = 0.7
        om0 = (1-derf((x-r1)/d1))/2.
        om0p =-exp(-((x-r1)/d1)**2)/sqrt(PI)/d1

       else if(ans1.eq.'H2')then
C
C OMEGA = 
C       = 
C       = om0    + om2 * P2
C  
C       P2 = -1/2 * (1-3*COS^2 THETA) 
C test benchmark case : check normalization

      c4 =0
c  SUN ! to be input via namelist
c      rc1 = 0.715
c  ior!
c      rc1 = 0.80
c      dd1 = 0.03 
c      dd1 = 0.02 

      omegaeq = 1.
      omegac  = oco
      c2 =   rc2 
      delta   = (1.e0 + derf( (x-rc1)/dd1  ) )/ 2.e0

c c2= -coefficient of cos(theta) ^2 at the surface DR:   Omega_surface=Omega_eq-c2*cos(theta)^2

      om2  = -2*delta*c2/3.d0 
      om2p = (-2*c2/3.d0)*(exp(-((x-rc1)/dd1)**2)/sqrt(PI)/dd1) 

      om0 = delta*(-c2/3.0+1-omegac) +omegac   
      om0p = (-c2/3.0+1-omegac)*(exp(-((x-rc1)/dd1)**2)/sqrt(PI)/dd1)


c K giant model

        else if(ans1.eq.'K2')then
c arcturus !!!!! 
C OMEGA = 
C       = 
C       = om0    + om2 * P2
C  
C       P2 = -1/2 * (1-3*COS^2 THETA) 

c mimic the leonid omega

      om2=0.01/x-0.25
      om2p=-0.01/x/x

      om0  = 0.22/x+1.   
      om0p = -0.22/x/x


      else if(ans1.eq.'H4')then
C
C omega is given by the following expression taken from helioseismology
C OMEGA = omegac + delta*(omegaeq-omegac) - delta*C2*cos(theta)^2 -delta*C4*cos(theta)^4 
C       = OMEGA0 + omega2 * cos(theta)^2 + omega4 *cos(theta)^4
C       = om0    + om2*cos(theta)^2  + om4*cos(theta)^4
C              
      delta   = (1.e0 + derf( (x-rc1)/dd1  ) )/ 2.e0
      om0  =  (omegac + delta * ( omegaeq - omegac))/omegaeq
      om0p =  (exp(-((x-rc1)/dd1)**2)/sqrt(PI)/dd1)*(omegaeq - omegac)/omegaeq
      om2  = -delta*c2/omegaeq
      om2p = -(exp(-((x-rc1)/dd1)**2)/sqrt(PI)/dd1)*c2/omegaeq
C
      om4  = -delta*c4/omegaeq
      om4p =  -(exp(-((x-rc1)/dd1)**2)/sqrt(PI)/dd1)*c4/omegaeq

C compute bench DR with alternative formalism with no P2: the result is the same as with flag H2, just uncomment the lines!
c
c      dd1=0.02
c      rc1=0.7      
c      c4 =0
c
c      omegaeq = 1.
c      omegac  = oco
c      c2 =   rc2 
c
c      delta   = (1.e0 + derf( (x-rc1)/dd1  ) )/ 2.e0
c
c      om0  =  (omegac + delta * ( omegaeq - omegac))/omegaeq
c      om0p =  (exp(-((x-rc1)/dd1)**2)/sqrt(PI)/dd1)*
c     &         (omegaeq - omegac)/omegaeq
c
c      om2  = -delta*(c2)/omegaeq
c      om2p = -(exp(-((x-rc1)/dd1)**2)/sqrt(PI)/dd1)*(c2)/omegaeq
c      om4  =  0 
c      om4p =  0 
c
      else if(ans1.eq.'H5')then !PURE LATITUDINAL DEPENDENCE
C
      om0  =  2.78
      om0  =  1.
      om0p =  0
      om2  =  om0*(-0.13)
      om2p =  0
      om4  =  om0*(-0.16)
      om4p =  0

      else if(ans1.eq.'H6')then  ! SUBSURFACE SHEAR

      ome0 = 435.
      omeq = 452.5
      rt   = 0.69
      rc   = 0.71
      a2   = -61.0
      a4   = -73.5
      omet = 0.05
      omec = 0.05
      omes = 0.05
c      omes = 0.03
      rs   = 0.95
      b0   = 437.
      b4   = -1445./2
      r    =  x 
c
c      b0   = s6
c      b4   = s4
c      rs   = 0.7

       om0 = (b0 + ome0 + omeq - b0*r + 
     -    (-((r - rc)*(7*a2 + 3*a4 + 
     -            35*(b0 - ome0 + omeq - b0*rs))*
     -          derf((2*(r - rc))/omec)) + 
     -       (7*a2 + 3*a4 + 35*(b0 - ome0 + omeq - b0*rc))*
     -        (r - rs)*derf((2*(r - rs))/omes) + 
     -       (7*a2 + 3*a4)*(-rc + rs)*derf((2*(r - rt))/omet))/(35.*(rc - rs)))/2.

      om0p = (-b0 + ((4*(7*a2 + 3*a4 + 
     -            35*(b0 - ome0 + omeq - b0*rc))*(r - rs))/
     -        (dexp((4*(r - rs)**2)/omes**2)*omes*Sqrt(Pi)) + 
     -       (4*(7*a2 + 3*a4)*(-rc + rs))/
     -        (dexp((4*(r - rt)**2)/omet**2)*omet*Sqrt(Pi)) - 
     -       (4*(r - rc)*
     -          (7*a2 + 3*a4 + 35*(b0 - ome0 + omeq - b0*rs))
     -          )/
     -        (dexp((4*(r - rc)**2)/omec**2)*omec*Sqrt(Pi)) - 
     -       (7*a2 + 3*a4 + 35*(b0 - ome0 + omeq - b0*rs))*
     -        derf((2*(r - rc))/omec) + 
     -       (7*a2 + 3*a4 + 35*(b0 - ome0 + omeq - b0*rc))*
     -        derf((2*(r - rs))/omes))/(35.*(rc - rs)))/2.

      om2 = (a2*(1 + derf((2*(r - rt))/omet)))/2.
      om2p= (2*a2)/(dexp((4*(r - rt)**2)/omet**2)*omet*Sqrt(Pi))

      om4=(-(b4*(r - rc)*(1 - rs)*derf((2*(r - rc))/omec)) + 
     -    b4*(1 - rc)*(r - rs)*derf((2*(r - rs))/omes) + 
     -    (rc - rs)*(a4 + b4 - b4*r + a4*derf((2*(r - rt))/omet)))/(2.*(rc - rs))

      om4p=((-4*b4*(r - rc)*(1 - rs))/(dexp((4*(r - rc)**2)/omec**2)*omec*Sqrt(Pi)) + 
     -    (4*b4*(1 - rc)*(r - rs))/(dexp((4*(r - rs)**2)/omes**2)*omes*Sqrt(Pi)) + 
     -    (-b4 + (4*a4)/(dexp((4*(r - rt)**2)/omet**2)*omet*Sqrt(Pi)))*(rc - rs) - 
     -    b4*(1 - rs)*derf((2*(r - rc))/omec) + b4*(1 - rc)*derf((2*(r - rs))/omes))/
     -  (2.*(rc - rs))


        om0   = om0  /omeq
        om0p  = om0p /omeq
        om2   = om2  /omeq
        om2p  = om2p /omeq
        om4   = om4  /omeq
        om4p  = om4p /omeq


      else if(ans1.eq.'R2')then

      om0  = -3.*sqrt(3.)*(1-x**2)**2/8.
      om0p = -(3.*sqrt(3.)/8)*2*(1-x**2)*(-2*x)

      om2=0
      om2p=0
      om4 =0
      om4p =0


      endif
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE ALT(x,a1,a2)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C derf(X-X1) is one for X> X1 and -1 for X< X1   
C 
C by definition we have  derf(x) = 2/sqrt(pi) * int exp(-t^2) t=0..x
C                        derf(x)' =2/sqrt(pi) * exp(-x^2) 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT REAL(A-H,O-Z)

      include 'cio'
      REAL xa1,xa2,xa3,xb,xda1,xda2
      COMMON/apar/xa1,xa2,xa3,xb,xda1,xda2

      REAL x_in, bct, c3,mmm,SR,gd, aqu, flg
      COMMON/ppar/x_in,bct,c3,mmm,SR,gd, aqu, flg

      character*52 dir
      CHARACTER*8 ANS1,ANS2,ANS3,ANS4
      COMMON/var3/ANS1,ANS2,ANS3,ANS4,dir


      if(ANS1.eq.'NS') then

         a1 = (1+derf((x-xa1)/xda1))/2.e0
         a2 = exp(-(x-xa1)**2/xda1/xda1)/sqrt(pi)/xda1

c         a1 = 1.0
c         a2 = 0.0

      else if(ANS1.eq.'R86') then

C   model (i) R86

        d_xi = xda1 
        x_xi = xa1 
        xi = (x-x_xi)/d_xi

         if(x.le.x_xi-d_xi)then
         a1=0
         a2=0
         else if(x.ge.x_xi-d_xi.and.x.le.x_xi+d_xi)then
         a1 = (-15+10*3*xi**2-3*5*xi**4)/16./d_xi
         a2 = (10*3*2*xi-3*5*4*xi**3)/16./d_xi/d_xi
         else if(x.gt.x_xi+d_xi)then
         a1=0
         a2=0
          endif
c      else if(ANS1.eq.'vega') then

c        d1 = 0.025
c        r1 = 0.95
c        a1 = (1-derf((x-r1)/d1))/2.
c        a2 =-exp(-((x-r1)/d1)**2)/sqrt(PI)/d1

c        if(xa2.ne.xa1)then

c         a1 =(1-derf((x-xa2)/xda2) )
c     &   *(1+derf((x-xa1)/xda1))/2.e0/2.e0
c
c
c       a2=((-2 + derf((x - xa1)/xda1))/(dexp((x - xa2)**2/xda2**2)*xda2) + 
c     -    derf((x - xa2)/xda2)/(dexp((x - xa1)**2/xda1**2)*xda1))/(2.*Sqrt(Pi))

c
       else if(ANS1.eq.'AQ1') then
C    simplest alpha2 model 
            a1=1
            a2=0
      else if(ANS1.eq.'AQ2')then 
C     SK' AN 291, p271, 1969
            if(x.lt.xa1)then
             a1=0
             a2=0  
            else
             a1 = (x-xa1)*(x*3.-xa1-2.)*(27./4.)
             a2 = (6.*x -4*xa1 -2.)*27./((1-xa1)**3.0)/4.
            endif

       else if(ANS1.eq.'Rp72')then
c ROBERTS 72 alphaomega model  1 eq. 4.3
          a1= 1
          a2= 0
       else if(ANS1.eq.'Rm72')then
c ROBERTS 72 alphaomega model  1 eq. 4.3
          a1= -1
          a2= 0

      else if(ANS1.eq.'R72'.or.ANS1.eq.'R72M')then
c ROBERTS 72 alphaomega model  2 eq. 4.5 and 5.6

          a1= 729*x**8*(1-x**2)**2/16.0
          a2= 729*(8*x**7*(1-x**2)**2+x**8*2*(1-x**2)*(-2*x))/16.0

      else if(ANS1.eq.'R72B')then
c ROBERTS 72 alphaomega model 2 eq. 4.5

          a1 = 24.*sqrt(3.)*x**2*(1-x)**2
          a2 = 24.*sqrt(3.)*(2*x*(1-x)**2-x**2*2*(1-x))

        else if(ANS1.eq.'ST74') then
        d1 = 0.075
        r1 = 0.7
        a1 = (1-derf((x-r1)/d1))/2.
        a2 =-exp(-((x-r1)/d1)**2)/sqrt(PI)/d1

      ELSE if(ANS1.eq.'H2'.or.ANS1.eq.'H4'.or.ANS1.eq.'K2'
     +     .or.ANS1.eq.'H6'.or.ANS1.eq.'vega'.or.ans1.eq.'ekeri')then
c three steps....
c       xbi = 2.e0*xb/(xb-1.e0)
c      a1 =(1-xbi+derf((x-xa1)/xda1 ))*(1-derf((x-xa3)/xda2))
c     &   *(1+derf((x-xa2)/xda1)) / (2-xbi)/2.e0/2.e0

c      a2 = (1+derf((x-xa2)/xda1))*exp(-((x-xa1)/xda1)**2 )
c     &     *(1-derf((x-xa3)/xda2))/(2-xbi)/2.0/sqrt(pi)/xda1 
c     &     -(1-xbi+derf((x-xa1)/xda1 ))*exp(-((x-xa3)/xda2)**2 )
c     &     *(1+derf((x-xa2)/xda1)) / (2-xbi)/2./sqrt(pi)/xda2
c     &     +(1-xbi+derf((x-xa1)/xda1 ))*exp(-((x-xa2)/xda1)**2 )
c     &     *(1-derf((x-xa3)/xda2) )/ (2-xbi)/2./sqrt(pi)/xda1
c

c  two  error functions
c
c        if(xa2.ne.xa1)then

         a1 =(1-derf((x-xa2)/xda2) )
     &   *(1+derf((x-xa1)/xda1))/2.e0/2.e0
c
c
         a2 = (1-derf((x-xa2)/xda2))
     &         *exp(-(x-xa1)**2/xda1/xda1)/sqrt(pi)/xda1/2
     &        -(1+derf((x-xa1)/xda1))
     &         *exp(-(x-xa2)**2/xda2/xda2)/sqrt(pi)/xda2/2.
c
c       a2=((-2 + derf((x - xa1)/xda1))/(dexp((x - xa2)**2/xda2**2)*xda2) + 
c     -    derf((x - xa2)/xda2)/(dexp((x - xa1)**2/xda1**2)*xda1))/(2.*Sqrt(Pi))
c
c          else
c
c  one erf (alpha not zero after a given xa)
c
c         a1 = (1+derf((x-xa1)/xda1))/2.e0
c         a2 = exp(-(x-xa1)**2/xda1/xda1)/sqrt(pi)/xda1
c
c         endif
c
c
c divide by 3srt(3)/4 to normalize alpha  when sin^2 cos functional form is used

        if(c3.ne.0)then  
        a1 = a1 /0.3849
        a2 = a2 /0.3849
        endif         

         else if(ANS1.eq.'SK69'.or.ANS1.eq.'STM') then
        
c         xa1=0.9
         xda1=0.075

         a1 = (1+derf((x-xa1)/xda1))/2.e0
         a2 = exp(-(x-xa1)**2/xda1/xda1)/sqrt(pi)/xda1

c         a1 = a1*0.11
c         a2 = a2*0.11


      ELSE if(ANS1.eq.'BR')then

         a1 = sqrt(3.)*24*x**2*(1-x)**2
         a2 = sqrt(3.)*24*(2*x*(1-x)**2-x**2*(1-x)*2)

      ELSE if(ANS1.eq.'DJ')then
         a1 = (1+derf((x-xa1)/xda1))/2.e0
         a2 = exp(-(x-xa1)**2/xda1/xda1)/sqrt(pi)/xda1

      ELSE if(ANS1.eq.'H5')then

      a1=1
      a2=0

      ELSE if(ANS1.eq.'corona')then
c
         a1 =(1-derf((x-xa2)/xda2) )
     &   *(1+derf((x-xa1)/xda1))/2.e0/2.e0
c
         a2 = (1-derf((x-xa2)/xda2))
     &         *exp(-(x-xa1)**2/xda1/xda1)/sqrt(pi)/xda1/2
     &        -(1+derf((x-xa1)/xda1))
     &         *exp(-(x-xa2)**2/xda2/xda2)/sqrt(pi)/xda2/2.

c     charbonneau and macgregor 2001 o-b star model
      ELSE IF (ans1.EQ.'C01') THEN

      delta = 1. - derf( (x-xa1)/xda1 )
      a1 = derf(2.*x/xda1)*delta/2.
      a2 = ((2.*EXP(-(2.*x/xda1)**2))/(xda1*SQRT(PI)))*delta  - (EXP(-((x-xa1)/xda1)**2 )/(SQRT(PI)*xda1))*derf(2.*x/xda1) 
      
c      WRITE(444,'(F7.3,5X,F7.3)') x, a1

c     JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
      


       else 

        write(*,*) 'not available!!', ans1

      stop
      endif



      return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE eta_turb(x, e1, e2, e3)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   It gives eta_turb= e1, and derivatives eta' = e2, eta'' = e3 
C   edr is the ratio eta_c/eta_top
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT REAL(A-H,O-Z)

      include 'cio'

      REAL x, x1,x2,xd
      REAL e1,e2,e3

      REAL edr,xe1,xde1,hd
      COMMON/epar/edr,xe1,xde1,hd

      REAL s0,s2,s4,s6,A2P,A4P,xm
      COMMON/PSI/s0,s2,s4,s6,A2P,A4P,xm

      character*53 dir
      CHARACTER*8 ANS1,ANS2,ANS3,ANS4
      COMMON/var3/ANS1,ANS2,ANS3,ANS4,dir

      x1 = xe1
      xd = xde1 

c     ETA = ETA_C+(1/2)*(ETA_T-ETA_C)*(...)
c
c      PSI = 1/((exp(x-x1)/xd)+1)
c      e1 = edr*PSI +1-PSI
c      x3= 0.8e0
c      xd3=0.1
c
      ett=1.e0
      et2= 1.e0 !0.01

      e1 =   edr+(et2-edr)*(1+derf((x-x1)/xd))/2.

      e2 =  (et2-edr)*exp(-((x-x1)/xd)**2)/sqrt(PI)/xd 

      e3 =  -2*(et2-edr)*(x-x1)*exp(-((x-x1)/xd)**2)/sqrt(PI)/xd**3


      if(ANS1.eq.'NS')then
c REVERSE !!!
       e1 = et2+(edr-et2)*(1+derf((x-x1)/xd))/2.
       e2 =  (edr-et2)*exp(-((x-x1)/xd)**2)/sqrt(PI)/xd 
       e3 =  -2*(edr-et2)*(x-x1)*exp(-((x-x1)/xd)**2)/sqrt(PI)/xd**3


c      if(x .ge. 0.9) e1=e1*10.
c      if(x .ge. 0.9) e2=e2*10.

c     charbonneau and macgregor 2001 o-b star model
      ELSE IF (ans1.EQ.'C01') THEN

      et2 = 1.
   
      delta = 1. - derf( (x-x1)/xd )
      e1 = edr + delta*(et2-edr)/2.
      e2 = -(EXP( -((x-x1)/xd)**2 )/(SQRT(PI)*xd))*(et2-edr)
      e3 =(2.*(x-x1)/(SQRT(PI)*(xd)**3))*(EXP( -((x-x1)/xd)**2) )*(et2-edr)
    
c      WRITE(445,'(F7.3,5X,F7.3)') x, e1

      endif

      return
      END

