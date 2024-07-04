c  

      subroutine stream(x, ax, axp,  bx, bxp, fx)
      IMPLICIT REAL(A-H,O-Z)

      ! STREAM FUNCTION AND DENSITY / normalized
      ! OUTPUT 

      include 'cio' 

      REAL x_in, bct, c3, mmm,SR,rotp,gd, aqu, flg
      COMMON/ppar/x_in,bct,c3,mmm,SR,rotp,gd, aqu, flg

      REAL xa1,xa2,xa3,xb,xda1,xda2
      COMMON/apar/xa1,xa2,xa3,xb,xda1,xda2
    
      REAL s0,s2,s4,s6,A2P,A4P,xm
      COMMON/PSI/s0,s2,s4,s6,A2P,A4P,xm

      REAL ur, ut
      REAL rho     
      REAL xi, xo, xd, y1, y2
      REAL NORM, M1,M2,M3,M4
      
      character*52 dir
      CHARACTER*8 ANS1,ANS2,ANS3,ANS4
      COMMON/var3/ANS1,ANS2,ANS3,ANS4,dir
 
      REAL X,  me
c 
c
C STREAM FUNCTION DEFINITION AND NORMALIZATION  (one cell in latitude)
C
C       P2  = (3 cos(theta)**2 - 1) / 2
C
C       P12   = 3 cos(theta) sin(theta)  = - Diff (P2, theta) 
C     
C       Define the stream function as:        Psi = fx * P12     
C
C       Up = - Curl (0,0,Psi) = (3 fx (1 - 3 cos^2) /x , 3 (fx x)'/x  sin cos, 0) =  (ax * (1- 3 cos^2) , bx sin cos, 0) = input code!
C       
C       so that Div (Up) = 0
C
C       But mass conservation implies  Div ( rho Up ) = 0
C   
C       thus :    
C    
C       Up  ->  Up/rho   
c   
c       Ur =  -(3*cos(x)^2-1) fx / x   =>  ax = fx/x/rho       see eq.10a
c       
c       Utheta = (x fx)'/x  cos sin     =>  bx = (x fx )'  / x / rho   see eq.10b
c
C       so that Div ( rho Up ) = 0
C
c       Bullard-Gellard formalism:
c
c       Curl Curl (s2 P2, 0, 0) = (- 3 fx/x^2 * (3 cos^2 -1 ), -3 fx'/x * sin cos, 0)  
c
c NEW PSI
c         goto 3232


       if(x.ge.xb)then   !

        sigma=s2
        x0=s0

        fx=(1 - dexp(-(x - xb)**2/sigma**2))*(-1 + x)*x**2

        ax=-(((-1 + dexp((x - xb)**2/sigma**2))*(-1 + x)*(1 - x0)**gd)/
     -    (dexp(((-1 + x)*(1 + x - 2*xb))/sigma**2)*
     -      (-1 + dexp((-1 + xb)**2/sigma**2))*(1/x - x0)**gd))

        axp=((1 - x0)**gd*(-((-1 + dexp((x - xb)**2/sigma**2))*gd*sigma**2) + 
     -      x*((-1 + dexp((x - xb)**2/sigma**2))*(1 + gd)*sigma**2 + 2*xb - 
     -         x*(2 + (-1 + dexp((x - xb)**2/sigma**2))*sigma**2*x0 + 2*x**2*x0 + 
     -            2*(1 + x0)*xb - 2*x*(1 + x0 + x0*xb)))))/
     -  (dexp(((-1 + x)*(1 + x - 2*xb))/sigma**2)*
     -    (-1 + dexp((-1 + xb)**2/sigma**2))*sigma**2*x*(1/x - x0)**gd*(-1 + x*x0))


        bx=-((dexp((1 - xb)**2/sigma**2 - (x - xb)**2/sigma**2)*(1 - x0)**gd*
     -      (2*sigma**2 - 3*sigma**2*x - 2*x**2 + 2*x**3 + 
     -        dexp((x - xb)**2/sigma**2)*sigma**2*(-2 + 3*x) + 2*x*xb - 2*x**2*xb))
     -     /((-sigma**2 + dexp((1 - xb)**2/sigma**2)*sigma**2)*(1/x - x0)**gd))

       bxp= ((1 - x0)**gd*(-2*(-1 + dexp((x - xb)**2/sigma**2))*gd*sigma**4 + 
     -      4*x**6*x0 + sigma**2*x*
     -       (3*(-1 + dexp((x - xb)**2/sigma**2))*(1 + gd)*sigma**2 + 
     -         2*(3 + gd)*xb) - 4*x**5*(1 + x0 + 2*x0*xb) + 
     -      x**2*(-2*(4 + gd)*sigma**2 - 
     -         3*(-1 + dexp((x - xb)**2/sigma**2))*sigma**4*x0 - 
     -         2*sigma**2*(5 + gd + 3*x0)*xb + 4*xb**2) + 
     -      2*x**3*(-2*xb*(2 + xb + x0*xb) + 
     -         sigma**2*(6 + gd + 4*x0 + 5*x0*xb)) + 
     -      4*x**4*(1 - 3*sigma**2*x0 + xb*(2 + x0*(2 + xb)))))/
     -  (dexp(((-1 + x)*(1 + x - 2*xb))/sigma**2)*
     -    (-1 + dexp((-1 + xb)**2/sigma**2))*sigma**4*x*(1/x - x0)**gd*(-1 + x*x0))


c IMPORTANT you multiply by 2 because you then divide by sin_theta cos_theta 
           fx=fx*2
           ax=ax*2
           axp=axp*2
           bx = bx*2
           bxp=bxp*2

c          if(x.gt.0.95) ax=0
c          if(x.gt.0.95) axp=0
c          if(x.gt.0.95) bx=0
c          if(x.gt.0.95) bxp=0


       else 
        fx=0
        ax=0
        axp=0
        bx=0
        bxp=0
       endif

      return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
      include 'b2func'
      include 'b3func'
      include 'b4func'
      include 'b5func'
      include 'c1s2func'
      include 'c3s2func'      
      include 'cdfunc'




