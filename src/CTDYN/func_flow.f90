module func_flow 

  use kind_parameters
  use cio

  implicit none

  private

  public :: stream

contains 

  subroutine stream(x, ax, axp,  bx, bxp, fx)
    !Stream function and density / normalized
    !output 
    !
    !STREAM FUNCTION DEFINITION AND NORMALIZATION  (one cell in latitude)
    !
    !P2  = (3 cos(theta)**2 - 1) / 2
    !
    !P12   = 3 cos(theta) sin(theta)  = - Diff (P2, theta) 
    !
    !Define the stream function as:        Psi = fx * P12     
    !
    !Up = - Curl (0,0,Psi) = (3 fx (1 - 3 cos^2) /x , 3 (fx x)'/x  sin cos, 0) =  (ax * (1- 3 cos^2) , bx sin cos, 0) = input code!
    !
    !so that Div (Up) = 0
    !
    !But mass conservation implies  Div ( rho Up ) = 0
    !
    !thus :    
    !
    !Up  ->  Up/rho   
    ! 
    !Ur =  -(3*cos(x)^2-1) fx / x   =>  ax = fx/x/rho       see eq.10a
    !
    !Utheta = (x fx)'/x  cos sin     =>  bx = (x fx )'  / x / rho   see eq.10b
    !
    !so that Div ( rho Up ) = 0
    !
    !Bullard-Gellard formalism:
  
    !Curl Curl (s2 P2, 0, 0) = (- 3 fx/x^2 * (3 cos^2 -1 ), -3 fx'/x * sin cos, 0)  
  
    real(dp) :: ur, ut, ax, axp, bx, bxp, fx, sigma, x0
    real(dp) :: rho     
    real(dp) :: xi, xo, xd, y1, y2
    real(dp) :: norm, m1,m2,m3,m4
    real(dp) :: x, me
  
    if (x .ge. xb) then   
      sigma=s2
      x0=s0
  
      fx=(1 - dexp(-(x - xb)**2/sigma**2))*(-1d0 + x)*x**2
  
      ax=-(((-1 + dexp((x - xb)**2/sigma**2))*(-1d0 + x)*(1d0 - x0)**gd)/ &
        &    (dexp(((-1d0 + x)*(1d0 + x - 2d0*xb))/sigma**2)*  &
        &      (-1d0 + dexp((-1d0 + xb)**2/sigma**2))*(1d0/x - x0)**gd)) 
  
      axp=((1 - x0)**gd*(-((-1d0 + dexp((x - xb)**2/sigma**2))*gd*sigma**2) + & 
        &      x*((-1d0 + dexp((x - xb)**2d0/sigma**2))*(1d0 + gd)*sigma**2 + 2d0*xb - &
        &         x*(2d0 + (-1d0 + dexp((x - xb)**2d0/sigma**2))*sigma**2*x0 + 2d0*x**2*x0 + & 
        &            2d0*(1d0 + x0)*xb - 2d0*x*(1 + x0 + x0*xb)))))/ &
        &  (dexp(((-1d0 + x)*(1d0 + x - 2d0*xb))/sigma**2)* &
        &    (-1d0 + dexp((-1d0 + xb)**2/sigma**2))*sigma**2*x*(1d0/x - x0)**gd*(-1d0 + x*x0))
  
  
      bx=-((dexp((1d0 - xb)**2/sigma**2 - (x - xb)**2/sigma**2)*(1d0 - x0)**gd*  &
        &      (2d0*sigma**2 - 3d0*sigma**2*x - 2d0*x**2 + 2d0*x**3 + &
        &        dexp((x - xb)**2d0/sigma**2)*sigma**2*(-2d0 + 3d0*x) + 2d0*x*xb - 2d0*x**2*xb)) &
        &     /((-sigma**2 + dexp((1d0 - xb)**2/sigma**2)*sigma**2)*(1d0/x - x0)**gd))
  
      bxp= ((1d0 - x0)**gd*(-2d0*(-1d0 + dexp((x - xb)**2/sigma**2))*gd*sigma**4 + &
        &      4d0*x**6*x0 + sigma**2*x* &
        &       (3d0*(-1d0 + dexp((x - xb)**2/sigma**2))*(1d0 + gd)*sigma**2 + &
        &         2d0*(3d0 + gd)*xb) - 4d0*x**5*(1d0 + x0 + 2d0*x0*xb) + &
        &      x**2*(-2d0*(4d0 + gd)*sigma**2 - &
        &         3d0*(-1d0 + dexp((x - xb)**2/sigma**2))*sigma**4*x0 - &
        &         2d0*sigma**2*(5d0 + gd + 3d0*x0)*xb + 4d0*xb**2) + &
        &      2d0*x**3*(-2d0*xb*(2d0 + xb + x0*xb) + &
        &         sigma**2*(6d0 + gd + 4d0*x0 + 5d0*x0*xb)) + & 
        &      4d0*x**4*(1d0 - 3d0*sigma**2*x0 + xb*(2d0 + x0*(2d0 + xb)))))/ &
        &  (dexp(((-1d0 + x)*(1d0 + x - 2d0*xb))/sigma**2)* &
        &    (-1d0 + dexp((-1d0 + xb)**2d0/sigma**2))*sigma**4*x*(1d0/x - x0)**gd*(-1d0 + x*x0))
  
      ! important you multiply by 2 because you then divide by sin_theta cos_theta 
      fx=fx*2d0
      ax=ax*2d0
      axp=axp*2d0
      bx = bx*2d0
      bxp=bxp*2d0
  
    else 
      fx=0d0
      ax=0d0
      axp=0d0
      bx=0d0
      bxp=0d0
    endif
    return
  end subroutine

end module func_flow
  
