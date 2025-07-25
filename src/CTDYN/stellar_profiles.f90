module stellar_profiles

  use kind_parameters
  use cio

  implicit none

  private

  public :: rot, alt, eta_turb

contains

  subroutine rot(x, om0, om0p, om2, om2p, om4, om4p)
  
    real(dp) :: x, x1, x2, x3, xd ,xd2, xd3, &
            & om0, om2, om0p, om2p, om4, om4p, a1, a2, xmax, &
            & omc, omeq, hzn, c2, c4, omegac, oms, omegaeq, &
            & a4, b0, b4, d1, delta, dw, facn, ome0, &
            & omec, omes, omet, r, r1, rc, rs, rt, xi, xw, &
            & xx1, xx2, xxd1, xxd2 
         
    
    !--------- definition of omega see Dikpati & Chcarbonneau apj99.
                 
    facn = 1d0 
    omc = 432.8e0*facn
    omeq = 460.7e0*facn
    hzn = 1.0d0
  
    c2 = 62.69e0*hzn*2.0d0*pi*facn
    c4 = 67.13e0*hzn*2.0d0*pi
    omegac = 2.0e0*pi*omc*hzn
    omegaeq = 2.0e0*pi*omeq*hzn
  
    om2  = 0.0e0
    om2p = 0.0e0
    om4  = 0.0e0
    om4p = 0.0e0
  
    if (regime.eq.'r86') then
    ! eq (9) pag 91, rad86 an paper, appendix tables
    !-----------  if regime = r86  then we test Radler '86  an paper --------
      dw = xda2
      xw = xa2 
      xi = (x-xw)/dw
  
      if (x.le.xw-dw) then
        om0 = -1d0
        om0p = 0d0
      else if (x.gt.xw-dw.and.x.le.xw+dw) then
          om0 = -(1d0-3d0*xi/2d0+xi**3/2d0)/2.d0
          om0p =-(-3d0+3d0*xi**2 )/4.d0/dw
      else 
          om0 = 0d0
          om0p = 0d0
      endif
    else if(regime.eq.'ns')then
      ! pns model
      om0  = 1.d0+x**2
      om0p = 2d0*x
  
      ! for cylindrical rotation law
      om0  = 1d0+2d0*x**2/3d0
      om0p = 4d0*x/3d0
      om2  = -2d0*x**2/3d0 
      om2p = -4d0*x/3d0 
      om4  = 0d0
      om4p = 0d0
  
    else if (regime .eq. 'dj') then
      ! Dudley and James
      om0  = x*sin(x*pi)
      om0p = sin(pi*x) + cos(pi*x)*pi*x
  
    else if(regime.eq.'vega') then
      om0 = x
      om0p = 1.d0
  
    else if(regime.eq.'rp72'.or.regime.eq.'rm72'.or.regime.eq.'ekeri') then
      ! Roberts 72 alpha omega model 1 eq.4.3
      om0 = (x-x_in)
      om0p = 1.d0
  
    else if(regime.eq.'corona') then
      xx1=0.4d0
      xx2=0.9d0
      xxd1=0.1d0
      xxd2=0.02d0
      om0 = 1d0
      om0p= (-2d0*(1d0 + erf((x - xx1)/xxd1)))/(exp((x - xx2)**2/xxd2**2)*sqrt(pi)*xxd2) + & 
         & (2d0*(1d0 - erf((x - xx2)/xxd2)))/(exp((x - xx1)**2/xxd1**2)*sqrt(pi)*xxd1)
  
    else if (regime .eq. 'r72' .or. regime .eq. 'r72m' ) then
      ! Roberts 72 alpha omega model 2 eq.4.5  and model 5.6
      om0 = -19683.d0*(1-x**2)**5/40960.d0
      om0p = -19683.0d0*5d0*(1d0-x**2)**4*(-2d0*x)/40960.d0
  
    else if (regime .eq. 'r72b') then
      ! Roberts 72 alpha omega model eq.4.5
      om0 = -3.d0*sqrt(3.d0)*(1d0-x**2)**2/8.d0
      om0p = 3.d0*sqrt(3.d0)*(1d0-x**2)*4d0*x/8.d0
  
    else if(regime.eq.'br') then
      !  Roberts 72 alpha omega model 5.2  (bragisnki)
      om0 = -(3.d0*sqrt(3.d0)/8d0)*((1d0-x**2)**2d0)
      om0p = (3.d0*sqrt(3.d0)/8d0)*(+2d0*x*(1d0-x**2)*2d0)
  
    else if(regime.eq.'sk69'.or.regime.eq.'stm') then
      ! model 1 Steenbeck-Krause sk69 paper, alpha omega
      d1 = 0.075d0
      r1 = xc1
      om0 = (1d0-derf((x-r1)/d1))/2.d0
      om0p =-exp(-((x-r1)/d1)**2)/sqrt(pi)/d1
  
    else if(regime.eq.'st74') then
      !  Stix 74 paper, alpha omega
      d1 = 0.075d0
      r1 = 0.7d0
      om0 = (1d0-derf((x-r1)/d1))/2.d0
      om0p =-exp(-((x-r1)/d1)**2)/sqrt(pi)/d1
  
    else if (regime.eq.'h2') then
  
      c4 = 0d0
      omegaeq = 1.d0
      omegac  = oco
      delta   = (1.e0 + derf( (x-xc1)/dd1  ) )/ 2.e0
      c2 = c2_h2 
      ! c2= -coefficient of cos(theta) ^2 at the surface dr:   omega_surface=omega_eq-c2*cos(theta)^2
      om2  = -2d0*delta*c2/3.d0 
      om2p = (-2d0*c2/3.d0)*(exp(-((x-xc1)/dd1)**2)/sqrt(pi)/dd1) 
      om0 = delta*(-c2/3.0+1d0-omegac) + omegac   
      om0p = (-c2/3.0d0+1d0-omegac)*(exp(-((x-xc1)/dd1)**2)/sqrt(pi)/dd1)
  
    else if(regime.eq.'k2')then
      ! K giant model
      ! arcturus !!!!! 
      ! mimic the leonid omega
      om2=0.01d0/x-0.25d0
      om2p=-0.01d0/x/x
      om0  = 0.22d0/x+1.d0   
      om0p = -0.22d0/x/x
    
    else if(regime.eq.'h4')then
      ! omega is given by the following expression taken from helioseismology
      ! omega = omegac + delta*(omegaeq-omegac) - delta*c2*cos(theta)^2 -delta*c4*cos(theta)^4 
      !       = omega0 + omega2 * cos(theta)^2 + omega4 *cos(theta)^4
      !       = om0    + om2*cos(theta)^2  + om4*cos(theta)^4
      delta   = (1.e0 + derf( (x-xc1)/dd1  ) )/ 2.e0
      om0  =  (omegac + delta * ( omegaeq - omegac))/omegaeq
      om0p =  (exp(-((x-xc1)/dd1)**2)/sqrt(pi)/dd1)*(omegaeq - omegac)/omegaeq
      om2  = -delta*c2/omegaeq
      om2p = -(exp(-((x-xc1)/dd1)**2)/sqrt(pi)/dd1)*c2/omegaeq
      om4  = -delta*c4/omegaeq
      om4p =  -(exp(-((x-xc1)/dd1)**2)/sqrt(pi)/dd1)*c4/omegaeq
  
    else if(regime.eq.'h5')then !pure latitudinal dependence
      om0  =  2.78d0
      om0  =  1.d0
      om0p =  0d0
      om2  =  om0*(-0.13d0)
      om2p =  0d0
      om4  =  om0*(-0.16d0)
      om4p =  0d0
  
    else if(regime.eq.'h6')then  ! subsurface shear
      ome0 = 435.d0
      omeq = 452.5d0
      rt   = 0.69d0
      rc   = 0.71d0
      a2   = -61.0d0
      a4   = -73.5d0
      omet = 0.05d0
      omec = 0.05d0
      omes = 0.05d0
      rs   = 0.95d0
      b0   = 437.d0
      b4   = -1445.d0/2d0
      r    =  x 
  
      om0 = (b0 + ome0 + omeq - b0*r + & 
         &    (-((r - rc)*(7d0*a2 + 3d0*a4 + &
         &            35d0*(b0 - ome0 + omeq - b0*rs))* &
         &          derf((2d0*(r - rc))/omec)) + &
         &       (7d0*a2 + 3d0*a4 + 35d0*(b0 - ome0 + omeq - b0*rc))* &
         &        (r - rs)*derf((2d0*(r - rs))/omes) + & 
         &       (7d0*a2 + 3d0*a4)*(-rc + rs)*derf((2d0*(r - rt))/omet))/(35.d0*(rc - rs)))/2.d0
  
      om0p = (-b0 + ((4d0*(7d0*a2 + 3d0*a4 + &
        &            35d0*(b0 - ome0 + omeq - b0*rc))*(r - rs))/    &
        &        (dexp((4d0*(r - rs)**2)/omes**2)*omes*sqrt(pi)) +  &
        &       (4d0*(7d0*a2 + 3d0*a4)*(-rc + rs))/                     & 
        &        (dexp((4d0*(r - rt)**2)/omet**2)*omet*sqrt(pi)) -  &
        &       (4d0*(r - rc)*                                      &
        &          (7d0*a2 + 3d0*a4 + 35d0*(b0 - ome0 + omeq - b0*rs))  &
        &          )/                                             &
        &        (dexp((4d0*(r - rc)**2)/omec**2)*omec*sqrt(pi)) -  &
        &       (7d0*a2 + 3d0*a4 + 35d0*(b0 - ome0 + omeq - b0*rs))*    &
        &        derf((2*(r - rc))/omec) +                        &
        &       (7d0*a2 + 3d0*a4 + 35d0*(b0 - ome0 + omeq - b0*rc))*    &
        &        derf((2d0*(r - rs))/omes))/(35.d0*(rc - rs)))/2.d0     
  
      om2 = (a2*(1d0 + derf((2d0*(r - rt))/omet)))/2.d0
      om2p= (2d0*a2)/(dexp((4d0*(r - rt)**2)/omet**2)*omet*sqrt(pi))
  
      om4=(-(b4*(r - rc)*(1d0 - rs)*derf((2d0*(r - rc))/omec)) + & 
        &    b4*(1d0 - rc)*(r - rs)*derf((2d0*(r - rs))/omes) + & 
        &    (rc - rs)*(a4 + b4 - b4*r + a4*derf((2d0*(r - rt))/omet)))/(2.d0*(rc - rs))
  
      om4p=((-4d0*b4*(r - rc)*(1d0 - rs))/(dexp((4d0*(r - rc)**2)/omec**2)*omec*sqrt(pi)) + & 
         &    (4d0*b4*(1d0 - rc)*(r - rs))/(dexp((4d0*(r - rs)**2)/omes**2)*omes*sqrt(pi)) + &
         &    (-b4 + (4d0*a4)/(dexp((4d0*(r - rt)**2)/omet**2)*omet*sqrt(pi)))*(rc - rs) - &
         &    b4*(1d0 - rs)*derf((2d0*(r - rc))/omec) + b4*(1d0 - rc)*derf((2d0*(r - rs))/omes))/ &
         &  (2.d0*(rc - rs))
  
      om0   = om0  /omeq
      om0p  = om0p /omeq
      om2   = om2  /omeq
      om2p  = om2p /omeq
      om4   = om4  /omeq
      om4p  = om4p /omeq
  
    else if(regime.eq.'r2')then
      om0  = -3.d0*sqrt(3.d0)*(1d0-x**2)**2/8.d0
      om0p = -(3.d0*sqrt(3.d0)/8d0)*2d0*(1d0-x**2)*(-2d0*x)
      om2  = 0d0
      om2p = 0d0
      om4  = 0d0
      om4p = 0d0
    endif
    return
  end
  
  !-------------------------------------------------------------------------
  
  subroutine alt(x,a1,a2)
  
    !-------------------------------------------------------------------------
    !
    ! derf(x-x1) is one for x> x1 and -1 for x< x1   
    ! 
    ! by definition we have  derf(x) = 2/sqrt(pi) * int exp(-t^2) t=0..x
    !                        derf(x)' =2/sqrt(pi) * exp(-x^2) 
    !
    !-------------------------------------------------------------------------
  
    real(dp) :: x, a1, a2, d1, d_xi, delta, r1, x_xi, xi
  
    if (regime.eq.'ns') then
       a1 = (1d0+derf((x-xa1)/xda1))/2.e0
       a2 = exp(-(x-xa1)**2/xda1/xda1)/sqrt(pi)/xda1
  
    else if (regime.eq.'r86') then
       d_xi = xda1 
       x_xi = xa1 
       xi = (x-x_xi)/d_xi
  
       if (x .le. x_xi-d_xi) then
         a1=0d0
         a2=0d0
       else if (x.ge.x_xi-d_xi .and. x.le.x_xi+d_xi) then
         a1 = (-15d0+10d0*3d0*xi**2-3d0*5d0*xi**4)/16.d0/d_xi
         a2 = (10d0*3d0*2d0*xi-3d0*5d0*4d0*xi**3)/16.d0/d_xi/d_xi
       else if (x .gt. x_xi+d_xi) then
         a1=0d0
         a2=0d0
       endif
  
    else if (regime .eq. 'aq1') then
      ! simplest alpha2 model 
      a1=1d0
      a2=0d0
    else if(regime.eq.'aq2')then 
      ! sk' an 291, p271, 1969
      if(x.lt.xa1)then
        a1=0d0
        a2=0d0  
      else
        a1 = (x-xa1)*(x*3.d0-xa1-2.d0)*(27.d0/4.d0)
        a2 = (6.d0*x -4d0*xa1 -2.d0)*27.d0/((1d0-xa1)**3.0d0)/4.d0
      endif
  
    else if(regime.eq.'rp72')then
      ! roberts 72 alphaomega model  1 eq. 4.3
      a1= 1d0
      a2= 0d0
  
    else if(regime.eq.'rm72')then
      ! roberts 72 alphaomega model  1 eq. 4.3
      a1= -1d0
      a2= 0d0
  
    else if(regime.eq.'r72'.or.regime.eq.'r72m')then
      ! roberts 72 alphaomega model  2 eq. 4.5 and 5.6
      a1= 729d0*x**8*(1d0-x**2)**2/16.0d0
      a2= 729d0*(8d0*x**7*(1d0-x**2)**2+x**8*2d0*(1d0-x**2)*(-2d0*x))/16.0d0
  
    else if(regime.eq.'r72b')then
      ! roberts 72 alphaomega model 2 eq. 4.5
      a1 = 24.d0*sqrt(3.d0)*x**2*(1d0-x)**2
      a2 = 24.d0*sqrt(3.d0)*(2d0*x*(1d0-x)**2-x**2*2d0*(1d0-x))
  
    else if(regime.eq.'st74') then
      d1 = 0.075d0
      r1 = 0.7d0
      a1 = (1d0-derf((x-r1)/d1))/2.d0
      a2 =-exp(-((x-r1)/d1)**2)/sqrt(pi)/d1
  
    else if(regime.eq.'h2'.or.regime.eq.'h4'.or.regime.eq.'k2' &
       & .or.regime.eq.'h6'.or.regime.eq.'vega'.or.regime.eq.'ekeri')then
  
       a1 =(1d0-derf((x-xa2)/xda2) ) &
         &   *(1d0+derf((x-xa1)/xda1))/2.e0/2.e0
  
       a2 = (1d0-derf((x-xa2)/xda2)) &
         &         *exp(-(x-xa1)**2/xda1/xda1)/sqrt(pi)/xda1/2 &
         &        -(1d0+derf((x-xa1)/xda1)) &
         &         *exp(-(x-xa2)**2/xda2/xda2)/sqrt(pi)/xda2/2.
       ! divide by 3srt(3)/4 to normalize alpha  when sin^2 cos functional form is used
       if (c3 .ne. 0) then  
         a1 = a1 /0.3849d0
         a2 = a2 /0.3849d0
       endif         
  
    else if(regime.eq.'sk69'.or.regime.eq.'stm') then
      xda1=0.075d0
      a1 = (1d0+derf((x-xa1)/xda1))/2.e0
      a2 = exp(-(x-xa1)**2/xda1/xda1)/sqrt(pi)/xda1
  
    else if(regime.eq.'br')then
      a1 = sqrt(3.d0)*24d0*x**2*(1d0-x)**2
      a2 = sqrt(3.d0)*24d0*(2d0*x*(1d0-x)**2-x**2*(1d0-x)*2d0)
  
    else if(regime.eq.'dj')then
      a1 = (1d0+derf((x-xa1)/xda1))/2.e0
      a2 = exp(-(x-xa1)**2/xda1/xda1)/sqrt(pi)/xda1
  
    else if(regime.eq.'h5')then
      a1=1d0
      a2=0d0
  
    else if (regime .eq. 'corona') then
      a1 =(1d0-derf((x-xa2)/xda2) ) &
        &   *(1d0+derf((x-xa1)/xda1))/2.e0/2.e0
      a2 = (1d0-derf((x-xa2)/xda2)) &
        &         *exp(-(x-xa1)**2/xda1/xda1)/sqrt(pi)/xda1/2d0 &
        &        -(1d0+derf((x-xa1)/xda1)) &
        &         *exp(-(x-xa2)**2/xda2/xda2)/sqrt(pi)/xda2/2.d0
      ! charbonneau and macgregor 2001 o-b star model
    else if (regime.eq.'c01') then
      delta = 1.d0 - derf( (x-xa1)/xda1 )
      a1 = derf(2.d0*x/xda1)*delta/2.d0
      a2 = ((2.d0*exp(-(2.d0*x/xda1)**2))/(xda1*sqrt(pi)))*delta  &
        & - (exp(-((x-xa1)/xda1)**2 )/(sqrt(pi)*xda1))*derf(2.d0*x/xda1) 
    else 
      write(*,*) 'Configuration not available:', regime
      stop
    endif
    return
  end subroutine
  
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  
  subroutine eta_turb(x, e1, e2, e3)
  
    !---------------------------------------------------------------------------
    !   it gives eta_turb= e1, and derivatives eta' = e2, eta'' = e3 
    !   edr is the ratio eta_c/eta_top
    !---------------------------------------------------------------------------
  
    real(dp) :: x, x1, x2, xd
    real(dp) :: e1, e2, e3
    real(dp) :: delta, et2, ett
  
    x1 = xe1
    xd = xde1 
    ett=1.e0
    et2= 1.e0 !0.01
  
    e1 =   edr+(et2-edr)*(1d0+derf((x-x1)/xd))/2.d0
    e2 =  (et2-edr)*exp(-((x-x1)/xd)**2)/sqrt(pi)/xd 
    e3 =  -2d0*(et2-edr)*(x-x1)*exp(-((x-x1)/xd)**2)/sqrt(pi)/xd**3
  
    if (regime.eq.'ns') then
      e1 = et2+(edr-et2)*(1d0+derf((x-x1)/xd))/2.d0
      e2 =  (edr-et2)*exp(-((x-x1)/xd)**2)/sqrt(pi)/xd 
      e3 =  -2d0*(edr-et2)*(x-x1)*exp(-((x-x1)/xd)**2)/sqrt(pi)/xd**3
  
    !Charbonneau and MacGregor 2001 o-b star model
    else if (regime.eq.'c01') then
      et2 = 1.d0
      delta = 1.d0 - derf( (x-x1)/xd )
      e1 = edr + delta*(et2-edr)/2.d0
      e2 = -(exp( -((x-x1)/xd)**2 )/(sqrt(pi)*xd))*(et2-edr)
      e3 = (2.d0*(x-x1)/(sqrt(pi)*(xd)**3))*(exp( -((x-x1)/xd)**2) )*(et2-edr)
    endif
    return
  end subroutine
  
end module stellar_profiles 
