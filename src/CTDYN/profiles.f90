subroutine rot(x, om0, om0p, om2, om2p, om4, om4p)

  use cio
  implicit none

  real :: x, x1, x2, x3, xd ,xd2, xd3, &
          & om0, om2, om0p, om2p, om4, om4p, a1, a2, xmax, &
          & omc, omeq, hzn, c2, c4, omegac, oms, xbi, omegaeq, &
          & a4, b0, b4, d1, delta, dw, facn, ome0, &
          & omec, omes, omet, r, r1, rc, rs, xi, xw, &
          & xx1, xx2, xxd1, xxd2 
       
  
  !--------- definition of omega see Dikpati & Chcarbonneau apj99.
               
  facn = 1 
  omc = 432.8e0*facn
  omeq = 460.7e0*facn
  hzn = 1.0

  c2 = +62.69e0*hzn*2.0*pi*facn
  c4 = +67.13e0*hzn*2.0*pi
  omegac = 2.0e0*pi*omc*hzn
  omegaeq = 2.0e0*pi*omeq*hzn

  om2  = 0.0e0
  om2p = 0.0e0
  om4  = 0.0e0
  om4p = 0.0e0

  if (ans1.eq.'r86') then
  ! eq (9) pag 91, rad86 an paper, appendix tables
  !-----------  if ans1 = r86  then we test radler '86  an paper --------
    dw = xda2
    xw = xa2 
    xi = (x-xw)/dw

    if (x.le.xw-dw) then
      om0 = -1
      om0p = 0
    else if (x.gt.xw-dw.and.x.le.xw+dw) then
        om0 = -(1-3*xi/2+xi**3/2)/2.
        om0p =-(-3+3*xi**2 )/4./dw
    else 
        om0 = 0
        om0p = 0
    endif
  else if(ans1.eq.'ns')then
    ! pns model
    om0  = 1.+x**2
    om0p = 2*x

    ! for cylindrical rotation law
    om0  = 1+2*x**2/3
    om0p = 4*x/3 
    om2  = -2*x**2/3 
    om2p = -4*x/3 
    om4  = 0.0e+00
    om4p = 0.0e+00

  else if (ans1 .eq. 'dj') then
    ! dudley and james
    om0  = x*sin(x*pi)
    om0p = sin(pi*x) + cos(pi*x)*pi*x

  else if(ans1.eq.'vega') then
    om0 = x
    om0p = 1.

  else if(ans1.eq.'rp72'.or.ans1.eq.'rm72'.or.ans1.eq.'ekeri') then
    ! roberts 72 alpha omega model 1 eq.4.3
    om0 = (x-x_in)
    om0p = 1.

  else if(ans1.eq.'corona') then
    xx1=0.4
    xx2=0.9
    xxd1=0.1
    xxd2=0.02
    om0 = 1
    om0p= (-2*(1 + erf((x - xx1)/xxd1)))/(exp((x - xx2)**2/xxd2**2)*sqrt(pi)*xxd2) + & 
       & (2*(1 - erf((x - xx2)/xxd2)))/(exp((x - xx1)**2/xxd1**2)*sqrt(pi)*xxd1)

  else if (ans1 .eq. 'r72' .or. ans1 .eq. 'r72m' ) then
    ! roberts 72 alpha omega model 2 eq.4.5  and model 5.6
    om0 = -19683.*(1-x**2)**5/40960.
    om0p = -19683.0*5*(1-x**2)**4*(-2*x)/40960.

  else if (ans1 .eq. 'r72b') then
    ! roberts 72 alpha omega model eq.4.5
    om0 = -3.*sqrt(3.)*(1-x**2)**2/8.
    om0p = 3.*sqrt(3.)*(1-x**2)*4*x/8.

  else if(ans1.eq.'br') then
    !  roberts 72 alpha omega model 5.2  (bragisnki)
    om0 = -(3.*sqrt(3.)/8)*((1-x**2)**2)
    om0p = (3.*sqrt(3.)/8)*(+2*x*(1-x**2)*2)

  else if(ans1.eq.'sk69'.or.ans1.eq.'stm') then
    ! model 1 steenbeck-krause sk69 paper, alpha omega
    d1 = 0.075
    r1 = rc1
    om0 = (1-derf((x-r1)/d1))/2.
    om0p =-exp(-((x-r1)/d1)**2)/sqrt(pi)/d1

  else if(ans1.eq.'st74') then
    !  stix 74 paper, alpha omega
    d1 = 0.075
    r1 = 0.7
    om0 = (1-derf((x-r1)/d1))/2.
    om0p =-exp(-((x-r1)/d1)**2)/sqrt(pi)/d1

  else if (ans1.eq.'h2') then

    c4 =0
    omegaeq = 1.
    omegac  = oco
    c2 =   rc2 
    delta   = (1.e0 + derf( (x-rc1)/dd1  ) )/ 2.e0
    ! c2= -coefficient of cos(theta) ^2 at the surface dr:   omega_surface=omega_eq-c2*cos(theta)^2
    om2  = -2*delta*c2/3.d0 
    om2p = (-2*c2/3.d0)*(exp(-((x-rc1)/dd1)**2)/sqrt(pi)/dd1) 
    om0 = delta*(-c2/3.0+1-omegac) +omegac   
    om0p = (-c2/3.0+1-omegac)*(exp(-((x-rc1)/dd1)**2)/sqrt(pi)/dd1)

  else if(ans1.eq.'k2')then
    ! k giant model
    ! arcturus !!!!! 
    ! mimic the leonid omega
    om2=0.01/x-0.25
    om2p=-0.01/x/x
    om0  = 0.22/x+1.   
    om0p = -0.22/x/x
  
  else if(ans1.eq.'h4')then
    ! omega is given by the following expression taken from helioseismology
    ! omega = omegac + delta*(omegaeq-omegac) - delta*c2*cos(theta)^2 -delta*c4*cos(theta)^4 
    !       = omega0 + omega2 * cos(theta)^2 + omega4 *cos(theta)^4
    !       = om0    + om2*cos(theta)^2  + om4*cos(theta)^4
    delta   = (1.e0 + derf( (x-rc1)/dd1  ) )/ 2.e0
    om0  =  (omegac + delta * ( omegaeq - omegac))/omegaeq
    om0p =  (exp(-((x-rc1)/dd1)**2)/sqrt(pi)/dd1)*(omegaeq - omegac)/omegaeq
    om2  = -delta*c2/omegaeq
    om2p = -(exp(-((x-rc1)/dd1)**2)/sqrt(pi)/dd1)*c2/omegaeq
    om4  = -delta*c4/omegaeq
    om4p =  -(exp(-((x-rc1)/dd1)**2)/sqrt(pi)/dd1)*c4/omegaeq

  else if(ans1.eq.'h5')then !pure latitudinal dependence
    om0  =  2.78
    om0  =  1.
    om0p =  0
    om2  =  om0*(-0.13)
    om2p =  0
    om4  =  om0*(-0.16)
    om4p =  0

  else if(ans1.eq.'h6')then  ! subsurface shear
    ome0 = 435.
    omeq = 452.5
    rt   = 0.69
    rc   = 0.71
    a2   = -61.0
    a4   = -73.5
    omet = 0.05
    omec = 0.05
    omes = 0.05
    rs   = 0.95
    b0   = 437.
    b4   = -1445./2
    r    =  x 

    om0 = (b0 + ome0 + omeq - b0*r + & 
       &    (-((r - rc)*(7*a2 + 3*a4 + &
       &            35*(b0 - ome0 + omeq - b0*rs))* &
       &          derf((2*(r - rc))/omec)) + &
       &       (7*a2 + 3*a4 + 35*(b0 - ome0 + omeq - b0*rc))* &
       &        (r - rs)*derf((2*(r - rs))/omes) + & 
       &       (7*a2 + 3*a4)*(-rc + rs)*derf((2*(r - rt))/omet))/(35.*(rc - rs)))/2.

    om0p = (-b0 + ((4*(7*a2 + 3*a4 + &
      &            35*(b0 - ome0 + omeq - b0*rc))*(r - rs))/    &
      &        (dexp((4*(r - rs)**2)/omes**2)*omes*sqrt(pi)) +  &
      &       (4*(7*a2 + 3*a4)*(-rc + rs))/                     & 
      &        (dexp((4*(r - rt)**2)/omet**2)*omet*sqrt(pi)) -  &
      &       (4*(r - rc)*                                      &
      &          (7*a2 + 3*a4 + 35*(b0 - ome0 + omeq - b0*rs))  &
      &          )/                                             &
      &        (dexp((4*(r - rc)**2)/omec**2)*omec*sqrt(pi)) -  &
      &       (7*a2 + 3*a4 + 35*(b0 - ome0 + omeq - b0*rs))*    &
      &        derf((2*(r - rc))/omec) +                        &
      &       (7*a2 + 3*a4 + 35*(b0 - ome0 + omeq - b0*rc))*    &
      &        derf((2*(r - rs))/omes))/(35.*(rc - rs)))/2.     

    om2 = (a2*(1 + derf((2*(r - rt))/omet)))/2.
    om2p= (2*a2)/(dexp((4*(r - rt)**2)/omet**2)*omet*sqrt(pi))

    om4=(-(b4*(r - rc)*(1 - rs)*derf((2*(r - rc))/omec)) + & 
      &    b4*(1 - rc)*(r - rs)*derf((2*(r - rs))/omes) + & 
      &    (rc - rs)*(a4 + b4 - b4*r + a4*derf((2*(r - rt))/omet)))/(2.*(rc - rs))

    om4p=((-4*b4*(r - rc)*(1 - rs))/(dexp((4*(r - rc)**2)/omec**2)*omec*sqrt(pi)) + & 
       &    (4*b4*(1 - rc)*(r - rs))/(dexp((4*(r - rs)**2)/omes**2)*omes*sqrt(pi)) + &
       &    (-b4 + (4*a4)/(dexp((4*(r - rt)**2)/omet**2)*omet*sqrt(pi)))*(rc - rs) - &
       &    b4*(1 - rs)*derf((2*(r - rc))/omec) + b4*(1 - rc)*derf((2*(r - rs))/omes))/ &
       &  (2.*(rc - rs))

    om0   = om0  /omeq
    om0p  = om0p /omeq
    om2   = om2  /omeq
    om2p  = om2p /omeq
    om4   = om4  /omeq
    om4p  = om4p /omeq

  else if(ans1.eq.'r2')then
    om0  = -3.*sqrt(3.)*(1-x**2)**2/8.
    om0p = -(3.*sqrt(3.)/8)*2*(1-x**2)*(-2*x)
    om2=0
    om2p=0
    om4 =0
    om4p =0
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

  use cio
  implicit none

  real :: x, a1, a2, d1, d_xi, delta, r1, x_xi, xi

  if (ans1.eq.'ns') then
     a1 = (1+derf((x-xa1)/xda1))/2.e0
     a2 = exp(-(x-xa1)**2/xda1/xda1)/sqrt(pi)/xda1

  else if (ans1.eq.'r86') then
     d_xi = xda1 
     x_xi = xa1 
     xi = (x-x_xi)/d_xi

     if (x .le. x_xi-d_xi) then
       a1=0
       a2=0
     else if (x.ge.x_xi-d_xi .and. x.le.x_xi+d_xi) then
       a1 = (-15+10*3*xi**2-3*5*xi**4)/16./d_xi
       a2 = (10*3*2*xi-3*5*4*xi**3)/16./d_xi/d_xi
     else if (x .gt. x_xi+d_xi) then
       a1=0
       a2=0
     endif

  else if (ans1 .eq. 'aq1') then
    ! simplest alpha2 model 
    a1=1
    a2=0
  else if(ans1.eq.'aq2')then 
    ! sk' an 291, p271, 1969
    if(x.lt.xa1)then
      a1=0
      a2=0  
    else
      a1 = (x-xa1)*(x*3.-xa1-2.)*(27./4.)
      a2 = (6.*x -4*xa1 -2.)*27./((1-xa1)**3.0)/4.
    endif

  else if(ans1.eq.'rp72')then
    ! roberts 72 alphaomega model  1 eq. 4.3
    a1= 1
    a2= 0

  else if(ans1.eq.'rm72')then
    ! roberts 72 alphaomega model  1 eq. 4.3
    a1= -1
    a2= 0

  else if(ans1.eq.'r72'.or.ans1.eq.'r72m')then
    ! roberts 72 alphaomega model  2 eq. 4.5 and 5.6
    a1= 729*x**8*(1-x**2)**2/16.0
    a2= 729*(8*x**7*(1-x**2)**2+x**8*2*(1-x**2)*(-2*x))/16.0

  else if(ans1.eq.'r72b')then
    ! roberts 72 alphaomega model 2 eq. 4.5
    a1 = 24.*sqrt(3.)*x**2*(1-x)**2
    a2 = 24.*sqrt(3.)*(2*x*(1-x)**2-x**2*2*(1-x))

  else if(ans1.eq.'st74') then
    d1 = 0.075
    r1 = 0.7
    a1 = (1-derf((x-r1)/d1))/2.
    a2 =-exp(-((x-r1)/d1)**2)/sqrt(pi)/d1

  else if(ans1.eq.'h2'.or.ans1.eq.'h4'.or.ans1.eq.'k2' &
     & .or.ans1.eq.'h6'.or.ans1.eq.'vega'.or.ans1.eq.'ekeri')then

     a1 =(1-derf((x-xa2)/xda2) ) &
       &   *(1+derf((x-xa1)/xda1))/2.e0/2.e0

     a2 = (1-derf((x-xa2)/xda2)) &
       &         *exp(-(x-xa1)**2/xda1/xda1)/sqrt(pi)/xda1/2 &
       &        -(1+derf((x-xa1)/xda1)) &
       &         *exp(-(x-xa2)**2/xda2/xda2)/sqrt(pi)/xda2/2.
     ! divide by 3srt(3)/4 to normalize alpha  when sin^2 cos functional form is used
     if (c3 .ne. 0) then  
       a1 = a1 /0.3849
       a2 = a2 /0.3849
     endif         

  else if(ans1.eq.'sk69'.or.ans1.eq.'stm') then
    xda1=0.075
    a1 = (1+derf((x-xa1)/xda1))/2.e0
    a2 = exp(-(x-xa1)**2/xda1/xda1)/sqrt(pi)/xda1

  else if(ans1.eq.'br')then
    a1 = sqrt(3.)*24*x**2*(1-x)**2
    a2 = sqrt(3.)*24*(2*x*(1-x)**2-x**2*(1-x)*2)

  else if(ans1.eq.'dj')then
    a1 = (1+derf((x-xa1)/xda1))/2.e0
    a2 = exp(-(x-xa1)**2/xda1/xda1)/sqrt(pi)/xda1

  else if(ans1.eq.'h5')then
    a1=1
    a2=0

  else if (ans1 .eq. 'corona') then
    a1 =(1-derf((x-xa2)/xda2) ) &
      &   *(1+derf((x-xa1)/xda1))/2.e0/2.e0
    a2 = (1-derf((x-xa2)/xda2)) &
      &         *exp(-(x-xa1)**2/xda1/xda1)/sqrt(pi)/xda1/2 &
      &        -(1+derf((x-xa1)/xda1)) &
      &         *exp(-(x-xa2)**2/xda2/xda2)/sqrt(pi)/xda2/2.
    ! charbonneau and macgregor 2001 o-b star model
  else if (ans1.eq.'c01') then
    delta = 1. - derf( (x-xa1)/xda1 )
    a1 = derf(2.*x/xda1)*delta/2.
    a2 = ((2.*exp(-(2.*x/xda1)**2))/(xda1*sqrt(pi)))*delta  - (exp(-((x-xa1)/xda1)**2 )/(sqrt(pi)*xda1))*derf(2.*x/xda1) 
  else 
    write(*,*) 'Configuration not available:', ans1
    stop
  endif
  return
end

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine eta_turb(x, e1, e2, e3)

  !---------------------------------------------------------------------------
  !   it gives eta_turb= e1, and derivatives eta' = e2, eta'' = e3 
  !   edr is the ratio eta_c/eta_top
  !---------------------------------------------------------------------------

  use cio
  implicit none

  real :: x, x1, x2, xd
  real :: e1, e2, e3
  real :: delta, et2, ett

  x1 = xe1
  xd = xde1 
  ett=1.e0
  et2= 1.e0 !0.01

  e1 =   edr+(et2-edr)*(1+derf((x-x1)/xd))/2.
  e2 =  (et2-edr)*exp(-((x-x1)/xd)**2)/sqrt(pi)/xd 
  e3 =  -2*(et2-edr)*(x-x1)*exp(-((x-x1)/xd)**2)/sqrt(pi)/xd**3

  if (ans1.eq.'ns') then
    e1 = et2+(edr-et2)*(1+derf((x-x1)/xd))/2.
    e2 =  (edr-et2)*exp(-((x-x1)/xd)**2)/sqrt(pi)/xd 
    e3 =  -2*(edr-et2)*(x-x1)*exp(-((x-x1)/xd)**2)/sqrt(pi)/xd**3

  !charbonneau and macgregor 2001 o-b star model
  else if (ans1.eq.'c01') then
    et2 = 1.
    delta = 1. - derf( (x-x1)/xd )
    e1 = edr + delta*(et2-edr)/2.
    e2 = -(exp( -((x-x1)/xd)**2 )/(sqrt(pi)*xd))*(et2-edr)
    e3 =(2.*(x-x1)/(sqrt(pi)*(xd)**3))*(exp( -((x-x1)/xd)**2) )*(et2-edr)
  endif
  return
end

