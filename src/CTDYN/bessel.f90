module bessel 

  use kind_parameters 
  use cio

  implicit none

contains
  
  real(dp) function chebev(a,b,c,m,x)
    integer :: m
    real(dp) :: a,b,x,c(m)

    integer :: j
    real(dp) :: d,dd,sv,y,y2
  
    if ((x-a)*(x-b) .gt. 0.) write (*,*) 'x not in range in chebev'
    d=0.
    dd=0.
    y=(2.*x-a-b)/(b-a)
    y2=2.*y
    do j=m,2,-1
      sv=d
      d=y2*d-dd+c(j)
      dd=sv
    enddo
    chebev=y*d-dd+0.5*c(1)
    return
  end

  subroutine bess(x,j,corr)
  
    ! enter beta (=x) and j,  and output the correction at r/R=1
  
    real(dp) :: x, corr, cord, corn, rj, &
         & rj1, rjm1, rjp, rjp1, rjpm1, rjpz, &
         & rjpz1, rjz, rjz1, ry, ry1, rym1, &
         & ryp, ryp1, rypm, rypm1, rypz, rypz1, ryz, &
         & ryz1, xmu, xnu, xnu1, z
    integer :: j
  
    ! gamma is fixed so that at about z=zeta_r  so that the field lines open!
  
    z=zeta_r
  
    if (x.eq.0) then 
      corr= 0
      gam = 0
  
    else 
      if (x.gt.0) then
        xnu=1.0*j+1/2.
        xnu1=1.0*j+3./2.
        xmu = 1.0*j-1./2.
        ! bessel function at r=R, argument is beta!
        call bessjy(x,xnu,rj,ry,rjp,ryp)
        call bessjy(x,xnu1,rj1,ry1,rjp1,ryp1)
        ! bessel function at r=z*R, argument is beta*r
        call bessjy(x*z,xmu,rjm1,rym1,rjpm1,rypm1)
        call bessjy(x*z,xnu,rjz,ryz,rjpz,rypz)
        call bessjy(x*z,xnu1,rjz1,ryz1,rjpz1,rypz1)
        ! derivative of psi and phi are zero at the normalized radius z=r/R
        gam = (-z*x*rjm1+rjz+z*x*rjz1)/(z*x*rym1-ryz-z*x*ryz1) 
        corn = (ry1*(-z*x*rjm1+rjz+z*x*rjz1)+ rj1*(z*x*rym1-ryz-z*x*ryz1))
        cord = (ry*(-z*x*rjm1+rjz+z*x*rjz1)+ rj*(z*x*rym1-ryz-z*x*ryz1))
        corr = x*corn/cord-(j*2.d0 +1)
        
      else 
        x=-x
        xnu=1.0*j+1/2.
        xnu1=1.0*j+3./2.
        xmu = 1.0*j-1./2.
        ! bessel function at r=r, argument is beta!
        call bessjy(x,xnu,rj,ry,rjp,ryp)
        call bessjy(x,xnu1,rj1,ry1,rjp1,ryp1)
        ! bessel function at r=z*r, argument is beta*r
        call bessjy(x*z,xmu,rjm1,rym1,rjpm1,rypm1)
        call bessjy(x*z,xnu,rjz,ryz,rjpz,rypz)
        call bessjy(x*z,xnu1,rjz1,ryz1,rjpz1,rypz1)
        ! derivative of psi and phi are zero at the normalized radius z=r/r ??? ???
        gam = -(-z*x*rjm1+rjz+z*x*rjz1)/(z*x*rym1-ryz-z*x*ryz1) 
        corn = -(ry1*(-z*x*rjm1+rjz+z*x*rjz1)+ rj1*(z*x*rym1-ryz-z*x*ryz1))
        cord = -(ry*(-z*x*rjm1+rjz+z*x*rjz1)+ rj*(z*x*rym1-ryz-z*x*ryz1))
        corr = x*corn/cord-(j*2.d0 +1)
      endif
    endif
    return
  end
  
  subroutine bessjy(x,xnu,rj,ry,rjp,ryp)
  
    real(dp) :: rj,rjp,ry,ryp,x,xnu
    integer, parameter :: maxit=10000
    real(dp), parameter :: eps=1.e-10, fpmin=1.e-30, xmin=2.
    integer :: i,isign,l,nl
    real(dp) :: a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,dr,e, &
       &f,fact,fact2,fact3,ff,gam,gam1,gam2,gammi,gampl,h,p,pimu,pimu2,q, &
       &r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,rymu,rymup,rytemp,sum,sum1, &
       &temp,w,x2,xi,xi2,xmu,xmu2
  
    if (x.le.0. .or. xnu.lt.0.) then
      write (*,*) 'bad arguments in bessjy'
    endif 
    if (x.lt.xmin) then
      nl=int(xnu+.5d0)
    else
      nl=max(0,int(xnu-x+1.5d0))
    endif
    xmu=xnu-nl
    xmu2=xmu*xmu
    xi=1.d0/x
    xi2=2.d0*xi
    w=xi2/pi
    isign=1
    h=xnu*xi
    if(h.lt.fpmin)h=fpmin
    b=xi2*xnu
    d=0.d0
    c=h
    i = 1 
    do while (i .le. maxit .and. abs(del-1.d0) .gt. eps)
      b=b+xi2
      d=b-d
      if (abs(d) .lt. fpmin) d=fpmin
      c=b-1.d0/c
      if (abs(c) .lt. fpmin) c=fpmin
      d=1.d0/d
      del=c*d
      h=del*h
      if (d .lt. 0.d0) isign=-isign
      i = i+1
    enddo
    if (abs(del-1.d0) .gt.eps) then
        write (*,*) 'x too large in bessjy; try asymptotic expansion'
    endif 
    rjl=isign*fpmin
    rjpl=h*rjl
    rjl1=rjl
    rjp1=rjpl
    fact=xnu*xi
    do l=nl,1,-1
      rjtemp=fact*rjl+rjpl
      fact=fact-xi
      rjpl=fact*rjtemp-rjl
      rjl=rjtemp
    enddo
    if (rjl .eq. 0.d0) rjl=eps
    f=rjpl/rjl
    if (x.lt.xmin) then
      x2=.5d0*x
      pimu=pi*xmu
      if (abs(pimu) .lt. eps) then
        fact=1.d0
      else
        fact=pimu/sin(pimu)
      endif
      d=-log(x2)
      e=xmu*d
      if (abs(e) .lt. eps )then
        fact2=1.d0
      else
        fact2=sinh(e)/e
      endif
      call beschb(xmu,gam1,gam2,gampl,gammi)
      ff=2.d0/pi*fact*(gam1*cosh(e)+gam2*fact2*d)
      e=exp(e)
      p=e/(gampl*pi)
      q=1.d0/(e*pi*gammi)
      pimu2=0.5d0*pimu
      if (abs(pimu2) .lt. eps) then
        fact3=1.d0
      else
        fact3=sin(pimu2)/pimu2
      endif
      r=pi*pimu2*fact3*fact3
      c=1.d0
      d=-x2*x2
      sum=ff+r*q
      sum1=p
      i = 1
      do while (i .le. maxit .and. abs(del) .gt. (1.d0+abs(sum))*eps)
        ff=(i*ff+p+q)/(i*i-xmu2)
        c=c*d/i
        p=p/(i-xmu)
        q=q/(i+xmu)
        del=c*(ff+r*q)
        sum=sum+del
        del1=c*p-i*del
        sum1=sum1+del1
        i = i + 1
      enddo
      if(abs(del).gt.(1.d0+abs(sum))*eps) then
          write (*,*) 'bessy series failed to converge'
      endif
      rymu=-sum
      ry1=-sum1*xi2
      rymup=xmu*xi*rymu-ry1
      rjmu=w/(rymup-f*rymu)
    else
      a=.25d0-xmu2
      p=-.5d0*xi
      q=1.d0
      br=2.d0*x
      bi=2.d0
      fact=a*xi/(p*p+q*q)
      cr=br+q*fact
      ci=bi+p*fact
      den=br*br+bi*bi
      dr=br/den
      di=-bi/den
      dlr=cr*dr-ci*di
      dli=cr*di+ci*dr
      temp=p*dlr-q*dli
      q=p*dli+q*dlr
      p=temp
      i=2
      do while (i .le. maxit .and. abs(dlr-1.d0)+abs(dli) .gt. eps)
        a=a+2*(i-1)
        bi=bi+2.d0
        dr=a*dr+br
        di=a*di+bi
        if(abs(dr)+abs(di).lt.fpmin)dr=fpmin
        fact=a/(cr*cr+ci*ci)
        cr=br+cr*fact
        ci=bi-ci*fact
        if(abs(cr)+abs(ci).lt.fpmin)cr=fpmin
        den=dr*dr+di*di
        dr=dr/den
        di=-di/den
        dlr=cr*dr-ci*di
        dli=cr*di+ci*dr
        temp=p*dlr-q*dli
        q=p*dli+q*dlr
        p=temp
        i = i+1
      enddo
      if (abs(dlr-1.d0)+abs(dli) .gt. eps) then
        write (*,*) 'cf2 failed in bessjy'
      endif
      gam=(p-f)/q
      rjmu=sqrt(w/((p-f)*gam+q))
      rjmu=sign(rjmu,rjl)
      rymu=rjmu*gam
      rymup=rymu*(p+q/gam)
      ry1=xmu*xi*rymu-rymup
    endif
    fact=rjmu/rjl
    rj=rjl1*fact
    rjp=rjp1*fact
    do i=1,nl
      rytemp=(xmu+i)*xi2*ry1-rymu
      rymu=ry1
      ry1=rytemp
    enddo
    ry=rymu
    ryp=xnu*xi*rymu-ry1
    return
  end
  
  subroutine beschb(x,gam1,gam2,gampl,gammi)
    real(dp) :: gam1,gam2,gammi,gampl,x
    integer, parameter :: nuse1=5, nuse2=5
    real(dp) :: xx,c1(7),c2(8) 
    save c1,c2
    data c1/-1.142022680371168d0,6.5165112670737d-3,3.087090173086d-4, &
         &-3.4706269649d-6,6.9437664d-9,3.67795d-11,-1.356d-13/ 
    data c2/1.843740587300905d0,-7.68528408447867d-2, &
         &1.2719271366546d-3,-4.9717367042d-6,-3.31261198d-8,2.423096d-10, &
         &-1.702d-13,-1.49d-15/
  
    xx=8.d0*x*x-1.d0
    gam1=chebev(-1.d0,1.d0,c1,nuse1,xx)
    gam2=chebev(-1.d0,1.d0,c2,nuse2,xx)
    gampl=gam2-x*gam1
    gammi=gam2+x*gam1
    return
  end

end module bessel
