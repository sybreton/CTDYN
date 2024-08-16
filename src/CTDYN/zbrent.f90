!  (c) Copr. 1986-92 Numerical Recipes Software <&)1.

subroutine zbr(x1,x2,tol,zbrent)
  
  use cio

  implicit none 

  integer :: itmax
  real :: zbrent,tol,x1,x2,func,eps,x3
  parameter (itmax=100,eps=3.e-8)
  integer :: iter
  real :: a, b, c, d, e, fa, fb, fc 
  real :: p, q, r, s, tol1
  
  a=x1
  b=x2
  it=1
  call dynamo(a,fa) 
  it=2
  call dynamo(b,fb)
  ! --------------------------------------------------------
  if ((fa.gt.0. .and. fb.gt.0.) .or. (fa.lt.0. .and. fb.lt.0.)) then
    write (*,*) 'Root must be bracketed for zbrent'
    stop
  endif
  ! --------------------------------------------------------
  c=b
  fc=fb
  do iter=1, itmax
    it=iter+2
    if ((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.)) then
      c=a
      fc=fa
      d=b-a
      e=d
    endif
    if (abs(fc).lt.abs(fb)) then
      a=b
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
    endif
    tol1=2.*eps*abs(b)+0.5*tol
    xm=.5*(c-b)
    if(abs(xm).le.tol1 .or. fb.eq.0.)then
      write(*,*) 'tolerance threshold met at iteration', it
      zbrent=b
      return
    endif
    if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
      s=fb/fa
      if(a.eq.c) then
        p=2.*xm*s
        q=1.-s
      else
        q=fa/fc
        r=fb/fc
        p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
        q=(q-1.)*(r-1.)*(s-1.)
      endif
      if(p.gt.0.) q=-q
      p=abs(p)
      if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
        e=d
        d=p/q
      else
        d=xm
        e=d
      endif
    else
      d=xm
      e=d
    endif
    a=b
    fa=fb
    if(abs(d) .gt. tol1) then
      b=b+d
    else
      b=b+sign(tol1,xm)
    endif
    call dynamo(b,fb)
    if (abs(fb) .le. 1.0e-4) then
      zbrent=b
      write(*,*) 'fb below 1e-4 at iteration number', iter
      write(*,*)
      return
     endif
  enddo
  
  ! Case where the loop terminated without finishing
  write(*,*) 'zbrent exceeding maximum iterations'
  zbrent=b

  return
end
