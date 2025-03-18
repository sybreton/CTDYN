!  (c) Copr. 1986-92 Numerical Recipes Software <&)1.

module zbrent

  use kind_parameters
  use cio
  use dyna

  implicit none 

  private

  public :: zbr

contains 

  subroutine zbr(x1, x2, tol, root)
    
    ! Arguments
    real(dp) :: root, tol, x1, x2

    ! Local variables
    real(dp) :: eta, period, imag
    real(dp) :: x3
    integer :: iter
    real(dp) :: a, b, c, d, e, fa, fb, fc 
    real(dp) :: p, q, r, s, tol1
    integer, parameter :: itmax=100
    real(dp), parameter :: eps=3.e-8 
    
    a=x1
    b=x2
    it=1
    call dynamo(a, fa, imag, eta, period) 
    it=2
    call dynamo(b, fb, imag, eta, period)
    ! --------------------------------------------------------
    if ((fa.gt.0.d0 .and. fb.gt.0.d0) .or. (fa.lt.0.d0 .and. fb.lt.0.d0)) then
      write (*,*) 'Root must be bracketed for zbrent'
      stop
    endif
    ! --------------------------------------------------------
    c=b
    fc=fb
    do iter=1, itmax
      it=iter+2
      if ((fb.gt.0.d0 .and. fc.gt.0.d0).or.(fb.lt.0.d0 .and. fc.lt.0.d0)) then
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
      tol1=2.*eps*abs(b)+0.5d0*tol
      xm=.5d0*(c-b)
      if(abs(xm).le.tol1 .or. fb.eq.0.d0) then
        write(*,*) 'Tolerance threshold met at iteration', it
        root=b
        return
      endif
      if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
        s=fb/fa
        if(a.eq.c) then
          p=2.d0*xm*s
          q=1.d0-s
        else
          q=fa/fc
          r=fb/fc
          p=s*(2.d0*xm*q*(q-r)-(b-a)*(r-1.d0))
          q=(q-1.d0)*(r-1.d0)*(s-1.d0)
        endif
        if(p.gt.0.d0) q=-q
        p=abs(p)
        if(2.d0*p .lt. min(3.d0*xm*q-abs(tol1*q),abs(e*q))) then
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
      call dynamo(b, fb, imag, eta, period)
      if (abs(fb) .le. 1.0e-4) then
        root=b
        write(*,*) 'fb below 1e-4 at iteration number', iter
        return
       endif
    enddo
    
    ! Case where the loop terminated without finishing
    write(*,*) 'Brent algorithm exceeding maximum iterations'
    root=b
  
    return
  end

end module zbrent 
