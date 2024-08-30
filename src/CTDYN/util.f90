module util

  use kind_parameters 

  implicit none

  private

  public :: hunt, sort2

contains

  subroutine hunt(xx, n, x, jlo)
  
    ! Arguments
    integer :: jlo, n
    real(dp) :: x, xx(n)

    ! Local variables
    integer :: inc, jhi, jm
    logical :: ascnd
  
    ascnd = xx(n).gt.xx(1)
    if (jlo .le. 0 .or. jlo .gt. n) then
      jlo=0
      jhi=n+1
      goto 3
    endif
    inc=1
    if (x.ge.xx(jlo).eqv.ascnd) then
  1   jhi=jlo+inc
      if(jhi.gt.n)then
        jhi=n+1
      else if (x.ge.xx(jhi).eqv.ascnd) then
        jlo=jhi
        inc=inc+inc
        goto 1
      endif
    else
      jhi=jlo
  2   jlo=jhi-inc
      if (jlo.lt.1) then
        jlo=0
      else if (x.lt.xx(jlo).eqv.ascnd) then
        jhi=jlo
        inc=inc+inc
        goto 2
      endif
    endif
  3 if(jhi-jlo.eq.1) then 
      return
    endif 
    jm=(jhi+jlo)/2
    if (x.gt.xx(jm).eqv.ascnd) then
      jlo=jm
    else
      jhi=jm
    endif
    goto 3
  end subroutine

  subroutine sort2(n,arr,brr)

    integer :: n
    real(dp) :: arr(n), brr(n)
  
    integer, parameter :: m=7
    integer, parameter :: nstack=500
    integer :: i, ir, j, jstack, k, l, istack(nstack)
    real(dp) :: a, b, temp
  
    jstack=0
    l=1
    ir=n
  1 if (ir-l.lt.m) then
      do j=l+1,ir
        a=arr(j)
        b=brr(j)
        do i=j-1,l,-1
          if (arr(i).le.a) goto 2
          arr(i+1)=arr(i)
          brr(i+1)=brr(i)
        enddo
        i=l-1
  2     arr(i+1)=a
        brr(i+1)=b
      enddo
      if(jstack.eq.0)return
      ir=istack(jstack)
      l=istack(jstack-1)
      jstack=jstack-2
    else
      k=(l+ir)/2
      temp=arr(k)
      arr(k)=arr(l+1)
      arr(l+1)=temp
      temp=brr(k)
      brr(k)=brr(l+1)
      brr(l+1)=temp
      if(arr(l).gt.arr(ir))then
        temp=arr(l)
        arr(l)=arr(ir)
        arr(ir)=temp
        temp=brr(l)
        brr(l)=brr(ir)
        brr(ir)=temp
      endif
      if (arr(l+1).gt.arr(ir)) then
        temp=arr(l+1)
        arr(l+1)=arr(ir)
        arr(ir)=temp
        temp=brr(l+1)
        brr(l+1)=brr(ir)
        brr(ir)=temp
      endif
      if(arr(l).gt.arr(l+1))then
        temp=arr(l)
        arr(l)=arr(l+1)
        arr(l+1)=temp
        temp=brr(l)
        brr(l)=brr(l+1)
        brr(l+1)=temp
      endif
      i=l+1
      j=ir
      a=arr(l+1)
      b=brr(l+1)
  3   i=i+1
      if(arr(i).lt.a) goto 3
  4   j=j-1
      if (arr(j).gt.a) goto 4
      if (j.lt.i) goto 5
      temp=arr(i)
      arr(i)=arr(j)
      arr(j)=temp
      temp=brr(i)
      brr(i)=brr(j)
      brr(j)=temp
      goto 3
  5   arr(l+1)=arr(j)
      arr(j)=a
      brr(l+1)=brr(j)
      brr(j)=b
      jstack=jstack+2
      if (jstack .gt. nstack) write (*,*) 'nstack too small in sort2'
      if(ir-i+1.ge.j-l)then
        istack(jstack)=ir
        istack(jstack-1)=i
        ir=j-1
      else
        istack(jstack)=j-1
        istack(jstack-1)=l
        l=i
      endif
    endif
    goto 1
  end subroutine
  
  function factrl(n)

    integer :: n
    real(dp) :: factrl
    integer :: j, ntop
    real(dp) :: a(33),gammln
    save ntop,a
    data ntop,a(1)/0,1./
  
    if (n.lt.0) then
      write (*,*) 'Negative factorial in factrl'
    else if (n.le.ntop) then
      factrl=a(n+1)
    else if (n.le.32) then
      do j=ntop+1,n
        a(j+1)=j*a(j)
      enddo
      ntop=n
      factrl=a(n+1)
    else
      factrl=exp(gammln(n+1.))
    endif
    return
  end function
  
  function gammln(xx)

    real(dp) :: gammln, xx
    integer :: j
    real(dp) :: ser,stp,tmp,x,y,cof(6)
    save cof,stp
    data cof,stp/76.18009172947146d0,-86.50532032941677d0, &
         & 24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
         & -.5395239384953d-5,2.5066282746310005d0/
  
    x=xx
    y=x
    tmp=x+5.5d0
    tmp=(x+0.5d0)*log(tmp)-tmp
    ser=1.000000000190015d0
    do j=1,6
      y=y+1.d0
      ser=ser+cof(j)/y
    enddo
    gammln=tmp+log(stp*ser/x)
    return
  end function

end module util
