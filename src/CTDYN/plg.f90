module plg

  use kind_parameters 

  implicit none

  private

  public :: plgndr

contains

  real (dp) function plgndr(l,m,x)
  
    integer :: l, m
    real(dp) :: x

    integer :: i, ll
    real(dp) :: fact, pll, pmm, pmmp1, somx2
  
    if (m.lt.0 .or. m.gt.l .or. abs(x).gt.1.d0) write (*,*) 'bad arguments in plgndr'
    pmm=1.d0
    if(m.gt.0) then
      somx2=sqrt((1.d0-x)*(1.d0+x))
      fact=1.d0
      do i=1,m
        pmm=-pmm*fact*somx2
        fact=fact+2.d0
      enddo
    endif
    if(l.eq.m) then
      plgndr=pmm
    else
      pmmp1=x*(2d0*m+1d0)*pmm
      if(l.eq.m+1) then
        plgndr=pmmp1
      else
        do ll=m+2,l
          pll=(x*(2d0*ll-1d0)*pmmp1-(ll+m-1)*pmm)/(ll-m)
          pmm=pmmp1
          pmmp1=pll
        enddo
        plgndr=pll
      endif
    endif
    return
  end function

end module plg
