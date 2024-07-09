function plgndr(l,m,x)

  implicit none 
  integer :: l,m
  real :: plgndr,x
  integer :: i,ll
  real :: fact,pll,pmm,pmmp1,somx2

  if (m.lt.0 .or. m.gt.l .or. abs(x).gt.1.) write (*,*) 'bad arguments in plgndr'
  pmm=1.
  if(m.gt.0) then
    somx2=sqrt((1.-x)*(1.+x))
    fact=1.
    do i=1,m
      pmm=-pmm*fact*somx2
      fact=fact+2.
    enddo
  endif
  if(l.eq.m) then
    plgndr=pmm
  else
    pmmp1=x*(2*m+1)*pmm
    if(l.eq.m+1) then
      plgndr=pmmp1
    else
      do ll=m+2,l
        pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
        pmm=pmmp1
        pmmp1=pll
      enddo
      plgndr=pll
    endif
  endif
  return
 end
