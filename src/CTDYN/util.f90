module util

  use kind_parameters 
  use stdlib_sorting, only: sort_index

  implicit none

  private

  public :: hunt, sort2

contains

  integer function hunt (xx, n, x) result (jlo)
    !> Determine jlo such as xx(jl0) <= x < x(jl0+1).
  
    ! Arguments
    integer(i4) :: n
    real(dp) :: x, xx(n)

    ! Local variables
    integer(i4) :: inc, jhi, jm
    logical :: ascnd
  
    ascnd = xx(n).gt.xx(1)
    if (jlo .le. 0 .or. jlo .gt. n) then
      jlo=0
      jhi=n+1
    else 
      inc=1
      if (x.ge.xx(jlo).eqv.ascnd) then
        do while ((x.ge.xx(jhi).eqv.ascnd) .and. (jhi.lt.n+1))
          jhi = jlo+inc
          jlo = jhi
          inc = inc+inc
        enddo
        if (jhi.gt.n) then
          jhi = n+1
        endif
      else
        jhi=jlo
        do while ((x.lt.xx(jlo).eqv.ascnd) .and. (jlo.gt.0))
          jlo = jhi-inc
          jhi = jlo
          inc = inc+inc
        enddo
        if (jlo.lt.1) then
          jlo = 0
        endif
      endif
    endif
    do while (.not.jhi-jlo.eq.1)
      jm = (jhi+jlo)/2
      if (x.gt.xx(jm).eqv.ascnd) then
        jlo = jm
      else
        jhi = jm
      endif
    enddo 
    return 
  end function hunt

  subroutine sort2(n,arr,brr)
    ! Arguments
    integer(i4) :: n
    real(dp) :: arr(n), brr(n)

    ! Internal variables
    integer(i4) :: index(n)

    call sort_index (arr, index)
    brr = brr(index)
  end subroutine sort2
  
end module util
