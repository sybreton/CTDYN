module b3func

  use kind_parameters

  implicit none

contains

  real(dp) function b3m4 (n)
    integer :: n
    real(dp) :: cgret
    cgret = -n ** 2 * (n + 0.1E1) / (0.2E1 * n + 0.1E1) / (0.2E1 * n & 
    &- 0.1E1) * (n - 0.2E1) * (n - 0.1E1) / (0.2E1 * n - 0.5E1) / (0.2E1 &
    & * n - 0.3E1)
    b3m4 = cgret
    return
  end
  
  real(dp) function b3p4 (n)
    integer :: n
    real(dp) :: cgret
    cgret = n * (n + 0.1E1) ** 2 / (0.2E1 * n + 0.1E1) / (0.2E1 * n + &
    & 0.3E1) * (n + 0.2E1) * (n + 0.3E1) / (0.2E1 * n + 0.5E1) / (0.2E1 & 
    &* n + 0.7E1)
    b3p4 = cgret
    return
  end
  
  real(dp) function b3m2 (n)
    integer :: n
    real(dp) :: cgret
    cgret = -(n + 0.1E1) * n ** 2 * (0.2E1 * n ** 2 - 0.3E1 * n - 0.8E1) &
    &/ (0.2E1 * n + 0.3E1) / (0.4E1 * n ** 2 - 0.1E1) / (0.2E1 * n - 0.5E1)
    b3m2 = cgret
    return
  end
  
  real(dp) function b3p2 (n)
    integer :: n
    real(dp) :: cgret
    cgret = (n + 0.1E1) ** 2 * n * (0.2E1 * n ** 2 + 0.7E1 * n - 0.3E1) &
    &/ (0.2E1 * n + 0.7E1) / (0.2E1 * n + 0.3E1) / (0.4E1 * n ** 2 - 0.1E1)
    b3p2 = cgret
    return
  end
  
  real(dp) function b30 (n)
    integer :: n
    real(dp) :: cgret
    cgret = 0.3E1 * (n + 0.1E1) * n * (n ** 2 + n - 0.3E1) / (0.2E1 * &
    & n - 0.1E1) / (0.2E1 * n + 0.5E1) / (0.4E1 * n ** 2 - 0.9E1)
    b30 = cgret
    return
  end

end module b3func
