module b4func

  use kind_parameters 

  implicit none

contains

  real(dp) function b4p5 (n)
    integer :: n
    real(dp) :: cgret
    cgret = (n + 0.4E1) * n * (n + 0.1E1) ** 2 / (0.2E1 * n + 0.1E1) & 
      &/ (0.2E1 * n + 0.3E1) * (n + 0.2E1) * (n + 0.3E1) / (0.2E1 * n + 0.5E1) &
      & / (0.2E1 * n + 0.7E1) / (0.2E1 * n + 0.9E1)
    b4p5 = cgret
    return
  end
  
  real(dp) function b4p3 (n)
    integer :: n
    real(dp) :: cgret
    cgret = (n + 0.2E1) * n * (n + 0.1E1) ** 2 * (0.3E1 * n ** 2 + 0.13E2 &
      & * n - 0.6E1) / (0.2E1 * n + 0.9E1) / (0.2E1 * n + 0.5E1) / (0.2E1 &
      & * n + 0.3E1) / (0.4E1 * n ** 2 - 0.1E1)
    b4p3 = cgret
    return
  end
  
  real(dp) function b4p1 (n)
    integer :: n
    real(dp) :: cgret
    cgret = (0.2E1 * n ** 4 + 0.12E2 * n ** 3 + 0.4E1 * n ** 2 - 0.42E2 &
      & * n + 0.9E1) * (n + 0.1E1) * n / (0.2E1 * n - 0.3E1) / (0.2E1 * & 
      &n + 0.5E1) / (0.2E1 * n + 0.7E1) / (0.4E1 * n ** 2 - 0.1E1)
    b4p1 = cgret
    return
  end
  
  real(dp) function b4m1 (n)
    integer :: n
    real(dp) :: cgret
    cgret = -n * (n + 0.1E1) * (0.2E1 * n ** 4 - 0.4E1 * n ** 3 - 0.20E2 &
      & * n ** 2 + 0.22E2 * n + 0.45E2) / (0.2E1 * n + 0.1E1) / (0.4E1 & 
      &* n ** 2 - 0.9E1) / (0.4E1 * n ** 2 - 0.25E2)
    b4m1 = cgret
    return
  end
  
  real(dp) function b4m3 (n)
    integer :: n
    real(dp) :: cgret
    cgret = -n ** 2 * (0.3E1 * n ** 2 - 0.7E1 * n - 0.16E2) / (0.2E1 &
      &* n - 0.7E1) * (n ** 2 - 0.1E1) / (0.4E1 * n ** 2 - 0.1E1) / (0.4E1 &
      & * n ** 2 - 0.9E1)
    b4m3 = cgret
    return
  end
  
  real(dp) function b4m5 (n)
    integer :: n
    real(dp) :: cgret
    cgret = -n ** 2 * (n - 0.2E1) / (0.2E1 * n - 0.5E1) / (0.2E1 * n &
      &- 0.3E1) * (n - 0.3E1) / (0.2E1 * n - 0.7E1) * (n ** 2 - 0.1E1) / ( &
      &0.4E1 * n ** 2 - 0.1E1)
    b4m5 = cgret
    return
  end

end module b4func
