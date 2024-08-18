module b5func

  use kind_parameters 

  implicit none

contains

  real(dp) function b5m6 (n)
    integer :: n
    real(dp) :: cgret
    cgret = -n ** 2 * (n - 0.2E1) / (0.2E1 * n - 0.5E1) / (0.2E1 * n & 
      &- 0.3E1) * (n - 0.4E1) * (n - 0.3E1) / (0.2E1 * n - 0.9E1) / (0.2E1 &
      & * n - 0.7E1) * (n ** 2 - 0.1E1) / (0.4E1 * n ** 2 - 0.1E1)
    b5m6 = cgret
    return
  end
  
  real(dp) function b5p6 (n)
    integer :: n
    real(dp) :: cgret
    cgret = (n + 0.4E1) * (n + 0.5E1) / (0.2E1 * n + 0.9E1) / (0.2E1 &
       &* n + 0.11E2) * n * (n + 0.1E1) ** 2 / (0.2E1 * n + 0.1E1) / (0.2E1 &
       & * n + 0.3E1) * (n + 0.2E1) * (n + 0.3E1) / (0.2E1 * n + 0.5E1) / ( &
       &0.2E1 * n + 0.7E1)
    b5p6 = cgret
    return
  end
  
  real(dp) function b5m4 (n)
    integer :: n
    real(dp) :: cgret
    cgret = -(0.4E1 * n ** 2 - 0.13E2 * n - 0.27E2) * (n - 0.2E1) * n &
      & ** 2 / (0.4E1 * n ** 2 - 0.1E1) / (0.2E1 * n - 0.9E1) / (0.2E1 * n &
      & - 0.5E1) * (n ** 2 - 0.1E1) / (0.4E1 * n ** 2 - 0.9E1) 
    b5m4 = cgret
    return
  end
  
  real(dp) function b5p4 (n)
    integer :: n
    real(dp) :: cgret
    cgret = (n + 0.3E1) * (n + 0.2E1) * n * (n + 0.1E1) ** 2 * (0.4E1 &
      & * n ** 2 + 0.21E2 * n - 0.10E2) / (0.2E1 * n + 0.11E2) / (0.2E1 * & 
      &n + 0.5E1) / (0.2E1 * n + 0.3E1) / (0.2E1 * n + 0.7E1) / (0.4E1 * n &
      & ** 2 - 0.1E1)
    b5p4 = cgret
    return
  end
  
  real(dp) function b5m2 (n)
    integer :: n
    real(dp) :: cgret
    cgret = -0.5E1 * n ** 2 * (n + 0.1E1) * (n ** 4 - 0.3E1 * n ** 3 &
      &- 0.11E2 * n ** 2 + 0.21E2 * n + 0.37E2) / (0.2E1 * n + 0.3E1) / (0.4E1 &
      & * n ** 2 - 0.1E1) / (0.2E1 * n - 0.7E1) / (0.4E1 * n ** 2 - 0.25E2)
    b5m2 = cgret
    return
  end
  
  real(dp) function b5p2 (n)
    integer :: n
    real(dp) :: cgret
    cgret = 0.5E1 * (n + 0.1E1) ** 2 * n * (n ** 4 + 0.7E1 * n ** 3 + &
      & 0.4E1 * n ** 2 - 0.30E2 * n + 0.9E1) / (0.2E1 * n + 0.9E1) / (0.4E1 &
      & * n ** 2 - 0.1E1) / (0.2E1 * n + 0.7E1) / (0.4E1 * n ** 2 - 0.9E1 )
    b5p2 = cgret
    return
  end
  
  real(dp) function b50 (n)
    integer :: n
    real(dp) :: cgret
    cgret = 0.5E1 * n * (n + 0.1E1) * (0.2E1 * n ** 4 + 0.4E1 * n ** 3 &
      & - 0.20E2 * n ** 2 - 0.22E2 * n + 0.45E2) / (0.2E1 * n + 0.7E1) / & 
      &(0.4E1 * n ** 2 - 0.9E1) / (0.2E1 * n - 0.1E1) / (0.4E1 * n ** 2 - &
      &0.25E2)
    b50 = cgret
    return
  end

end module b5func
