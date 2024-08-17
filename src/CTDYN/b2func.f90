module b2func 

  use cio 

  implicit none

contains 

  real(dp) function b0p1(m)
    integer :: m
    real(dp) :: cgret
  
    cgret = m*(m+1.e0)/(2*m+1.e0)
    b0p1= cgret
    return
  end 
  
  real(dp) function b0m1(m)
    integer :: m
    real(dp) :: cgret
  
    cgret = -m*(m+1.e0)/(2*m+1.e0)
    b0m1= cgret
    return
  end 
  
  real(dp) function b1p2 (n)
    integer :: n
    real(dp) :: cgret
  
    cgret = n*(1.0+n)*(1-abs(mmm)+n)*(2-abs(mmm)+n)/(6.0+19.0*n+16.0*n**2+4*n**3)
    b1p2 = cgret
    return
  end
  
  real(dp) function b1m2 (n)
    integer :: n
    real(dp) :: cgret
  
    cgret = -n ** 2 * (n + 0.1e1) / (0.2e1 * n + 0.1e1) / (0.2e1 * n - 0.1e1)
    b1m2 = cgret
    return
  end
  
  real(dp) function b10 (n)
    integer :: n
    real(dp) :: cgret
  
    cgret = (n * (n + 0.1e1)  -3.0*mmm**2) / (0.2e1 * n  - 0.1e1) / (0.2e1 * n + 0.3e1)
    b10 = cgret
    return
  end
  
  real(dp) function b2p3 (m)
    integer :: m
    real(dp) :: cgret
  
    cgret = m * (m + 0.1e1) ** 2 / (0.2e1 * m + 0.1e1) / (0.2e1 * m + 0.3e1) * (m + 0.2e1) / (0.2e1 * m + 0.5e1)
    b2p3 = cgret
    return
  end
  
  real(dp) function b2m3 (n)
    integer :: n
    real(dp) :: cgret
  
    cgret = -n ** 2 / (0.2e1 * n - 0.3e1) * (n ** 2 - 0.1e1) / (0.4e1 * n ** 2 - 0.1e1)
    b2m3 = cgret
    return
  end
  
  real(dp) function b2p1 (n)
    integer :: n
    real(dp) :: cgret
  
    cgret = (n + 0.1e1) * n * (n ** 2 + 0.3e1 * n - 0.1e1) / (0.4e1 * n ** 2 - 0.1e1) / (0.2e1 * n + 0.5e1)
    b2p1 = cgret
    return
  end
  
  real(dp) function b2m1 (n)
    integer :: n
    real(dp) ::cgret
  
    cgret = -(n + 0.1e1) * n * (n ** 2 - n - 0.3e1) / (0.2e1 * n + 0.1e1) / (0.4e1 * n ** 2 - 0.9e1)
    b2m1 = cgret
    return
  end

end module b2func
