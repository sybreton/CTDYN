module c1s2func

  implicit none

contains

  real function c1s2p3 (n)
    integer n
    real cgret
    cgret = -n * (n + 0.1E1) * (n + 0.2E1) / (0.2E1 * n + 0.1E1) / (0
      &.2E1 * n + 0.3E1) / (0.2E1 * n + 0.5E1)
    c1s2p3 = cgret
    return
  end
  
  real function c1s2m3 (n)
    integer n
    real cgret
    cgret = -n * (n + 0.1E1) * (n - 0.1E1) / (0.2E1 * n - 0.1E1) / (0
      &.2E1 * n + 0.1E1) / (0.2E1 * n - 0.3E1)
    c1s2m3 = cgret
    return
  end
  
  real function c1s2p1 (n)
    integer n
    real cgret
    cgret = n * (n ** 2 + 0.2E1 * n + 0.1E1) / (0.2E1 * n + 0.5E1) / 
      &(0.4E1 * n ** 2 - 0.1E1)
    c1s2p1 = cgret
    return
  end
  
  real function c1s2m1 (n)
    integer n
    real cgret
    cgret = (n + 0.1E1) * n ** 2 / (0.2E1 * n + 0.1E1) / (0.4E1 * n *
      &* 2 - 0.9E1)
    c1s2m1 = cgret
    return
  end

end module c1s2func
