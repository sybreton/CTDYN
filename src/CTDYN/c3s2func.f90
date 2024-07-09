module c3s2func

  implicit none

contains

  real function c3s2p5 (n)
    integer n
    real cgret
    cgret = -(n + 0.3E1) * (n + 0.4E1) * n * (n + 0.1E1) * (n + 0.2E1
       &) / (0.2E1 * n + 0.1E1) / (0.2E1 * n + 0.3E1) / (0.2E1 * n + 0.5E1)
       & / (0.2E1 * n + 0.7E1) / (0.2E1 * n + 0.9E1)
    c3s2p5 = cgret
    return
  end
  
  real function c3s2p3 (n)
    integer n
    real cgret
    cgret = -(n + 0.2E1) * (n + 0.1E1) * n * (n ** 2 + 0.4E1 * n - 0.
       &6E1) / (0.4E1 * n ** 2 - 0.1E1) / (0.2E1 * n + 0.9E1) / (0.2E1 * n 
       &+ 0.5E1) / (0.2E1 * n + 0.3E1)
    c3s2p3 = cgret
    return
  end
  
  real function c3s2p1 (n)
    integer n
    real cgret
    cgret = (0.2E1 * n ** 4 + 0.8E1 * n ** 3 + n ** 2 - 0.14E2 * n - 
       &0.9E1) * n / (0.2E1 * n - 0.3E1) / (0.2E1 * n + 0.7E1) / (0.2E1 * n
       & + 0.5E1) / (0.4E1 * n ** 2 - 0.1E1)
    c3s2p1 = cgret
    return
  end
  
  real function c3s2m1 (n)
    integer n
    real cgret
    cgret = (n + 0.1E1) * n ** 2 * (0.2E1 * n ** 2 - 0.11E2) / (0.2E1
       & * n + 0.1E1) / (0.4E1 * n ** 2 - 0.25E2) / (0.4E1 * n ** 2 - 0.9E1
       &)
    c3s2m1 = cgret
    return
  end
  
  real function c3s2m3 (n)
    integer n
    real cgret
    cgret = -n * (n ** 2 - 0.2E1 * n - 0.9E1) / (0.2E1 * n - 0.7E1) *
       & (n ** 2 - 0.1E1) / (0.4E1 * n ** 2 - 0.1E1) / (0.4E1 * n ** 2 - 0.
       &9E1)
    c3s2m3 = cgret
    return
  end
  
  real function c3s2m5 (n)
    integer n
    real cgret
    cgret = -n / (0.2E1 * n - 0.3E1) * (n - 0.3E1) * (n - 0.2E1) / (0
       &.2E1 * n - 0.7E1) / (0.2E1 * n - 0.5E1) * (n ** 2 - 0.1E1) / (0.4E1
       & * n ** 2 - 0.1E1)
    c3s2m5 = cgret
    return
  end

end module c3s2func
