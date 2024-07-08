real function d1p2 (n)
  integer n
  real cgret
  real x_in, bct, c3,mmm
  common/ppar/x_in,bct,c3,mmm
  cgret = (n-1.0)*n*(n-abs(mmm)+1)*(n-abs(mmm)+2)/(6+19.0*n+16.0*n**2+4.0*n**3) 
  d1p2 = cgret
  return
end

real function d1m2 (n)
  integer n
  real cgret
  real x_in, bct, c3,mmm
  common/ppar/x_in,bct,c3,mmm
  cgret = -(n+1.0)*(n+2.0)*(n+abs(mmm)-1)*(n+abs(mmm))/(1.0-n-4.0*n**2+4.0*n**3)
  d1m2 = cgret
  return
end

real function d10 (n)
  integer n
  real cgret
  real x_in, bct, c3,mmm
  common/ppar/x_in,bct,c3,mmm
  cgret = -0.3e1 * (n + 0.2e1) * (n - 0.1e1) * (n*(n+1)-3*mmm**2)/ (0.2e1 * n - 0.1e1) 
     &/ (0.2e1 * n + 0.3e1)/n/(n+1.0)
  d10 = cgret
  return
end

real function d3p4 (n)
  integer n
  real cgret
  cgret = n / (0.2e1 * n + 0.1e1) / (0.2e1 * n + 0.3e1) * (n + 0.2e
     &1) * (n + 0.3e1) / (0.2e1 * n + 0.5e1) / (0.2e1 * n + 0.7e1) * (n *
     &* 2 - 0.1e1)
  d3p4 = cgret
  return
end

real function d3m4 (n)
  integer n
  real cgret
  cgret = -n / (0.2e1 * n - 0.5e1) / (0.2e1 * n - 0.3e1) * (n ** 2 
     &- 0.1e1) * (n ** 2 - 0.4e1) / (0.4e1 * n ** 2 - 0.1e1)
  d3m4 = cgret
  return
end

real function d3p2 (n)
  integer n
  real cgret
  cgret = (0.2e1 * n ** 2 + 0.3e1 * n - 0.17e2) * n / (0.2e1 * n + 
     &0.7e1) / (0.2e1 * n + 0.3e1) * (n ** 2 - 0.1e1) / (0.4e1 * n ** 2 -
     & 0.1e1)
  d3p2 = cgret
  return
end

real function d3m2 (n)
  integer n
  real cgret
  cgret = -(0.2e1 * n ** 2 + n - 0.18e2) * (n + 0.2e1) * (n + 0.1e1
     &) * n / (0.2e1 * n + 0.3e1) / (0.2e1 * n - 0.5e1) / (0.4e1 * n ** 2
     & - 0.1e1)
  d3m2 = cgret
  return
end

real function d30 (n)
  integer n
  real cgret
  cgret = -0.9e1 * (n ** 2 + n - 0.5e1) * (n - 0.1e1) * (n + 0.2e1)
     & / (0.2e1 * n + 0.5e1) / (0.2e1 * n - 0.1e1) / (0.4e1 * n ** 2 - 0.
     &9e1)
  d30 = cgret
  return
end

real function d5p6 (n)
  integer n
  real cgret
  cgret = n / (0.2e1 * n + 0.1e1) / (0.2e1 * n + 0.3e1) * (n + 0.2e
     &1) * (n + 0.3e1) / (0.2e1 * n + 0.5e1) / (0.2e1 * n + 0.7e1) * (n *
     &* 2 - 0.1e1) * (n + 0.4e1) * (n + 0.5e1) / (0.2e1 * n + 0.9e1) / (0
     &.2e1 * n + 0.11e2)
  d5p6 = cgret
  return
end

real function d5m6 (n)
  integer n
  real cgret
  cgret = -(n - 0.4e1) * (n - 0.3e1) / (0.2e1 * n - 0.9e1) / (0.2e1
     & * n - 0.7e1) * n / (0.2e1 * n - 0.5e1) / (0.2e1 * n - 0.3e1) * (n 
     &** 2 - 0.1e1) * (n ** 2 - 0.4e1) / (0.4e1 * n ** 2 - 0.1e1)
  d5m6 = cgret
  return
end

real function d5p4 (n)
  integer n
  real cgret
  cgret = n * (n + 0.2e1) * (n + 0.3e1) * (n ** 2 - 0.1e1) * (0.4e1
     & * n ** 2 + 0.17e2 * n - 0.32e2) / (0.2e1 * n + 0.11e2) / (0.2e1 * 
     &n + 0.5e1) / (0.2e1 * n + 0.3e1) / (0.2e1 * n + 0.7e1) / (0.4e1 * n
     & ** 2 - 0.1e1)
  d5p4 = cgret
  return
end

real function d5m4 (n)
  integer n
  real cgret
  cgret = -(0.4e1 * n ** 6 - 0.9e1 * n ** 5 - 0.65e2 * n ** 4 + 0.4
     &5e2 * n ** 3 + 0.241e3 * n ** 2 - 0.36e2 * n - 0.180e3) * n / (0.2e
     &1 * n - 0.9e1) / (0.2e1 * n - 0.5e1) / (0.4e1 * n ** 2 - 0.1e1) / (
     &0.4e1 * n ** 2 - 0.9e1)
  d5m4 = cgret
  return
end

real function d5p2 (n)
  integer n
  real cgret
  cgret = 0.5e1 * n * (n ** 6 + 0.3e1 * n ** 5 - 0.23e2 * n ** 4 - 
     &0.45e2 * n ** 3 + 0.139e3 * n ** 2 + 0.42e2 * n - 0.117e3) / (0.2e1
     & * n + 0.9e1) / (0.4e1 * n ** 2 - 0.9e1) / (0.2e1 * n + 0.7e1) / (0
     &.4e1 * n ** 2 - 0.1e1)
  d5p2 = cgret
  return
end

real function d5m2 (n)
  integer n
  real cgret
  cgret = -0.5e1 * n * (n ** 6 + 0.4e1 * n ** 5 - 0.20e2 * n ** 4 -
     & 0.80e2 * n ** 3 + 0.64e2 * n ** 2 + 0.391e3 * n + 0.270e3) / (0.2e
     &1 * n + 0.3e1) / (0.2e1 * n - 0.7e1) / (0.4e1 * n ** 2 - 0.1e1) / (
     &0.4e1 * n ** 2 - 0.25e2)
  d5m2 = cgret
  return
end

real function d50 (n)
  integer n
  real cgret
  cgret = -0.15e2 * (n + 0.2e1) * (0.2e1 * n ** 5 + 0.2e1 * n ** 4 
     &- 0.32e2 * n ** 3 - 0.2e1 * n ** 2 + 0.135e3 * n - 0.105e3) / (0.2e
     &1 * n - 0.1e1) / (0.2e1 * n + 0.7e1) / (0.4e1 * n ** 2 - 0.25e2) / 
     &(0.4e1 * n ** 2 - 0.9e1)
  d50 = cgret
  return
end

real function c2p2 (n)
  integer n
  real cgret
  real x_in, bct, c3,mmm
  common/ppar/x_in,bct,c3,mmm
  cgret = n*(1.e0-abs(mmm)+n)*(2.e0 -abs(mmm)+n) / (6.e0 +19.e0*n+16.e0*n**2+4*n**3) 
  c2p2 = cgret
  return
end

real function c2m2 (n)
  integer n
  real cgret
  real x_in, bct, c3,mmm
  common/ppar/x_in,bct,c3,mmm
  cgret  = (1.0+n)*(-1.0+abs(mmm)+n)*(abs(mmm)+n)/(1.0-n-4.0*n**2+4.0*n**3)
  c2m2 = cgret
  return
end

real function c20 (n)
  integer n
  real cgret
  real x_in, bct, c3,mmm
  common/ppar/x_in,bct,c3,mmm
  cgret = ((n*(n+1)-3*mmm**2)*2*(n*(n+1.)-3)/n/(n+1.)/(2*n-1.)/(2*n+3.) +1)/3.
  cgret = (0.2e1 * n ** 2 + 0.2e1 * n - 0.3e1) / (0.2e1 * n - 0.1e1
     &) / (0.2e1 * n + 0.3e1)
  c20 = cgret
  return
end

real function c3p3 (n)
  integer n
  real cgret
  cgret = n * (n + 0.1e1) * (n + 0.2e1) / (0.2e1 * n + 0.1e1) / (0.
     &2e1 * n + 0.3e1) / (0.2e1 * n + 0.5e1)
  c3p3 = cgret
  return
end

real function c3m3 (n)
  integer n
  real cgret
  cgret = n * (n + 0.1e1) * (n - 0.1e1) / (0.2e1 * n - 0.1e1) / (0.
     &2e1 * n + 0.1e1) / (0.2e1 * n - 0.3e1)
  c3m3 = cgret
  return
end

real function c3p1 (n)
  integer n
  real cgret
  cgret = 0.3e1 * n * (n ** 2 + 0.2e1 * n - 0.2e1) / (0.2e1 * n + 0
     &.5e1) / (0.2e1 * n - 0.1e1) / (0.2e1 * n + 0.1e1)
  c3p1 = cgret
  return
end

real function c3m1 (n)
  integer n
  real cgret
  cgret = 0.3e1 * (n ** 2 - 0.3e1) * (n + 0.1e1) / (0.2e1 * n + 0.3
     &e1) / (0.2e1 * n - 0.3e1) / (0.2e1 * n + 0.1e1)
  c3m1 = cgret
  return
end

real function c4p4 (n)
  integer n
  real cgret
  cgret = n * (n + 0.1e1) * (n + 0.2e1) / (0.2e1 * n + 0.1e1) / (0.
     &2e1 * n + 0.3e1) / (0.2e1 * n + 0.5e1) * (n + 0.3e1) / (0.2e1 * n +
     & 0.7e1)
  c4p4 = cgret
  return
end

real function c4m4 (n)
  integer n
  real cgret
  cgret = n * (n + 0.1e1) * (n - 0.1e1) / (0.2e1 * n - 0.1e1) / (0.
     &2e1 * n + 0.1e1) / (0.2e1 * n - 0.3e1) * (n - 0.2e1) / (0.2e1 * n -
     & 0.5e1)
  c4m4 = cgret
  return
end

real function c4p2 (n)
  integer n
  real cgret
  cgret = 0.2e1 * (0.2e1 * n ** 2 + 0.6e1 * n - 0.5e1) * n * (n + 0
     &.1e1) / (0.2e1 * n + 0.7e1) / (0.2e1 * n + 0.3e1) / (0.4e1 * n ** 2
     & - 0.1e1)
  c4p2 = cgret
  return
end

real function c4m2 (n)
  integer n
  real cgret
  cgret = 0.2e1 * (0.2e1 * n ** 2 - 0.2e1 * n - 0.9e1) * n * (n + 0
     &.1e1) / (0.2e1 * n + 0.3e1) / (0.2e1 * n - 0.5e1) / (0.4e1 * n ** 2
     & - 0.1e1)
  c4m2 = cgret
  return
end

real function c40 (n)
  integer n
  real cgret
  cgret = 0.3e1 * (0.2e1 * n ** 4 + 0.4e1 * n ** 3 - 0.10e2 * n ** 
     &2 - 0.12e2 * n + 0.15e2) / (0.2e1 * n - 0.1e1) / (0.2e1 * n + 0.5e1
     &) / (0.4e1 * n ** 2 - 0.9e1)
  c40 = cgret
  return
end

real function c6p6 (n)
  integer n
  real cgret
  cgret = n * (n + 0.1e1) * (n + 0.2e1) / (0.2e1 * n + 0.1e1) / (0.
     &2e1 * n + 0.3e1) / (0.2e1 * n + 0.5e1) * (n + 0.3e1) / (0.2e1 * n +
     & 0.7e1) * (n + 0.4e1) * (n + 0.5e1) / (0.2e1 * n + 0.9e1) / (0.2e1 
     &* n + 0.11e2)
  c6p6 = cgret
  return
end

real function c6m6 (n)
  integer n
  real cgret
  cgret = n * (n + 0.1e1) * (n - 0.1e1) / (0.2e1 * n - 0.1e1) / (0.
     &2e1 * n + 0.1e1) / (0.2e1 * n - 0.3e1) * (n - 0.2e1) / (0.2e1 * n -
     & 0.5e1) * (n - 0.4e1) * (n - 0.3e1) / (0.2e1 * n - 0.9e1) / (0.2e1 
     &* n - 0.7e1)
  c6m6 = cgret
  return
end

real function c6p4 (n)
  integer n
  real cgret
  cgret = 0.3e1 * (n + 0.1e1) * n * (n + 0.2e1) * (n + 0.3e1) * (0.
     &2e1 * n ** 2 + 0.10e2 * n - 0.7e1) / (0.4e1 * n ** 2 - 0.1e1) / (0.
     &2e1 * n + 0.11e2) / (0.2e1 * n + 0.7e1) / (0.2e1 * n + 0.5e1) / (0.
     &2e1 * n + 0.3e1)
  c6p4 = cgret
  return
end

real function c6m4 (n)
  integer n
  real cgret
  cgret = 0.3e1 * (n - 0.2e1) * n * (0.2e1 * n ** 2 - 0.6e1 * n - 0
     &.15e2) / (0.2e1 * n - 0.9e1) / (0.2e1 * n - 0.5e1) * (n ** 2 - 0.1e
     &1) / (0.4e1 * n ** 2 - 0.1e1) / (0.4e1 * n ** 2 - 0.9e1)
  c6m4 = cgret
  return
end

real function c6p2 (n)
  integer n
  real cgret
  cgret = 0.15e2 * (n ** 4 + 0.6e1 * n ** 3 - n ** 2 - 0.30e2 * n +
     & 0.21e2) * n * (n + 0.1e1) / (0.4e1 * n ** 2 - 0.9e1) / (0.2e1 * n 
     &+ 0.9e1) / (0.2e1 * n + 0.7e1) / (0.4e1 * n ** 2 - 0.1e1)
  c6p2 = cgret
  return
end

real function c6m2 (n)
  integer n
  real cgret
  cgret = 0.5e1 * (n + 0.1e1) * n * (n ** 4 - 0.2e1 * n ** 3 - 0.13
     &e2 * n ** 2 + 0.14e2 * n + 0.45e2) / (0.2e1 * n + 0.3e1) / (0.2e1 *
     & n - 0.7e1) / (0.4e1 * n ** 2 - 0.25e2) / (0.4e1 * n ** 2 - 0.1e1)
  c6m2 = cgret
  return
end

real function c60 (n)
  integer n
  real cgret
  cgret = 0.5e1 * (0.4e1 * n ** 6 + 0.12e2 * n ** 5 - 0.50e2 * n **
     & 4 - 0.120e3 * n ** 3 + 0.208e3 * n ** 2 + 0.270e3 * n - 0.315e3) /
     & (0.4e1 * n ** 2 - 0.9e1) / (0.2e1 * n - 0.1e1) / (0.2e1 * n + 0.7e
     &1) / (0.4e1 * n ** 2 - 0.25e2)
  c60 = cgret
  return
end

real function  abg(x)
  real x
  real get
  if(x.ne.0)then
    get = abs(x)/x
  else
    get = 0
  endif
  abg = get
  return
end
