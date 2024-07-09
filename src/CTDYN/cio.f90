module cio

  integer :: np ! mesh points
  integer :: nt ! matrix dimension 
  integer :: na !     na = 1,3,5,... radial order in an and an
  integer :: nb !     nb = 2,4,6,... radial order in bn and bn
  
  parameter(na = 49)     ! harmonics
  parameter(np = 50)     ! mesh points 
  
  parameter(nb = na+1)           ! harmonics b number
  parameter(nt = nb*np)          ! global matrix dimension
  
  integer :: lwork
  parameter(lwork=540000)      
  
  integer :: nueg
  parameter(nueg=8)
  
  integer :: n_theta
  parameter(n_theta=301) 
  
  integer :: n_time
  parameter(n_time=300)
  
  real :: pi
  parameter(pi=3.14159265359)

end module cio
