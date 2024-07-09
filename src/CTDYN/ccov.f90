module ccov

  use cio 
  implicit none 

  integer, dimension (nt) :: indeg     ! index the eigenvalues      
  real, dimension (nt) :: inde     ! index the eigenvalues      
  real, dimension (nt, nt) :: vr
  complex*16, dimension (nt, nt) :: cvr
  complex*16, dimension (nt)  :: ww
  real, dimension (nt) :: wr, wi  ! eigenvalues
  common/vec1/vr,indeg,wr,wi
  
end module ccov 
