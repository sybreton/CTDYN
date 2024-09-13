module ccov

  use kind_parameters 
  use cio 

  implicit none 

  integer, dimension (nt) :: indeg     ! index the eigenvalues      
  real(dp), dimension (nt) :: inde     ! index the eigenvalues      
  real(dp), dimension (nt, nt) :: vr
  complex(dp), dimension (nt, nt) :: cvr
  complex(dp), dimension (nt)  :: ww
  real(dp), dimension (nt) :: wr, wi  ! eigenvalues
  
end module ccov 
