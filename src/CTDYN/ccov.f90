module ccov

  use cio 

  implicit none 

  integer, dimension (nt) :: indeg     ! index the eigenvalues      
  real(dp), dimension (nt) :: inde     ! index the eigenvalues      
  real(dp), dimension (nt, nt) :: vr
  complex(qp), dimension (nt, nt) :: cvr
  complex(qp), dimension (nt)  :: ww
  real(dp), dimension (nt) :: wr, wi  ! eigenvalues
  
end module ccov 
