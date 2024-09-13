module cvect

  use kind_parameters
  use cio
  implicit none

  !########## definition of matrix arrays
  
  complex(dp) :: ala(np,nb)              ! alpha and gamma
  complex(dp) :: gaa(np,nb)              
  complex(dp) :: bea(np,nb)              ! array n-dep beta 
  complex(dp) :: alb(np,nb)              ! alpha and gamma
  complex(dp) :: gab(np,nb)              
  complex(dp) :: beb(np,nb)               ! array n-dep beta 
  
  complex(dp) :: fanp2(np, nb)
  complex(dp) :: fanm2(np, nb)
  complex(dp) :: fanp2p(np,nb)
  complex(dp) :: fanm2p(np,nb)
  complex(dp) :: fbnp2(np, nb)
  complex(dp) :: fbnm2(np, nb)
  complex(dp) :: fbnp2p(np,nb)
  complex(dp) :: fbnm2p(np,nb)
  
  complex(dp) :: fanp4(np, nb)
  complex(dp) :: fanm4(np, nb)
  complex(dp) :: fanp4p(np,nb)
  complex(dp) :: fanm4p(np,nb)
  complex(dp) :: fbnp4(np, nb)
  complex(dp) :: fbnm4(np, nb)
  complex(dp) :: fbnp4p(np,nb)
  complex(dp) :: fbnm4p(np,nb)
  
  complex(dp) :: fanp6(np, nb)
  complex(dp) :: fanm6(np, nb)
  complex(dp) :: fanp6p(np,nb)
  complex(dp) :: fanm6p(np,nb)
  complex(dp) :: fbnp6(np, nb)
  complex(dp) :: fbnm6(np, nb)
  complex(dp) :: fbnp6p(np,nb)
  complex(dp) :: fbnm6p(np,nb)
  
  complex(dp) :: omep1(np, nb)
  complex(dp) :: omem1(np, nb)
  complex(dp) :: omep1p(np,nb)
  complex(dp) :: omem1p(np,nb)
  complex(dp) :: omep3(np, nb)
  complex(dp) :: omem3(np, nb)
  complex(dp) :: omep3p(np,nb)
  complex(dp) :: omem3p(np,nb)
  complex(dp) :: omep5(np, nb)
  complex(dp) :: omem5(np, nb)
  complex(dp) :: omep5p(np,nb)
  complex(dp) :: omem5p(np,nb)
  
  complex(dp) :: alqm1(np,nb)
  complex(dp) :: alqp1(np,nb)
  complex(dp) :: gaqm1(np,nb)
  complex(dp) :: gaqp1(np,nb)
  complex(dp) :: beqm1(np,nb)
  complex(dp) :: beqp1(np,nb)
  
  complex(dp) :: alqm3(np,nb)
  complex(dp) :: alqp3(np,nb)
  complex(dp) :: gaqm3(np,nb)
  complex(dp) :: gaqp3(np,nb)
  complex(dp) :: beqm3(np,nb)
  complex(dp) :: beqp3(np,nb)
  
end module cvect  


