module cvect

  use cio
  implicit none

  !########## definition of matrix arrays
  
  complex(qp) :: ala(np,nb)              ! alpha and gamma
  complex(qp) :: gaa(np,nb)              
  complex(qp) :: bea(np,nb)              ! array n-dep beta 
  complex(qp) :: alb(np,nb)              ! alpha and gamma
  complex(qp) :: gab(np,nb)              
  complex(qp) :: beb(np,nb)               ! array n-dep beta 
  
  complex(qp) :: fanp2(np, nb)
  complex(qp) :: fanm2(np, nb)
  complex(qp) :: fanp2p(np,nb)
  complex(qp) :: fanm2p(np,nb)
  complex(qp) :: fbnp2(np, nb)
  complex(qp) :: fbnm2(np, nb)
  complex(qp) :: fbnp2p(np,nb)
  complex(qp) :: fbnm2p(np,nb)
  
  complex(qp) :: fanp4(np, nb)
  complex(qp) :: fanm4(np, nb)
  complex(qp) :: fanp4p(np,nb)
  complex(qp) :: fanm4p(np,nb)
  complex(qp) :: fbnp4(np, nb)
  complex(qp) :: fbnm4(np, nb)
  complex(qp) :: fbnp4p(np,nb)
  complex(qp) :: fbnm4p(np,nb)
  
  complex(qp) :: fanp6(np, nb)
  complex(qp) :: fanm6(np, nb)
  complex(qp) :: fanp6p(np,nb)
  complex(qp) :: fanm6p(np,nb)
  complex(qp) :: fbnp6(np, nb)
  complex(qp) :: fbnm6(np, nb)
  complex(qp) :: fbnp6p(np,nb)
  complex(qp) :: fbnm6p(np,nb)
  
  complex(qp) :: omep1(np, nb)
  complex(qp) :: omem1(np, nb)
  complex(qp) :: omep1p(np,nb)
  complex(qp) :: omem1p(np,nb)
  complex(qp) :: omep3(np, nb)
  complex(qp) :: omem3(np, nb)
  complex(qp) :: omep3p(np,nb)
  complex(qp) :: omem3p(np,nb)
  complex(qp) :: omep5(np, nb)
  complex(qp) :: omem5(np, nb)
  complex(qp) :: omep5p(np,nb)
  complex(qp) :: omem5p(np,nb)
  
  complex(qp) :: alqm1(np,nb)
  complex(qp) :: alqp1(np,nb)
  complex(qp) :: gaqm1(np,nb)
  complex(qp) :: gaqp1(np,nb)
  complex(qp) :: beqm1(np,nb)
  complex(qp) :: beqp1(np,nb)
  
  complex(qp) :: alqm3(np,nb)
  complex(qp) :: alqp3(np,nb)
  complex(qp) :: gaqm3(np,nb)
  complex(qp) :: gaqp3(np,nb)
  complex(qp) :: beqm3(np,nb)
  complex(qp) :: beqp3(np,nb)
  
end module cvect  


