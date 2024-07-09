module cvect

  use cio
  implicit none

  !########## definition of matrix arrays
  
  complex*16 :: ala(np,nb)              ! alpha and gamma
  complex*16 :: gaa(np,nb)              
  complex*16 :: bea(np,nb)              ! array n-dep beta 
  complex*16 :: alb(np,nb)              ! alpha and gamma
  complex*16 :: gab(np,nb)              
  complex*16 :: beb(np,nb)               ! array n-dep beta 
  
  complex*16 :: fanp2(np, nb)
  complex*16 :: fanm2(np, nb)
  complex*16 :: fanp2p(np,nb)
  complex*16 :: fanm2p(np,nb)
  complex*16 :: fbnp2(np, nb)
  complex*16 :: fbnm2(np, nb)
  complex*16 :: fbnp2p(np,nb)
  complex*16 :: fbnm2p(np,nb)
  
  complex*16 :: fanp4(np, nb)
  complex*16 :: fanm4(np, nb)
  complex*16 :: fanp4p(np,nb)
  complex*16 :: fanm4p(np,nb)
  complex*16 :: fbnp4(np, nb)
  complex*16 :: fbnm4(np, nb)
  complex*16 :: fbnp4p(np,nb)
  complex*16 :: fbnm4p(np,nb)
  
  complex*16 :: fanp6(np, nb)
  complex*16 :: fanm6(np, nb)
  complex*16 :: fanp6p(np,nb)
  complex*16 :: fanm6p(np,nb)
  complex*16 :: fbnp6(np, nb)
  complex*16 :: fbnm6(np, nb)
  complex*16 :: fbnp6p(np,nb)
  complex*16 :: fbnm6p(np,nb)
  
  complex*16 :: omep1(np, nb)
  complex*16 :: omem1(np, nb)
  complex*16 :: omep1p(np,nb)
  complex*16 :: omem1p(np,nb)
  complex*16 :: omep3(np, nb)
  complex*16 :: omem3(np, nb)
  complex*16 :: omep3p(np,nb)
  complex*16 :: omem3p(np,nb)
  complex*16 :: omep5(np, nb)
  complex*16 :: omem5(np, nb)
  complex*16 :: omep5p(np,nb)
  complex*16 :: omem5p(np,nb)
  
  complex*16 :: alqm1(np,nb)
  complex*16 :: alqp1(np,nb)
  complex*16 :: gaqm1(np,nb)
  complex*16 :: gaqp1(np,nb)
  complex*16 :: beqm1(np,nb)
  complex*16 :: beqp1(np,nb)
  
  complex*16 :: alqm3(np,nb)
  complex*16 :: alqp3(np,nb)
  complex*16 :: gaqm3(np,nb)
  complex*16 :: gaqp3(np,nb)
  complex*16 :: beqm3(np,nb)
  complex*16 :: beqp3(np,nb)
  
end module cvect  


