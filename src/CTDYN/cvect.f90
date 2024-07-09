module cvect

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
  
  ! coefficients of cos(theta)^n * \partial (p1n * sin(theta) / \partial theta
  ! n=1
  real :: b1m2,b1p2,b10
  external b1m2,b1p2,b10      
  ! n=2
  real :: b2m3,b2p3,b2m1,b2p1
  external b2m3,b2p3,b2m1,b2p1
  ! n=3
  real :: b3m4,b3p4,b3m2,b3p2,b30
  external b3m4,b3p4,b3m2,b3p2,b30
  ! n=4
  real :: b4m5,b4p5,b4m3,b4p3,b4m1,b4p1
  external b4m5,b4p5,b4m3,b4p3,b4m1,b4p1
  
  real :: b5m6,b5p6,b5m4,b5p4,b5m2,b5p2,b50
  external b5m6,b5p6,b5m4,b5p4,b5m2,b5p2,b50  

  ! coefficient functions of cos(theta)^n *p1n   
  real :: c2m2,c2p2,c20
  external c2m2,c2p2,c20
  
  real :: c3m3,c3p3,c3m1,c3p1
  external c3m3,c3p3,c3m1,c3p1
  
  real :: c4m4,c4p4,c4m2,c4p2,c40
  external c4m4,c4p4,c4m2,c4p2,c40
   
  real :: c6m6,c6p6,c6m4,c6p4,c6m2,c6p2,c60
  external c6m6,c6p6,c6m4,c6p4,c6m2,c6p2,c60
  
  ! coefficient functions of cos *sin^2 p1n
  real :: c1s2p3,c1s2m3,c1s2p1,c1s2m1
  external c1s2p3,c1s2m3,c1s2p1,c1s2m1
  
  real :: c3s2p5,c3s2m5,c3s2p3,c3s2m3,c3s2p1,c3s2m1
  external c3s2p5,c3s2m5,c3s2p3,c3s2m3,c3s2p1,c3s2m1
  
  real :: d1p2,d1m2,d10
  external d1p2,d1m2,d10
  
  real :: d3p4,d3m4,d3p2,d3m2,d30
  external d3p4,d3m4,d3p2,d3m2,d30
  
  real :: d5p6,d5m6,d5p4,d5m4,d5p2,d5m2
  external d5p6,d5m6,d5p4,d5m4,d5p2,d5m2

end module cvect  


