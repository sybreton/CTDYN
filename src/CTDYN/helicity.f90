! compute the helicity if beta not zero
! real part
esz1=0
! img part
esz2=0

! compute  boundary term 
! real part
bsz1=0
! img part
bsz2=0

if (beta.ne.0) then 
  do k2 = nsp,nep,2          !plm loop 
    xnu=k2*1.d0+1.0/2.0
    call bess(beta, k2, ffc)
    vc(np+2,k2,1)  = (4*vc(np+1,k2,1)-vc(np,k2,1))/(2*hh*(k2+1+ffc)+3)  
    vc(np+2,k2,2)  = (4*vc(np+1,k2,2)-vc(np,k2,2))/(2*hh*(k2+1+ffc)+3)    
    do j =np+2,np+2+nf    ! radial loop, perform the integral to obtain helicity and bulk term of the total energy
      x = x1+(j-1)*hh
      ! load the a_n vectors: use bad approx for integral to stay simple... 
      call bessjy(x*beta,xnu,rj,ry,rjp,ryp)
      call bessjy(beta,xnu,rj1,ry1,rjp1,ryp1)
      esz1 = esz1+k2*(k2+1)*(vc(np+2,k2,1)*(gam*ry+rj)/sqrt(x)/(gam*ry1+rj1))**2*x**2*hh
      esz2 = esz1+k2*(k2+1)*(vc(np+2,k2,2)*(gam*ry+rj)/sqrt(x)/(gam*ry1+rj1))**2*x**2*hh
   enddo ! radial
 ! boundary term eq. 19 berger
 enddo ! plm   
   !surface  term (nota che il rapporto dei gamma dovrebbe essere =1) ricontrolla normalization
   bsz1 = bsz1+k2*(k2+1) 
endif
heli= 2*beta*(esz1+esz2)
write(22,*) beta, heli


