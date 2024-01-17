
C compute the helicity if beta not zero

! real part
      ESZ1=0
! img part
      ESZ2=0

C compute  boundary term 
 
! real part
      BSZ1=0
! img part
      BSZ2=0

      if(beta.ne.0)then 
      DO   k2 = nsp,nep,2          !Plm loop 

      xnu=k2*1.d0+1.0/2.0
      call bess(beta, k2, ffc)
      VC(NP+2,k2,1)  = (4*VC(NP+1,k2,1)-VC(NP,k2,1))/(2*hh*(k2+1+ffc)+3)  
      VC(NP+2,k2,2)  = (4*VC(NP+1,k2,2)-VC(NP,k2,2))/(2*hh*(k2+1+ffc)+3)    


      DO J =NP+2,NP+2+NF    ! radial loop, perform the integral to obtain helicity and bulk term of the total energy
      x = x1+(J-1)*hh

               ! LOAD THE a_n vectors: use bad approx for integral to stay simple... 

                   call bessjy(x*beta,xnu,rj,ry,rjp,ryp)
                   call bessjy(beta,xnu,rj1,ry1,rjp1,ryp1)

                   ESZ1 = ESZ1+k2*(k2+1)*(VC(NP+2,k2,1)*(gam*ry+rj)/sqrt(x)/(gam*ry1+rj1))**2*x**2*hh
                   ESZ2 = ESZ1+k2*(k2+1)*(VC(NP+2,k2,2)*(gam*ry+rj)/sqrt(x)/(gam*ry1+rj1))**2*x**2*hh

         ENDDO ! Radial

        ! boundary term eq. 19 berger

c         print*, x, k2
c         stop

         ENDDO ! Plm   

         !surface  term (nota che il rapporto dei gamma dovrebbe essere =1) RICONTROLLA NORMALIZATION

         BSZ1 = BSZ1+k2*(k2+1) 

         endif

         heli= 2*beta*(ESZ1+ESZ2)

         write(22,*) beta, heli


