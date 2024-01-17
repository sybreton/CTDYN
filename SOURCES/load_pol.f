
      theta1 = 0.
      theta2 = 3.1415926
      theta = theta1-theta2/(n_theta)

      DO K1=1, n_theta   ! theta-loop
       theta = theta + theta2/(n_theta)
       DO J =1,NP+2+NF    ! radial loop
        x = x1+(J-1)*hh
               ! LOAD THE a_n vectors for DIPOLES: THIS DEFINES THE POLOIDAL (potential)
        DO   k2 = nsp,nep,2          !Plm loop 

          call bess(beta, k2, ffc)

           VC(NP+2,k2,1)   = (4*VC(NP+1,k2,1)-VC(NP,k2,1))/(2*hh*(k2+1+ffc)+3)  
           VC(NP+2,k2,2)   = (4*VC(NP+1,k2,2)-VC(NP,k2,2))/(2*hh*(k2+1+ffc)+3)    

             if(J.ge.(NP+2))then
              if(beta.ne.0)then
               ! check the x infront of the definition as for the beta=0 case
                xnu=k2*1.d0+1.0/2.0
                call bessjy(x*beta,xnu,rj,ry,rjp,ryp)
                call bessjy(beta,xnu,rj1,ry1,rjp1,ryp1)
                SZ1 = VC(NP+2,k2,1)*(-plgndr(k2,1, cos(theta)))*sin(theta)*(gam*ry+rj)/sqrt(x)/(gam*ry1+rj1)
                SZ2 = VC(NP+2,k2,2)*(-plgndr(k2,1, cos(theta)))*sin(theta)*(gam*ry+rj)/sqrt(x)/(gam*ry1+rj1)
              else
              !potential solution
               SZ1 = VC(NP+2,k2,1)*(-plgndr(k2,1, cos(theta)))*sin(theta)/(x**(k2+1))
               SZ2 = VC(NP+2,k2,2)*(-plgndr(k2,1, cos(theta)))*sin(theta)/(x**(k2+1))
              endif
             else ! interior 
                SZ1 = x*VC(J,k2,1)*(-plgndr(k2,1, cos(theta)))*sin(theta)
                SZ2 = x*VC(J,k2,2)*(-plgndr(k2,1, cos(theta)))*sin(theta) 
             endif

!load array
                APR(J,k1)  = APR(J,k1)+SZ1
                API(J,k1)  = API(J,k1)+SZ2   

         ENDDO ! Plm
         ENDDO ! Radial   
         ENDDO ! theta

