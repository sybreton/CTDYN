
      theta1 = 0.
      theta2 = 3.1415926
      theta = theta1-theta2/(n_theta)

      DO K1=1, n_theta   ! theta-loop
        theta = theta + theta2/(n_theta)
      DO J =1,NP+2+NF    ! radial loop
      x = x1+(J-1)*hh
               ! LOAD THE b_n vectors for FF-solution or vacuum 
      DO   k2 = nst,net,2          

      if(beta.ne.0)then
                 call bess(betb, k2, ffc)

          if(betb.gt.10) ffc=betb

                 VC(NP+2,k2,1) = (4*VC(NP+1,k2,1)-VC(NP,k2,1))/(2*hh*(k2+1+ffc)+3) 
                 VC(NP+2,k2,2) = (4*VC(NP+1,k2,2)-VC(NP,k2,2))/(2*hh*(k2+1+ffc)+3) 
      endif

                  if(J.ge.(NP+2))then
                   if(beta.ne.0)then
                    xnu=k2*1.d0+1.0/2.0
                    call bessjy(x*betb,xnu,rj,ry,rjp,ryp)
                    call bessjy(betb,xnu,rj1,ry1,rjp1,ryp1)

          if(betb.gt.10) ffc=betb

                    SZ1 = VC(NP+2,k2,1)*(-plgndr(k2,1, cos(theta)))*(gam*ry+rj)/sqrt(x)/(gam*ry1+rj1)
                    SZ2 = VC(NP+2,k2,2)*(-plgndr(k2,1, cos(theta)))*(gam*ry+rj)/sqrt(x)/(gam*ry1+rj1)
                   endif
                  else
                   SZ1=0
                   SZ2=0 
                  endif

                 BPHI(J,k1)  = BPHI(J,k1)+SZ1
                 IPHI(J,k1)  = IPHI(J,k1)+SZ2   

         ENDDO ! Plm
         ENDDO ! Radial   
         ENDDO ! theta


