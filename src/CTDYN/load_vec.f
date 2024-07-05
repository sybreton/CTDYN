
C load A-vectors

      DO   k2 = nsp,nep,2 
             call bess(beta, k2, ffc)
             VC(NP+2,k2,1) = (4*VC(NP+1,k2,1)-VC(NP,k2,1))/(2*hh*(k2+1+ffc)+3)  
             VC(NP+2,k2,2) = (4*VC(NP+1,k2,2)-VC(NP,k2,2))/(2*hh*(k2+1+ffc)+3)    
      DO J =2,NP+2+NF
      x = x1+(J-1)*hh
           aar = VC(J, k2, 1)
           aai = VC(J, k2, 2)
             if(J.ge.(NP+2))then
              if(beta .ne. 0)then

               xnu=1.0*k2+1/2.
               call bessjy(x*beta,xnu,rj,ry,rjp,ryp)
               call bessjy(beta,xnu,rj1,ry1,rjp1,ryp1)

               aar = VC(NP+2,k2,1)*(gam*ry+rj)/sqrt(x)/(gam*ry1+rj1)
               aai = VC(NP+2,k2,2)*(gam*ry+rj)/sqrt(x)/(gam*ry1+rj1)
              else 
! check exponent
               aar = VC(NP+2,k2,1)/x**(k2+1)
               aai = VC(NP+2,k2,2)/x**(k2+1)
              endif
             endif
            write(16,'(i2,e15.6,e15.6,e15.6)')  k2, x, aar,aai 
         ENDDO 
         ENDDO 

C load B-vectors

      DO   k2 = nst,net,2 

          if(beta .ne. 0)then 
           call bess(betb, k2, ffc)
           VC(NP+2,k2,1) = (4*VC(NP+1,k2,1)-VC(NP,k2,1))/(2*hh*(k2+1+ffc)+3)
           VC(NP+2,k2,2) = (4*VC(NP+1,k2,2)-VC(NP,k2,2))/(2*hh*(k2+1+ffc)+3)
          endif 

      DO J =2,NP+2+NF
      x = x1+(J-1)*hh
           bbr = VC(J, k2, 1)
           bbi = VC(J, k2, 2)
             if(J.ge.(NP+2))then
              if (beta .ne. 0)then

               xnu=1.0*k2+1/2.
               call bessjy(x*betb,xnu,rj,ry,rjp,ryp)
               call bessjy(betb,xnu,rj1,ry1,rjp1,ryp1)

               bbr = VC(NP+2,k2,1)*(gam*ry+rj)/sqrt(x)/(gam*ry1+rj1)
               bbi = VC(NP+2,k2,2)*(gam*ry+rj)/sqrt(x)/(gam*ry1+rj1)
             else 
               bbr = VC(NP+2,k2,1)/x**(k2+1)
               bbi = VC(NP+2,k2,2)/x**(k2+1)
             endif
            endif
            write(15,'(i2,e15.6,e15.6,e15.6)')  k2, x, bbr,bbi 
         ENDDO 
         ENDDO 

C surface energy

      etor = 0
    

      DO   k2 = nst,net,2 

          if(beta .ne. 0)then 
           call bess(betb, k2, ffc)
           VC(NP+2,k2,1) = (4*VC(NP+1,k2,1)-VC(NP,k2,1))/(2*hh*(k2+1+ffc)+3)
           VC(NP+2,k2,2) = (4*VC(NP+1,k2,2)-VC(NP,k2,2))/(2*hh*(k2+1+ffc)+3)
          endif 

           if(mm .ne. 0) stop

c          normt=(factrl(k2+mm+1)/factrl(k2-mm-1))/(2*k2+1.0)

         rosym = k2*(k2+1.0)/(2.0*k2+1.0)

        etor= etor+(VC(NP+2,k2,1)**2+VC(NP+2,k2,2)**2)*rosym

c        etor= etor+(VC(NP+1,k2,1)**2+VC(NP+1,k2,2)**2)*rosym

c         print*, k2, etor, VC(NP+1,k2,1), VC(NP+1,k2,2)  


         ENDDO 

c        print*, VC(NP+2,nst,1), nst
c        stop

C surface energy

      epol = 0
     
      DO   k2 = nsp,nep,2 
 
           call bess(beta, k2, ffc)
           VC(NP+2,k2,1) = (4*VC(NP+1,k2,1)-VC(NP,k2,1))/(2*hh*(k2+1+ffc)+3)
           VC(NP+2,k2,2) = (4*VC(NP+1,k2,2)-VC(NP,k2,2))/(2*hh*(k2+1+ffc)+3)

c          normt=(factrl(k2+mm+1)/factrl(k2-mm-1))/(2*k2+1.0)

         rosym = k2*(k2+1.0)/(2.0*k2+1.0)

         dd1= (VC(NP+2,k2,1)-VC(NP+1,k2,1))/hh
 
c aggiungi parte complessa...
 
        epol= epol+(VC(NP+2,k2,1)**2+VC(NP+2,k2,2)**2 +dd1**2+2*vc(np+2,k2,1)*dd1 )*rosym
     .        +(VC(NP+2,k2,1)**2+VC(NP+2,k2,2)**2)*k2*(k2+1)*rosym

         ENDDO 

        etep= etor/epol
        etet=etor/(epol+etor)

        print*, 'Et/Ep', etep
        print*, 'Et/Etot', etet

      
c          print*, etor, epol


c        write(18,*) etor/epol, etor/(epol+etor)



