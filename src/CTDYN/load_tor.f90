theta1 = 0.
theta2 = 3.1415926
theta = theta1-theta2/(n_theta)

do k1=1, n_theta   ! theta-loop
  theta = theta + theta2/(n_theta)
  do j =1,np+2+nf    ! radial loop
    x = x1+(j-1)*hh
    ! load the b_n vectors for ff-solution or vacuum 
    do k2 = nst,net,2          
      if(beta.ne.0)then
        call bess(betb, k2, ffc)
        if (betb.gt.10) ffc=betb
        vc(np+2,k2,1) = (4*vc(np+1,k2,1)-vc(np,k2,1))/(2*hh*(k2+1+ffc)+3) 
        vc(np+2,k2,2) = (4*vc(np+1,k2,2)-vc(np,k2,2))/(2*hh*(k2+1+ffc)+3) 
      endif

      if (j.ge.(np+2)) then
        if (beta.ne.0) then
          xnu=k2*1.d0+1.0/2.0
          call bessjy(x*betb,xnu,rj,ry,rjp,ryp)
          call bessjy(betb,xnu,rj1,ry1,rjp1,ryp1)
          if(betb.gt.10) ffc=betb
          sz1 = vc(np+2,k2,1)*(-plgndr(k2,1, cos(theta)))*(gam*ry+rj)/sqrt(x)/(gam*ry1+rj1)
          sz2 = vc(np+2,k2,2)*(-plgndr(k2,1, cos(theta)))*(gam*ry+rj)/sqrt(x)/(gam*ry1+rj1)
        endif
      else
        sz1=0
        sz2=0 
      endif
      bphi(j,k1)  = bphi(j,k1)+sz1
      iphi(j,k1)  = iphi(j,k1)+sz2   

    enddo ! plm
  enddo ! radial   
enddo ! theta


