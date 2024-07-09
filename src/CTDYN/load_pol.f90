theta1 = 0.
theta2 = 3.1415926
theta = theta1-theta2/(n_theta)

do k1=1, n_theta   ! theta-loop
  theta = theta + theta2/(n_theta)
  do j =1,np+2+nf    ! radial loop
    x = x1+(j-1)*hh
    ! load the a_n vectors for dipoles: this defines the poloidal (potential)
    do k2 = nsp,nep,2          !plm loop 

    call bess(beta, k2, ffc)
    vc(np+2,k2,1)   = (4*vc(np+1,k2,1)-vc(np,k2,1))/(2*hh*(k2+1+ffc)+3)  
    vc(np+2,k2,2)   = (4*vc(np+1,k2,2)-vc(np,k2,2))/(2*hh*(k2+1+ffc)+3)    
    if(j.ge.(np+2))then
      if(beta.ne.0)then
        ! check the x infront of the definition as for the beta=0 case
        xnu=k2*1.d0+1.0/2.0
        call bessjy(x*beta,xnu,rj,ry,rjp,ryp)
        call bessjy(beta,xnu,rj1,ry1,rjp1,ryp1)
        sz1 = vc(np+2,k2,1)*(-plgndr(k2,1, cos(theta)))*sin(theta)*(gam*ry+rj)/sqrt(x)/(gam*ry1+rj1)
        sz2 = vc(np+2,k2,2)*(-plgndr(k2,1, cos(theta)))*sin(theta)*(gam*ry+rj)/sqrt(x)/(gam*ry1+rj1)
      else
        !potential solution
        sz1 = vc(np+2,k2,1)*(-plgndr(k2,1, cos(theta)))*sin(theta)/(x**(k2+1))
        sz2 = vc(np+2,k2,2)*(-plgndr(k2,1, cos(theta)))*sin(theta)/(x**(k2+1))
      endif
    else ! interior 
      sz1 = x*vc(j,k2,1)*(-plgndr(k2,1, cos(theta)))*sin(theta)
      sz2 = x*vc(j,k2,2)*(-plgndr(k2,1, cos(theta)))*sin(theta) 
    endif

    !load array
    apr(j,k1)  = apr(j,k1)+sz1
    api(j,k1)  = api(j,k1)+sz2   

    enddo ! plm
  enddo ! radial   
enddo ! theta

