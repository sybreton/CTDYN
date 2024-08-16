! load a-vectors
do k2 = nsp,nep,2 
  call bess(beta, k2, ffc)
  vc(np+2,k2,1) = (4*vc(np+1,k2,1)-vc(np,k2,1))/(2*hh*(k2+1+ffc)+3)  
  vc(np+2,k2,2) = (4*vc(np+1,k2,2)-vc(np,k2,2))/(2*hh*(k2+1+ffc)+3)    
  do j =2,np+2+nf
    x = x1+(j-1)*hh
    aar = vc(j, k2, 1)
    aai = vc(j, k2, 2)
    if (j.ge.(np+2)) then
      if (beta .ne. 0) then
        xnu=1.0*k2+1/2.
        call bessjy(x*beta,xnu,rj,ry,rjp,ryp)
        call bessjy(beta,xnu,rj1,ry1,rjp1,ryp1)
        aar = vc(np+2,k2,1)*(gam*ry+rj)/sqrt(x)/(gam*ry1+rj1)
        aai = vc(np+2,k2,2)*(gam*ry+rj)/sqrt(x)/(gam*ry1+rj1)
      else 
        ! check exponent
        aar = vc(np+2,k2,1)/x**(k2+1)
        aai = vc(np+2,k2,2)/x**(k2+1)
      endif
    endif
  enddo 
enddo 

! load b-vectors
do k2 = nst,net,2 
  if (beta .ne. 0) then 
     call bess(betb, k2, ffc)
     vc(np+2,k2,1) = (4*vc(np+1,k2,1)-vc(np,k2,1))/(2*hh*(k2+1+ffc)+3)
     vc(np+2,k2,2) = (4*vc(np+1,k2,2)-vc(np,k2,2))/(2*hh*(k2+1+ffc)+3)
  endif 
  do j=2, np+2+nf
    x = x1+(j-1)*hh
    bbr = vc(j, k2, 1)
    bbi = vc(j, k2, 2)
    if (j.ge.(np+2)) then
      if (beta .ne. 0) then
        xnu=1.0*k2+1/2.
        call bessjy(x*betb,xnu,rj,ry,rjp,ryp)
        call bessjy(betb,xnu,rj1,ry1,rjp1,ryp1)

        bbr = vc(np+2,k2,1)*(gam*ry+rj)/sqrt(x)/(gam*ry1+rj1)
        bbi = vc(np+2,k2,2)*(gam*ry+rj)/sqrt(x)/(gam*ry1+rj1)
      else 
        bbr = vc(np+2,k2,1)/x**(k2+1)
        bbi = vc(np+2,k2,2)/x**(k2+1)
      endif
    endif
  enddo 
enddo 

! surface energy
etor = 0
do k2 = nst,net,2 
  if (beta .ne. 0) then 
    call bess(betb, k2, ffc)
    vc(np+2,k2,1) = (4*vc(np+1,k2,1)-vc(np,k2,1))/(2*hh*(k2+1+ffc)+3)
    vc(np+2,k2,2) = (4*vc(np+1,k2,2)-vc(np,k2,2))/(2*hh*(k2+1+ffc)+3)
  endif 
  if (mm .ne. 0) stop
    rosym = k2*(k2+1.0)/(2.0*k2+1.0)
    etor= etor+(vc(np+2,k2,1)**2+vc(np+2,k2,2)**2)*rosym
enddo 

! surface energy
epol = 0
do k2=nsp, nep, 2 
  call bess(beta, k2, ffc)
  vc(np+2,k2,1) = (4*vc(np+1,k2,1)-vc(np,k2,1))/(2*hh*(k2+1+ffc)+3)
  vc(np+2,k2,2) = (4*vc(np+1,k2,2)-vc(np,k2,2))/(2*hh*(k2+1+ffc)+3)
  rosym = k2*(k2+1.0)/(2.0*k2+1.0)
  dd1= (vc(np+2,k2,1)-vc(np+1,k2,1))/hh
  epol= epol+(vc(np+2,k2,1)**2+vc(np+2,k2,2)**2 +dd1**2+2*vc(np+2,k2,1)*dd1 )*rosym &
        & +(vc(np+2,k2,1)**2+vc(np+2,k2,2)**2)*k2*(k2+1)*rosym
enddo 

etep= etor/epol
etet=etor/(epol+etor)

print*, 'et/ep', etep
print*, 'et/etot', etet

