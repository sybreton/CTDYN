write(bfeld4,2011) trim(dir)//'/pfld.',ii,'.t',jj, q, mm  
open(61,status='unknown',file=adjustl(bfeld4))
write(61,'(a)') 'Radius'
write(61,'(a)') 'theta'
write(61,'(1x,I4,2x,I4)')NP+2+NF,n_theta
write(61,'(3x,f9.5,2x,f9.5,2x,f9.5,2x,f9.5)')x1,xf,theta1,theta2

do j=1, n_theta
  theta1 = 0.
  theta2 = 3.1415926
  theta = theta1-theta2/(n_theta)
  do I = 1,np+2+nf
    tc= jj*pi/nj
    write(61,'(e15.6)')(apr(i,j)*cos(abg(imag)*tc)+api(i,j)*sin(abg(imag)*tc))*x*sin(theta)
  enddo
enddo

close(61)
