write(bfeld3,2011) trim(dir)//'/tfld.',ii,'.t', jj, q, mm  
open(60,status='unknown',file=adjustl(bfeld3))
write(60,'(a)') 'Radius'
write(60,'(a)') 'theta'
write(60,'(1x,I4,2x,I4)') np+2+nf, n_theta
write(60,'(3x,f9.5,2x,f9.5,2x,f9.5,2x,f9.5)') x1, xf, theta1, theta2
do j=1, n_theta
  do I = 1, np+2+nf
  tc= jj*pi/nj
  write(60,'(e15.6)') bphi(i,j)*cos(abg(imag)*tc)+iphi(i,j)*sin(abg(imag)*tc)
  enddo
enddo
close(60)
