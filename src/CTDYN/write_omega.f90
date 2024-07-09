
!-------------------------------------------------------------
!
!   write omega(r,theta)/omega_equator  and stream function
!
!-------------------------------------------------------------

write(50,'(a)') 'Radius'
write(50,'(a)') 'theta'
write(50,'(1x,I4,2x,I4)') NP+2,n_theta
write(50,'(1x,f9.5,2x,f9.5,2x,f9.5,2x,f9.5)')x1,x2,theta1,theta2
write(51,'(a)') 'Radius'
write(51,'(a)') 'theta'
write(51,'(1x,I4,2x,I4)') NP+2,n_theta
write(51,'(1x,f9.5,2x,f9.5,2x,f9.5,2x,f9.5)')x1,x2,theta1,theta2
write(52,'(a)') 'Radius'
write(52,'(a)') 'theta'
write(52,'(1x,I4,2x,I4)') NP+2,n_theta
write(52,'(1x,f9.5,2x,f9.5,2x,f9.5,2x,f9.5)')x1,x2,theta1,theta2

do J=1, N_theta
  do I=1, NP+2
    write(50,'(e15.6)') OME(I,J)
    write(51,'(e15.6)') SFU(I,J)
    write(52,'(e15.6)') UTE(I,J)
  enddo
enddo

close(50)
close(51)
close(52)

