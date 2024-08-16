subroutine write_omega (np, n_theta, x1, x2, theta1, theta2, &
                        ome, sfu, ute)
  !-------------------------------------------------------------
  !
  !   write omega(r,theta)/omega_equator and stream function
  !
  !-------------------------------------------------------------
  integer :: np, n_theta
  real :: x1, x2, theta1, theta2
  real :: ome(np+2,n_theta)             ! omega
  real :: sfu(np+2,n_theta)             ! stream function 
  real :: ute(np+2,n_theta)             ! utheta

  integer :: i, j
  
  write(50,'(a)') 'Radius'
  write(50,'(a)') 'theta'
  write(50,'(1x,I4,2x,I4)') np+2,n_theta
  write(50,'(1x,f9.5,2x,f9.5,2x,f9.5,2x,f9.5)') x1, x2, theta1, theta2
  write(51,'(a)') 'Radius'
  write(51,'(a)') 'theta'
  write(51,'(1x,I4,2x,I4)') np+2, n_theta
  write(51,'(1x,f9.5,2x,f9.5,2x,f9.5,2x,f9.5)') x1, x2, theta1, theta2
  write(52,'(a)') 'Radius'
  write(52,'(a)') 'theta'
  write(52,'(1x,I4,2x,I4)')  np+2, n_theta
  write(52,'(1x,f9.5,2x,f9.5,2x,f9.5,2x,f9.5)') x1, x2, theta1, theta2
  
  do j=1, n_theta
    do i=1, np+2
      write(50,'(e15.6)') ome(i,j)
      write(51,'(e15.6)') sfu(i,j)
      write(52,'(e15.6)') ute(i,j)
    enddo
  enddo
  
  close(50)
  close(51)
  close(52)

end 
