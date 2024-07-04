
        write(bfeld4,2011) trim(dir)//'/pfld.',ii,'.t',jj, q, mm  
        open(61,status='unknown',file=adjustl(bfeld4))
        write(61,'(a)') 'Radius'
        write(61,'(a)') 'theta'
        write(61,'(1x,I4,2x,I4)')NP+2+NF,n_theta
        write(61,'(3x,f9.5,2x,f9.5,2x,f9.5,2x,f9.5)')x1,xf,theta1,theta2

        DO J =1,N_theta
        theta1 = 0.
        theta2 = 3.1415926
        theta = theta1-theta2/(n_theta)

          DO I = 1,NP+2+NF
          tc= jj*PI/nj
c          ts= sgn*jj*PI/nj
c          write(61,'(e15.6)')(APR(I,J)*cos(tc)+API(I,J)*sin(ts))*x*sin(theta)
c          write(61,'(e15.6)')(APR(I,J)*cos(abg(imag)*tc)+API(I,J)*sin(abg(imag)*ts))*x*sin(theta)
          write(61,'(e15.6)')(APR(I,J)*cos(abg(imag)*tc)+API(I,J)*sin(abg(imag)*tc))*x*sin(theta)
          ENDDO
        ENDDO

        close(61)
