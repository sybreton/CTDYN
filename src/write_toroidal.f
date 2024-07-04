

        write(bfeld3,2011) trim(dir)//'/tfld.',ii,'.t', jj, q, mm  
        open(60,status='unknown',file=adjustl(bfeld3))
        write(60,'(a)') 'Radius'
        write(60,'(a)') 'theta'
        write(60,'(1x,I4,2x,I4)') NP+2+NF,n_theta
        write(60,'(3x,f9.5,2x,f9.5,2x,f9.5,2x,f9.5)')x1,xf,theta1,theta2
        DO J =1,N_theta
          DO I = 1,NP+2+NF
          tc= jj*PI/nj
c          ts= sgn*jj*PI/nj
c          write(60,'(e15.6)')BPHI(I,J)*cos(abg(imag)*tc)+IPHI(I,J)*sin(abg(imag)*ts)
          write(60,'(e15.6)')BPHI(I,J)*cos(abg(imag)*tc)+IPHI(I,J)*sin(abg(imag)*tc)
          ENDDO
        ENDDO
        close(60)
C
