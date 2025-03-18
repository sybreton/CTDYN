! 
!  IMPORTANT: the eigenfunctions are always complex conjugate variable. If the eigenvalue has a positive Im part, 
!  be sure you keep the correct conjugate eigenvector
! 
!
!  IMPORTANT:
!
! The conventions for the vector potential and scalar potentials are in the Maple ws called currenthelicity.mw (and similar)
! and on the mathematica notebook notation.nb
!
! In spherical symmetry in order to make contact with the standard representation used in papers:
!
! Br= (1/r sintheta) \partial_\theta [A(r,theta) sintheta]
! 
! Btheta =(-1/r) \partial_r [A r]
!
! Bphi = B
!
! you must remember that, according to the Radler formalism of the Psi and Phi functions we have:
!    
!      B = -\partial_\theta P0n Psi_n = -P1n Psi_n
! and
!      A = -\partial_\theta P0n Phi_n = -P1n Phi_n
!
!    where we used the identity: \partial P0n = P1n to  make contact with the notations in the papers 
!    (note that the above identity differs by a sign from what you find in Rudiger's book in the appendix)
!
!     
module write_outputs 

  use kind_parameters
  use cio
  use cdfunc
  use bessel
  use util
  use func_flow
  use stellar_profiles
  use plg

  implicit none 
  private 

  public :: writefield 

  character(*), parameter :: fmt_2010 = "(a220,i6.6,a2,i2.2)" 
  character(*), parameter :: fmt_2011 = "(a220,i6.6,a2,i2.2,a2,i2.2)" 
  integer, parameter :: nft=200   
  integer :: mm 

contains 

  subroutine compute_pol (nsp, nep, nf, nft, vc, apr, api, &
                          ffc, hh, x1)
    ! -----------------------------------------------------------------------    
    ! > Compute the external extrapolation for the poloidal field.   
    ! -----------------------------------------------------------------------    

    ! arguments
    integer :: nsp, nep, nf, nft
    real(dp) :: vc(np+2,nb,2)                 ! eigenvector
    real(dp) :: apr(np+2+nft, n_theta)         ! b poloidal (potential)
    real(dp) :: api(np+2+nft, n_theta)         ! b poloidal (potential)
    real(dp) :: ffc, hh, x1

    ! local variables
    real(dp) :: x, xnu, theta, sz1, sz2
    integer :: k1, k2, j
    real(dp) :: rj, ry, rjp, ryp
    real(dp) :: rj1, ry1, rjp1, ryp1

    do k1=1, n_theta   ! theta-loop
      theta = theta1 + k1*dt
      do j=1, np+2+nf    ! radial loop
        x = x1+(j-1)*hh
        do k2=nsp, nep, 2          !plm loop 

          call bess(beta, k2, ffc)
          vc(np+2,k2,1) = (4*vc(np+1,k2,1)-vc(np,k2,1))/(2*hh*(k2+1+ffc)+3)
          vc(np+2,k2,2) = (4*vc(np+1,k2,2)-vc(np,k2,2))/(2*hh*(k2+1+ffc)+3)
          if (j.ge.(np+2)) then
            if (beta.ne.0) then
              ! check the x infront of the definition as for the beta=0 case
              xnu=k2*1.d0+1.0d0/2.0d0
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
          apr(j,k1) = apr(j,k1)+sz1
          api(j,k1) = api(j,k1)+sz2

        enddo ! plm
      enddo ! radial   
    enddo ! theta

  end subroutine

  subroutine compute_tor (nst, net, nf, nft, vc, bphi, iphi, ffc, &
                          hh, x1)
    ! -----------------------------------------------------------------
    ! > Compute the external extrapolation for the toroidal field.
    ! -----------------------------------------------------------------
  
    ! arguments
    integer :: nst, net, nf, nft
    real(dp) :: vc(np+2,nb,2)                 ! eigenvector
    real(dp) :: bphi(np+2+nft, n_theta)         ! b poloidal (potential)
    real(dp) :: iphi(np+2+nft, n_theta)         ! b poloidal (potential)
    real(dp) :: ffc, hh, x1
  
    ! local variables
    real(dp) :: x, xnu, theta, sz1, sz2
    integer :: k1, k2, j
    real(dp) :: rj, ry, rjp, ryp
    real(dp) :: rj1, ry1, rjp1, ryp1
  
    do k1=1, n_theta   ! theta-loop
      theta = theta1 + k1*dt
      do j=1, np+2+nf    ! radial loop
        x = x1+(j-1)*hh
        ! load the b_n vectors for ff-solution or vacuum 
        do k2=nst, net, 2
          if(beta.ne.0)then
            call bess(betb, k2, ffc)
            if (betb.gt.10d0) ffc=betb
            vc(np+2,k2,1) = (4d0*vc(np+1,k2,1)-vc(np,k2,1))/(2d0*hh*(k2+1d0+ffc)+3d0)
            vc(np+2,k2,2) = (4d0*vc(np+1,k2,2)-vc(np,k2,2))/(2d0*hh*(k2+1d0+ffc)+3d0)
          endif
  
          if (j.ge.(np+2)) then
            if (beta.ne.0) then
              xnu=k2*1.d0+1.0d0/2.0d0
              call bessjy(x*betb,xnu,rj,ry,rjp,ryp)
              call bessjy(betb,xnu,rj1,ry1,rjp1,ryp1)
              if(betb.gt.10d0) ffc=betb
              sz1 = vc(np+2,k2,1)*(-plgndr(k2,1, cos(theta)))*(gam*ry+rj)/sqrt(x)/(gam*ry1+rj1)
              sz2 = vc(np+2,k2,2)*(-plgndr(k2,1, cos(theta)))*(gam*ry+rj)/sqrt(x)/(gam*ry1+rj1)
            endif
          else
            sz1=0d0
            sz2=0d0
          endif
          bphi(j,k1) = bphi(j,k1)+sz1
          iphi(j,k1) = iphi(j,k1)+sz2
  
        enddo ! plm
      enddo ! radial   
    enddo ! theta
  
  end subroutine

  subroutine compute_ab_vector (vc, aai, aar, bbi, bbr, &
                                ffc, x1, hh, nsp, nep, &
                                nst, net, nf)
    ! -----------------------------------------------------------------
    ! Compute A and B vectors and surface energy
  
    ! The presence of aai and aar would make sense only if they were
    ! vectors with the same size as vc.
    ! -----------------------------------------------------------------
  
    ! arguments
    real(dp) :: vc(np+2,nb,2)                 ! eigenvector
    real(dp) :: aai, aar, bbi, bbr
    real(dp) :: ffc, x1, hh
    integer :: nsp, nep, nst, net, nf, mm
  
    ! local variables
    real(dp) :: etor, epol
    real(dp) :: x, xnu, theta, rosym, dd1
    integer :: k2, j
    real(dp) :: rj, ry, rjp, ryp
    real(dp) :: rj1, ry1, rjp1, ryp1
  
    ! load a-vectors
    do k2=nsp, nep, 2
      call bess(beta, k2, ffc)
      vc(np+2,k2,1) = (4d0*vc(np+1,k2,1)-vc(np,k2,1))/(2d0*hh*(k2+1d0+ffc)+3d0)
      vc(np+2,k2,2) = (4d0*vc(np+1,k2,2)-vc(np,k2,2))/(2d0*hh*(k2+1d0+ffc)+3d0)
      do j=2, np+2+nf
        if (j.lt.(np+2)) then
          aar = vc(j, k2, 1)
          aai = vc(j, k2, 2)
        else if (j.ge.(np+2)) then
          x = x1+(j-1)*hh
          if (beta .ne. 0) then
            xnu=1.0d0*k2+1d0/2.d0
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
         vc(np+2,k2,1) = (4d0*vc(np+1,k2,1)-vc(np,k2,1))/(2d0*hh*(k2+1d0+ffc)+3d0)
         vc(np+2,k2,2) = (4d0*vc(np+1,k2,2)-vc(np,k2,2))/(2d0*hh*(k2+1d0+ffc)+3d0)
      endif
      do j=2, np+2+nf
        if (j.lt.(np+2)) then
          bbr = vc(j, k2, 1)
          bbi = vc(j, k2, 2)
        else if (j.ge.(np+2)) then
          x = x1+(j-1)*hh
          if (beta .ne. 0) then
            xnu=1.0d0*k2+0.5d0
            call bessjy(x*betb, xnu,rj, ry, rjp, ryp)
            call bessjy(betb, xnu, rj1, ry1, rjp1, ryp1)
  
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
        vc(np+2,k2,1) = (4d0*vc(np+1,k2,1)-vc(np,k2,1))/(2d0*hh*(k2+1d0+ffc)+3d0)
        vc(np+2,k2,2) = (4d0*vc(np+1,k2,2)-vc(np,k2,2))/(2d0*hh*(k2+1d0+ffc)+3d0)
      endif
      if (mm .ne. 0) stop
      rosym = k2*(k2+1.0d0)/(2.0d0*k2+1.0d0)
      etor = etor+(vc(np+2,k2,1)**2+vc(np+2,k2,2)**2)*rosym
    enddo
  
    ! surface energy
    epol = 0
    do k2=nsp, nep, 2
      call bess(beta, k2, ffc)
      vc(np+2,k2,1) = (4d0*vc(np+1,k2,1)-vc(np,k2,1))/(2d0*hh*(k2+1d0+ffc)+3d0)
      vc(np+2,k2,2) = (4d0*vc(np+1,k2,2)-vc(np,k2,2))/(2d0*hh*(k2+1d0+ffc)+3d0)
      rosym = k2*(k2+1.0d0)/(2.0d0*k2+1.0d0)
      dd1 = (vc(np+2,k2,1)-vc(np+1,k2,1))/hh
      epol = epol+(vc(np+2,k2,1)**2+vc(np+2,k2,2)**2 +dd1**2+2d0*vc(np+2,k2,1)*dd1 )*rosym &
            & +(vc(np+2,k2,1)**2+vc(np+2,k2,2)**2)*k2*(k2+1)*rosym
    enddo
  
    etep = etor/epol
    etet = etor/(epol+etor)
  
    write (*,*) 'Et/Ep', etep
    write (*,*) 'Et/Etot', etet
  
  end subroutine

  subroutine write_toroidal (bfeld3, imag, nf, nj, jj, q, &
                             theta, x1, xf, tc, bphi, iphi)

    ! -----------------------------------------------------------------
    !> Write toroidal field to output file. 
    ! -----------------------------------------------------------------
  
    character(len=512) :: bfeld3
    character(len=2) :: q
    integer :: nf, nj, jj
    real(dp) :: imag, theta, x1, xf, tc
    real(dp) :: bphi(np+2+nft, n_theta)        ! toroidal real
    real(dp) :: iphi(np+2+nft, n_theta)        ! toroidal im
  
    integer :: i, j
  
    write(bfeld3, '(a220,i6.6,a2,i2.2,a2,i2.2)') trim(dir)//'/tfld.',ii,'.t', jj, q, mm
    open(60, status='unknown', file=adjustl(bfeld3))
    write(60, '(a)') 'Radius'
    write(60, '(a)') 'theta'
    write(60, '(1x,I4,2x,I4)') np+2+nf, n_theta
    write(60, '(3x,f9.5,2x,f9.5,2x,f9.5,2x,f9.5)') x1, xf, theta1, theta2
    do j=1, n_theta
      do i=1, np+2+nf
        tc= jj*pi/nj
        write(60,'(e15.6)') bphi(i,j)*cos(abg(imag)*tc)+iphi(i,j)*sin(abg(imag)*tc)
      enddo
    enddo
    close(60)
  end subroutine
  
  subroutine write_poloidal (bfeld4, imag, nf, nj, jj, q, &
                             theta, hh, x1, xf, tc, apr, api) 
 
    ! -----------------------------------------------------------------
    !> Write poloidal field to output file. 
    ! -----------------------------------------------------------------
  
    ! Arguments
    character(len=512) :: bfeld4
    character(len=2) :: q
    integer :: nf, nj, jj
    real(dp) :: imag, theta, hh, x1, xf, tc
    real(dp) :: apr(np+2+nft, n_theta)         ! b poloidal (potential)
    real(dp) :: api(np+2+nft, n_theta)         ! b poloidal (potential)
  
    ! Local variables
    integer :: i, j
    real(dp) :: x
  
    write(bfeld4, fmt_2011) trim(dir)//'/pfld.', ii, '.t', jj, q, mm
    open(61,status='unknown', file=adjustl(bfeld4))
    write(61,'(a)') 'Radius'
    write(61,'(a)') 'theta'
    write(61,'(1x,I4,2x,I4)') np+2+nf, n_theta
    write(61,'(3x,f9.5,2x,f9.5,2x,f9.5,2x,f9.5)') x1, xf, theta1, theta2
  
    do j=1, n_theta
      theta = theta1 + j*dt
      do i=1, np+2+nf
        tc = jj*pi/nj
        x = x1+(i-1d0)*hh
        write(61,'(e15.6)') (apr(i,j)*cos(abg(imag)*tc)+api(i,j)*sin(abg(imag)*tc))*x*sin(theta)
      enddo
    enddo
  
    close(61)
  end subroutine
  
  subroutine write_omega (np, n_theta, x1, x2, theta1, theta2, &
                          ome, sfu, ute)
    !-------------------------------------------------------------
    !
    ! > Write omega(r,theta)/omega_equator and stream function
    !
    !-------------------------------------------------------------
  
    integer :: np, n_theta
    real(dp) :: x1, x2, theta1, theta2
    real(dp) :: ome(np+2,n_theta)             ! omega
    real(dp) :: sfu(np+2,n_theta)             ! stream function 
    real(dp) :: ute(np+2,n_theta)             ! utheta
  
    integer :: i, j

    open(50,status='unknown',file=trim(dir)//'/omega.dat')
    open(51,status='unknown',file=trim(dir)//'/stream.dat')
    open(52,status='unknown',file=trim(dir)//'/utheta.dat')
  
    write(50,'(a)') 'Radius'
    write(50,'(a)') 'theta'
    write(50,'(1x,I4,2x,I4)') np+2, n_theta
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
  end subroutine 

  subroutine write_vect (bfeld2, q) 
    !-------------------------------------------------------------
    ! TODO
    !-------------------------------------------------------------

    character(len=512) :: bfeld2
    character(len=2) :: q

    write(bfeld2, fmt_2010) trim(dir)//'/vect.', ii, q, mm 
    open(13, status='unknown', file=adjustl(bfeld2))
    close(13)

  end subroutine


  subroutine time_butterfly (bfeld1, bfeld6, bfeld7, bfeld8, &
                             imag, nf, nft, x1, x2, hh, bphi, iphi, &
                             brr, bri, apr, api, q) 
    !-------------------------------------------------------------
    !  > Butterfly diagram block      
    !  Time evolution butterfly diagram 
    !-------------------------------------------------------------

    character(len=512) :: bfeld1, bfeld6, bfeld7, bfeld8
    character(len=2) :: q
    integer :: nf, nft
    real(dp) :: imag, x1, x2, hh
    real(dp) :: bphi(np+2+nft,n_theta)        ! toroidal real
    real(dp) :: iphi(np+2+nft,n_theta)        ! toroidal im
    real(dp) :: brr(np+2+nft,n_theta)         ! radial real
    real(dp) :: bri(np+2+nft,n_theta)         ! radial im
    real(dp) :: apr(np+2+nft,n_theta)         ! b poloidal (potential)
    real(dp) :: api(np+2+nft,n_theta)         ! b poloidal (potential)


    ! local arguments
    real(dp), parameter :: t_fin = 12.d0
    real(dp) :: time, ratio, theta
    real(dp) :: x
    real(dp) :: xr(np+2+nft)                
    integer :: bin(2), bax(2), rin(2), rax(2)
    integer :: k1, k2, j, jo, jbt, jbo
    real(dp) :: bphi7, brs, apt, bpt, bptkm, aptkm, &
            bptkp, aptkp, bptjm, aptjm, bptjp, &
            aptjp, dbtheta, datheta, dbtheta2, &
            datheta2, dax, dbx, dax2, chel, chel2


    write(bfeld1, fmt_2010) trim(dir)//'/butf.',ii,q,mm 
    open(18, status='unknown',file=adjustl(bfeld1))      
    write(bfeld6, fmt_2010) trim(dir)//'/jcin.',ii,q,mm 
    open(30, status='unknown',file=adjustl(bfeld6))      
    write(bfeld7, fmt_2010) trim(dir)//'/jcou.',ii,q,mm 
    open(31, status='unknown',file=adjustl(bfeld7))      
    write(bfeld8, fmt_2010) trim(dir)//'/brbp.',ii,q,mm 
    open(32, status='unknown',file=adjustl(bfeld8))      
  
    write(18,'(6x,i4)') n_theta-2
    write(30,'(6x,i4)') n_theta-2
    write(31,'(6x,i4)') n_theta-2
    write(32,'(6x,i4)') n_theta-2
  
    ! location of max bphi
    bax=maxloc(sqrt(bphi**2+iphi**2))
    bin=minloc(sqrt(bphi**2+iphi**2))
    rax=maxloc(sqrt(apr**2+api**2))
    rin=minloc(sqrt(apr**2+api**2))
  
    ratio= (sqrt(bphi(bax(1), bax(2))**2 + iphi(bax(1), bax(2))**2) &
           & /sqrt(apr(rax(1), rax(2))**2 + api(rax(1), rax(2))**2) )**(-1)
  
    write(*,'(a,f10.4)') 'max poloidal / max toroidal  ',  ratio
    write(*,'(a,f10.4)') 'Bphi max location', bax(1)*(x2-x1)/(np+1) + x1
    write(*,'(a,f10.4)') 'a max location', rax(1)*(x2-x1)/(np+1) + x1
  
    do k2=2, n_theta-1
      theta = theta1 + k2*dt
      write(18,'(e12.5)') theta
      write(30,'(e12.5)') theta
      write(31,'(e12.5)') theta
      write(32,'(e12.5)') theta
    enddo
  
    do j=1, np+2+nf        ! radial loop
      x = x1+(j-1)*hh
      xr(j) = x
    enddo
  
    jbt = hunt (xr, np+2+nf, xbt)
    jbo = hunt (xr, np+2+nf, xbo)
  
    time = -t_fin/n_time
    do k1=1, n_time
      time = time + t_fin/n_time
      write(18,'(e12.5)') time
      write(30,'(e12.5)') time
      write(31,'(e12.5)') time
      write(32,'(e12.5)') time
      do k2=2, n_theta-1
        theta = theta1 + k2*dt
        bphi7 = bphi(jbt,k2)*cos(abg(imag)*time)+iphi(jbt,k2)*sin(abg(imag)*time)
        j = jbt
        brs  = brr(np,k2)*cos(abg(imag)*time)+bri(np,k2)*sin(abg(imag)*time)
        apt    = apr(j,k2)*cos(abg(imag)*time)+api(j,k2)*sin(abg(imag)*time)
        bpt    = bphi(j,k2)*cos(abg(imag)*time)+iphi(j,k2)*sin(abg(imag)*time)
        bptkm  = bphi(j,k2-1)*cos(abg(imag)*time)+iphi(j,k2-1)*sin(abg(imag)*time)
        aptkm  = apr(j,k2-1)*cos(abg(imag)*time)+api(j,k2-1)*sin(abg(imag)*time)
        bptkp  = bphi(j,k2+1)*cos(abg(imag)*time)+iphi(j,k2+1)*sin(abg(imag)*time)
        aptkp  = apr(j,k2+1)*cos(abg(imag)*time)+api(j,k2+1)*sin(abg(imag)*time)
        bptjm  = bphi(j-1,k2)*cos(abg(imag)*time)+iphi(j-1,k2)*sin(abg(imag)*time)
        aptjm  = apr(j-1,k2)*cos(abg(imag)*time)+api(j-1,k2)*sin(abg(imag)*time)
        bptjp  = bphi(j+1,k2)*cos(abg(imag)*time)+iphi(j+1,k2)*sin(abg(imag)*time)
        aptjp  = apr(j+1,k2)*cos(abg(imag)*time)+api(j+1,k2)*sin(abg(imag)*time)
        ! first derivative wrt theta
        dbtheta =(bptkp-bptkm)/dt/2 
        datheta =(aptkp-aptkm)/dt/2
        ! second derivative wrt theta
        datheta2 = (aptkp-2*apt+aptkm)/dt/dt
        dbtheta2 = (bptkp-2*bpt+bptkm)/dt/dt
        ! first derivative wrt x
        dax= (apt-aptjm)/hh
        dbx= (bpt-bptjm)/hh
        ! second derivative wrt x
        dax2=(aptjp-2*apt+aptjm)/hh/hh
        ! to understand this formula for the cross helicity, just open the mathematica notebook chelicity2.nb
        chel=(2*apt*bpt/sin(theta)**2 + apt*dbtheta*cos(theta)/sin(theta) & 
         & +datheta*dbtheta -bpt*datheta2)/xr(j)/xr(j) &
         & +(-bpt*dax+apt*dbx)/xr(j)+dax*dbx -bpt*dax2
  
        jo = jbo
        apt    = apr(jo,k2)*cos(abg(imag)*time)+api(jo,k2)*sin(abg(imag)*time)
        bpt    = bphi(jo,k2)*cos(abg(imag)*time)+iphi(jo,k2)*sin(abg(imag)*time)
        bptkm  = bphi(jo,k2-1)*cos(abg(imag)*time)+iphi(jo,k2-1)*sin(abg(imag)*time)
        aptkm  = apr(jo,k2-1)*cos(abg(imag)*time)+api(jo,k2-1)*sin(abg(imag)*time)
        bptkp  = bphi(jo,k2+1)*cos(abg(imag)*time)+iphi(jo,k2+1)*sin(abg(imag)*time)
        aptkp  = apr(jo,k2+1)*cos(abg(imag)*time)+api(jo,k2+1)*sin(abg(imag)*time)
        bptjm  = bphi(jo-1,k2)*cos(abg(imag)*time)+iphi(jo-1,k2)*sin(abg(imag)*time)
        aptjm  = apr(jo-1,k2)*cos(abg(imag)*time)+api(jo-1,k2)*sin(abg(imag)*time)
        bptjp  = bphi(jo+1,k2)*cos(abg(imag)*time)+iphi(jo+1,k2)*sin(abg(imag)*time)
        aptjp  = apr(jo+1,k2)*cos(abg(imag)*time)+api(jo+1,k2)*sin(abg(imag)*time)
        ! first derivative wrt theta
        dbtheta =(bptkp-bptkm)/dt/2 
        datheta =(aptkp-aptkm)/dt/2
        ! second derivative wrt theta
        datheta2 = (aptkp-2*apt+aptkm)/dt/dt
        dbtheta2 = (bptkp-2*bpt+bptkm)/dt/dt
        ! first derivative wrt x
        dax= (apt-aptjm)/hh
        dbx= (bpt-bptjm)/hh
        ! second derivative wrt x
        dax2=(aptjp-2*apt+aptjm)/hh/hh
        ! to understand this formula for the cross helicity, just open the mathematica notebook chelicity2.nb
        chel2=(2*apt*bpt/sin(theta)**2 + apt*dbtheta*cos(theta)/sin(theta) & 
         & +datheta*dbtheta -bpt*datheta2)/xr(jo)/xr(jo) &
         & +(-bpt*dax+apt*dbx)/xr(jo)+dax*dbx -bpt*dax2

        write(18,'(e12.5)')  bphi7
        write(30,'(e12.5)') chel
        write(31,'(e12.5)')  chel2
        write(32,'(e12.5)') brs !* bphi7
      enddo
    enddo
  
    close(18)
    close(30)
    close(31)
    close(32)

  end subroutine 

  subroutine radial_butterfly (bfeld9, bfeld10, imag, nf, nft, x1, x2, &
                               hh, bphi, iphi, brr, bri, apr, api, q) 

    !----------------------------------------------------------
    ! > Radial butterfly diagram n and s 
    !----------------------------------------------------------

    character(len=512) :: bfeld9, bfeld10
    character(len=2) :: q
    integer :: nf, nft
    real(dp) :: imag, x1, x2, hh
    real(dp) :: bphi(np+2+nft,n_theta)        ! toroidal real
    real(dp) :: iphi(np+2+nft,n_theta)        ! toroidal im
    real(dp) :: brr(np+2+nft,n_theta)         ! radial real
    real(dp) :: bri(np+2+nft,n_theta)         ! radial im
    real(dp) :: apr(np+2+nft,n_theta)         ! b poloidal (potential)
    real(dp) :: api(np+2+nft,n_theta)         ! b poloidal (potential)


    ! local arguments
    real(dp), parameter :: t_fin = 12.d0
    real(dp) :: ang(n_theta)
    real(dp) :: time, theta, zq1, zq2
    real(dp) :: x, rnor, cja
    real(dp) :: xr(np+2+nft)                
    integer :: bin(2), bax(2), rin(2), rax(2)
    integer :: k1, k2, kb, km, j, jo, ko, jbt, jbo
    real(dp) :: bphi7, brs, apt, bpt, bptkm, aptkm, &
            bptkp, aptkp, bptjm, aptjm, aptjm2, bptjp, &
            aptjp, dbtheta, datheta, dbtheta2, &
            datheta2, dax, dbx, dax2, chels, cheln
    integer, parameter  :: nag = 2
    real(dp), parameter :: abt = 1.047d0
    real(dp), parameter :: abm = 2.094d0
    
    bax=maxloc(abs(apr))
    bin=minloc(abs(apr))
    zq1= bphi(bax(1), bax(2))
    bax=maxloc(abs(api))
    bin=minloc(abs(api))
    zq2=iphi(bax(1), bax(2))
    rnor=((zq1+zq2)/2.d0)**2
  
    write(bfeld9, fmt_2010) trim(dir)//'/brtn.',ii,q,mm 
    open(33,status='unknown',file=adjustl(bfeld9))      
    write(bfeld10, fmt_2010) trim(dir)//'/brts.',ii,q,mm 
    open(34,status='unknown',file=adjustl(bfeld10))      
  
    write(33,'(6x,i4)') np+2+nf-2
    write(34,'(6x,i4)') np+2+nf-2
  
    do j=2, np+2+nf-1       ! radial loop
      x = x1+(j-1d0)*hh
      write(33,'(e12.5)') x
      write(34,'(e12.5)') x
    enddo
  
    do k2=1, n_theta
      theta = theta1 + k2*dt
      ang(k2) = theta
    enddo
  
    kb = hunt (ang, n_theta, abt)
    km = hunt (ang, n_theta, abm)
  
    time = -t_fin/n_time
    do k1=1, n_time
      time = time + t_fin/n_time
      write(33,'(e12.5)') time
      write(34,'(e12.5)') time
      do j=3, np+2+nf-1       ! radial loop
         cja = 0d0
         do k2=kb-nag, kb+nag  ! average over angle 
           ko = k2
           apt    = apr(j,ko)*cos(abg(imag)*time)+api(j,ko)*sin(abg(imag)*time)
           bpt    = bphi(j,ko)*cos(abg(imag)*time)+iphi(j,ko)*sin(abg(imag)*time)
           bptkm  = bphi(j,ko-1)*cos(abg(imag)*time)+iphi(j,ko-1)*sin(abg(imag)*time)
           aptkm  = apr(j,ko-1)*cos(abg(imag)*time)+api(j,ko-1)*sin(abg(imag)*time)
           bptkp  = bphi(j,ko+1)*cos(abg(imag)*time)+iphi(j,ko+1)*sin(abg(imag)*time)
           aptkp  = apr(j,ko+1)*cos(abg(imag)*time)+api(j,ko+1)*sin(abg(imag)*time)
           bptjm  = bphi(j-1,ko)*cos(abg(imag)*time)+iphi(j-1,ko)*sin(abg(imag)*time)
           aptjm  = apr(j-1,ko)*cos(abg(imag)*time)+api(j-1,ko)*sin(abg(imag)*time)
           aptjm2  = apr(j-2,ko)*cos(abg(imag)*time)+api(j-2,ko)*sin(abg(imag)*time)
           bptjp  = bphi(j+1,ko)*cos(abg(imag)*time)+iphi(j+1,ko)*sin(abg(imag)*time)
           aptjp  = apr(j+1,ko)*cos(abg(imag)*time)+api(j+1,ko)*sin(abg(imag)*time)
           ! first derivative wrt theta
           dbtheta =(bptkp-bptkm)/dt/2 
           datheta =(aptkp-aptkm)/dt/2
           ! second derivative wrt theta
           datheta2 = (aptkp-2*apt+aptkm)/dt/dt
           dbtheta2 = (bptkp-2*bpt+bptkm)/dt/dt
           ! first derivative wrt x
           dax = (aptjp-aptjm)/hh/2
           dbx = (bptjp-bptjm)/hh/2
           ! second derivative wrt x 
           if (j.ne.np+2) then
             dax2=(aptjp-2d0*apt+aptjm)/hh/hh
           else
             dax2=(apt-2d0*aptjm+aptjm2)/hh/hh
           endif
  
           cheln=(2d0*apt*bpt/sin(ang(ko))**2 + apt*dbtheta*cos(ang(ko))/sin(ang(ko)) &
             & +datheta*dbtheta -bpt*datheta2)/xr(j)/xr(j) &
             & +(-bpt*dax+apt*dbx)/xr(j)+dax*dbx -bpt*dax2
  
           cja=cja+cheln   
         enddo
         write(33,'(e12.5)')  cja/rnor/2.d0/nag
       enddo
  
       do j=3, np+2+nf-1       ! radial loop
         cja = 0d0
         do k2=km-nag,km+nag  ! average over angle 
           ko = k2
           apt    = apr(j,ko)*cos(abg(imag)*time)+api(j,ko)*sin(abg(imag)*time)
           bpt    = bphi(j,ko)*cos(abg(imag)*time)+iphi(j,ko)*sin(abg(imag)*time)
           bptkm  = bphi(j,ko-1)*cos(abg(imag)*time)+iphi(j,ko-1)*sin(abg(imag)*time)
           aptkm  = apr(j,ko-1)*cos(abg(imag)*time)+api(j,ko-1)*sin(abg(imag)*time)
           bptkp  = bphi(j,ko+1)*cos(abg(imag)*time)+iphi(j,ko+1)*sin(abg(imag)*time)
           aptkp  = apr(j,ko+1)*cos(abg(imag)*time)+api(j,ko+1)*sin(abg(imag)*time)
           bptjm  = bphi(j-1,ko)*cos(abg(imag)*time)+iphi(j-1,ko)*sin(abg(imag)*time)
           aptjm  = apr(j-1,ko)*cos(abg(imag)*time)+api(j-1,ko)*sin(abg(imag)*time)
           aptjm2  = apr(j-2,ko)*cos(abg(imag)*time)+api(j-2,ko)*sin(abg(imag)*time)
           bptjp  = bphi(j+1,ko)*cos(abg(imag)*time)+iphi(j+1,ko)*sin(abg(imag)*time)
           aptjp  = apr(j+1,ko)*cos(abg(imag)*time)+api(j+1,ko)*sin(abg(imag)*time)
  
           ! first derivative wrt theta
           dbtheta =(bptkp-bptkm)/dt/2d0 
           datheta =(aptkp-aptkm)/dt/2d0
           ! second derivative wrt theta
           datheta2 = (aptkp-2d0*apt+aptkm)/dt/dt
           dbtheta2 = (bptkp-2d0*bpt+bptkm)/dt/dt
           ! first derivative wrt x
           dax= (aptjp-aptjm)/hh/2d0
           dbx= (bptjp-bptjm)/hh/2d0
           ! second derivative wrt x
  
           if (j.ne.np+2) then
             dax2=(aptjp-2d0*apt+aptjm)/hh/hh
           else
             dax2=(apt-2d0*aptjm+aptjm2)/hh/hh
           endif
  
           chels=(2d0*apt*bpt/sin(ang(ko))**2 + apt*dbtheta*cos(ang(ko))/sin(ang(ko)) &
             & +datheta*dbtheta -bpt*datheta2)/xr(j)/xr(j) &
             & +(-bpt*dax+apt*dbx)/xr(j)+dax*dbx -bpt*dax2
  
           cja=cja+chels
         enddo
         write(34,'(e12.5)')  cja/rnor/2.d0/nag
       enddo
    enddo
  
    close(33)
    close(34)

  end subroutine 

  subroutine compute_interior_solution (nst, net, nsp, nep, hh, x1, &
                                        bphi, iphi, brr, bri, apr, api, ffc, &
                                        sfu, ute, vc, ax, axp, bx, bxp, fx, & 
                                        om0, om0p, om2, om2p, om4, om4p, ome) 

    !--------------------------------------------------------
    ! > Compute interior solution for toroidal and poloidal 
    ! fields.
    !--------------------------------------------------------

    ! Arguments
    integer :: nst, net, nsp, nep 
    real(dp) :: hh, x1, ffc
    real(dp) :: om0, om0p, om2, om2p, om4, om4p
    real(dp) :: ax, axp, bx, bxp, fx

    real(dp) :: bphi(np+2+nft,n_theta)        ! toroidal real
    real(dp) :: iphi(np+2+nft,n_theta)        ! toroidal im
    real(dp) :: brr(np+2+nft,n_theta)         ! radial real
    real(dp) :: bri(np+2+nft,n_theta)         ! radial im
    real(dp) :: apr(np+2+nft,n_theta)         ! b poloidal (potential)
    real(dp) :: api(np+2+nft,n_theta)         ! b poloidal (potential)
    real(dp) :: sfu(np+2,n_theta)             ! stream function 
    real(dp) :: ome(np+2,n_theta)             ! omega
    real(dp) :: ute(np+2,n_theta)             ! utheta
    real(dp) :: vc(np+2,nb,2)                 ! eigenvector

    ! Local variables
    real(dp) :: theta, x
    integer :: j, k1, k2 

    bphi=0d0
    iphi=0d0
    brr=0d0
    bri=0d0
    apr=0d0
    api=0d0
    ! imp: interior solution
    ! imp: initialize to zero all the other quantities
    do k1=1, n_theta    ! theta loop
      theta = theta1 + k1*dt
      do j=1, np+2         ! radial loop
        x = x1+(j-1d0)*hh
        !  load the bphi_n  vectors in b_n p1(cos(theta))_n
        do k2=nst, net, 2    ! legendre polynomial  l-loop
          if (mm.eq.0) then
            !Note: in spherical  symmetry Btoroidal = - \partial_\theta
            !T(r,theta)  where  T=\sum_n P^0_n(cos(theta))  and P^0_n are the associate
            !Legendre functions of the first kind (see maple and mathematical definition
            !with the correct branch cut.) 
            !Therefore we can use the (Ruediger corrected) fundamental identity:  
            !\partial_\theta P^0_n(cos(theta)), theta) = (plus) P^1_n(cos(theta))
            bphi(j,k1)=bphi(j,k1)+vc(j,k2,1)*(-plgndr(k2,1,cos(theta)))
            iphi(j,k1)=iphi(j,k1)+vc(j,k2,2)*(-plgndr(k2,1,cos(theta)))
          else
            ! mm not zero, but should be checked! 
            write (*,*) 'check please m not zero'
            stop 
            bphi(j,k1)=bphi(j,k1)+vc(j,k2,1)*(- plgndr(k2,1,cos(theta))*(-cos(theta)/sin(theta))+k2*(k2+1)*plgndr(k2,0,cos(theta)))
            iphi(j,k1)=iphi(j,k1)+vc(j,k2,2)*(- plgndr(k2,1,cos(theta))*(-cos(theta)/sin(theta))+k2*(k2+1)*plgndr(k2,0,cos(theta)))
          endif 
        enddo
        ! load the a_n vectors at the boundary...
        do k2=nsp, nep, 2   ! Legendre polynomial even l-loop         
          !
          !  obtain B_r in the interior !!!
          !  to get B_r  use second fundamental identity: (see mathematica notebook notation.nb)
          !
          ! B_r = -(1/sintheta) \partial_\theta (sintheta\partial_\theta p0n) = 
          !     = -(1/sintheta) \partia_\theta [ p1n sintheta ] = +n*(n+1) p0n
          ! 
          !  so that  (note a minus sign in the definition of Br and Bphi) 
          !   
          !  B_r = (sum over n) + n(n+1) P^0_n / x
          vc(np+2,k2,1) = (4d0*vc(np+1,k2,1)-vc(np,k2,1))/(2d0*hh*(k2+1+ffc)+3d0)  
          vc(np+2,k2,2) = (4d0*vc(np+1,k2,2)-vc(np,k2,2))/(2d0*hh*(k2+1+ffc)+3d0)    
          brr(j,k1) = brr(j,k1)+(k2*(k2+1)/x)*vc(j,k2,1)*(+plgndr(k2,0,cos(theta)))
          bri(j,k1) = bri(j,k1)+(k2*(k2+1)/x)*vc(j,k2,2)*(+plgndr(k2,0,cos(theta)))
  
        enddo ! legendre polynomial even l-loop  k2
  
        call rot(x, om0, om0p, om2, om2p, om4, om4p)
        call stream(x, ax, axp, bx, bxp, fx)
        if (regime.eq.'h4' .or. regime.eq.'h6' .or. regime.eq.'h5') then 
          ome(j,k1) = om0 + om2*cos(theta)**2 + om4*cos(theta)**4
        else
          ome(j,k1) = om0 + om2 *( -(1d0-3d0*cos(theta)**2)/2.0d0 )
        endif
        sfu(j,k1) = fx*sin(theta)*cos(theta)
        ute(j,k1) = bx*sin(theta)*cos(theta) 
  
      enddo   ! close radial  
    enddo   ! close theta

  end subroutine compute_interior_solution
  
  subroutine writefield (imag)

    !------------------------------------------------------
    ! > Write the fields computed when solving the eigenvalue
    ! problem. 
    !------------------------------------------------------

    ! Arguments
    real(dp) :: imag ! imaginary part of the cycle eigenvalue

    ! Local variables
    real(dp) :: aai, aar, ax, axp, bbi, bbr, bpt, bx, bxp, &
            & dbx, ffc, fx, h2, hh, &
            & om0, om0p, om2, om2p, om4, om4p, &
            & tc, theta, x, x1, x2, xf
    integer :: i, i2, i3, jj, k, &
               & nep, net, nf, nsp, nst
  
    character(len=2) :: q
    character(len=512) :: bfeld1, bfeld2, bfeld3, bfeld4, bfeld5, &
                          bfeld6, bfeld7, bfeld8, bfeld9, bfeld10
  
    real(dp) :: bphi(np+2+nft,n_theta)        ! toroidal real
    real(dp) :: iphi(np+2+nft,n_theta)        ! toroidal im
    real(dp) :: brr(np+2+nft,n_theta)         ! radial real
    real(dp) :: bri(np+2+nft,n_theta)         ! radial im
    real(dp) :: apr(np+2+nft,n_theta)         ! b poloidal (potential)
    real(dp) :: api(np+2+nft,n_theta)         ! b poloidal (potential)
    real(dp) :: ome(np+2,n_theta)             ! omega
    real(dp) :: sfu(np+2,n_theta)             ! stream function 
    real(dp) :: ute(np+2,n_theta)             ! utheta
    real(dp) :: vc(np+2,nb,2)                 ! eigenvector
    real(dp) :: xr(np+2+nft)                
  
    integer :: bin(2), bax(2), rin(2), rax(2)
    
    ! initialize
    mm = int (mmm)
    vc = 0d0
    !  important mmm can only be le 1 in this subroutine !!!! 
    if (mm.ge.2) write (*, *) 'Be careful that m is larger than 1 !' 
    x1 = x_in           !inner boundary
    x2 = 1.0d0          !outer bound
    hh =(x2-x1)/(np+1)  !stepsize: radial accuracy parm.
    h2 = hh/2.e0
  
    if (mod(mm,2).eq.0) then 
      if (degree.eq.'d') then
        q='.a'
        nst=2
        net=nb
        nsp=1
        nep=na
      else
        q='.s'
        nst=1
        net=na
        nsp=2
        nep=nb
      endif
    else if (mod(mm,2).eq.1) then
      if (degree.eq.'d') then
        q='.s'
        nst=2
        net=nb
        nsp=1
        nep=na
      else
        q='.a'
        nst=1
        net=na
        nsp=2
        nep=nb
      endif
    endif
  
    if (flg.eq.1) then
      do i2=1, nb
        k = 1          
        do i3 = i2, nt, nb
          k = k + 1   
          vc(k,i2,1) = real(cvr(i3,int(indeg(nt))))
          vc(k,i2,2) = aimag(cvr(i3,int(indeg(nt))))
        enddo
      enddo
    else   ! flg = 0 (real matrix inversion)
      do i2=1, nb
        k = 1          
        do i3=i2, nt, nb
          k = k + 1   
          vc(k,i2,1) = vr(i3,indeg(nt))
          ! it should be -  (after long check with lapack libraries)
          vc(k,i2,2) = -vr(i3,indeg(nt-1))
        enddo
      enddo
    endif   ! flg mode

    call write_vect (bfeld2, q)
    call compute_interior_solution (nst, net, nsp, nep, hh, x1, &
                                    bphi, iphi, brr, bri, apr, api, ffc, &
                                    sfu, ute, vc, ax, axp, bx, bxp, fx, &
                                    om0, om0p, om2, om2p, om4, om4p, ome)  
    call write_omega (np, n_theta, x1, x2, theta1, theta2, &
                      ome, sfu, ute)
  
    ! determine the outer mesh point equal to zeta_r+1
    if (zeta_r .lt. 1.1d0) zeta_r = 1.1d0 
    xf = zeta_r
    nf = int((xf-1.0d0)/hh)
    if(nf .ge. nft) then
      write (*, *) 'too small nft'
      stop
    endif
  
    ! the exterior solution is different if ffree is not zero (pure force-free)
    ! parameter, or beta is zero or not zero (potential vs helmholtz extrapolation)    
    call compute_pol (nsp, nep, nf, nft, vc, apr, api, &
                      ffc, hh, x1)
    call compute_tor (nst, net, nf, nft, vc, bphi, iphi, ffc, &
                      hh, x1)
    call compute_ab_vector (vc, aai, aar, bbi, bbr, &
                            ffc, x1, hh, nsp, nep, &
                            nst, net, nf)
  
    ! nj is the number of time slices
    do jj=1, nj
      call write_toroidal (bfeld3, imag, nf, nj, jj, q, &
                           theta, x1, xf, tc, bphi, iphi)
      call write_poloidal (bfeld4, imag, nf, nj, jj, q, &
                           theta, hh, x1, xf, tc, apr, api) 
    enddo

    call time_butterfly (bfeld2, bfeld6, bfeld7, bfeld8, &
                         imag, nf, nft, x1, x2, hh, bphi, iphi, &
                         brr, bri, apr, api, q) 
    call radial_butterfly (bfeld9, bfeld10, imag, nf, nft, x1, x2, &
                           hh, bphi, iphi, brr, bri, apr, api, q) 
  
  end subroutine writefield

end module write_outputs
