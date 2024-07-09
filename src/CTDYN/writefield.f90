! 
!  IMPORTANT: the eigenfunctions are always complex conjugate variable. If the eigenvalue has a positive Im part, 
!  be sure you keep the correct conjugate eigenvector
! 
!
!  IMPORTANT:
!
! The convetions for the vector potential and scalar potentials are in the Maple ws called currenthelicity.mw (and similar)
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

subroutine writefield

  implicit none

  integer :: ii,it                 ! iteration counter
  common/ipar/ii,it

  real :: xa1,xa2,xa3,xb,xda1,xda2
  common/apar/xa1,xa2,xa3,xb,xda1,xda2
  real :: edr,xe1,xde1
  common/epar/edr,xe1,xde1,hd
  real :: x_in, bct, c3, mmm,sr,rotp,gd, aqu, flg
  common/ppar/x_in,bct,c3,mmm,sr,rotp,gd, aqu, flg

  real :: s0,s2,s4,s6, a2p, a4p,xm
  common/psi/s0,s2,s4,s6,a2p,a4p,xm

  character*102 :: dir 
  character*8 :: ans1,ans2,ans3,ans4 
  common/var3/ans1,ans2,ans3,ans4, dir

  character*2 :: jobvr,jobvl
  common/lap/jobvr,jobvl

  character*30 :: inp
  character*43 :: version, ver
  character*82 :: fel,fes,fel2,fes2,fen
  character*2 :: qq, q

  real :: imag
  common/part/vtu,rt,imag,co,c_u,beta,ffree,betb,etep,etet,xbt,xbo

  character*512 :: bfeld1,bfeld2,bfeld3,bfeld4,bfeld5, &
                   & bfeld6,bfeld7,bfeld8,bfeld9,bfeld10

  include 'cio'
  include 'ccov'

  common/vec/cvr
  common/parker/gam,zeta_r,ratio

  parameter(nft = 200)   

  real :: ang(n_theta)              

  real :: bphi(np+2+nft,n_theta), &        ! toroidal real
     &     iphi(np+2+nft,n_theta), &        ! toroidal im
     &     brr(np+2+nft,n_theta), &        ! radial real
     &     bri(np+2+nft,n_theta) , &        ! radial im
     &     apr(np+2+nft,n_theta), &        ! b poloidal (potential)
     &     api(np+2+nft,n_theta) , &        ! b poloidal (potential)
     &     ome(np+2,n_theta) , &        ! omega
     &     sfu(np+2,n_theta) , &        ! stream function 
     &     ute(np+2,n_theta) , &        ! utheta
     &     vc(np+2,nb,2), &         ! eigenvector
     &     xr(np+2+nft)                

  external plgndr
  integer :: bin(2), bax(2), rin(2), rax(2)

  ! initialize

  vc=0
  mm = int(mmm)
  if(mm.gt.2)print*, 'not really checked!!!!' 
  !  important mmm can only be le 1 in this subroutine !!!! 
  x1 = x_in           !inner boundary
  x2 = 1.0            !outer bound
  hh =(x2-x1)/(np+1)  !stepsize: radial accuracy parm.
  h2 = hh/2.e0

  if(mod(mm,2).eq.0)then 
    if(ans2.eq.'d')then
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
  else if(mod(mm,2).eq.1)then
    if(ans2.eq.'d')then
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
  jj = 1

  open(50,status='unknown',file=trim(dir)//'/omega.dat')
  open(51,status='unknown',file=trim(dir)//'/stream.dat')
  open(52,status='unknown',file=trim(dir)//'/utheta.dat')
  write(bfeld2,2010) trim(dir)//'/vect.',ii,q,mm 
  open(13,status='unknown',file=adjustl(bfeld2))

  if(flg.eq.1)then
    do i2 =1, nb
      k = 1          
      do i3 = i2, nt, nb
        k = k + 1   
        vc(k,i2,1) = real(cvr(i3,int(indeg(nt))))
        vc(k,i2,2) = aimag(cvr(i3,int(indeg(nt))))
      enddo
    enddo
  else   ! flg = 0 (real matrix inversion)
    do i2 =1,nb
      k = 1          
      do i3 = i2,nt,nb
        k = k + 1   
        vc(k,i2,1) = vr(i3,indeg(nt))
        ! it should be -  (after long check with lapack libraries)
        vc(k,i2,2) = -vr(i3,indeg(nt-1))
      enddo
    enddo
  endif   ! flg mode
  close(13)

  theta1 = 0.
  theta2 = 3.1415926
  theta = theta1-theta2/(n_theta)
  bphi=0
  iphi=0
  brr=0
  bri=0
  apr=0
  api=0
  adum=0

  ! imp: interior solution
  ! imp: initialize to zero all the other quantities

  do k1=1, n_theta    ! theta loop
    theta = theta + theta2/(n_theta)
    ! only for testing
    !       bdum(k1)=cos(theta)
    do j =1,np+2         ! radial loop
      x = x1+(j-1)*hh
      !  load the bphi_n  vectors in b_n p1(cos(theta))_n
      do k2 = nst, net,2    ! legendre polynomial  l-loop
        if (mm.eq.0) then
          !note: in spherical  symmetry btoroidal = - \partial_\theta
          !t(r,theta)  where  t=\sum_n p^0_n(cos(theta))  and p^0_n are the associate
          !legendre functions of the first kind (see maple and mathematical definition
          !with the correct branch cut.) 
          !Therefore we can use the (ruediger corrected) fundamental identity:  
          !\partial_\theta p^0_n(cos(theta)), theta) = (plus) p^1_n(cos(theta))
          bphi(j,k1)=bphi(j,k1)+vc(j,k2,1)*(-plgndr(k2,1,cos(theta)))
          iphi(j,k1)=iphi(j,k1)+vc(j,k2,2)*(-plgndr(k2,1,cos(theta)))
        else
          ! mm not zero, but should be checked! 
          print*, 'check please m not zero'
          stop 
          bphi(j,k1)=bphi(j,k1)+vc(j,k2,1)*(- plgndr(k2,1,cos(theta))*(-cos(theta)/sin(theta))+k2*(k2+1)*plgndr(k2,0,cos(theta)))
          iphi(j,k1)=iphi(j,k1)+vc(j,k2,2)*(- plgndr(k2,1,cos(theta))*(-cos(theta)/sin(theta))+k2*(k2+1)*plgndr(k2,0,cos(theta)))
        endif 
      enddo
      ! load the a_n vectors at the boundary...
      do k2 = nsp,nep,2   ! legendre polynomial even l-loop         
        !
        !  obtain B_r in the interior !!!
        !  to get B_r  use second fundamental identity: (see mathematica notebook notation.nb)
        !
        ! B_r = -(1/sintheta) \partial_\theta (sintheta\partial_\theta p0n) = 
        !     = -(1/sintheta) \partia_\theta [ p1n sintheta ] = +n*(n+1) p0n
        ! 
        !  so that  (note a minus sign in the definition of br and bphi) 
        !   
        !  B_r = (sum over n) + n(n+1) p0n / x

        vc(np+2,k2,1)   = (4*vc(np+1,k2,1)-vc(np,k2,1))/(2*hh*(k2+1+ffc)+3)  
        vc(np+2,k2,2)   = (4*vc(np+1,k2,2)-vc(np,k2,2))/(2*hh*(k2+1+ffc)+3)    

        brr(j,k1)   = brr(j,k1)+(k2*(k2+1)/x)*vc(j,k2,1)*(+plgndr(k2,0, cos(theta)))
        bri(j,k1)   = bri(j,k1)+(k2*(k2+1)/x)*vc(j,k2,2)*(+plgndr(k2,0, cos(theta)))

      enddo ! legendre polynomial even l-loop  k2

      call rot(x, om0, om0p, om2, om2p, om4, om4p)
      call stream(x,ax,axp,bx,bxp,fx)
      if (ans1.eq.'h4' .or. ans1.eq.'h6' .or. ans1.eq.'h5') then 
        ome(j,k1) = om0 + om2*cos(theta)**2 + om4*cos(theta)**4
      else
        ome(j,k1) = om0 + om2 *( -(1-3*cos(theta)**2)/2.0 )
      endif
      sfu(j,k1) = fx*sin(theta)*cos(theta)
      ute(j,k1) = bx*sin(theta)*cos(theta) 

    enddo   ! close radial  
  enddo   ! close theta

  include 'write_omega.f'

  ! determine the outer mesh point equal to zeta_r+1

  xf=zeta_r+0.5
  xf=zeta_r+0.3
  nf=int((xf-1.0)/hh)
  if(nf .ge. nft) then
    print*, 'too small nft'
    stop
  endif

  ! the exterior solution is different if ffree is not zero (pure force-free)
  ! parameter, or beta is zero or not zero (potential vs helmholtz extrapolation)    
    
  include 'load_pol.f'   
  include 'load_tor.f'
  include 'load_vec.f'

  !  ---------- write toroidal  retor.dat ed itor
  !  nj is the number of time slices
  nj = 8

  do jj =1,nj
    include 'write_toroidal.f'
    include 'write_poloidal.f'
  enddo

  !  butterfly diagram block      
  !  time evolution butterfly diagram 
  write(bfeld1,2010) trim(dir)//'/butf.',ii,q,mm 
  open(18,status='unknown',file=adjustl(bfeld1))      

  write(bfeld6,2010) trim(dir)//'/jcin.',ii,q,mm 
  open(30,status='unknown',file=adjustl(bfeld6))      

  write(bfeld7,2010) trim(dir)//'/jcou.',ii,q,mm 
  open(31,status='unknown',file=adjustl(bfeld7))      

  write(bfeld8,2010) trim(dir)//'/brbp.',ii,q,mm 
  open(32,status='unknown',file=adjustl(bfeld8))      

  t_in = 0.
  t_fin = 12.
  time = -t_fin/n_time
  theta1 = 0.
  theta2 = 3.1415926
  theta = theta1-theta2/(n_theta)

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

  write(*, '(a,f10.4)') 'max poloidal / max toroidal  ',  ratio
  write(*,'(a,f10.4)') 'bphi max location', bax(1)*(x2-x1)/(np+1) +x1
  write(*,'(a,f10.4)') 'a max location',   rax(1)*(x2-x1)/(np+1) +x1

  do k2 =2,n_theta-1
    theta = theta + theta2/(n_theta)
    write(18,'(e12.5)') theta
    write(30,'(e12.5)') theta
    write(31,'(e12.5)') theta
    write(32,'(e12.5)') theta
  enddo

  x=0
  do j =1,np+2 + nf        ! radial loop
    x = x1+(j-1)*hh
    xr(j) = x
  enddo

  call hunt(xr,np+2+nf, xbt, jbt)
  call hunt(xr,np+2+nf, xbo, jbo)

  do k1=1,n_time
    time = time + t_fin/n_time
    write(18,'(e12.5)') time
    write(30,'(e12.5)') time
    write(31,'(e12.5)') time
    write(32,'(e12.5)') time
    theta=0
    do k2 =2, n_theta-1
      theta = theta + theta2/(n_theta)
      bphi7= bphi(jbt,k2)*cos(abg(imag)*time)+iphi(jbt,k2)*sin(abg(imag)*time)
      write(18,'(e12.5)')  bphi7
      dt=theta2/n_theta
      j=jbt
      brs  = brr(np,k2)*cos(abg(imag)*time)+bri(np,k2)*sin(abg(imag)*time)
      write(32,'(e12.5)') brs !* bphi7
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
      write(30,'(e12.5)') chel

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
      write(31,'(e12.5)')  chel2
    enddo
  enddo

2010   format(a220,i6.6,a2,i2.2) 
2011   format(a220,i6.6,a2,i2.2,a2,i2.2)
           
  close(18)
  close(30)
  close(31)
  close(32)

!----------------------------------------------------------
! radial butterfly diagram n and s 
!----------------------------------------------------------

  bax=maxloc(abs(bphi))
  bin=minloc(abs(bphi))
  zq1= bphi(bax(1), bax(2))
  bax=maxloc(abs(iphi))
  bin=minloc(abs(iphi))
  zq2=iphi(bax(1), bax(2))
  rnor=((zq1+zq2)/2.)**2
  
  bax=maxloc(abs(apr))
  bin=minloc(abs(apr))
  zq1= bphi(bax(1), bax(2))
  bax=maxloc(abs(api))
  bin=minloc(abs(api))
  zq2=iphi(bax(1), bax(2))
  rnor=((zq1+zq2)/2.)**2

  write(bfeld9,2010) trim(dir)//'/brtn.',ii,q,mm 
  open(33,status='unknown',file=adjustl(bfeld9))      
  write(bfeld10,2010) trim(dir)//'/brts.',ii,q,mm 
  open(34,status='unknown',file=adjustl(bfeld10))      

  t_in = 0.
  t_fin = 12.
  time = -t_fin/n_time

  write(33,'(6x,i4)') np+2+nf-2
  write(34,'(6x,i4)') np+2+nf-2

  x=0
  do j =2,np+2+nf-1       ! radial loop
    x = x1+(j-1)*hh
    write(33,'(e12.5)') x
    write(34,'(e12.5)') x
  enddo

  theta=3.1415926/2.0
  theta=0
  do k2 =1,n_theta
    theta = theta + theta2/n_theta
    ang(k2)=theta
  enddo

  nag=2
  abt=1.047
  abm=2.094

  call hunt(ang,n_theta, abt, kb)
  call hunt(ang,n_theta, abm, km)

  dt=theta2/n_theta

  do k1=1,n_time
    time = time + t_fin/n_time
    write(33,'(e12.5)') time
    write(34,'(e12.5)') time
    do j =3,np+2+nf-1       ! radial loop
       cja=0
       do k2=kb-nag,kb+nag  ! average over angle 
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
         dax= (aptjp-aptjm)/hh/2
         dbx= (bptjp-bptjm)/hh/2
         ! second derivative wrt x 
         if (j.ne.np+2) then
           dax2=(aptjp-2*apt+aptjm)/hh/hh
         else
           dax2=(apt-2*aptjm+aptjm2)/hh/hh
         endif

         cheln=(2*apt*bpt/sin(ang(ko))**2 + apt*dbtheta*cos(ang(ko))/sin(ang(ko)) &
           & +datheta*dbtheta -bpt*datheta2)/xr(j)/xr(j) &
           & +(-bpt*dax+apt*dbx)/xr(j)+dax*dbx -bpt*dax2

         cja=cja+cheln   
       enddo
       write(33,'(e12.5)')  cja/rnor/2./nag
     enddo

     do j =3,np+2+nf-1       ! radial loop
       cja=0
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
         dbtheta =(bptkp-bptkm)/dt/2 
         datheta =(aptkp-aptkm)/dt/2
         ! second derivative wrt theta
         datheta2 = (aptkp-2*apt+aptkm)/dt/dt
         dbtheta2 = (bptkp-2*bpt+bptkm)/dt/dt
         ! first derivative wrt x
         dax= (aptjp-aptjm)/hh/2
         dbx= (bptjp-bptjm)/hh/2
         ! second derivative wrt x

         if (j.ne.np+2) then
           dax2=(aptjp-2*apt+aptjm)/hh/hh
         else
           dax2=(apt-2*aptjm+aptjm2)/hh/hh
         endif

         chels=(2*apt*bpt/sin(ang(ko))**2 + apt*dbtheta*cos(ang(ko))/sin(ang(ko)) &
           & +datheta*dbtheta -bpt*datheta2)/xr(j)/xr(j) &
           & +(-bpt*dax+apt*dbx)/xr(j)+dax*dbx -bpt*dax2

         cja=cja+chels
       enddo
       write(34,'(e12.5)')  cja/rnor/2./nag
     enddo
  enddo

  close(33)
  close(34)

end


