C 
C  IMPORTANT: the eigenfunctions are always complex conjugate variable. If the eigenvalue has a positive Im part, 
C  be sure you keep the correct conjugate eigenvector
C 
C
C  IMPORTANT:
C
C The convetions for the vector potential and scalar potentials are in the Maple ws called currenthelicity.mw (and similar)
C and on the mathematica notebook notation.nb
C
C In spherical symmetry in order to make contact with the standard representation used in papers:
C
C Br= (1/r sintheta) \partial_\theta [A(r,theta) sintheta]
C 
C Btheta =(-1/r) \partial_r [A r]
C
C Bphi = B
C
C you must remember that, according to the Radler formalism of the Psi and Phi functions we have:
C    
C      B = -\partial_\theta P0n Psi_n = -P1n Psi_n
C and
C      A = -\partial_\theta P0n Phi_n = -P1n Phi_n
C
C    where we used the identity: \partial P0n = P1n to  make contact with the notations in the papers 
C    (note that the above identity differs by a sign from what you find in Rudiger's book in the appendix)
C
C
      subroutine writefield

      IMPLICIT REAL(A-H,O-Z)

      INTEGER II,IT                 ! ITERATION COUNTER
      COMMON/ipar/II,IT

      REAL xa1,xa2,xa3,xb,xda1,xda2
      COMMON/apar/xa1,xa2,xa3,xb,xda1,xda2
      REAL edr,xe1,xde1
      COMMON/epar/edr,xe1,xde1,hd
      REAL x_in, bct, c3, mmm,SR,rotp,gd, aqu, flg
      COMMON/ppar/x_in,bct,c3,mmm,SR,rotp,gd, aqu, flg

      REAL s0,s2,s4,s6, A2P, A4P,xm
      COMMON/PSI/s0,s2,s4,s6,A2P,A4P,xm

      CHARACTER*102 dir 
c      CHARACTER*72 dir 
      CHARACTER*8 ANS1,ANS2,ANS3,ANS4 
      COMMON/var3/ANS1,ANS2,ANS3,ANS4, dir

      CHARACTER*2 JOBVR,JOBVL
      COMMON/lap/JOBVR,JOBVL
 
      CHARACTER*30 inp
      CHARACTER*43 version, ver
      CHARACTER*82 fel,fes,fel2,fes2,fen
      CHARACTER*2 qq, q

c      REAL vtu,rt,imag,co,c_u,beta,ffree
      REAL imag
      COMMON/part/vtu,rt,imag,CO,C_U,beta,ffree,betb,etep,etet,xbt,xbo

c      CHARACTER*52 bfeld1,bfeld2,bfeld3,bfeld4,bfeld5
c     & ,bfeld6,bfeld7,bfeld8,bfeld9,bfeld10

      CHARACTER*512 bfeld1,bfeld2,bfeld3,bfeld4,bfeld5
     & ,bfeld6,bfeld7,bfeld8,bfeld9,bfeld10

      include 'cio'
      include 'ccov'

      COMMON/vec/CVR
      COMMON/parker/gam,zeta_r,ratio

c      PARAMETER(NFT = 400)   
      PARAMETER(NFT = 200)   
c      PARAMETER(NFT = 500)   
c      PARAMETER(NFT = 1200)   

      REAL ang(N_THETA)              

      REAL BPHI(NP+2+NFT,N_THETA),         ! toroidal real
     &     IPHI(NP+2+NFT,N_THETA),         ! toroidal im
     &     BRR(NP+2+NFT,N_theta),         ! radial real
     &     BRI(NP+2+NFT,N_THETA) ,         ! radial im
     &     APR(NP+2+NFT,N_theta),         ! B poloidal (potential)
     &     API(NP+2+NFT,N_THETA) ,         ! B poloidal (potential)
     &     OME(NP+2,N_THETA) ,         ! OMEGA
     &     SFU(NP+2,N_THETA) ,         ! stream function 
     &     UTE(NP+2,N_THETA) ,         ! utheta
     &     VC(NP+2,NB,2)    ,          ! eigenvector
     &     XR(NP+2+NFT)                


      EXTERNAL PLGNDR
c      real bin(2), bax(2), rin(2), rax(2)
      INTEGER bin(2), bax(2), rin(2), rax(2)

c initialize

       VC=0

      mm = int(mmm)

      if(mm.gt.2)print*, 'NOT REALLY CHECKED!!!!' 

CCC   IMPORTANT mmm can only be LE 1 in this subroutine !!!! 

      x1 = x_in           !inner boundary
      x2 = 1.0            !outer bound
      hh =(x2-x1)/(NP+1)  !stepsize: radial accuracy parm.
      h2 = hh/2.e0

       if(mod(mm,2).eq.0)then 
          if(ans2.eq.'D')then
            q='.A'
            nst=2
            net=nb
            nsp=1
            nep=na
          else
            q='.S'
            nst=1
            net=na
            nsp=2
            nep=nb
          endif
        else if(mod(mm,2).eq.1)then
          if(ans2.eq.'D')then
           q='.S'
           nst=2
           net=nb
           nsp=1
           nep=na
          else
           q='.A'
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

c      if(flg.eq.1)then

c        write(13,'(1x,I4,2x,I4)')    NP, NB
c        write(13,'(1x,f9.5,2x,f9.5)') x1,x2
c        do i3 =1, NT
c        write(13,'(e15.5)') REAL(CVR(i3,INT(INDEG(NT))))    ! WRITE FIRST EIGENVECTOR
c        enddo
c        do i3 =1, NT
c        write(13,'(e15.5)') AIMAG(CVR(i3,INT(INDEG(NT))))   ! WRITE SECOND / IM part of first
c        enddo

c         write(*,*) flg

      if(flg.eq.1)then
        do i2 =1, NB
         k = 1          
         do i3 = i2, NT, NB
          k = k + 1   
          VC(k,i2,1) = REAL(CVR(i3,INT(INDEG(NT))))
          VC(k,i2,2) = AIMAG(CVR(i3,INT(INDEG(NT))))
         enddo
        enddo

      else   ! flg = 0 (real matrix inversion)

        do i2 =1,NB
        k = 1          
        do i3 = i2,NT,NB
        k = k + 1   
        VC(k,i2,1) = VR(i3,INDEG(NT))
c it should be -  (after long check with lapack libraries)
           VC(k,i2,2) = -VR(i3,INDEG(NT-1))
c        write(13,'(I4, 3e15.5)') i2, x1+(k+1)*hh, VC(k,i2,1) , VC(k,i2,2)

        enddo
        enddo

       endif   ! flg mode

       close(13)
c       print*, 'imag', imag 
c
      theta1 = 0.
      theta2 = 3.1415926
      theta = theta1-theta2/(n_theta)
      BPHI=0
      IPHI=0
      BRR=0
      BRI=0
      APR=0
      API=0
      ADUM=0


C IMP: INTERIOR SOLUTION

C IMP: initialize to zero all the other quantities

      DO K1=1, n_theta    ! theta loop

        theta = theta + theta2/(n_theta)

c only for testing
c       BDUM(k1)=cos(theta)

       DO J =1,NP+2         ! radial loop
       x = x1+(J-1)*hh

c  LOAD THE BPHI_n  vectors in b_n P1(COS(THETA))_n

        DO k2 = nst, net,2    ! Legendre Polynomial  l-loop
         if(mm.eq.0)then
C
C note: in spherical  symmetry Btoroidal = - \partial_\theta T(r,theta)  
C 
C where  T=\sum_n P^0_n(cos(theta))  and P^0_n are the associate legendre functions of the first kind (see Maple and Mathematical definition
C with the correct branch cut.) 
C 
C therefore we can use the (ruediger corrected) FUNDAMENTAL identity:  \partial_\theta P^0_n(cos(theta)), theta) = (PLUS) P^1_n(cos(theta))
C
         BPHI(J,k1)=BPHI(J,k1)+VC(J,k2,1)*(-plgndr(k2,1,cos(theta)))
         IPHI(j,k1)=IPHI(j,k1)+VC(J,k2,2)*(-plgndr(k2,1,cos(theta)))
C
         else
c mm not zero, but should be checked! 
        print*, 'check please m not zero'
        stop 
 
         BPHI(J,k1)=BPHI(J,k1)+VC(J,k2,1)*(- plgndr(k2,1,cos(theta))*(-cos(theta)/sin(theta))+k2*(k2+1)*plgndr(k2,0,cos(theta)))
         IPHI(j,k1)=IPHI(j,k1)+VC(J,k2,2)*(- plgndr(k2,1,cos(theta))*(-cos(theta)/sin(theta))+k2*(k2+1)*plgndr(k2,0,cos(theta)))

         endif 

CCC  cambia questa BC

c         BPHI(1,k1) = 0 ! change with pc 
c         IPHI(1,k1) = 0 ! change with pc
        ENDDO
               ! LOAD THE a_n vectors at the boundary...
        DO k2 = nsp,nep,2   ! Legendre Polynomial even l-loop         
C
C  obtain Bradial IN THE INTERIOR !!!
C  to get Bradial  USE second FUNDAMENTAL identity: (see mathematica notebook notation.nb)
C
C Br = -(1/sintheta) \partial_\theta (sintheta\partial_\theta P0n) = 
C    = -(1/sintheta) \partia_\theta [ P1n sintheta ] = +n*(n+1) P0n
C 
C  so that  (note a minus sign in the definition of Br and Bphi) 
C   
C  Bradial = (sum over n) + n(n+1) P0n / x
C     
C          VC(1,k2,1)      =  0
C          VC(1,k2,2)      =  0
C

          VC(NP+2,k2,1)   = (4*VC(NP+1,k2,1)-VC(NP,k2,1))/(2*hh*(k2+1+ffc)+3)  
          VC(NP+2,k2,2)   = (4*VC(NP+1,k2,2)-VC(NP,k2,2))/(2*hh*(k2+1+ffc)+3)    

          BRR(J,k1)   = BRR(J,k1)+(k2*(k2+1)/x)*VC(J,k2,1)*(+plgndr(k2,0, cos(theta)))
          BRI(J,k1)   = BRI(J,k1)+(k2*(k2+1)/x)*VC(J,k2,2)*(+plgndr(k2,0, cos(theta)))

        ENDDO ! Legendre Polynomial even l-loop  K2

            call rot(x, om0, om0p, om2, om2p, om4, om4p)
            call stream(x,ax,axp,bx,bxp,fx)
            if(ans1.eq.'H4'.or.ans1.eq.'H6'.or.ans1.eq.'H5')then 
            OME(J,K1) = om0 + om2*cos(theta)**2 + om4*cos(theta)**4
            else
            OME(J,K1) = om0 + om2 *( -(1-3*cos(theta)**2)/2.0 )
            endif
            SFU(J,K1) = fx*sin(theta)*cos(theta)
            UTE(J,K1) = bx*sin(theta)*cos(theta) 

         ENDDO   ! close radial  
         ENDDO   ! close theta

        INCLUDE 'write_omega.f'

C determine the outer mesh point equal to zeta_r+1

      xf=zeta_r+0.5
      xf=zeta_r+0.3
      NF=int((xf-1.0)/hh)
      if(NF .ge. NFT) then
      print*, 'too small NFT'
      stop
      endif
C
C the exterior solution is different if ffree is not zero (pure force-free)
C parameter, or beta is zero or not zero (potential vs helmholtz extrapolation)    
C
    
          INCLUDE 'load_pol.f'   
          INCLUDE 'load_tor.f'
          INCLUDE 'load_vec.f'


C  ---------- write toroidal  retor.dat ed itor
c  nj is the number of time slices
      nj = 8

      do jj =1,nj

       INCLUDE 'write_toroidal.f'
       INCLUDE 'write_poloidal.f'

      enddo

C  BUTTERFLY DIAGRAM BLOCK      
C  time evolution butterfly diagram 

      write(bfeld1,2010) trim(dir)//'/butf.',ii,q,mm 
      open(18,status='unknown',file=adjustl(bfeld1))      

      write(bfeld6,2010) trim(dir)//'/jcin.',ii,q,mm 
      open(30,status='unknown',file=adjustl(bfeld6))      

      write(bfeld7,2010) trim(dir)//'/jcou.',ii,q,mm 
      open(31,status='unknown',file=adjustl(bfeld7))      

      write(bfeld8,2010) trim(dir)//'/brbp.',ii,q,mm 
      open(32,status='unknown',file=adjustl(bfeld8))      

C
      t_in = 0.
      t_fin = 12.
      time = -t_fin/n_time
C
      theta1 = 0.
      theta2 = 3.1415926
      theta = theta1-theta2/(n_theta)

      write(18,'(6x,I4)') n_theta-2
      write(30,'(6x,I4)') n_theta-2
      write(31,'(6x,I4)') n_theta-2
      write(32,'(6x,I4)') n_theta-2

C     location of max bphi

        bax=maxloc(sqrt(bphi**2+Iphi**2))
        bin=minloc(sqrt(bphi**2+Iphi**2))
        rax=maxloc(sqrt(apr**2+api**2))
        rin=minloc(sqrt(apr**2+api**2))

c        write(*, '(a,f10.4)') 'max poloidal / max toroidal  ' , 
c     .  (sqrt(bphi(bax(1), bax(2))**2 + iphi(bax(1), bax(2))**2)
c     .  /sqrt(apr(rax(1), rax(2))**2 + api(rax(1), rax(2))**2) )**(-1)

         ratio= (sqrt(bphi(bax(1), bax(2))**2 + iphi(bax(1), bax(2))**2)
     .   /sqrt(apr(rax(1), rax(2))**2 + api(rax(1), rax(2))**2) )**(-1)

        write(*, '(a,f10.4)') 'max poloidal / max toroidal  ',  ratio

        write(*,'(a,f10.4)') 'Bphi max location', bax(1)*(x2-x1)/(np+1) +x1
        write(*,'(a,f10.4)') 'A max location',   rax(1)*(x2-x1)/(np+1) +x1

      do k2 =2,n_theta-1
        theta = theta + theta2/(n_theta)
        write(18,'(e12.5)') theta
        write(30,'(e12.5)') theta
        write(31,'(e12.5)') theta
        write(32,'(e12.5)') theta
      enddo

C TEST

       x=0
c       DO J =1,NP+2 + NFT        ! radial loop
       DO J =1,NP+2 + NF        ! radial loop
       x = x1+(J-1)*hh
       XR(J) = x
       enddo

       call hunt(xr,NP+2+NF, xbt, jbt)

       call hunt(xr,NP+2+NF, xbo, jbo)

      Do K1=1,n_time
       time = time + t_fin/n_time
       write(18,'(e12.5)') time
       write(30,'(e12.5)') time
       write(31,'(e12.5)') time
       write(32,'(e12.5)') time
       theta=0
       do k2 =2,n_theta-1

       theta = theta + theta2/(n_theta)

         bphi7= BPHI(jbt,k2)*cos(abg(imag)*time)+IPHI(jbt,k2)*sin(abg(imag)*time)

         write(18,'(e12.5)')  bphi7

         dt=theta2/n_theta

         J=jbt

        brs  = BRR(NP,k2)*cos(abg(imag)*time)+BRI(NP,K2)*sin(abg(imag)*time)

        write(32,'(e12.5)') brs !* bphi7
        apt    = APR(J,k2)*cos(abg(imag)*time)+API(J,K2)*sin(abg(imag)*time)
        bpt    = BPHI(J,k2)*cos(abg(imag)*time)+IPHI(J,k2)*sin(abg(imag)*time)
        bptkm  = BPHI(J,k2-1)*cos(abg(imag)*time)+IPHI(J,k2-1)*sin(abg(imag)*time)
        aptkm  = APR(J,k2-1)*cos(abg(imag)*time)+API(J,K2-1)*sin(abg(imag)*time)
        bptkp  = BPHI(J,k2+1)*cos(abg(imag)*time)+IPHI(J,k2+1)*sin(abg(imag)*time)
        aptkp  = APR(J,k2+1)*cos(abg(imag)*time)+API(J,K2+1)*sin(abg(imag)*time)
        bptjm  = BPHI(J-1,k2)*cos(abg(imag)*time)+IPHI(J-1,k2)*sin(abg(imag)*time)
        aptjm  = APR(J-1,k2)*cos(abg(imag)*time)+API(J-1,K2)*sin(abg(imag)*time)
        bptjp  = BPHI(J+1,k2)*cos(abg(imag)*time)+IPHI(J+1,k2)*sin(abg(imag)*time)
        aptjp  = APR(J+1,k2)*cos(abg(imag)*time)+API(J+1,K2)*sin(abg(imag)*time)
c first derivative wrt theta
        dbtheta =(bptkp-bptkm)/dt/2 
        datheta =(aptkp-aptkm)/dt/2
c second derivative wrt theta
        datheta2 = (aptkp-2*apt+aptkm)/dt/dt
        dbtheta2 = (bptkp-2*bpt+bptkm)/dt/dt
c first derivative wrt x
        dax= (apt-aptjm)/hh
        dbx= (bpt-bptjm)/hh
c second derivative wrt x
        dax2=(aptjp-2*apt+aptjm)/hh/hh
c to understand this formula for the cross helicity, just open the mathematica notebook chelicity2.nb
        chel=(2*apt*bpt/sin(theta)**2 + apt*dbtheta*cos(theta)/sin(theta) 
     . +datheta*dbtheta -bpt*datheta2)/xr(J)/xr(J)
     . +(-bpt*dax+apt*dbx)/xr(J)+dax*dbx -bpt*dax2

         write(30,'(e12.5)') chel

        JO = jbo
        apt    = APR(JO,k2)*cos(abg(imag)*time)+API(JO,K2)*sin(abg(imag)*time)
        bpt    = BPHI(JO,k2)*cos(abg(imag)*time)+IPHI(JO,k2)*sin(abg(imag)*time)
        bptkm  = BPHI(JO,k2-1)*cos(abg(imag)*time)+IPHI(JO,k2-1)*sin(abg(imag)*time)
        aptkm  = APR(JO,k2-1)*cos(abg(imag)*time)+API(JO,K2-1)*sin(abg(imag)*time)
        bptkp  = BPHI(JO,k2+1)*cos(abg(imag)*time)+IPHI(JO,k2+1)*sin(abg(imag)*time)
        aptkp  = APR(JO,k2+1)*cos(abg(imag)*time)+API(JO,K2+1)*sin(abg(imag)*time)
        bptjm  = BPHI(JO-1,k2)*cos(abg(imag)*time)+IPHI(JO-1,k2)*sin(abg(imag)*time)
        aptjm  = APR(JO-1,k2)*cos(abg(imag)*time)+API(JO-1,K2)*sin(abg(imag)*time)
        bptjp  = BPHI(JO+1,k2)*cos(abg(imag)*time)+IPHI(JO+1,k2)*sin(abg(imag)*time)
        aptjp  = APR(JO+1,k2)*cos(abg(imag)*time)+API(JO+1,K2)*sin(abg(imag)*time)
c first derivative wrt theta
        dbtheta =(bptkp-bptkm)/dt/2 
        datheta =(aptkp-aptkm)/dt/2
c second derivative wrt theta
        datheta2 = (aptkp-2*apt+aptkm)/dt/dt
        dbtheta2 = (bptkp-2*bpt+bptkm)/dt/dt
c first derivative wrt x
        dax= (apt-aptjm)/hh
        dbx= (bpt-bptjm)/hh
c second derivative wrt x
        dax2=(aptjp-2*apt+aptjm)/hh/hh
c to understand this formula for the cross helicity, just open the mathematica notebook chelicity2.nb
        chel2=(2*apt*bpt/sin(theta)**2 + apt*dbtheta*cos(theta)/sin(theta) 
     . +datheta*dbtheta -bpt*datheta2)/xr(JO)/xr(JO)
     . +(-bpt*dax+apt*dbx)/xr(JO)+dax*dbx -bpt*dax2
         write(31,'(e12.5)')  chel2
       enddo
      enddo

 2010   format(a220,i6.6,a2,i2.2) 
 2011   format(a220,i6.6,a2,i2.2,a2,i2.2)
         
      close(18)
      close(30)
      close(31)
      close(32)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C radial butterfly diagram N and S 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C   normalization
c       rnor=1
c       k2=50
c       DO J =2,NP+2+NF-1       ! radial loop
c      rnor= rnor+ BPHI(J, K2)**2 + IPHI(J,K2)**2
c       enddo
c      print*, rnor
 

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

      write(33,'(6x,I4)') NP+2+NF-2
      write(34,'(6x,I4)') NP+2+NF-2

       x=0
       DO J =2,NP+2+NF-1       ! radial loop
        x = x1+(J-1)*hh
        write(33,'(e12.5)') x
        write(34,'(e12.5)') x
       enddo

        theta=3.1415926/2.0
        theta=0
        do k2 =1,n_theta
        theta = theta + theta2/n_theta
        ang(k2)=theta
        enddo


c  average on theta 
        nag=2


C +/- 30 degrees latitude

       abt=1.047
       abm=2.094


c       abt=.07
c       abm=2.094


       call hunt(ang,n_theta, abt, kb)

       call hunt(ang,n_theta, abm, km)

       dt=theta2/n_theta

       Do K1=1,n_time
       time = time + t_fin/n_time
       write(33,'(e12.5)') time
       write(34,'(e12.5)') time
       
       DO J =3,NP+2+NF-1       ! radial loop

       cja=0
       DO k2=kb-nag,kb+nag  ! average over angle 
        ko = k2
        apt    = APR(J,ko)*cos(abg(imag)*time)+API(J,ko)*sin(abg(imag)*time)
        bpt    = BPHI(J,ko)*cos(abg(imag)*time)+IPHI(J,ko)*sin(abg(imag)*time)
        bptkm  = BPHI(J,ko-1)*cos(abg(imag)*time)+IPHI(J,ko-1)*sin(abg(imag)*time)
        aptkm  = APR(J,ko-1)*cos(abg(imag)*time)+API(J,ko-1)*sin(abg(imag)*time)
        bptkp  = BPHI(J,ko+1)*cos(abg(imag)*time)+IPHI(J,ko+1)*sin(abg(imag)*time)
        aptkp  = APR(J,ko+1)*cos(abg(imag)*time)+API(J,Ko+1)*sin(abg(imag)*time)
        bptjm  = BPHI(J-1,ko)*cos(abg(imag)*time)+IPHI(J-1,ko)*sin(abg(imag)*time)
        aptjm  = APR(J-1,ko)*cos(abg(imag)*time)+API(J-1,ko)*sin(abg(imag)*time)
        aptjm2  = APR(J-2,ko)*cos(abg(imag)*time)+API(J-2,ko)*sin(abg(imag)*time)
        bptjp  = BPHI(J+1,ko)*cos(abg(imag)*time)+IPHI(J+1,ko)*sin(abg(imag)*time)
        aptjp  = APR(J+1,ko)*cos(abg(imag)*time)+API(J+1,ko)*sin(abg(imag)*time)

c first derivative wrt theta
        dbtheta =(bptkp-bptkm)/dt/2 
        datheta =(aptkp-aptkm)/dt/2
c second derivative wrt theta
        datheta2 = (aptkp-2*apt+aptkm)/dt/dt
        dbtheta2 = (bptkp-2*bpt+bptkm)/dt/dt
c first derivative wrt x
        dax= (aptjp-aptjm)/hh/2
        dbx= (bptjp-bptjm)/hh/2

c second derivative wrt x 
 
        if(j.ne.NP+2)then
        dax2=(aptjp-2*apt+aptjm)/hh/hh
        else
c        dax2=0
        dax2=(apt-2*aptjm+aptjm2)/hh/hh
        endif

c        dbx2=(bptjp-2*bpt+bptjm)/hh/hh
c to understand this formula for the cross helicity, just open the mathematica notebook notation.nb

        chelN=(2*apt*bpt/sin(ang(ko))**2 + apt*dbtheta*cos(ang(ko))/sin(ang(ko))
     . +datheta*dbtheta -bpt*datheta2)/xr(J)/xr(J)
     . +(-bpt*dax+apt*dbx)/xr(J)+dax*dbx -bpt*dax2

        cja=cja+chelN   

c         write(33,'(e12.5)')  chelN/rnor
c         write(33,'(e12.5)')  apt*bpt
c         write(19, '(7e14.5)') xr(J), apt, bpt, dax, dbx, dax2, dbx2
       ENDDO

         write(33,'(e12.5)')  cja/rnor/2./nag

       ENDDO


       DO J =3,NP+2+NF-1       ! radial loop
       cja=0
       DO k2=km-nag,km+nag  ! average over angle 
        ko = k2
        apt    = APR(J,ko)*cos(abg(imag)*time)+API(J,ko)*sin(abg(imag)*time)
        bpt    = BPHI(J,ko)*cos(abg(imag)*time)+IPHI(J,ko)*sin(abg(imag)*time)
        bptkm  = BPHI(J,ko-1)*cos(abg(imag)*time)+IPHI(J,ko-1)*sin(abg(imag)*time)
        aptkm  = APR(J,ko-1)*cos(abg(imag)*time)+API(J,ko-1)*sin(abg(imag)*time)
        bptkp  = BPHI(J,ko+1)*cos(abg(imag)*time)+IPHI(J,ko+1)*sin(abg(imag)*time)
        aptkp  = APR(J,ko+1)*cos(abg(imag)*time)+API(J,Ko+1)*sin(abg(imag)*time)
        bptjm  = BPHI(J-1,ko)*cos(abg(imag)*time)+IPHI(J-1,ko)*sin(abg(imag)*time)
        aptjm  = APR(J-1,ko)*cos(abg(imag)*time)+API(J-1,ko)*sin(abg(imag)*time)
        aptjm2  = APR(J-2,ko)*cos(abg(imag)*time)+API(J-2,ko)*sin(abg(imag)*time)
        bptjp  = BPHI(J+1,ko)*cos(abg(imag)*time)+IPHI(J+1,ko)*sin(abg(imag)*time)
        aptjp  = APR(J+1,ko)*cos(abg(imag)*time)+API(J+1,ko)*sin(abg(imag)*time)

c first derivative wrt theta
        dbtheta =(bptkp-bptkm)/dt/2 
        datheta =(aptkp-aptkm)/dt/2
c second derivative wrt theta
        datheta2 = (aptkp-2*apt+aptkm)/dt/dt
        dbtheta2 = (bptkp-2*bpt+bptkm)/dt/dt
c first derivative wrt x
        dax= (aptjp-aptjm)/hh/2
        dbx= (bptjp-bptjm)/hh/2
c second derivative wrt x

        if(j.ne.NP+2)then
        dax2=(aptjp-2*apt+aptjm)/hh/hh
        else
c        dax2=0
        dax2=(apt-2*aptjm+aptjm2)/hh/hh
        endif

c        dax2=(aptjp-2*apt+aptjm)/hh/hh
c to understand this formula for the cross helicity, just open the mathematica notebook notation.nb
        chelS=(2*apt*bpt/sin(ang(ko))**2 + apt*dbtheta*cos(ang(ko))/sin(ang(ko))
     . +datheta*dbtheta -bpt*datheta2)/xr(J)/xr(J)
     . +(-bpt*dax+apt*dbx)/xr(J)+dax*dbx -bpt*dax2

        cja=cja+chelS
       ENDDO


         write(34,'(e12.5)')  cja/rnor/2./nag
c         write(34,'(e12.5)')  apt*bpt
c         write(19, '(7e14.5)') xr(J), apt, bpt, dax, dbx, dax2, dbx2
       ENDDO

c      stop

      ENDDO

      close(33)
      close(34)

      END


