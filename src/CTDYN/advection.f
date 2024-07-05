C     LINEAR DYNAMO CODE CTDYN 
C                     DOUBLE PRECISION/CHANGE INTRINSIC es: erf -> derf 
C
C     one cell   s0 = -1/2    s2 =  3/2     s4=0    s6=0  
C     two cells  s0 =  3/8    s2 = -15/4    s4=35/8 
C
C
C   the induction equation is written Dt B = rot ( v X B "+" alpha B) -rot eta rot B
C
C   therefore the turbulent ELF is E = + alpha B
C   
C   this alpha is "minus" the Raedler AN 307 89 1986 paper (R86)
C   
C   The B field is written  as  
C
C   B \propto exp (  \nu \tau)
C   
C   and in the R86 paper instead we have 
C   
C   B \propto exp (im \phi + (\lambda - i \Omega) \tau )      
C
C   therefore Im(nu) =  - \Omega
C
C  As a test example you have for  n = 14 m = 30
C
C    A0 Comega = 20      -> Calpha = 4.51   \nu = 0
C    A1 Comega = 20      -> Calpha = 5.01   \nu = 7.43 = - \Omega R86
C    

      PROGRAM ADVECTION
      IMPLICIT REAL(A-H,O-Z)

      include 'cio'

      REAL imag
      COMMON/part/vtu,rt,imag,CO,c_u,beta,ffree,betb,etep,etet,xbt,xbo
      
      REAL CA(10,10,4)
      REAL REG(10),IEG(10)
      REAL eep
      COMMON/eira/REG,IEG,eep
      INTEGER II,IT
      COMMON/ipar/II,IT

      REAL xa1,xa2,xa3,xb,xda1,xda2
      COMMON/apar/xa1,xa2,xa3,xb,xda1,xda2

      REAL edr,xe1,xde1
      COMMON/epar/edr,xe1,xde1,hd

      REAL x_in, bct, c3, mmm,SR,rotp,gd, aqu, flg
      COMMON/ppar/x_in,bct,c3,mmm,SR,rotp, gd,aqu, flg

      REAL dd1,rc1,rc2,oco
      COMMON/dpar/dd1,rc1,rc2,oco
      
      REAL s0,s2,s4,s6, A2P, A4P,xm
      COMMON/PSI/s0,s2,s4,s6,A2P,A4P,xm

      character(len=128)    ::  dir
      CHARACTER*8 ANS1,ANS2,ANS3,ANS4
      COMMON/var3/ANS1,ANS2,ANS3,ANS4, dir

      CHARACTER*2 JOBVR,JOBVL
      COMMON/lap/JOBVR,JOBVL
      
      CHARACTER*30 inp
      CHARACTER*52 mach, ddt
      CHARACTER*43 version, ver

      character(len=5)     ::  char_m  ! converti m in char
      character(len=128)    ::  char_fm ! file input_Am
      character(len=128)    ::  char_cm ! file crtal_Am
      character(len=128)    ::  char_cf ! file crtaf_Am
      character(len=128)    ::  fom     ! crea un format dinamico 
      character(len=5)      ::  qq      ! _Am
      EXTERNAL version
      INTEGER mm

      integer               ::  iome,nso
      integer               ::  ialo,nsa
      integer               ::  ires,nsr

      COMMON/parker/gam,zeta_r,ratio

      include 'cint'

      CALL getarg(1,inp)
      OPEN(31,file=inp,status='old')
      READ(31,NAMP)
      CLOSE(31)
       !converti mm in char e crea stringa _Am o _Sm
       mm=int(mmm)
       write(char_m,'(I2)')  mm
       if(mod(mm,2).eq.0)then 
          if(ans2.eq.'D')then
            qq=trim('_A'//adjustl(char_m))
          else
            qq=trim('_S'//adjustl(char_m))
          endif
        else if(mod(mm,2).eq.1)then
          if(ans2.eq.'D')then
            qq=trim('_S'//adjustl(char_m))
          else
            qq=trim('_A'//adjustl(char_m))
          endif
        endif

       char_fm=trim(trim(dir)//'/input'//qq)
       char_cm=trim(trim(dir)//'/crtal'//qq)
       char_cf=trim(trim(dir)//'/crtaf'//qq)

       open(34,status='unknown',file=char_fm)
       write(34,NAMP)
       close(34) 
       open(35,status='unknown',file=char_cm)

      ! MAIN ff-beta parameter
      beta=beta_i
      eep=-99
      II=0

      JOBVR='N'
      ANS3='CR'

      WRITE(35,'(a,2I4)')  '#N c_alpha,c_omega,R,beta,etor,epol,z=', np, na+1 

      DO 333 iome = 0,nso
      if(CM_I.ne.0)then
      astep=(abs(CM_F/CM_I))**(1./nso) 
      CO = astep**iome*CM_I

      else
      CO = CM_I + iome*(CM_F-CM_I)/(nso+1) 
      endif

      c_u=RM_I+RM_F*co**xm

      II=II+1    
      CALL ZBR(AL_I,AL_F,accu,critical) 
      WRITE(35,'(I4,16e12.4)') II, critical,CO,C_U,beta,etep,etet,zeta_r,(REG(I),ABS(IEG(I)),I=1,1), eep
      if(ANS4.eq.'V')then
      JOBVR='V'
      CALL DYNAMO(critical,rate)
      JOBVR='N'
      ENDIF
      
 333  continue

      close(35)
      
      end

