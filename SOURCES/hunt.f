
      SUBROUTINE hunt(xx,n,x,jlo)

      INTEGER jlo,n
      REAL x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd

      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
	jlo=0
	jhi=n+1
	goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1	jhi=jlo+inc
	if(jhi.gt.n)then
	  jhi=n+1
	else if(x.ge.xx(jhi).eqv.ascnd)then
	  jlo=jhi
	  inc=inc+inc
	  goto 1
	endif
      else
	jhi=jlo
2	jlo=jhi-inc
	if(jlo.lt.1)then
	  jlo=0
	else if(x.lt.xx(jlo).eqv.ascnd)then
	  jhi=jlo
	  inc=inc+inc
          goto 2
	endif
      endif
3     if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
	jlo=jm
      else
	jhi=jm
      endif
      goto 3
      END



      FUNCTION factrl(n)
      INTEGER n
      REAL factrl
CU    USES gammln
      INTEGER j,ntop
      REAL a(33),gammln
      SAVE ntop,a
      DATA ntop,a(1)/0,1./
      if (n.lt.0) then
        pause 'negative factorial in factrl'
      else if (n.le.ntop) then
        factrl=a(n+1)
      else if (n.le.32) then
        do 11 j=ntop+1,n
          a(j+1)=j*a(j)
11      continue
        ntop=n
        factrl=a(n+1)
      else
        factrl=exp(gammln(n+1.))
      endif
      return
      END


      FUNCTION gammln(xx)
      REAL gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
