c--------------------------------------------------------
      subroutine lepeos(t,mu,den,eps,pre,ent)
c--------------------------------------------------------
      implicit none
      real*8 t,mu,den,eps,pre,ent,pi2
      parameter (pi2=9.8696044)
 
      eps = (mu**4 + 2.d0*pi2*(t*mu)**2 + 
     &         7.d0/15.d0*pi2*pi2*t**4)/(4.d0*pi2)
      pre = eps/3.d0
      ent = t*(mu*mu + 7.d0/15.d0*pi2*t*t)/3.d0
      den = mu*(mu*mu + pi2*t*t)/(3.d0*pi2) 

      return
      end


c--------------------------------------------------------
      subroutine JELdpe(thc,mhc,muhc,den,pres,enrd,xns,entr)
c--------------------------------------------------------
      implicit none
      real*8 thc,mhc,muhc,den,pres,enrd,xns,entr
      real*8 f, muot, g, opf, opg, sqopf, sqopfa, FMUOT
      real*8 pred, sump, pi2, dd, a
      real*8 dsumdf, dsumdg, dpdf, dpdg
      integer mm,nn
      parameter (mm=2,nn=2,pi2=9.8696044,a=0.442)
      real*8 pmn(0:mm,0:nn)
      data pmn / 5.2381, 11.1751, 5.918,
     &          12.4991, 25.3687, 12.4945,
     &          8.36058, 15.6453, 6.8253/ 
      save pmn

      muot = (muhc-mhc)/thc

      if (muot.gt.-50.d0) then

        f = a*FMUOT(muot)
        g = (thc/mhc)*dsqrt(1.d0+f)
        opf = 1.d0 + f
        opg = 1.d0 + g
        sqopf = dsqrt(opf)
        sqopfa = dsqrt(1.d0+f/a)

        pred = (f*g*(g*opg)**1.5)/(opf**(mm+1)*opg**(nn))
 
        sump =  pmn(0,0) + pmn(0,1)*g + pmn(0,2)*g*g
     &       + f*(pmn(1,0) + pmn(1,1)*g + pmn(1,2)*g*g)
     &       + f*f*(pmn(2,0) + pmn(2,1)*g + pmn(2,2)*g*g)

        dsumdf = (pmn(1,0) + pmn(1,1)*g + pmn(1,2)*g*g)
     &       + 2.d0*f*(pmn(2,0) + pmn(2,1)*g + pmn(2,2)*g*g)

        dsumdg = (pmn(0,1) + pmn(1,1)*f + pmn(2,1)*f*f)
     &       + 2.d0*g*(pmn(0,2) + pmn(1,2)*f + pmn(2,2)*f*f)

        dpdf = (pred/f-(mm+1)*pred/opf)*sump + pred*dsumdf
        dpdg = (2.5*pred/g-(float(nn)-1.5)*pred/opg)*sump + pred*dsumdg

        dd = mhc**3/pi2 

        pres = mhc*dd*pred*sump 
        den  = dd*(dpdf+g*dpdg/(2.0*opf))*(f*sqopf/(g*sqopfa))
        entr = -muot*den + dd*sqopf*dpdg 
        enrd = thc*entr - pres + muhc*den
        xns = (enrd - 3.d0*pres)/mhc

      else
      
        pres = 0.
        den  = 0.
        entr = 0.
        enrd = 0.
        xns  = 0.
       
      endif


      return 
      end



c--------------------------------------------------
      FUNCTION FMUOT(x)
c--------------------------------------------------
c   This function solves the implicit equation
c  x = 2*sqrt(1+f)+log( (sqrt(1+f)-1) / (sqrt(1+f)+1) )
c
c     INPUT : x      ;    OUTPUT: f
c--------------------------------------------------
      implicit none
      real*8 x,f,logf,fmuot
      integer i,iload,nt,j
      parameter (nt=400)
      real*8 xa(nt),lfa(nt),xmin,dx,a,b

      save iload,xa,lfa,xmin,dx
      data iload/0/

      if (iload.eq.0) then
        open(unit=20,file='FMUOT.TAB')
        do i=1,nt
          read(20,*) xa(i),lfa(i)
        enddo
        close(20)
        xmin=xa(1)
        dx=xa(2)-xa(1)
        iload=1
      endif

      if (x.lt.-9.9) then
        f = 1.d0/(dsinh(1.d0-0.5d0*x))**2
      elseif (x.gt.9.9) then
        f = 1.d0+0.25d0*x*x
      else
        j = int((x-xmin)/dx)+1
        if(j.lt.1.or.j.gt.399) then
         print*,j,x,dx
         pause
        endif
        a = (x-xa(j))/dx
        b = 1.d0 - a
        logf = lfa(j)*b+lfa(j+1)*a
        f=10.d0**logf
      endif

      fmuot=f

      return
      end

c--------------------------------------------------
      FUNCTION FHBOSE(psi)
c--------------------------------------------------
c   This function solves the implicit equation
c              psi = f/(1+f)-log(1+f)
c
c     INPUT : psi      ;    OUTPUT: f
c--------------------------------------------------
      implicit none
      real*8 psi,fun,dfun,fhbose,f

      if (psi.gt.0.) then
        print*,'FHBOSE: PSI=',psi
        pause
      elseif (psi.gt.-1.0) then
        f=sqrt(-2.0*psi)
      else
        f=exp(1.0-psi)
      endif

 10   fun = f/(1.d0+f)-dlog(1.d0+f)-psi
      if (dabs(fun).gt.1.d-10) then
        dfun = -f/(1.0+f)**2
        f=f-fun/dfun
        goto 10
      endif

      fhbose=f

      return
      end


c--------------------------------------------------
      SUBROUTINE BOSE(T,PSI,MK,PRES,EPS,RHO,AKTH)
c--------------------------------------------------
      IMPLICIT NONE
      REAL*8 COEFF(5,5),A, PI2
      INTEGER M,N, I, J, K
      PARAMETER (A=1.029, M=4, N=5, PI2=9.8696044)
      REAL*8 T,PSI,MK,PRES,EPS,RHO,AKTH,SA,H,FHBOSE
      REAL*8 FRONT,FRONTN,FAC,XX,YY,ZZ,ZZN

      DATA (COEFF(1,J),J=1,5)/1.68134,6.8507,10.8537,7.81843,2.16461/
      DATA (COEFF(2,J),J=1,5)/6.72536,27.4028,43.4148,31.2737,8.65844/
      DATA (COEFF(3,J),J=1,5)/8.49651,35.6058,57.7134,42.3593,11.8199/
      DATA (COEFF(4,J),J=1,5)/3.45614,15.1152,25.5254,19.2745,5.51757/
      DATA (COEFF(5,J),J=1,5)/0.,0.,0.,0.,0./
      SAVE COEFF

      SA=DSQRT(A)
      H=SA*FHBOSE(PSI)

      PRES=0.D0
      EPS=0.D0
      RHO=0.D0

      FAC = (T*(T+1.D0))**1.5D0/((1.D0+H)**M*(1.D0+T)**(N-1))/(2.D0*PI2)
      FRONT = T*FAC
      FRONTN = FAC*(H+SA)**2/((1.D0+H))

      DO I=1,M
        XX=0.D0
        YY=0.D0
        IF(H.EQ.0.D0)THEN
          IF(I.EQ.1)XX=1.D0
          IF(I.EQ.2)YY=1.D0
        ELSE
          ZZ=57.D0 - DFLOAT(M-I)*DLOG(H)
          ZZN=ZZ-DLOG(H)
          IF(ZZ.GT.0.D0)XX=H**(I-1)
          IF(ZZN.GT.0.D0)YY=H**(I-2)
        ENDIF

        DO K=1,N
          PRES=PRES+COEFF(I,K)*XX*T**(K-1)
          EPS=EPS+COEFF(I,K)*XX*T**(K-1)*
     &     (0.5D0+DFLOAT(K)+(2.5D0-DFLOAT(N))*T/(1.D0+T))
          RHO=RHO+YY*(T**(K-1))*(-DFLOAT(I)*COEFF(I+1,K)
     &     +DFLOAT(M+1-I)*COEFF(I,K))

        ENDDO
      ENDDO

      AKTH = MK**2*FRONT*(EPS-3.D0*PRES)
      PRES = MK**4*FRONT*PRES
      EPS  = MK**4*(FRONT*EPS+FRONTN*RHO)
      RHO  = MK**3*FRONTN*RHO

      RETURN
      END


c--------------------------------------------------
      FUNCTION FUN2(psi)
c--------------------------------------------------
c   This function solves the implicit equation
c              psi = 1-y+log(y)
c
c     INPUT : psi      ;    OUTPUT: y
c--------------------------------------------------
      implicit none
      real*8 psi,fun,dfun,fun2,y

      if (psi.gt.0.) then
        print*,'FUN2: PSI=',psi
        pause
      elseif (psi.lt.-1.0) then
        y=dexp(psi-1.d0)
      else
        y=1.0
      endif

 10   fun = 1.d0-y+dlog(y)-psi
      if (dabs(fun).gt.1.d-10) then
        dfun = -1.d0+1.d0/y
        y=y-fun/dfun
        goto 10
      endif

      fun2=y

      return
      end

