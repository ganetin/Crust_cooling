C cgs units for P and rho
C nb in fm^-3
CCC Ep=938.272046,En=939.565379
      IMPLICIT NONE
 3334 FORMAT (1X,21E16.8)
      real*8 dattab(17)
      real*8 nb(500),P(500),rho(500),Ye(500),Ym(500)
      real*8 Yp(500),Yn(500),Yl(500),Ysm(500),Ys0(500)
      real*8 Ysp(500),mp(500),mn(500),ml(500),msm(500)
      real*8 ms0(500),msp(500)
      real*8 JUNK,ZERO,SEVEN
      integer*8 K,KDAT,I

      ZERO=0.d0
      JUNK=42.d0
      SEVEN=0.7d0

      K=0
      OPEN (3,file='BHF_UVIX_EOS_Cat_Fe.dat',status='old')
!crho_t(i),xnut,nbar_t(i),
!     2         yelect_t(i),ymuon_t(i),yneutr_t(i),yprot_t(i),
!     3         ylambda_t(i),ysminus_t(i),yszero_t(i),ysplus_t(i)
c*** 7 lines of comments in input file
      DO 90 I=1,7
 90      READ (3,*)
 7    K=K+1
 71   READ (3,*,END=8,ERR=71) (dattab(i),i=1,17)
      rho(K)=dattab(1)
      P(K)=dattab(2)
      nb(K)=dattab(3)
      Ye(K)=dattab(4)
      Ym(K)=dattab(5)
      Yn(K)=dattab(6)
      Yp(K)=dattab(7)
      Yl(K)=dattab(8)
      Ysm(K)=dattab(9)
      Ys0(K)=dattab(10)
      Ysp(K)=dattab(11)
      mp(K)=dattab(12)
      mn(K)=dattab(13)
      ml(K)=dattab(14)
      msm(K)=dattab(15)
      ms0(K)=dattab(16)
      msp(K)=dattab(17)
      write(*,*)K
      GOTO 7
 8    CLOSE(3)
      KDAT=K-1
c NSCool Input Core
!   Rho	   Press	    nbar	    Ye	    Ymu	    Yn	    Yp	    Yla	    Ysm	    Ys0	    Ysp	Yxim	Yxi0	   mstp	   mstn	   mstla	   mstsm	   msts0	   mstsp	mstxim	mstxi0
!  g/cm3	 dyne/cm2	   #/fm3	  [A_cell]	   [A_ion]	    [Z]	

      DO K=1,KDAT
         write(90,3334)RHO(K),P(K),nb(K),Ye(K),Ym(K),Yn(K),Yp(K),
     1   Yl(K),Ysm(K),Ys0(K),Ysp(K),ZERO,ZERO,
     1   mp(K),mn(K),ml(K),msm(K),ms0(K),msp(K),ZERO,ZERO
      ENDDO 
!c Input
!      integer KEOS
!      real*8 nb
!c Output
!      real*8 P,RHO,Ye,Ymu,Yp,Yn,meffp,meffn,Zcell,Acell,Aion,Zclust
!      real*8 RLG,PLG,JUNK,ZERO
!      real*8 nbcore,nbcrust
!      real*8 XND,XNC            
!      integer*8 I
!       ZERO=0.D0
!       KEOS=21
!       nb=0.14

!       if (KEOS.eq.19) nbcore=.0885
!       if (KEOS.eq.20) nbcore=.0854
!       if (KEOS.eq.21) nbcore=.0809

!       DO I=1,300
!       nb=10**(-8+(I-1)/(300.-1.)*(0.+8.)) !nbcore+(I-1)/(100.-1.)*(2.-nbcore)
!       call BSkRofN(KEOS,nb,RHO)
!       RLG=log10(RHO)
!       write(*,*) nb,rlg
!       call BSKfit(KEOS,RLG,PLG,JUNK,JUNK)
!       P=10**PLG



!       if (nb.ge.nbcore) then 
!         call FRACORE(KEOS,nb,Ye,Ymu)
!         Yp=Ye+Ymu
!         Yn=1.-Yp
!         call EFFMASS(KEOS,nb,Yp,meffp,meffn)
!         write(90,*)RHO,P,nb,Ye,Ymu,Yn,Yp,ZERO,ZERO
!     *   ,ZERO,ZERO,meffp,meffn,ZERO,ZERO,ZERO,ZERO
!       endif   

!       if (nb.lt.nbcore) then 
!         call INCRUST(KEOS,nb,XND,XNC,Aion,Acell,Zcell,Zclust,
!     *  JUNK,JUNK,JUNK,JUNK,JUNK,JUNK,JUNK,JUNK)
!         write(90,*)RHO,P,nb, Acell,Aion,Zcell,Zclust
!       endif  


 

      END


