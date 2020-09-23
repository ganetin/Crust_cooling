C***  PROGRAM NSPRHO (INPUT=1, OUTPUT=5)
C
C version A - old types of data file (EOS) - k,n,ro,P in CGS units
C
C     PURPOSE:
C
C        TO CALCULATE A MODEL OF A NEUTRON STAR.
C
C        TO OBTAIN, FOR A GIVEN CENTRAL PRESSURE, THE FOLLOWING
C        PARAMETERS WHICH CHARACTERIZE THE NEUTRON STAR MODEL:
C        (1)  THE RADIUS OF CONFIGURATION - R,
C        (2)  THE (GRAVITATIONAL) MASS - M,
C        (3)  THE TOTAL NUMBER OF BARYONS - A,
C        (4)  THE PROPER MASS - MPROP,
C        (5)  THE SURFACE REDSHIFT - Z,
C        (6)  THE MOMENT OF INERTIA FOR SLOW AND RIGID ROTATION - I.
C
C     REMARKS:
C
C        IN ORDER TO OBTAIN THESE PARAMETERS WITH
C        A HIGH PRECISION, 'NSTARP' USES ENERGY PER
C
C        WITHIN THIS PROGRAM PRESSURE IS EXPRESSED IN
C        UNITS OF MATTER DENSITY (G/CM**3), THAT IS
C        P=PRESSURE/C0**2, WHERE C0 IS THE SPEED OF LIGHT.
C
C        THE RELATIVE ERRORS IN THE PARAMETERS CHARACTERIZING
C        THE MODEL DO NOT EXCEED THE VALUE  K*TOL, WHERE
C        K IS A NUMBER OF INTEGRATION STEPS AND TOL IS
C        SPECIFIED BY THE USER (TOL IS AN INPUT PARAMETER).
C        INDEED, THEY ARE USUALLY MUCH SMALLER THAN  K*TOL,
C        BEING OF THE ORDER OF K*TOL/100 (SUCH A VALUE IS
C        PRINTED).
C
C     INPUT/OUTPUT UNITS:
C
C        file1 = IS USED TO READ A FILE WITH THE (TABULATED)
C            EQUATION OF STATE.
C            THIS ARRAY SHOULD BE PREPARED AS FOLLOWS:
C            - THE FIRST AND SECOND LINE IS A DUMMY TEXT DESCRIBING THE DATA,
C            - NEXT LINES SHOULD CONTAIN THE TABULATED EQUATION
C            OF STATE, EACH LINE CONSISTING OF:
C                       NDAT, RHODAT, PDAT,
C            WHERE KK IS A DUMMY INTEGER, NDAT IS A BARYON
C            DENSITY, EPSDAT IS ENERGY PER ONE BARYON AND
C            PDAT IS CORRESPONDING PRESSURE (in cgs units !).
C
C        5 = IS USED TO READ THE INPUT DATA (PCENTR, TOL, JPRINT)
C            AND TO PRINT THE RESULTS OF INTEGRATION, I.E.
C            PARAMETERS DESCRIBING THE NEUTRON STAR MODEL.
C
C     INPUT PARAMETERS:
C
C        PC      = PRESSURE (DIVIDED BY C0**2) IN THE CENTER OF THE STAR.
C        TOL     = RELATIVE TOLERANCE OF THE INTEGRATED FUNCTIONS
C                  AT ONE INTEGRATION STEP.
C
C     VERSION: 15-MAY-82
C
C     WARNING:
C
C        VARIABLES BEGINNING WITH  I, L, M, N  ARE DECLARED AS
C        IMPLICIT DOUBLE PRECISION, INSTEAD OF BEING INTEGERES.
C
C--------------------------------------------------------------------
C
C      W NSP ZMIENNA NIEZALEZNA JEST CISNIENIE - P
c      file1 - dane rownanie stanu (opis powyzej)
c      file2 - wyniki calkowite (nie w tabelce)
C      NSP3 LICZY DLA ZESTAWU DANYCH PC W ZBIORZE file3
C      WYNIKI ZAPISUJE TEZ W POSTACI TABELKI:
C       PC,ROC,NC,R,M,I,Z,RCORE,MCORE
C      W ZBIORZE file4
c      do IBM/PC
C
C            WERSJA 21.01.1998
C
C--------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,I,L,M,N,O-Z)
      double precision lnbmin,lnbmax,dlnxnb,lnxnb
      DIMENSION Y(10),PDAT(9999),RHODAT(9999),NDAT(9999),TEXT(140),
     1   F1(10),F2(10),F3(10),F4(10),F5(10),F6(10),F7(10),YDUM(10)
      character*30 file1,file2,file3,file4,file7
      EXTERNAL DIF9

      COMMON /AAA/ JFLAG
      COMMON /BBB/ NDAT,RHODAT,PDAT,KDAT
      COMMON /CCC/ PCORE,ROCORE
      COMMON /DDD/ C0,MNEUTR

      DATA ZZ1,ZZ2 /1H-,1H*/
c      DATA GC2 /7.4250581D-29/
      DATA GC2 /7.4256039440906D-29/
c      DATA MSOL /1.989D+33/
      DATA MSOL /1.9885D+33/
      DATA RKM /1.D+05/

      MNEUTR=1.67492
c      C0=2.9979D0
      C0=2.99792458D0
      PI=DACOS(-1.D0)
c...   max ilosc krokow w kutcie
      kmax=20000
c      write (*,*) 'zbior eos ? '
c      read (*,3333) file1
c      write (*,*) 'zbior danych cisnien ? '
c      read (*,3333) file3
c      write (*,*) 'zbior wynikowy (tabelka) ? '
c      read (*,3333) file4
c      write (*,*) 'zbior wynikowy (wszystko) ? '
c      read (*,3333) file2

c    name of the datafile with datafilenames

      file7='nsprho.dat'
c

 3333 format (a30)

      open (7,file=file7,status='old')
      read (7,3333) file1
c      read (7,3333) file3
      read (7,3333) file4
      read (7,3333) file2
c read min and max pressure (log) and step in logs
      read (7,*) lnbmin,lnbmax,dlnxnb
c
      READ (7,*) TOL,JPRINT
      close (7)

      open (1,file=file1,status='old')
c      open (3,file=file3,status='old')
      open (2,file=file2,status='new')
      open (4,file=file4,status='new')

c      WRITE (*,6000)
C****************************************************************
c      READ (1,6001) (TEXT(J),J=1,8)
c	read (1,*)
c	read (1,*)
c      write (*,*) ' krok 1'
c      READ (1,*)
c      READ (1,*)

c      write (*,*) ' krok 1'
c      READ (1,6006) PCORE,ROCORE
c      READ (1,*) PCORE,ROCORE
c      write (*,*) ' krok 2'
c      write (90,*) PCORE, ROCORE
c
c      ROCORE    = ROCORE*1.7827E12*197.33
c      PCORE     = PCORE*1.6022E33*197.33
c

      PCORE=0.
      ROCORE=0.
      PCORE=PCORE/C0/C0*1.E-20


      K=0
    7 K=K+1
c
c  NOWA KOLEJNOSC KOLUMN
c
c      READ (1,*,END=8) K1,NDAT(K),RHODAT(K),PDAT(K)
      READ (1,*,END=8) RHODAT(K),PDAT(K),NDAT(K)
      WRITE (*,*) RHODAT(K),PDAT(K),NDAT(K)
 
c      READ (1,*,END=8)  NDAT(K),RHODAT(K),PDAT(K)
c      READ (1,*,END=8)  NDAT(K),PDAT(K),RHODAT(K)

      K1=K
      PDAT(K)=PDAT(K)/C0/C0*1.E-20
cc
c      RHODAT(K) = RHODAT(K)*1.7827E12*197.33
c      PDAT(K)   = PDAT(K)*1.6022E33*197.33

cc
      GOTO 7
    8 CLOSE (1)
      KDAT=K-1
c      WRITE (4,6005) (TEXT(J),J=1,20)
c      WRITE (2,6005) (TEXT(J),J=1,20)
c     if JPRINT=1 write whole configuration to file2
      IF (JPRINT.EQ.1) WRITE (2,6019)
 6019 FORMAT ( 3X,'K',6X,'DP/PC',6X,'R[KM]',10X,'P',9X,'RO',9X,'N',9X,
     X'M/MSOL')
      WRITE (4,7001)
      WRITE (*,7001)
c      close (4)
c      close (2)

      lnxnb=lnbmin

 50   PC=10.d0**lnxnb

c 20   READ (3,*,END=25) PC
      PC=PC/C0/C0*1.E-20
c      WRITE (*,6007) PC
      CALL STATE (NC,ROC,PC)
      P=PC
      R=1.D0
      DP=-1.D0
      M=0.D0
      PHI=0.D0
      A=0.D0
      MCORE1=0.D0
      RCORE1=0.D0
      ACORE1=0.D0
      LCORE1=0.D0
      MCORE2=0.D0
      RCORE2=0.D0
      ACORE2=0.D0
      OMEGA=1.D0
      L=0.D0
      MPROP=0.D0
      POLD=PC
      ROOLD=ROC
      JFLAG=0
      JN=9
      Y(1)=R
      Y(2)=M
      Y(3)=PHI
      Y(4)=A
      Y(5)=OMEGA
      Y(6)=L
      Y(7)=MPROP
      Y(8)=1.D0
      Y(9)=1.D0
      K=0
C***
C* 11 ROLD=R
C*    REND=ROLD+1.D+04
C***
   12 ROLD=R
      DPOLD=DP
      POLD=P
      ROOLD=RO
C*    IF (DR.GT.REND-R)   DR=REND-R
      CALL RK78 (JFL,JN,P,DP,Y,TOL,YDUM,F1,F2,F3,F4,F5,F6,F7,DIF9)
      CALL STATE (NNEW,RONEW,P)
      IF (JFLAG.EQ.1)   GOTO 13
      K=K+1
      N=NNEW
      RO=RONEW
      R=Y(1)
      M=Y(2)
      PHI=Y(3)
      A=Y(4)
      OMEGA=Y(5)
      L=Y(6)
      MPROP=Y(7)
      CHP=(P+RO)/N
      CONST=CHP*DEXP(PHI)
c
c     if JPRINT=1 write whole configuration to file2

      IF (JPRINT.EQ.1)   WRITE (2,6009) K,DP/PC,R/RKM,P,RO,N,M/MSOL
      IF (JPRINT.EQ.1)   WRITE (*,6009) K,DP/PC,R/RKM,P,RO,N,M/MSOL
c
      IF (JPRINT.EQ.2)   WRITE (*,6017) K,DP/PC,R/RKM,P,CHP,CONST
      CALL CORE1 (POLD,P,ROOLD,RO,R,M,A,L,RCORE1,MCORE1,ACORE1,LCORE1,
     1                                    RCORE2,MCORE2,ACORE2,LCORE2)
c      if (K.gt.kmax) go to 20
      if (K.gt.kmax) go to 50
      GOTO 12
C***
C*    IF (R.LT.REND)   GOTO 12
C*    DP=DPOLD
C*    WRITE (*,6018) K,R/RKM,N,RO,P*(C0*1.D+10)**2,M/MSOL,A*1.D-19
C*    GOTO 11
C***
   13 R=ROLD
c
c   if JPRINT=1 do not write parameter of the star and exit
c    NO LOOP
      IF (JPRINT.EQ.1) GOTO 25

C*    WRITE (*,6018) K,R/RKM,N,RO,P*(C0*1.D+10)**2,M/MSOL,A*1.D-19
      DELTA=1.D0-2.D0*GC2*M/R
      Z=1.D0/DSQRT(DELTA)-1.D0
      OMEGA=OMEGA*DSQRT(DELTA)*DEXP(-PHI)
      OMEGAI=OMEGA+2.D0*L/R**3
      I=L/OMEGAI
      ICORE1=LCORE1/OMEGAI
      ICORE2=LCORE2/OMEGAI
      Q1=ICORE1/I
      Q2=ICORE2/I
c      call opn2 (file2)
c      call opn4 (file4)
c      WRITE (2,6013) (ZZ2,J=1,70)
c      WRITE (2,6013) (ZZ1,J=1,70)
c      WRITE (2,*) '     '
c      WRITE (2,*) '.................NSTAR P3.................'
c      WRITE (2,*) '     '
c      WRITE (2,*) '-------THE EQUATION OF STATE:'
c      WRITE (2,6005) (TEXT(J),J=1,10)
c      WRITE (2,*) '-------THE CENTER:'
      WRITE (2,6003) PC,ROC,NC
c      WRITE (2,*) '-------THE CORE(S):'
c      WRITE (2,6010) PCORE,ROCORE,RCORE1/RKM,RCORE2/RKM,
c     1 MCORE1/MSOL,MCORE2/MSOL,ACORE1*1.D-19,ACORE2*1.D-19,
c     2 ICORE1*1.D-30/GC2*1.D-15,ICORE2*1.D-30/GC2*1.D-15,Q1,Q2
c      WRITE (2,*) '-------THE SURFACE:'
      WRITE (2,6004) R,R/RKM,M,M/MSOL,A,A*1.D-19,I,
     2 I*1.D-30/GC2*1.D-15,DELTA,Z,MPROP,MPROP/MSOL
c      WRITE (2,*) '-------THE ACCURACY:'
      WRITE (2,6014) K,TOL,TOL*K/100.D0
c      WRITE (2,6013) (ZZ1,J=1,70)
c      WRITE (2,6013) (ZZ2,J=1,70)
c      WRITE (4,7000) PC*C0*C0*1.E20,ROC,NC,R/RKM,M/MSOL,A*1.D-19,
c     #MPROP/MSOL,I*1.D-30/GC2*1.D-15,Z
c      WRITE (4,7002) PC*C0*C0*1.E20,ROC,NC,R/RKM,M/MSOL,A*1.D-19,
c     #MPROP/MSOL,I*1.D-30/GC2*1.D-15,Z,RCORE1/R,MCORE1/M,
c     #RCORE2/R,MCORE2/M
      WRITE (4,7003) PC*C0*C0*1.E20,ROC,NC,R/RKM,M/MSOL,A*1.D-19,
     #MPROP/MSOL,I*1.D-30/GC2*1.D-15,Z,
     #RCORE1/RKM,MCORE1/MSOL
c     #1.-RCORE1/R,1.-MCORE1/M,
c     #1.-RCORE2/R,1.-MCORE2/M
      WRITE (*,7000) PC*C0*C0*1.E20,ROC,NC,R/RKM,M/MSOL,A*1.D-19,
     #MPROP/MSOL,I*1.D-30/GC2*1.D-15,Z
c     #,Z,RCORE1/R,MCORE1/M
c      close (2)
c      close (4)
      lnxnb=lnxnb+dlnxnb
      if (lnxnb.le.lnbmax) go to 50
c      GO TO 20
 25   CLOSE (3)
      close (2)
      close (4)
      STOP
 6000 FORMAT (1X/1X,'TYPE  TOL, JPRINT >>')
 6001 FORMAT (20A7)
 6002 FORMAT (E20.10,I10)
 6003 FORMAT (1X/1X,'PC     =',E20.10/
     1           1X,'ROC    =',E20.10/
     2           1X,'NC     =',E20.10/)
 6004 FORMAT (1X/1X,'R      =',E20.10,10X,F15.9,5X,'KM'/
     1           1X,'M      =',E20.10,10X,F15.9,5X,'MSOL'/
     2           1X,'A      =',E20.10,10X,F15.9,5X,'*1.E+58'/
     3           1X,'I      =',E20.10,10X,F15.9,5X,'*1.E+45'/
     4           1X,'DELTA  =',E20.10/
     1           1X,'Z      =',E20.10/
     2           1X,'MPROP  =',E20.10,10X,F15.9,5X,'MSOL'/)
 6005 FORMAT (1X,20A7/)
 6006 FORMAT (2E20.10)
 6007 FORMAT ('  PC =',G15.8)
 6009 FORMAT (1X,I4,3X,2(E10.4,1X),2X,3(E10.4,1X),2X,E10.4)
 6010 FORMAT (1X/1X,'PCORE  =',E19.9,5X,'ROCORE =',E19.9//
     1           1X,'RCORE  =',F15.9,9X,'RCORE  =',F15.9,5X,'KM'/
     2           1X,'MCORE  =',F15.9,9X,'MCORE  =',F15.9,5X,'MSOL'/
     3           1X,'ACORE  =',F15.9,9X,'ACORE  =',F15.9,5X,'*1.E+58'/
     4           1X,'ICORE  =',F15.9,9X,'ICORE  =',F15.9,5X,'*1.E+45'/
     5           1X,'Q      =',F15.9,9X,'Q      =',F15.9/)
 6011 FORMAT (1X,I3,3E20.12)
 6013 FORMAT (1X,70A1)
 6014 FORMAT (1X/1X,'K    =',3X,I5,10X,'TOL =',E10.3//
     1           1X,'CONSERVATIVE ESTIMATE OF REL.ERR.:  ',E10.3/)
 6017 FORMAT (1X,I4,2E11.4,3X,2E16.8,E18.10)
 6018 FORMAT (1X,I4,2X,F6.3,2X,3E12.4,2X,2E12.4)
 7000 FORMAT (1X,9G13.6)
 7002 FORMAT (1X,13G13.6)
 7003 FORMAT (1X,11G17.8)
 7001 FORMAT ('#',6X,'PC',12X,'ROC',10X,'NC',9X,'R(KM)',7X,'M/MSOL',7X,
     #'A(E58)',5x,'MPROP/MSOL',5x,'I(E45)',8X,'Z')
      END
C--------------------------------------------------------------------
      SUBROUTINE CORE1 (POLD,P,ROOLD,RO,R,M,A,L,RCORE1,MCORE1,ACORE1,
     1                  LCORE1,RCORE2,MCORE2,ACORE2,LCORE2)
      IMPLICIT DOUBLE PRECISION (A-H,I,L,M,N,O-Z)
      COMMON /CCC/ PCORE,ROCORE
      IF ((POLD.LT.PCORE).OR.(P.GT.PCORE))   GOTO 10
      RCORE1=R
      MCORE1=M
      ACORE1=A
      LCORE1=L
   10 IF ((ROOLD.LT.ROCORE).OR.(RO.GT.ROCORE))   GOTO 11
      RCORE2=R
      MCORE2=M
      ACORE2=A
      LCORE2=L
   11 RETURN
      END
C--------------------------------------------------------------------
      SUBROUTINE DIF9 (P,Y,F)
      IMPLICIT DOUBLE PRECISION (A-H,I,L,M,N,O-Z)
      DIMENSION Y(10),F(10)
      COMMON /AAA/ JFLAG
      COMMON /CCC/ PCORE,ROCORE
      COMMON /DDD/ C0,MNEUTR
      DATA PI4,GC2 /12.5663706143592D0,7.4250581D-29/
      R=Y(1)
      M=Y(2)
      PHI=Y(3)
      A=Y(4)
      OMEGA=Y(5)
      L=Y(6)
      MPROP=Y(7)
      CALL STATE (N,RO,P)
      DELTA=1.D0-2.D0*GC2*M/R
      IF (DELTA.LE.0.)   JFLAG=1
      IF (DELTA.LE.0.)   DELTA=-DELTA
      SQ=DSQRT(DELTA)
      DMDR=PI4*R**2*RO
      DPHIDR=GC2*(M/R**2+PI4*R*P)/DELTA
      DPDR=-(RO+P)*DPHIDR
      DNDR=PI4*R**2*N/SQ
      DOMEDR=6.D0/R**4*L*DEXP(PHI)/SQ
      DLDR=2.D0/3.D0*GC2*PI4*OMEGA*R**4*(RO+P)*
     1 DEXP(-PHI)/SQ
      DMPDR=PI4*R**2*RO/SQ
      F(1)=1.D0/DPDR
      F(2)=DMDR/DPDR
      F(3)=DPHIDR/DPDR
      F(4)=DNDR/DPDR
      F(5)=DOMEDR/DPDR
      F(6)=DLDR/DPDR
      F(7)=DMPDR/DPDR
C***
      F(8)=1.D0
      F(9)=1.D0
      IF (P.LT.PCORE)   F(8)=0.D0
      IF (RO.LT.ROCORE)   F(9)=0.D0
C***
      RETURN
      END
C--------------------------------------------------------------------
      SUBROUTINE STATE (N,RO,P)
      IMPLICIT DOUBLE PRECISION (A-H,I,L,M,N,O-Z)
      DIMENSION PDAT(9999),RHODAT(9999),NDAT(9999)
      COMMON /AAA/ JFLAG
      COMMON /BBB/ NDAT,RHODAT,PDAT,KDAT
      COMMON /CCC/ PCORE,ROCORE
      COMMON /DDD/ C0,MNEUTR
      KDAT1=KDAT-1
      IF (P.GT.PDAT(1))   GOTO 1
      JFLAG=1
      GOTO 10
    1 DO 3 K=1,KDAT1
      IF ((P.GT.PDAT(K)).AND.(P.LE.PDAT(K+1)))   GOTO 4
    3 CONTINUE
      WRITE (*,6005)
      STOP 111
    4 PHELP=(DLOG(P)-DLOG(PDAT(K)))/(DLOG(PDAT(K+1))-DLOG(PDAT(K)))
      N=NDAT(K)*DEXP((DLOG(NDAT(K+1))-DLOG(NDAT(K)))*PHELP)
      B=(DLOG(PDAT(K+1))-DLOG(PDAT(K)))/(DLOG(NDAT(K+1))-DLOG(NDAT(K)))
      RO=RHODAT(K)*DEXP((DLOG(RHODAT(K+1))-DLOG(RHODAT(K)))*PHELP)
c      RHO=RHODAT(K)+(P/N-PDAT(K)/NDAT(K))/(B-1.D0)*C0**2*1.D-19
c      RO=(EPS/C0**2*1.D+19+MNEUTR*1.D+15)*N
   10 RETURN
 6005 FORMAT(1X//1X,'(INPUT) P IS OUTSIDE THE RANGE OF PDAT(K)'//)
      END
C--------------------------------------------------------------------
      SUBROUTINE RK78 (IFLAG,N,T,DT,X,TOL,XDUM,F1,F2,F3,F4,F5,F6,F7,DER)
C
C     CHANGES FOR NNS -- PARAMETER IFLAG1 IS INTRODUCED.
C
C     VERSION: 1-SEP-81
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),XDUM(N),F1(N),F2(N),F3(N),F4(N),F5(N),F6(N),F7(N)
      DIMENSION CH(13),ALPH(13)
      COMMON /AAA/ IFLAG1
      DATA IFIRST/1/
C
c      IF (IFIRST.NE.1)  GOTO 70
C
      CH(1)=0.D0
      CH(2)=0.D0
      CH(3)=0.D0
      CH(4)=0.D0
      CH(5)=0.D0
      CH(11)=0.D0
      CH(6)=34.D0/105.D0
      CH(7)=9.D0/35.D0
      CH(8)=CH(7)
      CH(9)=9.D0/280.D0
      CH(10)=CH(9)
      CH(12)=41.D0/840.D0
      CH(13)=CH(12)
      ALPH(1)=0.D0
      ALPH(2)=2.D0/27.D0
      ALPH(3)=1.D0/9.D0
      ALPH(4)=1.D0/6.D0
      ALPH(5)=5.D0/12.D0
      ALPH(6)=0.5D0
      ALPH(7)=5.D0/6.D0
      ALPH(8)=1.D0/6.D0
      ALPH(9)=2.D0/3.D0
      ALPH(10)=1.D0/3.D0
      ALPH(11)=1.D0
      ALPH(12)=0.D0
      ALPH(13)=1.D0
      B2 1=2.D0/27.D0
      B3 1=1.D0/36.D0
      B4 1=1.D0/24.D0
      B5 1=5.D0/12.D0
      B6 1=0.05D0
      B7 1=-25.D0/108.D0
      B8 1=31.D0/300.D0
      B9 1=2.D0
      B10 1=-91.D0/108.D0
      B11 1=2383.D0/4100.D0
      B12 1=3.D0/205.D0
      B13 1=-1777.D0/4100.D0
      B3 2=1.D0/12.D0
      B4 3=1.D0/8.D0
      B5 3=-25.D0/16.D0
      B5 4=-B5 3
      B6 4=0.25D0
      B7 4=125.D0/108.D0
      B9 4=-53.D0/6.D0
      B10 4=23.D0/108.D0
      B11 4=-341.D0/164.D0
      B13 4=B11 4
      B6 5=0.2D0
      B7 5=-65.D0/27.D0
      B8 5=61.D0/225.D0
      B9 5=704.D0/45.D0
      B10 5=-976.D0/135.D0
      B11 5=4496.D0/1025.D0
      B13 5=B11 5
      B7 6=125.D0/54.D0
      B8 6=-2.D0/9.D0
      B9 6=-107.D0/9.D0
      B10 6=311.D0/54.D0
      B11 6=-301.D0/82.D0
      B12 6=-6.D0/41.D0
      B13 6=-289.D0/82.D0
      B8 7=13.D0/900.D0
      B9 7=67.D0/90.D0
      B10 7=-19.D0/60.D0
      B11 7=2133.D0/4100.D0
      B12 7=-3.D0/205.D0
      B13 7=2193.D0/4100.D0
      B9 8=3.D0
      B10 8=17.D0/6.D0
      B11 8=45.D0/82.D0
      B12 8=-3.D0/41.D0
      B13 8=51.D0/82.D0
      B10 9=-1.D0/12.D0
      B11 9=45.D0/164.D0
      B12 9=3.D0/41.D0
      B13 9=33.D0/164.D0
      B11 10=18.D0/41.D0
      B12 10=6.D0/41.D0
      B13 10=12.D0/41.D0
      B13 12=1.D0
      IFIRST=0
C
   70 IFLAG=0
C
   71 IF (IFLAG.EQ.1)   IFLAG1=0
      CALL DER (T,X,F1)
      T2=T+ALPH(2)*DT
      T3=T+ALPH(3)*DT
      T4=T+ALPH(4)*DT
      T5=T+ALPH(5)*DT
      T6=T+ALPH(6)*DT
      T7=T+ALPH(7)*DT
      T8=T+ALPH(8)*DT
      T9=T+ALPH(9)*DT
      T10=T+ALPH(10)*DT
      T11=T+DT
      T12=T
      T13=T+DT
      DO 72 I=1,N
   72 XDUM(I)=X(I) + DT*B2 1*F1(I)
      CALL DER (T2,XDUM,F2)
      DO 73 I=1,N
   73 XDUM(I)=X(I) + DT*(B3 1*F1(I)+B3 2*F2(I))
      CALL DER (T3,XDUM,F3)
      DO 74 I=1,N
   74 XDUM(I)=X(I)+DT*(B4 1*F1(I)+B4 3*F3(I))
      CALL DER (T4,XDUM,F4)
      DO 75 I=1,N
   75 XDUM(I)=X(I)+DT*(B5 1*F1(I)+B5 3*F3(I)+B5 4*F4(I))
      CALL DER (T5,XDUM,F5)
      DO 76 I=1,N
   76 XDUM(I)=X(I)+DT*(B6 1*F1(I)+B6 4*F4(I)+B6 5*F5(I))
      CALL DER (T6,XDUM,F6)
      DO 77 I=1,N
   77 XDUM(I)=X(I)+DT*(B7 1*F1(I)+B7 4*F4(I)+B7 5*F5(I)+B7 6*F6(I))
      CALL DER (T7,XDUM,F7)
      DO 78 I=1,N
   78 XDUM(I)=X(I)+DT*(B8 1*F1(I)+B8 5*F5(I)+B8 6*F6(I)+B8 7*F7(I))
      CALL DER (T8,XDUM,F2)
      DO 79 I=1,N
   79 XDUM(I)=X(I)+DT*(B9 1*F1(I)+B9 4*F4(I)+B9 5*F5(I)+B9 6*F6(I)
     C +B9 7*F7(I)+B9 8*F2(I))
      CALL DER (T9,XDUM,F3)
      DO 80 I=1,N
      X4=F4(I)
      X5=F5(I)
      X6=F6(I)
      X7=F7(I)
      X8=F2(I)
      X9=F3(I)
      F2(I)=CH(6)*X6+CH(7)*X7+CH(8)*X8+CH(9)*X9
      XDUM(I)=X(I)+DT*(B10 1*F1(I)+B10 4*X4+B10 5*X5+B10 6*X6+B10 7*X7 +
     C B10 8*X8+B10 9*X9)
      F4(I)=B11 1*F1(I)+B11 4*X4+B11 5*X5+B11 6*X6+B11 7*X7+B11 8*X8
     C +B11 9*X9
      F5(I)=B12 1*F1(I)+B12 6*X6+B12 7*X7+B12 8*X8+B12 9*X9
   80 F6(I)=B13 1*F1(I)+B13 4*X4+B13 5*X5+B13 6*X6+B13 7*X7+
     C B13 8*X8+B13 9*X9
      CALL DER (T10,XDUM,F3)
      DO 81 I=1,N
      XDUM(I)=X(I)+DT*(F4(I)+B11 10*F3(I))
   81 F1(I)=XDUM(I)
      CALL DER (T11,XDUM,F4)
      DO 82 I=1,N
   82 XDUM(I)=X(I)+DT*(F5(I)+B12 10*F3(I))
      CALL DER (T12,XDUM,F5)
      DO 83 I=1,N
      XDUM(I)=X(I)+DT*(F6(I)+B13 10*F3(I)+B13 12*F5(I))
   83 F7(I)=XDUM(I)
      CALL DER (T13,XDUM,F6)
      DO 84 I=1,N
   84 XDUM(I)=X(I)+DT*(CH(10)*F3(I)+CH(11)*F4(I)+CH(12)*F5(I)
     C + CH(13)*F6(I)+F2(I))
C
C     ...THE ESTIMATION OF ERRORS...
C
      DOUBLE=2.D0
      DO 85 I=1,N
      ERR1=DABS(F7(I)-F1(I))
      ERR2=DABS(F7(I)-XDUM(I))
      ERR3=DABS(F1(I)-XDUM(I))
      ERR=DMAX1(ERR1,ERR2,ERR3)
      TOLABS=DABS(TOL*XDUM(I))
      IF (XDUM(I).EQ.0.D0)   GOTO 85
      IF (ERR.GT.0.03125D0*TOLABS)   DOUBLE=1.D0
      IF (ERR.LT.TOLABS)   GOTO 85
      DT=DT/2.D0
      IFLAG=1
      GOTO 71
   85 CONTINUE
      IF (IFLAG.EQ.1)   DOUBLE=1.D0
C
C     ...STEP WAS EVALUATED...
C
      T=T+DT
      DT=DT*DOUBLE
      DO 86 I=1,N
   86 X(I)=XDUM(I)
      RETURN
      END
C--------------------------------------------------------------------
c..
c    Program otwiera zbior o nr. 2 i nazwie f2 i ustawia sie
c    na jego koncu
c..
        subroutine opn2 (f2)
        character*1 b
        character*30 f2
        close(unit=2)
        open(unit=2,status='old',file=f2
     $ ,form='formatted',err=5)
1       read(2,2,end=3)b
2       format(a1)
        goto 1
3       backspace 2
6       return
5       close(unit=2)
        open(unit=2,status='new',file=f2
     $ ,form='formatted')
        goto 6
        end
c..
c    Program otwiera zbior o nr. 4 i nazwie f4 i ustawia sie
c    na jego koncu
c..
        subroutine opn4 (f4)
        character*1 b
        character*30 f4
        close(unit=4)
        open(unit=4,status='old',file=f4
     $ ,form='formatted',err=5)
1       read(4,2,end=3)b
2       format(a1)
        goto 1
3       backspace 4
6       return
5       close(unit=4)
        open(unit=4,status='new',file=f4
     $ ,form='formatted')
        goto 6
        end
