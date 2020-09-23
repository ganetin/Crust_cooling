C cgs units for P and rho
C nb in fm^-3
CCC Ep=938.272046,En=939.565379
      IMPLICIT NONE
 3334 FORMAT (1X,17G16.8)
c Input
      integer KEOS,K,i,KDAT
      real*8 nb(1000),Pres(1000),rho(1000),mass(1000),rad(1000)
      real*8 phi(1000),PLG
      double precision dattab(8),JUNK,XN

      KEOS=21

      open (3,file='Prof_BSK21_Cat_Fe_1.4_Ioffe_extended_ori.dat')
c 1          ,status='old')
         K=0
 7       K=K+1
 71      READ (3,*,END=8,ERR=71) (dattab(i),i=1,8)
         mass(K)=dattab(1)
         rad(K)=dattab(2)*1e3
         pres(K)=dattab(3)
         rho(K)=dattab(4)
         phi(K)=dattab(5)
         GOTO 7
 8       CLOSE (1)
         KDAT=K-1
      close(3)
      
      do i=1,KDAT
       call BSkNofR(KEOS,rho(i),XN)
       nb(i)=XN
      enddo

      write(91,*)
      write(91,*) '    EOS file :  '
      write(91,*)
      write(91,10002)'step','radius ','baryon#  ','density   ',
     x               '  pressure  ','encl. mass   ','phi   ',
     x               'encl. bar. mass '
      write(91,10002) '   ','  (m)  ','(#/fm3)  ','(g/cm3)   ',
     x               '  (dyn/cm2) ',' (sol. mass)  ','     ',
     x               ' (sol. mass)   '
      write(91,*)
     
      do i=1,KDAT
      write(91,10000)i-1,rad(i),nb(i),rho(i),pres(i),mass(i),phi(i),42.
      enddo
10000 format(i6,0pf15.6,1p1e15.6,1e15.6,1e15.5,1e18.9,1e15.6,1e18.9)
10002 format(a6,a15,a15,a15,a15,a18,a15,a18)

c Output
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


 
!c NSCool Input Core
!c   Rho        Press       nbar       Ye         Ymu        Yn         Yp         Yla        Ysm        Ys0        Ysp       mstp       mstn       mstla      mstsm      msts0      mstsp
!c NSCool Input Crust
!c   Rho        Press       nbar       [A_cell]    [A_ion]     [Z]
!      ENDDO

      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Fitting functions reproducing EOSs BSk19,20,21.
* Last change: 17.02.2013.
* Contact: Alexander Potekhin <palex@astro.ioffe.ru>
* List of subroutines:
*   BSKfit returns pressure P and log.derivative d ln P / d ln\rho
*     at a given mass density \rho (at \rho > 10^6 g/cc);
*   BSKfitH returns \rho (at \rho > 10^6 g/cc) as a function of
*     H1 = e^H-1, where H is the dimensinoless pseudo-enthalpy;
*   BSkRofN returns mass density \rho as a function of n_b;
*   BSkNofR returns baryon number density n_b as a function of \rho;
*   FRACORE returns fractions of electrons Y_e and muons Y_\mu
*     as functions of n_b in the stellar core;
*   INCRUST returns the number of bound nucleons in a nucleus A,
*     the total number of nucleons per one nucleus A',
*     nuclear-shape parameters a_p, a_n, C_p, C_n, 
*     as well as additional parameters n_{0p} and n_{0n} - the heights
*     of the bumps of proton and neutron densities near the center of
*     a Wigner-Seitz cell, and the ratios xnuc and xnucn
*     of the effective proton and neutron radii of these bumps,
*     respectively, to the cell radius.
 
      subroutine BSKfit(KEOS,RLG,PLG,CHIrho,Gamma)
* Fit of P(\rho) for BSk19,20,21 at rho > 10^6 g/cc
* average error 1.%; max.error 3.7% at RLG=9.51
* Input: KEOS - number of BSK model (19,20,or 21)
*        RLG = lg(rho[g/cc])
* Output: PLG = lg(P[dyn/cm^2]),
*         CHIrho = d log P / d log rho
*         Gamma = d log P / d log n_b
*                                                          Version 12/11
      implicit double precision (A-H), double precision (O-Z)
      parameter(c2=2.99792458d10**2)
      dimension A(23,3)
      save A
      data A/
* BSk19:
     *  3.916,7.701,.00858,.22114,3.269,11.964,13.349,1.3683,
     *    3.254,-12.953,0.9237,6.20,14.383,16.693,-1.0514,2.486,
     *    15.362,.085,11.68,6.23,-.029,14.19,20.1,
* BSk20:
     *  4.078,7.587,.00839,.21695,3.614,11.942,13.751,1.3373,
     *    3.606,-22.996,1.6229,4.88,14.274,23.560,-1.5564,2.095,
     *    15.294,.084,11.67,6.36,-.042,14.18,14.8,
* BSk21:
     *  4.857,6.981,.00706,.19351,4.085,12.065,10.521,1.5905,
     *    4.104,-28.726,2.0845,4.89,14.302,22.881,-1.7690,0.989,
     *    15.313,.091,11.65,4.68,-.086,14.15,10.0/
      if (RLG.lt.-1.) stop'BSKfit: too low RLG'
      if (RLG.gt.16.5) stop'BSKfit: too high RLG'
      K=KEOS-18
      if (K.lt.1.or.K.gt.3) stop'BSKfit: unknown EOS'
      X=RLG
      D4=1.d0+A(4,K)*X
      P01=(A(1,K)+A(2,K)*X+A(3,K)*X**3)/D4
      F1=FERMI(A(5,K)*(X-A(6,K)))
      P02=A(7,K)+A(8,K)*X
      F2=FERMI(A(9,K)*(A(6,K)-X))
      P03=A(10,K)+A(11,K)*X
      F3=FERMI(A(12,K)*(A(13,K)-X))
      P04=A(14,K)+A(15,K)*X
      F4=FERMI(A(16,K)*(A(17,K)-X))
      G05=(X-A(19,K))*A(20,K)
      P05=A(18,K)/(1.d0+G05**2)
      G06=(X-A(22,K))*A(23,K)
      P06=A(21,K)/(1.d0+G06**2)
      PLG=P01*F1+P02*F2+P03*F3+P04*F4+P05+P06
      CHIrho=F1/D4* ! d\zeta/dX=d ln P/d ln rho
     *  ((A(2,K)-A(1,K)*A(4,K)+3.*A(3,K)*X**2+
     +    2.*A(3,K)*A(4,K)*X**3)/D4-
     -  A(5,K)*(A(1,K)+A(2,K)*X+A(3,K)*X**3)*(1.d0-F1)) +
     +  F2*(A(8,K)+A(9,K)*(1.d0-F2)*P02) +
     +  F3*(A(11,K)+A(12,K)*(1.d0-F3)*P03)+
     +  F4*(A(15,K)+A(16,K)*(1.d0-F4)*P04)-
     -  2.d0*(P05*G05*A(20,K)/(1.d0+G05**2)+
     +        P06*G06*A(23,K)/(1.d0+G06**2))
      Gamma=(1.d0+10.d0**(PLG-RLG)/c2)*CHIrho
      return
      end

      subroutine BSKfitH(KEOS,H1,RLG)
* Fit of rho(H) for BSk19,20,21 at rho > 10^6 g/cc (i.e.,e^H-1 > 5.9e-5)
* average error 1.3%; max.error < 5% (at both inner-crust boundaries)
* Input: KEOS - number of BSK model (19,20,or 21)
*        H1 = e^H-1, where H = pseudo-enthalpy [Eq.(7) of Haens.Pot.04]
* Output: RLG = lg(rho[g/cc])
*                                                          Version 02/12
      implicit double precision (A-H), double precision (O-Z)
      save A
      dimension A(10,3)
      data A/
     *  93.650,36.893,15.450,.672,61.240,68.97,.292,5.2,.48,6.8, ! BSk19
     *  81.631,31.583,15.310,.594,58.890,56.74,.449,4.5,.58,7.5, ! BSk20
     *  63.150,23.484,15.226,.571,54.174,37.15,.596,3.6,.51,10.4/! BSk21
      X=dlog10(H1)
      K=KEOS-18
      if (K.lt.1.or.K.gt.3) stop'BSKfitH: unknown EOS'
      RLG1=2.367+21.84*H1**.1843/(1.+.7*H1)
      RLG3=(A(3,K)+A(4,K)*X)
      G=(A(5,K)*H1)**A(10,K)
      RLG2=(A(1,K)+A(2,K)*X+G*RLG3)/(1.+A(6,K)*H1+G)
      Y=84.*(X+1.99)
      FY=FERMI(Y)
      RLG=RLG1*FY+RLG2*dim(1.d0,FY)+FERMI(A(8,K)*(A(9,K)-X))*A(7,K)
      return
      end

      subroutine BSkRofN(KEOS,XN,RHO)
*                                                          Version 02/12
* Fit of (rho(n_b)-m_0*n_b) consistent with BSKfit
*   for BSk19,20,21 at rho > 10^6 g/cc (i.e., n_b > 6.e-10 fm^{-3})
* Criterion is F = \rho/(m_0 n_b)-1 = e^H - 1 - P/(m_0 n_b c^2);
* average error in F is 1.7%; max.error 4.2% at rho=10^6 g/cc = min(rho)
* Input: KEOS - number of BSK model (19,20,or 21)
*        XN = n_b - baryon number density in fm^{-3}
* Output: RHO - mass-energy density in [g/cc]
      implicit double precision (A-H), double precision (O-Z)
      save A
      dimension A(8,3)
      data A/
     *  .259,2.30,.0339,.1527,3.170,1.085E-05,.6165,3.50, ! BSk19
     *  .632,2.71,.0352,.383,3.165,1.087E-05,.6167,3.51,  ! BSk20
     *  3.85,3.19,.0436,1.99,3.210,1.075E-05,.6154,3.44/  ! BSk21
      K=KEOS-18
      if (K.lt.1.or.K.gt.3) stop'BSkRofN: unknown EOS'
      F1=(A(1,K)*XN**A(2,K)+A(3,K)*dsqrt(XN))/(1.d0+A(4,K)*XN)**2
      F2=XN/(A(6,K)+XN**A(7,K)*A(8,K))
      XLG=dlog10(XN)
      FRM=FERMI(1.1*(XLG+A(5,K)))
      F=FRM*F2+(1.d0-FRM)*F1
      RHO=XN*1.66d15*(1.d0+F)
      return
      end

      subroutine BSkNofR(KEOS,RHO,XN)
*                                                          Version 02/12
* Fit of (rho-m_0*n_b(rho)) consistent with BSKfit at rho > 10^6 g/cc.
* Criterion is F = \rho/(m_0 n_b)-1 = e^H - 1 - P/(m_0 n_b c^2);
* average error in F is 1.8% for BSk19, 1.5% for BSk20 and 21,
*   max.error in F is 3.9% for BSk19,20 and 3.7% for BSk21 (at min\rho)
* Input: KEOS - number of BSK model (19,20,or 21)
*        RHO - mass-energy density in [g/cc]
* Output: XN = n_b - baryon number density in fm^{-3}
      implicit double precision (A-H), double precision (O-Z)
      save A
      dimension A(9,3)
      data A/
* BSk19:
     *  .378,1.28,17.20,3.844,3.778,12.0741,1.071E-05,3.287,.613,
* BSk20:
     *  .152,1.02,10.26,3.691,2.586,12.0570,1.067E-05,3.255,.6123,
* BSk21:
     *  .085,.802,16.35,3.634,2.931,12.0190,1.050E-05,3.187,.611/
      K=KEOS-18
      if (K.lt.1.or.K.gt.3) stop'BSkNofR: unknown EOS'
      X=RHO/1.66d15
      F1=(A(1,K)*X**A(2,K)+A(3,K)*X**A(4,K))/(1.d0+A(5,K)*X)**3
      F2=X/(A(7,K)+A(8,K)*X**A(9,K))
      RLG=dlog10(RHO)
      FRM=FERMI(RLG-A(6,K))
      F=FRM*F2+(1.-FRM)*F1
      XN=X/(1.d0+F)     
      return
      end

      function FERMI(X)
      implicit double precision (A-H), double precision (O-Z)
      if (X.gt.40.d0) then
         F=0.d0
         goto 50
      endif
      if (X.lt.-40.d0) then
         F=1.d0
         goto 50
      endif
      F=1.d0/(dexp(X)+1.d0)
   50 FERMI=F
      return
      end

      subroutine FRACORE(KEOS,XN,Ye,Ymu)
*                                                          Version 05/12
* Fit of particle fractions in the core.
* Input: KEOS - number of BSK model (19,20,or 21)
*        XN = n_b - baryon number density in fm^{-3}
* Output: Ye - number of electrons divided by number of nucleons
*         Ymu - number of muons divided by number of nucleons
* NB: proton fraction is Yp=Ye+Ymu, neutron fraction is Yn=1-Yp.
* average error 0.00011 (Ye, BSk19), 0.00013 (Ymu, BSk19),
*               0.00010 (Ye, BSk20), 0.00012 (Ymu, BSk20),
*               0.00023 (Ye, BSk21), 0.00026 (Ymu, BSk21),
*     max.error 0.00024 (Ye, BSk19), 0.00034 (Ymu, BSk19),
*               0.00020 (Ye, BSk20), 0.00070 (Ymu, BSk20),
*               0.00042 (Ye, BSk21), 0.00058 (Ymu, BSk21),
      implicit double precision (A-H), double precision (O-Z)
      save A
      dimension A(12,3)
      data A/
* BSk19:
     *  -.0157,.9063,0.,26.97,106.5,4.82,
     *  -.0315,.25,0.,12.42,72.4,19.5,
* BSk20:
     *  -.0078,.745,.508,22.888,.449,.00323,
     *  -.0364,.2748,.2603,12.99,.0767,.00413,
* BSk21:
     *  .00575,.4983,9.673,16.31,38.364,0.,
     *  -.0365,.247,11.49,24.55,48.544,0./
      K=KEOS-18
      if (K.lt.1.or.K.gt.3) stop'FRACORE: unknown EOS'
      Ye=(A(1,K)+A(2,K)*XN+A(3,K)*XN**4)/
     /  (1.+A(4,K)*XN*dsqrt(XN)+A(5,K)*XN**4)*
     *  exp(-A(6,K)*XN**5)
      if (Ye.lt.0.) Ye=0.
      Ymu=(A(7,K)+A(8,K)*XN+A(9,K)*XN**4)/
     /  (1.+A(10,K)*XN*dsqrt(XN)+A(11,K)*XN**4)*
     *  exp(-A(12,K)*XN**5)
      if (Ymu.lt.0.) Ymu=0.
      return
      end

      subroutine INCRUST(KEOS,XN,XND,XNC,CMI,CMI1,Zcell,Zclust,
     *  ap,an,Cp,Cn,Dn_p,Dn_n,xnuc,xnucn)
*                                                       Version 17.02.13
* Fit of the nuclear-shape parameters.
* Input: KEOS - number of BSK model (19,20,or 21)
*        XN = n_b - baryon number density in fm^{-3}
* Output: XND - neutron-drip number density of baryons in fm^{-3} 
*              = the inner/outer crust interface (Pearson et al. 2012),
*         XNC - number density of baryons in fm^{-3} at the crust/core
*              interface (Pearson et al. 2012),
*         CMI - average number of bound nucleons in a nucleus,
*         CMI1 - average total number of nucleons per 1 nucleus,
*         Zcell - total number of protons per a nucleus,
*         Zclust - number of clusterized protons in a nucleus,
*         a_q,C_q,Dn_q (q=p,n) - shape parameters (Onsi et al. 2008),
*         xnuc, xnucn - proton and neutron equivalent size parameters,
*          defined as (Gnedin et al. 2003) xnuc = sqrt{5<r^2>/3}
* maximal fractional errors [measured for n_b up to max(n_b[fm^{-3}])]:
* BSk:          19                    20                   21   
* ap:     1.1%  [.065]            1.2%  [.066]           0.34% [.062]
* an:  5.9%[.082]/1.0%[.047] 6.0%[.081]/0.85%[.046]       2.0% [.065]
* Cp: 0.43%[.048]/3%[.0815]    1%[.055]/3.3%[.08]    1%[.052]/9%[.078]
* Cn: 1%[.053]/3%[.0815] 1%[.056]/3%[.0795] 0.7%[.061]/2%[.075]/7%[.082]
* Dn_p: 1%[.056]/7%[.0745]     1%[.059]/9%[.0815]    2%[.056]/9%[.078]
* Dn_n: 0.7% [.065]/7%[.081] 0.6%[.066]/7%[.079]   0.4%[.065]/7%[.079]
* xnuc:    5% [.0815]               4% [.080]      2.5%[.065]/6.3%[.078]
* xnucn:  4.6% [.086]               4% [.084]              2.5% [.082]
* CMI:  0.8%[.063]/7%[.081]  0.7%[.063]/6%[.079]   0.8%[.059]/8%[.079]
* CMI1: 0.6%[.054]/4%[.070]        1.5% [0.1]      0.6%[.065]/2%[.079]
* The listed max.errors are reached near the lower and/or higher ends
* of the XN range. The second (larger) values (after slash) are probably
* due to inaccuracies in the numerical data.
* Typical errors are 10 times smaller than the first (lower) max.errors.
      implicit double precision (A-H), double precision (O-Z)
* Fit parameters for shape parameters a_q,C_q,xnuc_q,dn_q (q=p,n):
      save PARap,PARan,PARCp,PARCn,PARxp,PARxn,PARdp,PARdn,PARNi,A1,LZ
      dimension PARap(3,3),PARan(4,3),PARCp(4,3),PARCn(4,3),
     *  PARxp(4,3),PARxn(4,3),PARdp(4,3),PARdn(3,3),
     *  PARNi(6,3),A1(6,3),LZ(2,3)
      data PARap/
     *  .4377,4.360,1084., ! BSk19
     *  .4353,4.440,1154., ! BSk20
     *  .4316,4.704,1253./ ! BSk21
      data PARan/
     *  .639,1.461,.457,1137., ! BSk19
     *  .632,1.98,0.514,1122., ! BSk20
     *  .636,5.32,0.739,624./  ! BSk21
      data PARCp/
     *  5.500,11.7,.643,472., ! BSk19
     *  5.493,12.8,.636,484., ! BSk20
     *  5.457,14.2,.601,566./ ! BSk21
      data PARCn/
     *  5.714,14.05,.642,182., ! BSk19
     *  5.714,16.3,0.645,175., ! BSk20
     *  5.728,22.2,0.663,144./ ! BSk21
      data PARxp/
     *  .1120,2.06,.633,507., ! BSk19
     *  .1094,2.04,.613,509., ! BSk20
     *  .1045,2.09,.586,513./ ! BSk21
      data PARxn/
     *  .122,2.27,.618,193., ! BSk19
     *  .119,2.30,.603,182., ! BSk20
     *  .114,2.56,.595,107./ ! BSk21
      data PARdp/
     *  .05509,.1589,  0., .4917, ! BSk19
     *  .05382,.1400,.01715,.566, ! BSk20
     *  .05273,.1107,.02218,.6872/ ! BSk21
      data PARdn/
     *  .10336,.0772,1.129, ! BSk19
     *  .10283,.0825,1.189, ! BSk20
     *  .10085,.0942,1.279/ ! BSk21
      data PARNi/ ! for number of bound neutrons
     *  93.0,11.90,1.490,.334,5.05, 8.40, ! BSk19
     *  92.8,12.95,1.493,.354,7.57, 9.30, ! BSk20
     *  92.3,13.80,1.625,.3874,13.8,10.8/ ! BSk21
      data A1/
     *  134.7,183.7,308.7,.3814,-.00058,.00049, ! BSk19
     *  134.7,188.2,275.6,.4346,.00163,.00149,  ! BSk20
     *  132.6,187.6,229.2,.5202,.00637,.00151/  ! BSk21
      data LZ/19,16, 20,19, 27,17/
      K=KEOS-18
      if (KEOS.eq.19) then
         XND=2.63464d-4
         XNC=.0885
      elseif (KEOS.eq.20) then
         XND=2.62873d-4
         XNC=.0854
      elseif (KEOS.eq.21) then
         XND=2.57541d-4
         XNC=.0809
      else
         stop'INCRUST: unknown EOS'
      endif
      X=XN/XND
      Xlg=dmax1(dlog10(X),0.d0)
      XD=XN/XNC
      Zcell=40.d0
      Zclust=Zcell*dexp(-XD**LZ(1,K))*dim(1.d0,XD**LZ(2,K))
      P3Xi=PARNi(3,K)*Xlg
      CNI=(PARNi(1,K)+PARNi(2,K)*Xlg+P3Xi**3*dsqrt(P3Xi))/
     /  (1.d0+(PARNi(4,K)*Xlg)**PARNi(5,K))*
     *  (dim(1.d0,XD**PARNi(6,K)))
      CMI=Zclust+CNI
      CMI1=(A1(1,K)+A1(2,K)*Xlg+A1(3,K)*Xlg**2)/(1+(A1(4,K)*Xlg)**4)*
     *  (1+A1(5,K)*X)*dim(1.d0,(A1(6,K)*X)**2)
      ap=(PARap(1,K)+PARap(2,K)*XN)/(1.d0-PARap(3,K)*XN**3)
      an=(PARan(1,K)+PARan(2,K)*XN**PARan(3,K))/(1.d0-PARan(4,K)*XN**3)
      Cp=(PARCp(1,K)+PARCp(2,K)*XN**PARCp(3,K))/(1.d0-PARCp(4,K)*XN**3)
      Cn=(PARCn(1,K)+PARCn(2,K)*XN**PARCn(3,K))/(1.d0-PARCn(4,K)*XN**3)
      xnuc=(PARxp(1,K)+PARxp(2,K)*XN**PARxp(3,K))/
     /  (1.d0-PARxp(4,K)*XN**3)
      xnucn=(PARxn(1,K)+PARxn(2,K)*XN**PARxn(3,K))/
     /  (1.d0-PARxn(4,K)*XN**3)
      Dn_p=PARdp(1,K)-PARdp(2,K)*XN/(PARdp(3,K)+XN**PARdp(4,K))
      Dn_p=Dn_p*dim(1.d0,XD**9)
      Dn_n=PARdn(1,K)+PARdn(2,K)*XN**.5-PARdn(3,K)*XN
      Dn_n=Dn_n*dim(1.d0,XD**16)
      return
      end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Effective nucleon masses.
      subroutine EFFMASS(KEOS,XN,Yp,EFMp,EFMn)
*                                                       Version 11.04.13
*                                                     Corrected 19.03.14
* Proton and neutron effective masses in the nuclear matter.
* Calculation according to Eq.(A10) of Chamel, Goriely, & Pearson (2009)
*     with parameters from Table IV of Goriely, Chamel, & Pearson (2010)
* The code by A.Y.Potekhin <palex@astro.ioffe.ru>
* Input: KEOS - number of BSK model (19,20,or 21)
*        XN = n_b - baryon number density in fm^{-3}
*        Yp - proton fraction
* Output: EFMp, EFMn - ratios of the proton and neutron effective masses
*        to the masses of an isolated proton and neutron, respectively.
      implicit double precision (A-H), double precision (O-Z)
      save AT1,AT2X2,AT4,AT5,AX1,AX4,AX5
      dimension AT1(3),AT2X2(3),AT4(3),AT5(3),AX1(3),AX4(3),AX5(3)
      parameter (HbarC=197.3269718) ! \hbar*c [Mev*fm]
      parameter (Ep=938.272046,En=939.565379) ! m_p*c^2, m_n*c^2 [MeV]
      data AT1/403.072,438.219,396.131/,
     *  AT2X2/-1055.55,-1147.64,-1390.38/,
     *  AT4/-60.,-100.,-100./
     *  AT5/-90.,-120.,-150./,
     *  AX1/-0.137960,-0.392047,0.0648452/, ! corrected 19.03.14
     *  AX4/-6.,-3.,2./,
     *  AX5/-13.,-11.,-11./
      K=KEOS-18
      Yn=dim(1.d0,Yp)
      if (K.eq.1) then
         BETA=1.d0/3.d0
      elseif (K.eq.2) then
         BETA=1.d0/6.d0
      elseif (K.eq.3) then
         BETA=.5d0
      else
        stop'EFFMASS: unknown EOS'
      endif
      GAMMA=1.d0/12.d0
      DFMp=AT1(K)*((1.d0+.5d0*AX1(K))-(.5d0+AX1(K))*Yp)+
     +  AT2X2(K)*(.5d0+Yp)+
     +  AT4(K)*((1.d0+.5d0*AX4(K))-(.5d0+AX4(K))*Yp)*XN**BETA+
     +  AT5(K)*((1.d0+.5d0*AX5(K))+(.5d0+AX5(K))*Yp)*XN**GAMMA
      DFMn=AT1(K)*((1.d0+.5d0*AX1(K))-(.5d0+AX1(K))*Yn)+
     +  AT2X2(K)*(.5d0+Yn)+
     +  AT4(K)*((1.d0+.5d0*AX4(K))-(.5d0+AX4(K))*Yn)*XN**BETA+
     +  AT5(K)*((1.d0+.5d0*AX5(K))+(.5d0+AX5(K))*Yn)*XN**GAMMA
      EFMp=1.d0/(1.d0+.5d0*DFMp*XN*Ep/HbarC**2)
      EFMn=1.d0/(1.d0+.5d0*DFMn*XN*En/HbarC**2)
      return
      end
