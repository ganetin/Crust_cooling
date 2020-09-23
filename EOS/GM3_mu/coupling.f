      SUBROUTINE FIND_CONS ( RHO0,EA,ESYM,MS0,COMPRES,XSIGMA )
C------------------------------------------
C              Hyperon version
C------------------------------------------
      IMPLICIT NONE
      include 'const.inc'
      integer inpf
      parameter (inpf=1)
      REAL*8 RHO0,EA,ES0,KF0,MS0,ESYM,XSIGMA
      REAL*8 PHI,FPHI,DFPHI,D2UDPHI,fac1,fac2,fac3,COMPRES
      REAL*8 FUN,RHOS0
      REAL*8 MS2, MW2, MR2, GS, GW, GR, PI2, HC
      PARAMETER (PI2=9.8696044,  HC=197.33)
      REAL*8 C_S,C_W,C_R,BFAC,CFAC,GSK,GWK,GRK, W0, C3
      COMMON /COUPLING/ GS,GW,GR,BFAC,CFAC,GSK,GWK,GRK, C3

C------------------------------------------
C              Constants
C------------------------------------------
      MS2 = MS*MS
      MW2 = MW*MW
      MR2 = MR*MR
      MS0  = MS0*MN

      if (inpf.eq.1) then
        write(*,*) ' ------------- GM Lagrangian  ------------------- '
        PHI   = MN*(1.D0-MS0/MN)              ! GM
        FPHI  = 1.D0                          ! GM
        DFPHI = 0.D0                          ! GM
      elseif(inpf.eq.2) then
        write(*,*) ' ------------- ZM Lagrangian  ------------------- '
        PHI   = MN*(MN/MS0-1.D0)                 ! ZM
        FPHI  = 1.D0/(1.D0+PHI/MN)**2            ! ZM
        DFPHI = -2.D0/(1.D0+PHI/MN)**3/MN        ! ZM
      endif

C------------------------------------------
C------------------------------------------

      KF0 = (1.5D0*PI2*RHO0)**(1.D0/3.D0)
      ES0   = DSQRT(KF0*KF0+MS0*MS0) 
      RHOS0 = (MS0/PI2)*(KF0*ES0-MS0*MS0*DLOG((ES0+KF0)/MS0))

      D2UDPHI = (6.D0/PHI**2)*
     &        ( RHO0*(0.5*ES0+EA-MN)+RHOS0*(PHI*FPHI+0.5D0*MS0) )

      W0 = MN-EA-ES0
cc      C3 = 0.00278
      C3 = 0.d0

      C_W = W0/(RHO0-C3*W0**3)
      C_R = (ESYM-KF0*KF0/(6.D0*ES0))*(8.D0/RHO0) 

      fac3 = 3.d0*(RHO0/ES0-RHOS0/MS0)
      fac2 = COMPRES/(9.D0*HC) + (EA+ES0-MN) - KF0*KF0/(3.D0*ES0)
      fac1 = RHO0*(MS0/ES0)**2/fac2 - fac3

      fun = -FPHI*FPHI*fac1 + RHOS0*DFPHI - D2UDPHI

      C_S = 1.D0/FUN

C------------------------------------------
C------------------------------------------

      GS = DSQRT(MS2*C_S)
      GW = DSQRT(MW2*C_W)
      GR = DSQRT(MR2*C_R)

      CFAC = ( D2UDPHI + 2.D0/C_S - 2.D0*RHOS0*FPHI/PHI ) / PHI**2
      BFAC = ( D2UDPHI - 3.D0*CFAC*PHI*PHI ) / (2.D0*MN*PHI)

      GSK = XSIGMA*GS
      GWK = (-28.0/HC + XSIGMA*PHI)/(C_W*RHO0)*GW
      GRK = XSIGMA*GR

      print*, -bfac*mn*gs**3, cfac*gs**4, C3*gw**4

      RETURN
      END

