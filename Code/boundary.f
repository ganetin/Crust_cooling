c
c      File cleaned up on Oct. 10, 2007
c

       function fteff(Tb,ifteff,etabound,bfield,istep,time,
     1        Ts1,Ts2,  Z,A,Rho,debug)
c WARNING : Tb must be the local value, not the redshifted one !
c           fteff is also non red-shifted !
c           Ts1 and Ts2 ARE red-shifted !
       implicit real*8 (a-h,k-z)
        if (debug.ge.1.) print *,'Entering fteff'
        if (ifteff.eq.0) then
         fteff=fteff_table(Tb,Ts1,Ts2)
        else if (ifteff.eq.1) then
         fteff=fteff_GPE(Tb)
        else if (ifteff.eq.2) then
         fteff=fteff_NT(Tb)
        else if (ifteff.eq.3) then
         fteff=fteff_acc(Tb,etabound)
        else if (ifteff.eq.4) then
         fteff=fteff_field_iron(Tb,bfield)
        else if (ifteff.eq.5) then
         fteff=fteff_ZARho(Tb,Z,A,Rho)
        else if (ifteff.eq.6) then
         fteff=Tb
        else if (ifteff.eq.7) then
         fteff=fteff_fullacc(Tb)
        else if (ifteff.eq.8) then
         fteff=fteff_iron(Tb)
        else if (ifteff.eq.9) then
         fteff=fteff_ANS(Tb,etabound)
        else if (ifteff.eq.10) then
         fteff=fteff_HZ90(Tb,etabound)
        else if (ifteff.eq.21) then
         icompo=1
         deltaMM=etabound
         fteff=fteff_Bez_DMM(Tb, deltaMM,icompo)
        else if (ifteff.eq.22) then
         icompo=2
         deltaMM=etabound
         fteff=fteff_Bez_DMM(Tb, deltaMM,icompo)
        else if (ifteff.eq.23) then
         icompo=3
         deltaMM=etabound
         fteff=fteff_Bez_DMM(Tb, deltaMM,icompo)
        else
         write(*,*)'Houston: ifteff out of range',ifteff
        end if

        if (debug.ge.1.) print *,'Done'

      return
       end
c *********************************************************************
c *********************************************************************

      function fteff_ZARho(Tb,Z,A,Rho)
c FORTIN 4/02/11 Correction of the numerical factor
c Check vs. fig 1 of HA84
c This is the Hernquist & Applegate: ApJ 287, 244 (1984), analytical
c boundary condition for an arbitrary density (their Equ (3.12).

c WARNING : Tb must be the local value, not the redshifted one !
     
      implicit real*8 (a-h,k-z)
      common/gravity/gs14,compactness
      save print

      if (print.eq.1.) goto 10
       print *,'-----------------------------------------------------'
       print *,'Using Hernquist & Applegate boundary condition at  '
       print *,'  Rho =',Rho
       print *,'-----------------------------------------------------'
       print=1.
 10   continue

c      fteff_ZARho = gs14**0.25 * 
c     1              (A**3/Z**4 / 4.25d30 *Tb**6.5/Rho**2)**0.25
      fteff_ZARho = gs14**0.25 * 
     1              (A**3/Z**4 / 4.93d16 *Tb**6.5/Rho**2)**0.25
      end

c**********************************************************************
c**********************************************************************

      function ftb_ZARho(Teff,Z,A,Rho)
c FORTIN 9/02/11 Correction of the numerical factor
c Check vs. fig 1 of HA84
c This is the Hernquist & Applegate: ApJ 287, 244 (1984), analytical
c boundary condition for an arbitrary density (their Equ (3.12).

c WARNING : Tb must be the local value, not the redshifted one !
     
      implicit real*8 (a-h,k-z)
      common/gravity/gs14,compactness
      save print

      if (print.eq.1.) goto 10
       print *,'-----------------------------------------------------'
       print *,'Tb(Teff) ' 
       print *,'Using Hernquist & Applegate boundary condition at  '
       print *,'  Rho =',Rho
       print *,'-----------------------------------------------------'
       print=1.
 10   continue

      ftb_ZARho = (Teff**4.*Rho**2*4.93d16*Z**4/
     1             A**3/gs14)**(2./13.)
      end

C************************************************************************
C************************************************************************


       function fteff_table(Tb,Ts1,Ts2)
c checked on nov. 6 1990 and July 2006
c WARNING : Tb must be the local value, not the redshifted one !
c
c Ts1 & Ts2 are two auxiliary temperatures from the envelope calculations,
c They could be, e.g., the Max and min surface T, or anything wanted !
c Ts1 and Ts2 are NOT used for the cooling, they are only informative.
c
       implicit real*8 (a-h,k-z)
       character*100 f_TbTs
       common/bound/Tb_acc0,f_TbTs
       parameter(jtpmax=100)
       dimension tempb(0:jtpmax),temps(0:jtpmax),t2(0:jtpmax),
     1           temps1(0:jtpmax),t21(0:jtpmax),
     2           temps2(0:jtpmax),t22(0:jtpmax)
       common/gravity/gs14,compactness
       save read,tempb,temps,t2
        if (read.eq.1.) goto 999
c ************ read temps and tempb ************************************
         print *,'-------------------------------------------------'
         print *,'Using envelope boundary condition from table'
         print *,'     ',f_TbTs
         print *,'WARNING:   No gs14 correction applied ! '
         print *,'-------------------------------------------------'
        open(unit=15,file=f_TbTs,status='old')
         read(15,*)jmax
         jmax=jmax+1
         if (jmax.gt.jtpmax) then
          print '(a20,i5)',' jtpmax =',jtpmax
          pause 'function_fteff_table: jmax > jtpmax '
         end if
         do j=1,3
          read(15,*)
         end do
         do j=1,jmax
c---------------------
c This assumes Ts1 & Ts2 are defined in the table:
c          read(15,*)temps(j),tempb(j),temps1(j),temps2(j)
c---------------------
c If not, then use this:
          read(15,*)temps(j),tempb(j)
          temps1(j)=tempb(j)
          temps2(j)=tempb(j)
c---------------------
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c HERE DANY: this version assumes the gs14^(1/4) effect is already 
c            included in the table !
c          temps(j)=temps(j)*gs14**0.25
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         end do
        close(unit=15,status='keep')
        do j=jmax+1,jtpmax
         tempb(j)=tempb(jmax)+
     1            float(j-jmax)*(tempb(jmax)-tempb(jmax-1))
         temps(j)=temps(jmax)+
     1            float(j-jmax)*(temps(jmax)-temps(jmax-1))
         temps1(j)=temps1(jmax)+
     1            float(j-jmax)*(temps1(jmax)-temps1(jmax-1))
         temps2(j)=temps2(jmax)+
     1            float(j-jmax)*(temps2(jmax)-temps2(jmax-1))
        end do
        dt0=(temps(2)-temps(1))/(tempb(2)-tempb(1))
        dt1=(temps(jtpmax)-temps(jtpmax-1))/
     1      (tempb(jtpmax)-tempb(jtpmax-1))
        call spline1(tempb,temps,jtpmax,dt0,dt1,t2)
        dt0=(temps1(2)-temps1(1))/(tempb(2)-tempb(1))
        dt1=(temps1(jtpmax)-temps1(jtpmax-1))/
     1      (tempb(jtpmax)-tempb(jtpmax-1))
        call spline1(tempb,temps1,jtpmax,dt0,dt1,t21)
        dt0=(temps2(2)-temps2(1))/(tempb(2)-tempb(1))
        dt1=(temps2(jtpmax)-temps2(jtpmax-1))/
     1      (tempb(jtpmax)-tempb(jtpmax-1))
        call spline1(tempb,temps2,jtpmax,dt0,dt1,t22)
        read=1.
c ***************************************************************************
999     continue

        call splint1(tempb,temps ,t2 ,jtpmax,Tb,Ts )
        call splint1(tempb,temps1,t21,jtpmax,Tb,Ts1)
        call splint1(tempb,temps2,t22,jtpmax,Tb,Ts2)
       
        fteff_table=Ts

       return
      end
c *********************************************************************
c *********************************************************************

      function fteff_NT(Tb)

c this is the Nomoto & Tsuruta boundary condition at rho=1.e10 g/cm3

c WARNING : Tb must be the local value, not the redshifted one !
     
      implicit real*8 (a-h,k-z)
      common/gravity/gs14,compactness
      dimension ltb(11),lts(11)

      data lts/5.698,5.764,5.943,6.103,6.228,
     1         6.363,6.493,6.633,6.773,6.918,7.068/
      data ltb/7.67,7.75,8.00,8.25,8.50,
     2         8.75,9.00,9.25,9.50,9.75,10.00/

      save lts,ltb
      save print
c*********************
      if (print.eq.1.) goto 10
       print *,'-------------------------------------------------'
       print *,'Using Nomoto & Tsuruta boundary condition'
       print *,'-------------------------------------------------'
       print=1.
 10   continue
c*********************

      logtb=dlog10(Tb)

      if (logtb.le.ltb(1)) then
       i1=1
      else if (logtb.ge.ltb(11)) then
       i1=10
      else
       do i=1,10
        if ((logtb.ge.ltb(i)).and.(logtb.le.ltb(i+1))) then
         i1=i
         goto 100
        end if
       end do
100    continue
      end if

      i2=i1+1
      deltb=ltb(i2)-ltb(i1)
      w1=(ltb(i2)-logtb)/deltb
      w2=1.-w1
      logts=lts(i1)*w1+lts(i2)*w2
      logts=logts+0.25*dlog10(gs14)

      fteff_NT=10.**logts      
 
      end

c *********************************************************************
c *********************************************************************

      function fteff_GPE(Tb)

c this is the Gudmundsson, Pethick & Epstein, ApJ 259, L19 (1982)
c boundary condition at rho=1.e10 g/cm3

c WARNING : Tb must be the local value, not the redshifted one !
     
      implicit real*8 (a-h,k-z)
      common/gravity/gs14,compactness
      save print

      if (print.eq.1.) goto 10
       print *,'-------------------------------------------'
       print *,'Using Gundmundsson et al boundary condition'
       print *,'-------------------------------------------'
       print=1.
 10   continue

      Tb8=Tb/1.e8
      fteff_GPE=0.87e6*gs14**.25*(Tb8)**0.55
     
      end

c *********************************************************************
c *********************************************************************

      function fteff_acc(Tb,etabound)

c this is the boundary condition for accreted envelope
c From Potekhin, Chabrier & Yakovlev, A&A 323, 415 (1999)
c CHECKED on MARCH 13, 2001
c WARNING : Tb must be the local value, not the redshifted one !
     
      implicit real*8 (a-h,k-z)
      common/gravity/gs14,compactness
      save print
c*****
      if (print.eq.1.) goto 10
       print *,'------------------------------------------'
       print *,'Using accreted envelope boundary condition'
       print '(a15,1p1e10.2)','with   eta =',etabound
       print *,'------------------------------------------'
       print=1.
 10   continue
c*****
      write(*,*)"Are you sure you want to use this one?"
      write(*,*)"Better use fteff_ANS"

      Tb9=Tb/1.d9
      Ts=sqrt(7.0d0*Tb9*sqrt(gs14))
      z=Tb9-Ts/1.d3
      t4_iron=gs14*((7.0d0*z)**2.25+(z/3.0d0)**1.25)
      t4_wacc=gs14*(18.1d0*Tb9)**2.42
      if (etabound.gt.1.d-30) then
       a=(1.2d0+(5.3d-6/etabound)**0.38)*Tb9**(5./3.)
       t4_acc=(a*t4_iron+t4_wacc)/(a+1.0d0)
      else
       t4_acc=t4_iron
      end if

      fteff_acc=t4_acc**0.25*1.d6
     
      end

c *********************************************************************
c *********************************************************************

      function fteff_fullacc(Tb)

c this is the boundary condition for a fully accreted envelope
c From Potekhin et al. A&A 594 (2003)
c Equation A2
c WARNING : Tb must be the local value, not the redshifted one !
     
      implicit real*8 (a-h,k-z)
      common/gravity/gs14,compactness
      save print
c*****
      if (print.eq.1.) goto 10
       print *,'------------------------------------------'
       print *,'Using fully accreted envelope boundary condition'
       print *,'------------------------------------------'
       print=1.
 10   continue
c*****
      Tb9=Tb/1.d9
      zeta=Tb9-1d-3*gs14**0.25*(7.*Tb9)**0.5
      T4_fe=gs14*((7.*zeta)**2.25+(0.33*zeta)**1.25)

      T4_acc=(gs14*(18.1d0*Tb9)**2.42*
     1 (0.447+0.075*log10(Tb)/(1.+(6.2*Tb9)**4.))+3.2*Tb9**1.67*T4_fe)
     1 /(1.+3.2*Tb9**1.67)

      fteff_fullacc=t4_acc**0.25*1.d6
     
      end

c *********************************************************************
c *********************************************************************

      function fteff_ANS(Tb,etabound)

c this is the boundary condition for a partly accreted envelope
c From Potekhin et al. ApJ 594 (2003)
c Equations A7
c WARNING : Tb must be the local value, not the redshifted one !
     
      implicit real*8 (a-h,k-z)
      common/gravity/gs14,compactness
      save print
c*****
      if (print.eq.1.) goto 10
       print *,'------------------------------------------'
       print *,'Using partly accreted envelope boundary condition'
       print *,'Delta M/M=',DMM/gs14**2.
       print *,'------------------------------------------'
       print=1.
 10   continue
c*****
      Tb9=Tb/1.d9
      zeta=Tb9-1d-3*gs14**0.25*(7.*Tb9)**0.5
      T4_fe=gs14*((7.*zeta)**2.25+(0.33*zeta)**1.25)
      T4_acc=(gs14*(18.1d0*Tb9)**2.42*
     1 (0.447+0.075*log10(Tb)/(1.+(6.2*Tb9)**4.))+3.2*Tb9**1.67*T4_fe)
     1 /(1.+3.2*Tb9**1.67)

      xhi=dlog10(1e6*etabound)*(-1)
      gamma=1./(1.+3.8*(0.1*xhi)**9.)/(1.+0.171*xhi**(7./2.)*Tb9)
      
      fteff_ANS=(gamma*t4_acc+(1.-gamma)*t4_fe)**0.25*1e6
      end

c *********************************************************************
c *********************************************************************

c *********************************************************************
c *********************************************************************

      function fteff_iron(Tb)

c this is the boundary condition for an iron envelope envelope
c From Potekhin et al. ApJ 594 (2003)
c Equation A1
c WARNING : Tb must be the local value, not the redshifted one !
     
      implicit real*8 (a-h,k-z)
      common/gravity/gs14,compactness
      save print
c*****
      if (print.eq.1.) goto 10
       print *,'------------------------------------------'
       print *,'Using an iron envelope boundary condition'
       print *,'------------------------------------------'
       print=1.
 10   continue
c*****

      Tb9=Tb/1.d9
      zeta=Tb9-1d-3*gs14**0.25*(7.*Tb9)**0.5
      t4_fe=gs14*((7.*zeta)**2.25+(1./3.*zeta)**1.25)
      fteff_iron=t4_fe**0.25*1.d6
      end


c *********************************************************************
c *********************************************************************

      function fteff_field_iron(Tb,bfield)

c this is the boundary condition for iron envelope with B field
c From Potekhin & Yakovlev, A&A 374 (2001), p. 213
c CHECKED on June 5, 2002 
c WARNING : Tb must be the local value, not the redshifted one !
c WARNING: bfield is B at the magnetic pole
c WARNING: good only for Tb > 1e7 K
     
      implicit real*8 (a-h,k-z)
      common/gravity/gs14,compactness
      save print
c*****
      if (print.eq.1.) goto 10
       print *,'------------------------------------------'
       print *,'Using iron envelope with magnetic field:'
       print *,'(a15,1p1e10.2)','B_pole =',bfield
       print *,'------------------------------------------'
       print=1.
 10   continue
c*****
      etabound=0.d0
      f_zero=fteff_acc(Tb,etabound)

      B12=bfield/1.d12
      T9=Tb/1.d9
      beta=0.074d0*dsqrt(B12)/T9**0.45
      a1=5059.d0*T9**0.75/
     1   (1.d0+20.4d0*T9**0.5+138.d0*T9**1.5+1102.d0*T9**2)**0.5
      a2=1484.d0*T9**0.75/
     1   (1.d0+90.d0*T9**1.5+125.d0*T9**2)**0.5
      a3=5530.d0/dsqrt(1.d0-0.4d0*compactness)*T9**0.75/
     1   (1.d0+8.16d0*T9**0.5+107.8d0*T9**1.5+560.d0*T9**2)**0.5

      ratio=(1.d0+a1*beta**2+a2*beta**3+0.007d0*a3*beta**4)/
     1      (1.d0+a3*beta**2)

      fteff_field_iron=f_zero*ratio**0.25
     
      end

c *********************************************************************
c *********************************************************************

      function fteff_HZ90(Tb,etabound)

c this is the boundary condition for a partly accreted envelope
c From Yakovlev et al. AA 417 (2004)
c Equations 5a-d
c WARNING : Tb must be the local value, not the redshifted one !
     
      implicit real*8 (a-h,k-z)
      common/gravity/gs14,compactness
      save print
      dimension ahz90(7)
      data ahz90/2.87,0.534,2.203,1.349,8.328,0.01609,0.1378/

c*****
      if (print.eq.1.) goto 10
       print *,'------------------------------------------'
       print *,'Using accreted envelope boundary condition HZ90'
       print *,'------------------------------------------'
       print=1.
 10   continue
c*****
      Tb9=Tb/1.d9
      T_non=gs14**0.25*ahz90(1)*Tb9**ahz90(2)
 
      cfac=(1.+ahz90(3)*Tb9**ahz90(4))/(1.+ahz90(5)*Tb9**1.6)
      T_full=6.9*gs14**0.25*Cfac**0.25*Tb9**0.62

      p=-ahz90(6)/etabound**ahz90(7)
      gamma=(1.+120.*Tb9)**p
      
      fteff_HZ90=(gamma*t_full**4.+(1.-gamma)*t_non**4.)**0.25*1e6
      end

C************************************************************************
C************************************************************************

C************************************************************************
C Models of two-components envelope
C Beznogov et al. MNRAS 2016
C Ts(Tb,deltaM/M)
C************************************************************************
      function fteff_Bez_DMM(Tb, deltaMM,icompo)
      implicit real*8 (a-h,k-z)
      common/gravity/gs14,compactness
      PARAMETER (imax=30)
      real*8 lTst(imax),lTbt(imax)

      lTb=dlog10(Tb)
      do i=1,imax
      lTst(i)=4.+(i-1.)/(imax-1)*(8.-4.)
      lTbt(i)=dlog10(Tb_Bez_DMM(10**lTst(i), deltaMM,icompo))
      enddo

      call polint(lTbt,lTst,imax,lTb,lTs,dy)
      write(90,*)lTb,lTs,dy
      fteff_Bez_DMM=10.**lTs
      end
C************************************************************************

C************************************************************************
C Tb(Ts,deltaM/M)
C************************************************************************

      function Tb_Bez_DMM(Ts, deltaMM,icompo)
      implicit real*8 (a-h,k-z)
      common/gravity/gs14,compactness

      rhostar=rho_DMM(deltaMM,icompo)
      Tb_Bez_DMM=Tb_Bez(Ts, rhostar,icompo)
      end
C************************************************************************


C************************************************************************
C Equation A1 Tb(Ts,rho*)
C************************************************************************
      function Tb_Bez(Ts, rhostar,icompo)
      implicit real*8 (a-h,k-z)
      dimension p1(14),p2(11),p3(14)
      common/gravity/gs14,compactness
      Ts6=Ts/1d6
      Y=Ts6*(2.4271d0/gs14)**0.25d0

      Ymin=0.32d0
      Ymax=2.9d0

c H/He
      if (icompo.eq.1) then
      data p1/3.150d0,1.546d0,0.3225d0,1.132d0,1.621d0,
     1       1.083d0,7.734d0,1.894d0,2.335d5,7.071d0,
     1       5.202d0,10.01d0,2.007d0,0.4703d0       /
      fact=1.75d0
      
      if ((Y.ge.Ymin).and.(Y.le.Ymax)) then
      f1=p1(1)*Y**p1(2)*dsqrt(1.0d0+p1(3)*Y**p1(4))	
      f4=p1(5)*Y**p1(6)*dsqrt(1.0d0+p1(7)*Y**p1(8))	
      f2=p1(9)*Y**p1(10)/(1.0d0-p1(11)*Y+p1(12)*Y*Y)**2.d0
      f3=p1(13)*(Y**(-1.0d0*p1(14)))
      f5=-0.3d0   
      Tb7=f4+(f1-f4)*(1.0d0+(rhostar/f2)**f3)**f5

      else if (Y.lt.Ymin) then
      f1=p1(1)*Ymin**p1(2)*dsqrt(1.0d0+p1(3)*Ymin**p1(4))	
      f4=p1(5)*Ymin**p1(6)*dsqrt(1.0d0+p1(7)*Ymin**p1(8))	
      f2=p1(9)*Ymin**p1(10)/(1.0d0-p1(11)*Ymin+p1(12)*Ymin*Ymin)**2.d0
      f3=p1(13)*(Ymin**(-1.0*p1(14)))
      f5=-0.3d0   
      Tb7=(f4+(f1-f4)*(1.0+(rhostar/f2)**f3)**f5)
     1     *10.d0**(fact*(dlog10(Y)-dlog10(Ymin)))

      else

      f1=p1(1)*Ymax**p1(2)*dsqrt(1.0d0+p1(3)*Ymax**p1(4))	
      f4=p1(5)*Ymax**p1(6)*dsqrt(1.0d0+p1(7)*Ymax**p1(8))	
      f2=p1(9)*Ymax**p1(10)/(1.0d0-p1(11)*Ymax+p1(12)*Ymax*Ymax)**2.d0
      f3=p1(13)*(Ymax**(-1.0d0*p1(14)))
      f5=-0.3d0   
      
      Tb7=(f4+(f1-f4)*(1.0d0+(rhostar/f2)**f3)**f5)
     1     *10.d0**(fact*(dlog10(Y)-dlog10(Ymax)))

      endif

c He/C
      else if (icompo.eq.2) then
      data p2/5.386d0,0.1027d0,1.719d0,3.872d0,0.1344d0,
     1      1.759d0,1.056d5,1.881d0,3.680d0,3.857d0,
     1      1.102d0/
      fact=1.8d0

      if ((Y.ge.Ymin).and.(Y.le.Ymax)) then
      f1=p2(1)*Y**(p2(2)*dlog10(Y)+p2(3))
      f4=p2(4)*Y**(p2(5)*dlog10(Y)+p2(6))
      f2=p2(7)*Y**(p2(8)*dlog10(Y)*dlog10(Y)+p2(9))
      f3=p2(10)*dsqrt(Y/(Y*Y+p2(11)*p2(11)))
      f5=-0.2d0  
      Tb7=f4+(f1-f4)*(1.0d0+(rhostar/f2)**f3)**f5

      else if (Y.lt.Ymin) then
      f1=p2(1)*Ymin**(p2(2)*dlog10(Ymin)+p2(3))
      f4=p2(4)*Ymin**(p2(5)*dlog10(Ymin)+p2(6))
      f2=p2(7)*Ymin**(p2(8)*dlog10(Ymin)*dlog10(Ymin)+p2(9))
      f3=p2(10)*dsqrt(Ymin/(Ymin*Ymin+p2(11)*p2(11)))
      f5=-0.2d0  
      Tb7=(f4+(f1-f4)*(1.0d0+(rhostar/f2)**f3)**f5)
     1     *10.d0**(fact*(dlog10(Y)-dlog10(Ymin)))

      else
      f1=p2(1)*Ymax**(p2(2)*dlog10(Ymax)+p2(3))
      f4=p2(4)*Ymax**(p2(5)*dlog10(Ymax)+p2(6))
      f2=p2(7)*Ymax**(p2(8)*dlog10(Ymax)*dlog10(Ymax)+p2(9))
      f3=p2(10)*dsqrt(Ymax/(Ymax*Ymax+p2(11)*p2(11)))
      f5=-0.2d0   
      Tb7=(f4+(f1-f4)*(1.0d0+(rhostar/f2)**f3)**f5)
     1     *10.d0**(fact*(dlog10(Y)-dlog10(Ymax)))

      endif

c C/Fe
      else if (icompo.eq.3) then
      data p3/0.175,0.4004,53.91,1.96,5.208,
     1      1.651,0.03245,0.005433,36570.,1.824,
     1      3.942,2.035,1.861,0.02637/
      fact=1.8d0

      if ((Y.ge.Ymin).and.(Y.le.Ymax)) then
      f1=p3(1)*Y**(p3(2)*(-1.))*
     1   (p3(3)*Y*Y+p3(4)*Y*Y*Y*Y-1.)
      f4=p3(5)*Y**p3(6)*
     1   (1.+p3(7)*Y*Y-p3(8)*Y*Y*Y*Y)
      f2=p3(9)*Y**(p3(11)-p3(10)*dlog10(Y)*dlog10(Y))
      f3=p3(12)*dsqrt(1./(Y*Y+p3(13)*p3(13)))
     1   *(1.-p3(14)*Y*Y)
      f5=-0.4   
      Tb7=(f4+(f1-f4)*(1.0+(rhostar/f2)**f3)**f5)

      else if (Y.lt.Ymin) then
      f1=p3(1)*Ymin**(p3(2)*(-1.))*
     1  (p3(3)*Ymin*Ymin+p3(4)*Ymin*Ymin*Ymin*Ymin-1.)
      f4=p3(5)*Ymin**p3(6)*
     1  (1.+p3(7)*Ymin*Ymin-p3(8)*Ymin*Ymin*Ymin*Ymin)
      f2=p3(9)*Ymin**(p3(11)-p3(10)*dlog10(Ymin)*dlog10(Ymin))
      f3=p3(12)*dsqrt(1./(Ymin*Ymin+p3(13)*p3(13)))
     1  *(1.-p3(14)*Ymin*Ymin)
      f5=-0.4   
      Tb7=(f4+(f1-f4)*(1.0+(rhostar/f2)**f3)**f5)
     1     *10.**(fact*(dlog10(Y)-dlog10(Ymin)))

      else
      f1=p3(1)*Ymax**(p3(2)*(-1.))*
     1  (p3(3)*Ymax*Ymax+p3(4)*Ymax*Ymax*Ymax*Ymax-1.)
      f4=p3(5)*Ymax**p3(6)*
     1  (1.+p3(7)*Ymax*Ymax-p3(8)*Ymax*Ymax*Ymax*Ymax)
      f2=p3(9)*Ymax**(p3(11)-p3(10)*dlog10(Ymax)*dlog10(Ymax))
      f3=p3(12)*dsqrt(1./(Ymax*Ymax+p3(13)*p3(13)))*
     1  (1.-p3(14)*Ymax*Ymax)
      f5=-0.4    
      Tb7=(f4+(f1-f4)*(1.0+(rhostar/f2)**f3)**f5)
     1     *10.**(fact*(dlog10(Y)-dlog10(Ymax)))

      endif

C Houston we have a problem      
      else
       write(*,*)'icompo out of range',icompo

      endif

      Tb_Bez=Tb7*1d7

      end
C************************************************************************


C************************************************************************
C Equation 17
C************************************************************************
      function DMM_rho(rhostar,icompo)
      implicit real*8 (a-h,k-z)
      common/gravity/gs14,compactness
      if (icompo.eq.1) then
      ZAratio=1.d0
      else if ((icompo.eq.2).or.(icompo.eq.3)) then
      ZAratio=0.5d0
      endif
      ksi = 0.01009*(rhostar*ZAratio)**(1.0/3.0)
      DeltaMM1=ksi*dsqrt(1.0+ksi**2.)*(2./3*ksi**2.-1.0)+
     2 dlog(ksi+dsqrt(1.0+ksi**2.))
      DMM_rho=DeltaMM1*1.51e-11/gs14**2.
      end
C************************************************************************

C************************************************************************
c Inversion of equation 17
C************************************************************************
      function rho_DMM(DeltaMM,icompo)
      implicit real*8 (a-h,k-z)
      common/gravity/gs14,compactness

      eps = 1.0e-6
      left =1.0
      right =4.0e9

      rho = 0.5*(left+right)
      DL = DMM_rho(left,icompo)
      DR = DMM_rho(right,icompo)
      if (DeltaMM .gt. DR) then
       write(*,*)"DeltaMM_Rho_eff error: too high DeltaMM value."
      else if (DeltaMM .lt. DL) then
       write(*,*)"DeltaMM_Rho_eff error: too low DeltaMM value."
      else
       DMM =DMM_rho(rho,icompo)
       do while (((2.0*(right-left)/(right+left)).gt.eps).and.
     2   (abs((DMM-DeltaMM)/DeltaMM).gt.eps)) 
        if (DMM > DeltaMM) then
         right = 0.5*(left+right)
        else
 	 left = 0.5*(left+right)
        endif
       rho = 0.5*(left+right)
       DMM = DMM_rho(rho,icompo)
       enddo
      endif 
      rho_DMM=rho
      end
C************************************************************************



C************************************************************************
      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      DOUBLE PRECISION dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=30)
      INTEGER i,m,ns
      DOUBLE PRECISION den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.d0)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software *!^3]87J,.
C************************************************************************




C************************************************************************

      SUBROUTINE SPLINT1(XA,YA,Y2A,IN,X,Y)

C************************************************************************
C Given the arrays XA and YA of length IN, which tabulate a function    *
C (with the XA(i)'s in order), and given the array Y2A, which is the    *
C output from SPLINE, calculate the cubic-spline interpolated value     *
C Y for a given value X.                                                *
C                                                                       *
C From NUMERICAL RECIPES, p.89.                                         *
C************************************************************************

      IMPLICIT REAL*8 (A-H,L-Z)
      DIMENSION XA(IN),YA(IN),Y2A(IN)

      KLO=1
      KHI=IN
1     IF (KHI-KLO.GT.1) THEN
       K=(KHI+KLO)/2
        IF(XA(K).GT.X) THEN
         KHI=K
        ELSE
         KLO=K
        END IF
       GOTO 1
      END IF
      
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.) PAUSE 'Bad XA input'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO) + B*YA(KHI)
     1  + ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.

      RETURN

      END

C************************************************************************
C************************************************************************
C************************************************************************

      SUBROUTINE SPLINE1(X,Y,IN,YP1,YPN,Y2)

C************************************************************************
C Given the arrays X and Y of length IN, which tabulate a function      *
C (with the XA(i)'s in order), and given values YP1 and YPN of the      *
C first derivative of the interpolating function at points 1 and N,     *
C this subroutine returns an array Y2 of length IN which contains the   *
C second derivatives of the interpolating function at the tabulated     *
C points X(i)'s.                                                        *
C                                                                       *
C From NUMERICAL RECIPES, p.88                                          *
C************************************************************************

      IMPLICIT REAL*8 (A-H,L-Z)
      PARAMETER (JMAX=100)
      DIMENSION X(IN),Y(IN),Y2(IN),U(JMAX)

      IF (YP1.GT.1.E30) THEN
       Y2(1)=0.
       U(1)=0.
      ELSE
       Y2(1)=-0.5
       U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO I=2,IN-1
       SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
       P=SIG*Y2(I-1)+2.
       Y2(I)=(SIG-1.)/P
       U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     1      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
      END DO

      IF (YPN.GT.1.E30) THEN
       QN=0.
       UN=0.
      ELSE
       QN=0.5
       UN=(3./(X(IN)-X(IN-1)))*(YPN-(Y(IN)-Y(IN-1))/(X(IN)-X(IN-1)))
      END IF
      Y2(IN)=(UN-QN*U(IN-1))/(QN*Y2(IN-1)+1.)
      
      DO K=IN-1,1,-1
       Y2(K)=Y2(K)*Y2(K+1)+U(K)
      END DO

      RETURN
      END   

