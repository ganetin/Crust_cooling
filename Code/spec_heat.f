      subroutine specheat_hfb (i,t,rho,acell,aion,zion,cv,
     1       cvneutron,cvproton,cvelectron,cvmuon,
     2       cvlambda,cvsminus,cvszero,cvsplus,
     3       cvquark,cvions)


c***********************************************************************
c  Includes the Levenfish-Yakovlev coefficients for pairing suppression
c***********************************************************************

      implicit real*8(a-h,k-z)
      INCLUDE 'size.inc.f'
      INCLUDE 'rho_limits.inc.f'
c      parameter (pi = 3.14159265)
        
      INCLUDE 'profile_star.inc.f'
      INCLUDE 'profile_comp.inc.f'
      INCLUDE 'pairing.inc.f'
      INCLUDE 'pairing_quark_PU2002.inc.f'
      INCLUDE 'spec_heat.inc.f'
      INCLUDE 'fermi.inc.f'
      INCLUDE 'Physical_Constants.inc.f'
      REAL*8 AY(1),AZION(1),ACMI(1)
      INTEGER*4 LIQSOL
c ***************************

      fexp(x)=dexp(max(x,-7.d2))


      u_1s0(t)=dsqrt(1.-t)*(1.456-0.157/dsqrt(t)+1.764/t)
      r_1s0(u)=(0.4186+dsqrt(1.007**2+(0.5010*u)**2))**2.5*
     1         dexp(1.456-dsqrt(1.456**2+u**2))
      u_3p2(t)=dsqrt(1.-t)*(5.596+8.424/t)
      r_3p2(u)=(0.6893+dsqrt(0.790**2+(0.03983*u)**2))**2*
     1         fexp(1.934-dsqrt(1.934**2+u**2/(16.*pi)))

c ***************************
c ****** get Cv-ions :
     
      if (rho .lt. rhocore) then
	if ((icvel_nodeg.eq.1).and.(rho.lt.1e10)) then
	cvions=0.d0
	else
       call cvion(t,rho,acell,aion,zion,cv_ions)
       cv_nrot=0.d0
       cv_nvib=0.d0
       cvions=cv_ions+cv_nrot+cv_nvib
	endif
      else
       cvions=0.
      end if
c ****** get Cv-electrons :

c      if ((rho .lt. rhodrip).and.(icvel_nodeg.eq.1)) then
c       call cvelec(t,rho,acell,aion,zion,i,cvelectron)
	if ((icvel_nodeg.eq.1).and.(rho.lt.1e10)) then
	cvelectron=0.d0
      else
       cvelectron=cve(i)*t
      end if

c ***** get Cv-muons :
      cvmuon=cvm(i)*t


!c ****** get Cv-neutrons :

!c HERE DANY cccccccccccccccccccccccc
!c       raise=1.1d0
!       raise=1.0d0
!cccccccccccccccccccccccccccccccccccc
!      kboltz=1.380d-16
!      MeV2erg=1.602177d-6

!      if ((rho.ge.rhodrip).and.(rho.le.rhocore).and.(istrange.eq.0)) 
!     1 then
!            raise=1.0d0
!            if (t .gt. 20*tcn(i)) then
!               r=0.d0  
!               t0=1.d-10
!            else if (t .gt. tcn(i)) then
!                u=0.001
!                r=r_1s0(u)
!                t0=1.d-10
!            else
!                t0=t/tcn(i)
!                if (i.le.isf) then
!                     u=u_3p2(t0)
!                     r=r_3p2(u)
!                 else
!                     u=u_1s0(t0)
!                     r=r_1s0(u)
!                 endif
!            endif
!cc Include the new specific heat
!cc density of free neutrons (at t=0)
!            nn0=kfn(i)**3/(3.*pi**2)
!	        if (acell.eq.0.d0) then
!	          nnm=bar(i)
!	        else
!		      nnm=(acell-zion)/acell*bar(i)
!		    endif
!            if ((kbol*t).gt.5.5d0) then
!              nn=nnm
!            else
!              nn=nn0+kbol*t/5.5*(nnm-nn0)	
!            endif
!c            write(*,*)i,nn0,nnm,acell,zion,bar(i)
!c            stop
!            nn=nn0
!cc calculate the FErmi energy for the free neutrons
!            mn=939.56d0*mstn(i)
!            kfnn=(3.*pi**2*nn)**0.33333
!c            enfern=0.5*hbc**2/mn*kfn(i)**2
!            enfern=0.5*hbc**2/mn*kfnn**2
!cc temperature effect in cvn (beyond the linear term)
!            vartemp=pi*t*kbol/enfern
!            fcvn=1.-7./40.*vartemp**2-155./896.*vartemp**4
!cc classial neutron specific heat in erg.cm-3.K-1
!            cvnclass = 1.5*nn*kbol*MeV*1.d39
!cc smooth transition from classical gas to quantal gas            
!             if (vartemp.lt.0.5d0) then
!                  xop=1.
!             else if (vartemp.gt.1.5d0) then
!                  xop=0.d0
!             else
!          	   xop=1./(1.+dexp(+5.*(vartemp-1.)))
!             endif
!c          write(*,*)'transition quantique classique prise en compte'
!            if ((sfn1s0.eq.301.).or.(sfn1s0.eq.302.)) then
!cc modification of the function R
!               hfb_f1=(1.d0+dexp(-hfb_a0(i)*hfb_a1(i)*hfb_d0(i)
!     1                /hfb_a3(i)))/(1.d0+dexp((t/(MeV2erg/kboltz)
!     2                -hfb_a0(i)*hfb_a1(i)*hfb_d0(i))/hfb_a3(i)))
!               hfb_f2=(1.d0+dexp(-hfb_a0(i)*hfb_a2(i)*hfb_d0(i)/
!     1                hfb_a3(i)))/(1.d0+dexp((t/(MeV2erg/kboltz)
!     2                -hfb_a0(i)*hfb_a2(i)*hfb_d0(i))/hfb_a3(i)))
!               hfb_r=r*hfb_f1+(1-hfb_f2)
!c		     write(*,'(i4,2d14.4)')i,cvn(i),cvnclass
!c		     stop
!cc calculate the specific heat
!               cvneutron =cvn(i)*hfb_r*t*fcvn*xop+cvnclass*(1.-xop)
!c		     write(*,'(i4,3d14.4)')i,cvn(i)*t,cvnclass,cvneutron
!             else if (sfn1s0.eq.300.) then
!               cvneutron =cvn(i)*t*fcvn*xop+cvnclass*(1.-xop)
!             else
!               cvneutron =cvn(i)*t*r       
!             endif
!cc end of new specific heat
!      else if ((rho.lt.rhodrip).and.(istrange.eq.0)) then
!             cvneutron=0.d0     
!      else if ((rho.gt.rhocore).and.(istrange.eq.0)) then
!            if (t .lt. tcn(i))then
!             t0=min(0.999999999999d0,t/tcn(i))
!                  if (i.le.isf) then
!                   u=u_3p2(t0)
!                   r=r_3p2(u)
!                  else
!                   u=u_1s0(t0)
!                   r=r_1s0(u)
!                  end if
!            else
!             r=1.d0
!            end if
!            cvneutron =cvn(i)*t*r
!c corresponds to Shternin et al. 2007 with normal neutrons everywhere
!       if ((sfn1s0.eq.400).and.(sfn3p2.eq.400)) then
!               cvneutron =cvn(i)*t
!       endif
!c corresponds to Shternin et al. 2007 with normal neutrons in the core
!       if (sfn1s0.eq.401.and.(sfn3p2.eq.401)) then
!          if (rho.gt.rhocore) then
!               cvneutron =cvn(i)*t
!          endif
!       endif
!        end if
!      cvneutron=cvn(i)*t

!c here cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Milano specific heat:
!c      if ((rho.ge.rhodrip).and.(i.gt.isf)) then
!c       call cvn_milano(t,rho,cvneutron)
!c      end if
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c ****** get Cv-neutrons :

      if ((rho.ge.rhodrip).and.(istrange.eq.0)) then
c HERE DANY cccccccccccccccccccccccc
c       raise=1.1d0
       raise=1.0d0
cccccccccccccccccccccccccccccccccccc
       if (t .lt. raise*tcn(i))then
        t0=min(0.999999999999d0,t/tcn(i))
        if (i.le.isf) then
         u=u_3p2(t0)
         r=r_3p2(u)
        else
         u=u_1s0(t0)
         r=r_1s0(u)
        end if
        if (t .gt. tcn(i)) then
         w1=(raise*tcn(i)-t)/((raise-1.)*tcn(i))
         w2=1.-w1
         r=w1*r+w2*1.0
        end if
       else
        r=1.d0
       end if
       cvneutron =cvn(i)*t*r
      else
       cvneutron=0.d0
      end if
c here cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Milano specific heat:
c      if ((rho.ge.rhodrip).and.(i.gt.isf)) then
c       call cvn_milano(t,rho,cvneutron)
c      end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c ****** get Cv-protons :
  
      if ((rho.ge.rhocore).and.(istrange.eq.0)) then
       raise=1.1d0
       if(t .lt. raise*tcp(i))then
        t0=min(0.999999999999d0,t/tcp(i))
        u=u_1s0(t0)
        r=r_1s0(u)
        if (t .gt. tcp(i)) then
         w1=(raise*tcp(i)-t)/((raise-1.)*tcp(i))
         w2=1.-w1
         r=w1*r+w2*1.0
        end if
       else
        r=1.d0
       end if
       cvproton  =cvp(i)*t*r
c corresponds to Shternin et al. 2007 with normal protons everywhere
       if ((sfp1s0.eq.400.).or.(sfp1s0.eq.401.)) then
               cvproton =cvp(i)*t
       endif
      else
       cvproton=0.d0
      end if

c ****** get Cv-lambdas :
  
      if ((rho.ge.rhocore).and.(istrange.eq.0)) then
       raise=1.1
       if(t .lt. raise*tcla(i))then
        t0=min(0.999999999999d0,t/tcla(i))
        u=u_1s0(t0)
        r=r_1s0(u)
        if (t .gt. tcla(i)) then
         w1=(raise*tcla(i)-t)/((raise-1.)*tcla(i))
         w2=1.-w1
         r=w1*r+w2*1.0
        end if
       else
        r=1.d0
       end if
       cvlambda  =cvla(i)*t*r
      else
       cvlambda=0.d0
      end if

c ****** get Cv-Sigma- :
  
      if ((rho.ge.rhocore).and.(istrange.eq.0)) then
       raise=1.1
       if(t .lt. raise*tcsm(i))then
        t0=min(0.999999999999d0,t/tcsm(i))
        u=u_1s0(t0)
        r=r_1s0(u)
        if (t .gt. tcsm(i)) then
         w1=(raise*tcsm(i)-t)/((raise-1.)*tcsm(i))
         w2=1.-w1
         r=w1*r+w2*1.0
        end if
       else
        r=1.d0
       end if
       cvsminus  =cvsm(i)*t*r
      else
       cvsminus=0.d0
      end if

c ****** get Cv-Sigma0 :
  
      if ((rho.ge.rhocore).and.(istrange.eq.0)) then
       raise=1.1
       if(t .lt. raise*tcs0(i))then
        t0=min(0.999999999999d0,t/tcs0(i))
        u=u_1s0(t0)
        r=r_1s0(u)
        if (t .gt. tcs0(i)) then
         w1=(raise*tcs0(i)-t)/((raise-1.)*tcs0(i))
         w2=1.-w1
         r=w1*r+w2*1.0
        end if
       else
        r=1.d0
       end if
       cvszero  =cvs0(i)*t*r
      else
       cvszero=0.d0
      end if

c ****** get Cv-Sigma+ :
  
      if ((rho.ge.rhocore).and.(istrange.eq.0)) then
       raise=1.1
       if(t .lt. raise*tcsp(i))then
        t0=min(0.999999999999d0,t/tcsp(i))
        u=u_1s0(t0)
        r=r_1s0(u)
        if (t .gt. tcsp(i)) then
         w1=(raise*tcsp(i)-t)/((raise-1.)*tcsp(i))
         w2=1.-w1
         r=w1*r+w2*1.0
        end if
       else
        r=1.d0
       end if
       cvsplus  =cvsp(i)*t*r
      else
       cvsplus=0.d0
      end if

c ****** get Cv-photons :

c      cvphot= 4.*7.56e-15*t**3
      cvphot=0.d0

c ***** get Cv-quarks:

      if (rho.ge.rhocore) then
       if (c_cv_str.le.1.0) then
        call specheat_quark(i,t,cvquark)
       else
        cvquark=c_cv_str*(T/1.d9)
       end if
      else
       cvquark=0.d0
      end if

c ***** get Cv-EIP when rho.le.1e10
	if ((icvel_nodeg.eq.1).and.(rho.lt.1e10)) then
	iNMIX=1
	AY(1)=1.
	AZion(1)=zion
	ACMI(1)=aion
	TEMP=T/3.1577e5
	iKRAD=0
	call MELANGE8(iNMIX,AY,AZion,ACMI,RHO,TEMP,iKRAD,
     *   DENS,Zmean,CMImean,Z2mean,GAMImean,CHI,TPT,LIQSOL,
     *   PnkT,UNkT,SNk,CV,CHIR,CHIT)	
	cv=cv*kb*Na*rho/Aion
	else
c ****** Total Cv :
      cvelectron=cvelectron*fhad(i)
      cvmuon    =cvmuon    *fhad(i)
      cvproton  =cvproton  *fhad(i)
      cvneutron =cvneutron *fhad(i)
      cvlambda  =cvlambda  *fhad(i)
      cvsminus  =cvsminus  *fhad(i)
      cvszero   =cvszero   *fhad(i)
      cvsplus   =cvsplus   *fhad(i)
      cvquark   =cvquark*(1.d0-fhad(i))

      cv=cvions+
     2   cvelectron+cvmuon+
     3   cvproton+cvneutron+
     4   cvlambda+cvsminus+cvszero+cvsplus+
     5   cvphot+cvquark

	endif

c Ofengeim
!      if (rcvcrust.ne.1.) then
!      if (rho.lt.rhocore) cv=cv*rcvcr
!      if (rho.ge.rhocore) then
!       c0cv=1.12e20
!       acv=2.68
!       bcv=1.14
!       ccv=0.0170
!       gamcv=2.11
!       funcJ=xrho**(1./3.)*(1.-compactness)**(-0.5)*
!     1       (1.+ccv*xrho**(gamcv-1.))**((gamcv-1./3.)/(gamcv-1.))/
!     1       (1.-bcv*compactness)**(0.5)
!       cvoftot=c0cv*acv*radius**3.*(T/1e9)*funcJ !-> OK. checked vs fig 4
!       cv=cvoftot/dvcore
!       write(*,*)'cvoftot',cvoftot
!      endif
!      endif
      
      end


c***********************************************************************
c                        *Original routine*
c***********************************************************************
      subroutine specheat (i,t,rho,aion,zion,cv,
     1       cvneutron,cvproton,cvelectron,cvmuon,
     2       cvlambda,cvsminus,cvszero,cvsplus,
     3       cvquark,cvions)

c ****** checked on March 30, 1993

c***********************************************************************
c  Includes the Levenfish-Yakovlev coefficients for pairing suppression
c***********************************************************************

      implicit real*8(a-h,k-z)
      INCLUDE 'size.inc.f'
      INCLUDE 'rho_limits.inc.f'
      parameter (pi = 3.14159265)
        
      INCLUDE 'profile_star.inc.f'
      INCLUDE 'profile_comp.inc.f'
      INCLUDE 'pairing.inc.f'
      INCLUDE 'pairing_quark_PU2002.inc.f'
      INCLUDE 'spec_heat.inc.f'

c ***************************

      fexp(x)=dexp(max(x,-7.d2))
      u_1s0(t)=dsqrt(1.-t)*(1.456-0.157/dsqrt(t)+1.764/t)
      r_1s0(u)=(0.4186+dsqrt(1.007**2+(0.5010*u)**2))**2.5*
     1         fexp(1.456-dsqrt(1.456**2+u**2))

      u_3p2(t)=dsqrt(1.-t)*(5.596+8.424/t)
      r_3p2(u)=(0.6893+dsqrt(0.790**2+(0.03983*u)**2))**2*
     1         fexp(1.934-dsqrt(1.934**2+u**2/(16.*pi)))

c ***************************

c ****** get Cv-ions :
     
      if (rho .lt. rhocore) then
c       call cvion(t,rho,aion,zion,cv_ions)
       call cvion(t,rho,acell,aion,zion,cv_ions)

c       call cvnrot(t,rho,aion,zion,cv_nrot)
       cv_nrot=0.d0
c       call cvnvib(t,i,rho,aion,zion,cv_nvib)
       cv_nvib=0.d0
       cvions=cv_ions+cv_nrot+cv_nvib
      else
       cvions=0.
      end if

c ****** get Cv-electrons :

      if ((rho .lt. rhodrip).and.(icvel_nodeg.eq.1)) then
       call cvelec(t,rho,aion,zion,i,cvelectron)
      else
       cvelectron=cve(i)*t
      end if

c ***** get Cv-muons :

      cvmuon=cvm(i)*t

c ****** get Cv-neutrons :

      if ((rho.ge.rhodrip).and.(istrange.eq.0)) then
c HERE DANY cccccccccccccccccccccccc
c       raise=1.1d0
       raise=1.0d0
cccccccccccccccccccccccccccccccccccc
       if (t .lt. raise*tcn(i))then
        t0=min(0.999999999999d0,t/tcn(i))
        if (i.le.isf) then
         u=u_3p2(t0)
         r=r_3p2(u)
        else
         u=u_1s0(t0)
         r=r_1s0(u)
        end if
        if (t .gt. tcn(i)) then
         w1=(raise*tcn(i)-t)/((raise-1.)*tcn(i))
         w2=1.-w1
         r=w1*r+w2*1.0
        end if
       else
        r=1.d0
       end if
       cvneutron =cvn(i)*t*r
      else
       cvneutron=0.d0
      end if
c here cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Milano specific heat:
c      if ((rho.ge.rhodrip).and.(i.gt.isf)) then
c       call cvn_milano(t,rho,cvneutron)
c      end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c ****** get Cv-protons :
  
      if ((rho.ge.rhocore).and.(istrange.eq.0)) then
       raise=1.1d0
       if(t .lt. raise*tcp(i))then
        t0=min(0.999999999999d0,t/tcp(i))
        u=u_1s0(t0)
        r=r_1s0(u)
        if (t .gt. tcp(i)) then
         w1=(raise*tcp(i)-t)/((raise-1.)*tcp(i))
         w2=1.-w1
         r=w1*r+w2*1.0
        end if
       else
        r=1.d0
       end if
       cvproton  =cvp(i)*t*r
      else
       cvproton=0.d0
      end if

c ****** get Cv-lambdas :
  
      if ((rho.ge.rhocore).and.(istrange.eq.0)) then
       raise=1.1
       if(t .lt. raise*tcla(i))then
        t0=min(0.999999999999d0,t/tcla(i))
        u=u_1s0(t0)
        r=r_1s0(u)
        if (t .gt. tcla(i)) then
         w1=(raise*tcla(i)-t)/((raise-1.)*tcla(i))
         w2=1.-w1
         r=w1*r+w2*1.0
        end if
       else
        r=1.d0
       end if
       cvlambda  =cvla(i)*t*r
      else
       cvlambda=0.d0
      end if

c ****** get Cv-Sigma- :
  
      if ((rho.ge.rhocore).and.(istrange.eq.0)) then
       raise=1.1
       if(t .lt. raise*tcsm(i))then
        t0=min(0.999999999999d0,t/tcsm(i))
        u=u_1s0(t0)
        r=r_1s0(u)
        if (t .gt. tcsm(i)) then
         w1=(raise*tcsm(i)-t)/((raise-1.)*tcsm(i))
         w2=1.-w1
         r=w1*r+w2*1.0
        end if
       else
        r=1.d0
       end if
       cvsminus  =cvsm(i)*t*r
      else
       cvsminus=0.d0
      end if

c ****** get Cv-Sigma0 :
  
      if ((rho.ge.rhocore).and.(istrange.eq.0)) then
       raise=1.1
       if(t .lt. raise*tcs0(i))then
        t0=min(0.999999999999d0,t/tcs0(i))
        u=u_1s0(t0)
        r=r_1s0(u)
        if (t .gt. tcs0(i)) then
         w1=(raise*tcs0(i)-t)/((raise-1.)*tcs0(i))
         w2=1.-w1
         r=w1*r+w2*1.0
        end if
       else
        r=1.d0
       end if
       cvszero  =cvs0(i)*t*r
      else
       cvszero=0.d0
      end if

c ****** get Cv-Sigma+ :
  
      if ((rho.ge.rhocore).and.(istrange.eq.0)) then
       raise=1.1
       if(t .lt. raise*tcsp(i))then
        t0=min(0.999999999999d0,t/tcsp(i))
        u=u_1s0(t0)
        r=r_1s0(u)
        if (t .gt. tcsp(i)) then
         w1=(raise*tcsp(i)-t)/((raise-1.)*tcsp(i))
         w2=1.-w1
         r=w1*r+w2*1.0
        end if
       else
        r=1.d0
       end if
       cvsplus  =cvsp(i)*t*r
      else
       cvsplus=0.d0
      end if

c ****** get Cv-photons :

c      cvphot= 4.*7.56e-15*t**3
      cvphot=0.d0

c ***** get Cv-quarks:

      if (rho.ge.rhocore) then
       if (c_cv_str.le.1.0) then
        call specheat_quark(i,t,cvquark)
       else
        cvquark=c_cv_str*(T/1.d9)
       end if
      else
       cvquark=0.d0
      end if

c ****** Total Cv :
      cvelectron=cvelectron*fhad(i)
      cvmuon    =cvmuon    *fhad(i)
      cvproton  =cvproton  *fhad(i)
      cvneutron =cvneutron *fhad(i)
      cvlambda  =cvlambda  *fhad(i)
      cvsminus  =cvsminus  *fhad(i)
      cvszero   =cvszero   *fhad(i)
      cvsplus   =cvsplus   *fhad(i)
      cvquark   =cvquark*(1.d0-fhad(i))

      cv=cvions+
     2   cvelectron+cvmuon+
     3   cvproton+cvneutron+
     4   cvlambda+cvsminus+cvszero+cvsplus+
     5   cvphot+cvquark

      return

      end

c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
      subroutine specheat_quark(i,T,cv)
c     CHECKED on March 3, 2002
       implicit real*8(a-h,k-z)
       INCLUDE 'size.inc.f'
       INCLUDE 'profile_star.inc.f'
       INCLUDE 'profile_comp.inc.f'
       INCLUDE 'spec_heat.inc.f'
       INCLUDE 'pairing.inc.f'
       INCLUDE 'pairing_quark_PU2002.inc.f'

       parameter (pi=3.1415926535d0)
       parameter (arad=7.56d-15)
c ***************************
       fexp(x)=dexp(max(x,-7.d2))
       u_1s0(t)=dsqrt(1.-t)*(1.456-0.157/dsqrt(t)+1.764/t)
       r_1s0(u)=(0.4186+dsqrt(1.007**2+(0.5010*u)**2))**2.5*
     1          fexp(1.456-dsqrt(1.456**2+u**2))
c ***************************
        raise=1.2d0

        cv_e=cve(i) *T
        cv_u=cvqu(i)*T              ! WARNING:
        cv_d=cvqd(i)*T              ! must still multiply by 3
        cv_s=cvqs(i)*T              ! to take into account color 

c ***************************
        nbar=bar(i)
        r_CFL=1.d0
        r_3SC=1.d0
        r_2SC=1.d0
        r_u  =1.d0
        r_d  =1.d0
        r_s  =1.d0
c *******************************************************************
c NO PAIRING:
         if (normal.eq.1.d0) then
          cv = cv_e + 3.d0*(cv_u+cv_d+cv_s)
c *******************************************************************
c CFL PHASE:
         else if (nbar.ge.nb_CFL) then
          if (T.le.raise*Tc_CFL(i)) then
           tt=min(0.999999999999d0,T/Tc_CFL(i))
           u=u_1s0(tt)
           if      (type_CFL.eq.0 .d0) then
            r_CFL=r_1s0(u)
            if (T .ge. Tc_CFL(i)) then
             w1=(raise*Tc_CFL(i)-T)/((raise-1.d0)*Tc_CFL(i))
             w2=1.d0-w1
             r_CFL=w1*r_CFL+w2*1.d0
            end if
           else if (type_CFL.eq.1.d0) then 
            r_CFL=tt**2
           else if (type_CFL.eq.2.d0) then
            r_CFL=tt
           else 
            pause 'Spec_heat: type_CFL wrong !'
           end if
          end if
          cv_photons=3.d0*arad* T**3
          cv = 3.d0*r_CFL*(cv_u+cv_d+cv_s) + cv_photons
c          print '(1p2e20.3)',t,cv
c *******************************************************************
c 3SC PHASE:
         else if ((nbar.ge.nb_3SC).and.
     1            (nbar.lt.nb_CFL)) then 
          if (T.le.raise*Tc_3SC(i)) then
           tt=min(0.999999999999d0,T/Tc_3SC(i))
           u=u_1s0(tt)
           if      (type_3SC.eq.0.d0) then
            r_3SC=r_1s0(u)
            if (T .gt. Tc_3SC(i)) then
             w1=(raise*Tc_3SC(i)-T)/((raise-1.d0)*Tc_3SC(i))
             w2=1.d0-w1
             r_3SC=w1*r_3SC+w2*1.d0
            end if
           else if (type_3SC.eq.1.d0) then
            r_3SC=tt**2
           else if (type_3SC.eq.2.d0) then
            r_3SC=tt
           else
            pause 'Spec_heat: type_3SC wrong !'
           end if
          end if
          cv = cv_e + 3.d0*r_3SC*(cv_u+cv_d+cv_s)
c          print '(i5,0p1f8.5,1p2e12.3,5x,1p2e12.3)',
c     1       i,nbar, Tc_3SC(i),r_3_SC,3.d0*r_3SC*(cv_u+cv_d+cv_s),cv
c *******************************************************************
c 2SC PHASE:
         else
c Quarks u & d, colors 1 & 2:
          if (T.le.raise*Tc_2SC(i)) then
           tt=min(0.999999999999d0,T/Tc_2SC(i))
           u=u_1s0(tt)
           if      (type_2SC.eq.0.d0) then
            r_2SC=r_1s0(u)
            if (T .gt. Tc_2SC(i)) then
             w1=(raise*Tc_2SC(i)-T)/((raise-1.d0)*Tc_2SC(i))
             w2=1.d0-w1
             r_2SC=w1*r_2SC+w2*1.d0
            end if
           else if (type_2SC.eq.1.d0) then
            r_2SC=tt**2
           else if (type_2SC.eq.2.d0) then
            r_2SC=tt
           else
            pause 'Spec_heat: type_2SC wrong !'
           end if
          end if
c Quarks u color 3:
          if (T.le.raise*Tc_u(i)) then
           tt=min(0.999999999999d0,T/Tc_u(i))
           u=u_1s0(tt)
           if      (type_u.eq.0.d0) then
            r_u=r_1s0(u)
            if (T .gt. Tc_u(i)) then
             w1=(raise*Tc_u(i)-T)/((raise-1.d0)*Tc_u(i))
             w2=1.d0-w1
             r_u=w1*r_u+w2*1.d0
            end if
           else if (type_u.eq.1.d0) then
            r_u=tt**2
           else if (type_u.eq.2.d0) then
            r_u=tt
           else
            pause 'Spec_heat: type_u wrong !'
           end if
          end if
c Quarks d color 3:
          if (T.le.raise*Tc_d(i)) then
           tt=min(0.999999999999d0,T/Tc_d(i))
           u=u_1s0(tt)
           if      (type_d.eq.0.d0) then
            r_d=r_1s0(u)
            if (T .gt. Tc_d(i)) then
             w1=(raise*Tc_d(i)-T)/((raise-1.d0)*Tc_d(i))
             w2=1.d0-w1
             r_d=w1*r_d+w2*1.d0
            end if
           else if (type_d.eq.1.d0) then
            r_d=tt**2
           else if (type_d.eq.2.d0) then
            r_d=tt
           else
            pause 'Spec_heat: type_d wrong !'
           end if
          end if
c Quarks s with 3 colors:
          if (T.le.raise*Tc_s(i)) then
           tt=min(0.999999999999d0,T/Tc_s(i))
           u=u_1s0(tt)
           if      (type_s.eq.0.d0) then
            r_s=r_1s0(u)
            if (T .gt. Tc_s(i)) then
             w1=(raise*Tc_s(i)-T)/((raise-1.d0)*Tc_s(i))
             w2=1.d0-w1
             r_s=w1*r_s+w2*1.d0
            end if
           else if (type_s.eq.1.d0) then
            r_s=tt**2
           else if (type_s.eq.2.d0) then
            r_s=tt
           else
            pause 'Spec_heat: type_s wrong !'
           end if
          end if
c *********************************************************************
          cv = cv_e + 2.d0*r_2SC*(cv_u+cv_d) +
     1         r_u*cv_u + r_d*cv_d + 3.d0*r_s*cv_s
         end if                           
c *********************************************************************
c OLD SUBROUTINE FROM USOV's PAPER:
c        nb=nbar
c        cv_old=2.5d20*(nb/0.16)**(2./3.)*(T/1.d9)
c        cv=cv_old
cc        print *,cv,cv_old
c *********************************************************************
       return
      end
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
      subroutine cvion (T,rho,A,A1,Z,cv)
c      A=A_cell  A1=A_ion
       implicit real*8(a-h,k-z)
       common /control_cvion/icvion
c HERE DANY ccccccccccccccccccc
        if (icvion.eq.1) then
c         ! Original version
         call cvion_1(T,rho,A,Z,cv)
        else if (icvion.eq.2) then
c         ! Modified original version
         call cvion_2(T,rho,A,A1,Z,cv)
        else if ((icvion.eq.3).or.(icvion.eq.4)) then
         call cvion_3(T,rho,A,A1,Z,cv)
        else
         pause 'Subroutine cvion: icvion undefined !'
        end if
       return
      end
c *********************************************************************
c *********************************************************************
c     ORIGINAL VERSION from my PhD thesis:
c     This only inputs a=A_cell and sets A_ion=a1=3*z
      subroutine cvion_1 (t,rho,a,z,cv)

c *** checked on oct. 22 1990 *******

      implicit real*8(a-h,k-z)
      parameter (rhodrip=4.3e11)
     
      dimension cv0(0:14)
      save bcv,ccv,dcv,hcv,cte,cv0
     
      data bcv,ccv,dcv,hcv,cte/.95043,.18956,-.81487,3225,1.417e2/
      data cv0/0.0,2.956,2.829,2.633,2.389,2.118,1.840,1.572,
     1             1.323,1.102,0.909,0.745,0.609,0.496,0.404/

     
      gamma=2.273e5*z**2*(rho/a)**(1./3.)/t
      if (rho .ge. rhodrip) then
       a1=3.*z
      else
       a1=a
      endif
      nionkb=1.38e-16*6.022e23*rho/a
      delta=1./t*z*dsqrt(rho/(a1*a))*6.022e23
      if (gamma .le. .1) then
       cv=1.5*nionkb
       return
      else if (gamma .le. .2) then
       cv1=1.5*nionkb
       cv2=nionkb*(.75*bcv*gamma**.25+1.25*ccv/gamma**.25+dcv+1.5)
       cv=(gamma-.1)/.1*cv2+(.2-gamma)/.1*cv1
       return
      else if (gamma .le. 178.) then
       cv=nionkb*(.75*bcv*gamma**.25+1.25*ccv/gamma**.25+dcv+1.5)
       return
      else if ((gamma.le.210.).and.(delta.ge.1.e19)) then
       cv1=nionkb*(1.5+3.*hcv/gamma**2+1.5)
       cv0(0)=1.5+3.*hcv/gamma**2+1.5
       i1=int(delta*2.e-20)
       cv2=nionkb*(cv0(i1)+(delta*2.e-20-i1)*(cv0(i1+1)-cv0(i1)))
       cv=(gamma-178.)/32.*cv2+(210.-gamma)/32.*cv1
       return
      else if (delta .le. 1e19) then
       cv=nionkb*(1.5+3.*hcv/gamma**2+1.5)
       return
      else if ((delta .gt. 1e19) .and. (delta .lt. 7e20)) then
       cv0(0)=1.5+3.*hcv/gamma**2+1.5
       i1=int(delta*2.e-20)
       cv=nionkb*(cv0(i1)+(delta*2.e-20-i1)*(cv0(i1+1)-cv0(i1)))
       return
      else
       delta1=delta*1.e-20
       cv=nionkb*cte/delta1**3
       return
      end if
     
      return
     
      end

c *********************************************************************
c *********************************************************************
c     MODIFIED ORIGINAL VERSION (from my PhD thesis)
c     This gets A_cell=a and also A_ion=a1 from input
      subroutine cvion_2 (t,rho,a,a1,z,cv)
     
c *** checked on oct. 22 1990 *******

      implicit real*8(a-h,k-z)
      parameter (rhodrip=4.3e11)
     
      dimension cv0(0:14)
      save bcv,ccv,dcv,hcv,cte,cv0
     
      data bcv,ccv,dcv,hcv,cte/.95043,.18956,-.81487,3225,1.417e2/
      data cv0/0.0,2.956,2.829,2.633,2.389,2.118,1.840,1.572,
     1             1.323,1.102,0.909,0.745,0.609,0.496,0.404/

     
      gamma=2.273e5*z**2*(rho/a)**(1./3.)/t
c      if (rho .ge. rhodrip) then
c       a1=3.*z
c      else
c       a1=a
c      endif
      nionkb=1.38e-16*6.022e23*rho/a
      delta=1./t*z*dsqrt(rho/(a1*a))*6.022e23
      if (gamma .le. .1) then
       cv=1.5*nionkb
       return
      else if (gamma .le. .2) then
       cv1=1.5*nionkb
       cv2=nionkb*(.75*bcv*gamma**.25+1.25*ccv/gamma**.25+dcv+1.5)
       cv=(gamma-.1)/.1*cv2+(.2-gamma)/.1*cv1
       return
      else if (gamma .le. 178.) then
       cv=nionkb*(.75*bcv*gamma**.25+1.25*ccv/gamma**.25+dcv+1.5)
       return
      else if ((gamma.le.210.).and.(delta.ge.1.e19)) then
       cv1=nionkb*(1.5+3.*hcv/gamma**2+1.5)
       cv0(0)=1.5+3.*hcv/gamma**2+1.5
       i1=int(delta*2.e-20)
       cv2=nionkb*(cv0(i1)+(delta*2.e-20-i1)*(cv0(i1+1)-cv0(i1)))
       cv=(gamma-178.)/32.*cv2+(210.-gamma)/32.*cv1
       return
      else if (delta .le. 1e19) then
       cv=nionkb*(1.5+3.*hcv/gamma**2+1.5)
       return
      else if ((delta .gt. 1e19) .and. (delta .lt. 7e20)) then
       cv0(0)=1.5+3.*hcv/gamma**2+1.5
       i1=int(delta*2.e-20)
       cv=nionkb*(cv0(i1)+(delta*2.e-20-i1)*(cv0(i1+1)-cv0(i1)))
       return
      else
       delta1=delta*1.e-20
       cv=nionkb*cte/delta1**3
       return
      end if
     
      return
     
      end

c *********************************************************************
c *********************************************************************
c     NEW VERSION
      subroutine cvion_3 (t,rho,a,a1,z,cv)
c     a=A_cell & a1=A_ion
c *** Checked on March 16, 2010
       implicit real*8(a-h,k-z)
       INCLUDE 'rho_limits.inc.f'
       INCLUDE 'gamma_limits.inc.f'
       common /control_cvion/icvion
c       INCLUDE 'plasma.inc.f'
       parameter (pi=4.d0*atan(1.d0))
       parameter (Na=6.022d23,kb=1.380d-16,Mu=1.d0/Na)
       parameter (me=9.109d-28,e=4.803206d-10)
       parameter (hb=1.054572d-27,c=2.99792458d10)
       parameter (hbc=197.327d0,MeV=1.602177d-6)
       data bcv,ccv,dcv/.95043,.18956,-.81487/
       save bcv,ccv,dcv
        ni=Na*rho/A
        a_WS=(3.d0/(4.d0*pi*ni))**(1./3.)
        Gamma=(Z*e)**2/(a_WS*kb*T)
        if (Gamma.le.gammaliq) then
         call cv_liquid (T,rho,A,A1,Z,cv1)
c	write(90,*)'liquid', rho
c         From Slattery et al.: actually gives the same result as cv_liquid
c         at gamma >1 !
c         cv1=.75*bcv*gamma**.25+1.25*ccv/gamma**.25+dcv+1.5
         cv=cv1
        else if (Gamma.gt.gammaliq) then
         call cv_harmonic_crystal (T,rho,A,A1,Z,cv2h)
         call cv_anharmonic_crystal (T,rho,A,A1,Z,cv2ah)
         cv=cv2h+cv2ah
c	write(90,*)'solid', rho
c HERE DANY cccccccccccccccccccccccccccccccccccc
c         cv=cv2h
c         if (idany_cvion.ne.13) then
c          print *,'============================================='
c          print *,'   WARNING: no anharmonic effecs in Cv_ion'
c          print *,'============================================='
c          idany_cvion=13
c          read(5,*)
c          end if
cccccccccccccccccccccccccccccccccccccccccccccccc
         if (Gamma.le.gammacryst) then
          call cv_liquid (T,rho,A,A1,Z,cv1)
c	  write(90,*)'mixture', rho
c         From Slattery et al.: actually gives the same result as cv_liquid
c         at gamma >1 !
c          cv1=.75*bcv*gamma**.25+1.25*ccv/gamma**.25+dcv+1.5
          cv2=cv
          dg=gammacryst-gammaliq
          w1=(gammacryst-Gamma)/dg
          w2=1.d0-w1
          cv=w1*cv1+w2*cv2
         end if
        end if
        cv=cv*ni*kb
       return
      end
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
c     Baiko, Potekhin & Yakovlev, 2001PhRvE..64e7402B
      subroutine cv_harmonic_crystal (T,rho,A,A1,Z,cv)
c      A=A_cell A1=A_ion
c      cv is returned per ion in units of kb 
c      (i.e., must be multiplied by nion*kb !)
c      CHECKED on March 18, 2010 against 
c      Fig 2 of Potekhin, A. Y.; Chabrier, G. 2010CoPP...50...82P
       implicit real*8(a-h,k-z)
       parameter (pi=4.d0*atan(1.d0))
       parameter (Na=6.022d23,kb=1.380d-16,Mu=1.d0/Na)
       parameter (me=9.109d-28,e=4.803206d-10)
       parameter (hb=1.054572d-27,c=2.99792458d10)
       parameter (hbc=197.327d0,MeV=1.602177d-6)
       dimension aa(0:8),bb(0:8),al(0:8)
       data aa/ 1.d0 , 1.839d-1 , 5.93586d-1 , 5.4814d-3 , 5.01813d-4 ,
     x         0.d0 , 3.9247d-7 , 0.d0 , 5.8356d-11 /
       data bb/ 2.6166d2 , 0.d0 , 7.07997d0 , 0.d0 , 4.09484d-2 , 
     x         3.97355d-4 , 5.11148d-5  , 2.19749d-6 , 0.d0 /
       data al/0.d0 , 9.32446d-1 , 3.34547d-1 , 2.65764d-1 , 0.d0 , 
     x         0.d0 , 4.757014d-3 , 0.d0 , 4.7770935d-3 /
       save aa,bb,al
c *********************************************************************
       FA(x)=aa(0)+aa(1)*x+aa(2)*x**2+aa(3)*x**3+aa(4)*x**4+
     x       aa(5)*x**5+aa(6)*x**6+aa(7)*x**7+aa(8)*x**8
       FA1(x)=aa(1)+2.d0*aa(2)*x+3.d0*aa(3)*x**2+4.d0*aa(4)*x**3+
     x  5.d0*aa(5)*x**4+6.d0*aa(6)*x**5+7.d0*aa(7)*x**6+8.d0*aa(8)*x**7
       FA2(x)=2.d0*aa(2)+6.d0*aa(3)*x+12.d0*aa(4)*x**2+20.d0*aa(5)*x**3+
     x     30.d0*aa(6)*x**4+42.d0*aa(7)*x**5+56.d0*aa(8)*x**6
       FB(x)=bb(0)+bb(1)*x+bb(2)*x**2+bb(3)*x**3+bb(4)*x**4+
     x       bb(5)*x**5+bb(6)*x**6+bb(7)*x**7+
     x       al(6)*aa(6)*x**9+al(8)*aa(8)*x**11
       FB1(x)=bb(1)+2.d0*bb(2)*x+3.d0*bb(3)*x**2+4.d0*bb(4)*x**3+
     x     5.d0*bb(5)*x**4+6.d0*bb(6)*x**5+7.d0*bb(7)*x**6+
     x     9.d0*al(6)*aa(6)*x**8+11.d0*al(8)*aa(8)*x**10
       FB2(x)=2.d0*bb(2)+6.d0*bb(3)*x+12.d0*bb(4)*x**2+20.d0*bb(5)*x**3+
     x     30.d0*bb(6)*x**4+42.d0*bb(7)*x**5+
     x     72.d0*al(6)*aa(6)*x**7+110.d0*al(8)*aa(8)*x**9
c *********************************************************************
        ni=Na*rho/A
        om_p=dsqrt(4.d0*pi*ni*Z**2*e**2/(A1*mu))
        T_p=hb*om_p/kb
c ****
c      This is to use as input T_p
c      entry cv_harmonic_crystal_entry(T_p0,T,cv)
c        T_p=T_p0
c ****
        the=T_p/T
        N1=al(1)**2*the**2
        D1=(dexp(al(1)*the/2.d0)-dexp(-al(1)*the/2.d0))**2
        N2=al(2)**2*the**2
        D2=(dexp(al(2)*the/2.d0)-dexp(-al(2)*the/2.d0))**2
        N3=al(3)**2*the**2
        D3=(dexp(al(3)*the/2.d0)-dexp(-al(3)*the/2.d0))**2
        N=FA2(the)*FB(the)**2 - 2.d0*FA1(the)*FB1(the)*FB(the) +
     x    2.d0*FA(the)*FB1(the)**2 - FA(the)*FB(the)*FB2(the)
        D=FB(the)**3
        cv=N1/D1 + N2/D2 + N3/D3 + the**2*N/D
c        cv=ni*kb*cv
       return
      end
c *********************************************************************
c *********************************************************************
c     Potekhin, A. Y.; Chabrier, G. 2010CoPP...50...82P
c     From their Eq. 8
      subroutine cv_anharmonic_crystal (T,rho,A,A1,Z,cv)
c      A=A_cell A1=A_ion
c      cv is returned per ion in units of kb 
c      (i.e., must be multiplied by nion*kb !)
c      CHECKED on March 18, 2010 against 
c      Fig 2 of Potekhin, A. Y.; Chabrier, G. 2010CoPP...50...82P
       implicit real*8(a-h,k-z)
       parameter (pi=4.d0*atan(1.d0))
       parameter (Na=6.022d23,kb=1.380d-16,Mu=1.d0/Na)
       parameter (me=9.109d-28,e=4.803206d-10)
       parameter (hb=1.054572d-27,c=2.99792458d10)
       parameter (hbc=197.327d0,MeV=1.602177d-6)
       data aa1,aa2,aa3,c1/10.9d0 , 247.d0 , 1.765d5 , 0.0112d0 /
        ni=Na*rho/A
        aWS=(3.d0/(4.d0*pi*ni))**(1./3.)
        om_p=dsqrt(4.d0*pi*ni*Z**2*e**2/(A1*mu))
        T_p=hb*om_p/kb
        GT=(Z*e)**2/aWS/kb
c ****
c      This is to us as input Gamma*T and T_p
c      entry cv_anharmonic_crystal_entry(GT0,T_p0,T,cv)
c        GT=GT0
c        T_p=T_p0
c ****
        cv=T *( 
     x      (4.d0*aa3/GT**3)*T**2 +
     x      (3.d0*aa2/GT**2)*T + 
     x      (10.d0*c1*T_p**2*aa3+6.d0*aa1*GT**2)/(3.d0*GT**3) +
     x      (3.d0*c1*T_p**2*aa2/GT**2)/T +
     x      (4.d0*c1**2*T_p**4*aa3 + 6.d0*c1*T_p**2*aa1*GT**2)/
     x                                      (3.d0*GT**3)/T**2 +
     x      (2.d0*c1**2*T_p**4*aa2/GT**2)/T**3 + 
     x      (4.d0*c1**2*T_p**4*aa1/GT)/T**4
     x      ) * dexp(-c1*T_p**2/T**2)
       return
      end
c *********************************************************************
c *********************************************************************
c     Potekhin, A. Y.; Chabrier, G. 2000PhRvE..62.8554P
c     From their Eq. 17
      subroutine cv_liquid (T,rho,A,A1,Z,cv)
c      A=A_cell A1=A_ion
c      cv is returned per ion in units of kb 
c      (i.e., must be multiplied by nion*kb !)
c      CHECKED on March 21, 2010
       implicit real*8(a-h,k-z)
       parameter (pi=4.d0*atan(1.d0))
       parameter (Na=6.022d23,kb=1.380d-16,Mu=1.d0/Na)
       parameter (me=9.109d-28,e=4.803206d-10)
       parameter (hb=1.054572d-27,c=2.99792458d10)
       parameter (hbc=197.327d0,MeV=1.602177d-6)
       data AA1,AA2,AA3/-0.9070,0.62954,0.27710/
       data BB1,BB2,BB3,BB4/4.56d-3,211.6d0,1.0d-4,4.62d-3/
        ni=Na*rho/A
        aWS=(3.d0/(4.d0*pi*ni))**(1./3.)
        G=(Z*e)**2/aWS/kb/T
        om_p=dsqrt(4.d0*pi*ni*Z**2*e**2/(A1*mu))
        T_p=hb*om_p/kb
c ****
cc      This is to use as input Gamma and T_p
cc       (Notice that if T_p0=0 you only get the classical term !)
c      entry cv_liquid_entry(G0,T_p0,cv)
c        G=G0
c        T_p=T_p0
c ****  Classical term (Chabrier & Potekhin, 1998):
        cv=1.5d0 + 
     x     G**1.5/2.d0 *
     x     ( AA3*(G-1.d0)/(G+1.d0)**2 - AA1*AA2/(G+AA2)**1.5 ) +
     x     G**2 *
     x     ( BB3*(G**2-BB4)/(G**2+BB4)**2 - BB1*BB2/(G+BB2)**2  )
c ***   First quantum correction (hb^2):
c HERE DANY cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        cv_q=-T_p**2/T**2/12.d0
c         cv_q=0.d0
c        if (idany_cvion.ne.13) then
c         print *,'====================================================='
c         print *,'           CV_liq: no quantum correction   '
c         print *,'              (press <ENTER> to continue)'
c         print *,'====================================================='
c         idany_cvion=13
c         read(5,*)
c        end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ***
c        write(94,*)cv,cv_q,cv+cv_q  
        cv=cv+cv_q
       return
      end
c *********************************************************************
c *********************************************************************
      subroutine cvnrot (t,rho,a,z,cv)
     
c ****** to be checked

      implicit real*8(a-h,k-z)
      parameter (rhodrip=4.3e11)
      parameter (alpha=1.)
     
      if (rho .ge. rhodrip) then
       a1=3.*z
      else
       a1=a
      endif
      nionkb=1.38e-16*6.022e23*rho/a
      thetar=4.4e11/a1**(5./3.)/alpha     
      x=thetar/t

      if (x.ge.1.e-1)then
       jmax=16
       if(x.ge.1.e0)jmax=6
       fz=0.
       fz1=0.
       fz2=0.
       do j=0,jmax,2
        dz=(2.*j+1.)*dexp(-x*j*(j+1.))
        dz1=dz*j*(j+1.)
        dz2=dz1*j*(j+1.)
        fz=fz+dz
        fz1=fz1-dz1
        fz2=fz2+dz2
       end do
      else
       fz=1./(2.*x)+0.17
       fz1=-1./(2.*x**2)
       fz2=1./x**3
      end if

      cv=nionkb*x**2/fz*(fz2-fz1**2/fz)

      return
     
      end
     
c ****************************************************************************
c ****************************************************************************
     
      subroutine cvnvib (t,i,rho,a,z,cv)
     
c ****** to be checked

      implicit real*8(a-h,k-z)
      INCLUDE 'size.inc.f'

      dimension theta_vib(4,0:isize),c(4)
      common/cvn_vib/theta_vib
      save c
      data c/5.,15.,35.,7./
c ********************************
      fexp(x)=dexp(max(x,-7.d2))
c ********************************
      nionkb=1.38e-16*6.022e23*rho/a
      s0=1.
      s1=0.
      s2=0.
      do j=1,4
       ds=c(j)/fexp(theta_vib(j,i)/t)
       s0=s0+ds
       s1=s1+(theta_vib(j,i)/t)*ds
       s2=s2+(theta_vib(j,i)/t)**2*ds
      end do

      cv=nionkb*(s1/s0+(s0*s2-s1**2)/s0**2)

      return
     
      end
     
c ************************************************************************
c ************************************************************************

      subroutine cvelec (t,rho,a,z,ikeep,cv)

c *******  checked on sept 4  1991 **************
     
      implicit real*8(a-h,k-z)
      INCLUDE 'size.inc.f'
      parameter (pi = 3.14159625)
     
      dimension cd(0:3,0:3),cp(0:3,0:3),cu(0:3,0:3)
      dimension fkeep(0:isize)
      save cd,cp,cu,fkeep
     
      data ((cd(i,j),i=0,3),j=0,3)
     1         / 2.315472, 7.128660, 7.504998, 2.665350,
     2           7.837752,23.507934,23.311317, 7.987465,
     3           9.215560,26.834068,25.082745, 8.020509,
     4           3.693280,10.333176, 9.168960, 2.668248/
     
      data ((cp(i,j),i=0,3),j=0,3)
     1         / 2.315472, 6.748104, 6.564912, 2.132280,
     2           7.837752,21.439740,19.080088, 5.478100,
     3           9.215560,23.551504,19.015888, 4.679944,
     4           3.693280, 8.859868, 6.500712, 1.334124/
     
      data ((cu(i,j),i=0,3),j=0,3)
     1         / 3.473208,10.122156, 9.847368, 3.198420,
     2          16.121172,43.477194,37.852852,10.496830,
     3          23.971040,60.392810,47.782844,11.361074,
     4          11.079840,26.579604,19.502136, 4.002372/
     
      hb=1.054588e-27
      kb=1.380662e-16
      c=2.997924e10
      na=6.022045e23
      me=9.109e-28
     
      ne=na*rho*z/a
      pf=hb*(3.*pi**2*ne)**(1./3.)
      ef=dsqrt((me*c**2)**2+(pf*c)**2)-me*c**2
      tf=ef/kb
      xe=pf/(me*c)
      ae=xe**2/dsqrt(1.+xe**2)
      cvt=ne*kb**2*pi**2/(me*c**2)/ae
     
      t0=tf/50.
     
      if(t .le. .5*t0)then
       cv=cvt*t
       return
      end if
     
      nehat=ne/1.7595e30
     
c ***** calculate f at .9999*t
     
      t1=0.9999 * t/5.93e9

      if (fkeep(ikeep).eq.0.) then
       f1=1.
      else
       f1=fkeep(ikeep)
      end if

1111  f=abs(f1)
      g=t1*dsqrt(1.+f)

      sum1=0.
      do 100 j1=0,3
       do 101 j2=0,3
        sum1=sum1+cd(j1,j2)*f**j1*g**j2
101    continue
100   continue
      sum2=0.
      do 102 j1=1,3
       do 103 j2=0,3
        sum2=sum2+j1*cd(j1,j2)*f**(j1-1)*g**j2
103    continue
102   continue
      sum3=0.
      do 104 j1=0,3
       do 105 j2=1,3
        sum3=sum3+j2*cd(j1,j2)*f**j1*g**j2
105    continue
104   continue
     
      pf=1.+f
      pg=1.+g
      if(f.lt.1.)then
       nef=f*g**1.5*pf**(-4)*pg**(-1.5)*sum1
       nef1=(g**1.5*pf**(-4)/pg**1.5-
     1       3.25*f*g**1.5*pf**(-5)*pg**(-1.5)-
     2       .75*f*g**2.5*pf**(-5)*pg**(-2.5)) * sum1 +
     3       f*g**1.5*pf**(-4)*pg**(-1.5)*(sum2+0.5/pf*sum3)
      else
       nef=(f**.25/pf)**4*(g/pg)**1.5*sum1
       nef1=(f**.25/pf)**4*(g/pg)**1.5*
     1      ( (1./f-3.25/pf-0.75*g/(pf*pg))*sum1+sum2+0.5/pf*sum3 )
      end if
     
      f1=abs(f+(nehat-nef)/nef1)
     
      if ((abs(nef-nehat)/abs(nehat).gt.1.e-10).or.
     1    (abs(f1-f)/abs(f1).gt.1.e-12)) goto 1111

      fkeep(ikeep)=f1

c *****  calculate f at 1.0001*t
     
      t2=1.0001 * t/5.93e9
      f2=f1
     
1112  f=abs(f2)
      g=t2*dsqrt(1.+f)
     
      sum1=0.
      do 106 j1=0,3
       do 107 j2=0,3
        sum1=sum1+cd(j1,j2)*f**j1*g**j2
107    continue
106   continue
      sum2=0.
      do 108 j1=1,3
       do 109 j2=0,3
        sum2=sum2+j1*cd(j1,j2)*f**(j1-1)*g**j2
109    continue
108   continue
      sum3=0.
      do 110 j1=0,3
       do 111 j2=1,3
        sum3=sum3+j2*cd(j1,j2)*f**j1*g**j2
111    continue
110   continue
     
      pf=1.+f
      pg=1.+g
      if(f.lt.1)then
       nef=f*g**1.5*pf**(-4)*pg**(-1.5)*sum1
       nef1=(g**1.5*pf**(-4)/pg**1.5-
     1       3.25*f*g**1.5*pf**(-5)*pg**(-1.5)-
     2       .75*f*g**2.5*pf**(-5)*pg**(-2.5)) * sum1 +
     3       f*g**1.5*pf**(-4)*pg**(-1.5)*(sum2+0.5/pf*sum3)
      else
       nef=(f**.25/pf)**4*(g/pg)**1.5*sum1
       nef1=(f**.25/pf)**4*(g/pg)**1.5*
     1      ( (1./f-3.25/pf-0.75*g/pf/pg)*sum1+sum2+0.5/pf*sum3 )
      end if
     
     
      f2=abs(f+(nehat-nef)/nef1)
     
      if ((abs(nef-nehat)/abs(nehat).gt.1.e-10).or.
     1    (abs(f2-f)/abs(f2).gt.1.e-12)) goto 1112
     
c *****  calculate cv
     
      sumu2=0.
      g2=t2*dsqrt(1.+f2)
      do 120 j1=0,3
       do 121 j2=0,3
        sumu2=sumu2+cu(j1,j2)*f2**j1*g2**j2
121    continue
120   continue
      sumu1=0.
      g1=t1*dsqrt(1.+f1)
      do 122 j1=0,3
       do 123 j2=0,3
        sumu1=sumu1+cu(j1,j2)*f1**j1*g1**j2
123    continue
122   continue
      if(f1.lt.1)then
       u1=1.44e24*f1*g1**2.5/(1.+f1)**4/(1.+g1)**1.5*sumu1
       u2=1.44e24*f2*g2**2.5/(1.+f2)**4/(1.+g2)**1.5*sumu2
      else
       u1=1.44e24*(f1**.25/(1.+f1))**4*g1*(g1/(1.+g1))**1.5*sumu1
       u2=1.44e24*(f2**.25/(1.+f2))**4*g2*(g2/(1.+g2))**1.5*sumu2
      end if
     
      cv=(u2-u1)/(.0002*t)
     
      if (t .lt. 1.5*t0) then
       w1=(t-.5*t0)/t0
       w2=1.-w1
       cv=w1*cv+w2*cvt*t
      end if
     
      return
     
      end
     
      subroutine MELANGE8(NMIX,AY,AZion,ACMI,RHO,TEMP,KRAD,
     *   DENS,Zmean,CMImean,Z2mean,GAMImean,CHI,TPT,LIQSOL,
     *   PnkT,UNkT,SNk,CV,CHIR,CHIT)
*                                                       Version 11.09.08
*                                                          corr.26.12.09
* EOS of fully ionized electron-ion plasma mixture.     
* Previous version: MELANGE7 v.11.06.07 (in file eos2007.f).
* Limitations:
* (a) inapplicable in the regimes of
*      (1) bound-state formation,
*      (2) quantum liquid,
*      (3) presence of positrons;
* (b) for the case of a composition gradually depending on RHO or TEMP,
*  second-order functions (CV,CHIR,CHIT in output) should not be trusted
* Choice of the liquid or solid regime - criterion GAMI [because the
*     choice based on comparison of total (non-OCP) free energies can be
*     sometimes dangerous because of the fit uncertainties ("Local field
*     correction" in solid and quantum effects in liquid are unknown)].
* Input: NMIX - number of different elements;
*        AY - their partial number densities,
*        AZion and ACMI - their charge and mass numbers,
*        RHO - total mass density [g/cc]
*        TEMP - temperature [in a.u.=2Ryd=3.1577e5 K].
*        KRAD=1/0 to include/exclude radiation (photon) quantities
* NB: instead of RHO, a true input is CHI, defined below
*     Hence, disagreement between RHO and DENS is the fit error (<0.4%)
* Output: AY - rescaled so that to sum up to 1
*         DENS - electron number density [in a.u.=6.7483346e24 cm^{-3}]
*         Zmean=<Z>, CMImean=<A> - mean ion charge and mass numbers,
*         Z2mean=<Z^2> - mean-square ion charge number
*         GAMImean - effective ion-ion Coulomb coupling constant
*         CHI = mu_e/kT, where mu_e is the electron chem.potential
*         TPT - effective ionic quantum parameter (T_p/T)
*         LIQSOL=0/1 for liquid/solid
*         SNk - dimensionless entropy per 1 ion
*         UNkT - internal energy per kT per ion
*         PnkT - pressure / n_i kT, where n_i is the ion number density
*         CV - heat capacity per ion, div. by Boltzmann const.
*         CHIR - inverse compressibility -(d ln P / d ln V)_T ("\chi_r")
*         CHIT = (d ln P / d ln T)_V ("\chi_T")
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter(TINY=1.d-7,MAXY=10)
      dimension AY(*),AZion(*),ACMI(*)
      parameter (PI=3.141592653d0,AUM=1822.888d0, ! a.m.u./m_e
     *   GAMIMELT=175., ! OCP value of Gamma_i for melting
     *   RSIMELT=140., ! ion density parameter of quantum melting
     *   RAD=2.554d-7) ! Radiation constant (=4\sigma/c) (in a.u.)
      if (NMIX.gt.MAXY) stop'MELANGE8: NMIX > MAXY'
      if (KRAD.ne.0.and.KRAD.ne.1) stop'MELANGE8: wrong KRAD'
      Y=0.
      do IX=1,NMIX
         Y=Y+AY(IX)
      enddo
      if (dabs(Y-1.d0).gt.TINY) then
        do IX=1,NMIX
           AY(IX)=AY(IX)/Y
        enddo
         print*,'MELANGE8: partial densities (and derivatives)',
     *    ' are rescaled by factor',1./Y
      endif
      Zmean=0.
      Z2mean=0.
      Z52=0.
      Z53=0.
      Z73=0.
      Z321=0. ! corr.26.12.09
      CMImean=0.
      do IX=1,NMIX
         Zmean=Zmean+AY(IX)*AZion(IX)
         Z2mean=Z2mean+AY(IX)*AZion(IX)**2
         Z13=AZion(IX)**(1./3.)
         Z53=Z53+AY(IX)*Z13**5 
         Z73=Z73+AY(IX)*Z13**7
         Z52=Z52+AY(IX)*dsqrt(AZion(IX))**5
         Z321=Z321+AY(IX)*AZion(IX)*dsqrt(AZion(IX)+1.d0)**3 ! 26.12.09
         CMImean=CMImean+AY(IX)*ACMI(IX)
      enddo
* (0) Photons (=0 if KRAD=0):
      UINTRAD=RAD*TEMP**4*KRAD
      PRESSRAD=UINTRAD/3.
      CVRAD=4.*UINTRAD/TEMP
* (1) ideal electron gas (including relativity and degeneracy)  -----  *
      DENS=RHO/11.20587*Zmean/CMImean ! number density of electrons [au]
      call CHEMFIT(DENS,TEMP,CHI)
* NB: CHI can be used as true input instead of RHO or DENS
      call ELECT9(TEMP,CHI,DENS,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE,DDN)
* NB: at this point DENS is redefined (the difference can be ~0.1%)
      DTE=DENS*TEMP
      PRESSE=PEid*DTE ! P_e [a.u.]
      UINTE=UEid*DTE ! U_e / V [a.u.]
* (2) non-ideal Coulomb EIP  ----------------------------------------  *
      RS=(.75/PI/DENS)**.333333 ! r_s - electron density parameter
      RSI=RS*CMImean*Z73*AUM ! R_S - ion density parameter
      GAME=1./RS/TEMP ! electron Coulomb parameter Gamma_e
      GAMImean=Z53*GAME   ! effective Gamma_i - ion Coulomb parameter
      if (GAMImean.lt.GAMIMELT.or.RSI.lt.RSIMELT) then
         LIQSOL=0 ! liquid regime
      else
         LIQSOL=1 ! solid regime
      endif
      TPT2=0.
      do IX=1,NMIX
         TPT2=TPT2+AY(IX)*GAME**2/RS*3./(AUM*ACMI(IX))*AZion(IX)
* In the case of a mixture, this estimate of the quantum ion parameter
* is as crude as the linear mixing model for Fharm employed below.
      enddo
      if (LIQSOL.eq.0.and.TPT2.gt.18.) then
         print*,'MELANGE8: strong quantum effects in liquid!'
C         read(*,'(A)')
      endif
      TPT=dsqrt(TPT2) ! effective T_p/T - ion quantum parameter
* Calculate partial thermodynamic quantities and combine them together:
      UINT=UINTE
      PRESS=PRESSE
      CVtot=CVE*DENS
      Stot=SEid*DENS
      PDLT=PRESSE*CHITE ! d P_e[a.u.] / d ln T
      PDLR=PRESSE*CHIRE ! d P_e[a.u.] / d ln\rho
      DENSI=DENS/Zmean ! number density of all ions
      PRESSI=DENSI*TEMP ! ideal-ions total pressure (normalization)
* Add Coulomb+xc nonideal contributions, and ideal free energy:
      do IX=1,NMIX
        if (AY(IX).lt.TINY) goto 10 ! skip this species
         Zion=AZion(IX)
         CMI=ACMI(IX)
         GAMI=Zion**1.666667/RS/TEMP ! Gamma_i for given ion species
         DNI=DENSI*AY(IX) ! number density of ions of given type
         PRI=DNI*TEMP ! = ideal-ions partial pressure (normalization)
         call EOSFI8(LIQSOL,CMI,Zion,RS,GAMI,
     *     FC1,UC1,PC1,SC1,CV1,PDT1,PDR1,
     *     FC2,UC2,PC2,SC2,CV2,PDT2,PDR2)
* First-order TD functions:
         UINT=UINT+UC2*PRI ! internal energy density (e+i+Coul.)
         Stot=Stot+DNI*(SC2-dlog(AY(IX))) !entropy per unit volume[a.u.]
         PRESS=PRESS+PC2*PRI ! pressure (e+i+Coul.) [a.u.]
* Second-order functions (they take into account compositional changes):
         CVtot=CVtot+DENSI*CV2*AY(IX) ! C_V (e+i+Coul.)/ V
         PDLT=PDLT+PRI*PDT2 ! d P / d ln T
         PDLR=PDLR+PRI*PDR2 ! d P / d ln\rho
   10   continue
      enddo ! next IX
* Corrections to the linear mixing rule:
      call CORMIX(RS,GAME,Zmean,Z2mean,Z52,Z53,Z321,
     *  FMIX,UMIX,PMIX,CVMIX,PDTMIX,PDRMIX)
      UINT=UINT+UMIX*PRESSI
      Stot=Stot+DENSI*(UMIX-FMIX)
      PRESS=PRESS+PMIX*PRESSI
      CVtot=CVtot+DENSI*CVMIX
      PDLT=PDLT+PRESSI*PDTMIX
      PDLR=PDLR+PRESSI*PDRMIX
* First-order:
      PRESS=PRESS+PRESSRAD ! total P [a.u.]
      CVtot=CVtot+CVRAD
      Stot=Stot+CVRAD/3.
      PnkT=PRESS/PRESSI ! beta P / n_i
      UNkT=(UINT+UINTRAD)/PRESSI ! beta U / N_i
      SNk=Stot/DENSI ! S / N_i
* Second-order:
      CV=CVtot/DENSI ! C_V per ion
      CHIR=PDLR/PRESS ! d ln P / d ln\rho
      CHIT=(PDLT+4.*PRESSRAD)/PRESS ! d ln P / d ln T
      return
      end

      subroutine EOSFI8(LIQSOL,CMI,Zion,RS,GAMI,
     *  FC1,UC1,PC1,SC1,CV1,PDT1,PDR1,
     *  FC2,UC2,PC2,SC2,CV2,PDT2,PDR2)
*                                                       Version 16.09.08
* Non-ideal parts of thermodynamic functions in the fully ionized plasma
* Stems from EOSFI5 and EOSFI05 v.04.10.05
* Input: LIQSOL=0/1(liquid/solid), 
*        Zion,CMI - ion charge and mass numbers,
*        RS=r_s (electronic density parameter),
*        GAMI=Gamma_i (ion coupling),
* Output: FC1 and UC1 - non-ideal "ii+ie+ee" contribution to the 
*         free and internal energies (per ion per kT),
*         PC1 - analogous contribution to pressure divided by (n_i kT),
*         CV1 - "ii+ie+ee" heat capacity per ion [units of k]
*         PDT1=(1/N_i kT)*(d P_C/d ln T)_V
*         PDR1=(1/N_i kT)*(d P_C/d ln\rho)_T
* FC2,UC2,PC2,SC2,CV2 - analogous to FC1,UC1,PC1,SC1,CV1, but including
* the part corresponding to the ideal ion gas. This is useful for 
* preventing accuracy loss in some cases (e.g., when SC2 << SC1).
* FC2 does not take into account the entropy of mixing S_{mix}: in a
* mixture, S_{mix}/(N_i k) has to be added externally (see MELANGE8).
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter(TINY=1.d-20) ! prevents division by zero
      parameter (AUM=1822.888d0) ! a.m.u/m_e
      if (LIQSOL.ne.1.and.LIQSOL.ne.0) stop'EOSFI8: invalid LIQSOL'
      if (CMI.le..1) stop'EOSFI8: too small CMI'
      if (Zion.le..1) stop'EOSFI8: too small Zion'
      if (RS.le..0) stop'EOSFI8: invalid RS'
      if (GAMI.le..0) stop'EOSFI8: invalid GAMI'
      GAME=GAMI/Zion**1.666667
      call EXCOR7(RS,GAME,FXC,UXC,PXC,CVXC,SXC,PDTXC,PDRXC) ! "ee"("xc")
* Calculate "ii" part:
      COTPT=dsqrt(3./AUM/CMI)/Zion**(7./6.) ! auxiliary coefficient
      TPT=GAMI/dsqrt(RS)*COTPT              ! T_p/T
      FidION=1.5*dlog(TPT**2/GAMI)-1.323515
* 1.3235=1+0.5*ln(6/pi); FidION = F_{id.ion gas}/(N_i kT), but without
* the term x_i ln x_i = -S_{mix}/(N_i k).
      if (LIQSOL.eq.0) then                 ! liquid
         call FITION9(GAMI,
     *     FION,UION,PION,CVii,PDTii,PDRii)
         FItot=FION+FidION
         UItot=UION+1.5
         PItot=PION+1.d0
         CVItot=CVii+1.5d0
         SCItot=UItot-FItot
         PDTi=PDTii+1.d0
         PDRi=PDRii+1.d0
         FWK=TPT**2/24. ! Wigner-Kirkwood (quantum diffr.) term
         UWK=2.*FWK
         CVWK=-UWK
         PWK=FWK
         PDRWK=2.*PWK
         PDTWK=-PWK
      else                                  ! solid
         call FHARM8(GAMI,TPT,
     *     Fharm,Uharm,Pharm,CVharm,Sharm,PDTharm,PDRharm) ! harm."ii"
         call ANHARM8(GAMI,TPT,Fah,Uah,Pah,CVah,PDTah,PDRah) ! anharm.
         FItot=Fharm+Fah
         FION=FItot-FidION
         UItot=Uharm+Uah
         UION=UItot-1.5d0 ! minus 1.5=ideal-gas, in order to get "ii"
         PItot=Pharm+Pah
         PION=PItot-1.d0 ! minus 1=ideal-gas
         PDTi=PDTharm+PDTah
         PDRi=PDRharm+PDRah
         PDTii=PDTi-1.d0 ! minus 1=ideal-gas
         PDRii=PDRi-1.d0 ! minus 1=ideal-gas
         CVItot=CVharm+CVah
         SCItot=Sharm+Uah-Fah
         CVii=CVItot-1.5d0 ! minus 1.5=ideal-gas
         FWK=0.
         UWK=0.
         PWK=0.
         CVWK=0.
         PDRWK=0.
         PDTWK=0.
      endif
* Calculate "ie" part:
      if (LIQSOL.eq.1) then
         call FSCRsol8(RS,GAMI,Zion,TPT,
     *     FSCR,USCR,PSCR,S_SCR,CVSCR,PDTSCR,PDRSCR)
      else
         call FSCRliq8(RS,GAME,Zion,
     *     FSCR,USCR,PSCR,CVSCR,PDTSCR,PDRSCR)
         S_SCR=USCR-FSCR
      endif
* Total excess quantities ("ii"+"ie"+"ee", per ion):
      FC0=FSCR+Zion*FXC+FWK
      UC0=USCR+Zion*UXC+UWK
      PC0=PSCR+Zion*PXC+PWK
      SC0=S_SCR+Zion*SXC+FWK
      CV0=CVSCR+Zion*CVXC+CVWK
      PDT0=PDTSCR+Zion*PDTXC+PDTWK
      PDR0=PDRSCR+Zion*PDRXC+PDRWK
      FC1=FION+FC0
      UC1=UION+UC0
      PC1=PION+PC0
      SC1=(UION-FION)+SC0
      CV1=CVii+CV0
      PDT1=PDTii+PDT0
      PDR1=PDRii+PDR0
* Total excess + ideal-ion quantities
      FC2=FItot+FC0
      UC2=UItot+UC0
      PC2=PItot+PC0
      SC2=SCItot+SC0
      CV2=CVItot+CV0
      PDT2=PDTi+PDT0
      PDR2=PDRi+PDR0
      return
      end

* ==================  ELECTRON-ION COULOMB LIQUID  =================== *
      subroutine FITION9(GAMI,
     *     FION,UION,PION,CVii,PDTii,PDRii)
*                                                       Version 11.09.08
* Non-ideal contributions to thermodynamic functions of classical OCP,
*       corrected at small density for a mixture.
*   Stems from FITION00 v.24.05.00.
* Input: GAMI - ion coupling parameter
* Output: FION - ii free energy / N_i kT
*         UION - ii internal energy / N_i kT
*         PION - ii pressure / n_i kT
*         CVii - ii heat capacity / N_i k
*         PDTii = PION + d(PION)/d ln T = (1/N_i kT)*(d P_{ii}/d ln T)
*         PDRii = PION + d(PION)/d ln\rho
*   Parameters adjusted to Caillol (1999).
      implicit double precision (A-H),double precision (O-Z)
      save
      parameter (A1=-.907347d0,A2=.62849d0,C1=.004500d0,G1=170.0,
     *  C2=-8.4d-5,G2=.0037,SQ32=.8660254038d0) ! SQ32=sqrt(3)/2
      A3=-SQ32-A1/dsqrt(A2)
      F0=A1*(dsqrt(GAMI*(A2+GAMI))-
     -     A2*dlog(dsqrt(GAMI/A2)+dsqrt(1.+GAMI/A2)))+
     +     2.*A3*(dsqrt(GAMI)-datan(dsqrt(GAMI)))
      U0=dsqrt(GAMI)**3*(A1/dsqrt(A2+GAMI)+A3/(1.d0+GAMI))
*   This is the zeroth approximation. Correction:
      UION=U0+C1*GAMI**2/(G1+GAMI)+C2*GAMI**2/(G2+GAMI**2)
      FION=F0+C1*(GAMI-G1*dlog(1.d0+GAMI/G1))+
     +   C2/2.*dlog(1.d0+GAMI**2/G2)
      CVii=-0.5*dsqrt(GAMI)**3*(A1*A2/dsqrt(A2+GAMI)**3+
     +  A3*(1.d0-GAMI)/(1.d0+GAMI)**2) -
     -  GAMI**2*(C1*G1/(G1+GAMI)**2+C2*(G2-GAMI**2)/(G2+GAMI**2)**2)
      PION=UION/3.
      PDRii=(4.*UION-CVii)/9. ! p_{ii} + d p_{ii} / d ln\rho
      PDTii=CVii/3. ! p_{ii} + d p_{ii} / d ln T
      return
      end

      subroutine FSCRliq8(RS,GAME,Zion,
     *     FSCR,USCR,PSCR,CVSCR,PDTSCR,PDRSCR) ! fit to the el.-ion scr.
*                                                       Version 11.09.08
*                                                       cleaned 16.06.09
* Stems from FSCRliq7 v. 09.06.07. Included a check for RS=0.
*   INPUT: RS - density parameter, GAME - electron Coulomb parameter,
*          Zion - ion charge number,
*   OUTPUT: FSCR - screening (e-i) free energy per kT per 1 ion,
*           USCR - internal energy per kT per 1 ion (screen.contrib.)
*           PSCR - pressure divided by (n_i kT) (screen.contrib.)
*           CVSCR - heat capacity per 1 ion (screen.contrib.)
*           PDTSCR,PDRSCR = PSCR + d PSCR / d ln(T,\rho)
      implicit double precision(A-H),double precision(O-Z)
      save
      parameter(XRS=.0140047,TINY=1.d-19)
      if (RS.lt.0.) stop'FSCRliq8: RS < 0'
      if (RS.lt.TINY) then
         FSCR=0.
         USCR=0.
         PSCR=0.
         CVSCR=0.
         PDTSCR=0.
         PDRSCR=0.
         return
      endif
      SQG=sqrt(GAME)
      SQR=sqrt(RS)
      SQZ1=dsqrt(1.+Zion)
      SQZ=dsqrt(Zion)
      CDH0=Zion/1.73205 ! 1.73205=sqrt(3.)
      CDH=CDH0*(SQZ1**3-SQZ**3-1.)
      SQG=sqrt(GAME)
      ZLN=dlog(Zion)
      Z13=exp(ZLN/3.) ! Zion**(1./3.)
      X=XRS/RS ! relativity parameter
      CTF=Zion**2*.2513*(Z13-1.+.2/sqrt(Z13))
* Thomas-Fermi constant; .2513=(18/175)(12/\pi)^{2/3}
      P01=1.11*exp(.475*ZLN)
      P03=0.2+0.078*ZLN**2
      PTX=1.16+.08*ZLN
      TX=GAME**PTX
      TXDG=PTX*TX/GAME
      TXDGG=(PTX-1.)*TXDG/GAME
      TY1=1./(1.d-3*Zion**2+2.*GAME)
      TY1DG=-2.*TY1**2
      TY1DGG=-4.*TY1*TY1DG
      TY2=1.+6.*RS**2
      TY2DX=-12.*RS**2/X
      TY2DXX=-3.*TY2DX/X
      TY=RS**3/TY2*(1.+TY1)
      TYX=3./X+TY2DX/TY2
      TYDX=-TY*TYX
      TYDG=RS**3*TY1DG/TY2
      P1=(Zion-1.)/9.
      COR1=1.+P1*TY
      COR1DX=P1*TYDX
      COR1DG=P1*TYDG
      COR1DXX=P1*(TY*(3./X**2+(TY2DX/TY2)**2-TY2DXX/TY2)-TYDX*TYX)
      COR1DGG=P1*RS**3*TY1DGG/TY2
      COR1DXG=-P1*TYDG*TYX
      U0=.78*sqrt(GAME/Zion)*RS**3
      U0DX=-3.*U0/X
      U0DG=.5*U0/GAME
      U0DXX=-4.*U0DX/X
      U0DGG=-.5*U0DG/GAME
      U0DXG=-3.*U0DG/X
      D0DG=Zion**3
      D0=GAME*D0DG+21.*RS**3
      D0DX=-63.*RS**3/X
      D0DXX=252.*RS**3/X**2
      COR0=1.+U0/D0
      COR0DX=(U0DX-U0*D0DX/D0)/D0
      COR0DG=(U0DG-U0*D0DG/D0)/D0
      COR0DXX=(U0DXX-(2.*U0DX*D0DX+U0*D0DXX)/D0+2.*(D0DX/D0)**2)/D0
      COR0DGG=(U0DGG-2.*U0DG*D0DG/D0+2.*U0*(D0DG/D0)**2)/D0
      COR0DXG=(U0DXG-(U0DX*D0DG+U0DG*D0DX)/D0+2.*U0*D0DX*D0DG/D0**2)/D0
* Relativism:
      RELE=dsqrt(1.d0+X**2)
      Q1=.18/dsqrt(dsqrt(Zion))
      Q2=.2+.37/dsqrt(Zion)
      H1U=1.+X**2/5.
      H1D=1.+Q1*X+Q2*X**2
      H1=H1U/H1D
      H1X=.4*X/H1U-(Q1+2.*Q2*X)/H1D
      H1DX=H1*H1X
      H1DXX=H1DX*H1X+
     +  H1*(.4/H1U-(.4*X/H1U)**2-2.*Q2/H1D+((Q1+2.*Q2*X)/H1D)**2)
      UP=CDH*SQG+P01*CTF*TX*COR0*H1
      UPDX=P01*CTF*TX*(COR0DX*H1+COR0*H1DX)
      UPDG=.5*CDH/SQG+P01*CTF*(TXDG*COR0+TX*COR0DG)*H1
      UPDXX=P01*CTF*TX*(COR0DXX*H1+2.*COR0DX*H1DX+COR0*H1DXX)
      UPDGG=-.25*CDH/(SQG*GAME)+
     +  P01*CTF*(TXDGG*COR0+2.*TXDG*COR0DG+TX*COR0DGG)*H1
      UPDXG=P01*CTF*(TXDG*(COR0DX*H1+COR0*H1DX)+
     +  TX*(COR0DXG*H1+COR0DG*H1DX))
      DN1=P03*SQG+P01/RS*TX*COR1
      DN1DX=P01*TX*(COR1/XRS+COR1DX/RS)
      DN1DG=.5*P03/SQG+P01/RS*(TXDG*COR1+TX*COR1DG)
      DN1DXX=P01*TX/XRS*(2.*COR1DX+X*COR1DXX)
      DN1DGG=-.25*P03/(GAME*SQG)+
     +  P01/RS*(TXDGG*COR1+2.*TXDG*COR1DG+TX*COR1DGG)
      DN1DXG=P01*(TXDG*(COR1/XRS+COR1DX/RS)+TX*(COR1DG/XRS+COR1DXG/RS))
      DN=1.+DN1/RELE
      DNDX=DN1DX/RELE-X*DN1/RELE**3
      DNDXX=(DN1DXX-((2.*X*DN1DX+DN1)-3.*X**2*DN1/RELE**2)/RELE**2)/RELE
      DNDG=DN1DG/RELE
      DNDGG=DN1DGG/RELE
      DNDXG=DN1DXG/RELE-X*DN1DG/RELE**3
      FSCR=-UP/DN*GAME
      FX=(UP*DNDX/DN-UPDX)/DN
      FXDG=((UPDG*DNDX+UPDX*DNDG+UP*DNDXG-2.*UP*DNDX*DNDG/DN)/DN-
     -  UPDXG)/DN
      FDX=FX*GAME
      FG=(UP*DNDG/DN-UPDG)/DN
      FDG=FG*GAME-UP/DN
      FDGDH=SQG*DNDG/DN**2 ! d FDG / d CDH
      FDXX=((UP*DNDXX+2.*(UPDX*DNDX-UP*DNDX**2/DN))/DN-UPDXX)/DN*GAME
      FDGG=2.*FG+GAME*((2.*DNDG*(UPDG-UP*DNDG/DN)+UP*DNDGG)/DN-UPDGG)/DN
      FDXG=FX+GAME*FXDG
      USCR=GAME*FDG
      CVSCR=-GAME**2*FDGG
      PSCR=(X*FDX+GAME*FDG)/3.
      PDTSCR=-GAME**2*(X*FXDG+FDGG)/3.
      PDRSCR=(12.*PSCR+X**2*FDXX+2.*X*GAME*FDXG+GAME**2*FDGG)/9.
      return
      end

* ==============   SUBROUTINES FOR THE SOLID STATE   ================= *
      subroutine FSCRsol8(RS,GAMI,Zion,TPT,
     *     FSCR,USCR,PSCR,S_SCR,CVSCR,PDTSCR,PDRSCR)
*                                                       Version 28.05.08
* Fit to the el.-ion screening in bcc or fcc Coulomb solid
* Stems from FSCRsol8 v.09.06.07. Included a check for RS=0.
*   INPUT: RS - el. density parameter, GAMI - ion coupling parameter,
*          Zion - ion charge, TPT=T_p/T - ion quantum parameter
*   OUTPUT: FSCR - screening (e-i) free energy per kT per 1 ion,
*           USCR - internal energy per kT per 1 ion (screen.contrib.)
*           PSCR - pressure divided by (n_i kT) (screen.contrib.)
*           S_SCR - screening entropy contribution / (N_i k)
*           CVSCR - heat capacity per 1 ion (screen.contrib.)
*           PDTSCR,PDRSCR = PSCR + d PSCR / d ln(T,\rho)
      implicit double precision(A-H),double precision(O-Z)
      save
      dimension AP(4) ! parameters of the fit
      parameter (ENAT=2.7182818285d0,TINY=1.d-19)
      data AP/1.1866,.684,17.9,41.5/,PX/.205/ ! for bcc lattice
cc      data AP/1.1857,.663,17.1,40./,PX/.212/ ! for fcc lattice
      if (RS.lt.0.) stop'FSCRliq8: RS < 0'
      if (RS.lt.TINY) then
         FSCR=0.
         USCR=0.
         PSCR=0.
         S_SCR=0.
         CVSCR=0.
         PDTSCR=0.
         PDRSCR=0.
         return
      endif
      XSR=.0140047/RS ! relativity parameter
      Z13=Zion**.3333333
      P1=.00352*(1.-AP(1)/Zion**.267+.27/Zion)
      P2=1.+2.25/Z13*
     *(1.+AP(2)*Zion**5+.222*Zion**6)/(1.+.222*Zion**6)
      ZLN=dlog(Zion)
      Finf=sqrt(P2/XSR**2+1.)*Z13**2*P1 ! The TF limit
      FinfX=-P2/((P2+XSR**2)*XSR)
      FinfDX=Finf*FinfX
      FinfDXX=FinfDX*FinfX-FinfDX*(P2+3.*XSR**2)/((P2+XSR**2)*XSR)
      R1=AP(4)/(1.+ZLN)
      R2=.395*ZLN+.347/Zion/sqrt(Zion)
      R3=1./(1.+ZLN*sqrt(ZLN)*.01+.097/Zion**2)
      Q1U=R1+AP(3)*XSR**2
      Q1D=1.+R2*XSR**2
      Q1=Q1U/Q1D
      Q1X=2.*XSR*(AP(3)/Q1U-R2/Q1D)
      Q1XDX=Q1X/XSR+4.*XSR**2*((R2/Q1D)**2-(AP(3)/Q1U)**2)
      Q1DX=Q1*Q1X
      Q1DXX=Q1DX*Q1X+Q1*Q1XDX
C* Old (PC'2000) quantum suppression factor and its derivatives:
C      SUP=sqrt(1.+(PX*TPT)**2) ! suppression of the Gamma-dependence
C      SUPDX=(PX*TPT)**2/SUP*.5/XSR
C      SUPDG=(PX*TPT)**2/SUP/GAMI
C      SUPDXX=-SUPDX**2/SUP
C      SUPDGG=SUPDG*(1./GAMI-SUPDG/SUP)
C      SUPDXG=SUPDX*(2./GAMI-SUPDG/SUP)
* New quantum factor, in order to suppress CVSCR at TPT >> 1
      if (TPT.lt.6./PX) then
         Y0=(PX*TPT)**2
         Y0DX=Y0/XSR
         Y0DG=2.*Y0/GAMI
         Y0DXX=0.
         Y0DGG=Y0DG/GAMI
         Y0DXG=Y0DG/XSR
         Y1=dexp(Y0)
         Y1DX=Y1*Y0DX
         Y1DG=Y1*Y0DG
         Y1DXX=Y1*(Y0DX**2+Y0DXX)
         Y1DGG=Y1*(Y0DG**2+Y0DGG)
         Y1DXG=Y1*(Y0DX*Y0DG+Y0DXG)
         SA=1.d0+Y1
         SUPA=dlog(SA)
         SUPADX=Y1DX/SA
         SUPADG=Y1DG/SA
         SUPADXX=(Y1DXX-Y1DX**2/SA)/SA
         SUPADGG=(Y1DGG-Y1DG**2/SA)/SA
         SUPADXG=(Y1DXG-Y1DX*Y1DG/SA)/SA
         EM2=ENAT-2.d0
         SB=ENAT-EM2/Y1
         SUPB=dlog(SB)
         EM2Y1=EM2/(Y1**2*SB)
         SUPBDX=EM2Y1*Y1DX
         SUPBDG=EM2Y1*Y1DG
         SUPBDXX=EM2Y1*(Y1DXX-2.*Y1DX**2/Y1-Y1DX*SUPBDX)
         SUPBDGG=EM2Y1*(Y1DGG-2.*Y1DG**2/Y1-Y1DG*SUPBDG)
         SUPBDXG=EM2Y1*(Y1DXG-2.*Y1DX*Y1DG/Y1-Y1DG*SUPBDX)
         SUP=dsqrt(SUPA/SUPB)
         SUPX=.5*(SUPADX/SUPA-SUPBDX/SUPB)
         SUPDX=SUP*SUPX
         SUPG=.5*(SUPADG/SUPA-SUPBDG/SUPB)
         SUPDG=SUP*SUPG
         SUPDXX=SUPDX*SUPX+
     +     SUP*.5*(SUPADXX/SUPA-(SUPADX/SUPA)**2-
     -             SUPBDXX/SUPB+(SUPBDX/SUPB)**2)
         SUPDGG=SUPDG*SUPG+
     +     SUP*.5*(SUPADGG/SUPA-(SUPADG/SUPA)**2-
     -             SUPBDGG/SUPB+(SUPBDG/SUPB)**2)
         SUPDXG=SUPDX*SUPG+
     +     SUP*.5*((SUPADXG-SUPADX*SUPADG/SUPA)/SUPA-
     -             (SUPBDXG-SUPBDX*SUPBDG/SUPB)/SUPB)
      else
         SUP=PX*TPT
         SUPDX=.5*PX*TPT/XSR
         SUPDG=PX*TPT/GAMI
         SUPDXX=-.5*SUPDX/XSR
         SUPDGG=0.
         SUPDXG=SUPDX/GAMI
      endif
      GR3=(GAMI/SUP)**R3
      GR3X=-R3*SUPDX/SUP
      GR3DX=GR3*GR3X
      GR3DXX=GR3DX*GR3X-R3*GR3*(SUPDXX/SUP-(SUPDX/SUP)**2)
      GR3G=R3*(1./GAMI-SUPDG/SUP)
      GR3DG=GR3*GR3G
      GR3DGG=GR3DG*GR3G+GR3*R3*((SUPDG/SUP)**2-SUPDGG/SUP-1./GAMI**2)
      GR3DXG=GR3DG*GR3X+GR3*R3*(SUPDX*SUPDG/SUP**2-SUPDXG/SUP)
      W=1.+Q1/GR3
      WDX=Q1DX/GR3-Q1*GR3DX/GR3**2
      WDG=-Q1*GR3DG/GR3**2
      WDXX=Q1DXX/GR3-(2.*Q1DX*GR3DX+Q1*(GR3DXX-2.*GR3DX**2/GR3))/GR3**2
      WDGG=Q1*(2.*GR3DG**2/GR3-GR3DGG)/GR3**2
      WDXG=Q1DXG/GR3-(Q1DX*GR3DG+Q1*(GR3DXG-2.*GR3DX*GR3DG/GR3))/GR3**2
      FSCR=-GAMI*Finf*W
      FDX=-GAMI*(FinfDX*W+Finf*WDX)
      FDXX=-GAMI*(FinfDXX*W+2.*FinfDX*WDX+Finf*WDXX)
      FDG=-Finf*W-GAMI*Finf*WDG
      FDGG=-2.*Finf*WDG-GAMI*Finf*WDGG
      FDXG=-FinfDX*W-Finf*WDX-GAMI*(FinfDX*WDG+Finf*WDXG)
      S_SCR=-GAMI**2*Finf*WDG
      USCR=S_SCR+FSCR
      CVSCR=-GAMI**2*FDGG
      PSCR=(XSR*FDX+GAMI*FDG)/3.
      PDTSCR=GAMI**2*(XSR*Finf*(FinfX*WDG+WDXG)-FDGG)/3.
      PDRSCR=(12.*PSCR+XSR**2*FDXX+2.*XSR*GAMI*FDXG+
     +  GAMI**2*FDGG)/9.
      return
      end

      subroutine ANHARM8(GAMI,TPT,Fah,Uah,Pah,CVah,PDTah,PDRah)
* ANHARMONIC free energy                                Version 27.07.07
*                                                       cleaned 16.06.09
* Stems from ANHARM8b. Difference: AC=0., B1=.12 (.1217 - over accuracy)
* Input: GAMI - ionic Gamma, TPT=Tp/T - ionic quantum parameter
* Output: anharm.free en. Fah=F_{AH}/(N_i kT), internal energy Uah,
*   pressure Pah=P_{AH}/(n_i kT), specific heat CVah = C_{V,AH}/(N_i k),
*   PDTah = Pah + d Pah / d ln T, PDRah = Pah + d Pah / d ln\rho
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter(NM=3)
      dimension AA(NM)
      data AA/10.9,247.,1.765d5/ ! Farouki & Hamaguchi'93
      data B1/.12/ ! coeff.at \eta^2/\Gamma at T=0
      CK=B1/AA(1) ! fit coefficient
      TPT2=TPT**2
      TPT4=TPT2**2
      TQ=B1*TPT2/GAMI ! quantum dependence
      TK2=CK*TPT2
      SUP=dexp(-TK2) ! suppress.factor of class.anharmonicity
      Fah=0.
      Uah=0.
      Pah=0.
      CVah=0.
      PDTah=0.
      PDRah=0.
      SUPGN=SUP
      do N=1,NM
         CN=N
         SUPGN=SUPGN/GAMI ! SUP/Gamma^n
         ACN=AA(N)
         Fah=Fah-ACN/CN*SUPGN
         Uah=Uah+(ACN*(1.+2.*TK2/CN))*SUPGN
         PN=AA(N)/3.+TK2*AA(N)/CN
         Pah=Pah+PN*SUPGN
         CVah=CVah+((CN+1.)*AA(N)+(4.-2./CN)*AA(N)*TK2+
     +     4.*AA(N)*CK**2/CN*TPT4)*SUPGN
         PDTah=PDTah+(PN*(1.+CN+2.*TK2)-2./CN*AA(N)*TK2)*SUPGN
         PDRah=PDRah+(PN*(1.-CN/3.-TK2)+AA(N)/CN*TK2)*SUPGN
      enddo
      Fah=Fah-TQ
      Uah=Uah-TQ
      Pah=Pah-TQ/1.5
      PDRah=PDRah-TQ/4.5
      return
      end

      subroutine FHARM8(GAMI,TPT,
     *   Fharm,Uharm,Pharm,CVth,Sharm,PDTharm,PDRharm)
* Thermodynamic functions of a harmonic crystal, incl.stat.Coul.lattice
*                                                       Version 15.02.08
* Stems from FHARM7 v.09.06.07.
* Replaced FthCHA7 by the more accurate fit HLfit8.
* Input: GAMI - ionic Gamma, TPT=T_{p,i}/T
* Output: Fharm=F/(N_i T), Uharm=U/(N_i T), Pharm=P/(n_i T),
* CVth=C_V/N_i, Sharm=S/N_i
* PDTharm = Pharm + d Pharm / d ln T, PDRharm = Pharm + d Pharm/d ln\rho
      implicit double precision (A-H), double precision (O-Z)
      save
      data CM/.895929256d0/ ! Madelung
      call HLfit8(TPT,Fth,Uth,CVth,Sth,1)
      U0=-CM*GAMI ! perfect lattice
      Uharm=U0+Uth
      Fharm=U0+Fth
      Sharm=Sth
      Pharm=U0/3.+Uth/2.
      PDTharm=.5*CVth
      PDRharm=U0/2.25+.75*Uth-.25*CVth
      return
      end

      subroutine HLfit8(eta,F,U,CV,S,LATTICE)
*                                                       Version 03.12.08
* Fit to thermal part of the thermodynamic functions.
* Baiko, Potekhin, & Yakovlev (2001). Stems from HLfit v.20.03.01.
* Zero-point lattice quantum energy 1.5u_1\eta INCLUDED (unlike HLfit).
* Input: eta=Tp/T, LATTICE=1 for bcc, 2 for fcc
* Output: F and U (normalized to NkT),
*   CV and S (normalized to Nk) in the HL model for bcc Coulomb lattice
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter(EPS=1.d-5,TINY=1.d-99)
      if (LATTICE.eq.1) then ! bcc lattice
         CLM=-2.49389d0 ! 3*ln<\omega/\omega_p>
         U1=.5113875d0
         ALPHA=.265764d0
         BETA=.334547d0
         GAMMA=.932446d0
         A1=.1839d0
         A2=.593586d0
         A3=.0054814d0
         A4=5.01813d-4
         A6=3.9247d-7
         A8=5.8356d-11
         B0=261.66d0
         B2=7.07997d0
         B4=.0409484d0
         B5=.000397355d0
         B6=5.11148d-5
         B7=2.19749d-6
         C9=.004757014d0
         C11=.0047770935d0
      elseif (LATTICE.eq.2) then ! fcc lattice
         CLM=-2.45373d0
         U1=.513194d0
         ALPHA=.257591d0
         BETA=.365284d0
         GAMMA=.9167070d0
         A1=.0
         A2=.532535d0
         A3=.0
         A4=3.76545d-4
         A6=2.63013d-7
         A8=6.6318d-11
         B0=303.20d0
         B2=7.7255d0
         B4=.0439597d0
         B5=.000114295d0
         B6=5.63434d-5
         B7=1.36488d-6
         C9=.00492387d0
         C11=.00437506d0
      else
         stop'HLfit: unknown lattice type'
      endif
      if (eta.gt.1./EPS) then ! asymptote of Eq.(13) of BPY'01
         U=3./(C11*eta**3)
         F=-U/3.
         CV=4.*U
        goto 50
      elseif (eta.lt.EPS) then ! Eq.(17) of BPY'01
        if (eta.lt.TINY) stop'HLfit8: eta is too small'
         F=3.*dlog(eta)+CLM-1.5*U1*eta+eta**2/24. 
         U=3.-1.5*U1*eta+eta**2/12.
         CV=3.-eta**2/12.
         goto 50
      endif
      B9=A6*C9
      B11=A8*C11
      UP=1.+A1*eta+A2*eta**2+A3*eta**3+A4*eta**4+A6*eta**6+A8*eta**8
      DN=B0+B2*eta**2+B4*eta**4+B5*eta**5+B6*eta**6+
     +  B7*eta**7+B9*eta**9+B11*eta**11
      EA=dexp(-ALPHA*eta)
      EB=dexp(-BETA*eta)
      EG=dexp(-GAMMA*eta)
      F=dlog(1.d0-EA)+dlog(1.d0-EB)+dlog(1.-EG)-UP/DN ! F_{thermal}/NT
      UP1=A1+
     + 2.*A2*eta+3.*A3*eta**2+4.*A4*eta**3+6.*A6*eta**5+8.*A8*eta**7
      UP2=2.*A2+6.*A3*eta+12.*A4*eta**2+30.*A6*eta**4+56.*A8*eta**6
      DN1=2.*B2*eta+4.*B4*eta**3+5.*B5*eta**4+6.*B6*eta**5+
     +  7.*B7*eta**6+9.*B9*eta**8+11.*B11*eta**10.
      DN2=2.*B2+12.*B4*eta**2+20.*B5*eta**3+30.*B6*eta**4+
     +  42.*B7*eta**5+72.*B9*eta**7+110.*B11*eta**9
      U=ALPHA*EA/(1.d0-EA)+BETA*EB/(1.d0-EB)+GAMMA*EG/(1.d0-EG)-
     -  (UP1*DN-DN1*UP)/DN**2 ! int.en./NT/eta
      CV=ALPHA**2*EA/(1.d0-EA)**2+BETA**2*EB/(1.d0-EB)**2+
     +  GAMMA**2*EG/(1.d0-EG)**2+
     +  ((UP2*DN-DN2*UP)*DN-2.*(UP1*DN-DN1*UP)*DN1)/DN**3 ! cV/eta^2
      U=U*eta
      CV=CV*eta**2
   50 continue
      S=U-F
* Add zero-point lattice energy:
      E0=1.5*U1*eta
      U=U+E0
      F=F+E0
      return
      end

      subroutine CORMIX(RS,GAME,Zmean,Z2mean,Z52,Z53,Z321,
     *  FMIX,UMIX,PMIX,CVMIX,PDTMIX,PDRMIX)
*                                                       Version 02.07.09
* Correction to the linear mixing rule for moderate to small Gamma
* Input: RS=r_s (if RS=0, then OCP, otherwise EIP)
*        GAME=\Gamma_e
*        Zmean=<Z> (average Z of all ions, without electrons)
*        Z2mean=<Z^2>, Z52=<Z^2.5>, Z53=<Z^{5/3}>, Z321=<Z(Z+1)^1.5>
* Output: FMIX=\Delta f - corr.to the reduced free energy f=F/N_{ion}kT
*         UMIX=\Delta u - corr.to the reduced internal energy u
*         PMIX=\Delta u - corr.to the reduced pressure P=P/n_{ion}kT
*         CVMIX=\Delta c - corr.to the reduced heat capacity c_V
*         PDTMIX=(1/n_{ion}kT)d\Delta P / d ln T
*               = \Delta p +  d \Delta p / d ln T
*         PDRMIX=(1/n_{ion}kT)d\Delta P / d ln n_e
* (composition is assumed fixed: Zmean,Z2mean,Z52,Z53=constant)
      implicit double precision (A-H), double precision (O-Z)
      parameter (TINY=1.d-9)
      GAMImean=GAME*Z53
      if (RS.lt.TINY) then ! OCP
         Dif0=Z52-dsqrt(Z2mean**3/Zmean)
      else
         Dif0=Z321-dsqrt((Z2mean+Zmean)**3/Zmean)
      endif
      DifR=Dif0/Z52
      DifFDH=Dif0*GAME*sqrt(GAME/3.) ! F_DH - F_LM(DH)
      D=Z2mean/Zmean**2
      if (dabs(D-1.d0).lt.TINY) then ! no correction
         FMIX=0.
         UMIX=0.
         PMIX=0.
         CVMIX=0.
         PDTMIX=0.
         PDRMIX=0.
         return
      endif
      P3=D**(-0.2)
      D0=(2.6*DifR+14.*DifR**3)/(1.d0-P3)
      GP=D0*GAMImean**P3
      FMIX0=DifFDH/(1.+GP)
      Q=D**2*.0117
      R=1.5/P3-1.
      GQ=Q*GP
      FMIX=FMIX0/(1.+GQ)**R
      G=1.5-P3*GP/(1.+GP)-R*P3*GQ/(1.+GQ)
      UMIX=FMIX*G
      PMIX=UMIX/3.d0
      GDG=-P3**2*(GP/(1.d0+GP)**2+R*GQ/(1.d0+GQ)**2) ! d G /d ln Gamma
      UDG=UMIX*G+FMIX*GDG ! d u_mix /d ln Gamma
      CVMIX=UMIX-UDG
      PDTMIX=PMIX-UDG/3.
      PDRMIX=PMIX+UDG/9.
      return
      end

* ===================  IDEAL ELECTRON GAS  =========================== *
      subroutine ELECT9(TEMP,CHI,
     *  DENS,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE,DDN)
*                                                       Version 04.03.09
* Compared to ELECTRON v.06.07.00, this S/R is completely rewritten: 
*        numerical differentiation is avoided now.
* Compared to ELECT7 v.06.06.07,
*    - call BLIN7 is changed to call BLIN8,
*    - Sommerfeld expansion is used at chi > 14 i.o. 1.e4
*    - Sommerfeld expansion is corrected: introduced DeltaEF, D1 and D2
* Ideal electron-gas EOS.
* Input: TEMP - T [a.u.], CHI=\mu/T
* Output: DENS - electron number density n_e [a.u.],
*         FEid - free energy / N_e kT, UEid - internal energy / N_e kT,
*         PEid - pressure (P_e) / n_e kT, SEid - entropy / N_e k,
*         CVE - heat capacity / N_e k,
*         CHITE=(d ln P_e/d ln T)_V, CHIRE=(d ln P_e/d ln n_e)_T
*         DDN=(d ln n_e/d CHI)_T = (T/n_e) (d n_e/d\mu)_T
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter (BOHR=137.036,PI=3.141592653d0)
      parameter (PI2=PI**2,BOHR2=BOHR**2,BOHR3=BOHR2*BOHR) !cleaned 15/6
      TEMR=TEMP/BOHR2 ! T in rel.units (=T/mc^2)
      if (CHI.lt.14.d0) then ! use Fermi-Dirac integrals(14.d0:16.02.09)
         call BLIN8(TEMR,CHI,
     *     W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT,
     *     W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT,
     *     W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT,
     *     W0XXX,W0XTT,W0XXT)
         TPI=TEMR*dsqrt(2.d0*TEMR)/PI2 ! common pre-factor
         DENR=TPI*(W1*TEMR+W0)
         PR=TEMR*TPI/3.*(W2*TEMR+2.*W1)
         U=TEMR*TPI*(W2*TEMR+W1)
* (these are density, pressure, and internal energy in the rel.units)
         PEid=PR/(DENR*TEMR)
         UEid=U/(DENR*TEMR)
         FEid=CHI-PEid
         DENS=DENR*BOHR3 ! converts from rel.units to a.u.
         SEid=UEid-FEid
* derivatives over T at constant chi:
         dndT=TPI*(1.5*W0/TEMR+2.5*W1+W0DT+TEMR*W1DT) ! (d n_e/dT)_\chi
         dPdT=TPI/3.*(5.*W1+2.*TEMR*W1DT+3.5*TEMR*W2+TEMR**2*W2DT)!dP/dT
         dUdT=TPI*(2.5*W1+TEMR*W1DT+3.5*TEMR*W2+TEMR**2*W2DT)!dU/dT_\chi
* derivatives over chi at constant T:
         dndH=TPI*(W0DX+TEMR*W1DX) ! (d n_e/d\chi)_T
         DDN=dndH/DENR ! (d ln n_e/d\chi)_T
         dPdH=TPI/3.*TEMR*(2.*W1DX+TEMR*W2DX) ! (d P_e/d\chi)_T
         dUdH=TPI*TEMR*(W1DX+TEMR*W2DX) ! (d U_e/d\chi)_T
         CVE=(dUdT-dUdH*dndT/dndH)/DENR
         CHITE=TEMR/PR*(dPdT-dPdH*dndT/dndH)
         CHIRE=DENR/PR*dPdH/dndH ! (dndH*TEMR*PEid) ! DENS/PRE*dPdH/dndH
      else ! use the Sommerfeld expansion at very large CHI
         EF=CHI*TEMR ! Fermi energy in mc^2 - zeroth aprox. = CMU1
         DeltaEF=PI2*TEMR**2/6.d0*(1.d0+2.d0*EF*(2.d0+EF))/
     /     (EF*(1.d0+EF)*(2.d0+EF)) ! corr. [page 125]
         EF=EF+DeltaEF ! corrected Fermi energy (14.02.09)
         G=1.d0+EF ! electron Lorentz-factor
        if (EF.gt.1.d-5) then ! relativistic expansion (Yak.&Shal.'89)
           PF=dsqrt(G**2-1.d0) ! Fermi momentum [rel.un.=mc]
           F=(PF*(1.+2.d0*PF**2)*G-PF**3/.375d0-dlog(PF+G))/8.d0/PI2!F/V
           DF=-TEMR**2*PF*G/6.d0 ! thermal correction to F/V
           P=(PF*G*(PF**2/1.5d0-1.d0)+dlog(PF+G))/8.d0/PI2 ! P(T=0)
           DP=TEMR**2*PF*(PF**2+2.d0)/G/18.d0 ! thermal correction to P
           CVE=PI2*TEMR*G/PF**2
        else ! nonrelativistic limit
           PF=dsqrt(2.d0*EF)
           F=PF**5*0.1d0/PI2
           DF=-TEMR**2*PF/6.d0
           P=F/1.5d0
           DP=TEMR**2*PF/9.d0
           CVE=PI2*TEMR/EF/2.d0
        endif
         F=F+DF
         P=P+DP
         S=-2.d0*DF ! entropy per unit volume [rel.un.]
         U=F+S
         CHIRE=PF**5/(9.d0*PI2*P*G)
         CHITE=2.d0*DP/P
         DENR=PF**3/3.d0/PI2 ! n_e [rel.un.=\Compton^{-3}]
         DENS=DENR*BOHR3 ! conversion to a.u.(=\Bohr_radius^{-3})
         CMUDLN=PF**2/(3.d0*G)+ ! d\mu/d ln n
     +     PI2*TEMR**2*(1.d0+2.d0*(PF*G)**2)/(9.d0*(PF*G**2)**2) !T-cor.
         DDN=TEMR/CMUDLN ! d ln n / d\chi = (T/n)*(dn/d\mu)
         DT=DENR*TEMR
         PEid=P/DT
         UEid=U/DT
         FEid=F/DT
         SEid=S/DT
* Empirical corrections of 16.02.09:
         D1=DeltaEF/EF
         D2=D1*(4.d0-2.d0*(PF/G))
         CVE=CVE/(1.d0+D2)
         SEid=SEid/(1.d0+D1)
         CHITE=CHITE/(1.d0+D2)
      endif
      return
      end

* ==============  ELECTRON EXCHANGE AND CORRELATION   ================ *
      subroutine EXCOR7(RS,GAME,FXC,UXC,PXC,CVXC,SXC,PDTXC,PDRXC)
*                                                       Version 09.06.07
* Exchange-correlation contribution for the electron gas
* Stems from TANAKA1 v.03.03.96. Added derivatives.
* Input: RS - electron density parameter =electron-sphere radius in a.u.
*        GAME - electron Coulomb coupling parameter
* Output: FXC - excess free energy of e-liquid per kT per one electron
*             according to Tanaka & Ichimaru 85-87 and Ichimaru 93 FXC,
*         UXC - free and internal energy contr.[per 1 electron, kT]
*         PXC - pressure contribution divided by (n_e kT)
*         CVXC - heat capacity divided by N_e k
*         SXC - entropy divided by N_e k
*         PDTXC,PDRXC = PXC + d ln PXC / d ln(T,\rho)
      implicit double precision(A-H),double precision(O-Z)
      save
      THETA=.543*RS/GAME ! non-relativistic degeneracy parameter
      SQTH=dsqrt(THETA)
      THETA2=THETA**2
      THETA3=THETA2*THETA
      THETA4=THETA3*THETA
      if (THETA.gt..005) then
         CHT1=dcosh(1.d0/THETA)
         SHT1=dsinh(1.d0/THETA)
         CHT2=dcosh(1.d0/SQTH)
         SHT2=dsinh(1.d0/SQTH)
         T1=SHT1/CHT1 ! dtanh(1.d0/THETA)
         T2=SHT2/CHT2 ! dtanh(1./dsqrt(THETA))
         T1DH=-1./(THETA*CHT1)**2 ! d T1 / d\theta
         T1DHH=2./(THETA*CHT1)**3*(CHT1-SHT1/THETA)
         T2DH=-.5*SQTH/(THETA*CHT2)**2
         T2DHH=(.75*SQTH*CHT2-.5*SHT2)/(THETA*CHT2)**3
      else
         T1=1.
         T2=1.
         T1DH=0.
         T2DH=0.
         T1DHH=0.
         T2DHH=0.
      endif
      A0=.75+3.04363*THETA2-.09227*THETA3+1.7035*THETA4
      A0DH=6.08726*THETA-.27681*THETA2+6.814*THETA3
      A0DHH=6.08726-.55362*THETA+20.442*THETA2
      A1=1.+8.31051*THETA2+5.1105*THETA4
      A1DH=16.62102*THETA+20.442*THETA3
      A1DHH=16.62102+61.326*THETA2
      A=.610887*A0/A1*T1 ! HF fit of Perrot and Dharma-wardana
      AH=A0DH/A0-A1DH/A1+T1DH/T1
      ADH=A*AH
      ADHH=ADH*AH+A*(A0DHH/A0-(A0DH/A0)**2-A1DHH/A1+(A1DH/A1)**2+
     +  T1DHH/T1-(T1DH/T1)**2)
      B0=.341308+12.070873d0*THETA2+1.148889d0*THETA4
      B0DH=24.141746d0*THETA+4.595556d0*THETA3
      B0DHH=24.141746d0+13.786668d0*THETA2
      B1=1.+10.495346d0*THETA2+1.326623*THETA4
      B1DH=20.990692d0*THETA+5.306492*THETA3
      B1DHH=20.990692d0+15.919476d0*THETA2
      B=SQTH*T2*B0/B1
      BH=.5/THETA+T2DH/T2+B0DH/B0-B1DH/B1
      BDH=B*BH
      BDHH=BDH*BH+B*(-.5/THETA2+T2DHH/T2-(T2DH/T2)**2+
     +  B0DHH/B0-(B0DH/B0)**2-B1DHH/B1+(B1DH/B1)**2)
      D0=.614925+16.996055d0*THETA2+1.489056*THETA4
      D0DH=33.99211d0*THETA+5.956224d0*THETA3
      D0DHH=33.99211d0+17.868672d0*THETA2
      D1=1.+10.10935*THETA2+1.22184*THETA4
      D1DH=20.2187*THETA+4.88736*THETA3
      D1DHH=20.2187+14.66208*THETA2
      D=SQTH*T2*D0/D1
      DH=.5/THETA+T2DH/T2+D0DH/D0-D1DH/D1
      DDH=D*DH
      DDHH=DDH*DH+D*(-.5/THETA2+T2DHH/T2-(T2DH/T2)**2+
     +  D0DHH/D0-(D0DH/D0)**2-D1DHH/D1+(D1DH/D1)**2)
      E0=.539409+2.522206*THETA2+.178484*THETA4
      E0DH=5.044412*THETA+.713936*THETA3
      E0DHH=5.044412+2.141808*THETA2
      E1=1.+2.555501*THETA2+.146319*THETA4
      E1DH=5.111002*THETA+.585276*THETA3
      E1DHH=5.111002+1.755828*THETA2
      E=THETA*T1*E0/E1
      EH=1./THETA+T1DH/T1+E0DH/E0-E1DH/E1
      EDH=E*EH
      EDHH=EDH*EH+E*(T1DHH/T1-(T1DH/T1)**2+E0DHH/E0-(E0DH/E0)**2-
     -  E1DHH/E1+(E1DH/E1)**2-1./THETA2)
      EXP1TH=dexp(-1./THETA)
      C=(.872496+.025248*EXP1TH)*E
      CDH=.025248*EXP1TH/THETA2*E+C*EDH/E
      CDHH=.025248*EXP1TH/THETA2*(EDH+(1.-2.*THETA)/THETA2*E)+
     +  CDH*EDH/E+C*EDHH/E-C*(EDH/E)**2
      DISCR=dsqrt(4.*E-D**2)
      DIDH=.5/DISCR*(4.*EDH-2.*D*DDH)
      DIDHH=(-((2.*EDH-D*DDH)/DISCR)**2+2.*EDHH-DDH**2-D*DDHH)/DISCR
      S1=-C/E*GAME
      S1H=CDH/C-EDH/E
      S1DH=S1*S1H
      S1DHH=S1DH*S1H+S1*(CDHH/C-(CDH/C)**2-EDHH/E+(EDH/E)**2)
      S1DG=-C/E ! => S1DGG=0
      S1DHG=S1DG*(CDH/C-EDH/E)
      B2=B-C*D/E
      B2DH=BDH-(CDH*D+C*DDH)/E+C*D*EDH/E**2
      B2DHH=BDHH-(CDHH*D+2.*CDH*DDH+C*DDHH)/E+
     +  (2.*(CDH*D+C*DDH-C*D*EDH/E)*EDH+C*D*EDHH)/E**2
      SQGE=dsqrt(GAME)
      S2=-2./E*B2*SQGE
      S2H=B2DH/B2-EDH/E
      S2DH=S2*S2H
      S2DHH=S2DH*S2H+S2*(B2DHH/B2-(B2DH/B2)**2-EDHH/E+(EDH/E)**2)
      S2DG=.5*S2/GAME
      S2DGG=-.5*S2DG/GAME
      S2DHG=.5*S2DH/GAME
      R3=E*GAME+D*SQGE+1.
      R3DH=EDH*GAME+DDH*SQGE
      R3DHH=EDHH*GAME+DDHH*SQGE
      R3DG=E+.5*D/SQGE
      R3DGG=-.25*D/(GAME*SQGE)
      R3DHG=EDH+.5*DDH/SQGE
      B3=A-C/E
      B3DH=ADH-CDH/E+C*EDH/E**2
      B3DHH=ADHH-CDHH/E+(2.*CDH*EDH+C*EDHH)/E**2-2.*C*EDH**2/E**3
      C3=(D/E*B2-B3)/E ! =D*B2/E**2-B3/E
      C3DH=(DDH*B2+D*B2DH+B3*EDH)/E**2-2.*D*B2*EDH/E**3-B3DH/E
      C3DHH=(-B3DHH+
     +  (DDHH*B2+2.*DDH*B2DH+D*B2DHH+B3DH*EDH+B3*EDHH+B3DH*EDH)/E-
     -  2.*((DDH*B2+D*B2DH+B3*EDH+DDH*B2+D*B2DH)*EDH+D*B2*EDHH)/E**2+
     +  6.*D*B2*EDH**2/E**3)/E
      S3=C3*dlog(R3)
      S3DH=S3*C3DH/C3+C3*R3DH/R3
      S3DHH=(S3DH*C3DH+S3*C3DHH)/C3-S3*(C3DH/C3)**2+
     +  (C3DH*R3DH+C3*R3DHH)/R3-C3*(R3DH/R3)**2
      S3DG=C3*R3DG/R3
      S3DGG=C3*(R3DGG/R3-(R3DG/R3)**2)
      S3DHG=(C3DH*R3DG+C3*R3DHG)/R3-C3*R3DG*R3DH/R3**2
      B4=2.-D**2/E
      B4DH=EDH*(D/E)**2-2.*D*DDH/E
      B4DHH=EDHH*(D/E)**2+2.*EDH*(D/E)**2*(DDH/D-EDH/E)-
     -  2.*(DDH**2+D*DDHH)/E+2.*D*DDH*EDH/E**2
      C4=2.*E*SQGE+D
      C4DH=2.*EDH*SQGE+DDH
      C4DHH=2.*EDHH*SQGE+DDHH
      C4DG=E/SQGE
      C4DGG=-.5*E/(GAME*SQGE)
      C4DHG=EDH/SQGE
      S4A=2./E/DISCR
      S4AH=EDH/E+DIDH/DISCR
      S4ADH=-S4A*S4AH
      S4ADHH=-S4ADH*S4AH-
     -  S4A*(EDHH/E-(EDH/E)**2+DIDHH/DISCR-(DIDH/DISCR)**2)
      S4B=D*B3+B4*B2
      S4BDH=DDH*B3+D*B3DH+B4DH*B2+B4*B2DH
      S4BDHH=DDHH*B3+2.*DDH*B3DH+D*B3DHH+B4DHH*B2+2.*B4DH*B2DH+B4*B2DHH
      S4C=datan(C4/DISCR)-datan(D/DISCR)
      UP1=C4DH*DISCR-C4*DIDH
      DN1=DISCR**2+C4**2
      UP2=DDH*DISCR-D*DIDH
      DN2=DISCR**2+D**2
      S4CDH=UP1/DN1-UP2/DN2
      S4CDHH=(C4DHH*DISCR-C4*DIDHH)/DN1-
     -  UP1*2.*(DISCR*DIDH+C4*C4DH)/DN1**2-
     -  (DDHH*DISCR-D*DIDHH)/DN2+UP2*2.*(DISCR*DIDH+D*DDH)/DN2**2
      S4CDG=C4DG*DISCR/DN1
      S4CDGG=C4DGG*DISCR/DN1-2.*C4*DISCR*(C4DG/DN1)**2
      S4CDHG=(C4DHG*DISCR+C4DG*DIDH-
     -  C4DG*DISCR/DN1*2.*(DISCR*DIDH+C4*C4DH))/DN1
      S4=S4A*S4B*S4C
      S4DH=S4ADH*S4B*S4C+S4A*S4BDH*S4C+S4A*S4B*S4CDH
      S4DHH=S4ADHH*S4B*S4C+S4A*S4BDHH*S4C+S4A*S4B*S4CDHH+
     +  2.*(S4ADH*S4BDH*S4C+S4ADH*S4B*S4CDH+S4A*S4BDH*S4CDH)
      S4DG=S4A*S4B*S4CDG
      S4DGG=S4A*S4B*S4CDGG
      S4DHG=S4A*S4B*S4CDHG+S4CDG*(S4ADH*S4B+S4A*S4BDH)
      FXC=S1+S2+S3+S4
      FXCDH=S1DH+S2DH+S3DH+S4DH
      FXCDG=S1DG+S2DG+S3DG+S4DG
      FXCDHH=S1DHH+S2DHH+S3DHH+S4DHH
      FXCDGG=S2DGG+S3DGG+S4DGG
      FXCDHG=S1DHG+S2DHG+S3DHG+S4DHG
      PXC=(GAME*FXCDG-2.*THETA*FXCDH)/3.
      UXC=GAME*FXCDG-THETA*FXCDH
      SXC=(GAME*S2DG-S2+GAME*S3DG-S3+S4A*S4B*(GAME*S4CDG-S4C))-
     -  THETA*FXCDH
      if (dabs(SXC).lt.1.d-9*dabs(THETA*FXCDH)) SXC=0. ! accuracy loss
      CVXC=2.*THETA*(GAME*FXCDHG-FXCDH)-THETA**2*FXCDHH-GAME**2*FXCDGG
      if (dabs(CVXC).lt.1.d-9*dabs(GAME**2*FXCDGG)) CVXC=0. ! accuracy
      PDLH=THETA*(GAME*FXCDHG-2.*FXCDH-2.*THETA*FXCDHH)/3.
      PDLG=GAME*(FXCDG+GAME*FXCDGG-2.*THETA*FXCDHG)/3.
      PDRXC=PXC+(PDLG-2.*PDLH)/3.
      PDTXC=GAME*(THETA*FXCDHG-GAME*FXCDGG/3.)-
     -  THETA*(FXCDH/.75+THETA*FXCDHH/1.5)
      return
      end

* ======================  AUXILIARY SUBROUTINES   ==================== *
      subroutine FERINV7(F,N,X,XDF,XDFF) ! Inverse Fermi intergals
*                                                       Version 24.05.07
* X_q(f)=F^{-1}_q(f) : H.M.Antia 93 ApJS 84, 101
* q=N-1/2=-1/2,1/2,3/2,5/2 (N=0,1,2,3)
* Input: F - argument, N=q+1/2
* Output: X=X_q, XDF=dX/df, XDFF=d^2 X / df^2
* Relative error: N = 0     1      2      3
*        for X:    3.e-9, 4.2e-9, 2.3e-9, 6.2e-9
* jump at f=4:
*         for XDF: 6.e-7, 5.4e-7, 9.6e-8, 3.1e-7
*       for XDFF: 4.7e-5, 4.8e-5, 2.3e-6, 1.5e-6
      implicit double precision (A-H), double precision (O-Z)
      save
      dimension A(0:5,0:3),B(0:6,0:3),C(0:6,0:3),D(0:6,0:3),
     *  LA(0:3),LB(0:3),LD(0:3)
      data A/-1.570044577033d4,1.001958278442d4,-2.805343454951d3,
     *            4.121170498099d2,-3.174780572961d1,1.d0, ! X_{-1/2}
     *    1.999266880833d4,5.702479099336d3,6.610132843877d2,
     *            3.818838129486d1,1.d0,0., ! X_{1/2}
     *    1.715627994191d2,1.125926232897d2,2.056296753055d1,1.d0,0.,0.,
     *    2.138969250409d2,3.539903493971d1,1.d0,0.,0.,0./, ! X_{5/2}
     *  B/-2.782831558471d4,2.886114034012d4,-1.274243093149d4,
     *            3.063252215963d3,-4.225615045074d2,3.168918168284d1,
     *            -1.008561571363d0, ! X_{-1/2}
     *    1.771804140488d4,-2.014785161019d3,9.130355392717d1,
     *            -1.670718177489d0,0.,0.,0., ! X_{1/2}
     *    2.280653583157d2,1.193456203021d2,1.16774311354d1,
     *            -3.226808804038d-1,3.519268762788d-3,0.,0., ! X_{3/2}
     *    7.10854551271d2,9.873746988121d1,1.067755522895d0,
     *            -1.182798726503d-2,0.,0.,0./, ! X_{5/2}
     *  C/2.206779160034d-8,-1.437701234283d-6,6.103116850636d-5,
     *     -1.169411057416d-3,1.814141021608d-2,-9.588603457639d-2,1.d0,
     *  -1.277060388085d-2,7.187946804945d-2,-4.262314235106d-1,
     *      4.997559426872d-1,-1.285579118012d0,-3.930805454272d-1,1.d0,
     *  -6.321828169799d-3,-2.183147266896d-2,-1.05756279932d-1,
     *       -4.657944387545d-1,-5.951932864088d-1,3.6844711771d-1,1.d0,
     *  -3.312041011227d-2,1.315763372315d-1,-4.820942898296d-1,
     *       5.099038074944d-1,5.49561349863d-1,-1.498867562255d0,1.d0/,
     *  D/8.827116613576d-8,-5.750804196059d-6,2.429627688357d-4,
     *       -4.601959491394d-3,6.932122275919d-2,-3.217372489776d-1,
     *       3.124344749296d0, ! X_{-1/2}
     *    -9.745794806288d-3,5.485432756838d-2,-3.29946624326d-1,
     *       4.077841975923d-1,-1.145531476975d0,-6.067091689181d-2,0.,
     *    -4.381942605018d-3,-1.5132365041d-2,-7.850001283886d-2,
     *     -3.407561772612d-1,-5.074812565486d-1,-1.387107009074d-1,0.,
     *    -2.315515517515d-2,9.198776585252d-2,-3.835879295548d-1,
     *       5.415026856351d-1,-3.847241692193d-1,3.739781456585d-2,
     *       -3.008504449098d-2/, ! X_{5/2}
     *  LA/5,4,3,2/,LB/6,3,4,3/,LD/6,5,5,6/
      if (N.lt.0.or.N.gt.3) stop'FERINV7: Invalid subscript'
      if (F.le.0.) stop'FERINV7: Non-positive argument'
      if (F.lt.4.) then
         T=F
         UP=0.
         UP1=0.
         UP2=0.
         DOWN=0.
         DOWN1=0.
         DOWN2=0.
         do I=LA(N),0,-1
            UP=UP*T+A(I,N)
           if (I.ge.1) UP1=UP1*T+A(I,N)*I
           if (I.ge.2) UP2=UP2*T+A(I,N)*I*(I-1)
         enddo
         do I=LB(N),0,-1
            DOWN=DOWN*T+B(I,N)
           if (I.ge.1) DOWN1=DOWN1*T+B(I,N)*I
           if (I.ge.2) DOWN2=DOWN2*T+B(I,N)*I*(I-1)
         enddo
         X=dlog(T*UP/DOWN)
         XDF=1.d0/T+UP1/UP-DOWN1/DOWN
         XDFF=-1.d0/T**2+UP2/UP-(UP1/UP)**2-DOWN2/DOWN+(DOWN1/DOWN)**2
      else
         P=-1./(.5+N) ! = -1/(1+\nu) = power index
         T=F**P ! t - argument of the rational fraction
         T1=P*T/F ! dt/df
         T2=P*(P-1.)*T/F**2 ! d^2 t / df^2
         UP=0.
         UP1=0.
         UP2=0.
         DOWN=0.
         DOWN1=0.
         DOWN2=0.
         do I=6,0,-1
            UP=UP*T+C(I,N)
           if (I.ge.1) UP1=UP1*T+C(I,N)*I
           if (I.ge.2) UP2=UP2*T+C(I,N)*I*(I-1)
         enddo
         do I=LD(N),0,-1
            DOWN=DOWN*T+D(I,N)
           if (I.ge.1) DOWN1=DOWN1*T+D(I,N)*I
           if (I.ge.2) DOWN2=DOWN2*T+D(I,N)*I*(I-1)
         enddo
         R=UP/DOWN
         R1=(UP1-UP*DOWN1/DOWN)/DOWN ! dR/dt
         R2=(UP2-(2.*UP1*DOWN1+UP*DOWN2)/DOWN+2.*UP*(DOWN1/DOWN)**2)/
     /     DOWN
         X=R/T
         RT=(R1-R/T)/T
         XDF=T1*RT
         XDFF=T2*RT+T1**2*(R2-2.*RT)/T
      endif
      return
      end

      subroutine BLIN8(TEMP,CHI,
     *  W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT,
     *  W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT,
     *  W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT,
     *  W0XXX,W0XTT,W0XXT)
*                                                       Version 24.12.08
* Stems from BLIN7
* Differences - 3rd derivatives of W0 added; high-\chi expansion changed
* This is realization of the fit described in:
* REF.: G.Chabrier & A.Y.Potekhin, Phys.Rev.E 58, 4941 (1998)
* Input: TEMP=T/mc^2; CHI=(\mu-mc^2)/T
* Output: Wk - Fermi-Dirac integral of the order k+1/2
*         WkDX=dWk/dCHI, WkDT = dWk/dT, WkDXX=d^2 Wk / d CHI^2,
*         WkDTT=d^2 Wk / d T^2, WkDXT=d^2 Wk /dCHIdT,
*         W0XXX=d^3 W0 / d CHI^3, W0XTT=d^3 W0 /(d CHI d^2 T),
*         W0XXT=d^3 W0 /dCHI^2 dT
* Typical accuracy for W: a few times 1.e-4, maximum rel.error 0.002.
* Discontinuity at CHI=0.6: < 1.e-7 for W, ~ 1.e-6 for 2nd derivatives
* Discontinuity at CHI=14: a few times 1.e-5 for W, 0.005 for 2nd deriv.
      implicit double precision (A-H), double precision (O-Z)
      save
      dimension AC(5,0:2),AU(5,0:2),AX(5),AXI(5),AH(5),AV(5),AA(5,0:2)
      data AC/.37045057 d0, .41258437 d0,
     &        9.777982 d-2, 5.3734153 d-3, 3.8746281 d-5, ! c_i^0
     &        .39603109 d0, .69468795 d0, 
     &        .22322760 d0, 1.5262934 d-2, 1.3081939 d-4, ! c_i^1
     &        .76934619 d0, 1.7891437 d0, 
     &        .70754974 d0, 5.6755672 d-2, 5.5571480 d-4/ ! c_i^2
      data AU/.43139881 d0, 1.7597537 d0, 
     &        4.1044654 d0, 7.7467038 d0, 13.457678 d0, ! \chi_i^0
     &        .81763176 d0, 2.4723339 d0, 
     &        5.1160061 d0, 9.0441465 d0, 15.049882 d0, ! \chi_i^1
     &        1.2558461 d0, 3.2070406 d0, 
     &        6.1239082 d0, 10.316126 d0, 16.597079 d0/ ! \chi_i^2
      data AX/7.265351 d-2, .2694608 d0, 
     &        .533122 d0, .7868801 d0, .9569313 d0/ ! x_i
      data AXI/.26356032 d0, 1.4134031 d0, 
     &         3.5964258 d0, 7.0858100 d0, 12.640801 d0/ ! \xi_i
      data AH/3.818735 d-2, .1256732 d0, 
     &        .1986308 d0, .1976334 d0, .1065420 d0/ ! H_i
      data AV/.29505869 d0, .32064856 d0, 7.3915570 d-2, 
     &        3.6087389 d-3, 2.3369894 d-5/ ! \bar{V}_i
      parameter (PI=3.141592653d0,PI26=PI*PI/6.)
      dimension AM(0:2),AMDX(0:2),AMDT(0:2),
     *  AMDXX(0:2),AMDTT(0:2),AMDXT(0:2)
      data KRUN/0/
      if (TEMP.le.0.) stop'BLIN8: T < 0'
      if (CHI.lt.-80.) stop'BLIN8: CHI is too low'
      KRUN=KRUN+1
      if (KRUN.eq.1) then ! initialize
        do J=0,2
        do I=1,5
           AA(I,J)=dexp(-AU(I,J))
        enddo
        enddo
      endif
      if (CHI.lt.14..or.CHI*TEMP.lt..1) then
        do K=0,2
           W=0.
           WDX=0.
           WDT=0.
           WDXX=0.
           WDTT=0.
           WDXT=0.
           WDXXX=0.
           WDXTT=0.
           WDXXT=0.
          if (CHI.lt..6) then ! CHI < .6
             ECHI=dexp(-CHI)
            do I=1,5
               SQ=dsqrt(1.d0+AU(I,K)*TEMP/2.)
               DN=AA(I,K)+ECHI ! e^{-\chi_i}+e^{-\chi})
               W=W+AC(I,K)*SQ/DN
               WDX=WDX+AC(I,K)*SQ/DN**2
               WDT=WDT+AC(I,K)*AU(I,K)/(SQ*DN)
               WDXX=WDXX+AC(I,K)*SQ*(ECHI-AA(I,K))/DN**3
               WDTT=WDTT-AC(I,K)*AU(I,K)**2/(DN*SQ**3)
               WDXT=WDXT+AC(I,K)*AU(I,K)/(SQ*DN**2)
               WDXXX=WDXXX+AC(I,K)*SQ*
     *           (ECHI**2-4.*ECHI*AA(I,K)+AA(I,K)**2)/DN**4
               WDXTT=WDXTT-AC(I,K)*AU(I,K)**2/(DN**2*SQ**3)
               WDXXT=WDXXT+AC(I,K)*AU(I,K)*(ECHI-AA(I,K))/(SQ*DN**3)
            enddo
             WDX=WDX*ECHI
             WDT=WDT/4.
             WDXX=WDXX*ECHI
             WDTT=WDTT/16.
             WDXT=WDXT/4.*ECHI
             WDXXX=WDXXX*ECHI
             WDXTT=WDXTT*ECHI/16.
             WDXXT=WDXXT/4.*ECHI
          elseif (CHI.lt.14.) then ! .6 < CHI < 14
             SQCHI=dsqrt(CHI)
            do I=1,5
               CE=AX(I)-1.d0
               ECHI=dexp(CE*CHI)
               DE=1.d0+ECHI
               D=1.d0+AX(I)*CHI*TEMP/2.
               H=CHI**(K+1)*SQCHI*dsqrt(D)/DE
               HX=(K+1.5)/CHI+.25*AX(I)*TEMP/D-ECHI*CE/DE
               HDX=H*HX
               HXX=(K+1.5)/CHI**2+.125*(AX(I)*TEMP/D)**2+ECHI*(CE/DE)**2
               HDXX=HDX*HX-H*HXX
               HT=.25*AX(I)*CHI/D
               HDT=H*HT
               HDTT=-H*HT**2
               HTX=1./CHI-.5*AX(I)*TEMP/D
               HDXT=HDX*HT+HDT*HTX
               HDXXT=HDXX*HT+HDX*HT*HTX+HDXT*HTX+
     +           HDT*(.25*(AX(I)*TEMP/D)**2-1./CHI**2)
               HDXTT=HDXT*HT-HDX*.125*(AX(I)*CHI/D)**2+HDTT*HTX+
     +           HDT*.5*AX(I)*(TEMP*.5*AX(I)*CHI/D**2-1./D)
               HXXX=(2*K+3)/CHI**3+.125*(AX(I)*TEMP/D)**3-
     -           ECHI*(1.d0-ECHI)*(CE/DE)**3
               HDXXX=HDXX*HX-2.*HDX*HXX+H*HXXX
               XICHI=AXI(I)+CHI
               DXI=1.d0+XICHI*TEMP/2.
               V=XICHI**K*dsqrt(XICHI*DXI)
               VX=(K+.5)/XICHI+.25*TEMP/DXI
               VDX=V*VX
               VT=.25*XICHI/DXI
               VDT=V*VT
               VXX=(K+.5)/XICHI**2+.125*(TEMP/DXI)**2
               VDXX=VDX*VX-V*VXX
               VDXXX=VDXX*VX-2.*VDX*VXX+
     +           V*((2*K+1)/XICHI**3+.125*(TEMP/DXI)**3)
               VXXT=(1.-.5*TEMP*XICHI/DXI)/DXI
               VDTT=-V*VT**2
               VXT=1./XICHI-.5*TEMP/DXI
               VDXT=VDT*VXT+VDX*VT
               VDXXT=VDXT*VX+VDX*.25*VXXT-VDT*VXX-V*.25*TEMP/DXI*VXXT
               VDXTT=VDTT*VXT-VDT*.5*VXXT+VDXT*VT-
     -           VDX*.125*(XICHI/DXI)**2
               W=W+AH(I)*AX(I)**K*H+AV(I)*V
               WDX=WDX+AH(I)*AX(I)**K*HDX+AV(I)*VDX
               WDT=WDT+AH(I)*AX(I)**K*HDT+AV(I)*VDT
               WDXX=WDXX+AH(I)*AX(I)**K*HDXX+AV(I)*VDXX
               WDTT=WDTT+AH(I)*AX(I)**K*HDTT+AV(I)*VDTT
               WDXT=WDXT+AH(I)*AX(I)**K*HDXT+AV(I)*VDXT
               WDXXX=WDXXX+AH(I)*AX(I)**K*HDXXX+AV(I)*VDXXX
               WDXTT=WDXTT+AH(I)*AX(I)**K*HDXTT+AV(I)*VDXTT
               WDXXT=WDXXT+AH(I)*AX(I)**K*HDXXT+AV(I)*VDXXT
            enddo
          else ! CHI > 14, CHI*TEMP < 0.1: high-\chi,low-(T\chi) expans.
            do J=0,4 ! for nonrel.Fermi integrals from k+1/2 to k+4.5
               CNU=K+J+.5 ! nonrelativistic Fermi integral index \nu
               CHINU=CHI**(K+J)*dsqrt(CHI) ! \chi^\nu
               F=CHINU*(CHI/(CNU+1.)+PI26*CNU/CHI+ ! nonrel.Fermi
     +           .7*PI26**2*CNU*(CNU-1.)*(CNU-2.)/CHI**3)
               FDX=CHINU*(1.+PI26*CNU*(CNU-1.)/CHI**2+
     +           .7*PI26**2*CNU*(CNU-1.)*(CNU-2.)*(CNU-3.)/CHI**4)
               FDXX=CHINU/CHI*CNU*(1.+PI26*(CNU-1.)*(CNU-2.)/CHI**2+
     +           .7*PI26**2*(CNU-1.)*(CNU-2.)*(CNU-3.)*(CNU-4.)/CHI**4)
               FDXXX=CHINU/CHI**2*CNU*(CNU-1.)*
     *           (1.+PI26*(CNU-2.)*(CNU-3.)/CHI**2+
     +           .7*PI26**2*(CNU-2.)*(CNU-3.)*(CNU-4.)*(CNU-5.)/CHI**4)
              if (J.eq.0) then
                 W=F
                 WDX=FDX
                 WDXX=FDXX
                 WDXXX=FDXXX
              elseif (J.eq.1) then
                 C=.25*TEMP
                 W=W+C*F ! Fermi-Dirac, expressed through Fermi
                 WDX=WDX+C*FDX
                 WDXX=WDXX+C*FDXX
                 WDT=F/4.
                 WDXT=FDX/4.
                 WDTT=0.
                 WDXXX=WDXXX+C*FDXXX
                 WDXXT=FDXX/4.
                 WDXTT=0.
              else
                 C=-C/J*(2*J-3)/4.*TEMP
                 W=W+C*F
                 WDX=WDX+C*FDX
                 WDT=WDT+C*J/TEMP*F
                 WDXX=WDXX+C*FDXX
                 WDTT=WDTT+C*J*(J-1)/TEMP**2*F
                 WDXT=WDXT+C*J/TEMP*FDX
                 WDXXX=WDXXX+C*FDXXX
                 WDXTT=WDXTT+C*J*(J-1)/TEMP**2*FDX
                 WDXXT=WDXXT+C*J/TEMP*FDXX
              endif
            enddo ! next J
          endif
          if (K.eq.0) then
             W0=W
             W0DX=WDX
             W0DT=WDT
             W0DXX=WDXX
             W0DTT=WDTT
             W0DXT=WDXT
             W0XXX=WDXXX
             W0XTT=WDXTT
             W0XXT=WDXXT
          elseif (K.eq.1) then
             W1=W
             W1DX=WDX
             W1DT=WDT
             W1DXX=WDXX
             W1DTT=WDTT
             W1DXT=WDXT
          else
             W2=W
             W2DX=WDX
             W2DT=WDT
             W2DXX=WDXX
             W2DTT=WDTT
             W2DXT=WDXT
          endif
        enddo ! next K
*   ----------------------------------------------------------------   *
      else ! CHI > 14, CHI*TEMP > 0.1: general high-\chi expansion
         D=1.d0+CHI*TEMP/2.d0
         R=dsqrt(CHI*D)
         RX=.5d0/CHI+.25d0*TEMP/D
         RDX=R*RX
         RDT=.25d0*CHI**2/R
         RXX=-.5d0/CHI**2-.125d0*(TEMP/D)**2
         RDXX=RDX*RX+R*RXX
         RDTT=-.25d0*RDT*CHI/D
         RXT=.25d0/D-.125d0*CHI*TEMP/D**2
         RDXT=RDT*RX+R*RXT
         RXXX=1.d0/CHI**3+.125d0*(TEMP/D)**3
         RDXXX=RDXX*RX+2.d0*RDX*RXX+R*RXXX
         RXTT=-.25d0/D**2*CHI+.125d0*CHI**2*TEMP/D**3
         RDXTT=RDTT*RX+2.d0*RDT*RXT+R*RXTT
         RXXT=-RXT*TEMP/D
         RDXXT=RDXT*RX+RDX*RXT+RDT*RXX+R*RXXT
        do K=0,2
           DM=K+.5d0+(K+1.d0)*CHI*TEMP/2.d0
           AM(K)=CHI**K*DM/R
           FMX1=.5d0*(K+1.)*TEMP/DM
           FMX2=.25d0*TEMP/D
           FMX=(K-.5d0)/CHI+FMX1-FMX2
           AMDX(K)=AM(K)*FMX
           CKM=.5d0*(K+1.d0)/DM
           FMT1=CKM*CHI
           FMT2=.25d0*CHI/D
           FMT=FMT1-FMT2
           AMDT(K)=AM(K)*FMT
           FMXX=-(K-.5d0)/CHI**2-FMX1**2+2.d0*FMX2**2
           AMDXX(K)=AMDX(K)*FMX+AM(K)*FMXX
           FMTT=2.d0*FMT2**2-FMT1**2
           AMDTT(K)=AMDT(K)*FMT+AM(K)*FMTT
           AMDXT(K)=AMDX(K)*FMT+AM(K)*(CKM*(1.d0-CKM*CHI*TEMP)-
     -       .25d0/D+.125d0*CHI*TEMP/D**2)
          if (K.eq.0) then
             FMXXX=(2*K-1)/CHI**3+2.d0*FMX1**3-8.d0*FMX2**3
             AMDXXX=AMDXX(K)*FMX+2.d0*AMDX(K)*FMXX+AM(K)*FMXXX
             FMT1DX=CKM-TEMP*CHI*CKM**2
             FMT2DX=(.25d0-CHI*TEMP*.125d0/D)/D
             FMXT=FMT1DX-FMT2DX
             FMTTX=4.d0*FMT2*FMT2DX-2.d0*FMT1*FMT1DX
             AMDXTT=AMDXT(K)*FMT+AMDT(K)*FMXT+AMDX(K)*FMTT+AM(K)*FMTTX
             FMX1DT=CKM-CHI*TEMP*CKM**2
             FMX2DT=.25d0/D*(1.d0-.5d0*CHI*TEMP/D)
             FMXXT=4.d0*FMX2*FMX2DT-2.d0*FMX1*FMX1DT
             AMDXXT=AMDXT(K)*FMX+AMDX(K)*FMXT+AMDT(K)*FMXX+AM(K)*FMXXT
          endif
        enddo
         SQ2T=dsqrt(2.d0*TEMP)
           A=1.d0+CHI*TEMP+SQ2T*R
           ADX=TEMP+SQ2T*RDX
           ADT=CHI+R/SQ2T+SQ2T*RDT
           ADXX=SQ2T*RDXX
           ADTT=-R/SQ2T**3+2.d0/SQ2T*RDT+SQ2T*RDTT
           ADXT=1.d0+RDX/SQ2T+SQ2T*RDXT
           ADXTT=-RDX/SQ2T**3+2.d0/SQ2T*RDXT+SQ2T*RDXTT
           ADXXT=RDXX/SQ2T+SQ2T*RDXXT
           XT1=CHI+1.d0/TEMP
           Aln=dlog(A)
           FJ0=.5d0*XT1*R-Aln/SQ2T**3
           ASQ3=A*SQ2T**3
           ASQ3DX=ADX*SQ2T**3
           FJ0DX=.5d0*(R+XT1*RDX)-ADX/ASQ3
           FJ0DT=.5d0*(XT1*RDT-R/TEMP**2)-ADT/ASQ3+
     +       .75d0/(TEMP**2*SQ2T)*Aln
           FJ0DXX=RDX+.5d0*XT1*RDXX+(ADX/A)**2/SQ2T**3-ADXX/ASQ3
           FJ0DTT=R/TEMP**3-RDT/TEMP**2+.5d0*XT1*RDTT+
     +       3.d0/(ASQ3*TEMP)*ADT+
     +     (ADT/A)**2/SQ2T**3-ADTT/ASQ3-1.875d0/(TEMP**3*SQ2T)*Aln
           BXT=1.5d0/TEMP*ADX+ADX*ADT/A-ADXT
           BXXT=1.5d0/TEMP*ADXX+(ADXX*ADT+ADX*ADXT)/A-
     -       (ADX/A)**2*ADT-ADXXT
           FJ0DXT=.5d0*(RDT-RDX/TEMP**2+XT1*RDXT)+BXT/ASQ3
           FJ0XXX=RDXX*1.5d0+.5d0*XT1*RDXXX+
     +      (2.d0*ADX*(ADXX/A-(ADX/A)**2)-
     -      SQ2T*RDXXX+ADXX/ASQ3*ASQ3DX)/ASQ3
           FJ0XTT=RDX/TEMP**3-RDXT/TEMP**2+.5d0*(RDTT+XT1*RDXTT)+
     +      3.d0/TEMP*(ADXT-ADT/ASQ3*ASQ3DX)/ASQ3+
     +      (2.d0*ADT*(ADXT/A-ADT*ADX/A**2)-
     -      ADXTT+ADTT*ASQ3DX/ASQ3)/ASQ3-1.875d0/(TEMP**3*SQ2T)*ADX/A
           FJ0XXT=.5d0*(RDXT-RDXX/TEMP**2+RDXT+XT1*RDXXT)+
     +      (BXXT-BXT*ASQ3DX/ASQ3)/ASQ3
         W0=FJ0+PI26*AM(0)
         W0DX=FJ0DX+PI26*AMDX(0)
         W0DT=FJ0DT+PI26*AMDT(0)
         W0DXX=FJ0DXX+PI26*AMDXX(0)
         W0DTT=FJ0DTT+PI26*AMDTT(0)
         W0DXT=FJ0DXT+PI26*AMDXT(0)
         W0XXX=FJ0XXX+PI26*AMDXXX
         W0XTT=FJ0XTT+PI26*AMDXTT
         W0XXT=FJ0XXT+PI26*AMDXXT
           FJ1=(R**3/1.5d0-FJ0)/TEMP
           FJ1DX=(2.d0*R**2*RDX-FJ0DX)/TEMP
           FJ1DT=(2.d0*R**2*RDT-FJ0DT-FJ1)/TEMP
           FJ1DXX=(4.d0*R*RDX**2+2.d0*R**2*RDXX-FJ0DXX)/TEMP
           FJ1DTT=(4.d0*R*RDT**2+2.d0*R**2*RDTT-FJ0DTT-2.d0*FJ1DT)/TEMP
           FJ1DXT=(4.d0*R*RDX*RDT+2.d0*R**2*RDXT-FJ0DXT-FJ1DX)/TEMP
         W1=FJ1+PI26*AM(1)
         W1DX=FJ1DX+PI26*AMDX(1)
         W1DT=FJ1DT+PI26*AMDT(1)
         W1DXX=FJ1DXX+PI26*AMDXX(1)
         W1DTT=FJ1DTT+PI26*AMDTT(1)
         W1DXT=FJ1DXT+PI26*AMDXT(1)
           FJ2=(.5d0*CHI*R**3-1.25d0*FJ1)/TEMP
           FJ2DX=(.5d0*R**3+1.5d0*CHI*R**2*RDX-1.25d0*FJ1DX)/TEMP
           FJ2DT=(1.5d0*CHI*R**2*RDT-1.25d0*FJ1DT-FJ2)/TEMP
           FJ2DXX=(3.d0*R*RDX*(R+CHI*RDX)+1.5d0*CHI*R**2*RDXX-
     -       1.25d0*FJ1DXX)/TEMP
          FJ2DTT=(3.d0*CHI*R*(RDT**2+.5d0*R*RDTT)-
     -      1.25d0*FJ1DTT-2.d0*FJ2DT)/TEMP
           FJ2DXT=(1.5d0*R*RDT*(R+2.d0*CHI*RDX)+1.5d0*CHI*R**2*RDXT-
     -       1.25d0*FJ1DXT-FJ2DX)/TEMP
         W2=FJ2+PI26*AM(2)
         W2DX=FJ2DX+PI26*AMDX(2)
         W2DT=FJ2DT+PI26*AMDT(2)
         W2DXX=FJ2DXX+PI26*AMDXX(2)
         W2DTT=FJ2DTT+PI26*AMDTT(2)
         W2DXT=FJ2DXT+PI26*AMDXT(2)
      endif
      return
      end

      subroutine CHEMFIT(DENS,TEMP,CHI)
*                                                       Version 07.06.07
* This is merely an interface to CHEMFIT7 for compatibility purposes.
* Input:  DENS - electron density [a.u.=6.7483346e24 cm^{-3}],
* TEMP - temperature [a.u.=2Ryd=3.1577e5 K]
* Output: CHI=\mu/TEMP, where \mu - electron chem.pot.w/o rest-energy
      implicit double precision (A-H), double precision (O-Z)
      save
      DENR=DENS/2.5733806d6 ! n_e in rel.un.=\lambda_{Compton}^{-3}
      TEMR=TEMP/1.8778865d4 ! T in rel.un.=(mc^2/k)=5.93e9 K
      call CHEMFIT7(DENR,TEMR,CHI,CMU1,0,CMUDENR,CMUDT,CMUDTT)
      return
      end

      subroutine CHEMFIT7(DENR,TEMR,CHI,CMU1,KDERIV,
     *  CMUDENR,CMUDT,CMUDTT)
*                                                       Version 28.05.07
*                                                 "d0" inserted 23.03.09
* Fit to the chemical potential of free electron gas described in:
*     G.Chabrier & A.Y.Potekhin, Phys.Rev.E 58, 4941 (1998)
* Stems from CHEMFIT v.10.10.96. The main difference - derivatives.
*  All quantities are by default in relativistic units
* Input:  DENR - electron density, TEMR - temperature
*         KDERIV=0 if the derivatives are not required
* Output: CHI=CMU1/TEMR, where CMU1 = \mu-1 - chem.pot.w/o rest-energy
*         CMUDENR = (d\mu/d n_e)_T
*         CMUDT = (d\mu/dT)_V
*         CMUDTT = (d^2\mu/dT^2)_V
* CMUDENR,CMUDT, and CMUDTT =0 on output, if KREDIV=0
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter (PARA=1.612d0,PARB=6.192d0,PARC=.0944d0,
     *  PARF=5.535d0,PARG=.698d0)
      PF0=(29.6088132d0*DENR)**.33333333d0 ! Classical Fermi momentum
      if (PF0.gt.1.d-4) then
         TF=dsqrt(1.d0+PF0**2)-1.d0 ! Fermi temperature
      else
         TF=.5d0*PF0**2
      endif
      THETA=TEMR/TF
      THETA32=THETA*dsqrt(THETA)
      Q2=12.d0+8.d0/THETA32
      T1=dexp(-THETA) ! former ('96) 1/T
      U3=T1**2+PARA
      THETAC=THETA**PARC
      THETAG=THETA**PARG
      D3=PARB*THETAC*T1**2+PARF*THETAG
      Q3=1.365568127d0-U3/D3 ! 1.365...=2/\pi^{1/3}
      if (THETA.gt.1.d-5) then 
         Q1=1.5d0*T1/(1.d0-T1)
      else
         Q1=1.5d0/THETA
      endif
      SQT=dsqrt(TEMR)
      G=(1.d0+Q2*TEMR*Q3+Q1*SQT)*TEMR
      H=(1.d0+.5d0*TEMR/THETA)*(1.d0+Q2*TEMR)
      CT=1.d0+G/H
      F=.666666667d0/THETA32
      call FERINV7(F,1,X,XDF,XDFF)
      CHI=X ! non-relativistic result
     -  -    1.5d0*dlog(CT) ! Relativistic fit
      CMU1=TEMR*CHI ! Fit to chemical potential w/o mc^2
      if (KDERIV.eq.0) then ! DISMISS DERIVATIVES
         CMUDENR=0.
         CMUDT=0.
         CMUDTT=0.
         return
      endif
* CALCULATE DERIVATIVES:
* 1: derivatives of CHI over THETA and T
* (a): Non-relativistic result:
      THETA52=THETA32*THETA
      CHIDY=-XDF/THETA52 ! d\chi/d\theta
      CHIDYY=(XDFF/THETA**4-2.5d0*CHIDY)/THETA ! d^2\chi/d\theta^2
* (b): Relativistic corrections:
      if (THETA.gt.1.d-5) then 
         Q1D=-Q1/(1.d0-T1)
         Q1DD=-Q1D*(1.d0+T1)/(1.d0-T1)
      else
         Q1D=-1.5d0/THETA**2
         Q1DD=2.d0*Q1D/THETA
      endif
      Q2D=-12.d0/THETA52 ! d q_2 / d \theta
      Q2DD=30.d0/(THETA52*THETA) ! d^2 q_2 / d \theta^2
      U3D=-2.d0*T1**2
      D3D=PARF*PARG*THETAG/THETA+PARB*T1**2*THETAC*(PARC/THETA-2.d0)
      D3DD=PARF*PARG*(PARG-1.d0)*THETAG/THETA**2+
     +PARB*T1**2*THETAC*(PARC*(PARC-1.d0)/THETA**2-4.d0*PARC/THETA+4.d0)
      Q3D=(D3D*U3/D3-U3D)/D3
      Q3DD=(2.d0*U3D+(2.d0*U3D*D3D+U3*D3DD)/D3-2.d0*U3*(D3D/D3)**2)/D3
      GDY=TEMR*(Q1D*SQT+(Q2D*Q3+Q2*Q3D)*TEMR) ! dG/d\theta
      GDT=1.d0+1.5d0*Q1*SQT+2.d0*Q2*Q3*TEMR
      GDYY=TEMR*(Q1DD*SQT+(Q2DD*Q3+2.d0*Q2D*Q3D+Q2*Q3DD)*TEMR)
      GDTT=.75d0*Q1/SQT+2.d0*Q2*Q3
      GDYT=1.5d0*Q1D*SQT+2.d0*(Q2D*Q3+Q2*Q3D)*TEMR
      HDY=(-.5d0/THETA**2+Q2D+.5d0*(Q2D-Q2/THETA)/THETA*TEMR)*TEMR
      HDT=(.5d0+Q2*TEMR)/THETA+Q2
      HDYY=TEMR/THETA**3+Q2DD*TEMR+
     +  TEMR**2*(.5d0*Q2DD-Q2D/THETA+Q2/THETA**2)/THETA
      HDTT=Q2/THETA
      HDYT=Q2D*(1.d0+TEMR/THETA)-(.5d0+Q2*TEMR)/THETA**2
      CTY=GDY/G-HDY/H
      CTT=GDT/G-HDT/H
      GH=G/H
      CTDY=GH*CTY
      CTDT=GH*CTT
      CTDYY=CTDY*CTY+GH*(GDYY/G-(GDY/G)**2-HDYY/H+(HDY/H)**2)
      CTDTT=CTDT*CTT+GH*(GDTT/G-(GDT/G)**2-HDTT/H+(HDT/H)**2)
      CTDYT=CTDT*CTY+GH*(GDYT/G-GDY*GDT/G**2-HDYT/H+HDY*HDT/H**2)
      CHIDY=CHIDY-1.5d0*CTDY/CT
      CHIDT=-1.5d0*CTDT/CT
      CHIDYY=CHIDYY+1.5d0*((CTDY/CT)**2-CTDYY/CT)
      CHIDTT=1.5d0*((CTDT/CT)**2-CTDTT/CT)
      CHIDYT=1.5d0*(CTDY*CTDT/CT**2-CTDYT/CT)
      CMUDENR=-(THETA*PF0)**2/(3.d0*DENR*(1.d0+TF))*CHIDY
      CMUDT=CHI+THETA*CHIDY+TEMR*CHIDT
      CMUDTT=2.d0*(CHIDY/TF+CHIDT+THETA*CHIDYT)+
     +  THETA/TF*CHIDYY+TEMR*CHIDTT
      return
      end
