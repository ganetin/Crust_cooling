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

c ***************************

      fexp(x)=dexp(max(x,-7.d2))


      u_1s0(t)=dsqrt(1.-t)*(1.456-0.157/dsqrt(t)+1.764/t)
      r_1s0(u)=(0.4186+dsqrt(1.007**2+(0.5010*u)**2))**2.5*
     1         dexp(1.456-dsqrt(1.456**2+u**2))

c      f_1(t,a0i,a1i,a3i,d0i)=(1.d0+dexp(-a0i*a1i*d0i/a3i))/
c     1        (1.d0+dexp((t/(MeV2erg/kboltz)-a0i*a1i*d0i)/a3i))
c      f_2(t,a0i,a2i,a3i,d0i)=(1.d0+dexp(-a0i*a2i*d0i/a3i))/
c     1        (1.d0+dexp((t/(MeV2erg/kboltz)-a0i*a2i*d0i)/a3i))
c      h(u,t,a0i,a1i,a2i,a3i,d0i)=
c     1           r_1s0(u)*f_1(t,a0i,a1i,a3i,d0i)
c     2         +( 1.d0-f_2(t,a0i,a2i,a3i,d0i) )

      u_3p2(t)=dsqrt(1.-t)*(5.596+8.424/t)
      r_3p2(u)=(0.6893+dsqrt(0.790**2+(0.03983*u)**2))**2*
     1         fexp(1.934-dsqrt(1.934**2+u**2/(16.*pi)))

c      write(90,*)'specheat_hfb'
c ***************************
c ****** get Cv-ions :
     
      if (rho .lt. rhocore) then
       call cvion(t,rho,acell,aion,zion,cv_ions)
       cv_nrot=0.d0
       cv_nvib=0.d0
       cvions=cv_ions+cv_nrot+cv_nvib
      else
       cvions=0.
      end if

c ****** get Cv-electrons :

      if ((rho .lt. rhodrip).and.(icvel_nodeg.eq.1)) then
       call cvelec(t,rho,acell,aion,zion,i,cvelectron)
      else
       cvelectron=cve(i)*t
      end if

c ***** get Cv-muons :

      cvmuon=cvm(i)*t

c ****** get Cv-neutrons :

c HERE DANY cccccccccccccccccccccccc
c       raise=1.1d0
       raise=1.0d0
cccccccccccccccccccccccccccccccccccc
      kboltz=1.380d-16
      MeV2erg=1.602177d-6

      if ((rho.ge.rhodrip).and.(rho.le.rhocore).and.(istrange.eq.0)) 
     1 then
            raise=1.0d0
            if (t .gt. 20*tcn(i)) then
               r=0.d0  
               t0=1.d-10
            else if (t .gt. tcn(i)) then
                u=0.001
                r=r_1s0(u)
                t0=1.d-10
            else
                t0=t/tcn(i)
                if (i.le.isf) then
                     u=u_3p2(t0)
                     r=r_3p2(u)
                 else
                     u=u_1s0(t0)
                     r=r_1s0(u)
                 endif
            endif
cc Include the new specific heat
cc density of free neutrons (at t=0)
            nn0=kfn(i)**3/(3.*pi**2)
	        if (acell.eq.0.d0) then
	          nnm=bar(i)
	        else
		      nnm=(acell-zion)/acell*bar(i)
		    endif
            if ((kbol*t).gt.5.5d0) then
              nn=nnm
            else
              nn=nn0+kbol*t/5.5*(nnm-nn0)	
            endif
c            write(*,*)i,nn0,nnm,acell,zion,bar(i)
c            stop
            nn=nn0
cc calculate the FErmi energy for the free neutrons
            mn=939.56d0*mstn(i)
            kfnn=3.*pi**2*nn**0.33333
c            enfern=0.5*hbc**2/mn*kfn(i)**2
            enfern=0.5*hbc**2/mn*kfnn**2
cc temperature effect in cvn (beyond the linear term)
            vartemp=pi*t*kbol/enfern
            fcvn=1.-7./40.*vartemp**2-155./896.*vartemp**4
cc classial neutron specific heat in erg.cm-3.K-1
            cvnclass = 1.5*nn*kbol*MeV*1.d39
cc smooth transition from classical gas to quantal gas            
             if (vartemp.lt.0.5d0) then
                  xop=1.
             else if (vartemp.gt.1.5d0) then
                  xop=0.d0
             else
          	   xop=1./(1.+dexp(+5.*(vartemp-1.)))
             endif
c            write(*,*)i,mn,kfn(i),enfern,vartemp,fcvn*(1.-xop)
cc
c          write(*,*)'transition quantique classique prise en compte'
            if ((sfn1s0.eq.301.).or.(sfn1s0.eq.302.)) then
cc modification of the function R
               hfb_f1=(1.d0+dexp(-hfb_a0(i)*hfb_a1(i)*hfb_d0(i)
     1                /hfb_a3(i)))/(1.d0+dexp((t/(MeV2erg/kboltz)
     2                -hfb_a0(i)*hfb_a1(i)*hfb_d0(i))/hfb_a3(i)))
               hfb_f2=(1.d0+dexp(-hfb_a0(i)*hfb_a2(i)*hfb_d0(i)/
     1                hfb_a3(i)))/(1.d0+dexp((t/(MeV2erg/kboltz)
     2                -hfb_a0(i)*hfb_a2(i)*hfb_d0(i))/hfb_a3(i)))
               hfb_r=r*hfb_f1+(1-hfb_f2)
c		     write(*,'(i4,2d14.4)')i,cvn(i),cvnclass
c		     stop
cc calculate the specific heat
               cvneutron =cvn(i)*hfb_r*t*fcvn*xop+cvnclass*(1.-xop)
c		     write(*,'(i4,3d14.4)')i,cvn(i)*t,cvnclass,cvneutron
             else if (sfn1s0.eq.300.) then
               cvneutron =cvn(i)*t*fcvn*xop+cvnclass*(1.-xop)
             else
               cvneutron =cvn(i)*t*r             
             endif
cc end of new specific heat
      else if ((rho.lt.rhodrip).and.(istrange.eq.0)) then
             cvneutron=0.d0     
      else if ((rho.gt.rhocore).and.(istrange.eq.0)) then
            if (t .lt. tcn(i))then
             t0=min(0.999999999999d0,t/tcn(i))
                  if (i.le.isf) then
                   u=u_3p2(t0)
                   r=r_3p2(u)
                  else
                   u=u_1s0(t0)
                   r=r_1s0(u)
                  end if
            else
             r=1.d0
            end if
            cvneutron =cvn(i)*t*r
        end if

c**********************************
c***** No Superfluidity************
c            cvneutron =cvn(i)*t
c**********************************
c**********************************


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
       call cvion(t,rho,aion,zion,cv_ions)
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
c *********************************************************************
c *********************************************************************
c      subroutine cvion (t,rho,a,z,cv)
c     
c *** checked on oct. 22 1990 *******
c
c      implicit real*8(a-h,k-z)
c      parameter (rhodrip=4.3e11)
c     
c      dimension cv0(0:14)
c      save bcv,ccv,dcv,hcv,cte,cv0
c     
c      data bcv,ccv,dcv,hcv,cte/.95043,.18956,-.81487,3225,1.417e2/
c      data cv0/0.0,2.956,2.829,2.633,2.389,2.118,1.840,1.572,
c     1             1.323,1.102,0.909,0.745,0.609,0.496,0.404/
c
c     
c      gamma=2.273e5*z**2*(rho/a)**(1./3.)/t
c      if (rho .ge. rhodrip) then
c       a1=3.*z
c      else
c       a1=a
c      endif
c      nionkb=1.38e-16*6.022e23*rho/a
c      delta=1./t*z*dsqrt(rho/(a1*a))*6.022e23
c      if (gamma .le. .1) then
c       cv=1.5*nionkb
c       return
c      else if (gamma .le. .2) then
c       cv1=1.5*nionkb
c       cv2=nionkb*(.75*bcv*gamma**.25+1.25*ccv/gamma**.25+dcv+1.5)
c       cv=(gamma-.1)/.1*cv2+(.2-gamma)/.1*cv1
c       return
c      else if (gamma .le. 178.) then
c       cv=nionkb*(.75*bcv*gamma**.25+1.25*ccv/gamma**.25+dcv+1.5)
c       return
c      else if ((gamma.le.210.).and.(delta.ge.1.e19)) then
c       cv1=nionkb*(1.5+3.*hcv/gamma**2+1.5)
c       cv0(0)=1.5+3.*hcv/gamma**2+1.5
c       i1=int(delta*2.e-20)
c       cv2=nionkb*(cv0(i1)+(delta*2.e-20-i1)*(cv0(i1+1)-cv0(i1)))
c       cv=(gamma-178.)/32.*cv2+(210.-gamma)/32.*cv1
c       return
c      else if (delta .le. 1e19) then
c       cv=nionkb*(1.5+3.*hcv/gamma**2+1.5)
c       return
c      else if ((delta .gt. 1e19) .and. (delta .lt. 7e20)) then
c       cv0(0)=1.5+3.*hcv/gamma**2+1.5
c       i1=int(delta*2.e-20)
c       cv=nionkb*(cv0(i1)+(delta*2.e-20-i1)*(cv0(i1+1)-cv0(i1)))
c       return
c      else
c       delta1=delta*1.e-20
c       cv=nionkb*cte/delta1**3
c       return
c      end if
c     
c      return
c     
c      end
c
c ************************************************************************
c ************************************************************************

c *********************************************************************
c *********************************************************************
      subroutine cvion (T,rho,A,A1,Z,cv)
c      A=A_cell  A1=A_ion
       implicit real*8(a-h,k-z)
       common /control_cvion/icvion
c HERE DANY ccccccccccccccccccc
        if (icvion.eq.1) then
         ! Original version
         call cvion_1(T,rho,A,Z,cv)
            else if (icvion.eq.2) then
         ! Modified original version
         call cvion_2(T,rho,A,A1,Z,cv)
        else if ((icvion.eq.3).or.(icvion.eq.4)) then
         call cvion_3(T,rho,A,A1,Z,cv)
        else
         pause 'Subroutine cvion: icvion undefined !'
        end if
       return
      end
cc *********************************************************************
cc *********************************************************************
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
c         From Slattery et al.: actually gives the same result as cv_liquid
c         at gamma >1 !
c         cv1=.75*bcv*gamma**.25+1.25*ccv/gamma**.25+dcv+1.5
         cv=cv1
        else if (Gamma.gt.gammaliq) then
         call cv_harmonic_crystal (T,rho,A,A1,Z,cv2h)
         call cv_anharmonic_crystal (T,rho,A,A1,Z,cv2ah)
         cv=cv2h+cv2ah
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
cc *********************************************************************
cc *********************************************************************
cc *********************************************************************
cc *********************************************************************
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
c        cv_q=-T_p**2/T**2/12.d0
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
     
