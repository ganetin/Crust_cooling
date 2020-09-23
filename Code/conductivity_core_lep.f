c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
      subroutine con_core_lep(Temp,kf_e,kf_m,
     1                    kf_p ,mst_p0 ,Tc_p ,        ! proton
     2                    kf_sm,mst_sm0,Tc_sm,        ! sigma-
     3                    kf_sp,mst_sp0,Tc_sp,        ! sigma+
     4                    sigma_lep,lambda_e,lambda_m,debug,
     5                    nu_e_s,nu_e_l)
c *********************************************************************
c Checked on Sept 4-7, 2009
c *********************************************************************
c Calculates the thermal conductivity in the core from
c Shternin & Yakovlev, Phys. Rev. D75, 103004, 2007
c for electrons and muons contributions.
c No electrical conductivity is calculated !
c Sigma- and Sigma+ scattering is NOT implemented !
c *********************************************************************
       implicit real*8(a-h,k-z)
       parameter(pi=3.1415926535d0)
       parameter(hbar=1.0546d-27,c=2.99792d10,kb=1.3806d-16)
       parameter(mp=1.6726d-24,me=9.1095d-28,e=4.803d-10)
c *********************************************************************
c Pairing correction factors: [t=T/Tc y=u(t)]
c 1S0:
      u_1s0(t)=dsqrt(1.d0-t)*(1.456d0-0.157d0/dsqrt(t)+1.764d0/t)
c *********************************************************************
       if (debug.eq.1.2d0) print *,'Entering con_core_lep:',
     1                   ' T, kfeo=',T,kf_e
c *** In case there are no leptons:
       if (kf_e.eq.0.d0) then
        lambda_e  =0.d0
        lambda_m  =0.d0
        lambda_lep=0.d0
        sigma_lep =0.d0
        if (debug.eq.1.2d0) print *,'Exiting con_core_lep:',
     1           ' sigma_lep, lambda_lep=',sigma_lep,lambda_lep
        return
       end if
c *** In case there are no muons:
      if (kf_m.eq.0.d0) then
       muons=0.d0
      else
       muons=1.d0
      end if
c ***
       T8=Temp/1.d8
c THERMAL CONDUCTIVITY: ***********************************************
c Trick to automatically eliminate absent baryons:
c That's because they show up only through the phase space integrals
c and this comes to zero if mst=0 !
       mst_p=mst_p0
       if (kf_p .eq.0.d0) mst_p =0.d0
       mst_sm=mst_sm0
       if (kf_sm.eq.0.d0) mst_sm=0.d0
       mst_sp=mst_sp0
       if (kf_sp.eq.0.d0) mst_sp=0.d0
c **** Define Fermi momenta ratios (massively used):
       rkf_0 =1.68d0/kf_e
       rkf_m =kf_m  /kf_e
       rkf_p =kf_p  /kf_e
c **** Screening momenta ratios:
c      (kf_e/q_l)**3:
       rkf_e_ql3=1./(0.00929d0*
     x          (1.d0+rkf_m+2.83d0*mst_p*rkf_0*rkf_p))**1.5
c      (kf_e/q_t)**2:
       rkf_e_qt2=1./(0.00929d0*
     x          (1.d0+rkf_m**2+rkf_p**2))
c **** Longitudinal collisional frequencies:
       nu_ee_par=1.43d11 *          rkf_0   *rkf_e_ql3 * T8**2
       nu_em_par=nu_ee_par *muons
       nu_ep_par=1.15d12 * mst_p**2*rkf_0**2*rkf_e_ql3 * T8**2
       if (muons.eq.1.d0) then
        nu_mm_par=nu_ee_par /rkf_m
        nu_me_par=nu_mm_par
        nu_mp_par=nu_ep_par /rkf_m
       else
        nu_mm_par=0.d0
        nu_me_par=0.d0
        nu_mp_par=0.d0
       end if
c **** Transverse ("perpendicular") collisional frequencies:
       nu_ee_per=6.49d14 * rkf_e_qt2 * T8
       nu_em_per=nu_ee_per * rkf_m**2
       nu_ep_per=nu_ee_per * rkf_p**2
       nu_mm_per=nu_ee_per * rkf_m**3
       nu_me_per=nu_ee_per * rkf_m
       nu_mp_per=nu_ee_per * rkf_m*rkf_p**2
c **** Cross ("prime") collisional frequencies:
       nu_ee_pri=4.38d12 * rkf_0**(2./3.)*rkf_e_qt2**(1./3.)*
     x                     rkf_e_ql3**(2./3.) * T8**(5./3.)
       if (muons.ne.0.d0) then
        nu_em_pri=nu_ee_pri *rkf_m**2
        nu_mm_pri=nu_ee_pri *rkf_m**3
        nu_me_pri=nu_ee_pri /rkf_m
       else
        nu_em_pri=0.d0
        nu_mm_pri=0.d0
        nu_me_pri=0.d0
       end if
c **** Effect of pairing:
c Proton 1S0 pairing:
       if ((Temp.le.Tc_p).and.(kf_p.gt.0.d0)) then
        y=u_1s0(Temp/Tc_p)
        r=(kf_e**2+kf_m**2)/kf_p**2
        R_l_pri  =(r+1.d0)**(1./3.) /
     x            ((r+1.d0)**2-0.757d0*y+(0.50651d0*y)**2)**(1./6.)
         p1=0.48d0-0.17d0*r
         p3=((1.d0-p1)*54.d0/4.d0/pi**2/r)**2
        R_tot_per=p1*dexp(-0.14d0*y**2)+(1.d0-p1)/dsqrt(1.d0+p3*y**2)
c       This is R_p_par fron the paper:
c
c        R_p_par  =(1.d0
c     x   + (26.33d0*y**2+0.376d0*y**4)*dexp(-dsqrt((3.675d0)**2+y**2))
c     x   + 0.742d0*(dexp((1.673d0)**2-dsqrt((1.673d0)**2+y**2))-1.d0)
c     y    ) *dexp((1.361d0)**2-dsqrt((1.361d0)**2+y**2))
c
c            which is obviously wrong !
c       Here is my guess for the correct value:
c
c        R_p_par  =(1.d0
c     x   + (26.33d0*y**2+0.376d0*y**4)*dexp(-dsqrt((3.675d0)**2+y**2))
c     x   + 0.742d0*(dexp((1.673d0)-dsqrt((1.673d0)**2+y**2))-1.d0)
c     y    ) *dexp((1.361d0)-dsqrt((1.361d0)**2+y**2))
c
c       Here is Peter Shternin's new version of the fit 
c       (it gives essentially the same result as my guess):
c
        R_p_par=(0.998d0 + 
     x           (2.04d0 + 0.68d0*dsqrt(y) + 5.7d0*y**2 + 1.71d0*y**4 ) 
     y           * dexp(-1.04d0*y) ) * dexp(-dsqrt(1.23d0+y**2))
       else
        R_tot_per=1.d0
        R_p_par  =1.d0
        R_l_pri  =1.d0
       end if
c **** Adjust collisional frequencies for pairing and add them
       nu_e_par=nu_ee_par+nu_em_par+nu_ep_par*R_p_par
       nu_e_per=(nu_ee_per+nu_em_per+nu_ep_per)*R_tot_per
       nu_e_pri=nu_ee_pri*R_l_pri
       nu_m_par=nu_mm_par+nu_me_par+nu_mp_par*R_p_par
       nu_m_per=(nu_mm_per+nu_me_per+nu_mp_per)*R_tot_per
       nu_m_pri=nu_mm_pri*R_l_pri
       nu_em_pri=nu_em_pri*R_l_pri
       nu_me_pri=nu_me_pri*R_l_pri
c **** Total collisional frecuencies:
       nu_e=nu_e_par+nu_e_per+nu_e_pri
       nu_m=nu_m_par+nu_m_per+nu_m_pri
c **** Relaxation times:
       if (muons.ne.0.d0) then
        tau_e = (nu_m-nu_em_pri) / (nu_e*nu_m-nu_em_pri*nu_me_pri)
        tau_m = (nu_e-nu_me_pri) / (nu_e*nu_m-nu_em_pri*nu_me_pri)
       else
        tau_e=1.d0/nu_e
        tau_m=0.d0
       end if
c **** thermal conductivies:
       lambda_e = 1.70d24 * T8 *(1.d15*tau_e) * (kf_e/1.68d0)**2
       lambda_m = 1.70d24 * T8 *(1.d15*tau_m) * (kf_m/1.68d0)**2
       lambda_lep=lambda_e+lambda_m

c ELECTRICAL CONDUCTIVITY: ********************************************
        sigma_lep=0.d0
c *********************************************************************
        nu_e_l=nu_e
        nu_e_s=0.d0
c ****
       if (debug.eq.1.2d0) print *,'Exiting con_core_lep:',
     1          ' sigma_lep, lambda_lep=',sigma_lep,lambda_lep
c ****
       return
      end
c *********************************************************************
c *********************************************************************

