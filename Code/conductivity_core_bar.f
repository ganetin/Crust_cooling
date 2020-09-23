c *********************************************************************
c *********************************************************************
      subroutine con_core_bar(Temp,kf_e,kf_mu,
     1                    kf_p ,mst_p ,Tc_p ,        ! proton
     2                    kf_n ,mst_n ,Tc_n ,isfn,   ! neutron
     3                    kf_la,mst_la,Tc_la,        ! lambda
     4                    kf_sm,mst_sm,Tc_sm,        ! sigma-
     5                    kf_s0,mst_s0,Tc_s0,        ! sigma0
     6                    kf_sp,mst_sp,Tc_sp,        ! sigma+
     7                    sigma_bar,lambda_bar,debug,
     8                    nu_e_s,nu_e_l,icontrol)
c *********************************************************************
c Checked on 2 Septembre 2009 against Fig 1, 3, & 4                   *
c *********************************************************************
c Calculates the neutron thermal conductivity in the core from        *
c Baiko, Haensel & Yakovlev, 2001, A&A 374, p 151                     *
c                                                                     *
c Chemical potentials mu's in MeV                                     *
c Fermi momenta kf's in fm^-1                                         *
c *********************************************************************
       implicit real*8(a-h,k-z)
       parameter(pi=3.1415926535d0)
       parameter(hbar=1.0546d-27,c=2.99792d10,kb=1.3806d-16)
       parameter(mu=1.6726d-24,me=9.1095d-28,e=4.803d-10)
c *********************************************************************
c Pairing correction factors: [t=T/Tc y=u(t)]
c 1S0:
       u_1s0(t)=dsqrt(1.d0-t)*(1.456d0-0.157d0/dsqrt(t)+1.764d0/t)
       u_3p2(t)=dsqrt(1.d0-t)*(0.7893d0+1.188d0/t)
c *********************************************************************
        if (debug.eq.1.2d0) print *,'Entering con_core_bar:',
     1                   ' T, kfe=',T,kf_e
c *********************************************************************
c      Neutrons contribution:
c ***
        Sn1=14.57d0/kf_n**1.5 *
     1      (1.d0-0.0788d0*kf_n+0.0883d0*kf_n**2) /
     2      (1.d0-0.1114d0*kf_n)
        Sn2=7.880d0/kf_n**2 *
     1      (1.d0-0.2241d0*kf_n+0.2006d0*kf_n**2) /
     2      (1.d0-0.1742d0*kf_n)
        Sp1=0.8007d0*kf_p/kf_n**2 *
     1      (1.d0+31.28d0*kf_p-0.0004285d0*kf_p**2+
     1       26.85d0*kf_n+0.08012d0*kf_n**2) /      
     2      (1.d0-0.5898d0*kf_n+0.2368d0*kf_n**2+
     2       0.5838d0*kf_p**2+0.884d0*kf_n*kf_p)      
        Sp2=0.3830d0*kf_p**4/kf_n**5.5 *
     1      (1.d0+102.d0*kf_p+53.91d0*kf_n) /
     2      (1.d0-0.7087d0*kf_n+0.2537d0*kf_n**2+
     2       9.404*kf_p**2-1.589d0*kf_n*kf_p)
c ***
        u=kf_n-1.665d0
        Kn1=(0.4583d0+0.892d0*u**2-0.5497d0*u**3-0.06205d0*kf_p
     1      +0.04022d0*kf_p**2+0.2122d0*u*kf_p)
     2      / mst_n**2
        u=kf_n-1.556d0
        Kn2=(0.4891d0+1.111d0*u**2-0.2283d0*u**3+0.01589d0*kf_p
     1      -0.02099*kf_p**2+0.2773*u*kf_p)
     2      / mst_n**2
        u=kf_n-2.126d0
        Kp1=(0.04377d0+1.100d0*u**2+0.1180d0*u**3+0.1626d0*kf_p 
     1       +0.3871d0*u*kf_p-0.2990d0*u**4)
     2      / mst_p**2
        u=kf_n-2.116d0
        Kp2=(0.0001313d0+1.248d0*u**2+0.2403d0*u**3+0.3257d0*kf_p
     1       +0.5536d0*u*kf_p-0.3237d0*u**4+0.09786d0*u**2*kf_p)
     2      / mst_p**2
        if (icontrol.ge.2) then
         Kn1=1.d0
         Kn2=1.d0
         Kp1=1.d0
         Kp2=1.d0
        end if
c *** Pairing effects:
        if (Temp.le.Tc_p) then
         tau=Temp/Tc_p
         yp=u_1s0(tau)
        else
         yp=0.d0
        end if
        if (Temp.le.Tc_n) then
         tau=Temp/Tc_n
         if (isfn.eq.1) then
          yn=u_1s0(tau)
         else if (isfn.eq.3) then
          yn=u_3p2(tau)
         else
          pause 'Subroutine con_core_bar: isfn badly defined !'
         end if
        else
         yn=0.d0
        end if
        call con_core_bar_pairing_supr(yn,yp,Rn1,Rn2,Rp1,Rp2,RC)
c ***
        Snn=Sn2*Kn2*Rn2 + 3.0d0*Sn1*Kn1*(1.d0*Rn1-Rn2)
        Snp=Sp2*Kp2*Rp2 + 0.5d0*Sp1*Kp1*(3.d0*Rp1-Rp2)
c       When T<<Tc_p and/or Tc_n Snn and Snp may vanish, so better use:
        Snn = max(Snn,1.d-200)
        Snp = max(Snp,1.d-200)
        if (icontrol.eq.3) Snp=0.d0
        nu_nn=3.48d15*  mst_n**3    *(Temp/1.d8)**2 *Snn
        nu_np=3.48d15*mst_n*mst_p**2*(Temp/1.d8)**2 *Snp
        tau_n=RC/(nu_nn+nu_np)
c ***
        lambda_n=7.2d23* (Temp/1.d8) *
     x           RC**2/mst_n* 1.d15/(nu_nn+nu_np) * (kf_n/1.68)**3
c *********************************************************************
c      Lambda contribution:
       lambda_la=0.d0
c *********************************************************************
c      Sigma0 contribution:
       lambda_s0=0.d0
c *********************************************************************
       lambda_bar=lambda_n+lambda_la+lambda_s0
c *********************************************************************
c ELECTRICAL CONDUCTIVITY:
        sigma_bar=0.d0
c *********************************************************************
       if (debug.eq.1.2d0) print *,'Exiting con_core_bar:',
     1            ' sigma_bar, lambda_bar=',sigma_bar,lambda_bar
c ****
        nu_e_l=0.d0
        nu_e_s=0.d0
c ****
       return
      end
c *********************************************************************
c *********************************************************************
      subroutine con_core_bar_pairing_supr(yn,yp,Rn1,Rn2,Rp1,Rp2,RC)
       implicit real*8(a-h,k-z)
        if (yn.eq.0.d0) then
         Rn1=1.d0
         Rn2=1.d0
         RC =1.d0
        else
         Rn1=
     1       (2.d0/3.d0) * 
     1       (0.9468d0+dsqrt((0.0532d0)**2+0.5346d0*yn**2))**3 *
     1       dexp(0.377d0-dsqrt((0.377d0)**2+4.d0*yn**2)) +
     2       (1.d0/3.d0) *
     2       (1.d0+1.351d0*yn**2)**2 *
     2       dexp(0.169d0-dsqrt((0.169d0)**2+9.d0*yn**2))
         Rn2=
     1       0.5d0 *
     1       (0.6242d0+dsqrt((0.3758d0)**2+0.07198d0*yn**2))**3 *
     1       dexp(3.6724d0-dsqrt((3.6724d0)**2+4.d0*yn**2)) +
     2       0.5d0 *
     2       (1.d0+0.01211d0*yn**2)**9 *
     2       dexp(7.5351d0-dsqrt((7.5351d0)**2+9.d0*yn**2))
         RC=(0.647d0+dsqrt((0.353d0)**2+0.109d0*yn**2))**1.5 *
     x      dexp(1.39d0-dsqrt((1.39d0)**2+yn**2))
        end if
        if ((yn.eq.0.d0).and.(yp.eq.0.d0)) then
         Rp1=1.d0
         Rp2=1.d0
        else if ((yn.gt.0.d0).and.(yp.eq.0.d0)) then
         Rp1=(0.4459d0+dsqrt((0.5541d0)**2+0.03016d0*yn**2))**2 *
     1       dexp(2.1178d0-dsqrt((2.1178d0)**2+yn**2))
         Rp2=(0.801d0+dsqrt((0.199d0)**2+0.04645d0*yn**2))**2 *
     1       dexp(2.3569d0-dsqrt((2.3569d0)**2+yn**2))
        else if ((yn.eq.0.d0).and.(yp.gt.0.d0)) then
         Rp1=
     1       0.5d0 *
     1       (0.3695d0+dsqrt((0.6305)**2+0.01064d0*yp**2)) *
     1       dexp(2.4451d0-dsqrt((2.4451d0)**2+yp**2)) +
     2       0.5d0 *       
     2       (1.d0+0.1917*yp**2)**1.4 *
     2       dexp(4.6627d0-dsqrt((4.6627d0)**2+4.d0*yp**2))
         Rp2=
     1       0.0436d0 *       
     1       (dsqrt((3.345d0)**2+19.55*yp**2)-3.345d0) *
     1       dexp(2.0247d0-dsqrt((2.0247d0)**2+yp**2))
     2       + 0.0654d0 * dexp(8.992d0-dsqrt((8.992d0)**2+1.5d0*yp**2))
     3       + 0.8910d0 * dexp(9.627d0-dsqrt((9.627d0)**2+9.0d0*yp**2))
        else
         y_p=max(yn,yp)
         y_m=min(yn,yp)
         u_p=dsqrt(y_p**2+(1.485d0)**2)-1.485d0
         u_m=dsqrt(y_m**2+(1.485d0)**2)-1.485d0
         up =dsqrt(yp**2 +(1.485d0)**2)-1.485d0
         un =dsqrt(yn**2 +(1.485d0)**2)-1.485d0
         Rp1=
     1       dexp(-u_p-u_m) * (0.7751d0+0.4823d0*un+0.1124d0*up+
     1       0.04991d0*un**2+0.08513d0*un*up+0.01284d0*un**2*up)       
     2       + dexp(-2.d0*u_p) * (0.2249d0+0.3539d0*u_p-0.2189d0*u_m-
     2       0.6069d0*un*u_m+0.7362d0*up*u_p)       
         u_p=dsqrt(y_p**2+(1.761d0)**2)-1.761d0
         u_m=dsqrt(y_m**2+(1.761d0)**2)-1.761d0
         up =dsqrt(yp**2 +(1.761d0)**2)-1.761d0
         un =dsqrt(yn**2 +(1.761d0)**2)-1.761d0
         Rp2=
     1       dexp(-u_p-u_m)*(1.1032d0+0.8645d0*un+0.2042d0*up+
     1       0.07937d0*un**2+0.1451d0*un*up+0.01333d0*un**2*up)       
     2       +dexp(-2.d0*u_p)*(-0.1032d0-0.2340d0*u_p+0.06152d0*un*u_p+
     2       0.7533d0*un*u_m-1.007d0*up*u_p)       
        end if
       return
      end
c *********************************************************************
c *********************************************************************
