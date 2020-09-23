      subroutine neutrino(i,t,rho,a,z,qtot,
     1   qeebrem,qnpb,qplasma,qsynch,qbubble,qpair,qphoto,qbrem_nn,
     2   qmurca_nucl,qbrem_nucl,qmurca_hyp,qbrem_hyp,
     3   qdurca_np,qdurca_lap,qdurca_smn,qdurca_smla,qdurca_sms0,
     4   qfast,
     5   qdurca_q,qmurca_q,
     6   qpbf_n1s0,qpbf_n3p2,qpbf_p1s0,qpbf_q,qdurca_crust,m_dot,
     7   debug,naa)

c **** checked  august 21 1991 **************

      implicit real*8 (a-h,k-z)
      INCLUDE 'size.inc.f'
      INCLUDE 'rho_limits.inc.f'
      INCLUDE 'control_nu.inc.f'
      INCLUDE 'pairing.inc.f'
      INCLUDE 'fermi.inc.f'
      INCLUDE 'profile_comp.inc.f'
      INCLUDE 'profile_star.inc.f'
      INCLUDE 'mag_field.inc.f'
      common/gravity/gs14,compactness,xrho,radius,rcore,dvcore

      if (debug.ge.2.) print *,'Entering subroutine `neutrino'' ',
     2       ' T, rho, A, Z = ',t,rho,a,z
c *** ELECTRON-ELECTRON PAIR BREMSSTRAHLUNG:
      if ((rho.lt.rhocore).and.(inu_crust==1)) then
         mu_el=kfe(i)*197.
         call neebrem(T,mu_el,qeebrem)
      else
         qeebrem=0.0d0
      end if
c *** ELECTRON-ION PAIR BREMSSTRAHLUNG:
      if (inu_eion.eq.1) then
         if ((rho.lt.rhocore).and.(inu_crust==1)) then
            call npb_new(t,rho,qnpb)
         else
            qnpb=0.0d0
         end if
      else if (inu_eion.eq.2) then
         if ((rho.lt.rhocore).and.(inu_crust==1)) then
            call npb(t,rho,a,z,qnpb)
         else
            qnpb=0.0d0
         end if
      else
         if ((rho.lt.rhocore).and.(inu_crust==1)) then
            qnpb=0.0d0
            if (print_it.ne.1.) then
             print *,'No npb: Rho, Qnpb=',rho,qnpb
             print_it=1.
            end if
         else
            qnpb=0.0d0
         end if
      end if
c *** PLASMA NEUTRINO:
      if (inu_plasma.eq.1) then
         if ((rho.lt.rhocore).and.(inu_crust==1)) then
            call nplasma(t,rho,a,z,qplasma)
         else
            qplasma=0.0d0
         end if
      else if (inu_plasma.eq.-1) then
         if ((rho.lt.rhocore).and.(inu_crust==1)) then
            call nplasma_old(t,rho,a,z,qplasma)
         else
            qplasma=0.0d0
         end if
      else
         qplasma=0.0d0
      end if
c *** SYNCHROTRON NEUTRINO:
      if (inu_synch.eq.1) then
         if ((rho.lt.rhocore).and.(inu_crust==1)) then
            call nsynch(t,nbfield2(i),kfe(i),qsynch)
         end if
      else
         qsynch=0.0d0
      end if
c *** BUBBLE NEUTRINO:
      if (inu_bubble.eq.1) then
         if ((rho.lt.rhocore).and.(inu_crust==1)) then
            call nbub(i,t,rho,a,z,qbubble)
         else
            qbubble=0.0d0
         end if
      else
         qbubble=0.0d0
      end if
c *** NEUTRINO PAIR:
      if (inu_pair.eq.1) then
         if ((rho.lt.rhocore).and.(inu_crust==1)) then
            call npair(t,rho,a,z,qpair)
         else
            qpair=0.0d0
         end if
      else
         qpair=0.0d0
      end if
c **** PHOTO-NEUTRINO:
      if (inu_photo.eq.1) then
         if ((rho.lt.rhocore).and.(inu_crust==1)) then
            call nphoto(t,rho,a,z,qphoto)
         else
            qphoto=0.0d0
         end if
      else
          qphoto=0.0d0
      end if
c **** DURCA_CRUST:
      if (inu_durca_crust.eq.1) then
         if (m_dot.gt.0.) then
         if ((rho.lt.rhocore).and.(inu_crust==1)) then
            call durca_crust(i,t,qdurca_crust)
         else
            qdurca_crust=0.0d0
         end if
      else
          qdurca_crust=0.0d0
      end if
      end if
c *** NN-BREMSTRAHLUNG in the inner crust:
      if ((rho.lt.rhocore).and.(rho.ge.rhodrip).and.(inu_crust==1)) then
         call nubrem_crust_nn(i,t,v_ion(i),qbrem_nn)
      else
         qbrem_nn=0.d0
      end if
c *** URCA et al. PROCESSES:
      if (rho.ge.rhocore) then
         if (istrange.eq.0) then
            call numurca_nucl(i,t,qmurca_nucl)
            qmurca_nucl=qmurca_nucl*(1.d0+murca_increase)
            qmurca_nucl=qmurca_nucl*fhad(i)
            call nubrem_nucl(i,t,qbrem_nucl)
            qbrem_nucl=qbrem_nucl*(1.d0+murca_increase)
            qbrem_nucl=qbrem_nucl*fhad(i)
            call numurca_hyp(i,t,qmurca_hyp)
            qmurca_hyp=qmurca_hyp*fhad(i)
            call nubrem_hyp(i,t,qbrem_hyp)
            qbrem_hyp=qbrem_hyp*fhad(i)
            if (inu_durca.eq.1) then
             call nudurca_h(i,t,rho,qdurca_np,qdurca_lap,
     1                        qdurca_smn,qdurca_smla,qdurca_sms0)
             qdurca_np=qdurca_np*fhad(i)
             qdurca_lap=qdurca_lap*fhad(i)
             qdurca_smn=qdurca_smn*fhad(i)
             qdurca_smla=qdurca_smla*fhad(i)
             qdurca_sms0=qdurca_sms0*fhad(i)
            else
             qdurca_np=0.0d0
             qdurca_lap=0.0d0
             qdurca_smn=0.0d0
             qdurca_smla=0.0d0
             qdurca_sms0=0.0d0
            end if
c *** FAST neutrino emission:
            call nufast (i,t,rho,qfast)
            qfast=qfast*fhad(i)
c *** QUARK processes:
            call nudurca_q(i,t,rho,qdurca_q)
            call numurca_q(i,t,rho,qmurca_q)
            qdurca_q=qdurca_q*(1.d0-fhad(i))
            qmurca_q=qmurca_q*(1.d0-fhad(i))
            qstrange=0.d0
c *** STRANGE QUARK MATTER processes:
         else if (istrange.eq.1) then
            if (c_nu_str.le.1.d0) then
               call nu_strange(i,T,qstrange,debug)
            else
               qstrange=c_nu_str*(T/1.d9)**p_nu_str
            end if
            qmurca_nucl=0.0d0
            qbrem_nucl=0.0d0
            qmurca_hyp=0.0d0
            qbrem_hyp=0.0d0
            qdurca_np=0.0d0
            qdurca_lap=0.0d0
            qdurca_smn=0.0d0
            qdurca_smla=0.0d0
            qdurca_sms0=0.0d0
            qfast=0.0d0
            qdurca_q=0.0d0
            qmurca_q=0.0d0
         else
            print *,'neutrino: istrange not defined !'
            pause
         end if
      else
         qmurca_nucl=0.0d0
         qbrem_nucl=0.0d0
         qmurca_hyp=0.0d0
         qbrem_hyp=0.0d0
         qdurca_np=0.0d0
         qdurca_lap=0.0d0
         qdurca_smn=0.0d0
         qdurca_smla=0.0d0
         qdurca_sms0=0.0d0
         qfast=0.0d0
         qdurca_q=0.0d0
         qmurca_q=0.0d0
         qstrange=0.d0
      end if
c *** PBF PROCESSES:
      if (istrange.eq.0) then
c      Neutrons 1S0:
       if ((inu_n1s0_pbf.eq.1).and.(i.gt.isf)) then
          call nu_1s0_pbf(t,tcn(i),mstn(i),kfn(i),qpbf_n1s0)
          qpbf_n1s0=qpbf_n1s0*fhad(i)
       else
          qpbf_n1s0=0.0d0
       end if
c      Neutron 3P2:
       if ((inu_n3p2_pbf.eq.1).and.(i.le.isf)) then
          call nu_n3p2_B_pbf(t,tcn(i),mstn(i),kfn(i),qpbf_n3p2)
          qpbf_n3p2=qpbf_n3p2*fhad(i)
       else
          qpbf_n3p2=0.0d0
       end if
c      Protons:
       if (inu_p_pbf.eq.1) then
          call nu_1s0_pbf(t,tcp(i),mstp(i),kfp(i),qpbf_p1s0)
          qpbf_p1s0=qpbf_p1s0*fhad(i)
       else
          qpbf_p1s0=0.0d0
       end if
c      Quarks: TO BE INCLUDED !!!!!!!!
       qpbf_q=0.0d0
       qpbf_q=qpbf_q*(1.d0-fhad(i))
      else
       qpbf_n1s0=0.0d0
       qpbf_n3p2=0.0d0
       qpbf_p1s0=0.0d0
       qpbf_q=0.0d0
      end if
c *** ADDING EVERYTHING:
      qtot=
     1   qnpb+qplasma+qsynch+qbubble+qpair+qphoto+qbrem_nn+
c FORTIN Jan 26
c     1   qeebrem+qnpb+qplasma+qsynch+qbubble+qpair+qphoto+qbrem_nn+
c FORTIN Jan 26
     2   qmurca_nucl+qbrem_nucl+qmurca_hyp+qbrem_hyp+
     3   qdurca_np+qdurca_lap+qdurca_smn+qdurca_smla+qdurca_sms0+
     4   qfast+
     5   qdurca_q+qmurca_q+qstrange+
     6   qpbf_n1s0+qpbf_n3p2+qpbf_p1s0+qpbf_q+
     7   qdurca_crust

!       write(95,*)rho,qtot,qeebrem,qnpb,qplasma,qsynch,
!c                   1   2    3       4    5       6OFF
!     1   qbubble,    qpair,  qphoto,qbrem_nn,
!c          7OFF      8OFF     9OFF      10
!     2   qmurca_nucl,qbrem_nucl,
!c          11           12        

!c OFF
!c        qmurca_hyp,qbrem_hyp,
!c          11           12        13      14
!     3   qdurca_np,qdurca_lap,qdurca_smn,qdurca_smla,
!c          15           16        17      18
!     4   qdurca_sms0,qfast,qdurca_q,qmurca_q,qstrange,
!c          19           20    21      22      23
!     6   qpbf_n1s0,qpbf_n3p2,qpbf_p1s0,qpbf_q,qdurca_crust
!c          24           25    26        27n      28n


CC MODIF DEC 15
cc      qtot=qdurca_np+qmurca_nucl


c*****

c Ofengeim
cc      if (inu_crust.eq.0) then
cc      if (rho.lt.rhocore) qtot=0.
c DURCA
!      if (rho.ge.rhocore) then
!       q0qnu=1.96e27
!       aqnu=1.01
!       bqnu=1.14
!       cqnu=0.0028
!       gamqnu=2.48
!       nqnu=6
!       kqnu=2
!       pqnu=4
!       funcJ=xrho**(kqnu/3.)*(1.-compactness)**(pqnu/(-2.))
!     1   *(1.+cqnu*xrho**(gamqnu-1.))
!     1    **((pqnu*gamqnu-kqnu/3.)/(gamqnu-1.))
!     1       /(1.-bqnu*compactness)**(0.5)
!       qnuoftot=q0qnu*aqnu*radius**3.*(T/1e9)**nqnu*funcJ !-> OK. checked vs fig 4
!       qtot=qnuoftot/dvcore
!       write(*,*)'qnuoftot',qnuoftot,log10(qnuoftot)
!!       write(*,*)log10(qnuoftot)
!      endif
!c MURCA
!      if (rho.ge.rhocore) then
!       q0qnu=1.17e21
!       aqnu=1.12
!       bqnu=1.14
!       cqnu=0.006
!       gamqnu=2.45
!       nqnu=8
!       kqnu=2
!       pqnu=6
!       funcJ=xrho**(kqnu/3.)*(1.-compactness)**(pqnu/(-2.))
!     1   *(1.+cqnu*xrho**(gamqnu-1.))
!     1    **((pqnu*gamqnu-kqnu/3.)/(gamqnu-1.))
!     1       /(1.-bqnu*compactness)**(0.5)
!       qnuoftot=q0qnu*aqnu*radius**3.*(T/1e9)**nqnu*funcJ !-> OK. checked vs fig 4
!       qtot=qnuoftot/dvcore
!       write(*,*)'qnuoftot',qnuoftot,log10(qnuoftot)
!      endif


c      endif

      if (debug.ge.2.) print *,'Exiting subroutine `neutrino'' '
      return

      end







