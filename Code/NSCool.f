c *********************************************************************
c *********************************************************************
c DEBUGGING STRATEGY:
c debug.eq.0.9d0 : prints during reading of input parameter file `I*.dat'
c debug.ge.1.    : prints program steps
c debug.ge.2.    : prints program steps and entering and exiting main subroutines
c 
c
c
c debug.eq.-50. : print boundary condition details
c
c *********************************************************************
c *********************************************************************
c                             WARNING:
c
c          Index "i" for radial zones: my convention is that
c          - Temp is defined at ODD zone numbers
c          - Lum  is defined at EVEN zone numbers
c
c          i runs from imin=0 to imax
c          - imax MUST BE AN ODD NUMBER:
c            Temp is defined there and is used for the outer boundary condition
c          - icore, idrip and ienv should also be odd
c          - imin=0:
c            Lum is defined in the center of the star and Lum(0)=0
c *********************************************************************
c *********************************************************************
      program NSCool

      implicit real*8 (a-h,k-z)

      parameter (hbar=1.0547266d-27,e=4.803d-10,kb=1.380658d-16)
      parameter (g=6.67259d-8,c=2.99792458d10)
      parameter (msol=1.98892d33,lsol=3.828d33)
      parameter (year=3.1556925d7)
      parameter (pi=3.1415926535d0)

      INCLUDE 'size.inc.f'
      INCLUDE 'rho_limits.inc.f'
      INCLUDE 'gamma_limits.inc.f'
      INCLUDE 'files.inc.f'

      INCLUDE 'profile_star.inc.f'
      INCLUDE 'profile_comp.inc.f'
      INCLUDE 'pairing.inc.f'
      INCLUDE 'fermi.inc.f'
      INCLUDE 'spec_heat.inc.f'
      INCLUDE 'mag_field.inc.f'
      INCLUDE 'control_nu.inc.f'
      INCLUDE 'control_con.inc.f'
      INCLUDE 'control_heat.inc.f'
      INCLUDE 'accretion.inc.f'

      common/stuff/time,istep

c**** Define real*4 plot variable for pgplot::

c HERE DANY cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      INCLUDE 'Plot/Plot_def.inc.f'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c**** Input/Output:

      character*90 filename
      character*79 Model
      character*80  f_concryst
      character*100 f_Bound,f_Pairing,f_Neutrino,f_Conduct,
     x              f_Heat,f_Bfield,f_Accretion,f_Strange,
     x              f_Stestai
      dimension tprint(0:5000),mdot_tab(0:5000)
      dimension m_dot0_tab(0:20),t_acc0_tab(0:20)
      dimension t_acc1_tab(0:20),t_acc2_tab(0:20)
      common/print/pcore,model
      integer pscreen

c**** Heating:

c     Deep crustal heating:
      dimension q_deep_crust(0:isize)
      common/deep_crust/q_deep_crust
c     Hydrogen burning in the envelope
      dimension q_h_burning(0:isize)
      common/h_burning/q_h_burning
c     Old stuff: MAY NOT WORK !
      real*8 j_44,j_heat
      dimension qdeposit(0:isize)
      dimension mag_field(0:isize),j_heat(0:isize)
      common/others/mag_field,j_heat
      common/deposit/total_heat,t_dep,del_t_dep,qdeposit,i_dep

c**** Boundary condition:

      character*100 f_TbTs
      common/bound/tb_acc0,f_TbTs
      common/gravity/gs14,compactness,xrho,radius,rcore,dvcore

C***** Local variable used only in the main program

      dimension orad(0:isize),bar1(0:isize),obar(0:isize),
     2          rrho1(0:isize),orrho(0:isize)
      dimension ephi(0:isize),e2phi(0:isize),a2ephin(0:isize),
     2          dephi(0:isize)
      dimension temp(0:isize),otemp(0:isize),ntemp(0:isize)
      dimension ntemp1(0:isize),delt(0:isize),dtemp(0:isize)
      dimension lum(0:isize),olum(0:isize),nlum(0:isize)
      dimension dell(0:isize),dlum(0:isize)

      dimension lambda(0:isize),lambda1(0:isize)
      dimension kappa(0:isize),kappa1(0:isize)
      dimension qnu(0:isize),qnu1(0:isize)
      dimension qqq(0:isize),qqq1(0:isize)
      dimension qeebrem(0:isize),
     1          qnpb(0:isize),qplasma(0:isize),qsynch(0:isize),
     1          qbubble(0:isize),qpair(0:isize),qphoto(0:isize),
     2          qbrem_nn(0:isize),
     3          qmurca_nucl(0:isize),qbrem_nucl(0:isize),
     4          qmurca_hyp(0:isize),qbrem_hyp(0:isize),
     5          qdurca_np(0:isize),
     5          qdurca_lap(0:isize),qdurca_smn(0:isize),
     6          qdurca_smla(0:isize),qdurca_sms0(0:isize),
     7          qfast(0:isize),
     8          qdurca_q(0:isize),qmurca_q(0:isize),
     9          qpbf_n1s0(0:isize),qpbf_n3p2(0:isize),
     x          qpbf_p1s0(0:isize),qpbf_q(0:isize),
     x          qdurca_crust(0:isize)
      dimension heat(0:isize),heat1(0:isize)
      dimension cv(0:isize),cv1(0:isize)
      dimension cv_n(0:isize),cv_p(0:isize),cv_e(0:isize),cv_m(0:isize),
     1      cv_l(0:isize),cv_sm(0:isize),cv_s0(0:isize),cv_sp(0:isize),
     2      cv_q(0:isize),cv_ion(0:isize)


      dimension gamma(0:isize),cryst(0:isize)

      dimension radts(0:isize),tts(0:isize),rrhots(0:isize)
     1	,taucts(0:isize),taub(0:isize),taup(0:isize),
     2   rts(0:isize),taudts(0:isize)


c Accuracy checking:
c The equation nlum+fp*a2ephi*dtemp is very delicate when the star is
c isothermal: to get convinced that it is solved use:
c      real*16 nlum,ntemp,dtemp,ff
c If these are only real*8 then it may look like the equation is not
c solved properly. However it is, by the Henyey method as soon as
c dell & delt are small enough. That's why it is considered as solved
c either when it is numerically solved (maxff1<chff1) or ratiol<mratl

      dimension fp(0:isize),fq(0:isize),fr(0:isize)
      dimension fp1(0:isize),fq1(0:isize),fr1(0:isize)
      dimension dfp(0:isize),dfq(0:isize),dfr(0:isize)
      dimension fa(0:isize),fb(0:isize),fc(0:isize),ff(0:isize)
      dimension fj(0:isize),fk(0:isize)

      dimension f_gr_field(0:isize),g_gr_field(0:isize),
     1          h_gr_field(0:isize)

c Auxiliary variables:
      character*5 what



c *********************************************************************
c *********************************************************************
c **********************     LET'S GO !      **************************
c *********************************************************************
c *********************************************************************

      i_model=0

c *********************************************************************
c *********************************************************************
c ************   BEGINNING OF A NEW MODEL CALCULATION   ***************
c *********************************************************************
c *********************************************************************

 1234 continue

      i_model=i_model+1

      idt=1
      htot=0.0
      contraction=0.d0

c *********************************************************************
c *****************     GET INPUT MODEL FILES    **********************
c *********************************************************************

      if (i_model.eq.1) then
c ***  Choose between two input: **************************************
c      Ask for the input file *****************************************
c       write(6,*) 'Input master file ='
c       read(5,*)filename
c ***  Can add here the directory where "Cool_*.in" is:
c       filename='CasA/'//filename
c       filename='Defence/'//filename
c       filename='Cv2/'//filename
c       filename='Sym_e/APR/'//filename
c       filename='Model_1/'//filename
c      filename='Envelope/XTE/'//filename
c      filename='Acc_steady_state/APR+HZ/'//filename
c      filename='QPXRT/IGR_J17480-2446/'//filename
c      filename='QPXRT/MXB_1659-29/'//filename
      filename='DCH_HZ08.in'
c ***  Or define it completely here: **********************************
c      filename='Acc_steady_state/APR+HZ/1.4_APR+HZ_SFBa.in'
c       filename='Ter5/test.in'
c       filename='Ofengeim/BSk21_2.0_O.in'
c       filename='Beznogov/Isolated_Ioffe.in'
c       filename='OPUS15/DCH_BSK21.in'
c       write(6,*)'Using as input: ',filename
c**********************************************************************
       open(unit=15,file=filename,status='old')
      else if ( (i_acc.eq.3) .and. (i_model.le.i_modelmax) .and. 
     1         (i_model.gt.1)) then
          go to 9996
      else
       read(15,*,end=9997,err=9997)
      end if
      read(15,*,end=9997,err=9997)version
      write(*,*)version
      if (version.eq.'STR') then
       istrange=1
      else
       istrange=0
      end if
c *** BASIC MODEL FILES: **********************************************
      read(15,*,end=9997,err=9997)
      read(15,*,end=9997,err=9997)f_crusteos
      read(15,*,end=9997,err=9997)f_stareos
      read(15,*,end=9997,err=9997)f_profile
c *** OTHER MODEL FILES: **********************************************
      read(15,*,end=9997,err=9997)
      read(15,*,end=9997,err=9997)f_Structure
      read(15,*,end=9997,err=9997)f_Bound
      read(15,*,end=9997,err=9997)f_Pairing
      read(15,*,end=9997,err=9997)f_Neutrino
      read(15,*,end=9997,err=9997)f_Conduct
      read(15,*,end=9997,err=9997)f_Heat
      read(15,*,end=9997,err=9997)f_Bfield
      read(15,*,end=9997,err=9997)f_Accretion
      if (istrange.eq.1) then
       read(15,*,end=9997,err=9997)f_Strange
      end if
c *** OUTPUT FILES: ***************************************************
      read(15,*)                    
      read(15,*,end=9997,err=9997)f_i
      read(15,*,end=9997,err=9997)f_Teff
      read(15,*,end=9997,err=9997)f_Temp
      read(15,*,end=9997,err=9997)f_Star
      read(15,*,end=9997,err=9997)f_SteStao
c**********************************************************************
c *** PRINT ON THE SCREEN THE FILES:
c Notice: pscreen will be read from file "I.dat"
      pscreen=1
      if (pscreen.gt.0) then
       print *,'-------------------------------------------------'
       print *,'Here are the files: '
       write(6,*)version(1:istrlen(version))
       write(6,*)f_crusteos(1:istrlen(f_crusteos))
       write(6,*)f_stareos(1:istrlen(f_stareos))
       write(6,*)f_profile(1:istrlen(f_profile))
       write(6,*)f_structure(1:istrlen(f_structure))
       write(6,*)f_Bound(1:istrlen(f_Bound))
       write(6,*)f_Pairing(1:istrlen(f_Pairing))
       write(6,*)f_Neutrino(1:istrlen(f_Neutrino))
       write(6,*)f_Heat(1:istrlen(f_Heat))
       write(6,*)f_Bfield(1:istrlen(f_Bfield))
       write(6,*)f_Accretion(1:istrlen(f_Accretion))
       if (istrange.eq.1) write(6,*)f_Strange(1:istrlen(f_Strange))
       write(6,*)f_i(1:istrlen(f_i))
       write(6,*)f_Teff(1:istrlen(f_Teff))
       write(6,*)f_Temp(1:istrlen(f_Temp))
       write(6,*)f_Star(1:istrlen(f_Star))
       write(6,*)f_SteStao(1:istrlen(f_SteStao))
       print *,'-------------------------------------------------'
       print *,'LET''S GO !'
       print *,'-------------------------------------------------'
       print '(1a28,1a50)','****************************',
     1  '**************************************************'
      end if


cc Cas A test 
c       write(6,*) 'Time DURCA [yr]?'
c       read(5,*)time_exoy
c       time_exo=time_exoy*3.15576d7 
c       write(6,*) 'Rho DURCA ?'
c       read(5,*)rhoexor
cc

 9996 continue 

c *********************************************************************
c *****************     READ THE ABOVE FILES     **********************
c *********************************************************************

      INCLUDE 'NSCool_READ.inc.f'

c *********************************************************************
c *****************        INITIALIZATION        **********************
c *********************************************************************
      if (debug.ge.1.) print *,'Initializing'

c *********************************************************************
c *** Get the time independent pieces of physics: *********************
c *********************************************************************
c     get_core_chemistry MUST be called BEFORE get_crust_chemistry

      if (istrange.eq.1) then
       call grid_strange(idec,rhocore,rhodrip,rhoenv,rhosurf,Model)
       call get_core_chemistry_strange(debug)
      else
       call grid(idec,rhocore,rhodrip,rhoenv,rhosurf,Model)
       call get_core_chemistry
      end if
      call get_crust_chemistry(debug)
      call get_Fermi_momenta
      call get_effective_masses
      write(*,*)mstn(icore),bar(icore)
      call get_spec_heat_degenerate
      call get_Tc
      call get_degenerate_density




c *********************************************************************
c ***** Calculate the T-independent coefficients **********************
c *********************************************************************
      if (debug.ge.1.) print *,'Calculating T-independent coeff.'

      do i=0,imax
       ephi(i)=dexp(phi(i))
       e2phi(i)=ephi(i)**2
       a2ephin(i)=(4.d0*pi*rad(i)**2)**2*ephi(i)
      end do

      dephi(0)=0.d0
      do i=1,imax-1
       dephi(i)=(ephi(i+1)-ephi(i-1))/(rad(i+1)-rad(i-1))
      end do
      dephi(imax)=dephi(imax-1)

      radius=rad(imax)
      root=dsqrt(1.-2.d0*g*msol*emas(imax)/rad(imax)/c**2)
      factor=root/(4.d0*pi*rad(imax)**2)/6.022d23*1.d39
      constant=4.d0*pi*g*msol*emas(imax)*4.d0/3.d0*5.67d-5*
     1         e2phi(imax)/pres(imax)/root

      gs=g*msol*emas(imax)/rad(imax)**2/root
      gs14=gs/1.d14

      compactness=2.d0*g*msol*emas(imax)/rad(imax)/c**2
      xrho=msol*emas(imax)/(2.8e14*rad(imax)**3.)
      rcore=rad(icore)

c Definition of eta
      etabound=etabound*gs14**2.

c *********************************************************************
c *** Initialize some more stuff: *************************************
c *********************************************************************

      if (i_heat_deep_crust.eq.1) then
       call initialize_heating_deep_crust
      end if
      if (i_heat_deep_crust.eq.2) then
       call initialize_heating_deep_crust_mb
      end if      
      if (i_heat_deep_crust.eq.20) then
       call initialize_heating_deep_crust_bsk20
      end if
      if (i_heat_deep_crust.eq.21) then
       call initialize_heating_deep_crust_bsk21
      end if
      if (i_heat_deposit.eq.1) then
       call initialize_heating_deposit
      end if
      if (i_heat_convert.eq.1) then
       call initialize_heating_convert(MeV_neutron)
      end if
      if (i_heat_h.eq.1) then
       call initialize_heating_h_burning
      end if
      call get_beta(emas(imax),rad(imax),beta_rot,m_i)
c      if (ifield.ne.0) call initialize_dipole(i0,i1)
      if (i_acc.eq.3) then 
          m_dot0=mdot_tab(i_model)
	  write(90,*)m_dot0,i_model
      endif
	if(i_acc.ne.4) call initialize_accretion_rate


      if (i_heat_vortex_creep.eq.1) then
       call initialize_heating_vortex_creep(imax,rad,rrho,tcn,dvol,j_44)
      end if
c BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
c This is to initialize the Bfield safely 
c when Bfield will be eliminated:
      do i=imax,1,-2
       nbfield2(i)=bfield0
      end do
c BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
c      INCLUDE 'Bfield/Bfield_1.inc.f'
c BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

c *********************************************************************
c ********* Plot or Print out some physics ****************************
c *********************************************************************

c      INCLUDE 'Plot/Plot_Qnu.inc.f'

c *********************************************************************
c ********* Calculate the initial Temp and Lum profiles ***************
c *********************************************************************

      INCLUDE 'NSCool_INIT_TPROF.inc.f'

c BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
c      INCLUDE 'Bfield/Bfield_2.inc.f'
c BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

c *********************************************************************
c ****************** Open print out files *****************************
c *********************************************************************

      INCLUDE 'NSCool_OPEN.inc.f'
 

c *********************************************************************
c *********************************************************************
c ************************     COOLING     ****************************
c *********************************************************************
c *********************************************************************

c *********************************************************************
      time=time0       ! Initialize the time
      icycle=0         ! Initialize the counter for accretion cycles
c *********************************************************************

c *********************************************************************
c      itprint=0       ! To print out the initial T and L profiles
      itprint=1        ! To print out only at the required times
c *********************************************************************

c *********************************************************************
c     THIS IS THE MAIN TIME LOOP:
      do 9999 istep=1,istepmax
c *********************************************************************

        debug=0.
        if (istep.ge.istep_debug) debug=debug_keep
        if (debug.ge.1.) print *,'Going: istep=',istep

c *********************************************************************
2345  itrial=0         ! Branch back here in case:
                       !  - Too many iteration in Newton-Raphson
                       !  - Envelope boundary condition cannot be solved
                       !  - Temp has changed too much
c *********************************************************************

      ratiot=1.d-2
      ratiol=1.d-2

cc
      if (i_acc.eq.4) then

       if ( (time+dtime) .le. t_acc0_tab(1) ) then
	t_acc0=t_acc0_tab(1)
	t_acc2=t_acc2_tab(1)
	m_dot0=m_dot0_tab(1)
	t_burst=time+dtime-t_acc0
	icycle=0
      else if ( ( (time+dtime) .ge. t_acc0_tab(1) ) .and.
     1           ( (time+dtime) .le. ( t_acc0_tab(1)+t_acc2_tab(1) ) )
     2         ) then
	t_acc0=t_acc0_tab(1)
	t_acc2=t_acc2_tab(1)
	m_dot0=m_dot0_tab(1)
	t_burst=time+dtime-t_acc0
      call accretion_rate_multiple(time+dtime,dtime,m_dot)
	icycle=1
       else if ( ( (time+dtime) .ge. ( t_acc0_tab(1)+t_acc2_tab(1) ) ) 
     1   .and.   ( (time+dtime) .le. ( t_acc0_tab(2) ) )
     2         ) then
	t_acc0=t_acc0_tab(2)
	t_acc2=t_acc2_tab(2)
	m_dot0=m_dot0_tab(2)
	t_burst=time+dtime-t_acc0
      call accretion_rate_multiple(time+dtime,dtime,m_dot)
	icycle=0
       else if ( ( (time+dtime) .ge. t_acc0_tab(2) ) .and.
     1           ( (time+dtime) .le. ( t_acc0_tab(2)+t_acc2_tab(2) ) )
     2         ) then
	t_acc0=t_acc0_tab(2)
	t_acc2=t_acc2_tab(2)
	m_dot0=m_dot0_tab(2)
	t_burst=time+dtime-t_acc0
      call accretion_rate_multiple(time+dtime,dtime,m_dot)
	icycle=1
	else
	t_acc0=t_acc0_tab(2)
	t_acc2=t_acc2_tab(2)
	m_dot=0.d0
        icycle=1
	t_burst=time+dtime-t_acc0
	endif

	else
c *** Accretion rate: *************************************************
      icycle_old=icycle
      if ((i_acc.eq.1).or.(i_acc.eq.2)) then
       if (time+dtime.ge.t_acc0) then
        icycle=int((time+dtime-t_acc0)/t_acc1)+1
        t_burst= time+dtime -t_acc0 - float(icycle-1)*t_acc1
                 ! t_burst = time since beginning of burst
       else
        icycle=0
        t_burst=0.d0
       end if
      end if

      call accretion_rate(time+dtime,dtime,m_dot)
	endif
c *** Gravitational energy of the infalling accreted material *********
      if (ilgrav.eq.1) then
cc Prescription from Miralda-Escude et al. 90, equation (5)
cc Local Lgrav at the surface of the star
	Lgrav=1/2.* g *msol*emas(imax) /radius *m_dot /ephi(imax)
     1        /lsol*e2phi(imax-1)
      endif
c ***

      call accretion_velocity(m_dot)
c *********************************************************************

c ---------------------------------------------------------------------
      if (pscreen.ge.2) then
       print '(1a28,1a50)','****************************',
     1  '**************************************************'
       print '(2a10,1i5,1a53)','**********','step#=',istep,
     1  '**************************************************'
       print '(1a28,1a50)','****************************',
     1  '**************************************************'
c       read(5,*)        ! Wait for <ENTER> befor printing out
       print '(2(a10,1p1e10.3),a30,0p1f6.3)',
     1                'time =',(time+dtime)/year,
     2               'dtime =',dtime/year,
     3        'dtime/odtime =',dtime/odtime
       print *
       if  ((dtime.ge.t_acc2).and.(icycle.ne.0)) then
       print '(1a28,1a50)','************problem**********',
     1  '**************************************************'
       endif
       if (chtemp.eq.1.) then
        print '(a42,0p1f5.2,a9,1p1e9.2,a3,1p1e9.2)',
     1        'dtime limited by TEMP change, max_dtemp =',max_dtemp,
     2        'at rho=',rrho(icht),'T=',temp(icht)
       end if
       if (chstoke.eq.1.) then
        print '(a42,0p1f5.2,a9,1p1e9.2,a3,1p1e9.2)',
     1        'dtime limited by STOKE change, mdstoke =',mdstoke,
     2        'at rho=',rrho(ichs),'S=',stoke(ichs)
       end if
       if (chtrial.eq.1) then
        print '(a40)',
     1        '   dtime limited by ITRIAL'
       end if
      end if
c ---------------------------------------------------------------------

c *********************************************************************
c ***** Calculate ntemp & nlum for first guess ************************
c *********************************************************************

      if (debug.ge.1.) print *,'Guessing NLum & NTemp'
      coeff_int=0.8d0

      do i=1,imax,2
       ntemp(i)=temp(i)+coeff_int*(temp(i)-otemp(i))*dtime/odtime
      end do

      dtemp(0)=0.d0
      do i=2,imax-1,2
       dtemp(i)=(ntemp(i+1)-ntemp(i-1))/(debar(i)+debar(i+1))
      end do

      nlum(0)=0.
      do i=2,imax-1,2
       nlum(i)=lum(i)+coeff_int*(lum(i)-olum(i))*dtime/odtime
      end do
      do i=1,imax-2,2
       dlum(i)=(nlum(i+1)-nlum(i-1))/(debar(i)+debar(i+1))
      end do
c look +++++++++++++++++++++++++++++++++++++++++++++
      dlum(imax)=0.d0
c      dlum(imax)=dlum(imax-2)
c ++++++++++++++++++++++++++++++++++++++++++++++++++

c BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
c      INCLUDE 'Bfield/Bfield_3.inc.f'
c BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

c ---------------------------------------------------------------------
      if (pscreen.eq.3)then
       print '(1i3,1a5)', 0,     '-----'
      end if
c ---------------------------------------------------------------------

c *********************************************************************
c *********************************************************************
c ***** Branch here if new trial **************************************
c *********************************************************************
c *********************************************************************

c *********************************************************************
2000  itrial=itrial+1                 ! This is the Newton-Raphson loop
c *********************************************************************
c	write(91,*)'536',itrial,1013,ntemp(1013)
       if (itrial.eq.itrial_max+1)then
        tcut=dsqrt(scale_dt0)
        if (time.le.1.e5) tcut=dsqrt(scale_dt1)
        dtime=dtime/tcut
        goto 2345
       end if

      do i=0,imax
       rad(i)=orad(i)
       rrho(i)=orrho(i)
       rrho1(i)=rrho(i)
      end do   

      oteffective=teffective
c ---------------------------------------------------------------------
       if (pscreen.eq.3)then
        print '(1i3,1a5)', itrial,'-----'
       end if
c ---------------------------------------------------------------------

c *********************************************************************
c ***** Calculate the new density in inner envelope at ntemp **********
c *********************************************************************

      if (debug.ge.1.) print *,'Calculating new density at NTemp'
c	write(91,*)'562',1013,ntemp(1013)
      do i=imax-1,ienv+1,-2
       x=(dlog(rrho(i+1))-dlog(rrho(i)))/
     1   (dlog(rrho(i+1))-dlog(rrho(i-1)))
       y=(dlog(rrho(i))-dlog(rrho(i-1)))/    
     1   (dlog(rrho(i+1))-dlog(rrho(i-1)))
       ltemp=y*log(ntemp(i+1))+x*log(ntemp(i-1))
       ntemp(i)=dexp(ltemp)
c	write(91,*)'569',1013,ntemp(1013)
      end do
c	write(91,*)'571',1013,ntemp(1013)

      do i=ienv,imax
       if (ntemp(i).lt.tcon) then      
        call density(ntemp(i)/ephi(i),
     1               pres(i),a_cell(i),z_ion(i),rrho(i))
        rrho(i)=min(rrho(i),rhod(i))
        bar(i)=1.d-39 * 6.022d23*rrho(i)
        dr=debar(i)/rrho(i)*factor
        rad(i+1)=rad(i)+dr
        a2ephin(i)=(4.d0*pi*rad(i)**2)**2*ephi(i)
       end if
      end do

c *********************************************************************
c ***** Calculate the physical parameters at ntemp ********************
c *********************************************************************

      if (debug.ge.1.) print *,'Calculating physics at NTemp'

cc Cooling time scales cc
c	itmax=31
c	do it=1,itmax
c	T=10**(3./(itmax-1)*(it-1)+7.)
c	tauptest=0.d0
c	taubtest=0.d0
c	tts(:)=0.d0
c	radts(:)=0.d0
c	rts(:)=0.d0		
c	rrhots(:)=0.d0
c	taucts(:)=0.d0
c	taudts(:)=0.d0
c	taub(:)=0.d0
c	taup(:)=0.d0
cccccccccccccccccccccccccc
! Physical volume of the core
      dvcore=0
      do i=1,imax,2
       if(rrho(i).ge.rhocore) then
       dvcore=dvcore+(dvol(i)+dvol(i+1))
      end if 
      end do 

c      do i=1,imax,2
c       if(rrho(i).ge.rhocore) then
c       write(90,*)rrho(i),pres(i),bar(i),yelect(i),ymuon(i),yneutr(i)
c     4 ,yprot(i),mstp(i),mstn(i)
c a_cell(i),a_ion(i),z_ion(i)
c      end if 
c      end do 
c      stop
!        write(6,*) 'Temperature ?'
!        read(5,*)t

      do i=1,imax,2
c      do i=imax,1,-2


c***************Original temperature ******
      t=ntemp(i)/ephi(i)
c******************************************
c	t=1e7/ephi(i)

       d=rrho(i)
       a=a_cell(i)
       a1=a_ion(i)
       z=z_ion(i)
       call neutrino(i,t,d,a,z,qnu(i),
     1   qeebrem(i),qnpb(i),qplasma(i),qsynch(i),qbubble(i),
     1   qpair(i),qphoto(i),qbrem_nn(i),
     2   qmurca_nucl(i),qbrem_nucl(i),qmurca_hyp(i),qbrem_hyp(i),
     3   qdurca_np(i),qdurca_lap(i),
     3   qdurca_smn(i),qdurca_smla(i),qdurca_sms0(i),
     4   qfast(i),
     5   qdurca_q(i),qmurca_q(i),
     6   qpbf_n1s0(i),qpbf_n3p2(i),qpbf_p1s0(i),qpbf_q(i),
     7   qdurca_crust(i),m_dot/ephi(i),debug,naa)

cc*** results of the neutrino emissivities calculations ***c
c       write(*,*)'T =',t
CC       write(92,'(11E14.4)')rrho(i),qnu(i),
CC     1   qeebrem(i),qnpb(i),qplasma(i),qsynch(i),qbubble(i),
CC     1   qpair(i),qphoto(i),qbrem_nn(i),qpbf_n1s0(i)
c       write(90,'(11E13.3)')rrho(i),
c     1     qplasma(i),qeebrem(i),qnpb(i),qbrem_nn(i),
c     2     qpbf_n1s0(i),qpbf_p1s0(i),
c     3     qmurca_nucl(i),qbrem_nucl(i),
c     4     qfast(i),qnu(i)
c         write(90,'(5E14.4)')rrho(i),qnu(i),
c     1                        qpbf_n1s0(i),qpbf_n3p2(i),qpbf_p1s0(i)
CC         write(94,'(I5,3E14.4)')i,rrho(i),qnu(i),qdurca_crust(i)
c         write(92,'(2E14.4)')rrho(i),qnu(i)*e2phi(i)
!         write(92,*)rrho(i),qnu(i),qeebrem(i),qnpb(i),
!     1   qplasma(i),qsynch(i),qbubble(i),qpair(i),qphoto(i),qbrem_nn(i),
!     2   qmurca_nucl(i),qbrem_nucl(i),
!     4   qfast(i),
!     6   qpbf_n1s0(i),qpbf_n3p2(i),qpbf_p1s0(i),qpbf_q(i),
!     7   qdurca_crust(i)

c         write(92,*)rrho(i),qnu(i),qeebrem(i),qnpb(i),qbrem_nn(i),
c     2   qmurca_nucl(i),qbrem_nucl(i)
c         write(94,'(3E14.4)')rrho(i),yelect(i)*bar(i),ymuon(i)*bar(i)
cc*********************************************************c

c          m_dot=m_dot_av

       call heating(i,time,dtime,t,d,a,a1,z,tcp(i),
     1              m_dot/ephi(i),i_eleconb,
     2              ephi(i),dephi(i),heat(i))     
c          lq_av=lq_av+e2phi(i)*heat(i)*(dvol(i)+dvol(i+1))

       qqq(i)=qnu(i)-heat(i)
c	t=1e9
       call specheat_hfb(i,t,d,a,a1,z,cv(i),
     1      cv_n(i),cv_p(i),cv_e(i),cv_m(i),
     2      cv_l(i),cv_sm(i),cv_s0(i),cv_sp(i),
     3      cv_q(i),cv_ion(i))

cc*** results of the heat capacities calculations ***c
c       write(*,*)'T =',t
c      write(90,'(6E14.4)')rrho(i),cv(i),
c     1      cv_n(i),cv_p(i),cv_e(i),cv_ion(i)


!      write(96,*)rrho(i),cv(i),
!c                   1      2  
!     1      cv_n(i),cv_p(i),cv_e(i),cv_m(i),
!c             3       4      5        6
!     2      cv_l(i),cv_sm(i),cv_s0(i),cv_sp(i),
!c             7       8      9        10
!     3      cv_q(i),cv_ion(i)
!c             11       12      
c       write(91,'(i7,3E14.4)')i,rrho(i),cv(i)
c         write(90,'(i7,7E14.4)')i,rad(i),rrho(i),
c     1                    cv_n(i),cv(i),a,a1,z
c       write(90,'(4E14.4)')rrho(i),cv(i),cv_e(i),cv_ion(i)
c*********************************************************c
c       write(90,'(2E14.4)')rrho(i),cv(i)
ccccc
c	write(91,*)'658',i,t
       call conduct(i,t,d,a,a1,z,qimp,nbfield2(i),
     1              sig,lambda(i),debug,
     2              nu_e_s,nu_e_l)

cc*** results of the thermal conductivities calculations ***c
c       write(*,*)'T =',t
c	if (d.le.rrho(icore)) then
CC       write(91,'(5E14.4)')d,a,a1,z,lambda(i)
c	endif
!        write(91,'(2E14.4)')rrho(i),lambda(i)

ccc*********************************************************c
       call opacity(t,d,a,z,kappa(i))
       acd=7.56d-15*c/(3.d0*kappa(i)*d)
       fp(i)=(lambda(i)+4.d0*acd*t**3)*bar(i)/lsol
       fq(i)=bar(i)/cv(i)*lsol
       fr(i)=e2phi(i)*qqq(i)/cv(i)
     1       -ephi(i)*pres(i)/cv(i)
     2       *(dlog(rrho(i))-dlog(orrho(i)))/dtime*contraction

!        write(91,'(3E14.4)')rrho(i),cv(i),lambda(i)

	enddo

!      cv_core=0
!      do i=1,imax,2
!       if(rrho(i).ge.rhocore) then
!       cv_core=cv_core+cv(i)*(dvol(i)+dvol(i+1))
!      end if 
!      end do 
!       write(*,*)'cv_core',cv_core

!      qnu_core=0
!      lnu_core=0
!      do i=1,imax,2
!       if(rrho(i).ge.rhocore) then
!       qnu_core=qnu_core+qnu(i)*(dvol(i)+dvol(i+1))
!       lnu_core=lnu_core+e2phi(i)*qnu(i)*(dvol(i)+dvol(i+1))
!      end if 
!      end do 

c       write(*,*)'qnu_core',qnu_core,log10(qnu_core)
c       write(*,*)'qnu_core*e2phi',lnu_core,log10(lnu_core)
c       write(*,*)'rhocore', rhocore
c	stop

cc*********************** Cooling time scales ********************
c
cc Thickness of the inner crust
c	ric=-rad(idrip)+rad(imax)
c	ric=-rad(1)+rad(imax)
cc Number of zones for the calculations of the time scales
c	its=200
cc Thickness of the zones for the calculations of the time scales
c	lts=ric/its
cc Radius of the zones for the time scale calculations
c	do is=1,its+1
c	radts(is)=rad(idrip)+(is-1)*lts
c	radts(is)=rad(1)+(is-1)*lts
c	enddo
c
c	j=idrip
c	j=1
c	do is=1,its+1
cc index of the lower shell 
c	itsmin=j
c	do while ( rad(j).le.radts(is) ) 
c	j=j+2
cc since the cv are defined for even indices
c	enddo
cc index of the upper shell
c	itsmax=j
cc in the case of the outermost shell, its boundary correspond to the boundary of the outer layer
c	if (is.eq.its+1) itsmax=imax
c
c	if (itsmax.ne.itsmin) then
cc radius of the middle of the zone
c	radtsmid=(rad(itsmin)+rad(itsmax)) / 2.
c          x=( dlog( radtsmid )   - dlog( rad(itsmin) ) )/
c     1      ( dlog( rad(itsmax) )- dlog( rad(itsmin) ) )
c          y=( dlog( rad(itsmax) )- dlog( radtsmid ) )/
c     1      ( dlog( rad(itsmax) )- dlog( rad(itsmin) ) )
c
c	  logrrhomid   = y * dlog( rrho(itsmax) )
c     1                 + x * dlog( rrho(itsmin) )
c	  logcvmid     = y * dlog( cv(itsmax) )
c     1                 + x * dlog( cv(itsmin) )
c	  loglambdamid = y * dlog( lambda(itsmax) )
c     1                 + x * dlog( lambda(itsmin) )
c	  logqnumid    = y * dlog( qnu(itsmax) )
c     1                 + x * dlog( qnu(itsmin) )
c
c	  rrhomid=dexp(logrrhomid)
c	  cvmid=dexp(logcvmid)
c	  lambdamid=dexp(loglambdamid)
c	  qnumid=dexp(logqnumid)
cc	write(92,*)itsmin,itsmax
cc	write(92,'(3E14.4)') rrho(itsmin),rrhomid,rrho(itsmax)
cc	write(92,'(5E14.4)') cv(itsmin),cvmid,cv(itsmax)
cc	write(92,'(3E14.4)') lambda(itsmin),lambdamid,lambda(itsmax)
cc	write(92,'(3E14.4)') qnu(itsmin),qnumid,qnu(itsmax)
c	taudmid=cvmid*t/qnumid
c	taucmid=( lts )**2*cvmid/lambdamid
c	tts(is)=t
c	rrhots(is)=rrhomid
c	rts(is)=radtsmid
c	taucts(is)=taucmid
c	taudts(is)=taudmid
c	else
c	tts(is)=t
c	rrhots(is)=rrho(itsmin)
c	rts(is)=rad(itsmin)
c	taucts(is)=0
c	taudts(is)=cv(itsmax)*t/qnu(itsmax)
c 	endif
c	enddo
c
c	taup(its+1)=taucts(its+1)*pi**2./4.
c	taub(its+1)=1./pi**2.*(sqrt(taucts(its+1)*pi**2./4.))**2.
c
c	do is=its,1,-1
cc Pizzochero et al. 2002 formula, ApJ 569
c	taup(is)=taup(is+1)+taucts(is)*pi**2./4.
c	tausum=0.d0
c	do isum=its+1,is,-1
c	tausum=sqrt(taucts(isum)*pi**2./4.)+tausum
c	enddo
cc Brown et al. 1998 formula, ApJ L95
c	taub(is)=1./pi**2.*(tausum)**2.
c	enddo
c
c	do is=1,its+1
c	write(93,'(i7,5E14.4)')is,tts(is),rts(is),rrhots(is),taudts(is),
c     1              taub(is)
c	enddo
c
c	write(93,*)
c	enddo
cc**********************************************************************
c	stop


c *********************************************************************
c ***** Calculate the new density at (1-tinc)*ntemp *******************
c *********************************************************************

      if (debug.ge.1.) print *,'Calculating density at NTemp'''
      tinc=max(1.d-12,ratiot/1.d1)
      do i=1,imax,2
       ntemp1(i)=ntemp(i)*(1.d0-tinc)
      end do

      do i=imax-1,ienv+1,-2
       x=(dlog(rrho(i+1))-dlog(rrho(i)))/
     1   (dlog(rrho(i+1))-dlog(rrho(i-1)))
       y=(dlog(rrho(i))-dlog(rrho(i-1)))/    
     1   (dlog(rrho(i+1))-dlog(rrho(i-1)))
       ltemp=y*log(ntemp1(i+1))+x*log(ntemp1(i-1))
       ntemp1(i)=dexp(ltemp)
      end do

      do i=ienv,imax
       if (ntemp(i).lt.tcon) then      
        call density(ntemp1(i)/ephi(i),pres(i),a_ion(i),z_ion(i),
     1               rrho1(i))
        rrho1(i)=min(rrho1(i),rhod(i))
        bar1(i)=1.d-39 * 6.022d23*rrho1(i)
       end if
      end do

c *********************************************************************
c ***** Calculate the physical parameters at (1-tinc)*ntemp ***********
c *********************************************************************

      if (debug.ge.1.) print *,'Calculating physics at NTemp'''
      do i=1,imax,2
       t=ntemp1(i)/ephi(i)
       d=rrho1(i)
       a=a_cell(i)
       a1=a_ion(i)
       z=z_ion(i)
       call neutrino(i,t,d,a,z,qnu1(i),
     1          qn00,qn01,qn02,qn03,qn04,qn05,qn06,qn07,qn08,qn09,q10,
     2               qn11,qn12,qn13,qn14,qn15,qn16,qn17,qn18,qn19,q20,
     3               qn21,qn22,qn23,
     4               qn24,m_dot/ephi(i),debug,naa)
       call heating(i,time,dtime,t,d,a,a1,z,tcp(i),
     1              m_dot/ephi(i),i_eleconb,
     2              ephi(i),dephi(i),heat1(i))
       qqq1(i)=qnu1(i)-heat1(i)
c       call specheat(i,t,d,a,z,cv1(i),
c     1               xx1,xx2,xx3,xx4,xx5,xx6,xx7,xx8,xx9,xx0)

       call specheat_hfb(i,t,d,a,a1,z,cv1(i),
     1      xx1,xx2,xx3,xx4,xx5,xx6,xx7,xx8,xx9,xx0)
       if (debug.ge.1.) print '(a20,i5,1p2e12.3)',
     1      'Calling conduct',i,d,t
       call conduct(i,t,d,a,a1,z,qimp,nbfield2(i),
     1              sig,lambda1(i),debug,
     2                   nu_e_s,nu_e_l)
       if (debug.ge.1.) print *,'Done'
       call opacity(t,d,a,z,kappa1(i))
       acd1=7.56d-15*c/(3.d0*kappa1(i)*d)
       fp1(i)=(lambda1(i)+4.*acd1*t**3)*bar1(i)/lsol
       fq1(i)=bar1(i)/cv1(i)*lsol
       fr1(i)=e2phi(i)*qqq1(i)/cv1(i)
     1        -ephi(i)*pres(i)/cv1(i)
     2        *(dlog(rrho1(i))-dlog(orrho(i)))/dtime*contraction
      end do

c *********************************************************************
c ***** Calculate the derivatives of fp,fq & fr ***********************
c *********************************************************************

      if (debug.ge.1.) print *,'Calculating derivatives of FP, FQ, FR'
       do i=1,imax,2
        t =ntemp(i)
        t1=ntemp1(i)
        dfp(i)=(fp(i)-fp1(i))/(t-t1)
        dfq(i)=(fq(i)-fq1(i))/(t-t1)
        dfr(i)=(fr(i)-fr1(i))/(t-t1)
       end do

c *********************************************************************
c ***** Calculate ff **************************************************
c *********************************************************************

      if (debug.ge.1.) print *,'Calculating FF'
      ff(0)=0.d0
      do i=2,imax-1,2
       ff(i)=nlum(i)+.5d0*(fp(i-1)+fp(i+1))*a2ephin(i)*dtemp(i)
       ff(i-1)=fr(i-1)+fq(i-1)*dlum(i-1)+(ntemp(i-1)-temp(i-1))/dtime
      end do
c look +++++++++++++++++++++++++++++++++++++++++++++
      ff(imax)=0.d0
c      ff(imax)=fr(i-max)+fq(imax)*dlum(imax)+
c     1         (ntemp(imax)-temp(imax))/dtime
c ++++++++++++++++++++++++++++++++++++++++++++++++++

c *********************************************************************
c ***** Matrix inversion for Newton-Raphson method ********************
c *********************************************************************

      if (debug.ge.1.) print *,'Newton-Raphson'
      do i=2,imax-1,2
       fa(i)=.5d0*dfp(i+1)*a2ephin(i)*dtemp(i)+
     1       .5d0*(fp(i+1)+fp(i-1))*a2ephin(i)/(debar(i)+debar(i+1))
       fb(i)=.5d0*dfp(i-1)*a2ephin(i)*dtemp(i)-
     1       .5d0*(fp(i+1)+fp(i-1))*a2ephin(i)/(debar(i)+debar(i+1))
       fc(i)=1.d0
      end do

      do i=1,imax-2,2
       fa(i)=fq(i)/(debar(i)+debar(i+1))
       fb(i)=-fq(i)/(debar(i)+debar(i+1))
       fc(i)=dfr(i)+dfq(i)*dlum(i)+1.d0/dtime
      end do

      fk(1)=-ff(1)/fc(1)
      fj(1)=+fa(1)/fc(1)
      do i=2,imax-1
       fk(i)=-(ff(i)+fb(i)*fk(i-1))/(fc(i)-fb(i)*fj(i-1))
       fj(i)=fa(i)/(fc(i)-fb(i)*fj(i-1))
      end do

c *********************************************************************
c ***************** Boundary condition ********************************
c *********************************************************************

!cc Check cc Jan 25 2016
!      do i=1,30
!       tblog=6.5d0+(i-1.)*0.1d0
!       tb=10.d0**tblog
!       teff0=fteff(tb,ifteff,etabound,bf_r(imax),
!     1           istep,time,ts1,ts2,z_ion(imax),a_ion(imax),rrho(imax),
!     2             debug)
!       write(*,*)tblog,log10(teff0)
!      enddo
!      stop
	ibc=0
       		if(debug.eq.-50.) print *,'Boundary Condition'
       		if(debug.eq.-50.) print *,'Lgrav =', Lgrav
      if (debug.ge.1.) print *,'Boundary Condition'
      if ( (ifteff.ne.15)) then
cc No gravitational energy included
	if ((ilgrav.eq.0).or.(lgrav.eq.0.)) then
       		if(debug.eq.-50.) print *,'NoLgrav loop'
	       epsilon=1.d-8
	       precision=1.d-12
	       coeff=4.d0*pi*radius**2*5.67d-5*e2phi(imax)/lsol
	       lhs=nlum(imax-1)+fk(imax-1)+fj(imax-1)*ntemp(imax)
	       ntp=ntemp(imax)
	       tp0_keep=ntp
7654   	       tp0=ntp
	       teff0=fteff(tp0/ephi(imax),ifteff,etabound,bf_r(imax),
     1           istep,time,ts1,ts2,z_ion(imax),a_ion(imax),rrho(imax),
     2             debug)
       		if(debug.eq.-50.) print *,'Tb0, Te0 =',tp0,teff0
	       tp1=(1.d0+epsilon)*tp0
	       teff1=fteff(tp1/ephi(imax),ifteff,etabound,bf_r(imax),
     1            istep,time,ts1,ts2,z_ion(imax),a_ion(imax),rrho(imax),
     2             debug)
	       if(debug.eq.-50.) print *,'Tb1, Te1 =',tp1,teff1
               derivative=coeff*(teff1**4-teff0**4)/(epsilon*tp0)
	       derivative=-fj(imax-1)-derivative
	       if(debug.eq.-50.) print *,'Derivative =',derivative
	       function=lhs-fj(imax-1)*tp0-coeff*teff0**4
	       if(debug.eq.-50.) print *,'Function =',function
	       ntp=tp0-function/derivative
	       if(debug.eq.-50.) print *,'Del(Tp)/Tp =',abs(tp0-ntp)/tp0
	       if(debug.eq.-50.) print *,'------> New Tb =',ntp
	       if ((ntp.le.0.).or.(ntp.gt.1.e12)) then ! In case the method diverges
                                 ! restart iterations with shorter time step
        	tcut=dsqrt(scale_dt0)
        	if (time.le.1.e5) tcut=dsqrt(scale_dt1)
        	dtime=dtime/tcut
        	goto 2345
      		end if
       		if(abs(tp0-ntp)/tp0.gt.precision)goto 7654
	 else 
       		if(debug.eq.-50.) print *,'Lgrav loop'
cc No gravitational energy included
	       epsilon=1.d-8
	       precision=1.d-12
	       coeff=4.d0*pi*radius**2*5.67d-5*e2phi(imax)/lsol
	       lhs=nlum(imax-1)+fk(imax-1)+fj(imax-1)*ntemp(imax)
	       ntp=ntemp(imax)
c Initial guess of tp0.
		if (ifteff.ne.5) print *,'Problem with BC'
	       teff_ini=(Lgrav/coeff)**0.25
	       ntp=(teff_ini**4.*1e4**2*4.93d16*1**4/
     1             1**3/gs14)**(2./13.)
       		if(debug.eq.-50.) print *,'Tb0 ini, Teff ini =',ntp,teff_ini
cc
	       tp0_keep=ntp
7655   	       tp0=ntp
	       teff0=fteff(tp0/ephi(imax),ifteff,etabound,bf_r(imax),
     1           istep,time,ts1,ts2,z_ion(imax),a_ion(imax),rrho(imax),
     2             debug)
       		if(debug.eq.-50.) print *,'Tb0, Te0 =',tp0,teff0
	       tp1=(1.d0+epsilon)*tp0
	       teff1=fteff(tp1/ephi(imax),ifteff,etabound,bf_r(imax),
     1           istep,time,ts1,ts2,z_ion(imax),a_ion(imax),rrho(imax),
     2             debug)
	       if(debug.eq.-50.) print *,'Tb1, Te1 =',tp1,teff1
               derivative=coeff*(teff1**4-teff0**4)/(epsilon*tp0)
	       derivative=-fj(imax-1)-derivative
	       if(debug.eq.-50.) print *,'Derivative =',derivative
	       function=lhs-fj(imax-1)*tp0-coeff*teff0**4
	       function=lhs-fj(imax-1)*tp0-
     1   (coeff*teff0**4-Lgrav)
	       if(debug.eq.-50.) print *,'Function =',function
	       ntp=tp0-function/derivative
	       if(debug.eq.-50.) print *,'Del(Tp)/Tp =',abs(tp0-ntp)/tp0
	       if(debug.eq.-50.) print *,'------> New Tb =',ntp
	       if ((ntp.le.0.).or.(ntp.gt.1.e12)) then ! In case the method diverges
                                 ! restart iterations with shorter time step
        	tcut=dsqrt(scale_dt0)
        	if (time.le.1.e5) tcut=dsqrt(scale_dt1)
        	dtime=dtime/tcut
        	goto 2345
      		end if
       		if(abs(tp0-ntp)/tp0.gt.precision)goto 7655
	endif
     	 else
       ntp=tb_acc0            ! Fixed T_b for accretion
      end if


c *********************************************************************
c ****** Get ntemp & nlum *********************************************
c *********************************************************************

      if (debug.ge.1.) print *,'Getting NTemp & NLum'
      delt(imax)=ntp-ntemp(imax)
      do i=imax-2,1,-2
       dell(i+1)=fk(i+1)-fj(i+1)*delt(i+2)
       dell(i+1)=sign(min(2.d3*abs(nlum(i+1)),abs(dell(i+1))),
     1                                                dell(i+1))
       delt(i)=fk(i)-fj(i)*dell(i+1)
       delt(i) = sign(min(.5d0*ntemp(i),abs(delt(i))),delt(i))
      end do
c	write(91,*)'1044',1013,delt(1013),imax
      dell(0)=0.d0     ! This is the inner boundary condition !

c Check for matrix inversion: *****************************************
c      if (debug.ge.1.) print *,'Checking Matrix Inversion'
c      max_equ_t=0.d0
c      max_equ_l=0.d0
c      do i=imax-2,1,-2
c       equ_t=ff( i )+
c     1       fa( i )*dell(i+1)+fb( i )*dell(i-1)+fc( i )*delt( i )	
c       equ_t=abs(equ_t)/(abs(ntemp(i)+delt(i))+1.d-8)*dtime
c       equ_l=ff(i+1)+
c     1       fa(i+1)*delt(i+2)+fb(i+1)*delt( i )+fc(i+1)*dell(i+1)
c       equ_l=abs(equ_l)/(abs(nlum(i+1)+dell(i+1)+1.d-8))
c       if (equ_t.gt.max_equ_t) max_equ_t=equ_t
c       if (equ_l.gt.max_equ_l) max_equ_l=equ_l
c      end do
c *********************************************************************

      do i=0,imax-1,2
       nlum(i)=nlum(i)+dell(i)
       ntemp(i+1)=ntemp(i+1)+delt(i+1)
      end do
c	write(91,*)'1066',1013,ntemp(1013),delt(1013)
      dtemp(0)=0.d0
      do i=2,imax-1
       dtemp(i)=(ntemp(i+1)-ntemp(i-1))/(debar(i)+debar(i+1))
      end do
      do i=1,imax-2,2
       dlum(i)=(nlum(i+1)-nlum(i-1))/(debar(i)+debar(i+1))
      end do
c look +++++++++++++++++++++++++++++++++++++++++++++
      dlum(imax)=0.d0
c      dlum(imax)=dlum(imax-2)
c ++++++++++++++++++++++++++++++++++++++++++++++++++

c BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
c      INCLUDE 'Bfield/Bfield_3.inc.f'
c BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

c *********************************************************************
c ***** Analyze the results to see if it has converged ****************
c *********************************************************************

      if (debug.ge.1.) print *,'Analyzing Results'
      ratiot=0.d0
      ratiol=0.d0
      do i=1,imax-2,2
       ratl=abs(dell(i+1))/(abs(nlum(i+1))+1.d-12)
       if(ratl.gt.ratiol)then
        ratiol=ratl
        iratl=i+1
       end if
       ratt=abs(delt(i)/(ntemp(i)+1.d-30))
       if (ratt.gt.ratiot)then
        ratiot=ratt
        iratt=i
       end if
      end do

c BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
c This was used when doing dipolar magnetic field evolution.
c Not used any more.
      ratios=0.d0
      irats=1
c      INCLUDE 'Bfield/Bfield_4.inc.f'
c BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

c ---------------------------------------------------------------------
      if (pscreen.eq.3) then
       if (ifield.eq.2) then
        print '(3x,3(a10,1p1e10.3,a3,i4,a1))',
     1                                     'dT/T=',ratiot,'(',iratt,')',
     2                                     'dL/L=',ratiol,'(',iratl,')',
     3                                     'dS/S=',ratios,'(',irats,')'
       else
        print '(3x,2(a10,1p1e10.3,a3,i4,a1))',
     1                                     'dT/T=',ratiot,'(',iratt,')',
     2                                     'dL/L=',ratiol,'(',iratl,')'
       end if
      end if
c ---------------------------------------------------------------------


c *********************************************************************
c     Decide if converged or not:
      if ((ratiot.lt.mratt).and.(ratiol.lt.mratl).and.(ratios.lt.mrats))
     x then
        continue        ! Converged ! continue to next time step
       else 
        goto 2000       ! Not converged ! Go back for another iteration
       end if
c *********************************************************************

	if ((ilgrav.eq.0).or.(lgrav.eq.0.)) then
      luminosity=nlum(imax-1)/ephi(imax-1)**2
      sign_l=abs(nlum(imax-1))/nlum(imax-1)
	else
      luminosity=nlum(imax-1)/ephi(imax-1)**2+Lgrav/ephi(imax-1)**2
      sign_l=abs(luminosity)/luminosity
	endif	
      teffective=sign_l*
     2   (abs(luminosity)/(4.d0*pi*radius**2*5.67d-5))**.25d0*
     3   lsol**.25d0*ephi(imax-1)
c In case of accretion with fixed outer Tb, lum at the surface can become negative: 
c the surface is injecting heat into the star (from the surface nuclear burning). 

c BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
c      INCLUDE 'Bfield/Bfield_5.inc.f'
c BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

c *********************************************************************
c ******* PREPARATION TO CALCULATE THE NEW TIME STEP  *****************
c *********************************************************************
c
c This is a delicate part, based on experience and many trials and errors.
c It works pretty well, so avoid changing it !
c
c PHILOSOPHY OF TIME STEP CONTROL:
c
c (Time of step just finished is "time+dtime", not just time !)
c The new "dtime" will be "scale_dt*dtime" with "scale_dt" calculated below.
c Allows for 2 different "scale_dt": at early time, while relaxing from initial
c conditions, accuracy is not important and one can allow for larger timestep:
c "scale_dt0" and "scale_dt1" are read from the file 
c NUM_PARAM.dat in NSCool_READ.inc.f
c and are the maximum allowed relative increase in "dtime"
c *********************************************************************

      if (debug.ge.1.) print *,'Calculating New Time Step'

      scale_dt=scale_dt0
      if (time.le.1.e5) scale_dt=scale_dt1
	test=0.0005
	frac=1.1
c      if ( i_heat_h.eq.1 ) then
c      if ( ( time.gt.t_acc2*(1-test) ) .and. 
c     2     ( time.lt.t_acc2*(1+test) )  ) then
c	scale_dt=scale_dt1/frac
c	write(90,*)'frac',scale_dt,scale_dt1
c	endif
c	endif

c *********************************************************************
c TEMP variation: "max_dtemp" is the max. relative variation of T in the
c                 star and "icht" the zone where "max_dtemp" is obtained
      icht=0
      max_dtemp=0.d0
      do i=1,imax-2,2
       mdt=abs(temp(i)-ntemp(i))/ntemp(i)
       if (mdt.gt.max_dtemp)then
        max_dtemp=mdt
        icht=i
       end if
      end do
c *********************************************************************
c Check if "max_dtemp" is not too large. If "1+max_dtemp" exceeds "tvar",
c then "scale_dt" is reduced correspondingly:
c ("tvar" read from the file NUM_PARAM.dat in NSCool_READ.inc.f)
      chtemp=0.d0
      if((tvar-1.d0).lt.max_dtemp)then
       scale_dt=scale_dt*(tvar-1.)/max_dtemp
       if (time.lt.1.d5) then
        scale_dt=min(scale_dt,scale_dt1)
       else
        scale_dt=min(scale_dt,scale_dt0)
       end if
       chtemp=1.d0     ! this means the time-step is controlled b
      end if           ! a too large max_dtemp
c BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
       chstoke=0.d0
c       INCLUDE 'Bfield/Bfield_6.inc.f'
c BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
c *********************************************************************
c If "scale_dt" has been reduced too much, i.e. "max_dtemp" too large,
c then the present time-step is recalculated with a smaller "dtime"
c unless one is still at the first few time-steps:
c ("repeat" and "istart" are read from the file
c   NUM_PARAM.dat in NSCool_READ.inc.f)
      if ( (scale_dt.lt.repeat).and.(istep.gt.istart) ) then
       dtime=repeat*dtime 
c ---------------------------------------------------------------------
       if (pscreen.ge.2) then
        write(6,*)
        print '(1a40)','dtime too large, do it again'
        print '(1a13,1p1e12.3)','time=',(time+dtime)/year
        print '(1a13,1p1e12.3)','dtime=',dtime/year
        print '(1a13,1p1e12.3)','dtime ratio=',dtime/odtime
        write(6,*)
       end if
c ---------------------------------------------------------------------
       goto 2345
      end if
c In case convergence is reached in too many trials, scale_dt is reduced
      if ((itrial.gt.itrial_opt).and.(istep.gt.istart)) then
       chtrial=1.d0     ! this means the time-step is controlled by needing too many iterations
       olddt=scale_dt
       if (time.lt.1.d5) then
        scale_dt=scale_dt/
     x           scale_dt1**(float(itrial-itrial_opt)/2.d0)
        scale_dt=min(scale_dt,scale_dt1)
       else
        scale_dt=scale_dt/
     x           scale_dt0**((1.d0+float(itrial-itrial_opt))/2.d0)
        scale_dt=min(scale_dt,scale_dt0)
       end if
      else
       chtrial=0.d0
      end if
c Before setting the next time and time-step, stuff are printed out and updated:

c *********************************************************************
c ***** End of iterations
c *********************************************************************

      do 171 i=1,imax,2
       otemp(i)=temp(i)
       temp(i)=ntemp(i)
       orrho(i)=rrho(i)
       orad(i)=rad(i)       
       obar(i)=bar(i)
171   continue
       do 172 i=2,imax-1,2
       olum(i)=lum(i)
       lum(i)=nlum(i)
       orrho(i)=rrho(i)
       orad(i)=rad(i)       
       obar(i)=bar(i)
172   continue

c BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
c      INCLUDE 'Bfield/Bfield_7.inc.f'
c BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

c *********************************************************************
c Spin down: **********************************************************
c *********************************************************************
      call spin_down(period,dtime,beta_rot,bfield2(imax),nperiod)
c      dp_dt=(nperiod-period)/dtime
c      p_av=(nperiod+period)/2.d0            ! average P
c      omega=2.d0*pi/p_av**2
c      domega_dt=-2.d0*pi*dp_dt/p_av**2
c      energy_loss=m_i*omega*domega_dt       ! Spin-down power
      period=nperiod
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Stuff below, till the next +++++ line is only informative and not
c used in the calculations.
c
c *********************************************************************
c Crystalization density: *********************************************
c *********************************************************************
      if (debug.ge.1.) print *,'Calculating Crystalization Densities'
      do i=imax,icore+2,-2
       gamma(i)=2.273d5*z_ion(i)**2*(rrho(i)/a_cell(i))**(1.d0/3.d0)/
     1          (temp(i)/ephi(i))
       if (gamma(i).lt.gammaliq) then
        cryst(i)=-1.
       else if ( (gamma(i).ge.gammaliq  ) .and.
     1           (gamma(i).le.gammacryst)       ) then
        cryst(i)=0.d0
       else
        cryst(i)=1.d0
       end if
      end do
      do i=icore,0,-2
       gamma(i)=1.d0
      end do
      iliq=0
      icryst=0
      do i=imax,icore+2,-2
       if (cryst(i).eq.0.d0) then
        iliq=i
        goto 600
       end if
      end do
600   continue
      do i=imax,icore+2,-2
       if (cryst(i).eq.1.d0) then
        icryst=i
        goto 700
       end if
      end do
700   continue

c *********************************************************************
c ***** Calculate the neutrino luminosity and heating: ****************
c *********************************************************************
       if (debug.ge.1.) print *,
     1     'Calculating Total Neutrino Luminosity and Heating'
       lnu =0.d0
       lnu0=0.d0
       lh  =0.d0
       lh0 =0.d0
       do i=1,imax,2
        lnu =lnu +e2phi(i)* qnu(i)*(dvol(i)+dvol(i+1))
        lnu0=lnu0+ephi(i) * qnu(i)*(dvol(i)+dvol(i+1))   ! without energy red-shift
        lh  =lh  +e2phi(i)*heat(i)*(dvol(i)+dvol(i+1))
        lh0 =lh0 +ephi(i) *heat(i)*(dvol(i)+dvol(i+1))   ! without energy red-shift
       end do
       lnu=lnu/lsol
       lh =lh/lsol
       htot=htot+lh*dtime

c ***** CALCULATE THE INTEGRATED NEUTRINO LUMINOSITIES: ***************
c     Note: lnu_tot, calculated from qnu(i), is the garanteed total
c     neutrino luminosity. The other ones are only informative.
      lmurca_nucl=0.0d0
      lbrem_nucl =0.0d0
      lplasma    =0.0d0
      lnpb       =0.0d0
      lpbf_n1S0  =0.0d0
      lpbf_n3P2  =0.0d0
      lpbf_p1S0  =0.0d0
      do i=1,imax,2
       e2p=e2phi(i)
       lmurca_nucl=lmurca_nucl+qmurca_nucl(i)*(dvol(i)+dvol(i+1))*e2p
       lbrem_nucl =lbrem_nucl +qbrem_nucl(i) *(dvol(i)+dvol(i+1))*e2p
       lplasma    =lplasma    +qplasma(i)    *(dvol(i)+dvol(i+1))*e2p
       lnpb       =lnpb       +qnpb(i)       *(dvol(i)+dvol(i+1))*e2p
       lpbf_n1S0  =lpbf_n1S0  +qpbf_n1S0(i)  *(dvol(i)+dvol(i+1))*e2p
       lpbf_n3P2  =lpbf_n3P2  +qpbf_n3P2(i)  *(dvol(i)+dvol(i+1))*e2p
       lpbf_p1S0  =lpbf_p1S0  +qpbf_p1S0(i)  *(dvol(i)+dvol(i+1))*e2p
      end do

c ***** CALCULATE THE INTEGRATED SPECIFIC HEATS: **********************
c     cv_tot_all, calculated from cv(i), is the garanteed total
c     specific heat. The other ones are only informative.
      cv_core=0.d0
      cv_crust=0.d0
      cv_tot_all=0.d0
      cv_tot_ion=0.d0
      cv_tot_neu=0.d0
      cv_tot_pro=0.d0
      cv_tot_ele=0.d0
      cv_tot_muo=0.d0
      cv_tot_lam=0.d0
      cv_tot_sim=0.d0
      cv_tot_si0=0.d0
      cv_tot_sip=0.d0
      cv_tot_qrk=0.d0
      cv_phot=0.d0
      do i=1,imax,2
       cv_tot_all=cv_tot_all+cv(i)   *(dvol(i)+dvol(i+1))
       cv_tot_ion=cv_tot_ion+cv_ion(i)*(dvol(i)+dvol(i+1))
       cv_tot_neu=cv_tot_neu+cv_n(i)  *(dvol(i)+dvol(i+1))
       cv_tot_pro=cv_tot_pro+cv_p(i)  *(dvol(i)+dvol(i+1))
       cv_tot_ele=cv_tot_ele+cv_e(i)  *(dvol(i)+dvol(i+1))
       cv_tot_muo=cv_tot_muo+cv_m(i)  *(dvol(i)+dvol(i+1))
       cv_tot_lam=cv_tot_lam+cv_l(i)  *(dvol(i)+dvol(i+1))
       cv_tot_sim=cv_tot_sim+cv_sm(i) *(dvol(i)+dvol(i+1))
       cv_tot_si0=cv_tot_si0+cv_s0(i) *(dvol(i)+dvol(i+1))
       cv_tot_sip=cv_tot_sip+cv_sp(i) *(dvol(i)+dvol(i+1))
       cv_tot_qrk=cv_tot_qrk+cv_q(i)  *(dvol(i)+dvol(i+1))
       cv_phot= 4.*7.56e-15*(ntemp(i)/ephi(i))**3  *(dvol(i)+dvol(i+1))
      end do
      do i=1,icore+2
       cv_core=cv_core+cv(i)*(dvol(i)+dvol(i+1))
      end do
      do i=icore+2,imax
       cv_crust=cv_crust+cv(i)*(dvol(i)+dvol(i+1))
      end do

c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c *********************************************************************
c ******************** print out results ******************************
c *********************************************************************

cc
cc print temperature at rho=1e10 g.cm-3 cc
cc
c	rho10=1e10
c	j=imax
c	do while (rrho(j).le.rho10) 
c	j=j-2
c	enddo
c	i10=j
c          x=(dlog(rrho(i10))-dlog(rho10))/
c     1      (dlog(rrho(i10))-dlog(rrho(i10+2)))
c          y=(dlog(rho10)-dlog(rrho(i10+2)))/
c     1      (dlog(rrho(i10))-dlog(rrho(i10+2)))
c          logntemp10=y*dlog(ntemp(i10))+x*dlog(ntemp(i10+2))
c	 ntemp10=dexp(logntemp10)
c 603    format(i8,1p1g16.9,5x,1p7g16.9)
c	write(90,603) istep,(time+dtime)/year,
c     1               rrho(i10),rho10,rrho(i10+2),
c     2               ntemp(i10),ntemp10,ntemp(i10+2),
c     3               lum(imax-1)*lsol

      INCLUDE 'NSCool_PRINT.inc.f'

c ----------------------------------------------------------------------
      if (pscreen.ge.2) then
       print *
       print '(1a21,1i3,1a9)','Iteration finished:',
     1                    itrial,'trials.'
       print *

       print '(a8,1p1e12.3,a30,1p1e12.3,a10,i5)','Teff =',teffective,
     1                   'max delt/ntemp=',max_dtemp,'at i=',icht
       print *
       if ((iliq.eq.0).and.(icryst.eq.0)) then
        print *,'   Whole crust liquid'
       else if (icryst.eq.imax) then
        print *,'   Whole crust solid'
       else if ((iliq.eq.0).and.(icryst.ne.0)) then
        print '(22x,a10,1p1e12.3,a3,i3,a1)',
     1    'Rho_s=',rrho(icryst),'(',icryst,')'
       else if ((iliq.ne.0).and.(icryst.eq.0)) then
        print '(a10,1p1e12.3,a3,i3,a1)',
     1    'Rho_l=',rrho(iliq),'(',iliq,')'
       else
        print '(2(a10,1p1e12.3,a3,i3,a1))',
     1    'Rho_l=',rrho(iliq),'(',iliq,')',
     2    'Rho_s=',rrho(icryst),'(',icryst,')'
       end if
       print *
       print '(5x,3(a10,1p1e10.3))',
     1   'L(ph)=',luminosity*lsol*e2phi(imax-1),
     2   'L(nu) =',lnu*lsol,
     3   'L(H)  =',lh*lsol
       lh_red=lh
       print *
       if ((i_acc.eq.1).or.(i_acc.eq.2)) then
        print '(a12,1p1e10.2)','M_dot =',m_dot*3.15576d7/2.d33
        print '(a12,1i10,a12,0p1f8.2)',
     1       'Cycle number',icycle,'Phase =',t_burst/t_acc1
        print *
       end if
      end if
c ----------------------------------------------------------------------
      if (pscreen.eq.1) then
       if (i_acc.eq.0) then
         print '(i5,1p1e12.3,3x,1p1e12.3,3x,1p3e12.3)',
     1        istep,(time+dtime)/year,teffective,
     2        nlum(imax-1)*lsol,lnu*lsol,lh*lsol
       else
        if (icycle.ge.1) then
         if (t_burst.le.t_acc2) then
          t_check=t_burst/t_acc2
         else
          t_check=0.d0
         end if
        end if
       if (i_acc.ne.3) then
        print 
     x '(i7,1p1e22.15,1p1e14.6,1p3e12.3,
     x   1p1e15.3,i10,0p1f10.3,1p1e12.3,
c     x   5x,1p2e12.3)',
     x   5x,1p6e14.6)',
     1   istep,(time+dtime)/year,teffective,
     2   nlum(imax-1)*lsol,lnu*lsol,lh*lsol,
     3   odtime,icycle,t_check,m_dot*3.15576d7/2.d33,
     4   ntemp(imax)/ephi(imax),
     5   ntemp(idrip)/ephi(idrip),ntemp(1)/ephi(1),
     6   ephi(imax),ephi(idrip),ephi(1)
c         write(90,'(i7,1p1e22.15,1p1e14.6,i10,1p1e12.3,0p1f8.2,1i3)')
c     x   istep,(time+dtime)/year,teffective,icycle,
c     x   m_dot*3.15576d7/2.d33,t_burst/t_acc1,itrial

      end if

       end if

c          write(91,*)(time+dtime)/year,m_dot*3.15576d7/2.d33
      end if
 



c ----------------------------------------------------------------------------

c *********************************************************************
c ****************  CALCULATE THE NEW TIME STEP  **********************
c *********************************************************************
      time=time+dtime
      odtime=dtime
      dtime=min(scale_dt*dtime,dtlimit)

c *********************************************************************
c For accretion scenarios: time step must moreover be shortened 
c dramatically when a new outburst is approaching (to make sure it is  
c much shorter than the outburst duration, or rise time, or any relevant 
c time scale which has to be resolved by the code):

c Transient FRED ("Fast rise and exponential decay") *******************
      if (i_acc.eq.1) then
       time_step_cut=100.d0
       if (((time+3.*dtime).ge.t_acc0).and.(icycle.eq.0)) then
        timeleft=t_acc0-time
        dtime=max(timeleft/3.d0,t_acc2/time_step_cut)
       end if
       if (time.gt.t_acc0) then
        if (delt_acc/t_acc2.le.10.d0) then
         scale_dt0=1.05d0
        else
         scale_dt0=1.2d0
        end if
        t_next=t_acc0+float(icycle+1)*t_acc1
        if ((time+3.*dtime).ge.t_next) then
         timeleft=t_next-time
         dtime=max(timeleft/3.d0,t_acc2/time_step_cut)
        end if
       end if
       if (time.ge.t_acc0) dtlimit=t_acc1/20.d0
      end if

c Transient STEP ******************************************************
      if (i_acc.eq.2) then
       if (((time+2.*dtime).ge.t_acc0).and.(time.lt.t_acc0)) then
        timeleft=t_acc0-time
        dtime=max(timeleft/3.d0,time_step_min)
        day=86400.d0
c        print '(a30)','Approaching first burst:'
c        print '(2(a20,0p1f18.8))',
c     1        'Time left =',timeleft/day,'DTime =',dtime/day
       end if
       if (time.gt.t_acc0) then
        t_next=t_acc0 + float(icycle-1)*t_acc1 + t_acc1
        t_end =t_acc0 + float(icycle-1)*t_acc1 + t_acc2
c       Slow down if approaching next burst
        if ((time+dtime.gt.t_end).and.
     1      (time+2.*dtime.ge.t_next)) then
         timeleft=t_next-time
         dtime=max(timeleft/3.d0,time_step_min)
c        print '(a30)','Approaching next burst:'
c        print '(2(a20,0p1f18.8))',
c     1        'Time left =',timeleft/day,'DTime =',dtime/day
c       Slow down if approaching end of burst
c        else if ((time+dtime.le.t_end).and.
c     1           (time+2.*dtime.ge.t_end)) then
        else if ((time.le.t_end).and.
     1           (time+2.*dtime.ge.t_end)) then
c        else if (time+2.*dtime.ge.t_end) then
         timeleft=t_end-time
         dtime=max(timeleft/3.d0,time_step_min)
c        print '(a30)','Approaching end of burst:'
c        print '(2(a20,0p1f18.8))',
c     1        'Time left =',timeleft/day,'DTime =',dtime/day
        end if
       end if
       if (time.ge.t_acc0) dtlimit=t_acc1/20.d0
      end if

cc Multiple STEP ******************************************************
      if (i_acc.eq.4) then
	t_acc1=1e30
       if (((time+2.*dtime).ge.t_acc0).and.(time.lt.t_acc0)) then
        timeleft=t_acc0-time
        dtime=max(timeleft/3.d0,time_step_min)
        day=86400.d0
       end if
       if (time.gt.t_acc0) then
        t_next=t_acc0 + float(icycle-1)*t_acc1 + t_acc1
        t_end =t_acc0 + float(icycle-1)*t_acc1 + t_acc2
        if ((time+dtime.gt.t_end).and.
     1      (time+2.*dtime.ge.t_next)) then
         timeleft=t_next-time
         dtime=max(timeleft/3.d0,time_step_min)
        else if ((time.le.t_end).and.
     1           (time+2.*dtime.ge.t_end)) then
         timeleft=t_end-time
         dtime=max(timeleft/3.d0,time_step_min)
        end if
       end if
       if (time.ge.t_acc0) dtlimit=t_acc1/20.d0
      end if



c Heat Deposition *****************************************************
      if (i_heat_deposit.eq.1) then
       t_slow=max(dtime,1000.*del_t_dep)
       if (abs(t_dep-(time+dtime)).lt.t_slow) then
        timeleft=abs(t_dep-(time+dtime))
        dtime=max(timeleft/30.d0,del_t_dep/100.d0)
       end if
      end if
c *********************************************************************
      if (time/year.ge.timemax) goto 9998
      if ((sign_l*teffective).lt.tempmin) goto 9998
ccc Steady state of transiently accreting NS 08/2010
      if ( (i_acc.eq.3) .and. 
     1   ( abs(otemp(imax)-temp(imax)).le.5e-10 ) ) then
c	write(*,*)time,otemp(imax),temp(imax),otemp(imax)-temp(imax)
	goto 9998
	endif
ccc
c *********************************************************************
9999  continue     ! This is the `end do' for the main time integration
c *********************************************************************

c *********************************************************************
9998  continue     ! Jump here if time > timemax
c *********************************************************************

c Close the output file for the present model: ************************
      close(unit=10,status='keep')
      close(unit=19,status='keep')

c Go back to the beginning: do the next model listed in the "Input file"
      goto 1234

c *********************************************************************
 9997 continue    ! This is the real end of it !
c *********************************************************************

      close(unit=15,status='keep')
      close(unit=29,status='keep')
      print *,'Done !'

      end
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
