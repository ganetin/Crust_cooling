c--------------------------------------------------------
      SUBROUTINE T_NP(y,g,no)
c--------------------------------------------------------
c          Finite temperature
c--------------------------------------------------------
      implicit none
      include 'const.inc'
      integer i,no
      real*8 y(no),g(no)
      real*8 epsk,presk,nbar,eps1,pres1
      real*8 epsb,presb,entb,epsl,presl,entl,entnu,ente

      real*8 sig, ome, rho, us, duds
      real*8 xnu(6),xmu(6),efm(6),dnum(6)
      real*8 pre(6),ed(6),xns(6),entr(6)

      real*8 nb,ye,yp,yn,yl,ynu,t,s,p,e,muh,mun,mue,ym,mns,mup,mum
      common /EOS_IO/ nb,ye,yp,yn,yl,ynu,t,s,p,e,muh,mun,mue,ym,mns


C
C--------------------------------------------------------
C
C    1 -> e neutrino    (or sigma field)
C    2 -> mu-           (or omega field)
C    3 -> mu+           (or rho   field)
C    4 -> electron      (or baryon number)
C    5 -> neutron       (or charge)
C    6 -> proton        (or electron number)
C
C--------------------------------------------------------
C
      sig = y(1)
      ome = y(2)
      rho = y(3)
      mun = y(4) 
      muh = y(5)
      mue = y(6)
cc      mue = muh
      mum = muh
      mup = mun-muh
c      t   = y(7)

c---------------------------
c    chemical potential as det. by beta-equilibrium
c---------------------------

      xnu(1) = mue-muh
      xnu(2) = -mue
      xnu(3) = mum
      xnu(4) = mue
      xnu(5) = mun - gw*ome + 0.5d0*gr*rho
      xnu(6) = mup - gw*ome - 0.5d0*gr*rho
                                    
c---------------------------
c       effective masses
c---------------------------

      mns = mn*(1.d0 - gs*sig/mn)
      efm(5) = mns
      efm(6) = mns

c---------------------------
c       lepton thermodinamics
c---------------------------

      call lepeos(t,xnu(1),dnum(1),ed(1),pre(1),entr(1))
      call JELdpe(t,mel,xnu(2),dnum(2),pre(2),ed(2),xns(2),entr(2))
      call JELdpe(t,mmu,xnu(3),dnum(3),pre(3),ed(3),xns(3),entr(3))
      call JELdpe(t,mel,xnu(4),dnum(4),pre(4),ed(4),xns(4),entr(4))

      presl = 0.5*pre(1) + pre(4) + pre(3) + pre(2)
      epsl  = 0.5*ed(1)  + ed(4)  + ed(3) + ed(2)
      entl  = 0.5*entr(1)+ entr(4) + entr(3) + entr(2)

c---------------------------
c       baryon thermodynamics
c---------------------------
      do i=5,6
        call JELdpe(t,efm(i),xnu(i),dnum(i),pre(i),ed(i),xns(i),entr(i))
      enddo

c---------------------------
c      Thermodinamical variables 
c---------------------------

      epsk = 0.d0
      presk = 0.d0
      nbar = 0.d0
      entb = 0.d0
      do i=5,6
        epsk  = epsk  + ed(i)
        presk = presk + pre(i)
        nbar  = nbar  + dnum(i)
        entb  = entb  + entr(i)
      enddo

      us = b/3.d0*mn*(gs*sig)**3 + c/4.d0*(gs*sig)**4
      duds = gs*( b*mn + c*(gs*sig) )*(gs*sig)**2

      eps1 = us + 0.5d0*(ms*sig)**2 + 0.5d0*(mw*ome)**2 
     &        + 0.5d0*(mr*rho)**2 
      pres1 = - us - 0.5d0*(ms*sig)**2 + 0.5d0*(mw*ome)**2 
     &        + 0.5d0*(mr*rho)**2 

      epsb  = epsk  + eps1 
      presb = presk + pres1 
     
      ye  = dnum(4)/nbar
      ynu = 0.5d0*dnum(1)/nbar
      yp  = dnum(6)/nbar
      yn  = dnum(5)/nbar
cc
      ym  = dnum(3)/nbar
cc
c---------------------------
c      set of equations to be solved
c---------------------------

        g(1) = ms*ms*sig + duds   
     &       - gs*(xns(5)+xns(6))

        g(2) = mw*mw*ome 
     &       - gw*(dnum(5)+dnum(6))

        g(3) = mr*mr*rho 
     &       - gr*0.5d0*(dnum(6) - dnum(5))

        g(4) = 1.d0 - nbar/nb

        g(5) = ( dnum(6)+dnum(2)
     &       -  (dnum(4)+dnum(3)) )/nb

c        g(6) = ye  - dnum(4)/nbar
        g(6) = yl  - ye-ynu!-ym

c        g(7) =  s - (entb+entl)/nb

      p = presb+pre(4)
      e = epsb+ed(4)
      s = entb/nb+entl/nb
cccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccc
      mns =(1.d0 - gs*sig/mn)/mn
cccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccc

      return
      end


