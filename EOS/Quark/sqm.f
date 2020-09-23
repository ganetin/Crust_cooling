c CODE FROM J. L. ZDUNIK
c JUNE 2017
c......
c         program sqm.f (Strange Quark Matter)
c
c      this program is a copy of a code qmataut (automatically calculating
c       quark matter) written by JLZ in 1989.
c      It is slightly changed - calculates also the function 
c          H=log((rho+P)/n)-log((rho0+P0)/n0)   (P0=0)
c
c		23.07.1998
c
c...
c            qmataut.for
c     Program oblicza sklad materii kwarkowej w rownowadze
c     silnej i slabej zlozonej z kwarkow u,d,s.
c     Wyznacza punkt P=0
c     Tworzy rownanie stanu na zbiorze file2 w postaci
c         i, N, RO, P
c     Na zbiorze file4 parametry materii kwarkowej
c     Jednostki:
c               P - cisnienie/B
c               ro - gestosc energii*c**2/B
c               u - pot. chem. czastek/(hc**3*B)**.25
c               n - gestosc czastek*(hc/B)**.75
c               om - funkcja termodyn. omega/B
c     Parametry materii:
c               q - ms*c**2/(hc**3*B)**.25
c               alfac - oddzialywanie
c               B - stala bagu (zawarta w skalowaniu)
c
c          25.09.1989
c...
      implicit double precision (a-h,n-z)
      double precision ms
      character*30 file1, file2, file4
      common /par/ q,alfac,ren
      common /con/ pi
      common /outnew/ istnew,iernew
      data hc, xkbmev /197.32858d0, 8.61735d-11/
      data chp /1.6021892d33/
      data cl /2.9979d10/
      chro=chp/cl/cl
      pi=3.141592653589793238d0
c...wspolczynniki zmian jednostek na cgs (P,ro) fm**-3 (n) i MeV (u)
c       dla B0=60 MeV*fm**-3
      B0=60.
      cp=B0*chp
      cro=B0*chro
      cu=(hc**3*B0)**0.25d0
      cn=(B0/hc)**0.75d0
      write (*,6005)
      read (*,6006) file2
      write (*,6007)
      read (*,6006) file4
 6005 format ('  Nazwa zbioru wynikowego pelnego > ')
 6007 format ('  Nazwa zbioru wynikowego z eosem > ')
 6006 format (a30)     
      open (2,file=file2,status='new')
      open (4,file=file4,status='new')
c...   wczytanie parametrow materii kwarkowej
c      write (*,*) 'podaj q '
      write (*,*) ' strange quark mass m_s (MeV) = '
      read (*,*) ms
      q=ms/(hc**3*B0)**.25
      write (*,*) ' alfa_c = '
      read (*,*) alfac
c...   parametr ren=313./B**.25  dla B0=60.
      ren=2.136062d0
c...   minimalne i maksymalne cisnienie w EOSie
      pmin=1.d-10
      pmax=7.
c
c...   wartosci qp, dp okreslaja zmiany cisnien w rownaniu stanu
c      qp - wzgledna zmiana cisnienia w tabeli eosu
c      dp - max bezwzgledna zmiana
      qp=1.d0
      dp=0.1d0
c      write (2,2009) q,alfac,ren,B0
      write (2,2019) ms,q,alfac,ren,B0
      write (2,2012)
c      write (2,2010)
c      write (4,4010) q,alfac,ren
      write (4,4011) ms,q,alfac,ren
      write (4,4012)
 2000 format (4x,'I',6x,'nu',10x,'nd',10x,'ns',10x,'ne',10x,' N',
     x10x,'RO',10x,' P',10x,'MIU')
 2012 format (4x,'i',5x,'n_u',10x,'n_d',10x,'n_s',10x,'n_e',10x,'n_b',
     x10x,'rho',12x,'P',11x,'mu',12X,'H')
 2009 format (' Quark matter u,d,s  q =',g12.5,' alfa =',g12.5,
     x' ren =',g14.7,5x,'nuclear units for bag B0 =',g12.5)
 2019 format (' Quark matter u,d,s  ms =',g12.5, '(q=',f7.5,');'
     x,' alfa=',f7.5,
     x'; ren=',f9.6,5x,'m_s and nuclear units for B0 =',f5.1)
 2010 format (4x,'I',6x,'N',8x,'nu/3n',10x,'uu',7x,'nd/3n',10x,'ud',7x,
     x'ns/3n',10x,'us',10x,'ne',10x,'ue',10x,'u')
 4010 format (' Quark matter u,d,s  q =',g12.5,' alfa =',g12.5,
     x' ren =',g14.7)
 4011 format (' Quark matter u,d,s   ms =',g12.5, '(q=',f7.5,');'
     x,' alfa=',f7.5,
     x'; ren=',f9.6)
 4012 format (4x,'i',2x,'n_b [(B/hc)^0.75]',6x,'rho [B/c^2]',11x,
     x'P [B]',19x,'H')
c...    wartosci startowe potencjalu kwarku d (s)
      u04=13.15947d0/(1.-2.d0*alfac/pi)
      u4=u04
      du4=1.
      call zercis (u4,du4,pmin)
      ud=u4**0.25d0
      call matq3 (ud,nu,nd,ns,ne,uu,us,ue,p,ro,n,u)
      H0=log10((ro+p)/n)
      i=1
      H=log10((ro+p)/n)-H0
      write (4,4001) i,n,ro,p,H
      write (*,5001) i,n,ro,p,istnew
      write (2,2002) i,nu*cn,nd*cn,ns*cn,ne*cn,n*cn,ro*cro,p*cp,u*cu,H
c      write (2,2011) i,n*cn,nu/3./n,uu*cu,nd/n/3.,ud*cu,ns/n/3.,us*cu,
c     xne*cn,ue*cu,u*cu
      u41=u4+du4
      ud1=u41**0.25d0
      call matq3 (ud1,nu,nd,ns,ne,uu,us,ue,p1,ro,n,u)
      dpdu4=(p1-p)/du4
  10  i=i+1
      du4=qp*p/dpdu4
      if (qp*p.gt.dp) du4=dp/dpdu4
      u4=u4+du4
      ud=u4**0.25d0
      pold=p
      call matq3 (ud,nu,nd,ns,ne,uu,us,ue,p,ro,n,u)
      H=log10((ro+p)/n)-H0
      dpdu4=(p-pold)/du4
      write (4,4001) i,n,ro,p,H
      write (*,5001) i,n,ro,p,istnew
      write (2,2002) i,nu*cn,nd*cn,ns*cn,ne*cn,n*cn,ro*cro,p*cp,u*cu,H
c      write (2,2011) i,n*cn,nu/3./n,uu*cu,nd/n/3.,ud*cu,ns/n/3.,us*cu,
c     xne*cn,ue*cu,u*cu
      if (p.lt.pmax) go to 10
 4001 format (i5,4(g20.10))
 5001 format (i5,3(2x,g12.5),i6)
 2001 format (i5,8(g12.5))
 2002 format (i5,9(g13.5))
 2011 format (i5,10(g12.5))
   20 close (2)
      close (4)
      stop
      end
c...
c       Procedura znajduje zero cisnienia
c       tzn. takie u4 dla ktorego 0<p<pmin
c
      subroutine zercis (u4,du4,pmin)
      implicit double precision (a-h,n-z)
      external qeq
      common /par/ q,alfac,ren
      common /con/ pi
      du4=dabs(du4)
      ud=u4**0.25d0
      call matq3 (ud,nu,nd,ns,ne,uu,us,ue,p,ro,n,u)
      if (p.gt.0.d0) du4=-du4
   10 u4=u4+du4
      ud=u4**0.25d0
      pold=p
      call matq3 (ud,nu,nd,ns,ne,uu,us,ue,p,ro,n,u)
      if (p.lt.0.d0) go to 10
      u41=u4-du4
      u42=u4
      p1=pold
      p2=p
c...    znajdowanie zera P metoda siecznych
      u4=1.d10
      i=0
  20  u4old=u4
      i=i+1
      u4=(p1*u42-p2*u41)/(p1-p2)
      ud=u4**0.25d0
      call matq3 (ud,nu,nd,ns,ne,uu,us,ue,p,ro,n,u)
      if (dabs(p).lt.pmin) go to 30
      if (p*p1.lt.0.d0) then
        p2=p
	u42=u4
	else
	p1=p
	u41=u4
	endif
      go to 20
  30  du4=dabs(u4-u4old)
      write (*,100) i
  100 format (' p0 sieczne il. krokow =',i4)
      if (p.gt.0.d0) return
      write (*,*) ' p0 zmiana znaku P'
      i=0
  40  u4=u4+du4
      ud=u4**0.25d0
      i=i+1
      pold=p
      call matq3 (ud,nu,nd,ns,ne,uu,us,ue,p,ro,n,u)
      if (p.lt.0.d0) go to 40
      if (p.lt.pmin) return
      write (*,*) ' p0 zmiana znaku P - bisekcje'
      u41=u4-du4
      p1=pold
      u42=u4
      p2=p
  50  u4=0.5d0*(u41+u42)
      ud=u4**0.25d0
      call matq3 (ud,nu,nd,ns,ne,uu,us,ue,p,ro,n,u)
      if ((p.gt.0.d0).and.(p.lt.pmin)) go to 60
      if (p.lt.0.d0) then
         u41=u4
	 p1=p
	 else
	 u42=u4
	 p2=p
	 endif
      go to 50
   60 du4=dabs(u42-u4)
      return
      end
c...
      subroutine matq3 (ud,nu,nd,ns,ne,uu,us,ue,p,ro,n,u)
      implicit double precision (a-h,n-z)
      external qeq
      common /par/ q,alfac,ren
      common /con/ pi
      common /doqeq/ qalf,nd0,ns0,ud0
      common /outnew/ istep,iernew
      pi2=pi*pi
      qalf=1.d0-2.d0*alfac/pi
      q2=q*q
      epsu=1.d-12
c...   rownowaga termodynamiczna
      us=ud
      nd=ud**3/pi2*qalf
      us2=us*us
      qsq=dsqrt(us2-q2)
      qlo=dlog((us+qsq)/q)
      ns=(us2-q2)/pi2*(qsq-2.d0*alfac/pi*(us-3.d0*q*q/qsq*
     xdlog((us+qsq)/ren)))
      uust=(0.5d0*(nd+ns)*pi2/qalf)**0.3333333333333333333d0
      nd0=nd
      ns0=ns
      ud0=ud
      uu=xnewt (qeq,uust,epsu)
c...
      ue=ud-uu
      nu=uu**3/pi2*qalf
      omu=-nu*uu/4.d0
      omd=-nd*ud/4.d0
      oms=-(us*qsq*(us2-2.5d0*q2)+1.5d0*q2*q2*qlo
     x-2.d0*alfac/pi*(3.d0*(us*qsq-q2*dlog((us+qsq)/ren))**2-2.d0*
     x(us2-q2)**2-3.d0*q2*q2*(dlog(q/ren))**2))/4.d0/pi2
      ne=ue**3/3.d0/pi2
      ome=-ue*ne/4.d0
      omega=omu+omd+oms+ome
      p=-omega-1.d0
      ro=omega+uu*nu+ud*nd+us*ns+ue*ne+1.d0
      n=(nu+nd+ns)/3.d0
      u=(ro+p)/n
      return
      end
c...
      subroutine qeq (uu,f,df)
      implicit double precision (a-h,n-z)
      common /doqeq/ qalf,nd,ns,ud
      common /con/ pi
      pi2=pi*pi
      f=(2.d0*uu**3*qalf+(uu-ud)**3)/pi2-nd-ns
      df=3.d0*(2.d0*uu**2*qalf+(uu-ud)**2)/pi2
      return
      end
c...
c    Procedura funkcyjna xnewt - metoda Newtona startujaca z xstart
c...
      double precision function xnewt (func,xstart,eps)
      implicit double precision (a-h,o-z)
      common /outnew/ istep,ier
      imax=200
      ier=0
      xacc=xstart*eps
      x=xstart
      do 200 i=1,imax
        call func (x,f,df)
	dx=-f/df
	x=x+dx
	if (dabs(dx).lt.xacc) go to 300
 200    continue
      ier=-2
      istep=imax
      write (*,*) ' xnewt - za duzo krokow'
      xnewt=x
      return
 300  x1=x+xacc
      call func (x1,f1,df1)
      x2=x-xacc
      call func (x2,f2,df2)
      if (f1*f2.gt.0.d0) then
        ier=-1
	write (*,*) ' xnewt - zero poza przedzialem <x-xacc,x+xacc>'
	endif
      xnewt=x
      istep=i
      return
      end
