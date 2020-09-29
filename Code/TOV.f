      subroutine TOV(title,targetmass,rad_t,bar_t,rho_t,pres_t,
     1            emas_t,phi_t,imax)
     
      implicit real*8 (a-h,o-z)
      real*8 k1,k2,k3,k4,l1,l2,l3,l4,m1,m2,m3,m4,n1,n2,n3,n4
      character*150 title
      parameter (ipmax=10000)
     
      dimension e(0:50000), den(0:50000), r(0:50000), 
     x          p(0:50000), em(0:50000), baryden(0:50000),
     x          barynum(0:50000)
      dimension prho(ipmax), phi(0:50000)
      dimension rad_t(10000),bar_t(10000),rho_t(10000),pres_t(10000),
     1            emas_t(10000),phi_t(10000)

      common/interp/ i1,i2,i3,limit
      common/const/ pi,g
      common/ans/ bfunc, pmfunc, emfunc, pfunc, phfunc
      common/eos_dat/ peos(50000),eeos(50000),deos(50000),const
      common/prof/r,den,e,p,em,phi,rhoc,baryden,propmas,
     x delta,stepi,barynum
     
      logical output, stndrd, qseek,init,nmax
     
      pi=acos(-1.)
      g=1.484
     
C ***** Defines default values *******************************************
c      write(*,*)'target mass', targetmass
      stepi=.0005
      delta=.01
c      title='../EOS/APR_EOS_Cat_new.dat'
      limit=96
      output=.true.!.false.

      rhol=9.629e-2
c This correspond to a density of 1.6e14 g/cm3 for the core limit
      rhod=2.393e-4
c This correspond to a density of 4e11 g/cm3 for neutron drip
      init=.true.     

1     const = 1.d0
      c2=1.d0
c      targetmass=1.4
      nmax=.false.

      if (init) then
      rhoc=.4
      stndrd=.false.
      else 
c      stndrd=.true.
      stndrd=.false.
      endif
      if (stndrd) then
       rhoc=0.25
      endif
C ***** Initializes result arrays ****************************************
     
      do i=0,50000
      em(i)=0.
      den(i)=0.
      p(i)=0.
      r(i)=0.
      e(i)=0.
      phi(i)=0.
      end do
     
C ***** Waits for instructions *********************************************
     
c       output=.true.
ccc CGS units
        const = 1.0d0 / 1.989d18
        c2=9.d20

       call eos(title,c2)

2      iii=1
      rhoc=0.25
      omassold=0.42
      write(*,*)omassold
12    continue
      i1=2
      i2=2
      i3=2
     
C ***** Do the first step ********************************************
     
      den(0)=rhoc
      p(0)=pres(den(0))
      e(0)=ener(p(0))
      phi(0)=0.
      icore=0
      idrip=0
     
      k1=0.
      l1=0.
      m1=0.
      n1=0.
      p1=0.
      step=stepi
     
      call twostep(r(0)+step/2.,p(0)+k1/2.,em(0)+l1/2.,phi(0)+p1/2.)
      k2 = step *  pfunc
      l2 = step * emfunc
      m2 = step *  bfunc
      n2 = step * pmfunc
      p2 = step * phfunc
     
      call twostep(r(0)+step/2.,p(0)+k2/2.,em(0)+l2/2.,phi(0)+p2/2.)
      k3 = step *  pfunc
      l3 = step * emfunc
      m3 = step *  bfunc
      n3 = step * pmfunc
      p3 = step * phfunc
     
      call twostep(r(0)+step,p(0)+k3,em(0)+l3,phi(0)+p3)
      k4 = step *  pfunc
      l4 = step * emfunc
      m4 = step *  bfunc
      n4 = step * pmfunc
      p4 = step * phfunc
     
      p(1)  =  p(0) + (k1+2.*k2+2.*k3+k4)/6.
      em(1) = em(0) + (l1+2.*l2+2.*l3+l4)/6.
      r(1)  =  r(0) + step
      den(1)= rho(p(1))
      baryden(1) = (m1+2.*m2+2.*m3+m4)/6.
      barynum(1) =baryden(1)*2d33/1.67d-24
      propmas = (n1+2.*n2+2.*n3+n4)/6.
      e(1)  = ener(p(1))
      phi(1)   = phi(0)   + (p1+2.*p2+2.*p3+p4)/6.
     
C ***** Do the next step **********************************************
     
      do i=1,50000-1
     
      imax=i+1
     
      call twostep(r(i),p(i),em(i),phi(i))
     
      step = delta/(emfunc/em(i)-pfunc/p(i))
     
      k1 = step *  pfunc
      l1 = step * emfunc
      m1 = step *  bfunc
      n1 = step * pmfunc
      p1 = step * phfunc
     
      call twostep(r(i)+step/2.,p(i)+k1/2.,em(i)+l1/2.,phi(i)+p1/2.)
      k2 = step *  pfunc
      l2 = step * emfunc
      m2 = step *  bfunc
      n2 = step * pmfunc
      p2 = step * phfunc
     
      call twostep(r(i)+step/2.,p(i)+k2/2.,em(i)+l2/2.,phi(i)+p2/2.)
      k3 = step *  pfunc
      l3 = step * emfunc
      m3 = step *  bfunc
      n3 = step * pmfunc
      p3 = step * phfunc
     
      call twostep(r(i)+step,p(i)+k3,em(i)+l3,phi(i)+p3)
      k4 = step *  pfunc
      l4 = step * emfunc
      m4 = step *  bfunc
      n4 = step * pmfunc
      p4 = step * phfunc
     
      p(i+1)  = p(i) + (k1+2.*k2+2.*k3+k4)/6.
      em(i+1) = em(i) + (l1+2.*l2+2.*l3+l4)/6.
      r(i+1)  = r(i) + step
      den(i+1)= rho(p(i+1))
      baryden(i+1) = baryden(i) + (m1+2.*m2+2.*m3+m4)/6.
      barynum(i+1) =baryden(i+1)*2d33/1.67d-24
      propmas = propmas + (n1+2.*n2+2.*n3+n4)/6.
      e(i+1)  = ener(p(i+1))
      phi(i+1)   = phi(i)   + (p1+2.*p2+2.*p3+p4)/6.
     
      if ((icore .eq. 0).and.(den(i+1) .le. rhol)) icore=i
     
      if ((idrip .eq. 0).and.(den(i+1) .lt. rhod)) idrip=i
     
      if (den(i+1).lt.deos(limit-1)) goto 1100
     
      end do
     
      write(6,*) ' '
      write(6,*) ' The star surface has not been reached! '
      write(6,*) ' '
     
C ***** The calculation is finished ****************************************
     
1100  continue
     
C ***** Rescales phi ********************************************************
     
      const1=.5*dlog(1.-2.*g*em(imax)/r(imax))-phi(imax)
      do i=0,imax
        phi(i)=phi(i)+const1
      end do
     
C **********************************************************************
          
      iii=iii+1
      if (stndrd) then
ccccccccccccccccccccccccccccccccccccccccccccccc
        emax = em(imax)-targetmass
c        write(*,*)emax,rhoc,'here'

ccccccccccccccccccccccccccccccccccccccccccccccc
        accuracy=1.e-7
        if(.not.qseek(emax,rhoc,emold1,rhoold1,emold2,rhoold2,
     1                accuracy,iii-1))  then
        go to 12
        else
       if(output) call profile(title,imax,icore,idrip)
      write(6,*) ' '
      write(6,*) title
      write(6,*) ' '
      write(6,'('' rhoc:'', 1pe12.4, '' radius:'', 1pe12.4,/,
     x '' masses(g,b,p):'', 1p3e17.9)' )
     x    rhoc,r(imax),em(imax),baryden(imax),propmas
      write(*,*)'imax',imax
      do i=1,imax+1
      rad_t(i) =r(i-1)*1.e3*100
      bar_t(i) =den(i-1)
      rho_t(i) =e(i-1)*1.989e18
      pres_t(i)=p(i-1)*1.989e18*9e20
      emas_t(i)=em(i-1)
      phi_t(i) =phi(i-1)
      enddo
        endif
      else
       if ((em(imax)>omassold)) then
        if (nmax) go to 13
        write(*,*)omassold,em(imax),rhoc,stndrd,deos(1)
        omassold=em(imax)       
        rhoc=rhoc*1.005
        if (rhoc.ge.deos(1)) then
         rhoc=rhoc/1.005
         targetmass=omassold
         nmax=.true.
        endif
          go to 12      
        else 
13      targetmass=em(imax)
        write(*,*)'targetmass',targetmass,rhoc,deos(1)
        stndrd=.true.
        init=.false.
        go to 2
       endif
      endif
      
      if (init) then
      init=.false.
      go to 1
      endif
c      stop
      return
      end
     
c *************************************************************************
c *************************************************************************
     
      function ener(p)
     
      implicit real*8 (a-h,o-z)
     
      common/interp/ i1,i2,i3,limit
      common/eos_dat/ peos(50000),eeos(50000),deos(50000),const
     
1     ener=0.
     
      if ( p.ge.peos(i1) ) then
        ener = eeos(i1) * exp( log(p/peos(i1)) *
     x                         log(eeos(i1-1)/eeos(i1)) /
     x                         log(peos(i1-1)/peos(i1)) )
        return
      else
        i1=i1+1
      endif
     
      if(i1.gt.limit) then
        i1=limit
        return
      endif
     
      goto 1
     
      end
     
c ************************************************************************
c ************************************************************************
     
      function pres(d)
     
      implicit real*8 (a-h,o-z)
     
      common/interp/ i1,i2,i3,limit
      common/eos_dat/ peos(50000),eeos(50000),deos(50000),const
     
1     pres=0.
     
      if ( d.ge.deos(i2) ) then
        pres = peos(i2) * exp( log(d/deos(i2)) *
     x                         log(peos(i2-1)/peos(i2)) /
     x                         log(deos(i2-1)/deos(i2)) )
        return
      else
        i2=i2+1
      endif
     
      if(i2.gt.limit) then
        i2=limit
        return
      endif
     
      goto 1
     
      end
     
c **************************************************************************
c **************************************************************************
     
      function rho(p)
     
      implicit real*8 (a-h,o-z)
     
      common/interp/ i1,i2,i3,limit
      common/eos_dat/ peos(50000),eeos(50000),deos(50000),const
     
1     rho=0.
     
      if ( p.ge.peos(i3) ) then
        rho = deos(i3) * exp( log(p/peos(i3)) *
     x                         log(deos(i3-1)/deos(i3)) /
     x                         log(peos(i3-1)/peos(i3)) )
        return
      else
        i3=i3+1
      endif
     
      if(i3.gt.limit) then
        i3=limit
        return
      endif
     
      goto 1
     
      end
     
c *************************************************************************
c *************************************************************************
c *************************************************************************
     
      subroutine eos(title,c2)
     
      implicit real*8 (a-h,o-z)
      character*150 title
     
      common/eos_dat/ peos(50000),eeos(50000),deos(50000),const
      common/interp/ i1,i2,i3,limit
      write(*,*)title
      open (unit=15, file=title, status='old')
       itext=6
       limit=1000
       do i=1,itext
        read(15,*)
       end do
       do i=1,limit
        read(15,*,err=100,end=100) x1,x2,x3
c Check order of columns in the EOS file:
        if(i.eq.1) then
         if ((x3.le.10.).and.(x2.ge.1e30)) then
          ilist=1
         else if ((x1.le.10.).and.(x3.ge.1.e30))then
          ilist=2
         else
c          pause 'Check EOS column ordering !'
         end if
        end if
        if (ilist.eq.1) then
         eeos(i)=x1
         peos(i)=x2
         deos(i)=x3
        else if (ilist.eq.2) then
         eeos(i)=x2
         peos(i)=x3
         deos(i)=x1
        else
c         pause 'Check EOS column ordering !'
        end if
        peos(i)=peos(i) * const / c2
        eeos(i)=eeos(i) * const
       end do
 100  continue
      limit=i-1
      close(unit=15, status='keep')

c      do i=1,limit
c       write(6,10010) i,eeos(i)/const,peos(i)/const*c2,deos(i)
c      end do
10010 format(1i5,1p3e16.4)
     
      return
      end
     
c ************************************************************************
c ************************************************************************
c ************************************************************************
     
      subroutine twostep (r, p, em, phi )
     
      implicit real*8 (a-h,o-z)
      real*8 infunc, nodrag,  nonrel
     
      common/ans/ bfunc, pmfunc, emfunc, pfunc, phfunc
      common/const/ pi,g
     
      dens = rho(p)
      energy = ener(p)
      vol = 4.*pi*r*r*r
     
      bfunc = dens *4.*pi*r*r / sqrt( 1.-2.*g*em / r )
      bfunc = bfunc * 8.42e-4
     
c  this converts to solar masses from [#-km**3/fm**3]
c  I hope. this is the baryonic mass and is good for
c  computing the gravitational binding energy of a star.
     
      pmfunc = energy *4.*pi*r*r / sqrt( 1.-2.*g*em / r )
     
      pfunc=-g * (energy+p) * (em+vol*p) / ( r*(r-2.*em*g) )
     
      emfunc= 4.*pi*r*r * energy
     
      phfunc = g * (em+vol*p) / ( r*(r-2.*em*g) )
     
      return
     
      end
     
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
     
      subroutine profile(title,imax,icore,idrip)
     
      implicit real*8 (a-h,o-z)
      character*70 title
      character*2 error
     
      dimension e(0:50000), den(0:50000), r(0:50000), 
     x          p(0:50000), em(0:50000), baryden(0:50000),
     x          barynum(0:50000)
      dimension phi(0:50000)
     
      common/prof/r,den,e,p,em,phi,rhoc,baryden,propmas,
     x delta,stepi,barynum
     
      write(6,*) ' Profile being output. '
     
      open (unit=15,file='prof.dat')
      write(15,10001) 6, imax, icore, idrip
      write(15,*)
      write(15,*) '    EOS file :  ',title
      write(15,*)
      write(15,10002)'step','radius ','baryon#  ','density   ',
     x               '  pressure  ','encl. mass   ','phi   ',
     x               'encl. bar. mass ', 'encl. bar. num '
      write(15,10002) '   ','  (m)  ','(#/fm3)  ','(g/cm3)   ',
     x               '  (dyn/cm2) ',' (sol. mass)  ','     ',
     x               ' (sol. mass)   '
      write(15,*)
     
      do i=0,imax
       write(15,10000) i,r(i)*1.e3,den(i),e(i)*1.989e18,
     x                 p(i)*1.989e18*9e20,em(i),phi(i),
     x                 baryden(i),barynum(i)
      end do
     
      close (unit=15, status='keep')
     
10000 format(i6,0pf15.6,1p1e15.6,1e15.6,1e15.5,1e18.9,1e15.6,1e18.9
     x ,1e15.6,1e18.9)
10001 format(4i8)
10002 format(a6,a15,a15,a15,a15,a18,a15,a18,a18,a18)

      return
     
      end
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
c          Below is included content of old file "seek.for"
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************    
c **
      logical function seek(f1,x1,f2,x2,tol,i)
c **
      implicit real*8 (a-h,o-z)
c **
      seek=.true.
      if(abs(f1).lt.tol) return
c **
      seek=.false.
      if(i.eq.1) then
        xnew=.9*x1
      else
        xnew=(f2*x1-f1*x2)/(f2-f1)
      endif
c **
      f2=f1
      x2=x1
      x1=xnew
c **
      return
      end
c **
c **
      logical function qseek(f1,x1,f2,x2,f3,x3,tol,i)
c **
      implicit real*8 (a-h,o-z)
      logical seek
c **
      qseek=.true.
      if(abs(f1).lt.tol) return
c **
      qseek=.false.
      if (i.gt.2) then
c **
        a=(x1-x2)*f3+(x2-x3)*f1+(x3-x1)*f2
        a=-a/((x1-x2)*(x2-x3)*(x3-x1))
        b=(f1-f3)/(x1-x3)-(x1+x3)*a
        c=f2-a*x2*x2-b*x2
        d=b*b-4.0*a*c
        xnew=-0.5*b/a
c **
        if (d.lt.0.0) then
          if (abs(xnew-x1).lt.tol.or.abs(f1-f2).lt.tol) qseek=.true.
        else if (xnew.lt.x1) then
          xnew=xnew + 0.5*abs(sqrt(d)/a)
        else
          xnew=xnew - 0.5*abs(sqrt(d)/a)
        endif
c **
        f3=f2
        f2=f1
        x3=x2
        x2=x1
        x1=xnew
c **
      else
c **
        f3=f2
        x3=x2
        qseek=seek(f1,x1,f2,x2,tol,i)
        x1=max( x2*x2/x1, min( x1, x1*x1/x2 ) )
c **
      endif
c **
      return
      end
c **

