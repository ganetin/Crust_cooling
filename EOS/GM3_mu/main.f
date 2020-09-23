      implicit none
      integer i,j,k,neq,l,icore

      include 'const.inc'
      real*8 hbarc,tmev,nb0
      parameter (hbarc=197.33)

      real*8 nb,ye,yp,yn,yl,ynu,t,s,p,eps,muh,mun,mue,ym,mns,mtot
      common /EOS_IO/nb,ye,yp,yn,yl,ynu,t,s,p,eps,muh,mun,mue,ym,mns
      real*8 ynui,rho0,ylprev
      real*8 nb1,e1,mue1,mun1,mup1,p1,s1,munu
      real*8 dpdn, dpdt, dedn, dedt, dsdn, dsdt, cs2, gam1
      real*8 y(7),y0(7),var
      logical check,conv
      real*8 ptot(600),utot(600),rho(600)
      real*8 n0     
      EXTERNAL T_NP

      data y0 /0.07, 0.04, -0.006, 4.7, 0.1, 0.5, 0.01/

C---------------------------------------------
C---------------------------------------------

        k=0
        open (unit=15, file='EOS.dat')
        neq=6
        do l=1,7
           y(l)=y0(l)
         enddo
cc Nuclear input
         tmev= 1d-2
         Ynui=0.001d0
         Ylprev=0.25
        rho0=1d16
        i=1

        DO WHILE (rho0.gt.1.6d14)
        conv=.true.
        nb =1.5-0.01*(i-1.)
        if (nb.le.0.d0) then
        i=i-1
        go to 13
        endif
        do j=1,80001
        yl=Ylprev-(j-1)*0.00001
        t=tmev/hbarc
        CALL NEWT(y,neq,check,T_NP,conv)

        nb1 = nb
        p1   = hbarc*p   
        p1=p1*1.6022d33 ! in dyn/cm^2  
        mue1 = hbarc*mue
        mun1 = hbarc*(mun-4.7535)
        mup1 = mun1 - hbarc*muh
        munu = hbarc*(mue-muh)
        e1   = hbarc*eps
        s1   = s
        tmev = hbarc*t
        rho0=eps*3.5178d14 ! in g.cm^-3
        if (ynu.le.ynui) then 
        if (conv) then
        k=k+1
c   Rho	   Press	    nbar	    Ye	    Ymu	    Yn	    Yp	    Yla	    Ysm	    Ys0	    Ysp	Yxim	   mstp	   mstn	   mstla	   mstsm	   msts0	   mstsp	mstxim

        write(15,111)rho0, p1,  nb1, Ye,ym,yn,yp,0.d0,0.d0,0.d0,0.d0,
     1  0.d0,mns ,mns, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0
111    format(19(e13.6,2x))
        write(*,10012)i,rho0,ye,ynu,yl!,conv
10012 format(1i5,1p1e16.4,1p1e16.4,1p1e16.4,1p1e16.4)!,2x,4L)
        go to 12
       endif
       endif
      enddo
12      continue
       i=i+1
       Ylprev=Yl
      ENDDO
13        icore=k
       write(*,*)'EOS calculated'
       write(*,*)'Output : EOS.dat' 
       write(*,*)'Icore =',icore
      end


