      SUBROUTINE newt(x,n,check,func,conv)
      implicit none
      INTEGER n,NP,MAXITS
      PARAMETER(NP=10,MAXITS=500000)

      LOGICAL check,conv,sing
      REAL*8 x(n),fvec(NP),TOLF,TOLX,STPMX,EPS
      PARAMETER(TOLF=5.d-7,TOLX=1.d-9,STPMX=0.1,EPS=5.d-4)
CU    USES fmin1,lnsrch,lubksb,ludcmp
      INTEGER i,its,j,k,l
      REAL*8 f,fold,stpmax,sum,temp,test,h
      REAL*8 fjac(NP,NP),g(NP),p(NP),xold(NP),ftemp(NP),fmin
      EXTERNAL func

      check = .false.
      conv=.true.
c----------
c     Calculation of first f value
c----------
      CALL func(x,fvec,n)
      f=fmin(n,fvec)

      test=0.
      do 11 i=1,n
        if(dabs(fvec(i)).gt.test) test=dabs(fvec(i))
11    continue
      if(test.lt.1.d-2*TOLF)return

      sum=0.
      do 12 i=1,n
        sum=sum+x(i)*x(i)
12    continue
      stpmax=STPMX*max(sqrt(sum),dble(n))

      do 21 its=1,MAXITS

c----------
c     Calculation of Jacobian matrix
c----------

        do k=1,n
          temp=x(k)
          h=EPS*abs(temp)
          if(h.eq.0.)h=1.d-5
          x(k)=temp+h
          h=x(k)-temp
          CALL FUNC(x,ftemp,n)
          x(k)=temp
          do l=1,n
            fjac(l,k)=(ftemp(l)-fvec(l))/h
          enddo
        enddo

        do 14 i=1,n
          sum=0.
          do 13 j=1,n
            sum=sum+fjac(j,i)*fvec(j)
13        continue
          g(i)=sum
14      continue
        do 15 i=1,n
          xold(i)=x(i)
15      continue
        fold=f
        do 16 i=1,n
          p(i)=-fvec(i)
16      continue

        call lusolve(fjac,n,NP,p,sing)
        if (sing.eqv..false.) then
        call lnsrch(n,xold,fold,g,p,x,f,stpmax,check,fvec,func)

        test=0.
        do 17 i=1,n
          if(dabs(fvec(i)).gt.test) test=dabs(fvec(i))
17      continue
        if(test.lt.TOLF)then
          check=.false.
          return
        endif

        test=0.
        do 19 i=1,n
          if (x(i).eq.0.d0) then
            temp=dabs(xold(i))
          else  
            temp=dabs(1.d0-xold(i)/x(i))
          endif  
          if(temp.gt.test)test=temp
19      continue
        if(test.lt.TOLX) then
          check=.true.
          return
        endif
        else 
        go to 22
        endif
21    continue
c      print*,'MAXITS exceeded in newt',f
22      conv=.false.
c      do i=1,n
c        print*,i,x(i),fvec(i)
c      enddo
c      pause
      END



      FUNCTION fmin(n,fvec)
      implicit none
      INTEGER i,n
      REAL*8 fmin,fvec(n),sum
      
      sum=0.
      do i=1,n
        sum=sum+fvec(i)*fvec(i)
      enddo
      fmin=0.5*sum

      return
      END




      SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check,fvec,func)
      implicit none
      INTEGER n
      LOGICAL check
      REAL*8 f,fold,stpmax
      REAL*8 g(n),p(n),x(n),xold(n),fvec(n),ALF,TOLX,fmin
      PARAMETER (ALF=1.e-4,TOLX=1.e-6)
      EXTERNAL func
CU    USES func
      INTEGER i
      REAL*8 a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2
      REAL*8 slope,sum,temp,test,tmplam

      check=.false.
      sum=0.
      do 11 i=1,n
        sum=sum+p(i)*p(i)
11    continue
      sum=sqrt(sum)

      if(sum.gt.stpmax)then
        do 12 i=1,n
          p(i)=p(i)*stpmax/sum
12      continue
      endif

      slope=0.
      do 13 i=1,n
        slope=slope+g(i)*p(i)
13    continue

      test=0.
      do 14 i=1,n
        temp=abs(p(i))/max(abs(xold(i)),1.d0)
        if(temp.gt.test)test=temp
14    continue
      alamin=TOLX/test

      alam=1.
1     continue

        do 15 i=1,n
          x(i)=xold(i)+alam*p(i)
15      continue

        CALL func(x,fvec,n)
        f=fmin(n,fvec)

        if(alam.lt.alamin)then

          do 16 i=1,n
            x(i)=xold(i)+alamin*p(i)
16        continue

cx          check=.true.
          return

        else if(f.le.fold+ALF*alam*slope)then

          return

        else

          if(alam.eq.1.)then
            tmplam=-slope/(2.*(f-fold-slope))
          else
            rhs1=f-fold-alam*slope
            rhs2=f2-fold2-alam2*slope
            a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
            if(a.eq.0.)then
              tmplam=-slope/(2.*b)
            else
              disc=b*b-3.*a*slope
              if (disc.lt.0.d0) then
c                print*,'disc=',disc
c                print*,b*b,3.*a*slope
                disc=0.d0
              endif
              tmplam=(-b+sqrt(disc))/(3.*a)
            endif
            if(tmplam.gt..5*alam)tmplam=.5*alam
          endif

        endif

        alam2=alam
        f2=f
        fold2=fold
        alam=max(tmplam,.1*alam)

      goto 1

      END



      SUBROUTINE lusolve(a,n,np,b,sing)
      implicit none
      INTEGER n,np,NMAX
      REAL*8 a(np,np),TINY,b(n)
      PARAMETER (NMAX=20,TINY=1.0e-20)
      INTEGER i,imax,j,k,ii,ll,indx(NMAX)
      REAL*8 aamax,dum,sum,vv(NMAX)
      logical sing

      sing=.false.       
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.lt.TINY) then
c        print*,'singular matrix in ludcmp'
        sing=.true.
        endif
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue

      ii=0
      do 22 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 21 j=ii,i-1
            sum=sum-a(i,j)*b(j)
21        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
22    continue
      do 24 i=n,1,-1
        sum=b(i)
        do 23 j=i+1,n
          sum=sum-a(i,j)*b(j)
23      continue
        b(i)=sum/a(i,i)
24    continue

      return
      END
      
