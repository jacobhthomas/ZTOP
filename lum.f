      subroutine lumt3(x,xmu,nt,nf,deltas,deltac,s12,fin,fout)
      implicit double precision (a-h,o-z)
      dimension fin(-6:6),fout(-6:6),ftmp(-6:6)
      dimension fvec(1)
      parameter (n=3)
      parameter (cf=4d0/3d0)

      omx=1d0-x
      xjac=omx
      call brm48(fvec,1)

      y=x+xjac*fvec(1)
      atmp=deltac*s12/(xmu*xmu)
      xlt=dlog(atmp)
      omy=1d0-y
      xly=xlt+dlog(omy/y)

      ptqg=xly*pqg0(y)-pqg1(y)
      ptgq=xly*pgq0(y)-pgq1(y)

      CALL PFTOPDG(x/y,xmu,ftmp)

      do i=-6,6
         ftmp(i)=ftmp(i)*y/x
      enddo
      if(nt.eq.-1)then
         do i=1,6
            ftmp1=ftmp(-i)
            ftmp(-i)=ftmp(i)
            ftmp(i)=ftmp1
         enddo
      endif

      fsum=0d0
      do i=1,nf
         fsum=fsum+ftmp(-i)+ftmp(+i)
      enddo

      oms=1d0-deltas
      if((x.le.oms).and.(y.le.oms))then
         f34=gfun(oms,atmp)-gfun(x,atmp)
         f34=f34/(omx-deltas)

         ptqq=xly*pqq0(y)-pqq1(y)
         ptgg=xly*pgg0(y)-pgg1(y)

         do i=1,nf
            j=+i
            fout(j) = 
     .                2d0*cf*fin(j)*f34
     .              + (ftmp(j)*ptqq-fin(j)*2d0*cf*xly/omy)/y
     .              + ftmp(0)*ptqg/y
            fout(j) = fout(j)*xjac
            j=-i
            fout(j) = 
     .                2d0*cf*fin(j)*f34
     .              + (ftmp(j)*ptqq-fin(j)*2d0*cf*xly/omy)/y
     .              + ftmp(0)*ptqg/y
            fout(j) = fout(j)*xjac
         enddo
         
         fout(0) = 
     .             2d0*dble(n)*fin(0)*f34
     .           + (ftmp(0)*ptgg-fin(0)*2d0*dble(n)*xly/omy)/y
     .           + fsum*ptgq/y
         fout(0) = fout(0)*xjac
         
      else
         do i=1,nf
            fout(+i) = xjac*ftmp(0)*ptqg/y
            fout(-i) = fout(+i)
         enddo
         
         fout(0) = xjac*fsum*ptgq/y
         
      endif

      return
      end

      function gfun(y,a)
      implicit double precision (a-h,o-z)
      tmp=dlog(a*(1d0-y)/y)
      gfun=-0.5d0*tmp*tmp
      return
      end
