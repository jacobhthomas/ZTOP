c.
c. 2 -> 2 nlo cross section
c.
      subroutine ssig_nlo2(p,tot)
      implicit none
C  
C ARGUMENTS 
C  
      double precision p(0:3,5),tot
C  
C CONSTANTS
C  
      double precision pi
      parameter       (pi=3.14159265358979d0)
      double precision zeta2
      parameter       (zeta2=1.64493406684823d0)
C
C LOCAL VARIABLES 
C
      DOUBLE PRECISION XPQ1(-6:6),XPQ2(-6:6)
      double precision f1(-6:6),f2(-6:6)
      double precision f1t(-6:6),f2t(-6:6)
      double precision xl(4),xl1(4),xl2(4),sig(4)
      double precision col,sft,vrt,subtot,vrtnb(4)
      double precision s,t,u,beta,lam

      integer i
C
C EXTERNAL FUNCTIONS
C
      DOUBLE PRECISION dot,psi0,xlog,xlog2,dilog
C
C GLOBAL VARIABLES
C
      double precision deltas,deltac
      common/slice/deltas,deltac
      double precision hbarc2,gweak2,xm2,wm2
      common/parm/hbarc2,gweak2,xm2,wm2
      double precision xa,xb,scale
      common/toevent/  xa,xb,scale
      integer              icollide,ismear,igroup,iset,iq2
      common /to_collider2/icollide,ismear,igroup,iset,iq2
      integer      itop,iatop
      common/ttbar/itop,iatop
c
c Begin code
c
      s=2d0*dot(p(0,1),p(0,2))
      t=-2d0*dot(p(0,1),p(0,3))
      u=xm2-t-s
      beta=1d0-xm2/s
      lam=1d0/beta

      CALL PFTOPDG(xa,SCALE,XPQ1)
      CALL PFTOPDG(xb,SCALE,XPQ2)
      DO I=-6,6
         XPQ1(i)=XPQ1(i)/xa
         XPQ2(i)=XPQ2(i)/xb
      ENDDO
      if (icollide .le. 2) then ! f1 is proton, f2 is antiproton
         do i=-6,6
            f1(i)=xpq1(i)
            f2(i)=xpq2(-i)
         enddo
         call lumt3(xa,scale,1,5,deltas,deltac,s,f1,f1t)
         call lumt3(xb,scale,-1,5,deltas,deltac,s,f2,f2t)
      elseif (icollide .eq. 3) then !PP
         do i=-6,6
            f1(i)=xpq1(i)
            f2(i)=xpq2(i)
         enddo
         call lumt3(xa,scale,1,5,deltas,deltac,s,f1,f1t)
         call lumt3(xb,scale,1,5,deltas,deltac,s,f2,f2t)
      endif
c
      call setlum(f1,f2,xl)
      call setlum(f1t,f2,xl1)
      call setlum(f1,f2t,xl2)
c
      col = (4d0*xlog(deltas)+3d0)*xlog(s/scale**2)             ! 2 * initial
     &     + 3.5d0 - 2d0*zeta2 - xlog2(deltas) - xlog2(beta)  ! final...
     &     + 2d0*xlog(deltas)*xlog(beta) - xlog(deltac)*
     &     (2d0*xlog(deltas)+1.5d0-2d0*xlog(beta))
      col = col*4d0/3d0      ! C_F = 4/3

      vrt = xlog2(s/xm2)/4d0 + 5d0/4d0*xlog(s/xm2) + xlog(1d0-lam)*
     &     xlog(s/xm2) + 7d0 + (1d0+beta/2d0)*xlog(1d0-lam)
     &     +xlog2(1d0-lam)/2d0 - dilog(lam) - zeta2
      vrt = -vrt*8d0/3d0     ! -2*C_F = -8/3

c non-Born piece
c u d~ -> W+
      vrtnb(1) = t*u*(1d0-beta)*xlog(1d0-lam)/(s-wm2)/(s-wm2)/3d0
c d~ u -> W+,   t <-> u
      vrtnb(2) = vrtnb(1)
c u~ d -> W-
      vrtnb(3) = vrtnb(1)
c d u~ -> W-
      vrtnb(4) = vrtnb(2)

c      call setsoft(s,sft)
      sft = 6d0*xlog2(deltas) - 2d0*xlog(deltas) +2d0*xlog(deltas)*
     &     log(s/xm2) + (s+xm2)/(s-xm2)*xlog(s/xm2) - 2d0*dilog(beta)
     &     - xlog2(s/xm2)/2d0
      sft = sft*4d0/3d0      ! C_F = 4/3

      subtot = vrt+sft+col

      sig(1)=psi0(s,t,u,xm2,wm2)       ! u d~ -> W+
      sig(2)=psi0(s,u,t,xm2,wm2)       ! d~ u -> W+
      sig(3)=sig(1)                    ! u~ d -> W-
      sig(4)=sig(2)                    ! d u~ -> W-
c
      tot=0d0
      if(itop.eq.1)tot=tot+(xl(1)*subtot+xl1(1)+xl2(1))*sig(1)+xl(1)
     &     *vrtnb(1)+(xl(2)*subtot+xl1(2)+xl2(2))*sig(2)+xl(2)*vrtnb(2)
      if(iatop.eq.1)tot=tot+(xl(3)*subtot+xl1(3)+xl2(3))*sig(3)+xl(3)
     &     *vrtnb(3)+(xl(4)*subtot+xl1(4)+xl2(4))*sig(4)+xl(4)*vrtnb(4)
c
      return
      end

c. log
      function xlog(z)
      implicit double precision (a-h,o-z)
      xlog=dlog(dabs(z))
      return
      end

c. log squared with negative argument
      function xlog2(z)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)
      xlog2=dlog(dabs(z))**2
      if(z.lt.0d0)xlog2=xlog2-pi**2
      return
      end
