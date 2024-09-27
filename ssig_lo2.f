c.
c. 2 -> 2 lo cross section
c.
      subroutine ssig_lo2(p,tot)
      implicit none
C  
C ARGUMENTS 
C  
      double precision p(0:3,5),tot
C
C LOCAL VARIABLES 
C
      DOUBLE PRECISION XPQ1(-6:6),XPQ2(-6:6)
      double precision f1(-6:6),f2(-6:6)
      double precision xl(4),sig(4)
      double precision s,t,u

      integer i
C
C EXTERNAL FUNCTIONS
C
      DOUBLE PRECISION dot,psi0
C
C GLOBAL VARIABLES
C
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
      elseif (icollide .eq. 3) then !PP
         do i=-6,6
            f1(i)=xpq1(i)
            f2(i)=xpq2(i)
         enddo
      endif
c
      call setlum(f1,f2,xl)
c
      sig(1)=psi0(s,t,u,xm2,wm2)       ! u d~ -> W+
      sig(2)=psi0(s,u,t,xm2,wm2)       ! d~ u -> W+
      sig(3)=sig(1)                    ! u~ d -> W-
      sig(4)=sig(2)                    ! d u~ -> W-
c
      tot=0d0
      if(itop.eq.1) tot=tot+xl(1)*sig(1)+xl(2)*sig(2)
      if(iatop.eq.1)tot=tot+xl(3)*sig(3)+xl(4)*sig(4)
c
      return
      end


      subroutine setlum(f1,f2,f)
c Returns luminosity
      implicit none
c
      double precision f1(-6:6),f2(-6:6),f(4)
c     g=0,d=1,u=2,s=3,c=4,b=5,anti = -n

      double precision vud2,vus2,vub2,vcd2,vcs2,vcb2,vtd2,vts2,vtb2
      common/ckm/vud2,vus2,vub2,vcd2,vcs2,vcb2,vtd2,vts2,vtb2
c
c Begin code
c
c  (u,c) (d~,s~,b~) -> W+
      f(1)=vtb2*(f1(2)*f2(-1)*vud2+f1(2)*f2(-3)*vus2+f1(2)*f2(-5)*vub2
     &     +f1(4)*f2(-1)*vcd2+f1(4)*f2(-3)*vcs2+f1(4)*f2(-5)*vcb2)
c  d~ u -> W+
      f(2)=vtb2*(f2(2)*f1(-1)*vud2+f2(2)*f1(-3)*vus2+f2(2)*f1(-5)*vub2
     &     +f2(4)*f1(-1)*vcd2+f2(4)*f1(-3)*vcs2+f2(4)*f1(-5)*vcb2)
c  u~ d -> W-
      f(3)=vtb2*(f1(-2)*f2(1)*vud2+f1(-2)*f2(3)*vus2+f1(-2)*f2(5)*vub2
     &     +f1(-4)*f2(1)*vcd2+f1(-4)*f2(3)*vcs2+f1(-4)*f2(5)*vcb2)
c  d u~ -> W-
      f(4)=vtb2*(f2(-2)*f1(1)*vud2+f2(-2)*f1(3)*vus2+f2(-2)*f1(5)*vub2
     &     +f2(-4)*f1(1)*vcd2+f2(-4)*f1(3)*vcs2+f2(-4)*f1(5)*vcb2)
c
      return
      end


      double precision function psi0(s,t,u,xm2,wm2)
c Returns the LO matrix element squared, averaged over spin/color
      implicit none
      double precision s,t,u,xm2,wm2
      psi0=0.25d0*t*(t-xm2)/(s-wm2)/(s-wm2)
      return
      end
