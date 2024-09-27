c.
c. 2 -> 2 lo cross section
c.
c Modified 9/17/02 to allow arbitrary scaling via isclf/isclfb
c Modified 6/18/24 to loop over multiple PDF sets
c     
      subroutine tsig_lo2(p,tot)
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
      DOUBLE PRECISION XPQB1(-6:6),XPQB2(-6:6)
      double precision fb1(-6:6),fb2(-6:6)
      double precision xl(8),sig(8)
      double precision s,t,u

      double precision scale1,scale2,scaleb1,scaleb2

      double precision totpdf(0:58)
      integer iloop
      
      integer i
C
C EXTERNAL FUNCTIONS
C
      DOUBLE PRECISION tpsi0
      DOUBLE PRECISION dot,pt,alphas2
C
C GLOBAL VARIABLES
C
      double precision hbarc2,gweak2,xm2,wm2
      common/parm/hbarc2,gweak2,xm2,wm2
      integer      isclf,ialphas,iscale
      common/flags/isclf,ialphas,iscale
      integer       isclfb,iscaleb
      common/flagsb/isclfb,iscaleb
      double precision xa,xb,scale
      common/toevent/  xa,xb,scale
      integer              icollide,ismear,igroup,iset,iq2
      common /to_collider2/icollide,ismear,igroup,iset,iq2
      integer      itop,iatop
      common/ttbar/itop,iatop
      double precision bwgt
      common/btagging/ bwgt

      double precision ration(100)
      common/pdfratios/ration

      integer MXX, MXQ, MXF, MaxVal, MXPQX, MXPQXTOT
      PARAMETER (MXX = 201, MXQ = 40, MXF = 6, MaxVal=4)
      PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)
      parameter (MXPQXTOT = MXPQX*59)
      integer icallset
      double precision updmulti(MXPQXTOT)
      common / ctmulti / updmulti, icallset

c
c Begin code
c
      tot=0d0
      bwgt=0d0
c
      s=2d0*dot(p(0,1),p(0,2))
      t=-2d0*dot(p(0,1),p(0,3))
      u=xm2-t-s
c
c Set scale, pdf's
c
      if(iscale.gt.7) then
c Do not choose scale 8 unless cutting out 2 -> 2 events.
      if(iscale.eq.8) then
         scale1=pt(p(0,3))               ! M_jj^2, but for 2 -> 2, M_jj=0
         scale2=scale1                   ! so use pt_d
      elseif (iscale.eq.9) then
         scale1=dsqrt(-t)                ! -t = Q^2 out of p, Q^2+M_jj^2
         scale2=dsqrt(-u)                ! -u = Q^2 out of pbar, M_jj=0
      elseif (iscale.eq.10) then
         scale1=pt(p(0,3))               ! MT_jj^2, but for 2 -> 2,
         scale2=scale1                   !  MT_jj=pt_d
      elseif (iscale.eq.11) then
         scale1=dsqrt(-t+(pt(p(0,3)))**2) ! -t = Q^2 out of p, Q^2+MT_jj^2
         scale2=dsqrt(-u+(pt(p(0,3)))**2) ! -u = Q^2 out of pbar, MT_jj=pt_d
      elseif (iscale.eq.12) then
         scale1=pt(p(0,3))+pt(p(0,4))                 ! HT = sum PT
         scale2=scale1
      elseif (iscale.eq.13) then
         scale1=pt(p(0,3))+dsqrt((pt(p(0,4)))**2+xm2) ! HT~ = sum MT
         scale2=scale1
      endif
      else            ! normal scale choices
      if(iscale.eq.4) then
         scale1=dsqrt(xm2)
         scale2=scale1
      elseif (iscale.eq.5) then
         scale1=dsqrt(s)
         scale2=scale1
      elseif (iscale.eq.6) then
         scale1=dsqrt(-t)                ! -t = Q^2 out of p
         scale2=dsqrt(-u)                ! -u = Q^2 out of pbar
      elseif (iscale.eq.7) then
         scale1=dsqrt(-t+xm2)            ! -t = Q^2 out of p
         scale2=dsqrt(-u+xm2)            ! -u = Q^2 out of pbar
      elseif (iscale.eq.0) then
         scale1=pt(p(0,4))
         scale2=scale1
      elseif (iscale.eq.1) then
         scale1=dsqrt((pt(p(0,4)))**2+xm2)
         scale2=scale1
      elseif (iscale.eq.2) then
         scale1=pt(p(0,3))               ! Highest ET jet (only 1 here)
         scale2=scale1
      elseif (iscale.eq.3) then
         scale1=dsqrt(xm2+2d0*dot(p(0,4),p(0,3))) ! 3 is highest ET jet
         scale2=scale1
      endif
      endif
      if(isclf.ne.0) then
         scale1=scale1*(dble(isclf)/1d2)
         scale2=scale2*(dble(isclf)/1d2)
      endif
      if((scale1.lt.1d0).or.(scale2.lt.1d0)) return

      if(iscaleb.gt.7) then
      if(iscaleb.eq.8) then
         scaleb1=dsqrt(xm2)              ! M_tj^2, but for 2 -> 2, M_tj=mt
         scaleb2=scaleb1
      elseif (iscaleb.eq.9) then
         scaleb1=dsqrt(-t+xm2)           ! -t = Q^2 out of p, Q^2+M_tj^2
         scaleb2=dsqrt(-u+xm2)           ! -u = Q^2 out of pbar, M_tj=mt
      elseif (iscaleb.eq.10) then
         scaleb1=dsqrt((pt(p(0,4)))**2+xm2)  ! MT_tj^2, but for 2 -> 2,
         scaleb2=scaleb1                     !  MT_tj=MT_t
      elseif (iscaleb.eq.11) then
         scaleb1=dsqrt(-t+(pt(p(0,4)))**2+xm2) ! -t = Q^2 out of p, Q^2+MT_tj^2
         scaleb2=dsqrt(-u+(pt(p(0,4)))**2+xm2) ! -u = Q^2 out of pbar
      elseif (iscaleb.eq.12) then
         scaleb1=pt(p(0,3))+pt(p(0,4))                 ! HT = sum PT
         scaleb2=scaleb1
      elseif (iscaleb.eq.13) then
         scaleb1=pt(p(0,3))+dsqrt((pt(p(0,4)))**2+xm2) ! HT~ = sum MT
         scaleb2=scaleb1
      endif
      else            ! normal scale choices
      if(iscaleb.eq.4) then
         scaleb1=dsqrt(xm2)
         scaleb2=scaleb1
      elseif (iscaleb.eq.5) then
         scaleb1=dsqrt(s)
         scaleb2=scaleb1
      elseif (iscaleb.eq.6) then
         scaleb1=dsqrt(-t)                ! -t = Q^2 out of p
         scaleb2=dsqrt(-u)                ! -u = Q^2 out of pbar
      elseif (iscaleb.eq.7) then
         scaleb1=dsqrt(-t+xm2)            ! -t = Q^2 out of p
         scaleb2=dsqrt(-u+xm2)            ! -u = Q^2 out of pbar
      elseif (iscaleb.eq.0) then
         scaleb1=pt(p(0,4))
         scaleb2=scaleb1
      elseif (iscaleb.eq.1) then
         scaleb1=dsqrt((pt(p(0,4)))**2+xm2)
         scaleb2=scaleb1
      elseif (iscaleb.eq.2) then
         scaleb1=pt(p(0,3))               ! Highest ET jet (only 1 here)
         scaleb2=scaleb1
      elseif (iscaleb.eq.3) then
         scaleb1=dsqrt(xm2+2d0*dot(p(0,4),p(0,3))) ! 3 is highest ET jet
         scaleb2=scaleb1
      endif
      endif
      if(isclfb.ne.0) then
         scaleb1=scaleb1*(dble(isclfb)/1d2)
         scaleb2=scaleb2*(dble(isclfb)/1d2)
      endif
      if((scaleb1.lt.1d0).or.(scaleb2.lt.1d0)) return
c

cz Move to here to enable PDF loop
c      as1=alphas2(scale1)       ! NOT used here
c      as2=alphas2(scaleb1)

      sig(1)=tpsi0(s,t,u,xm2,wm2)       ! u  b -> W+
      sig(2)=tpsi0(u,t,s,xm2,wm2)       ! d~ b -> W+
      sig(3)=tpsi0(s,u,t,xm2,wm2)       ! b u  -> W+
      sig(4)=tpsi0(t,u,s,xm2,wm2)       ! b d~ -> W+
      sig(5)=sig(1)                     ! u~ b~ -> W-
      sig(6)=sig(2)                     ! d  b~ -> W-
      sig(7)=sig(3)                     ! b~ u~ -> W-
      sig(8)=sig(4)                     ! b~ d  -> W-
c
      
      do iloop = 0, 58
         tot = 0d0
         icallset = iloop
c     
c Get PDFs
c
      CALL PFTOPDG(xa,SCALE1,XPQ1)
      CALL PFTOPDG(xb,SCALEB1,XPQB1)
      CALL PFTOPDG(xb,SCALE2,XPQ2)
      CALL PFTOPDG(xa,SCALEB2,XPQB2)
      DO I=-6,6
         XPQ1(i)=XPQ1(i)/xa
         XPQB1(i)=XPQB1(i)/xb
         XPQ2(i)=XPQ2(i)/xb
         XPQB2(i)=XPQB2(i)/xa
      ENDDO
      if (icollide .le. 2) then ! f1 is proton, f2 is antiproton
         do i=-6,6
            f1(i)=xpq1(i)
            f2(i)=xpq2(-i)
            fb1(i)=xpqb1(-i)
            fb2(i)=xpqb2(i)
         enddo
      elseif (icollide .eq. 3) then !PP
         do i=-6,6
            f1(i)=xpq1(i)
            f2(i)=xpq2(i)
            fb1(i)=xpqb1(i)
            fb2(i)=xpqb2(i)
         enddo
      endif
c
      call setlumt(f1,f2,fb1,fb2,xl)
c
      if(itop.eq.1) tot=tot+xl(1)*sig(1)+xl(2)*sig(2)+xl(3)*sig(3)
     &     +xl(4)*sig(4)
      if(iatop.eq.1)tot=tot+xl(5)*sig(5)+xl(6)*sig(6)+xl(7)*sig(7)
     &     +xl(8)*sig(8)
c
      totpdf(iloop) = tot
      enddo
cz Set lumi ratios here
      do i = 0, 58
         ration(i+1) = totpdf(i)/totpdf(0)
      enddo
      tot = totpdf(0)
cz//
      return
      end


      subroutine setlumt(f1,f2,fb1,fb2,f)
c Returns luminosity
      implicit none
c
      double precision f1(-6:6),f2(-6:6),f(8)
      double precision fb1(-6:6),fb2(-6:6)
c     g=0,d=1,u=2,s=3,c=4,b=5,anti = -n

      double precision vud2,vus2,vub2,vcd2,vcs2,vcb2,vtd2,vts2,vtb2
      common/ckm/vud2,vus2,vub2,vcd2,vcs2,vcb2,vtd2,vts2,vtb2
c
c Begin code
c
c TEST code
c      f(1)=(f1(2)+f1(4))*fb2(5)
c      f(2)=(f1(-1)+f1(-3)+f1(-5))*fb2(5)
c      f(3)=(f2(2)+f2(4))*fb1(5)
c      f(4)=(f2(-1)+f2(-3)+f2(-5))*fb1(5)
c      f(5)=(f1(-2)+f1(-4))*fb2(-5)
c      f(6)=(f1(1)+f1(3)+f1(5))*fb2(-5)
c      f(7)=(f2(-2)+f2(-4))*fb1(-5)
c      f(8)=(f2(1)+f2(3)+f2(5))*fb1(-5)
c      return

c  (u,c) (d,s,b) -> W+
      f(1)=(f1(2)+f1(4)) * (vtd2*fb1(1)+vts2*fb1(3)+vtb2*fb1(5))
c  (d~,s~,b~) (d,s,b) -> W+
      f(2)=((vud2+vcd2)*f1(-1)+(vus2+vcs2)*f1(-3)+(vub2+vcb2)*f1(-5))
     &     * (vtd2*fb1(1)+vts2*fb1(3)+vtb2*fb1(5))
c  b  u  -> W+
      f(3)=(f2(2)+f2(4)) * (vtd2*fb2(1)+vts2*fb2(3)+vtb2*fb2(5))
c  b d~ -> W+
      f(4)=((vud2+vcd2)*f2(-1)+(vus2+vcs2)*f2(-3)+(vub2+vcb2)*f2(-5))
     &     * (vtd2*fb2(1)+vts2*fb2(3)+vtb2*fb2(5))
c
c  u~ b~ -> W-
      f(5)=(f1(-2)+f1(-4)) * (vtd2*fb1(-1)+vts2*fb1(-3)+vtb2*fb1(-5))
c  d b~ -> W-
      f(6)=((vud2+vcd2)*f1(1)+(vus2+vcs2)*f1(3)+(vub2+vcb2)*f1(5))
     &     * (vtd2*fb1(-1)+vts2*fb1(-3)+vtb2*fb1(-5))
c  b~ u~ -> W-
      f(7)=(f2(-2)+f2(-4)) * (vtd2*fb2(-1)+vts2*fb2(-3)+vtb2*fb2(-5))
c  b~ d -> W-
      f(8)=((vud2+vcd2)*f2(1)+(vus2+vcs2)*f2(3)+(vub2+vcb2)*f2(5))
     &     * (vtd2*fb2(-1)+vts2*fb2(-3)+vtb2*fb2(-5))
c
      return
      end

      double precision function tpsi0(s,t,u,xm2,wm2)
c Returns the LO matrix element squared, averaged over spin/color
      implicit none
      double precision s,t,u,xm2,wm2
      tpsi0=0.25d0*s*(s-xm2)/(t-wm2)/(t-wm2)
      return
      end
