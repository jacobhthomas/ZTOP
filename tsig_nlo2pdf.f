c.
c. 2 -> 2 nlo cross section
c.
c Modified 9/17/02 to allow arbitrary scaling via isclf/isclfb, and renorm/fac
c  scales treated separately.
c Modified 6/18/24 to loop over multiple PDF sets
c
      subroutine tsig_nlo2(p,tot)
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
      double precision fb1t(-6:6),fb2t(-6:6)
      DOUBLE PRECISION XPQB1(-6:6),XPQB2(-6:6)
      double precision fb1(-6:6),fb2(-6:6)
      double precision xl(8),sig(8),xl1(8),xl2(8)
      double precision coll1,coll2,sftl,vrtl,subtot1
      double precision colh1,colh2,sfth,vrth,subtot2,vrtnbh(8)
      double precision collt1,collt2,sftlt,vrtlt
      double precision colht1,colht2,sftht,vrtht
      double precision collu1,collu2,sftlu,vrtlu
      double precision colhu1,colhu2,sfthu,vrthu
      double precision s,t,u,beta,lamt,omlt,lamu,omlu
      double precision lam,oml

      double precision scale1,scale2,scaleb1,scaleb2
      double precision rscale1,rscale2,rscaleb1,rscaleb2
      double precision asl1,asl2,ash1,ash2

      double precision totpdf(0:58)
      integer iloop

      integer i
C
C EXTERNAL FUNCTIONS
C
      DOUBLE PRECISION tpsi0,xlog,xlog2,dilog
      DOUBLE PRECISION dot,pt,alphas2
C
C GLOBAL VARIABLES
C
      double precision deltas,deltac
      common/slice/deltas,deltac
      double precision hbarc2,gweak2,xm2,wm2
      common/parm/hbarc2,gweak2,xm2,wm2
      integer      isclf,ialphas,iscale
      common/flags/isclf,ialphas,iscale
      integer       isclfb,iscaleb
      common/flagsb/isclfb,iscaleb
      integer       irsclf,irscale,irsclfb,irscaleb
      common/rflags/irsclf,irscale,irsclfb,irscaleb
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
      beta=1d0-xm2/s
      lamt=t/(t-xm2)
      lamu=u/(u-xm2)
      omlt=1d0-lamt
      omlu=1d0-lamu
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
      endif   ! iscale > 7
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
c Get PDFs ! MOVED below to enable loop over PDFs
c
c
c Collinear *** Set scales ***
      colh1=(2d0*xlog(deltas) + 1.5d0)*xlog(s/scaleb1**2) *4d0/3d0
      colh2=(2d0*xlog(deltas) + 1.5d0)*xlog(s/scaleb2**2) *4d0/3d0
      coll1=( 3.5d0 - pi*pi/3d0 -xlog2(deltas) -xlog2(beta) +2d0*     ! final
     &     xlog(deltas)*xlog(beta) -xlog(deltac)*(2d0*xlog(deltas)+
     &     1.5d0-2d0*xlog(beta)) )*4d0/3d0
      coll2= coll1+(2d0*xlog(deltas)+1.5d0)*xlog(s/scale2**2)*4d0/3d0 ! initial
      coll1= coll1+(2d0*xlog(deltas)+1.5d0)*xlog(s/scale1**2)*4d0/3d0 ! initial


c Virtual
      vrtlt=( -xlog2(s/-t) -3d0*xlog(s/-t) -8d0 -2d0*zeta2 )*4d0/3d0
      vrtht=( -0.5d0*xlog2(s/xm2) -2.5d0*xlog(s/xm2) -2d0*xlog(omlt)*
     &     xlog(s/xm2) -6d0 -xlog(omlt)/lamt -xlog2(omlt)
     &     -2d0*xlog(omlt) +2d0*dilog(lamt) -2d0*zeta2 )*4d0/3d0
      vrtlu=( -xlog2(s/-u) -3d0*xlog(s/-u) -8d0 -2d0*zeta2 )*4d0/3d0
      vrthu=( -0.5d0*xlog2(s/xm2) -2.5d0*xlog(s/xm2) -2d0*xlog(omlu)*
     &     xlog(s/xm2) -6d0 -xlog(omlu)/lamu -xlog2(omlu)
     &     -2d0*xlog(omlu) +2d0*dilog(lamu) -2d0*zeta2 )*4d0/3d0

c non-Born piece
c  (u,c) (d,s,b) -> W+
      vrtnbh(1)=s*u*xm2/t*xlog(omlt)/(t-wm2)/(t-wm2)/3d0    ! CF/4
c  (d~,s~,b~) (d,s,b) -> W+, s <-> u
      vrtnbh(2)=vrtnbh(1)
c  b  u  -> W+, t <-> u
      vrtnbh(3)=s*t*xm2/u*xlog(omlu)/(u-wm2)/(u-wm2)/3d0    ! CF/4
c  b  d~ -> W+, (3) s <-> t
      vrtnbh(4)=vrtnbh(3)
c
c  u~ b~ -> W-
      vrtnbh(5)=vrtnbh(1)
c  d  b~ -> W-
      vrtnbh(6)=vrtnbh(2)
c  b~ u~ -> W-
      vrtnbh(7)=vrtnbh(3)
c  b~ d  -> W-
      vrtnbh(8)=vrtnbh(4)

c Soft
      sftlt=( 4d0*xlog2(deltas) +4d0*xlog(deltas)*xlog(-t/s/beta)
     &     +xlog2(-t/s/beta) +2d0*dilog(1d0+t/s/beta) )*4d0/3d0
      sftht=( 2d0*xlog2(deltas) -2d0*xlog(deltas) +2d0*xlog(deltas)*
     &     xlog((xm2-t)*(xm2-t)/xm2/s) +(s+xm2)/(s-xm2)*xlog(s/xm2)
     &     +xlog2(xm2/(xm2-t)) -0.5d0*xlog2(s/xm2) +2d0*dilog(t/xm2)
     &     -2d0*dilog(u/(s+u)) )*4d0/3d0 ! CF=4/3
      sftlu=( 4d0*xlog2(deltas) +4d0*xlog(deltas)*xlog(-u/s/beta)
     &     +xlog2(-u/s/beta) +2d0*dilog(1d0+u/s/beta) )*4d0/3d0
      sfthu=( 2d0*xlog2(deltas) -2d0*xlog(deltas) +2d0*xlog(deltas)*
     &     xlog((xm2-u)*(xm2-u)/xm2/s) +(s+xm2)/(s-xm2)*xlog(s/xm2)
     &     +xlog2(xm2/(xm2-u)) -0.5d0*xlog2(s/xm2) +2d0*dilog(u/xm2)
     &     -2d0*dilog(t/(s+t)) )*4d0/3d0 ! CF=4/3

c Alpha_s
cz Added renormalization rscale different
      if(irscale.gt.7) then
c Do not choose rscale 8 unless cutting out 2 -> 2 events.
      if(irscale.eq.8) then
         rscale1=pt(p(0,3))               ! M_jj^2, but for 2 -> 2, M_jj=0
         rscale2=rscale1                   ! so use pt_d
      elseif (irscale.eq.9) then
         rscale1=dsqrt(-t)                ! -t = Q^2 out of p, Q^2+M_jj^2
         rscale2=dsqrt(-u)                ! -u = Q^2 out of pbar, M_jj=0
      elseif (irscale.eq.10) then
         rscale1=pt(p(0,3))               ! MT_jj^2, but for 2 -> 2,
         rscale2=rscale1                   !  MT_jj=pt_d
      elseif (irscale.eq.11) then
         rscale1=dsqrt(-t+pt(p(0,3)))     ! -t = Q^2 out of p, Q^2+MT_jj^2
         rscale2=dsqrt(-u+pt(p(0,3)))     ! -u = Q^2 out of pbar, MT_jj=pt_d
      elseif (irscale.eq.12) then
         rscale1=pt(p(0,3))+pt(p(0,4))                 ! HT = sum PT
         rscale2=rscale1
      elseif (irscale.eq.13) then
         rscale1=pt(p(0,3))+dsqrt((pt(p(0,4)))**2+xm2) ! HT~ = sum MT
         rscale2=rscale1
      endif
      else            ! normal rscale choices
      if(irscale.eq.4) then
         rscale1=dsqrt(xm2)
         rscale2=rscale1
      elseif (irscale.eq.5) then
         rscale1=dsqrt(s)
         rscale2=rscale1
      elseif (irscale.eq.6) then
         rscale1=dsqrt(-t)                ! -t = Q^2 out of p
         rscale2=dsqrt(-u)                ! -u = Q^2 out of pbar
      elseif (irscale.eq.7) then
         rscale1=dsqrt(-t+xm2)            ! -t = Q^2 out of p
         rscale2=dsqrt(-u+xm2)            ! -u = Q^2 out of pbar
      elseif (irscale.eq.0) then
         rscale1=pt(p(0,4))
         rscale2=rscale1
      elseif (irscale.eq.1) then
         rscale1=dsqrt((pt(p(0,4)))**2+xm2)
         rscale2=rscale1
      elseif (irscale.eq.2) then
         rscale1=pt(p(0,3))               ! Highest ET jet (only 1 here)
         rscale2=rscale1
      elseif (irscale.eq.3) then
         rscale1=dsqrt(xm2+2d0*dot(p(0,4),p(0,3))) ! 3 is highest ET jet
         rscale2=rscale1
      endif
      endif   ! irscale > 7
      if(irsclf.ne.0) then
         rscale1=rscale1*(dble(irsclf)/1d2)
         rscale2=rscale2*(dble(irsclf)/1d2)
      endif

      if(irscaleb.gt.7) then
      if(irscaleb.eq.8) then
         rscaleb1=dsqrt(xm2)              ! M_tj^2, but for 2 -> 2, M_tj=mt
         rscaleb2=rscaleb1
      elseif (irscaleb.eq.9) then
         rscaleb1=dsqrt(-t+xm2)           ! -t = Q^2 out of p, Q^2+M_tj^2
         rscaleb2=dsqrt(-u+xm2)           ! -u = Q^2 out of pbar, M_tj=mt
      elseif (irscaleb.eq.10) then
         rscaleb1=dsqrt((pt(p(0,4)))**2+xm2)  ! MT_tj^2, but for 2 -> 2,
         rscaleb2=rscaleb1                     !  MT_tj=MT_t
      elseif (irscaleb.eq.11) then
         rscaleb1=dsqrt(-t+(pt(p(0,4)))**2+xm2) ! -t=Q^2 out of p,Q^2+MT_tj^2
         rscaleb2=dsqrt(-u+(pt(p(0,4)))**2+xm2) ! -u=Q^2 out of pbar
      elseif (irscaleb.eq.12) then
         rscaleb1=pt(p(0,3))+pt(p(0,4))                 ! HT = sum PT
         rscaleb2=rscaleb1
      elseif (irscaleb.eq.13) then
         rscaleb1=pt(p(0,3))+dsqrt((pt(p(0,4)))**2+xm2) ! HT~ = sum MT
         rscaleb2=rscaleb1
      endif
      else            ! normal rscale choices
      if(irscaleb.eq.4) then
         rscaleb1=dsqrt(xm2)
         rscaleb2=rscaleb1
      elseif (irscaleb.eq.5) then
         rscaleb1=dsqrt(s)
         rscaleb2=rscaleb1
      elseif (irscaleb.eq.6) then
         rscaleb1=dsqrt(-t)                ! -t = Q^2 out of p
         rscaleb2=dsqrt(-u)                ! -u = Q^2 out of pbar
      elseif (irscaleb.eq.7) then
         rscaleb1=dsqrt(-t+xm2)            ! -t = Q^2 out of p
         rscaleb2=dsqrt(-u+xm2)            ! -u = Q^2 out of pbar
      elseif (irscaleb.eq.0) then
         rscaleb1=pt(p(0,4))
         rscaleb2=rscaleb1
      elseif (irscaleb.eq.1) then
         rscaleb1=dsqrt((pt(p(0,4)))**2+xm2)
         rscaleb2=rscaleb1
      elseif (irscaleb.eq.2) then
         rscaleb1=pt(p(0,3))               ! Highest ET jet (only 1 here)
         rscaleb2=rscaleb1
      elseif (irscaleb.eq.3) then
         rscaleb1=dsqrt(xm2+2d0*dot(p(0,4),p(0,3))) ! 3 is highest ET jet
         rscaleb2=rscaleb1
      endif
      endif
      if(irsclfb.ne.0) then
         rscaleb1=rscaleb1*(dble(irsclfb)/1d2)
         rscaleb2=rscaleb2*(dble(irsclfb)/1d2)
      endif

      if(rscale1.gt.1d0) then
         asl1=alphas2(rscale1)
      else
         asl1=0d0
      endif
      if(rscale2.gt.1d0) then
         asl2=alphas2(rscale2)
      else
         asl2=0d0
      endif
      if(rscaleb1.gt.1d0) then
         ash1=alphas2(rscaleb1)
      else
         ash1=0d0
      endif
      if(rscaleb2.gt.1d0) then
         ash2=alphas2(rscaleb2)
      else
         ash2=0d0
      endif

c LO
      sig(1)=tpsi0(s,t,u,xm2,wm2)       ! u  b -> W+
      sig(2)=tpsi0(u,t,s,xm2,wm2)       ! d~ b -> W+
      sig(3)=tpsi0(s,u,t,xm2,wm2)       ! b u  -> W+
      sig(4)=tpsi0(t,u,s,xm2,wm2)       ! b d~ -> W+
      sig(5)=sig(1)                     ! u~ b~ -> W-
      sig(6)=sig(2)                     ! d  b~ -> W-
      sig(7)=sig(3)                     ! b~ u~ -> W-
      sig(8)=sig(4)                     ! b~ d  -> W-

c
      subtot1=asl1*(sftlt+vrtlt+coll1)+ash1*(sftht+vrtht+colh1)
      subtot2=asl2*(sftlu+vrtlu+coll2)+ash2*(sfthu+vrthu+colh2)
cz
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
         call lumt3(xa,scale1,1,5,deltas,deltac,s,f1,f1t)
         call lumt3(xb,scaleb1,-1,5,deltas,deltac,s,fb1,fb1t)
         call lumt3(xb,scale2,-1,5,deltas,deltac,s,f2,f2t)
         call lumt3(xa,scaleb2,1,5,deltas,deltac,s,fb2,fb2t)
      elseif (icollide .eq. 3) then !PP
         do i=-6,6
            f1(i)=xpq1(i)
            f2(i)=xpq2(i)
            fb1(i)=xpqb1(i)
            fb2(i)=xpqb2(i)
         enddo
         call lumt3(xa,scale1,1,5,deltas,deltac,s,f1,f1t)
         call lumt3(xb,scaleb1,1,5,deltas,deltac,s,fb1,fb1t)
         call lumt3(xb,scale2,1,5,deltas,deltac,s,f2,f2t)
         call lumt3(xa,scaleb2,1,5,deltas,deltac,s,fb2,fb2t)
      endif
c
      call setlumt(f1,f2,fb1,fb2,xl)
      call setlumt(f1t,f2t,fb1,fb2,xl1)              ! I think this is right
      call setlumt(f1,f2,fb1t,fb2t,xl2)
c
      if(itop.eq.1) tot=tot
     &     +(xl(1)*subtot1+xl1(1)*asl1+xl2(1)*ash1)*sig(1)
     &     +ash1*xl(1)*vrtnbh(1)
     &     +(xl(2)*subtot1+xl1(2)*asl1+xl2(2)*ash1)*sig(2)
     &     +ash1*xl(2)*vrtnbh(2)
     &     +(xl(3)*subtot2+xl1(3)*asl2+xl2(3)*ash2)*sig(3)
     &     +ash2*xl(3)*vrtnbh(3)
     &     +(xl(4)*subtot2+xl1(4)*asl2+xl2(4)*ash2)*sig(4)
     &     +ash2*xl(4)*vrtnbh(4)
      if(iatop.eq.1)tot=tot
     &     +(xl(5)*subtot1+xl1(5)*asl1+xl2(5)*ash1)*sig(5)
     &     +ash1*xl(5)*vrtnbh(5)
     &     +(xl(6)*subtot1+xl1(6)*asl1+xl2(6)*ash1)*sig(6)
     &     +ash1*xl(6)*vrtnbh(6)
     &     +(xl(7)*subtot2+xl1(7)*asl2+xl2(7)*ash2)*sig(7)
     &     +ash2*xl(7)*vrtnbh(7)
     &     +(xl(8)*subtot2+xl1(8)*asl2+xl2(8)*ash2)*sig(8)
     &     +ash2*xl(8)*vrtnbh(8)
c
      tot=tot/2d0/pi       ! rest of as/2/pi
c
      totpdf(iloop) = tot
      enddo
cz Set lumi ratios here
      do i = 0, 58
         ration(i+1) = totpdf(i)/totpdf(0)
      enddo
      tot = totpdf(0)
cz //
      return
      end

c. log
      double precision function xlog(z)
      implicit none
      double precision z
      xlog=dlog(dabs(z))
      return
      end

c. log squared with negative argument
      double precision function xlog2(z)
      implicit none
      double precision z,pi
      parameter (pi=3.14159265358979d0)
      xlog2=(dlog(dabs(z)))**2
      if(z.lt.0d0)xlog2=xlog2-pi**2
      return
      end
