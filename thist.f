c. ------------------- all histogramming stuff is here ----------------
c.
c. call histogram initializer
c. make changes to limits, etc. here
c.
      subroutine shist(pref)
      character*40 pref,fname
      character*2 tend
      character*1 etachar

      double precision s,ecm
      common/energy/s,ecm
      integer              icollide,ismear,igroup,iset,iq2
      common /to_collider2/icollide,ismear,igroup,iset,iq2
      real*8 etmin(3:8),etmax(3:8),etamin(3:8),etamax(3:8),dr(3:8,3:8)
      integer nocutflag,npart
      common /cuts/etmin,etmax,etamin,etamax,dr,nocutflag,npart

      etachar = char(int(icollide)+48)

      if (nocutflag .eq. 0) tend = 'nc'
      if (nocutflag .eq. 1) tend = 'dc'
      if (nocutflag .eq. 2) tend = 'jc'
      if (nocutflag .eq. 3) tend = 'tc'

      call ihist

      call strcat(pref,'spt_t.'//etachar//tend,fname)
      call bbook(1,fname,300,300d0,0d0)
      call strcat(pref,'seta_t.'//etachar//tend,fname)
      call bbook(2,fname,120,6d0,-6d0)

      call strcat(pref,'spt_f.'//etachar//tend,fname)
      call bbook(3,fname,300,300d0,0d0)
      call strcat(pref,'seta_f.'//etachar//tend,fname)
      call bbook(4,fname,120,6d0,-6d0)

      call strcat(pref,'spt_j1.'//etachar//tend,fname)
      call bbook(5,fname,300,300d0,0d0)
      call strcat(pref,'seta_j1.'//etachar//tend,fname)
      call bbook(6,fname,120,6d0,-6d0)

      call strcat(pref,'spt_j2.'//etachar//tend,fname)
      call bbook(7,fname,300,300d0,0d0)
      call strcat(pref,'seta_j2.'//etachar//tend,fname)
      call bbook(8,fname,120,6d0,-6d0)

      call strcat(pref,'spt_bb.'//etachar//tend,fname)
      call bbook(9,fname,300,300d0,0d0)
      call strcat(pref,'seta_bb.'//etachar//tend,fname)
      call bbook(10,fname,120,6d0,-6d0)

      call strcat(pref,'spt_b1.'//etachar//tend,fname)
      call bbook(11,fname,300,300d0,0d0)
      call strcat(pref,'seta_b1.'//etachar//tend,fname)
      call bbook(12,fname,120,6d0,-6d0)

      call strcat(pref,'drtf.'//etachar//tend,fname)
      call bbook(13,fname,100,10d0,0d0)
      call strcat(pref,'drtj1.'//etachar//tend,fname)
      call bbook(14,fname,100,10d0,0d0)
      call strcat(pref,'drtj2.'//etachar//tend,fname)
      call bbook(15,fname,100,10d0,0d0)
      call strcat(pref,'drtjb.'//etachar//tend,fname)
      call bbook(16,fname,100,10d0,0d0)

      call strcat(pref,'dphitf.'//etachar//tend,fname)
      call bbook(17,fname,180,180d0,0d0)
      call strcat(pref,'dphitj1.'//etachar//tend,fname)
      call bbook(18,fname,180,180d0,0d0)
      call strcat(pref,'dphitj2.'//etachar//tend,fname)
      call bbook(19,fname,180,180d0,0d0)
      call strcat(pref,'dphitjb.'//etachar//tend,fname)
      call bbook(20,fname,180,180d0,0d0)

      call strcat(pref,'detatf.'//etachar//tend,fname)
      call bbook(21,fname,80,8d0,0d0)
      call strcat(pref,'detatj1.'//etachar//tend,fname)
      call bbook(22,fname,80,8d0,0d0)
      call strcat(pref,'detatj2.'//etachar//tend,fname)
      call bbook(23,fname,80,8d0,0d0)
      call strcat(pref,'detatjb.'//etachar//tend,fname)
      call bbook(24,fname,80,8d0,0d0)

      call strcat(pref,'spt_b2.'//etachar//tend,fname)
      call bbook(25,fname,300,300d0,0d0)
      call strcat(pref,'seta_b2.'//etachar//tend,fname)
      call bbook(26,fname,120,6d0,-6d0)

      call strcat(pref,'spt_q1.'//etachar//tend,fname)
      call bbook(27,fname,300,300d0,0d0)
      call strcat(pref,'seta_q1.'//etachar//tend,fname)
      call bbook(28,fname,120,6d0,-6d0)

      call strcat(pref,'spt_q2.'//etachar//tend,fname)
      call bbook(29,fname,300,300d0,0d0)
      call strcat(pref,'seta_q2.'//etachar//tend,fname)
      call bbook(30,fname,120,6d0,-6d0)

c
      return
      end

c.
c. fill histograms
c.
      subroutine fhist(wgt,p)
      implicit double precision (a-h,o-z)
      double precision  pi              , to_deg
      parameter        (pi = 3.1415927d0, to_deg=180d0/pi)
      double precision wgt,p(0:3,5)
      logical plot,passcut
      common/plots/plot
      common/to_cuts/passcut
      common/energy/s,ecm
      common/parm/hbarc2,gweak2,xm2,wm2
      real*8 PS(0:3,5)
      integer btag(2),njets,jet(3),pjet(3)
      common /graphing/ps,btag,njets,jet,pjet
      double precision xmtb,xmtbx,cthtb,gam,beta,cthtj
      double precision xa,xb,scale
      common/toevent/  xa,xb,scale
      real*8 etmin(3:8),etmax(3:8),etamin(3:8),etamax(3:8),dr(3:8,3:8)
      integer nocutflag,npart
      common /cuts/etmin,etmax,etamin,etamax,dr,nocutflag,npart

      double precision bwgt
      common/btagging/ bwgt

      double precision dot,pt,eta,rap,cosbeam,r2,delta_phi

      double precision tet3,teta3,tet5,teta5
      logical passplot

      wtmp=wgt
      if(plot)then
cz Cuts for plotting were moved here.
         passplot=.false.
         if (nocutflag.eq.0) passplot=.true.

         tet3=pt(ps(0,3))
         teta3=eta(ps(0,3))
         tet5=pt(ps(0,5))
         teta5=eta(ps(0,5))
         if ((nocutflag.eq.1).and.(njets.eq.3)) then  ! need njets=3
            if(tet3.ge.tet5) then
               if((tet3.ge.etmin(3)).and.(abs(teta3).le.etamax(3)).and.
     &              (tet5.ge.etmin(5)).and.(abs(teta5).le.etamax(5)))
     &              passplot=.true.
               jet(1)=3   ! irrelevant if previous false
               jet(2)=5
            else
               if((tet5.ge.etmin(3)).and.(abs(teta5).le.etamax(3)).and.
     &              (tet3.ge.etmin(5)).and.(abs(teta3).le.etamax(5)))
     &              passplot=.true.
               jet(1)=5   ! irrelevant if previous false
               jet(2)=3
            endif
         endif

         if (nocutflag.ge.2) then
            jet(1)=0
            jet(2)=0
            if((tet3.ge.etmin(3)).and.(abs(teta3).le.etamax(3)))
     &           jet(1)=3
            if (njets.eq.3) then
               if((tet5.ge.etmin(3)).and.(abs(teta5).le.etamax(3)))then
                  if ((jet(1).eq.0).or.((jet(1).eq.3).and.
     &                 (tet5.gt.tet3))) then
                     jet(2)=jet(1)
                     jet(1)=5
                  else
                     jet(2)=5
                  endif
               endif
               if ((jet(1).eq.3).and.(tet5.ge.etmin(6)).and.
     &              (abs(teta5).le.etamax(6))) jet(2)=5
               if ((jet(1).eq.5).and.(tet3.ge.etmin(6)).and.
     &              (abs(teta3).le.etamax(6))) jet(2)=3
            endif
            if ((nocutflag.eq.2).and.(jet(1).ne.0)) passplot=.true.
            if ((nocutflag.eq.3).and.(jet(1).ne.0).and.
     &           ((jet(2).eq.0).or.(pt(ps(0,jet(2))).lt.etmin(5)).or.
     &         (abs(eta(ps(0,jet(2)))).gt.etamax(5)))) passplot=.true.
         endif



         if(passplot) then
c top
            call bfill(1,pt(ps(0,4)),wtmp)
            call bfill(2,eta(ps(0,4)),wtmp)

c Plot "LO d" jet always
            call bfill(3,pt(ps(0,3)),wtmp)
            call bfill(4,eta(ps(0,3)),wtmp)
            call bfill(13,dsqrt(r2(ps(0,4),ps(0,3))),wtmp)
            call bfill(17,to_deg*delta_phi(ps(0,4),ps(0,3)),wtmp)
            call bfill(21,abs(eta(ps(0,4))-eta(ps(0,3))),wtmp)
c Plot j1 always
            call bfill(5,pt(ps(0,jet(1))),wtmp)
            call bfill(6,eta(ps(0,jet(1))),wtmp)
            call bfill(14,dsqrt(r2(ps(0,4),ps(0,jet(1)))),wtmp)
            call bfill(18,to_deg*delta_phi(ps(0,4),ps(0,jet(1))),wtmp)
            call bfill(22,abs(eta(ps(0,4))-eta(ps(0,jet(1)))),wtmp)
c Plot j2 if there
            if (jet(2).ne.0) then
               call bfill(7,pt(ps(0,jet(2))),wtmp)
               call bfill(8,eta(ps(0,jet(2))),wtmp)
               call bfill(15,dsqrt(r2(ps(0,4),ps(0,jet(2)))),wtmp)
             call bfill(19,to_deg*delta_phi(ps(0,4),ps(0,jet(2))),wtmp)
               call bfill(23,abs(eta(ps(0,4))-eta(ps(0,jet(2)))),wtmp)
            endif
c Plot b~ if there
         if (bwgt.gt.1d-10) then
            if(njets.eq.2) then        ! b~ absorbed into jet
               if((pt(ps(0,3)).ge.etmin(6)).and.(abs(eta(ps(0,3))).le.
     &              etamax(6))) then
               call bfill(9,pt(ps(0,3)),wtmp*bwgt)
               call bfill(10,eta(ps(0,3)),wtmp*bwgt)
               call bfill(11,pt(ps(0,3)),wtmp*bwgt) ! b in highest pt jet
               call bfill(12,eta(ps(0,3)),wtmp*bwgt)
               call bfill(16,dsqrt(r2(ps(0,4),ps(0,3))),wtmp*bwgt)
               call bfill(20,to_deg*delta_phi(ps(0,4),ps(0,3)),
     &              wtmp*bwgt)
               call bfill(24,abs(eta(ps(0,4))-eta(ps(0,3))),wtmp*bwgt)
               endif
               call bfill(27,pt(ps(0,jet(1))),wtmp*(1d0-bwgt))     ! q1
               call bfill(28,eta(ps(0,jet(1))),wtmp*(1d0-bwgt))
            else
               if((pt(ps(0,5)).ge.etmin(6)).and.(abs(eta(ps(0,5))).le.
     &              etamax(6))) then
               call bfill(9,pt(ps(0,5)),wtmp*bwgt)
               call bfill(10,eta(ps(0,5)),wtmp*bwgt)
               call bfill(16,dsqrt(r2(ps(0,4),ps(0,5))),wtmp*bwgt)
               call bfill(20,to_deg*delta_phi(ps(0,4),ps(0,5)),
     &              wtmp*bwgt)
               call bfill(24,abs(eta(ps(0,4))-eta(ps(0,5))),wtmp*bwgt)
               endif
               if (jet(1).eq.5) then
                  call bfill(11,pt(ps(0,5)),wtmp*bwgt) ! b IN highest pt jet
                  call bfill(12,eta(ps(0,5)),wtmp*bwgt)
                  call bfill(27,pt(ps(0,jet(1))),wtmp*(1d0-bwgt)) ! q1
                  call bfill(28,eta(ps(0,jet(1))),wtmp*(1d0-bwgt))
                  if(jet(2).ne.0) then
                     call bfill(29,pt(ps(0,jet(2))),wtmp)         ! q2
                     call bfill(30,eta(ps(0,jet(2))),wtmp)
                  endif
               elseif (jet(2).eq.5) then
                  call bfill(25,pt(ps(0,5)),wtmp*bwgt) ! b in 2nd highest jet
                  call bfill(26,eta(ps(0,5)),wtmp*bwgt)
                  call bfill(29,pt(ps(0,5)),wtmp*(1d0-bwgt))      ! q2
                  call bfill(30,eta(ps(0,5)),wtmp*(1d0-bwgt))
                  call bfill(27,pt(ps(0,jet(1))),wtmp)            ! q1
                  call bfill(28,eta(ps(0,jet(1))),wtmp)
               else                                       ! No jet 5 seen.
                  call bfill(27,pt(ps(0,jet(1))),wtmp)    ! q1
                  call bfill(28,eta(ps(0,jet(1))),wtmp)
                  if(jet(2).ne.0) then                    ! Doesn't ever occur
                     call bfill(29,pt(ps(0,jet(2))),wtmp) ! q2
                     call bfill(30,eta(ps(0,jet(2))),wtmp)
                  endif
               endif
            endif
         else
c no b~
            call bfill(27,pt(ps(0,jet(1))),wtmp)       ! q1
            call bfill(28,eta(ps(0,jet(1))),wtmp)
            if(jet(2).ne.0) then
               call bfill(29,pt(ps(0,jet(2))),wtmp)    ! q2
               call bfill(30,eta(ps(0,jet(2))),wtmp)
            endif
         endif
         endif

      endif

      return
      end

c.
c. write histograms to their respective files
c.
      subroutine whist
      implicit double precision (a-h,o-z)
      integer        norder
      common/histwrt/norder
      fac=1d0
      face=1d0             ! this histogrammer automatically adjusts
      call bwrite(1,fac)
      call bwrite(2,face)
      call bwrite(3,fac)
      call bwrite(4,face)
      call bwrite(5,fac)
      call bwrite(6,face)
      call bwrite(13,face)
      call bwrite(14,face)
      call bwrite(17,fac)
      call bwrite(18,fac)
      call bwrite(21,face)
      call bwrite(22,face)
      call bwrite(27,fac)
      call bwrite(28,face)
      if(norder.ne.0) then
         call bwrite(7,fac)
         call bwrite(8,face)
         call bwrite(9,fac)
         call bwrite(10,face)
         call bwrite(11,fac)
         call bwrite(12,face)
         call bwrite(15,face)
         call bwrite(16,face)
         call bwrite(19,fac)
         call bwrite(20,fac)
         call bwrite(23,face)
         call bwrite(24,face)
         call bwrite(25,fac)
         call bwrite(26,face)
         call bwrite(29,fac)
         call bwrite(30,face)
      endif
c
      return
      end

      double precision function cosbeam(p)
      implicit none
c Calculates cos_theta(p) wrt. beam
      double precision p(0:3),denom
c
      denom = dsqrt(p(1)*p(1)+p(2)*p(2)+p(3)*p(3))
      if (denom.gt.0) then
         cosbeam=p(3)/denom
      else
         cosbeam=0d0
      endif
      return
      end
