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

      call strcat(pref,'spt_j1.'//etachar//tend,fname)
      call bbook(3,fname,300,300d0,0d0)
      call strcat(pref,'seta_j1.'//etachar//tend,fname)
      call bbook(4,fname,120,6d0,-6d0)

      call strcat(pref,'spt_j2.'//etachar//tend,fname)
      call bbook(5,fname,300,300d0,0d0)
      call strcat(pref,'seta_j2.'//etachar//tend,fname)
      call bbook(6,fname,120,6d0,-6d0)

      call strcat(pref,'spt_bb.'//etachar//tend,fname)
      call bbook(7,fname,300,300d0,0d0)
      call strcat(pref,'seta_bb.'//etachar//tend,fname)
      call bbook(8,fname,120,6d0,-6d0)

      call strcat(pref,'spt_b1.'//etachar//tend,fname)
      call bbook(9,fname,300,300d0,0d0)
      call strcat(pref,'seta_b1.'//etachar//tend,fname)
      call bbook(10,fname,120,6d0,-6d0)

      call strcat(pref,'m_tb.'//etachar//tend,fname)
      call bbook(11,fname,500,650d0,150d0)

      call strcat(pref,'HT.'//etachar//tend,fname)
      call bbook(12,fname,500,500d0,0d0)

      call strcat(pref,'spt_b2.'//etachar//tend,fname)
      call bbook(13,fname,300,300d0,0d0)
      call strcat(pref,'seta_b2.'//etachar//tend,fname)
      call bbook(14,fname,120,6d0,-6d0)

      call strcat(pref,'spt_q1.'//etachar//tend,fname)
      call bbook(15,fname,300,300d0,0d0)
      call strcat(pref,'seta_q1.'//etachar//tend,fname)
      call bbook(16,fname,120,6d0,-6d0)

      call strcat(pref,'spt_q2.'//etachar//tend,fname)
      call bbook(17,fname,300,300d0,0d0)
      call strcat(pref,'seta_q2.'//etachar//tend,fname)
      call bbook(18,fname,120,6d0,-6d0)

      return
      end
c.
c. fill histograms
c.
      subroutine fhist(wgt,p)
      implicit double precision (a-h,o-z)
      double precision wgt,p(0:3,5)
      logical plot,passcut
      common/plots/plot
      common/to_cuts/passcut
      common/energy/s,ecm
      common/parm/hbarc2,gweak2,xm2,wm2
      real*8 PS(0:3,5)
      integer btag(2),njets,jet(3),pjet(3)
      common /graphing/ps,btag,njets,jet,pjet
      double precision xmtb,xmtbx,cthtb,gam,beta
      double precision cthcmtb,ptmp(0:3,3)
      real*8 etmin(3:8),etmax(3:8),etamin(3:8),etamax(3:8),dr(3:8,3:8)
      integer nocutflag,nev
      common /cuts/etmin,etmax,etamin,etamax,dr,nocutflag,nev

      double precision dot,pt,eta,rap,cosbeam,et
      double precision HT

      wtmp=wgt
      if(plot)then
c top
         call bfill(1,pt(ps(0,4)),wtmp)
         call bfill(2,eta(ps(0,4)),wtmp)
c j1
         call bfill(3,pt(ps(0,jet(1))),wtmp)
         call bfill(4,eta(ps(0,jet(1))),wtmp)
c Plot j2 if there
         if (jet(2).ne.0) then
            call bfill(5,pt(ps(0,jet(2))),wtmp)
            call bfill(6,eta(ps(0,jet(2))),wtmp)
            if (jet(2).eq.3) then
c plot b2
               call bfill(13,pt(ps(0,jet(2))),wtmp)
               call bfill(14,eta(ps(0,jet(2))),wtmp)
            else
c plot q2
               call bfill(17,pt(ps(0,jet(2))),wtmp)
               call bfill(18,eta(ps(0,jet(2))),wtmp)
            endif
         endif
c plot bb if seen
         if ((jet(1).eq.3).or.(jet(2).eq.3)) then
            call bfill(7,pt(ps(0,3)),wtmp)
            call bfill(8,eta(ps(0,3)),wtmp)
         endif
c plot b1
         if (jet(1).eq.3) then
            call bfill(9,pt(ps(0,3)),wtmp)
            call bfill(10,eta(ps(0,3)),wtmp)
         else
c plot q1
            call bfill(15,pt(ps(0,jet(1))),wtmp)
            call bfill(16,eta(ps(0,jet(1))),wtmp)
         endif
c plot mtb (if b seen)
         if ((jet(1).eq.3).or.(jet(2).eq.3)) then
            do i=0,3
               ptmp(i,1)=ps(i,4)+ps(i,3)
            enddo
            xmtb=dsqrt(dot(ptmp(0,1),ptmp(0,1)))
            call bfill(11,xmtb,wtmp)
         endif
c plot "HT"
         ht = et(ps(0,4))+et(ps(0,3))
         call bfill(12,ht,wtmp)

      endif

 1000 return
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
      call bwrite(7,fac)
      call bwrite(8,face)
      call bwrite(9,fac)
      call bwrite(10,face)
      call bwrite(11,fac)
      call bwrite(12,fac)
      if (norder.gt.0) then
         call bwrite(5,fac)
         call bwrite(6,face)
         call bwrite(13,fac)
         call bwrite(14,face)
         call bwrite(15,fac)
         call bwrite(16,face)
         call bwrite(17,fac)
         call bwrite(18,face)
      endif
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
