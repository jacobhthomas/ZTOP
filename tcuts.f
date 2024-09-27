cz Changed 2/26/04 to avoid boundary problems, place final cuts in histogram.
cz Turns out to be irrelevant, but leaving it to avoid adding typos.
      SUBROUTINE SETCUTS
      IMPLICIT NONE
c
c     Constants
c
c
c     Local
c
      integer i,j
C     
C     GLOBAL
C
c      real*8         pmass(8)
c      common/to_mass/pmass
c      double precision      spole(maxinvar),swidth(maxinvar),bwjac
c      common/to_brietwigner/spole          ,swidth          ,bwjac
c      logical             pp
c      common/to_collider/ pp

      real*8 etmin(3:8),etmax(3:8),etamin(3:8),etamax(3:8),dr(3:8,3:8)
      integer nocutflag,nev
      common /cuts/etmin,etmax,etamin,etamax,dr,nocutflag,nev
      integer qte
      common /to_smearing/ qte

      integer              icollide,ismear,igroup,iset,iq2
      common /to_collider2/icollide,ismear,igroup,iset,iq2
c
c     Data
c
      data qte/1/     ! W -> q q_     W from t
c     *** user inputs ***
c     Particle               #3      4     5     6     7    8
      data etmin            /20d0 ,20d0, 20d0 , 0d0,  0d0, 0d0/
      data etmax            /0d0  , 0d0,  0d0 , 0d0,  0d0, 0d0/
      data etamin           /0d0  , 0d0,  0d0 , 0d0,  0d0, 0d0/
      data etamax           /2d0  , 2d0,  2d0 , 0d0,  0d0, 0d0/
c
      data (dr(i,3),i=4,7)  /      .7d0 ,.7d0 , 0d0 , 0d0/
      data (dr(i,4),i=5,7)  /            .7d0 , 0d0 , 0d0/
      data (dr(i,5),i=6,7)  /                   0d0 , 0d0/
      data (dr(i,6),i=7,7)  /                         0d0/
c
C-----
C  BEGIN CODE
C-----
c
c     Write out smearing info.
      if (ismear .eq. 0) then
         write(*,*) 'No smearing.'
      elseif(ismear .eq. 1) then
         write(*,*)'Smearing jets, but using correct top momentum.'
      elseif(ismear .eq. 2) then
         write(*,*)'Smearing jets and top momentum.'
      else
         write(*,*) 'Unknown smearing value, turning smearing off.'
     &        ,ismear
         ismear=0
      endif
c
      RETURN
      END


      LOGICAL FUNCTION PASSCUTS(PX)
C**************************************************************************
C     INPUT:
C            P(0:3,1)           MOMENTUM OF INCOMING PARTON
C            P(0:3,2)           MOMENTUM OF INCOMING PARTON
C            P(0:3,3)           MOMENTUM OF d
C            P(0:3,4)           MOMENTUM OF t
C            P(0:3,5)           MOMENTUM OF extra jet
C            COMMON/JETCUTS/   CUTS ON JETS
C     OUTPUT:
C            TRUE IF EVENTS PASSES ALL CUTS LISTED
C**************************************************************************
      IMPLICIT NONE
C
C     ARGUMENTS
C
c      REAL*8 PX(0:3,5)
      double precision px(0:3,5)

C
C     LOCAL
C
      LOGICAL FIRSTTIME
      integer i,j,nparticles,method
      real*8 tet3,tet4,tet5
      integer seej(4),seej3,seej4
C
C     EXTERNAL
C
      REAL*8 R2,DOT,ET,ETA,DJ,pt

C
C     GLOBAL
C
      double precision rcone,rsep
      common /jet/     rcone,rsep
      real*8 etmin(3:8),etmax(3:8),etamin(3:8),etamax(3:8),dr(3:8,3:8)
      integer nocutflag,npart
      common /cuts/etmin,etmax,etamin,etamax,dr,nocutflag,npart
      real*8 P(0:3,5)
      integer btag(2),njets,jet(3),pjet(3)
      common /graphing/p,btag,njets,jet,pjet
      integer              icollide,ismear,igroup,iset,iq2
      common /to_collider2/icollide,ismear,igroup,iset,iq2
      integer qte
      common /to_smearing/ qte
      DATA FIRSTTIME/.TRUE./
C-----
C  BEGIN CODE
C-----
      PASSCUTS=.TRUE.             !EVENT IS OK UNLESS OTHERWISE CHANGED
      nparticles=5
      IF (FIRSTTIME) THEN
         FIRSTTIME=.FALSE.
         if ((nocutflag .lt. 0) .or. (nocutflag .gt. 3)) then
            write (*,*) 'CUT FLAG MUST BE (0-3)'
            stop
         endif
         if (nocutflag .eq. 0) then
            write(*,'(a)') '*** NO CUTS ***'
         else
            write(*,'(a10,6i8)') 'Particle',(i,i=3,nparticles)
            write(*,'(a10,6f8.1)') 'Et >',(etmin(i),i=3,nparticles)
            write(*,'(a10,6f8.1)') 'Et <',(etmax(i),i=3,nparticles)
            write(*,'(a10,6f8.1)') 'Eta >',(etamin(i),i=3,nparticles)
            write(*,'(a10,6f8.1)') 'Eta <',(etamax(i),i=3,nparticles)
            do j=3,nparticles
               write(*,'(a10,6f8.1)') 'd R >',(-0.0,i=3,j),
     &              (dr(i,j),i=j+1,nparticles)
            enddo
            if (nocutflag .eq. 1) write(*,'(a)') '** See j1 and j2 **'
            if (nocutflag .eq. 2) write(*,'(a)') '** A least 1 jet **'
            if (nocutflag .eq. 3) write(*,'(a)') '** Only 1 jet **'
         endif
      ENDIF
c
      method=3        ! k_T algorithm with rcone=rsep, pt
      call jetfinder(px,p,rcone,rsep,njets,method) ! if njets=2, just returns p
c
      if (ismear.ne.0) then
         call smearing(p(0,3),p(0,3))                    ! jet
         if (ismear.eq.2) call smearing(p(0,4),p(0,4))   ! top
         if (njets.eq.3) call smearing(p(0,5),p(0,5))    ! jet
      endif
c
      btag(1) = 0
      btag(2) = 0
c
cz Changed to pt for paper 6/4/02.
      tet3 = pt(p(0,3))
      tet4 = pt(p(0,4))
      tet5 = pt(p(0,5))

      if(njets.eq.2) then
         jet(1)=3
         jet(2)=0
      else        ! njets = 3
         if(tet3.ge.tet5) then
            jet(1)=3
            jet(2)=5
         else
            jet(1)=5
            jet(2)=3
         endif
      endif

      if (nocutflag .eq. 0) return

cz Changed for comparison w/ madgraph 8/19/02
      jet(1)=0
      jet(2)=0
cz Added for 8/20/04 release, first check to see if top quark passes cuts
      if(etamax(4).gt.1d-5) then
         if((tet4.lt.etmin(4)).or.(abs(eta(p(0,4))).gt.etamax(4))) then
            passcuts=.false.
            return
         endif
      endif
cz Changed 9/25/02 to get multiple cut options.  Jet cuts are pt ordered.
cz   Cuts on jet 2 cannot be tighter than cuts on jet 1.
      seej(1)=0
      seej(2)=0
      seej(3)=0
      seej(4)=0
      if((tet3.ge.etmin(3)).and.
     &     (abs(eta(p(0,3))).le.etamax(3))) seej(1)=1 ! see 3 over cut 1
      if((njets.eq.3).and.(tet5.ge.etmin(3)).and.
     &     (abs(eta(p(0,5))).le.etamax(3))) seej(2)=1 ! see 5 over cut 1
      if((tet3.ge.etmin(5)).and.
     &     (abs(eta(p(0,3))).le.etamax(5))) seej(3)=1 ! see 3 over cut 2
      if((njets.eq.3).and.(tet5.ge.etmin(5)).and.
     &     (abs(eta(p(0,5))).le.etamax(5))) seej(4)=1 ! see 5 over cut 2
c See both
      if (nocutflag.eq.1) then
         if(njets.eq.2) then
            passcuts=.false.
            return
         endif
cz Changed 2/26/04 to avoid boundary problems, place final cuts in histogram.
         if((tet3.lt.(etmin(6)-1d0)).or.(abs(eta(p(0,3))).gt.
     &        (etamax(6)+0.25d0)).or.(tet5.lt.(etmin(6)-1d0)).or.
     &        (abs(eta(p(0,5))).gt.(etamax(6)+0.25d0))) then
            passcuts=.false.
            return
         endif
cz NOTE this must be redone in histogramming stage.
         if(tet3.ge.tet5) then
            jet(1)=3
            jet(2)=5
         else
            jet(1)=5
            jet(2)=3
         endif
         return
      endif

cz At least or exactly 1 jet here are the same cut:
      if((nocutflag.eq.2).or.(nocutflag.eq.3)) then
         if(((tet3.lt.(etmin(6)-1d0)).or.(abs(eta(p(0,3))).gt.
     &        (etamax(6)+0.25d0))).and.((tet5.lt.(etmin(6)-1d0)).or.
     &        (abs(eta(p(0,5))).gt.(etamax(6)+0.25d0)))) then
            passcuts=.false.
            return
         endif
cz NOTE this must be redone in histogramming stage.
         if(tet3.ge.tet5) then
            jet(1)=3
            jet(2)=5
         else
            jet(1)=5
            jet(2)=3
         endif
         return
      endif

      return
      END
