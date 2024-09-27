      program ztoptchan
c Copyright (C) 2001, by Zack Sullivan
c
c This program calculates t-channel single-top-quark production fully
c  differentially.  The b is taken massless, and spin is averaged over.
c
c This has been modified to adjust factorization and renormalization scales
c  separately.
c
      implicit double precision (a-h,o-z)
      character * 40 pref,fname,lpref
      character * 2 tend
      integer i,j,it1,it2

      common/flags/isclf,ialphas,iscale
      common/flagsb/isclfb,iscaleb
      common/rflags/irsclf,irscale,irsclfb,irscaleb
      common/parm/hbarc2,gweak2,xm2,wm2
      common/energy/s,ecm
      common/jet/rcone,rsep
      common/pdfzero/izerog,izeroq
      common/slice/deltas,deltac
      common/ckm/vud2,vus2,vub2,vcd2,vcs2,vcb2,vtd2,vts2,vtb2
      integer      itop,iatop
      common/ttbar/itop,iatop
      integer        norder
      common/histwrt/norder

      common/btagging/bwgt

      integer              icollide,ismear,igroup,iset,iq2
      common /to_collider2/icollide,ismear,igroup,iset,iq2
      real*8 etmin(3:8),etmax(3:8),etamin(3:8),etamax(3:8),dr(3:8,3:8)
      integer nocutflag,npart
      common /cuts/etmin,etmax,etamin,etamax,dr,nocutflag,npart

      external tlo2,tnlo2,tnlo3

c. write title
      write(*,*) ''
      write(*,198) '====== ======= ====== ======     ======='
      write(*,198) '    /     |    |    | |    |        |   '
      write(*,198) '   /      |    |    | |    |        |   '
      write(*,198) '  /       |    |    | |====  ----   |   '
      write(*,198) ' /        |    |    | |             |   '
      write(*,198) '======    |    ====== |             |   '
      write(*,*) ''
      write(*,198) 't-channel single-top-quark production  '
      write(*,198) 'Copyright 2001, by Zack Sullivan    '
 198  format(T20,A40)
      write(*,*) 'Version 1.1'
      write(*,*) ''
      write(*,*) 'Based on: 1. B.W. Harris, E. Laenen, L. Phaf, ',
     &     'Z. Sullivan, and S. Weinzierl,'
      write(*,*) '              Phys. Rev. D 66, 054024 (2002) ',
     &     '[arXiv:hep-ph/0207055].'
      write(*,*) '          2. Zack Sullivan, Phys. Rev. D 70, ',
     &     '114012 (2004)'
      write(*,*) '              [arXiv:hep-ph/0408049].'
      write(*,*) ''

c. input data
      read (*,*) norder      ! 0:lo   1:nlo   2:sum=lo+nlo
      read (*,*) iscale      ! 0:t pt 1:t transverse mass 2:j pt
c                        3:t-j invariant mass 4:mt 5:shat 6:Q^2 7:Q^2+mt^2
c   for active corrections ->  8:M_jj 9:Q^2+M_jj^2 10:MT_jj 11:Q^2+MT_jj^2
c isclf is also different
      read (*,*) isclf       ! 0:mu, all others are factors*100
      read (*,*) iscaleb     ! "  "
      read (*,*) isclfb      ! "  "
      read (*,*) irscale     ! "  " for renormalization scales_l
      read (*,*) irsclf      ! "  "
      read (*,*) irscaleb    ! "  " for renormalization scales_h
      read (*,*) irsclfb     ! "  "
      read (*,*) its1        ! vegas iterations LO 2 body
      read (*,*) npt1        ! vegas points     LO 2 body
      read (*,*) its2        ! vegas iterations NLO 2 body
      read (*,*) npt2        ! vegas points     NLO 2 body
      read (*,*) its3        ! vegas iterations NLO 3 body
      read (*,*) npt3        ! vegas points     NLO 3 body
      read (*,*) xm          ! top mass
      read (*,*) wm          ! W mass
      read (*,*) ecm         ! energy cm
      read (*,*) icollide    ! collider type
      read (*,*) iset        ! PDF set
      read (*,*) itop
      read (*,*) iatop
      read (*,*) rcone       ! cone size
      read (*,*) rsep        ! rsep
      read (*,*) nocutflag   ! 0: no cuts 1: detector
      read (*,*) ismear      ! 0: none 1: all 2: just b/jet
      read (*,*) id          ! 0: delta_s 1:delta_c 2:deltac=deltas/a
      read (*,*) imin
      read (*,*) imax
      read (*,*) afac
      read (*,*) pref        ! output file prefix string
      read (*,*) etmin(3),etmin(4),etmin(5)
      read (*,*) etamax(3),etamax(4),etamax(5)
      etmin(6)=5d0
      etamax(6)=4d0
c. output data
      write (*,*) 'Input parameters:'
      write (*,*) 'norder    = ',norder
      write (*,*) 'iscale_fl = ',iscale
      write (*,*) 'isclf     = ',isclf
      write (*,*) 'iscaleb_fh= ',iscaleb
      write (*,*) 'isclfb    = ',isclfb
      write (*,*) 'irscale_rl= ',irscale
      write (*,*) 'irsclf    = ',irsclf
      write (*,*) 'irscaleb_rh= ',irscaleb
      write (*,*) 'irsclfb   = ',irsclfb
      write (*,*) 'its1      = ',its1
      write (*,*) 'pts1      = ',npt1
      write (*,*) 'its2      = ',its2
      write (*,*) 'pts2      = ',npt2
      write (*,*) 'its3      = ',its3
      write (*,*) 'pts3      = ',npt3
      write (*,*) 'top mass  = ',xm
      write (*,*) 'w mass    = ',wm
      write (*,*) 'ecm       = ',ecm
      write (*,*) 'icollide  = ',icollide
      write (*,*) 'iset      = ',iset
      write (*,*) 'itop      = ',itop
      write (*,*) 'iatop     = ',iatop
      write (*,*) 'rcone     = ',rcone
      write (*,*) 'rsep      = ',rsep
      write (*,*) 'cutflag   = ',nocutflag
      write (*,*) 'smearing  = ',ismear
      write (*,*) 'id        = ',id
      write (*,*) 'imin      = ',imin
      write (*,*) 'imax      = ',imax
      write (*,*) 'afac      = ',afac
      write (*,*) 'pref      = ',pref


c. setup masses, energy, coupling, etc.
      xm2=xm*xm
      wm2=wm*wm

      call setcuts

      ecm=1000d0*ecm
      s=ecm*ecm
      write(*,*) ' '
      write(*,*) 'sqrt(S) = ',dsqrt(s)

      hbarc2 = 0.389379292d9   ! conversion to picobarns
c older value!!!      hbarc2 = 0.38937966d9   ! conversion to picobarns
      gweak2 = 8d0 * (1.16639d-5) * wm2 / dsqrt(2d0)
c      gweak2 = 8d0 * (1.16639d-5) * (80.4d0)**2 / dsqrt(2d0) ! Changed 3/4/02

c. CKM matrix
      vud=0.9751d0
      vud2=vud**2
      vus=0.2215d0
      vus2=vus**2
      vub=0.0035d0
      vub2=vub**2
      vcd=0.2210d0
      vcd2=vcd**2
      vcs=0.9743d0
      vcs2=vcs**2
      vcb=0.0410d0
      vcb2=vcb**2
c V_td = V_ts = 0, V_tb = 1  change here
c      vtd=0.009d0
      vtd=0d0
      vtd2=vtd**2
c      vts=0.04d0
      vts=0d0
      vts2=vts**2
      vtb=1d0
      vtb2=vtb**2
      write(*,*) ' '
      write(*,*) 'CKM matrix Vij'
      write(*,5) vud,vus,vub
      write(*,5) vcd,vcs,vcb
      write(*,5) vtd,vts,vtb
      write(*,*) ' '
 5    format(1x,3(2x,f6.4))


c. setup Lambda(nf=5) fron proton PDF
      if(iset.le.63) then
         it1=int(iset/10)
         it2=iset-10*it1
      else
         it1=6
         it2=iset-6000
      endif
      call setctq(it1,it2)

c. show what alphas this lambda gives at the Z mass
      zm=91.187d0
      zm2=zm*zm
      as2=alphas2(zm)
      write(*,*) 'alphas(m_z=91.187) = ',as2
      write(*,*) ' '

c. phase space slicing parameters

      if (nocutflag .eq. 0) tend = 'nc'
      if (nocutflag .eq. 1) tend = 'dc'
      if (nocutflag .eq. 2) tend = 'jc'
      if (nocutflag .eq. 3) tend = 'tc'

      call strcat(pref,'delta.out'//tend,fname)
      open(unit=11,file=fname,status='unknown')

      do 1000 i=imin,imax

         if(id.eq.0)then
            deltas=0.25d-5*2**i
            deltac=afac            ! 1d-05
         elseif(id.eq.1)then
            deltas=afac            ! 1d-2
            deltac=1d-6*2**i
         elseif(id.eq.2)then
            deltas = 0.25d-5*2**i
            deltac = deltas/afac
         else
            write(*,*) 'unknown id: ',id
            stop
         endif
      
         write(*,*) 'deltac = ',deltac
         write(*,*) 'deltas = ',deltas
         write(*,*) ' '

c. initialize random number generator
      call brm48i(40,0,0)
      write(*,*) ' '

c. set up histogram files
      j=istrl(pref)
      lpref=pref
      if (norder.eq.1) then
         lpref(j+1:j+1)=char(int(i/10)+48)
         lpref(j+2:j+2)=char(i-10*int(i/10)+48)
         lpref(j+3:j+3)='_'
      endif
      j=istrl(lpref)
      write (*,*) lpref
      call shist(lpref)


c. norder=0: lo
c.        1: nlo
c.        3: lo+nlo
c  tlo2, tnlo2,3 are vegas-like functions that call a function with momenta
c
      if(norder.eq.0)then

c. lo
         call vsup(4,npt1,its1,tlo2,aitlo2,sdtlo2,chi2tlo2)

         write(*,102) aitlo2,sdtlo2
         write(11,102) aitlo2,sdtlo2
 102     format(1x,e14.6,e11.4)
         tsum= aitlo2
         tdsum=sdtlo2
         call whist
         goto 1001

      elseif(norder.eq.1)then

c. nlo
         ialphas=2
         if(nocutflag.ne.1) then
            call vsup(4,npt2,its2,tnlo2,aitnlo2,sdtnlo2,chi2tnlo2)
         else
            aitnlo2=0d0
            sdtnlo2=0d0
            chi2tnlo2=0d0
         endif
         call vsup(7,npt3,its3,tnlo3,aitnlo3,sdtnlo3,chi2tnlo2)

         tsum=       aitnlo2   +aitnlo3
         tdsum=dsqrt(sdtnlo2**2+sdtnlo3**2)

         write(*,101)
     &deltas,deltac,aitnlo2,sdtnlo2,aitnlo3,sdtnlo3,tsum,tdsum
         write(11,101)
     &deltas,deltac,aitnlo2,sdtnlo2,aitnlo3,sdtnlo3,tsum,tdsum

c         goto 1000       !    ******* Remove in general *******

      elseif(norder.eq.2)then

c. lo+nlo
         ialphas=2
         if(nocutflag.ne.1) then
            call vsup(4,npt1,its1,tlo2,aitlo2,sdtlo2,chi2tlo2)
            call vsup(4,npt2,its2,tnlo2,aitnlo2,sdtnlo2,chi2tnlo2)
         else
            aitlo2=0d0
            sdtlo2=0d0
            chi2tlo2=0d0
            aitnlo2=0d0
            sdtnlo2=0d0
            chi2tnlo2=0d0
         endif
         aitnlo3=0d0
         sdtnlo3=0d0
         chi2tnlo3=0d0
         call vsup(7,npt3,its3,tnlo3,aitnlo3,sdtnlo3,chi2tnlo3)

         tsum=       aitlo2   +aitnlo2   +aitnlo3
         tdsum=dsqrt(sdtlo2**2+sdtnlo2**2+sdtnlo3**2)

         write(*,103)
     &deltas,deltac,aitlo2,sdtlo2,aitnlo2,sdtnlo2,aitnlo3,sdtnlo3,
     &        tsum,tdsum
         write(11,103)
     &deltas,deltac,aitlo2,sdtlo2,aitnlo2,sdtnlo2,aitnlo3,sdtnlo3,
     &        tsum,tdsum

      else
        write(*,*) 'unknown norder: ',norder
        stop
      endif

c. write out histogram files
      call whist

 1000 continue

 101  format(1x,2(e11.4),3(e14.6,e11.4))
 103  format(1x,2(e11.4),4(e14.6,e11.4))
 104  format(A15,E14.6,A4,E12.4,A3)

 1001 close(11)

      write(*,104) 'Cross section: ',tsum,'+-',tdsum,'pb'

      end
