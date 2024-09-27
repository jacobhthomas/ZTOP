      program ztopschan
c Copyright (C) 2001, by Zack Sullivan
c
c This program calculates s-channel single-top-quark production fully
c  differentially.  The b is taken massless, and spin is averaged over.
c
      implicit double precision (a-h,o-z)
      character * 40 pref,fname,lpref,lhaset
      character * 2 tend
      integer i,j,it1,it2

      common/flags/isclf,ialphas,iscale
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

      integer              icollide,ismear,igroup,iset,iq2
      common /to_collider2/icollide,ismear,igroup,iset,iq2
      real*8 etmin(3:8),etmax(3:8),etamin(3:8),etamax(3:8),dr(3:8,3:8)
      integer nocutflag,npart
      common /cuts/etmin,etmax,etamin,etamax,dr,nocutflag,npart

      external flo2,fnlo2,fnlo3

c. write title
      write(*,*) ''
      write(*,198) '====== ======= ====== ======      ======'
      write(*,198) '    /     |    |    | |    |      |     '
      write(*,198) '   /      |    |    | |    |      ======'
      write(*,198) '  /       |    |    | |====  ----      |'
      write(*,198) ' /        |    |    | |                |'
      write(*,198) '======    |    ====== |           ======'
      write(*,*) ''
      write(*,198) 's-channel single-top-quark production  '
      write(*,198) 'Copyright 2001, by Zack Sullivan    '
 198  format(T20,A40)
      write(*,*) 'Version 1.2'
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
      read (*,*) iscale      ! 0:t pt 1:t transverse mass 2:b pt
c                              3:t-b invariant mass 4:mt 5:shat 6:mw
      read (*,*) isclf       ! 0:mu   1:2*mu   2:mu/2
      read (*,*) its1        ! vegas iterations LO 2 body
      read (*,*) npt1        ! vegas points     LO 2 body
      read (*,*) its2        ! vegas iterations NLO 2 body
      read (*,*) npt2        ! vegas points     NLO 2 body
      read (*,*) its3        ! vegas iterations NLO 3 body
      read (*,*) npt3        ! vegas points     NLO 3 body
      read (*,*) xm          ! top mass (GeV)
      read (*,*) wm          ! W mass (GeV)
      read (*,*) ecm         ! energy cm (TeV)
      read (*,*) icollide    ! 1,2: p-pbar  3: p-p collider type
c      read (*,*) iset        ! PDF set
      read (*,*) lhaset,iset      ! PDF set
      read (*,*) itop        ! 1: top production
      read (*,*) iatop       ! 1: anti-top production
      read (*,*) rcone       ! cone size
      read (*,*) rsep        ! rsep
      read (*,*) nocutflag   ! 0: no cuts 1: 2j 2: >=1 j 3: 1j (nc, dc, jc, tc)
      read (*,*) ismear      ! 0: no smearing 1: smear jets 2: smear jets+top
      read (*,*) id          ! 0: delta_s 1:delta_c 2:deltac=deltas/a
      read (*,*) imin        ! imin
      read (*,*) imax        ! imax
      read (*,*) afac        ! deltac = deltas/a
      read (*,*) pref        ! output file prefix string
      read (*,*) etmin(3),etmin(4),etmin(5)     ! b, t, j_2 ptmin cuts
      read (*,*) etamax(3),etamax(4),etamax(5)  ! b, t, j_2 etamax cuts
      etmin(6)=5d0
      etamax(6)=4d0
c. output data
      write (*,*) 'Input parameters:'
      write (*,*) 'norder    = ',norder
      write (*,*) 'iscale    = ',iscale
      write (*,*) 'isclf     = ',isclf
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
c      write (*,*) 'iset      = ',iset
      write (*,*) 'lhaset iset  = ',lhaset, iset
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

      hbarc2 = 0.389379292d9   ! conversion to picobarns 2000
c older value!!!      hbarc2 = 0.38937966d9   ! conversion to picobarns
      gweak2 = 8d0 * (1.16639d-5) * wm2 / dsqrt(2d0)

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
c      vtd=0.009d0
      vtd=0d0
      vtd2=vtd**2
c     vts=0.04d0
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
c      if(iset.le.63) then
c         it1=int(iset/10)
c         it2=iset-10*it1
c      else
c         it1=6
c         it2=iset-6000
c      endif
c      call setctq(it1,it2)
      call InitPDFsetByName(lhaset)
      call InitPDF(iset)


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
            deltas=0.25d-4*2**i
            deltac=1d-05
         elseif(id.eq.1)then
            deltas=1d-2
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
      write (*,*) lpref
      call shist(lpref)


c. norder=0: lo
c.        1: nlo
c.        3: lo+nlo
c  flo2, fnlo2,3 are vegas-like functions that call a function with momenta
c
      if(norder.eq.0)then

c. lo
         call vsup(4,npt1,its1,flo2,ailo2,sdlo2,chi2lo2)

         write(*,102) ailo2,sdlo2
         write(11,102) ailo2,sdlo2
 102     format(1x,e14.6,e11.4)
         sum= ailo2
         dsum=sdlo2
         call whist
         goto 1001

      elseif(norder.eq.1)then

c. nlo
         ialphas=2
         if(nocutflag.ne.1) then
            call vsup(4,npt2,its2,fnlo2,ainlo2,sdnlo2,chi2nlo2)
         else
            ainlo2=0d0
            sdnlo2=0d0
            chi2nlo2=0d0
         endif
         call vsup(7,npt3,its3,fnlo3,ainlo3,sdnlo3,chi2nlo2)

         sum=       ainlo2   +ainlo3
         dsum=dsqrt(sdnlo2**2+sdnlo3**2)

         write(*,101)
     &deltas,deltac,ainlo2,sdnlo2,ainlo3,sdnlo3,sum,dsum
         write(11,101)
     &deltas,deltac,ainlo2,sdnlo2,ainlo3,sdnlo3,sum,dsum

      elseif(norder.eq.2)then

c. lo+nlo
         ialphas=2
         if(nocutflag.ne.1) then
            call vsup(4,npt1,its1,flo2,ailo2,sdlo2,chi2lo2)
            call vsup(4,npt2,its2,fnlo2,ainlo2,sdnlo2,chi2nlo2)
         else
            ailo2=0d0
            sdlo2=0d0
            chi2lo2=0d0
            ainlo2=0d0
            sdnlo2=0d0
            chi2nlo2=0d0
         endif
         call vsup(7,npt3,its3,fnlo3,ainlo3,sdnlo3,chi2nlo2)

         sum=       ailo2   +ainlo2   +ainlo3
         dsum=dsqrt(sdlo2**2+sdnlo2**2+sdnlo3**2)

         write(*,103)
     &deltas,deltac,ailo2,sdlo2,ainlo2,sdnlo2,ainlo3,sdnlo3,sum,dsum
         write(11,103)
     &deltas,deltac,ailo2,sdlo2,ainlo2,sdnlo2,ainlo3,sdnlo3,sum,dsum

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

      write(*,104) 'Cross section: ',sum,'+-',dsum,'pb'

      end
