      double precision function fnlo2(xx,vwgt)
c
c This function takes x and wgt values from Vegas, generates the momenta
c  and phase space.
c
      implicit none
C  
C CONSTANTS
C  
      double precision pi
      parameter       (pi=3.14159265358979d0)
C  
C ARGUMENTS 
C  
      double precision xx(10),vwgt
C
C LOCAL VARIABLES 
C
      double precision p(0:3,5),scale2

      double precision a,as,at,gmin,gmax,xat,y,ycm,tau,cth,sth
      double precision s12,beta,t,u,srs2,sphi,cphi
      double precision xjac,xnorm,wgt,sig
      integer i
C
C EXTERNAL FUNCTIONS
C
      LOGICAL PASSCUTS
      DOUBLE PRECISION ALPHAS2,pt,dot
C
C GLOBAL VARIABLES
C
      double precision deltas,deltac
      common/slice/deltas,deltac
      double precision hbarc2,gweak2,xm2,wm2
      common/parm/hbarc2,gweak2,xm2,wm2
      integer      isclf,ialphas,iscale
      common/flags/isclf,ialphas,iscale
      double precision s,ecm
      common/energy/s,ecm
      integer      itop,iatop
      common/ttbar/itop,iatop

      double precision xa,xb,scale
      common/toevent/  xa,xb,scale

      integer              icollide,ismear,igroup,iset,iq2
      common /to_collider2/icollide,ismear,igroup,iset,iq2
      real*8 PS(0:3,5)
      integer btag(2),njets,jet(3),pjet(3)
      common /graphing/ps,btag,njets,jet,pjet
c
c Begin code
c
c
c Define integration variables, momenta
c
c. tau integration
      a = 2d0
      at = 1d0 - a
      xat = 1d0 / at
      gmin = xat * (xm2/s)**at
      gmax = xat
      y = gmin + xx(1) * ( gmax - gmin )
      tau = ( y * at )**xat
      xjac = ( gmax - gmin ) * tau**a
c. ycm integration
      gmax = - 0.5d0 * dlog( tau )
      gmin = - gmax
      ycm = gmin + xx(2) * ( gmax - gmin )
      xjac = xjac * ( gmax - gmin )
c. xa, xb, s, and beta may now be calculated
      xa = dsqrt( tau ) * dexp( ycm )
      xb = tau / xa
      s12 = tau * s
      beta = 1d0 - xm2 / s12
c. partonic cms angle integration
      gmin = -1d0
      gmax = 1d0
      cth = gmin + xx(3) * ( gmax - gmin )
      xjac = xjac * ( gmax - gmin )
c.  t, u, and sin(theta) may now be calculated
      t = - 0.5d0 * s12 * beta * ( 1d0 - cth )
      u = xm2 - s12 - t
      sth = dsqrt( 1d0 - cth * cth )
      srs2=0.5d0*ecm
c. arbitrary overall angle phi
      sphi = dsin(2d0*pi*xx(4))
      cphi = dcos(2d0*pi*xx(4))
      xjac = xjac*2d0*pi
c. incoming parton 4-vectors
      p(0,1)=srs2*xa
      p(1,1)=0d0
      p(2,1)=0d0
      p(3,1)=p(0,1)
      p(0,2)=srs2*xb
      p(1,2)=0d0
      p(2,2)=0d0
      p(3,2)=-p(0,2)
c. outgoing parton 4-vectors
      p(0,3)=-0.5d0 * ( xa * u + xb * t ) / tau / ecm
      p(1,3)=-0.5d0*dsqrt(s12)*beta * sth * sphi
      p(2,3)=0.5d0*dsqrt(s12)*beta * sth * cphi
      p(3,3)=-0.5d0 * ( xa * u - xb * t ) / tau / ecm
      p(0,4)=p(0,1)+p(0,2)-p(0,3)
      p(1,4)=-p(1,3)
      p(2,4)=-p(2,3)
      p(3,4)=p(3,1)+p(3,2)-p(3,3)
      p(0,5)=0d0
      p(1,5)=0d0
      p(2,5)=0d0
      p(3,5)=0d0
c
c ***** Insert pre-selection here *****
c
      njets=2
      if (.not.passcuts(p)) then
         fnlo2=0d0
         return
      endif
c
c Set scale, pdf's
c
      if(iscale.eq.0) scale2=(pt(p(0,4)))**2
      if(iscale.eq.1) scale2=(pt(p(0,4)))**2+xm2
      if(iscale.eq.2) scale2=(pt(p(0,3)))**2
      if(iscale.eq.3) scale2=xm2+2d0*dot(p(0,4),p(0,3))
      if(iscale.eq.4) scale2=xm2
      if(iscale.eq.5) scale2=s12
      if(iscale.eq.6) scale2=wm2    ! To test vs. PYTHIA
      scale = dsqrt(scale2)
      if(isclf.eq.1) scale=scale*2d0
      if(isclf.eq.2) scale=scale/2d0

      as=alphas2(scale)
c
c Call matrix element
c
      call ssig_nlo2(p,sig)
      xnorm=hbarc2*beta*gweak2**2/64d0/pi**2/s12*xjac
      wgt=xnorm*sig*(as/2d0/pi)
      call fhist(wgt*vwgt,p)
      fnlo2=wgt
      return
      end
