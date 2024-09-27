      double precision function fnlo3(xx,vwgt)
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

      double precision a,as,at,gmin,gmax,xat,y,ycm,tau,cth,sth,sphi
      double precision cth1,sth1,cth2,sth2,theta,theta1,theta2,cphi
      double precision s12,beta,srs2,s45,smin45,smax45,p4z,p14,p24
      double precision p30,p33,p40,p43,b,gamma,b2,gamma2,p40_12,p4z_12
      double precision xjac,xnorm,wgt,sig,ps3
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
c. theta
      theta = pi * xx(3)
      cth = dcos( theta )
      sth = dsin( theta )
      xjac = xjac * pi
c. theta1
      theta1 = pi * xx(4)
      cth1 = dcos( theta1 )
      sth1 = dsin( theta1 )
      xjac = xjac * pi
c. theta2
      theta2 = 2d0 * pi * xx(5)
      cth2 = dcos( theta2 )
      sth2 = dsin( theta2 )
      xjac = xjac * 2d0 * pi
c. s45
      smin45 = xm2
      smax45 = s12
      s45 = smin45 + ( smax45 - smin45 ) * xx(6)
      xjac = xjac * (smax45 - smin45)
c. arbitrary overall angle phi
      sphi = dsin(2d0*pi*xx(7))
      cphi = dcos(2d0*pi*xx(7))
      xjac = xjac*2d0*pi
c. put together PS3
      ps3=xjac*sth*sth1*(s45-xm2)*(s12-s45)/s45/s12/1024d0/pi**4
      srs2=0.5d0*ecm
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
      b = ( xa - xb ) / ( xa + xb )
      gamma = 1d0 / dsqrt( 1d0 - b**2 )
      p30 = 0.5d0 * ( s12 - s45 ) / sqrt( s12 )
      p33 = p30
      p(0,3) = gamma * ( p30 + b * p33 * cth )
      p(1,3) = p33 * sth * cphi
      p(2,3) = p33 * sth * sphi
      p(3,3) = gamma * ( p33 * cth + b * p30 )
      b2 = ( s45 - s12 ) / ( s45 + s12 )
      gamma2 = 1d0 / sqrt( 1d0 - b2**2 )
      p40 = 0.5d0 * ( s45 + xm2 ) / sqrt( s45 )
      p4z = 0.5d0 * ( s45 - xm2 ) / sqrt( s45 ) * cth1
      p40_12 = gamma2 * ( p40 + b2 * p4z )
      p4z_12 = gamma2 * ( p4z + b2 * p40 )
      p14 = 0.5d0*(s45-xm2)*sth1*cth2*cth/sqrt(s45)+p4z_12*sth
      p24 = 0.5d0*(s45-xm2)*sth1*sth2/sqrt(s45)
      p(1,4) = p14*cphi - p24*sphi
      p(2,4) = p14*sphi + p24*cphi
      p4z_12 = -0.5d0*(s45-xm2)*sth1*cth2*sth/sqrt(s45)+p4z_12*cth
      p(0,4) = gamma * ( p40_12 + b * p4z_12 )
      p(3,4) = gamma * ( p4z_12 + b * p40_12 )
      p(0,5)=p(0,1)+p(0,2)-p(0,3)-p(0,4)
      p(1,5)=p(1,1)+p(1,2)-p(1,3)-p(1,4)
      p(2,5)=p(2,1)+p(2,2)-p(2,3)-p(2,4)
      p(3,5)=p(3,1)+p(3,2)-p(3,3)-p(3,4)
c
c ***** Insert pre-selection here *****
c
      njets=3
      if (.not.passcuts(p)) then
         fnlo3=0d0
         return
      endif
c
c Set scale, pdf's
c
      if(iscale.eq.0) scale2=(pt(p(0,4)))**2
      if(iscale.eq.1) scale2=(pt(p(0,4)))**2+xm2
      if(iscale.eq.2) scale2=(pt(p(0,3)))**2
cz      if(iscale.eq.3) scale2=xm2+2d0*dot(p(0,4),p(0,3))
cz Changed 6/20/02 to always use the b-jet (similar to Q^2, but not quite).
      if(iscale.eq.3) scale2=xm2+2d0*dot(ps(0,4),ps(0,3))
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
      call ssig_nlo3(p,sig)
      xnorm =-0.25d0*hbarc2*gweak2**2/s12  ! phi above multiplies this by 2*pi
      wgt=ps3*xnorm*sig*as
      call fhist(wgt*vwgt,p)
      fnlo3=wgt
      return
      end
