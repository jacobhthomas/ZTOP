      subroutine smearing(p,ps)
c************************************************************************
c     Simulates the detector resolution by smearing the energy 
c     according to a gaussian, and then reassigning the momentum so that
c     the direction and mass of the jet remains constant.
c************************************************************************
      implicit none
c
c     Constants
c
      double precision width
      parameter       (width = 0.8d0)
c
c     Arguments
c
      double precision p(0:3)      !Input 4 momenta
      double precision ps(0:3)     !Smeared 4 momenta
c
c     Local
c     
      integer j
      integer idum
      double precision fac
c      logical good
c
c     Global
c
c
c     External
c
      double precision dot
      real   gasdev
      data idum/-1/
      save idum
c-----
c  Begin Code
c-----
c
      fac=(1d0+dsqrt(width**2/(p(0))+.05**2)*gasdev(idum))
c
c     Smear the momentum
c         
      do j=1,3
         ps(j)=p(j)*fac
      enddo
c
c     Now set the mass
c
      ps(0) = dsqrt(dot(p(0),p(0))
     &           +ps(1)**2+ps(2)**2+ps(3)**2)
      return
      end

      function gasdev(idum)
      data iset/0/
      save iset,gset
      if (iset.eq.0) then
1       v1=2.*zran1(idum)-1.
        v2=2.*zran1(idum)-1.
        r=v1**2+v2**2
        if(r.ge.1.)go to 1
        fac=sqrt(-2.*log(r)/r)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      end

      function zran1(idum)
      dimension r(97)
      parameter (m1=259200,ia1=7141,ic1=54773,rm1=3.8580247e-6)
      parameter (m2=134456,ia2=8121,ic2=28411,rm2=7.4373773e-6)
      parameter (m3=243000,ia3=4561,ic3=51349)
      data iff /0/
      save iff,ix1,ix2,ix3,r
      if (idum.lt.0.or.iff.eq.0) then
        iff=1
        ix1=mod(ic1-idum,m1)
        ix1=mod(ia1*ix1+ic1,m1)
        ix2=mod(ix1,m2)
        ix1=mod(ia1*ix1+ic1,m1)
        ix3=mod(ix1,m3)
        do 11 j=1,97
          ix1=mod(ia1*ix1+ic1,m1)
          ix2=mod(ia2*ix2+ic2,m2)
          r(j)=(float(ix1)+float(ix2)*rm2)*rm1
11      continue
        idum=1
      endif
      ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ia2*ix2+ic2,m2)
      ix3=mod(ia3*ix3+ic3,m3)
      j=1+(97*ix3)/m3
      if(j.gt.97.or.j.lt.1)pause 'ran1'
      zran1=r(j)
      r(j)=(float(ix1)+float(ix2)*rm2)*rm1
      return
      end
