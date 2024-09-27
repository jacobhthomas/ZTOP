      function pqg0(z)
      implicit double precision (a-h,o-z)
      omz=1d0-z
      pqg0=0.5d0*(z*z+omz*omz)
      return
      end

      function pqg1(z)
      implicit double precision (a-h,o-z)
      pqg1=z*(z-1d0)
      return
      end

      function pgg0(z)
      implicit double precision (a-h,o-z)
      omz=1d0-z
      pgg0=6d0*(z/omz+omz/z+z*omz)
      return
      end

      function pgg1(z)
      implicit double precision (a-h,o-z)
      pgg1=0d0
      return
      end

      function pqq0(z)
      implicit double precision (a-h,o-z)
      pqq0=4d0*(1d0+z*z)/(3d0*(1d0-z))
      return
      end

      function pqq1(z)
      implicit double precision (a-h,o-z)
      pqq1=4d0*(z-1d0)/3d0
      return
      end

      function pgq0(z)
      implicit double precision (a-h,o-z)
      pgq0=pqq0(1d0-z)
      return
      end

      function pgq1(z)
      implicit double precision (a-h,o-z)
      pgq1=pqq1(1d0-z)
      return
      end

c. dot product
      function dot(p,q)
      implicit double precision (a-h,o-z)
      dimension p(0:3),q(0:3)
      dot=p(0)*q(0)-p(1)*q(1)-p(2)*q(2)-p(3)*q(3)
      return
      end

c. concatenate str1 and str2 into str
      subroutine strcat(str1,str2,str)
      character *(*) str1,str2,str
      l1=istrl(str1)
      l2=istrl(str2)
      l =len(str)
      if(l.lt.l1+l2) then
          write(*,*) 'error: l1+l2>l in strcat'
          write(*,*) 'l1=',l1,' str1=',str1
          write(*,*) 'l2=',l2,' str2=',str2
          write(*,*) 'l=',l
          stop
      endif
      if(l1.ne.0) str(1:l1)=str1(1:l1)
      if(l2.ne.0) str(l1+1:l1+l2)=str2(1:l2)
      if(l1+l2+1.le.l) str(l1+l2+1:l)= ' '
      end
c. returns the position of the last non-blank character in string
      function istrl(string)
      character * (*) string
      i = len(string)
      dowhile(i.gt.0.and.string(i:i).eq.' ')
         i=i-1
      enddo
      istrl = i
      end


      DOUBLE PRECISION FUNCTION R2(P1,P2)
c************************************************************************
c     Distance in eta,phi between two particles.
c************************************************************************
      IMPLICIT NONE
c
c     Arguments
c
      double precision p1(0:3),p2(0:3)
c
c     External
c
      double precision ETA,DELTA_PHI
      external eta,delta_phi
c-----
c  Begin Code
c-----
      R2 = (DELTA_PHI(P1,P2))**2+(ETA(p1)-ETA(p2))**2
      RETURN
      END

      DOUBLE PRECISION  FUNCTION eta(p)
c************************************************************************
c     Returns pseudo-rapidity of particle
c************************************************************************
      IMPLICIT NONE
c
c     Arguments
c
      double precision  p(0:3)
c
c     Local
c
      double precision pm
c
c     Global
c
      double precision s,msti,mu,lampp(3),als,x1,x2
      common/to_event/ s,msti,mu,lampp,als,x1,x2
c-----
c  Begin Code
c-----
      pm=dsqrt(p(1)**2+p(2)**2+p(3)**2)
      eta = .5d0*dlog((pm+p(3))/(pm-p(3)))
      end

      DOUBLE PRECISION  FUNCTION etalab(p)
c************************************************************************
c     Returns rapidity of particle in the lab frame (CM frame is rest)
c************************************************************************
      IMPLICIT NONE
c
c     Arguments
c
      double precision  p(0:3)
c
c     Local
c
      double precision pm
c
c     Global
c
      REAL*8              S,X1,X2,PSWGT,JAC
      COMMON /PHASESPACE/ S,X1,X2,PSWGT,JAC
c-----
c  Begin Code
c-----
      pm=dsqrt(p(1)**2+p(2)**2+p(3)**2)
      etalab = .5d0*dlog((pm+p(3))/(pm-p(3)))+.5*dlog(x1/x2) ! CM-> lab
      end

      DOUBLE PRECISION  FUNCTION rap(p)
c************************************************************************
c     Returns rapidity of particle in frame of p
c************************************************************************
      IMPLICIT NONE
c
c     Arguments
c
      double precision  p(0:3)
c
c     Local
c
      double precision pm
c-----
c  Begin Code
c-----
      if ((p(0)-p(3)).gt.1d-8) then
         rap = .5d0*dlog((p(0)+p(3))/(p(0)-p(3)))
      else
         rap = 999d0
      endif
      end

      DOUBLE PRECISION FUNCTION DELTA_PHI(P1, P2)
c************************************************************************
c     Returns separation in phi of two particles p1,p2
c************************************************************************
      IMPLICIT NONE
c
c     Arguments
c
      double precision p1(0:3),p2(0:3)
c
c     Local
c
      REAL*8 DENOM, TEMP
c-----
c  Begin Code
c-----
      DENOM = SQRT(P1(1)**2 + P1(2)**2) * SQRT(P2(1)**2 + P2(2)**2)
      TEMP = MAX(-0.99999999D0, (P1(1)*P2(1) + P1(2)*P2(2)) / DENOM)
      TEMP = MIN( 0.99999999D0, TEMP)
      DELTA_PHI = ACOS(TEMP)
      END

      double precision function pt(p)
c************************************************************************
c     Returns transverse momentum of particle
c************************************************************************
      IMPLICIT NONE
c
c     Arguments
c
      double precision p(0:3)
c-----
c  Begin Code
c-----
      pt = dsqrt(p(1)**2+p(2)**2)
      end

      double precision function et(p)
c************************************************************************
c     Returns transverse energy of particle
c************************************************************************
      IMPLICIT NONE
c
c     Arguments
c
      double precision p(0:3)
c
c     Local
c
      double precision pt
c-----
c  Begin Code
c-----
      pt = dsqrt(p(1)**2+p(2)**2)
      if (pt .gt. 0d0) then
         et = p(0)*pt/dsqrt(pt**2+p(3)**2)
      else
         et = 0d0
      endif
      end

      double precision function etlab(p)
c************************************************************************
c     Returns transverse energy of particle in LAB frame (P in CM frame)
c************************************************************************
      IMPLICIT NONE
c
c     Arguments
c
      double precision p(0:3)
c
c     Local
c
      double precision pt,beta,gamma,ep,pz
c
c     Global
c
      double precision s,msti,mu,lampp(3),als,x1,x2
      common/to_event/ s,msti,mu,lampp,als,x1,x2
c-----
c  Begin Code
c-----
      pt = dsqrt(p(1)**2+p(2)**2)
      beta = (x1-x2)/(x1+x2)
      gamma = 0.5d0*(x1+x2)/dsqrt(x1*x2)
      ep = gamma*(p(0)+beta*p(3))
      pz = gamma*(p(3)+beta*p(0))
      if (pt .gt. 0d0) then
         etlab = ep*pt/dsqrt(pt**2+pz**2)
      else
         etlab = 0d0
      endif
      end

      double precision function etcm(p)
c************************************************************************
c     Returns transverse energy of particle in CM frame (P in LAB frame)
c************************************************************************
      IMPLICIT NONE
c
c     Arguments
c
      double precision p(0:3)
c
c     Local
c
      double precision pt,beta,gamma,ep,pz
c
c     Global
c
      double precision s,msti,mu,lampp(3),als,x1,x2
      common/to_event/ s,msti,mu,lampp,als,x1,x2
c-----
c  Begin Code
c-----
      pt = dsqrt(p(1)**2+p(2)**2)
      beta = (x1-x2)/(x1+x2)
      gamma = 0.5d0*(x1+x2)/dsqrt(x1*x2)
      ep = gamma*(p(0)-beta*p(3))
      pz = gamma*(p(3)-beta*p(0))
      if (pt .gt. 0d0) then
         etcm = ep*pt/dsqrt(pt**2+pz**2)
      else
         etcm = 0d0
      endif
      end

      double precision function DJ(p1,p2)
c***************************************************************************
c     Uses Durham algorythm to calculate the y value for two partons
c***************************************************************************
      implicit none
c
c     Arguments
c
      double precision p1(0:3),p2(0:3)
c
c     Local
c
      double precision pm1,pm2,costh
      integer j
c-----
c  Begin Code
c-----
      pm1 = dsqrt(p1(1)**2+p1(2)**2+p1(3)**2)
      pm2 = dsqrt(p2(1)**2+p2(2)**2+p2(3)**2)
      if (pm1*pm2 .ne. 0d0) then
         costh = (p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3))/(pm1*pm2)
         dj = 2d0*min(p1(0)**2,p2(0)**2)*(1d0-costh)   !Durham
c         dj = 2d0*p1(0)*p2(0)*(1d0-costh)    !JADE
      else
         print*,'Warning 0 momentum in Durham algorythm'
         write(*,'(4e15.5)') (p1(j),j=0,3)
         write(*,'(4e15.5)') (p2(j),j=0,3)
         dj = 0d0
      endif
      end

      double precision function DJ2(p1,p2)
c***************************************************************************
c     Uses Lorentz
c***************************************************************************
      implicit none
c
c     Arguments
c
      double precision p1(0:3),p2(0:3)
c
c     Local
c
      integer j
c
c     External
c
      double precision dot
c-----
c  Begin Code
c-----
      dj2 = dot(p1,p1)+2d0*dot(p1,p2)+dot(p2,p2)
      return
      end


      double precision function phi(p)
c************************************************************************
c     Returns angle phi of momentum p, phi=0 == p_T=p_x
c************************************************************************
      IMPLICIT NONE
C  
C     Constants
C  
      double precision pi
      parameter       (pi=3.14159265358979d0)
      double precision pi2
      parameter       (pi2=1.57079632679490d0)
      double precision pit2
      parameter       (pit2=6.28318530717959d0)
c
c     Arguments
c
      double precision p(0:3)
c-----
c  Begin Code
c-----
      if(p(1).gt.0d0) then
         phi=datan(p(2)/p(1))
      elseif(p(1).lt.0d0) then
         phi=datan(p(2)/p(1))+pi
      else      ! p(1)=0
         phi=dsign(pi2,p(2))
      endif
      if(phi.lt.0d0)phi=phi+pit2
      return
      end
