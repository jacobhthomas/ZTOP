      subroutine jetfinder(pin,p,rcone,rsep,njets,method)
      implicit none
c
c NOTE: This is specialized for a maximum njets=3, and ONLY checks whether
c        the gluon is close to the b, since the b is massless, methods 1,2=3,4.
c
c This subroutine combines 'njets' momenta 'pin' within a size 'rcone,rsep'
c  and outputs new 'p' with a final 'njets'.
c
c Input:  pin(0:3,njets), rcone, rsep, method, njets  - last one is gluon
c
c Output: p(0:3,njets), njets'
c
c method: 1) k_T alogrithm with rcone == rsep, E_T used
c         2) k_T alogrithm with rcone, rsep both used, E_T used
c         3) k_T alogrithm with rcone == rsep, p_T used
c         4) k_T alogrithm with rcone, rsep both used, p_T used
c
      double precision emin
      parameter       (emin=1d-10)
c
c Arguments
c
      integer njets,method
      double precision rcone,rsep
      double precision pin(0:3,5),p(0:3,5)
c
c Local
c
      integer i,j,k
c      double precision di(njets),dij(njets,njets)
      double precision di(3),dij2(3,3),rc2,rs2,rm12,rm13,rm23
      double precision dm1,dm2,dm3,etas,phis,ets
c
c External
c
      double precision R2,ET,PT,ETA,PHI,RAP
c
c Begin code
c
      if ((njets.eq.2).or.(pin(0,5).lt.emin)) then  ! do nothing, or gluon soft
         njets=2
         do i=1,4
            do j=0,3
               p(j,i)=pin(j,i)
            enddo
         enddo
         return
      else                          ! return 1,2 same always
         do j=0,3
            p(j,1)=pin(j,1)
            p(j,2)=pin(j,2)
         enddo
      endif
c
      if ((method.eq.1).or.(method.eq.2)) then
         di(1)=et(pin(0,3))
         di(2)=et(pin(0,4))
         di(3)=et(pin(0,5))
      elseif ((method.eq.3).or.(method.eq.4)) then
         di(1)=pt(pin(0,3))
         di(2)=pt(pin(0,4))
         di(3)=pt(pin(0,5))
      else
 999     write (*,*) 'Method must be (1-4), used:',method
         stop
      endif
      rc2=rcone*rcone
      rs2=rsep*rsep
c comparing gluon (5) to top (4) or bottom (3) subtract 2 to get dij(n)
      dij2(1,3) = r2(pin(0,5),pin(0,3))
c      dij2(2,3) = r2(pin(0,5),pin(0,4))
c      dij2(1,2) = r2(pin(0,4),pin(0,3))
c
c 1) k_T alogrithm with rcone == rsep, E_T used
c 3) k_T alogrithm with rcone == rsep, p_T used
      if ((method.eq.1).or.(method.eq.3)) then
         if(dij2(1,3).lt.rc2) goto 100 ! (bg)
c         if((dij2(1,3).lt.rc2).and.(dij2(1,3).lt.dij2(2,3))) goto 100 ! (bg)
c         if((dij2(2,3).lt.rc2).and.(dij2(2,3).lt.dij2(1,3))) goto 200 ! (tg)
         goto 300                                                     ! nj=3
c 2) k_T alogrithm with rcone, rsep both used, E_T used
c 4) k_T alogrithm with rcone, rsep both used, p_T used
      elseif ((method.eq.2).or.(method.eq.4)) then
         rm13=dmin1(rsep,rcone*(di(1)+di(3))/dmax1(di(1),di(3)))
c         rm23=dmin1(rsep,rcone*(di(2)+di(3))/dmax1(di(2),di(3)))
c         rm12=dmin1(rsep,rcone*(di(1)+di(2))/dmax1(di(1),di(2)))
         if(dij2(1,3).lt.rm13) goto 100     ! (bg)
c         if((dij2(1,3).lt.rm13).and.(dij2(2,3).ge.rm23)) goto 100     ! (bg)
c         if((dij2(1,3).ge.rm13).and.(dij2(2,3).lt.rm23)) goto 200     ! (tg)
c         if((dij2(1,3).lt.rm13).and.(dij2(2,3).lt.rm23)) then
c            dm1=dmin1(di(1),di(2))*dij2(1,2)/rm12
c            dm2=dmin1(di(1),di(3))*dij2(1,3)/rm13
c            dm3=dmin1(di(2),di(3))*dij2(2,3)/rm23
c            if(dm2.le.dm3) goto 100                                   ! (bg)
c            goto 200                                                  ! (tg)
c         endif
         goto 300                                                     ! nj=3
      else
         goto 999    ! whoops
      endif
      write (*,*) 'Logic bug in routine jetfinder.'
      stop

c njets=2 (bg) combined in k_T alogrithm, t passed through
 100  njets=2
      do j=0,3
         p(j,4)=pin(j,4)               ! t untouched
      enddo
      ets=di(1)+di(3)
      if (ets.lt.emin) then
         write (*,*) 'ERROR: Using unweighted average in jetfinder.'
         etas=0.5d0*(rap(pin(0,3))+rap(pin(0,5)))
c         etas=0.5d0*(eta(pin(0,3))+eta(pin(0,5)))
         phis=0.5d0*(phi(pin(0,3))+phi(pin(0,5)))
         goto 110
      endif
      etas=(di(1)*rap(pin(0,3))+di(3)*rap(pin(0,5)))/ets
c      etas=(di(1)*eta(pin(0,3))+di(3)*eta(pin(0,5)))/ets
      phis=(di(1)*phi(pin(0,3))+di(3)*phi(pin(0,5)))/ets

c For now ets==pts, no masses - et==pt, eta==rap

 110  p(0,3)=ets*dcosh(etas)      ! with mass, ets -> dsqrt(pts**2+m^2)
      p(1,3)=ets*dcos(phis)       ! with mass, ets -> pts
      p(2,3)=ets*dsin(phis)       ! with mass, ets -> pts
      p(3,3)=ets*dsinh(etas)      ! with mass, ets -> dsqrt(pts**2+m^2)

      return

c njets=2 (tg) combined in k_T alogrithm, b passed through
 200  njets=2
c
      goto 300                         ! **** tg disabled for now
c
      do j=0,3
         p(j,3)=pin(j,3)               ! b untouched
      enddo
      ets=di(2)+di(3)
      if (ets.lt.emin) then
         write (*,*) 'ERROR: Using unweighted average in jetfinder.'
         etas=0.5d0*(eta(pin(0,4))+eta(pin(0,5)))
         phis=0.5d0*(phi(pin(0,4))+phi(pin(0,5)))
         goto 110
      endif
      etas=(di(2)*eta(pin(0,4))+di(3)*eta(pin(0,5)))/ets
      phis=(di(2)*phi(pin(0,4))+di(3)*phi(pin(0,5)))/ets
      
 210  p(0,4)=pin(0,4)
      p(1,4)=pin(1,4)
      p(2,4)=pin(2,4)
      p(3,4)=pin(3,4)

      return

c njets=3 just pass through
 300  njets=3
      do i=3,5
         do j=0,3
            p(j,i)=pin(j,i)
         enddo
      enddo
c
      return
      end
