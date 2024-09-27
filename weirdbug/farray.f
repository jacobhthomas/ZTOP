      program testfarray

      integer sets, xpdfmax
      parameter (sets = 59, xpdfmax = 13*(sets-1)+6)
      double precision ar(-6:xpdfmax)

      integer i

      call fillar(ar)
      
c      write (*,*) (ar(i), i= -6, xpdfmax)
      call writearn(ar(-4))
      
      end

      subroutine fillar(ar)
      integer sets, xpdfmax
      parameter (sets = 59, xpdfmax = 13*(sets-1)+6)
      double precision ar(-6:xpdfmax)
      
      integer i
      do i=-6, xpdfmax
         ar(i) = i - ((i+6)/13)*13
      enddo

      return
      end
      
      subroutine writearn(ar)
      integer sets, xpdfmax
      parameter (sets = 59, xpdfmax = 13*(sets-1)+6)
c      double precision ar(-6:xpdfmax)
      double precision ar(-6:6)

      write (*,*) ar(-6)
      write (*,*) ar(-5)
      write (*,*) ar(-4)
      write (*,*) ar(4)
      write (*,*) ar(5)
      write (*,*) ar(6)

      return
      end
      
